"""
Post-stitch tile boundary quality check.

Detects tile boundary misalignment in stitched images by measuring the
LOSS OF SHARPNESS in the overlap region. A properly stitched boundary blends
matching tissue features smoothly, preserving local sharpness. A misaligned
boundary causes GHOSTING — features from the two tiles appear in slightly
different positions and blend to a blurred result.

Metric: ``sharpness_ratio = Laplacian_variance(overlap) /
Laplacian_variance(adjacent_non_overlap)``. A ratio near 1.0 indicates a
healthy boundary. A ratio below ~0.7 indicates significant ghosting from
misalignment.

Typical usage from workflow/scripts/stitch.py::

    from kintsugi.qc.boundary_check import (
        check_boundary_quality,
        select_alternate_zplanes,
    )

    passed, flagged, scores = check_boundary_quality(
        stitched, result_df, tile_h, tile_w, threshold=0.7,
    )
    if not passed:
        for alt_z in select_alternate_zplanes(ref_zplane, n_zplanes):
            # re-run stitch_images on alt_z, re-check, pick best model
            ...
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
import pandas as pd
from scipy.ndimage import laplace


@dataclass
class BoundaryInfo:
    """Geometry of a single tile boundary in the stitched image."""

    orientation: Literal["vertical", "horizontal"]
    # Overlap region in stitched coordinates:
    #   vertical:   overlap_start/end are x-coordinates; extent is y-range
    #   horizontal: overlap_start/end are y-coordinates; extent is x-range
    overlap_start: int
    overlap_end: int
    extent_start: int
    extent_end: int
    tile_a: tuple[int, int]    # (row, col) of one tile
    tile_b: tuple[int, int]    # (row, col) of the adjacent tile

    @property
    def coord(self) -> int:
        """Midline of the overlap region (for reporting/plotting)."""
        return (self.overlap_start + self.overlap_end) // 2


def compute_boundary_positions(
    result_df: pd.DataFrame,
    tile_h: int,
    tile_w: int,
) -> list[BoundaryInfo]:
    """Derive tile boundary positions and overlap extents in stitched coords.

    Parameters
    ----------
    result_df : pd.DataFrame
        Stitch model DataFrame with columns ``row``, ``col``, ``y_pos``,
        ``x_pos``.
    tile_h, tile_w : int
        Tile height and width in pixels.

    Returns
    -------
    list[BoundaryInfo]
        One entry per tile boundary. Sorted by (orientation, tile_a).
    """
    if "row" not in result_df.columns or "col" not in result_df.columns:
        raise ValueError("result_df must have 'row' and 'col' columns")

    df = result_df.copy()
    df["x_pos2"] = df["x_pos"] - df["x_pos"].min()
    df["y_pos2"] = df["y_pos"] - df["y_pos"].min()

    by_rc = {(int(r.row), int(r.col)): r for r in df.itertuples(index=False)}

    boundaries: list[BoundaryInfo] = []

    for (row, col), tile in by_rc.items():
        # Right neighbor → vertical boundary
        right = by_rc.get((row, col + 1))
        if right is not None:
            x_ov_start = int(round(right.x_pos2))
            x_ov_end = int(round(tile.x_pos2 + tile_w))
            if x_ov_end > x_ov_start:
                y_start = int(round(max(tile.y_pos2, right.y_pos2)))
                y_end = int(round(min(tile.y_pos2 + tile_h, right.y_pos2 + tile_h)))
                if y_end > y_start:
                    boundaries.append(BoundaryInfo(
                        orientation="vertical",
                        overlap_start=x_ov_start,
                        overlap_end=x_ov_end,
                        extent_start=y_start,
                        extent_end=y_end,
                        tile_a=(row, col),
                        tile_b=(row, col + 1),
                    ))

        # Bottom neighbor → horizontal boundary
        below = by_rc.get((row + 1, col))
        if below is not None:
            y_ov_start = int(round(below.y_pos2))
            y_ov_end = int(round(tile.y_pos2 + tile_h))
            if y_ov_end > y_ov_start:
                x_start = int(round(max(tile.x_pos2, below.x_pos2)))
                x_end = int(round(min(tile.x_pos2 + tile_w, below.x_pos2 + tile_w)))
                if x_end > x_start:
                    boundaries.append(BoundaryInfo(
                        orientation="horizontal",
                        overlap_start=y_ov_start,
                        overlap_end=y_ov_end,
                        extent_start=x_start,
                        extent_end=x_end,
                        tile_a=(row, col),
                        tile_b=(row + 1, col),
                    ))

    boundaries.sort(key=lambda b: (b.orientation, b.tile_a))
    return boundaries


def _sharpness(img: np.ndarray) -> float:
    """Laplacian variance — a standard focus/sharpness metric.

    High for images with sharp structure, low for blurred/ghosted regions.
    """
    if img.size == 0:
        return 0.0
    lap = laplace(img.astype(np.float32, copy=False))
    return float(lap.var())


def measure_boundary_discontinuity(
    stitched: np.ndarray,
    boundaries: list[BoundaryInfo],
    **_ignored,
) -> dict[int, float]:
    """Measure sharpness ratio (overlap vs adjacent) at each boundary.

    For each tile boundary:
      1. Extract the overlap region (where both tiles contribute via blending)
      2. Extract two adjacent same-size regions just outside the overlap
         (where only one tile contributes)
      3. Compute Laplacian variance (sharpness) of each
      4. Score = sharpness(overlap) / mean(sharpness(adjacent))

    Interpretation:
      - score ≈ 1.0: overlap is as sharp as the interior → well aligned
      - score < 0.7: overlap is noticeably blurred/ghosted → misaligned
      - score < 0.5: severe ghosting — obvious visible artifact

    The metric is inherently scale-invariant: it normalizes by the local
    tissue sharpness, so dense vs sparse tissue gives comparable scores
    and detection is driven purely by blending degradation.

    Parameters
    ----------
    stitched : np.ndarray
        The stitched image (2D, any dtype).
    boundaries : list[BoundaryInfo]
        Output from :func:`compute_boundary_positions`.

    Returns
    -------
    dict[int, float]
        Mapping from boundary index to sharpness ratio. NaN if the overlap
        or adjacent regions are out of bounds or empty.
    """
    if stitched.ndim != 2:
        raise ValueError(f"stitched must be 2D, got shape {stitched.shape}")

    img = stitched.astype(np.float32, copy=False)
    h, w = img.shape

    scores: dict[int, float] = {}
    for idx, b in enumerate(boundaries):
        if b.orientation == "vertical":
            x_ov0 = max(0, b.overlap_start)
            x_ov1 = min(w, b.overlap_end)
            y0 = max(0, b.extent_start)
            y1 = min(h, b.extent_end)
            ov_w = x_ov1 - x_ov0
            if ov_w <= 0 or y1 - y0 <= 0:
                scores[idx] = float("nan")
                continue

            overlap = img[y0:y1, x_ov0:x_ov1]

            # Adjacent non-overlap regions (same width as overlap)
            left_x0 = max(0, x_ov0 - ov_w)
            left_x1 = x_ov0
            right_x0 = x_ov1
            right_x1 = min(w, x_ov1 + ov_w)
            left_adj = img[y0:y1, left_x0:left_x1]
            right_adj = img[y0:y1, right_x0:right_x1]
        else:  # horizontal
            y_ov0 = max(0, b.overlap_start)
            y_ov1 = min(h, b.overlap_end)
            x0 = max(0, b.extent_start)
            x1 = min(w, b.extent_end)
            ov_h = y_ov1 - y_ov0
            if ov_h <= 0 or x1 - x0 <= 0:
                scores[idx] = float("nan")
                continue

            overlap = img[y_ov0:y_ov1, x0:x1]

            top_y0 = max(0, y_ov0 - ov_h)
            top_y1 = y_ov0
            bot_y0 = y_ov1
            bot_y1 = min(h, y_ov1 + ov_h)
            left_adj = img[top_y0:top_y1, x0:x1]
            right_adj = img[bot_y0:bot_y1, x0:x1]

        # Require both adjacent regions to have at least half the overlap area
        min_adj_area = 0.5 * overlap.size
        if left_adj.size < min_adj_area and right_adj.size < min_adj_area:
            scores[idx] = float("nan")
            continue

        s_ov = _sharpness(overlap)
        s_adjs: list[float] = []
        if left_adj.size >= min_adj_area:
            s_adjs.append(_sharpness(left_adj))
        if right_adj.size >= min_adj_area:
            s_adjs.append(_sharpness(right_adj))
        if not s_adjs:
            scores[idx] = float("nan")
            continue

        s_adj = float(np.mean(s_adjs))
        if s_adj < 1.0:
            # Adjacent regions have no structure (background) — can't judge
            scores[idx] = float("nan")
            continue

        scores[idx] = s_ov / s_adj

    return scores


def check_boundary_quality(
    stitched: np.ndarray,
    result_df: pd.DataFrame,
    tile_h: int,
    tile_w: int,
    threshold: float = 0.7,
    **_ignored,
) -> tuple[bool, list[dict], dict[int, float]]:
    """Check tile boundaries for misalignment via sharpness degradation.

    A boundary is flagged when its overlap region is noticeably less sharp
    than the adjacent non-overlap regions, indicating that blending is
    averaging out ghosted (misaligned) features.

    Parameters
    ----------
    stitched : np.ndarray
        Stitched 2D image.
    result_df : pd.DataFrame
        Stitch model DataFrame (see :func:`compute_boundary_positions`).
    tile_h, tile_w : int
        Tile dimensions in pixels.
    threshold : float
        Sharpness ratio BELOW which a boundary is flagged. Default 0.7
        (overlap has < 70% of adjacent sharpness → significant ghosting).
        Conservative default catches clearly-visible misalignments without
        false alarms from natural blending-induced smoothing.

    Returns
    -------
    passed : bool
        True if no boundaries are flagged.
    flagged : list[dict]
        Flagged-boundary descriptors sorted by ascending score (worst first).
        Keys: ``index``, ``orientation``, ``coord``, ``tile_a``, ``tile_b``,
        ``score``.
    scores : dict[int, float]
        Per-boundary sharpness ratios (includes non-flagged boundaries).
    """
    boundaries = compute_boundary_positions(result_df, tile_h, tile_w)
    scores = measure_boundary_discontinuity(stitched, boundaries)

    flagged: list[dict] = []
    for idx, b in enumerate(boundaries):
        score = scores.get(idx, float("nan"))
        if not np.isnan(score) and score < threshold:
            flagged.append({
                "index": idx,
                "orientation": b.orientation,
                "coord": b.coord,
                "tile_a": b.tile_a,
                "tile_b": b.tile_b,
                "score": score,
            })

    flagged.sort(key=lambda d: d["score"])  # worst (lowest) first
    return (len(flagged) == 0), flagged, scores


def select_alternate_zplanes(
    ref_zplane: int,
    n_zplanes: int,
    max_alternates: int = 4,
) -> list[int]:
    """Choose alternate z-planes to try for the stitch model.

    The selection prioritizes z-planes spread across the stack to maximize
    tissue-feature diversity — if the reference z-plane has weak features at
    a boundary, another z-plane is likely to have stronger features there.

    Z-planes are 1-indexed (matching the on-disk ``_Z{N:03d}_`` naming).

    Parameters
    ----------
    ref_zplane : int
        The reference z-plane that was already tried (1-indexed).
    n_zplanes : int
        Total number of z-planes (1..n_zplanes inclusive).
    max_alternates : int
        Maximum number of alternates to return.

    Returns
    -------
    list[int]
        Ordered candidate z-planes (closest-spread first, then extremes).
        Never includes ``ref_zplane``. Never returns duplicates.
    """
    if n_zplanes < 2:
        return []

    quarter = max(1, n_zplanes // 4)
    candidates = [
        ref_zplane - quarter,
        ref_zplane + quarter,
        1,
        n_zplanes,
    ]

    seen: set[int] = {ref_zplane}
    out: list[int] = []
    for z in candidates:
        if z < 1 or z > n_zplanes:
            continue
        if z in seen:
            continue
        seen.add(z)
        out.append(z)
        if len(out) >= max_alternates:
            break

    return out
