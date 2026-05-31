"""Tests for kintsugi.qc.boundary_check."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
from scipy.ndimage import gaussian_filter

from kintsugi.qc.boundary_check import (
    check_boundary_quality,
    compute_boundary_positions,
    measure_boundary_discontinuity,
    select_alternate_zplanes,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_result_df(
    n_rows: int,
    n_cols: int,
    tile_h: int,
    tile_w: int,
    overlap: float,
) -> pd.DataFrame:
    """Build a result_df for an n_rows x n_cols tile grid with uniform overlap."""
    step_y = int(tile_h * (1 - overlap))
    step_x = int(tile_w * (1 - overlap))
    rows = []
    for r in range(n_rows):
        for c in range(n_cols):
            rows.append(
                {
                    "row": r,
                    "col": c,
                    "y_pos": r * step_y,
                    "x_pos": c * step_x,
                }
            )
    return pd.DataFrame(rows)


def _build_sharp_image(shape: tuple[int, int], seed: int = 0) -> np.ndarray:
    """Random image with fine-scale features (high Laplacian variance)."""
    rng = np.random.default_rng(seed)
    img = rng.integers(100, 500, size=shape, dtype=np.uint16)
    return img


# ---------------------------------------------------------------------------
# compute_boundary_positions
# ---------------------------------------------------------------------------


def test_compute_boundary_positions_3x3_count():
    """3x3 grid → 6 vertical + 6 horizontal = 12 boundaries."""
    df = _make_result_df(3, 3, tile_h=100, tile_w=100, overlap=0.3)
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=100)
    vertical = [b for b in boundaries if b.orientation == "vertical"]
    horizontal = [b for b in boundaries if b.orientation == "horizontal"]
    assert len(vertical) == 6
    assert len(horizontal) == 6


def test_compute_boundary_positions_overlap_extent():
    """Verify vertical boundary overlap extents."""
    # 1x2 grid, 100x100 tiles, 30% overlap → step_x=70
    # col 1 starts at x=70, col 0 ends at x=100 → overlap is [70, 100]
    df = _make_result_df(1, 2, tile_h=100, tile_w=100, overlap=0.3)
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=100)
    assert len(boundaries) == 1
    b = boundaries[0]
    assert b.orientation == "vertical"
    assert b.overlap_start == 70
    assert b.overlap_end == 100
    assert b.coord == 85  # midline
    assert b.extent_start == 0
    assert b.extent_end == 100
    assert b.tile_a == (0, 0)
    assert b.tile_b == (0, 1)


def test_compute_boundary_positions_horizontal_overlap():
    df = _make_result_df(2, 1, tile_h=100, tile_w=100, overlap=0.3)
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=100)
    assert len(boundaries) == 1
    b = boundaries[0]
    assert b.orientation == "horizontal"
    assert b.overlap_start == 70
    assert b.overlap_end == 100


def test_compute_boundary_positions_rejects_bad_df():
    df = pd.DataFrame({"foo": [1], "bar": [2]})
    with pytest.raises(ValueError, match="row.*col"):
        compute_boundary_positions(df, tile_h=100, tile_w=100)


# ---------------------------------------------------------------------------
# measure_boundary_discontinuity (sharpness ratio)
# ---------------------------------------------------------------------------


def test_measure_sharpness_ratio_aligned():
    """When overlap has the same sharpness as adjacent → ratio ≈ 1.0."""
    # Build a long strip where sharpness is uniform (same random pattern
    # across the whole image) — overlap and adjacent have identical stats
    rng = np.random.default_rng(42)
    # Full stitched image: 100 tall, 230 wide (covers overlap + adjacent)
    img = rng.integers(100, 500, size=(100, 230), dtype=np.uint16)

    # Boundary: overlap spans x=100..130 (30px wide)
    df = pd.DataFrame(
        [
            {"row": 0, "col": 0, "y_pos": 0, "x_pos": 0},
            {"row": 0, "col": 1, "y_pos": 0, "x_pos": 100},
        ]
    )
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=130)
    scores = measure_boundary_discontinuity(img, boundaries)
    valid = [s for s in scores.values() if not np.isnan(s)]
    assert len(valid) == 1
    # Uniform texture → overlap sharpness matches adjacent
    assert 0.7 <= valid[0] <= 1.3, f"aligned ratio should be near 1, got {valid[0]}"


def test_measure_sharpness_ratio_blurred_overlap():
    """When overlap is blurred (ghosting) → ratio substantially less than 1."""
    rng = np.random.default_rng(7)
    img = rng.integers(100, 500, size=(100, 230)).astype(np.float32)
    # Blur ONLY the overlap region (x=100..130) — simulates ghosting
    blurred_overlap = gaussian_filter(img[:, 100:130], sigma=2.0)
    img_mixed = img.copy()
    img_mixed[:, 100:130] = blurred_overlap
    img_mixed = img_mixed.astype(np.uint16)

    df = pd.DataFrame(
        [
            {"row": 0, "col": 0, "y_pos": 0, "x_pos": 0},
            {"row": 0, "col": 1, "y_pos": 0, "x_pos": 100},
        ]
    )
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=130)
    scores = measure_boundary_discontinuity(img_mixed, boundaries)
    valid = [s for s in scores.values() if not np.isnan(s)]
    assert len(valid) == 1
    # Blurred overlap → sharpness ratio well below 1
    assert valid[0] < 0.5, f"blurred overlap should be < 0.5, got {valid[0]}"


def test_measure_sharpness_background_returns_nan():
    """All-background overlap and adjacent → NaN."""
    img = np.zeros((100, 230), dtype=np.uint16)
    df = pd.DataFrame(
        [
            {"row": 0, "col": 0, "y_pos": 0, "x_pos": 0},
            {"row": 0, "col": 1, "y_pos": 0, "x_pos": 100},
        ]
    )
    boundaries = compute_boundary_positions(df, tile_h=100, tile_w=130)
    scores = measure_boundary_discontinuity(img, boundaries)
    assert np.isnan(scores[0])


# ---------------------------------------------------------------------------
# check_boundary_quality
# ---------------------------------------------------------------------------


def test_check_boundary_quality_passes_on_aligned():
    """Uniform texture → all boundaries pass."""
    rng = np.random.default_rng(1)
    img = rng.integers(100, 500, size=(100, 230), dtype=np.uint16)
    df = pd.DataFrame(
        [
            {"row": 0, "col": 0, "y_pos": 0, "x_pos": 0},
            {"row": 0, "col": 1, "y_pos": 0, "x_pos": 100},
        ]
    )
    passed, flagged, scores = check_boundary_quality(
        img,
        df,
        tile_h=100,
        tile_w=130,
        threshold=0.7,
    )
    assert passed, f"aligned should pass, flagged={flagged}"


def test_check_boundary_quality_flags_blurred():
    """Ghosted overlap → flagged."""
    rng = np.random.default_rng(7)
    img = rng.integers(100, 500, size=(100, 230)).astype(np.float32)
    blurred = gaussian_filter(img[:, 100:130], sigma=2.0)
    img_mixed = img.copy()
    img_mixed[:, 100:130] = blurred
    img_mixed = img_mixed.astype(np.uint16)

    df = pd.DataFrame(
        [
            {"row": 0, "col": 0, "y_pos": 0, "x_pos": 0},
            {"row": 0, "col": 1, "y_pos": 0, "x_pos": 100},
        ]
    )
    passed, flagged, scores = check_boundary_quality(
        img_mixed,
        df,
        tile_h=100,
        tile_w=130,
        threshold=0.7,
    )
    assert not passed
    assert len(flagged) == 1
    assert flagged[0]["orientation"] == "vertical"
    assert flagged[0]["tile_a"] == (0, 0)
    assert flagged[0]["tile_b"] == (0, 1)


def test_check_boundary_quality_worst_first():
    """Flagged boundaries sorted by ascending score (worst first)."""
    # 3 tiles with 30% overlap:
    #   tile_w=100, step=70 → tile 0: x=0-100, tile 1: x=70-170, tile 2: x=140-240
    #   Overlap 0↔1: x=[70, 100], Overlap 1↔2: x=[140, 170]
    rng = np.random.default_rng(3)
    img = rng.integers(100, 500, size=(100, 240)).astype(np.float32)
    # Boundary 1 (x=70-100) heavily blurred
    img[:, 70:100] = gaussian_filter(img[:, 70:100], sigma=3.0)
    # Boundary 2 (x=140-170) moderately blurred
    img[:, 140:170] = gaussian_filter(img[:, 140:170], sigma=1.5)
    img = img.astype(np.uint16)

    df = _make_result_df(1, 3, tile_h=100, tile_w=100, overlap=0.3)
    passed, flagged, scores = check_boundary_quality(
        img,
        df,
        tile_h=100,
        tile_w=100,
        threshold=0.9,
    )
    assert not passed
    assert len(flagged) == 2, f"expected 2 flagged, got {flagged}; scores={scores}"
    # Worst (heavier blur) should come first
    assert flagged[0]["score"] < flagged[1]["score"]


# ---------------------------------------------------------------------------
# select_alternate_zplanes
# ---------------------------------------------------------------------------


def test_select_alternate_zplanes_13_from_7():
    alternates = select_alternate_zplanes(ref_zplane=7, n_zplanes=13, max_alternates=4)
    assert alternates == [4, 10, 1, 13]


def test_select_alternate_zplanes_20_from_10():
    alternates = select_alternate_zplanes(ref_zplane=10, n_zplanes=20, max_alternates=4)
    assert alternates == [5, 15, 1, 20]


def test_select_alternate_zplanes_respects_max():
    alternates = select_alternate_zplanes(ref_zplane=7, n_zplanes=13, max_alternates=2)
    assert len(alternates) == 2
    assert alternates == [4, 10]


def test_select_alternate_zplanes_small_stack():
    """4 z-planes, ref=2 → [1, 3, 4]."""
    alternates = select_alternate_zplanes(ref_zplane=2, n_zplanes=4, max_alternates=4)
    assert alternates == [1, 3, 4]


def test_select_alternate_zplanes_single_plane():
    assert select_alternate_zplanes(ref_zplane=1, n_zplanes=1) == []


def test_select_alternate_zplanes_never_includes_ref():
    for ref in range(1, 14):
        alternates = select_alternate_zplanes(ref_zplane=ref, n_zplanes=13)
        assert ref not in alternates
