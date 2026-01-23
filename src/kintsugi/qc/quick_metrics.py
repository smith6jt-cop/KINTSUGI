"""
Fast quality control metrics for stitched images.

This module provides lightweight QC checks that can be run inline during
the stitching pipeline without significant performance impact (~50ms per image).

These checks detect common issues:
- Saturation: BaSiC correction failures causing clipped values
- Low signal: Blank or failed acquisitions
- Tile grids: Visible seams from stitching failures

Usage:
    from kintsugi.qc import compute_quick_qc_metrics

    stitched = stitch_with_blending(...)
    qc = compute_quick_qc_metrics(stitched)

    if qc['issues']:
        print(f"QC issues: {qc['issues']}")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass
class QuickQCResult:
    """Result of quick QC metrics computation."""

    issues: list[str] = field(default_factory=list)
    severity: Literal["none", "moderate", "severe"] = "none"
    saturation_pct: float = 0.0
    mean: float = 0.0
    std_mean_ratio: float = 0.0

    @property
    def passed(self) -> bool:
        """Return True if no issues detected."""
        return len(self.issues) == 0

    def to_dict(self) -> dict:
        """Convert to dictionary for JSON serialization."""
        return {
            "issues": self.issues,
            "severity": self.severity,
            "saturation_pct": self.saturation_pct,
            "mean": self.mean,
            "std_mean_ratio": self.std_mean_ratio,
        }


def compute_quick_qc_metrics(
    img: np.ndarray,
    saturation_threshold: float = 50.0,
    low_signal_threshold: float = 100.0,
    tile_grid_threshold: float = 0.8,
) -> QuickQCResult:
    """
    Compute fast QC metrics for a single stitched plane.

    This function is designed for inline use during the stitching pipeline.
    It runs in approximately 50ms and detects common quality issues that
    would cause problems in downstream processing.

    Parameters
    ----------
    img : np.ndarray
        2D stitched image (typically uint16).
    saturation_threshold : float, default 50.0
        Percentage of pixels at max value to trigger saturation flag.
    low_signal_threshold : float, default 100.0
        Minimum mean intensity below which low_signal is flagged.
    tile_grid_threshold : float, default 0.8
        Std/mean ratio above which tile_grid pattern is flagged.

    Returns
    -------
    QuickQCResult
        Dataclass with issues list, severity, and metric values.

    Examples
    --------
    >>> from kintsugi.qc import compute_quick_qc_metrics
    >>> qc = compute_quick_qc_metrics(stitched_image)
    >>> if qc.issues:
    ...     print(f"Issues detected: {qc.issues}")
    ...     print(f"Severity: {qc.severity}")
    """
    issues = []

    # Determine max value for saturation check
    if img.dtype == np.uint16:
        max_val = 65535
    elif img.dtype == np.uint8:
        max_val = 255
    elif np.issubdtype(img.dtype, np.integer):
        max_val = np.iinfo(img.dtype).max
    else:
        max_val = img.max() if img.max() > 0 else 1.0

    # Saturation check: >X% pixels at near-max value
    saturation_pct = 100.0 * np.sum(img > 0.98 * max_val) / img.size
    if saturation_pct > saturation_threshold:
        issues.append(f"saturation ({saturation_pct:.0f}%)")

    # Low signal check
    mean_val = float(img.mean())
    if mean_val < low_signal_threshold:
        issues.append(f"low_signal (mean={mean_val:.0f})")

    # Tile grid indicator (high std/mean ratio)
    if mean_val > 0:
        std_mean_ratio = float(img.std() / mean_val)
        if std_mean_ratio > tile_grid_threshold:
            issues.append(f"tile_grid (std/mean={std_mean_ratio:.2f})")
    else:
        std_mean_ratio = 0.0

    # Determine severity
    if not issues:
        severity: Literal["none", "moderate", "severe"] = "none"
    elif saturation_pct > 80 or mean_val < 50:
        severity = "severe"
    else:
        severity = "moderate"

    return QuickQCResult(
        issues=issues,
        severity=severity,
        saturation_pct=saturation_pct,
        mean=mean_val,
        std_mean_ratio=std_mean_ratio,
    )


def average_with_neighbors(
    prev_plane: np.ndarray | None,
    next_plane: np.ndarray | None,
    dtype: np.dtype,
    shape: tuple[int, int] | None = None,
) -> np.ndarray:
    """
    Average two z-planes for replacement of a bad plane.

    For edge planes where only one neighbor exists, averages with
    a black (zeros) image as specified in the design.

    Parameters
    ----------
    prev_plane : np.ndarray or None
        Previous z-plane (z-1), or None if at z=0.
    next_plane : np.ndarray or None
        Next z-plane (z+1), or None if at z=max.
    dtype : np.dtype
        Output data type.
    shape : tuple of int, optional
        Shape for zeros array if both neighbors are None.
        Required if both prev_plane and next_plane are None.

    Returns
    -------
    np.ndarray
        Averaged plane in the specified dtype.

    Examples
    --------
    >>> # Interior plane: average both neighbors
    >>> replacement = average_with_neighbors(z_minus_1, z_plus_1, np.uint16)

    >>> # Edge plane z=0: average with zeros + z_above
    >>> replacement = average_with_neighbors(None, z_plus_1, np.uint16)

    >>> # Edge plane z=max: average with z_below + zeros
    >>> replacement = average_with_neighbors(z_minus_1, None, np.uint16)
    """
    # Determine shape from available plane
    if prev_plane is not None:
        ref_shape = prev_plane.shape
    elif next_plane is not None:
        ref_shape = next_plane.shape
    elif shape is not None:
        ref_shape = shape
    else:
        raise ValueError("At least one neighbor plane or shape must be provided")

    # Create zeros for missing neighbors
    zeros = np.zeros(ref_shape, dtype=np.float32)

    z_below = prev_plane.astype(np.float32) if prev_plane is not None else zeros
    z_above = next_plane.astype(np.float32) if next_plane is not None else zeros

    # Average
    result = (z_below + z_above) / 2.0

    # Convert to output dtype with proper clipping
    if np.issubdtype(dtype, np.integer):
        dtype_info = np.iinfo(dtype)
        result = np.clip(result, dtype_info.min, dtype_info.max)

    return result.astype(dtype)


def fix_bad_zplanes_in_zarr(
    zarr_store,
    cycle: int,
    channel: int,
    flagged_zplanes: list[int],
    n_zplanes: int,
    verbose: bool = True,
) -> list[dict]:
    """
    Fix flagged z-planes in a zarr store by replacing with neighbor averages.

    This function is designed to be called AFTER all z-planes for a channel
    have been processed, allowing access to neighbors for averaging.

    Parameters
    ----------
    zarr_store : KintsugiZarr
        The zarr store containing stitched images.
    cycle : int
        Cycle number.
    channel : int
        Channel number.
    flagged_zplanes : list of int
        List of z-plane indices (1-based) that need replacement.
    n_zplanes : int
        Total number of z-planes.
    verbose : bool, default True
        Print status messages.

    Returns
    -------
    list of dict
        List of fix records with 'zplane', 'action', and 'neighbors_used' keys.

    Examples
    --------
    >>> from kintsugi.qc import fix_bad_zplanes_in_zarr
    >>> fixes = fix_bad_zplanes_in_zarr(zarr_store, cycle=1, channel=3,
    ...                                  flagged_zplanes=[1, 7], n_zplanes=15)
    """

    fixes = []

    for z in flagged_zplanes:
        # Determine neighbors
        z_below = z - 1 if z > 1 else None
        z_above = z + 1 if z < n_zplanes else None

        # Read neighbor planes from zarr
        prev_plane = None
        next_plane = None

        if z_below is not None:
            try:
                prev_plane = zarr_store.read_stitched(cycle, channel, z_below)
            except Exception:
                prev_plane = None

        if z_above is not None:
            try:
                next_plane = zarr_store.read_stitched(cycle, channel, z_above)
            except Exception:
                next_plane = None

        # Read current plane to get dtype and shape
        current = zarr_store.read_stitched(cycle, channel, z)
        dtype = current.dtype

        # Compute replacement
        replacement = average_with_neighbors(prev_plane, next_plane, dtype, shape=current.shape)

        # Write back to zarr
        zarr_store.write_stitched(cycle, channel, z, replacement)

        # Record the fix
        neighbors_used = []
        if prev_plane is not None:
            neighbors_used.append(f"z{z_below}")
        else:
            neighbors_used.append("zeros")
        if next_plane is not None:
            neighbors_used.append(f"z{z_above}")
        else:
            neighbors_used.append("zeros")

        fix_record = {
            "zplane": z,
            "action": "replaced_with_neighbor_avg",
            "neighbors_used": neighbors_used,
        }
        fixes.append(fix_record)

        if verbose:
            print(f"  [QC FIX] Cyc{cycle} CH{channel} Z{z:02d}: replaced with avg({', '.join(neighbors_used)})")

    return fixes


def fix_bad_zplanes_in_tiff(
    stitch_dir: str,
    cycle: int,
    channel: int,
    flagged_zplanes: list[int],
    n_zplanes: int,
    verbose: bool = True,
) -> list[dict]:
    """
    Fix flagged z-planes in TIFF files by replacing with neighbor averages.

    Parameters
    ----------
    stitch_dir : str
        Path to the stitch output directory.
    cycle : int
        Cycle number.
    channel : int
        Channel number.
    flagged_zplanes : list of int
        List of z-plane indices (1-based) that need replacement.
    n_zplanes : int
        Total number of z-planes.
    verbose : bool, default True
        Print status messages.

    Returns
    -------
    list of dict
        List of fix records.
    """
    import os

    from skimage.io import imread, imsave

    fixes = []
    ch_dir = os.path.join(stitch_dir, f"cyc{cycle:02d}", f"CH{channel}")

    for z in flagged_zplanes:
        # Determine neighbors
        z_below = z - 1 if z > 1 else None
        z_above = z + 1 if z < n_zplanes else None

        # Read neighbor planes from TIFF
        prev_plane = None
        next_plane = None

        if z_below is not None:
            prev_path = os.path.join(ch_dir, f"{z_below:02d}.tif")
            if os.path.exists(prev_path):
                prev_plane = imread(prev_path)

        if z_above is not None:
            next_path = os.path.join(ch_dir, f"{z_above:02d}.tif")
            if os.path.exists(next_path):
                next_plane = imread(next_path)

        # Read current plane to get dtype and shape
        current_path = os.path.join(ch_dir, f"{z:02d}.tif")
        current = imread(current_path)
        dtype = current.dtype

        # Compute replacement
        replacement = average_with_neighbors(prev_plane, next_plane, dtype, shape=current.shape)

        # Write back to TIFF
        imsave(current_path, replacement, check_contrast=False)

        # Record the fix
        neighbors_used = []
        if prev_plane is not None:
            neighbors_used.append(f"z{z_below}")
        else:
            neighbors_used.append("zeros")
        if next_plane is not None:
            neighbors_used.append(f"z{z_above}")
        else:
            neighbors_used.append("zeros")

        fix_record = {
            "zplane": z,
            "action": "replaced_with_neighbor_avg",
            "neighbors_used": neighbors_used,
        }
        fixes.append(fix_record)

        if verbose:
            print(f"  [QC FIX] Cyc{cycle} CH{channel} Z{z:02d}: replaced with avg({', '.join(neighbors_used)})")

    return fixes
