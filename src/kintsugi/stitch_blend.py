"""
Sigmoid taper blending for artifact-free image stitching.

This module provides smooth blending of overlapping tiles using sigmoid taper
functions, adapted from the DeepTile library (https://github.com/arjunrajlaboratory/DeepTile).

The sigmoid taper approach eliminates visible seams without blurring underlying
structure, unlike median/gaussian filtering approaches.

Usage
-----
>>> from kintsugi.stitch_blend import stitch_with_blending, generate_taper
>>>
>>> # Option 1: Stitch with blending in one call
>>> stitched = stitch_with_blending(tiles, result_df, blend=True, sigma=10)
>>>
>>> # Option 2: Replace existing smooth_tile_borders_2d calls
>>> # Old: stitched = smooth_tile_borders_2d(stitched, result_df, tile_shape)
>>> # New: Use stitch_with_blending() instead of manual placement + smoothing

References
----------
- DeepTile: https://github.com/arjunrajlaboratory/DeepTile
- Sigmoid blending for seamless mosaics
"""

import numpy as np
import pandas as pd


def generate_taper(
    tile_size: tuple[int, int], overlap: tuple[float, float], sigma: float = 10.0
) -> np.ndarray:
    """
    Generate sigmoid taper for blending tile overlaps.

    Creates a 2D weight mask that smoothly transitions from 1 at tile center
    to 0 at edges, with the transition width controlled by sigma.

    Parameters
    ----------
    tile_size : tuple of int
        Size of each tile as (height, width).
    overlap : tuple of float
        Fractions of tile_size to use for overlap in each dimension (y, x).
        Values should be between 0 and 1 (e.g., 0.3 for 30% overlap).
    sigma : float, default 10.0
        Sigma bandwidth parameter controlling transition sharpness.
        - Smaller values (5-10): sharper transitions, better for high-contrast images
        - Larger values (15-30): smoother transitions, better for noisy images

    Returns
    -------
    taper : np.ndarray
        2D taper array of shape tile_size with values in [0, 1].
        Values are 1 at center and smoothly decrease to ~0 at edges.

    Examples
    --------
    >>> taper = generate_taper((512, 512), (0.3, 0.3), sigma=10)
    >>> taper.shape
    (512, 512)
    >>> taper[256, 256]  # Center value
    1.0
    """
    height, width = tile_size
    overlap_y, overlap_x = overlap

    # Create distance arrays from center
    x = np.arange(width)
    x = np.abs(x - np.median(x))
    y = np.arange(height)
    y = np.abs(y - np.median(y))

    # Sigmoid taper: 1 at center, ~0 at edges
    # The transition occurs at the non-overlap region boundary
    taper_x = 1 / (1 + np.exp((x - width * (1 - overlap_x) / 2) / sigma))
    taper_y = 1 / (1 + np.exp((y - height * (1 - overlap_y) / 2) / sigma))

    # Combine into 2D taper
    taper = taper_y[:, None] * taper_x

    return taper


def stitch_with_blending(
    tiles: np.ndarray,
    result_df: pd.DataFrame,
    blend: bool = True,
    sigma: float = 10.0,
    overlap_fraction: tuple[float, float] | None = None,
    output_dtype: np.dtype | None = None,
) -> np.ndarray:
    """
    Stitch tiles into a seamless mosaic using sigmoid taper blending.

    This function replaces the manual tile placement + smooth_tile_borders_2d()
    workflow with a mathematically smooth blending approach that eliminates
    visible seams without blurring underlying structure.

    Parameters
    ----------
    tiles : np.ndarray
        Array of tiles with shape (n_tiles, height, width).
    result_df : pd.DataFrame
        Stitching result DataFrame from stitch_images() containing
        'x_pos', 'y_pos' or 'x_pos2', 'y_pos2' columns with tile positions.
    blend : bool, default True
        If True, use sigmoid taper blending in overlap regions.
        If False, tiles are simply placed (last tile wins in overlaps).
    sigma : float, default 10.0
        Sigma bandwidth parameter for sigmoid taper.
        - Smaller (5-10): sharper transitions
        - Larger (15-30): smoother transitions
    overlap_fraction : tuple of float, optional
        Overlap fraction (y, x) for taper generation.
        If None, automatically estimated from tile positions.
    output_dtype : np.dtype, optional
        Output data type. If None, uses input dtype.

    Returns
    -------
    stitched : np.ndarray
        Stitched image with seamless blending.

    Notes
    -----
    The sigmoid taper approach:
    1. Creates a 2D weight mask for each tile (high at center, low at edges)
    2. Multiplies each tile by its weight mask
    3. Accumulates weighted tiles and weight masks separately
    4. Divides by accumulated weights to normalize

    This ensures smooth transitions without visible seams, even with
    intensity variations between tiles.

    Examples
    --------
    >>> from kintsugi.stitch_blend import stitch_with_blending
    >>> from Kstitch.stitching import stitch_images
    >>>
    >>> # Get tile positions
    >>> result_df, props = stitch_images(tiles, rows, cols, ...)
    >>>
    >>> # Stitch with blending (replaces manual placement + smooth_tile_borders_2d)
    >>> stitched = stitch_with_blending(tiles, result_df, blend=True, sigma=10)
    """
    tile_h, tile_w = tiles.shape[1], tiles.shape[2]

    # Determine position column names
    if "y_pos2" in result_df.columns and "x_pos2" in result_df.columns:
        y_col, x_col = "y_pos2", "x_pos2"
    elif "y_pos" in result_df.columns and "x_pos" in result_df.columns:
        y_col, x_col = "y_pos", "x_pos"
        # Normalize positions to start from 0
        result_df = result_df.copy()
        result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
        result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()
        y_col, x_col = "y_pos2", "x_pos2"
    else:
        raise ValueError("result_df must contain position columns (x_pos/y_pos or x_pos2/y_pos2)")

    # Calculate output size
    out_h = int(result_df[y_col].max() + tile_h)
    out_w = int(result_df[x_col].max() + tile_w)

    # Determine output dtype
    if output_dtype is None:
        output_dtype = tiles.dtype

    if not blend:
        # Simple placement (no blending)
        stitched = np.zeros((out_h, out_w), dtype=output_dtype)
        for i, row in result_df.iterrows():
            y0, x0 = int(row[y_col]), int(row[x_col])
            stitched[y0 : y0 + tile_h, x0 : x0 + tile_w] = tiles[i]
        return stitched

    # Estimate overlap fraction if not provided
    if overlap_fraction is None:
        overlap_fraction = _estimate_overlap_fraction(result_df, tile_h, tile_w, y_col, x_col)

    # Generate taper
    taper = generate_taper((tile_h, tile_w), overlap_fraction, sigma)

    # Accumulate weighted tiles and weights
    stitched = np.zeros((out_h, out_w), dtype=np.float64)
    weight_sum = np.zeros((out_h, out_w), dtype=np.float64)

    for i, row in result_df.iterrows():
        y0, x0 = int(row[y_col]), int(row[x_col])
        y1, x1 = y0 + tile_h, x0 + tile_w

        # Handle tiles that extend beyond output bounds (shouldn't happen normally)
        y1_clip = min(y1, out_h)
        x1_clip = min(x1, out_w)
        tile_h_clip = y1_clip - y0
        tile_w_clip = x1_clip - x0

        # Get tile and taper slices
        tile_slice = tiles[i, :tile_h_clip, :tile_w_clip]
        taper_slice = taper[:tile_h_clip, :tile_w_clip]

        # Accumulate weighted values
        stitched[y0:y1_clip, x0:x1_clip] += tile_slice.astype(np.float64) * taper_slice
        weight_sum[y0:y1_clip, x0:x1_clip] += taper_slice

    # Normalize by accumulated weights (avoid division by zero)
    stitched = stitched / (weight_sum + 1e-10)

    # Convert back to original dtype
    if np.issubdtype(output_dtype, np.integer):
        dtype_info = np.iinfo(output_dtype)
        stitched = np.clip(stitched, dtype_info.min, dtype_info.max)

    return stitched.astype(output_dtype)


def _estimate_overlap_fraction(
    result_df: pd.DataFrame, tile_h: int, tile_w: int, y_col: str, x_col: str
) -> tuple[float, float]:
    """
    Estimate overlap fraction from tile positions.

    Parameters
    ----------
    result_df : pd.DataFrame
        Stitching result DataFrame with tile positions.
    tile_h, tile_w : int
        Tile dimensions.
    y_col, x_col : str
        Column names for y and x positions.

    Returns
    -------
    overlap_fraction : tuple of float
        Estimated overlap fraction (y, x).
    """
    positions = result_df[[y_col, x_col]].values

    if len(positions) < 2:
        # Single tile, no overlap
        return (0.1, 0.1)

    # Find typical step sizes between adjacent tiles
    y_steps = []
    x_steps = []

    for i in range(len(positions)):
        for j in range(i + 1, len(positions)):
            dy = abs(positions[j, 0] - positions[i, 0])
            dx = abs(positions[j, 1] - positions[i, 1])

            # Check if tiles are adjacent (one dimension small, one dimension ~tile size)
            if dy < tile_h * 0.1 and tile_w * 0.5 < dx < tile_w * 1.5:
                x_steps.append(dx)
            elif dx < tile_w * 0.1 and tile_h * 0.5 < dy < tile_h * 1.5:
                y_steps.append(dy)

    # Calculate overlap from step sizes
    if y_steps:
        y_step = np.median(y_steps)
        overlap_y = max(0.05, min(0.5, 1 - y_step / tile_h))
    else:
        overlap_y = 0.3  # Default 30% overlap

    if x_steps:
        x_step = np.median(x_steps)
        overlap_x = max(0.05, min(0.5, 1 - x_step / tile_w))
    else:
        overlap_x = 0.3  # Default 30% overlap

    return (overlap_y, overlap_x)


def smooth_tile_borders_sigmoid(
    stitched_plane: np.ndarray,
    result_df: pd.DataFrame,
    tile_shape: tuple[int, int],
    sigma: float = 10.0,
    overlap_fraction: tuple[float, float] | None = None,
) -> np.ndarray:
    """
    Re-stitch an already-stitched image using sigmoid blending.

    This is a drop-in replacement for smooth_tile_borders_2d() that uses
    sigmoid taper blending instead of median/gaussian filtering.

    Note: This requires re-accessing the original tiles. For better efficiency,
    use stitch_with_blending() directly instead of stitching then smoothing.

    Parameters
    ----------
    stitched_plane : np.ndarray
        Already stitched 2D image (used only for shape reference).
    result_df : pd.DataFrame
        Stitching result DataFrame with tile positions.
    tile_shape : tuple of int
        Shape of individual tiles (height, width).
    sigma : float, default 10.0
        Sigma bandwidth for sigmoid taper.
    overlap_fraction : tuple of float, optional
        Overlap fraction (y, x). Auto-estimated if None.

    Returns
    -------
    stitched : np.ndarray
        Re-stitched image with sigmoid blending.

    Raises
    ------
    NotImplementedError
        This function requires access to original tiles.
        Use stitch_with_blending() instead.
    """
    raise NotImplementedError(
        "smooth_tile_borders_sigmoid requires access to original tiles.\n"
        "Use stitch_with_blending(tiles, result_df, blend=True) instead of:\n"
        "  1. Manual tile placement\n"
        "  2. smooth_tile_borders_2d() or smooth_tile_borders_sigmoid()\n"
        "\n"
        "This approach is both faster and produces better results."
    )


# Convenience function for common use case
def stitch_corrected_tiles(
    corrected: np.ndarray,
    result_df: pd.DataFrame,
    sigma: float = 10.0,
    overlap_percentage: float = 30.0,
) -> np.ndarray:
    """
    Stitch illumination-corrected tiles with sigmoid blending.

    Convenience wrapper for the common KINTSUGI workflow of stitching
    BaSiC-corrected tiles.

    Parameters
    ----------
    corrected : np.ndarray
        Corrected tiles with shape (n_tiles, height, width).
    result_df : pd.DataFrame
        Stitching result DataFrame from stitch_images().
    sigma : float, default 10.0
        Sigma bandwidth for blending.
    overlap_percentage : float, default 30.0
        Overlap percentage (0-100) for taper calculation.

    Returns
    -------
    stitched : np.ndarray
        Seamlessly stitched image.

    Examples
    --------
    >>> from kintsugi.stitch_blend import stitch_corrected_tiles
    >>>
    >>> # After BaSiC correction and stitch_images()
    >>> stitched = stitch_corrected_tiles(corrected, result_df)
    """
    overlap_fraction = (overlap_percentage / 100.0, overlap_percentage / 100.0)
    return stitch_with_blending(
        corrected,
        result_df,
        blend=True,
        sigma=sigma,
        overlap_fraction=overlap_fraction,
    )


# GPU-accelerated version for large images
def stitch_with_blending_gpu(
    tiles: np.ndarray,
    result_df: pd.DataFrame,
    sigma: float = 10.0,
    overlap_fraction: tuple[float, float] | None = None,
    device_id: int = 0,
) -> np.ndarray:
    """
    GPU-accelerated stitching with sigmoid taper blending.

    Uses CuPy for faster processing of large images. Falls back to CPU
    if CuPy is not available.

    Parameters
    ----------
    tiles : np.ndarray
        Array of tiles with shape (n_tiles, height, width).
    result_df : pd.DataFrame
        Stitching result DataFrame with tile positions.
    sigma : float, default 10.0
        Sigma bandwidth for sigmoid taper.
    overlap_fraction : tuple of float, optional
        Overlap fraction (y, x). Auto-estimated if None.
    device_id : int, default 0
        GPU device ID to use.

    Returns
    -------
    stitched : np.ndarray
        Stitched image (always returned as numpy array).
    """
    try:
        import cupy as cp

        cp.cuda.Device(device_id).use()
        HAS_CUPY = True
    except ImportError:
        HAS_CUPY = False

    if not HAS_CUPY:
        # Fallback to CPU
        return stitch_with_blending(
            tiles, result_df, blend=True, sigma=sigma, overlap_fraction=overlap_fraction
        )

    tile_h, tile_w = tiles.shape[1], tiles.shape[2]

    # Determine position columns
    if "y_pos2" in result_df.columns and "x_pos2" in result_df.columns:
        y_col, x_col = "y_pos2", "x_pos2"
    else:
        result_df = result_df.copy()
        result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
        result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()
        y_col, x_col = "y_pos2", "x_pos2"

    # Calculate output size
    out_h = int(result_df[y_col].max() + tile_h)
    out_w = int(result_df[x_col].max() + tile_w)

    # Estimate overlap if not provided
    if overlap_fraction is None:
        overlap_fraction = _estimate_overlap_fraction(result_df, tile_h, tile_w, y_col, x_col)

    # Generate taper on GPU
    taper_cpu = generate_taper((tile_h, tile_w), overlap_fraction, sigma)
    taper_gpu = cp.asarray(taper_cpu)

    # Allocate output arrays on GPU
    stitched_gpu = cp.zeros((out_h, out_w), dtype=cp.float64)
    weight_sum_gpu = cp.zeros((out_h, out_w), dtype=cp.float64)

    # Process tiles
    for i, row in result_df.iterrows():
        y0, x0 = int(row[y_col]), int(row[x_col])
        y1, x1 = y0 + tile_h, x0 + tile_w

        # Clip to bounds
        y1_clip = min(y1, out_h)
        x1_clip = min(x1, out_w)
        tile_h_clip = y1_clip - y0
        tile_w_clip = x1_clip - x0

        # Transfer tile to GPU and accumulate
        tile_gpu = cp.asarray(tiles[i, :tile_h_clip, :tile_w_clip], dtype=cp.float64)
        taper_slice = taper_gpu[:tile_h_clip, :tile_w_clip]

        stitched_gpu[y0:y1_clip, x0:x1_clip] += tile_gpu * taper_slice
        weight_sum_gpu[y0:y1_clip, x0:x1_clip] += taper_slice

        # Free tile memory
        del tile_gpu

    # Normalize
    stitched_gpu = stitched_gpu / (weight_sum_gpu + 1e-10)

    # Transfer back to CPU
    stitched = cp.asnumpy(stitched_gpu)

    # Clean up GPU memory
    del stitched_gpu, weight_sum_gpu, taper_gpu
    cp.get_default_memory_pool().free_all_blocks()

    # Convert to original dtype
    if np.issubdtype(tiles.dtype, np.integer):
        dtype_info = np.iinfo(tiles.dtype)
        stitched = np.clip(stitched, dtype_info.min, dtype_info.max)

    return stitched.astype(tiles.dtype)
