"""
Patch-Based Denoising Algorithms

Pure Python implementations of patch-based denoising methods
including a lightweight BM3D-inspired algorithm and patch similarity methods.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Union

import numpy as np
from scipy import ndimage
from scipy.fft import dct, idct

if TYPE_CHECKING:
    import dask.array

logger = logging.getLogger("kintsugi.denoise.patch_based")

# Type alias
ArrayLike = Union[np.ndarray, "dask.array.Array"]


def _to_numpy(arr: ArrayLike) -> np.ndarray:
    """Convert to numpy, computing if dask."""
    if hasattr(arr, "compute"):
        return arr.compute()
    return np.asarray(arr)


def _extract_patches(
    image: np.ndarray,
    patch_size: int,
    step: int,
) -> tuple[np.ndarray, list]:
    """Extract overlapping patches from image."""
    h, w = image.shape
    patches = []
    positions = []

    for y in range(0, h - patch_size + 1, step):
        for x in range(0, w - patch_size + 1, step):
            patch = image[y : y + patch_size, x : x + patch_size]
            patches.append(patch)
            positions.append((y, x))

    return np.array(patches), positions


def _reconstruct_from_patches(
    patches: np.ndarray,
    positions: list,
    image_shape: tuple[int, int],
    patch_size: int,
) -> np.ndarray:
    """Reconstruct image from overlapping patches with averaging."""
    h, w = image_shape
    output = np.zeros((h, w), dtype=np.float64)
    weights = np.zeros((h, w), dtype=np.float64)

    for patch, (y, x) in zip(patches, positions, strict=False):
        output[y : y + patch_size, x : x + patch_size] += patch
        weights[y : y + patch_size, x : x + patch_size] += 1

    weights = np.maximum(weights, 1)
    return output / weights


def _dct_2d(block: np.ndarray) -> np.ndarray:
    """2D DCT transform."""
    return dct(dct(block.T, norm="ortho").T, norm="ortho")


def _idct_2d(block: np.ndarray) -> np.ndarray:
    """Inverse 2D DCT transform."""
    return idct(idct(block.T, norm="ortho").T, norm="ortho")


def _hard_threshold(coeffs: np.ndarray, threshold: float) -> np.ndarray:
    """Apply hard thresholding to coefficients."""
    return np.where(np.abs(coeffs) < threshold, 0, coeffs)


def _wiener_threshold(
    coeffs: np.ndarray,
    estimate_coeffs: np.ndarray,
    noise_var: float,
) -> np.ndarray:
    """Apply Wiener filtering to coefficients."""
    signal_var = np.maximum(estimate_coeffs**2 - noise_var, 0)
    wiener_weights = signal_var / (signal_var + noise_var + 1e-8)
    return coeffs * wiener_weights


def denoise_bm3d_lite(
    image: ArrayLike,
    sigma: float | None = None,
    patch_size: int = 8,
    search_window: int = 21,
    max_matches: int = 16,
    step: int = 3,
    hard_threshold_factor: float = 2.7,
) -> np.ndarray:
    """
    Lightweight BM3D-inspired denoising.

    This is a simplified implementation of BM3D (Block-Matching and 3D filtering)
    that provides good denoising while being more efficient than full BM3D.

    The algorithm:
    1. Extract overlapping patches
    2. For each patch, find similar patches within search window
    3. Stack similar patches into 3D group
    4. Apply 3D transform (2D DCT + 1D Haar along stack)
    5. Threshold coefficients
    6. Inverse transform and aggregate

    Parameters
    ----------
    image : array-like
        Input image
    sigma : float, optional
        Noise standard deviation. If None, estimated automatically.
    patch_size : int
        Size of patches
    search_window : int
        Size of search window for finding similar patches
    max_matches : int
        Maximum number of similar patches to use
    step : int
        Step size for patch extraction
    hard_threshold_factor : float
        Factor for hard threshold (threshold = factor * sigma)

    Returns
    -------
    np.ndarray
        Denoised image
    """
    img = _to_numpy(image).astype(np.float64)
    original_dtype = image.dtype

    # Estimate noise if not provided
    if sigma is None:
        from kintsugi.denoise.filters import estimate_noise_level

        sigma = estimate_noise_level(img, method="mad")
        logger.info(f"Estimated noise sigma: {sigma:.2f}")

    h, w = img.shape
    half_window = search_window // 2
    threshold = hard_threshold_factor * sigma

    # Extract reference patches
    patches, positions = _extract_patches(img, patch_size, step)
    n_patches = len(patches)

    # Pre-compute DCT of all patches for efficient similarity search
    patch_dcts = np.array([_dct_2d(p) for p in patches])

    # Output patches
    output_patches = np.zeros_like(patches)
    patch_weights = np.zeros(n_patches)

    for i, (ref_patch, (ref_y, ref_x)) in enumerate(zip(patches, positions, strict=False)):
        # Find similar patches within search window
        similar_indices = [i]
        ref_dct = patch_dcts[i]

        for j, (_other_patch, (other_y, other_x)) in enumerate(
            zip(patches, positions, strict=False)
        ):
            if i == j:
                continue

            # Check if within search window
            if abs(other_y - ref_y) <= half_window and abs(other_x - ref_x) <= half_window:

                # Compute similarity in DCT domain (faster than spatial)
                diff = np.sum((patch_dcts[j] - ref_dct) ** 2)
                if diff < (patch_size**2) * (sigma**2) * 2:
                    similar_indices.append(j)

                    if len(similar_indices) >= max_matches:
                        break

        # Stack similar patches into 3D group
        group = patches[similar_indices]
        n_group = len(group)

        if n_group < 2:
            # Not enough similar patches - use single patch DCT filtering
            dct_coeffs = _dct_2d(ref_patch)
            filtered = _hard_threshold(dct_coeffs, threshold)
            output_patches[i] = _idct_2d(filtered)
            patch_weights[i] = 1.0
            continue

        # 3D transform: 2D DCT on each patch + 1D Haar along stack
        group_dct = np.array([_dct_2d(p) for p in group])

        # 1D Haar transform along stack dimension (simplified)
        group_3d = np.zeros_like(group_dct)
        for level in range(int(np.log2(n_group)) + 1):
            step_size = 2**level
            if step_size >= n_group:
                break
            for k in range(0, n_group - step_size, 2 * step_size):
                group_3d[k] = (group_dct[k] + group_dct[k + step_size]) / np.sqrt(2)
                if k + step_size < n_group:
                    group_3d[k + step_size] = (group_dct[k] - group_dct[k + step_size]) / np.sqrt(2)
            group_dct = group_3d.copy()

        # Hard threshold in transform domain
        group_filtered = _hard_threshold(group_3d, threshold)

        # Inverse 1D Haar transform
        group_ihaar = group_filtered.copy()
        for level in range(int(np.log2(n_group)), -1, -1):
            step_size = 2**level
            if step_size >= n_group:
                continue
            for k in range(0, n_group - step_size, 2 * step_size):
                low = group_ihaar[k]
                high = group_ihaar[k + step_size] if k + step_size < n_group else 0
                group_ihaar[k] = (low + high) / np.sqrt(2)
                if k + step_size < n_group:
                    group_ihaar[k + step_size] = (low - high) / np.sqrt(2)

        # Inverse 2D DCT and aggregate to output
        denoised_group = np.array([_idct_2d(p) for p in group_ihaar])

        # Weight by number of non-zero coefficients (sparsity)
        n_nonzero = np.sum(group_filtered != 0)
        weight = 1.0 / (1.0 + n_nonzero / (n_group * patch_size**2))

        # Accumulate to output patches
        for idx, patch_idx in enumerate(similar_indices):
            output_patches[patch_idx] += denoised_group[idx] * weight
            patch_weights[patch_idx] += weight

    # Normalize patch weights
    patch_weights = np.maximum(patch_weights, 1e-8)
    output_patches = output_patches / patch_weights[:, np.newaxis, np.newaxis]

    # Reconstruct image
    result = _reconstruct_from_patches(output_patches, positions, (h, w), patch_size)

    # Clip to valid range
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        result = np.clip(result, info.min, info.max)

    return result.astype(original_dtype)


def denoise_patch_similarity(
    image: ArrayLike,
    patch_size: int = 7,
    search_radius: int = 10,
    h_param: float | None = None,
    fast_mode: bool = True,
) -> np.ndarray:
    """
    Non-local patch similarity denoising.

    A pure Python implementation of non-local means denoising using
    patch similarity for weighting.

    Parameters
    ----------
    image : array-like
        Input image
    patch_size : int
        Size of comparison patches
    search_radius : int
        Radius for searching similar patches
    h_param : float, optional
        Filtering strength. If None, estimated from noise level.
    fast_mode : bool
        Use faster but less accurate algorithm

    Returns
    -------
    np.ndarray
        Denoised image
    """
    img = _to_numpy(image).astype(np.float64)
    original_dtype = image.dtype

    h, w = img.shape
    half_patch = patch_size // 2

    # Estimate h from noise if not provided
    if h_param is None:
        from kintsugi.denoise.filters import estimate_noise_level

        sigma = estimate_noise_level(img, method="mad")
        h_param = sigma * 1.2

    # Pad image
    padded = np.pad(img, half_patch + search_radius, mode="reflect")

    # Output image
    output = np.zeros_like(img)

    if fast_mode:
        # Faster version using vectorized operations
        # Pre-compute patch sums of squares for efficiency
        patch_sq_sum = ndimage.uniform_filter(padded**2, size=patch_size) * (patch_size**2)
        patch_sum = ndimage.uniform_filter(padded, size=patch_size) * (patch_size**2)

        offset = half_patch + search_radius

        for y in range(h):
            for x in range(w):
                py, px = y + offset, x + offset

                # Extract search window
                y_start = py - search_radius
                y_end = py + search_radius + 1
                x_start = px - search_radius
                x_end = px + search_radius + 1

                # Get reference patch statistics
                ref_sq_sum = patch_sq_sum[py, px]
                ref_sum = patch_sum[py, px]

                # Get search window statistics
                window_sq_sum = patch_sq_sum[y_start:y_end, x_start:x_end]
                window_sum = patch_sum[y_start:y_end, x_start:x_end]

                # Approximate distance using sum statistics
                # d^2 = sum((a-b)^2) = sum(a^2) + sum(b^2) - 2*sum(a*b)
                # Approximated as: ref_sq_sum + window_sq_sum - 2 * ref_sum * window_sum / n
                distances = ref_sq_sum + window_sq_sum - 2 * ref_sum * window_sum / (patch_size**2)
                distances = np.maximum(distances, 0)

                # Compute weights
                weights = np.exp(-distances / (h_param**2 * patch_size**2))

                # Weighted average of center pixels in search window
                center_values = padded[
                    y_start + half_patch : y_end + half_patch,
                    x_start + half_patch : x_end + half_patch,
                ]

                output[y, x] = np.sum(weights * center_values) / np.sum(weights)

    else:
        # More accurate but slower version
        offset = half_patch + search_radius

        for y in range(h):
            for x in range(w):
                py, px = y + offset, x + offset

                # Extract reference patch
                ref_patch = padded[
                    py - half_patch : py + half_patch + 1, px - half_patch : px + half_patch + 1
                ]

                weighted_sum = 0.0
                weight_sum = 0.0

                # Search window
                for dy in range(-search_radius, search_radius + 1):
                    for dx in range(-search_radius, search_radius + 1):
                        ny, nx = py + dy, px + dx

                        # Extract comparison patch
                        comp_patch = padded[
                            ny - half_patch : ny + half_patch + 1,
                            nx - half_patch : nx + half_patch + 1,
                        ]

                        # Compute distance
                        dist = np.sum((ref_patch - comp_patch) ** 2)

                        # Weight
                        weight = np.exp(-dist / (h_param**2 * patch_size**2))

                        weighted_sum += weight * padded[ny, nx]
                        weight_sum += weight

                output[y, x] = weighted_sum / weight_sum

    # Clip to valid range
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        output = np.clip(output, info.min, info.max)

    return output.astype(original_dtype)


def denoise_svd_patch(
    image: ArrayLike,
    patch_size: int = 8,
    rank: int | None = None,
    stride: int = 4,
) -> np.ndarray:
    """
    SVD-based patch denoising.

    Collects patches, performs SVD, and reconstructs using top singular values.

    Parameters
    ----------
    image : array-like
        Input image
    patch_size : int
        Size of patches
    rank : int, optional
        Number of singular values to keep. If None, determined automatically.
    stride : int
        Stride for patch extraction

    Returns
    -------
    np.ndarray
        Denoised image
    """
    img = _to_numpy(image).astype(np.float64)
    original_dtype = image.dtype

    h, w = img.shape

    # Extract patches
    patches, positions = _extract_patches(img, patch_size, stride)
    n_patches = len(patches)

    # Reshape patches to matrix (each patch is a row)
    patch_matrix = patches.reshape(n_patches, -1)

    # Center the data
    mean = patch_matrix.mean(axis=0)
    centered = patch_matrix - mean

    # SVD
    U, S, Vt = np.linalg.svd(centered, full_matrices=False)

    # Determine rank if not specified
    if rank is None:
        # Use energy-based selection
        energy = np.cumsum(S**2) / np.sum(S**2)
        rank = np.argmax(energy > 0.95) + 1
        rank = max(rank, 3)  # Keep at least 3 components

    logger.info(f"SVD patch denoising using rank {rank}")

    # Truncate
    S_truncated = S.copy()
    S_truncated[rank:] = 0

    # Reconstruct
    reconstructed = U @ np.diag(S_truncated) @ Vt + mean

    # Reshape back to patches
    output_patches = reconstructed.reshape(n_patches, patch_size, patch_size)

    # Reconstruct image
    result = _reconstruct_from_patches(output_patches, positions, (h, w), patch_size)

    # Clip
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        result = np.clip(result, info.min, info.max)

    return result.astype(original_dtype)
