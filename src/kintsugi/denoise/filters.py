"""
Traditional Denoising Filters

Pure Python implementations of classical denoising algorithms
optimized for fluorescence microscopy images.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Union

import numpy as np
from scipy import ndimage
from scipy.ndimage import gaussian_filter, median_filter

if TYPE_CHECKING:
    import dask.array

logger = logging.getLogger("kintsugi.denoise.filters")

# Type alias for array-like
ArrayLike = Union[np.ndarray, "dask.array.Array"]


def _to_numpy(arr: ArrayLike) -> np.ndarray:
    """Convert array to numpy, computing if dask."""
    if hasattr(arr, "compute"):
        return arr.compute()
    return np.asarray(arr)


def _preserve_dtype(func):
    """Decorator to preserve input dtype in output."""

    def wrapper(image: ArrayLike, *args, **kwargs):
        original_dtype = image.dtype
        result = func(image, *args, **kwargs)
        # Clip and convert back to original dtype
        if np.issubdtype(original_dtype, np.integer):
            info = np.iinfo(original_dtype)
            result = np.clip(result, info.min, info.max)
        return result.astype(original_dtype)

    return wrapper


def estimate_noise_level(
    image: ArrayLike,
    method: str = "mad",
) -> float:
    """
    Estimate the noise level in an image.

    Parameters
    ----------
    image : array-like
        Input image (2D or 3D)
    method : str
        Estimation method:
        - "mad": Median Absolute Deviation (robust to outliers)
        - "std_background": Std of low-intensity regions
        - "laplacian": High-frequency content estimation
        - "wavelet": Wavelet-based estimation

    Returns
    -------
    float
        Estimated noise standard deviation
    """
    img = _to_numpy(image).astype(np.float64)

    if method == "mad":
        # Robust estimator using Laplacian for high-frequency content
        laplacian = ndimage.laplace(img)
        # MAD * 1.4826 approximates std for Gaussian noise
        noise_std = np.median(np.abs(laplacian - np.median(laplacian))) * 1.4826
        # Adjust for Laplacian operator
        noise_std /= np.sqrt(6)

    elif method == "std_background":
        # Use low-intensity regions (bottom 20%)
        threshold = np.percentile(img, 20)
        background = img[img < threshold]
        if len(background) > 100:
            noise_std = float(np.std(background))
        else:
            noise_std = float(np.std(img) * 0.1)

    elif method == "laplacian":
        # Estimate from Laplacian variance
        laplacian = ndimage.laplace(img)
        noise_std = np.sqrt(np.var(laplacian) / 6.0)

    elif method == "wavelet":
        # Simplified wavelet-based estimation using Haar-like filter
        # Approximates detail coefficients
        h_detail = img[:, 1:] - img[:, :-1]
        v_detail = img[1:, :] - img[:-1, :]
        # MAD of detail coefficients
        noise_std = (
            np.median(
                np.abs(
                    np.concatenate(
                        [
                            h_detail.ravel(),
                            v_detail.ravel(),
                        ]
                    )
                )
            )
            * 1.4826
            / np.sqrt(2)
        )

    else:
        raise ValueError(f"Unknown method: {method}")

    return float(noise_std)


@_preserve_dtype
def denoise_median(
    image: ArrayLike,
    size: int = 3,
    mode: str = "reflect",
) -> np.ndarray:
    """
    Apply median filter denoising.

    Effective for salt-and-pepper noise while preserving edges.

    Parameters
    ----------
    image : array-like
        Input image
    size : int
        Filter footprint size (size x size)
    mode : str
        How to handle boundaries ('reflect', 'constant', 'nearest', 'wrap')

    Returns
    -------
    np.ndarray
        Denoised image
    """
    img = _to_numpy(image).astype(np.float64)
    return median_filter(img, size=size, mode=mode)


@_preserve_dtype
def denoise_gaussian(
    image: ArrayLike,
    sigma: float = 1.0,
    truncate: float = 4.0,
) -> np.ndarray:
    """
    Apply Gaussian filter denoising.

    Simple smoothing filter, effective but may blur edges.

    Parameters
    ----------
    image : array-like
        Input image
    sigma : float
        Standard deviation for Gaussian kernel
    truncate : float
        Truncate the filter at this many standard deviations

    Returns
    -------
    np.ndarray
        Denoised image
    """
    img = _to_numpy(image).astype(np.float64)
    return gaussian_filter(img, sigma=sigma, truncate=truncate)


@_preserve_dtype
def denoise_bilateral(
    image: ArrayLike,
    sigma_spatial: float = 2.0,
    sigma_range: float | None = None,
    win_size: int | None = None,
) -> np.ndarray:
    """
    Apply bilateral filter denoising.

    Edge-preserving smoothing that considers both spatial distance
    and intensity difference.

    Parameters
    ----------
    image : array-like
        Input image
    sigma_spatial : float
        Spatial Gaussian standard deviation
    sigma_range : float, optional
        Range Gaussian standard deviation. If None, estimated from image.
    win_size : int, optional
        Window size. If None, computed from sigma_spatial.

    Returns
    -------
    np.ndarray
        Denoised image
    """
    try:
        from skimage.restoration import denoise_bilateral as skimage_bilateral

        img = _to_numpy(image).astype(np.float64)

        # Normalize to [0, 1] for skimage
        img_min, img_max = img.min(), img.max()
        if img_max > img_min:
            img_norm = (img - img_min) / (img_max - img_min)
        else:
            return img

        if sigma_range is None:
            sigma_range = 0.1

        if win_size is None:
            win_size = int(4 * sigma_spatial + 1) | 1  # Ensure odd

        result = skimage_bilateral(
            img_norm,
            sigma_color=sigma_range,
            sigma_spatial=sigma_spatial,
            win_size=win_size,
            channel_axis=None,
        )

        # Denormalize
        return result * (img_max - img_min) + img_min

    except ImportError:
        logger.warning("skimage not available, falling back to Gaussian")
        return denoise_gaussian(image, sigma=sigma_spatial)


@_preserve_dtype
def denoise_nlm(
    image: ArrayLike,
    patch_size: int = 5,
    patch_distance: int = 6,
    h: float | None = None,
    fast_mode: bool = True,
) -> np.ndarray:
    """
    Apply Non-Local Means denoising.

    Averages pixels based on patch similarity across the image.
    Excellent for preserving texture while removing noise.

    Parameters
    ----------
    image : array-like
        Input image
    patch_size : int
        Size of patches for similarity comparison
    patch_distance : int
        Maximum distance for patch search
    h : float, optional
        Filtering strength. If None, estimated from noise level.
    fast_mode : bool
        Use faster but less accurate algorithm

    Returns
    -------
    np.ndarray
        Denoised image
    """
    try:
        from skimage.restoration import denoise_nl_means

        img = _to_numpy(image).astype(np.float64)

        # Normalize for skimage
        img_min, img_max = img.min(), img.max()
        if img_max > img_min:
            img_norm = (img - img_min) / (img_max - img_min)
        else:
            return img

        if h is None:
            # Estimate h from noise level
            noise_std = estimate_noise_level(img_norm, method="mad")
            h = 1.15 * noise_std

        result = denoise_nl_means(
            img_norm,
            patch_size=patch_size,
            patch_distance=patch_distance,
            h=h,
            fast_mode=fast_mode,
            channel_axis=None,
        )

        return result * (img_max - img_min) + img_min

    except ImportError:
        logger.warning("skimage not available, falling back to bilateral")
        return denoise_bilateral(image, sigma_spatial=patch_size / 2)


def adaptive_denoise(
    image: ArrayLike,
    strength: str = "auto",
    preserve_edges: bool = True,
) -> np.ndarray:
    """
    Automatically select and apply appropriate denoising.

    Analyzes the image to determine noise level and selects
    the best denoising strategy.

    Parameters
    ----------
    image : array-like
        Input image
    strength : str
        Denoising strength: "light", "medium", "strong", or "auto"
    preserve_edges : bool
        Prioritize edge preservation

    Returns
    -------
    np.ndarray
        Denoised image with same dtype as input
    """
    img = _to_numpy(image)
    original_dtype = img.dtype

    # Estimate noise level
    noise_std = estimate_noise_level(img, method="mad")

    # Normalize noise estimate relative to image range
    img_range = float(np.percentile(img, 99) - np.percentile(img, 1))
    if img_range > 0:
        relative_noise = noise_std / img_range
    else:
        relative_noise = 0.01

    # Determine strength if auto
    if strength == "auto":
        if relative_noise < 0.02:
            strength = "light"
        elif relative_noise < 0.1:
            strength = "medium"
        else:
            strength = "strong"

    logger.info(
        f"Adaptive denoise: noise={noise_std:.2f}, relative={relative_noise:.4f}, strength={strength}"
    )

    # Apply denoising based on strength and edge preference
    if strength == "light":
        if preserve_edges:
            result = denoise_bilateral(img, sigma_spatial=1.0)
        else:
            result = denoise_gaussian(img, sigma=0.5)

    elif strength == "medium":
        if preserve_edges:
            result = denoise_nlm(img, patch_size=5, patch_distance=6)
        else:
            result = denoise_median(img, size=3)

    else:  # strong
        if preserve_edges:
            # Two-stage denoising
            result = denoise_nlm(img, patch_size=7, patch_distance=11)
            result = denoise_bilateral(result, sigma_spatial=1.5)
        else:
            result = denoise_median(img, size=5)
            result = denoise_gaussian(result, sigma=1.0)

    # Convert back to original dtype
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        result = np.clip(result, info.min, info.max)

    return result.astype(original_dtype)


def denoise_dask(
    image: dask.array.Array,
    method: str = "median",
    overlap: int = 16,
    **kwargs,
) -> dask.array.Array:
    """
    Apply denoising to a Dask array with proper overlap handling.

    Parameters
    ----------
    image : dask.array.Array
        Input Dask array
    method : str
        Denoising method: "median", "gaussian", "bilateral", "nlm", "adaptive"
    overlap : int
        Overlap size for block processing
    **kwargs
        Additional arguments passed to the denoising function

    Returns
    -------
    dask.array.Array
        Denoised Dask array
    """
    import dask.array as da

    method_map = {
        "median": denoise_median,
        "gaussian": denoise_gaussian,
        "bilateral": denoise_bilateral,
        "nlm": denoise_nlm,
        "adaptive": adaptive_denoise,
    }

    if method not in method_map:
        raise ValueError(f"Unknown method: {method}. Choose from {list(method_map.keys())}")

    denoise_func = method_map[method]

    def apply_denoise(block):
        return denoise_func(block, **kwargs)

    # Use map_overlap for proper boundary handling
    result = da.map_overlap(
        apply_denoise,
        image,
        depth=overlap,
        boundary="reflect",
        dtype=image.dtype,
    )

    return result
