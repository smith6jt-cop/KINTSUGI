"""
SMO (Silver Mountain Operator) Background Estimation Wrapper

Parameter-free background estimation for fluorescence microscopy.
SMO estimates the background distribution from a single image using only
two parameters (sigma and size), which are robust across a wide range of values.

This provides an alternative to blank-channel subtraction for channels where
blank channels are unavailable, unreliable, or where blank subtraction
produces artifacts.

Reference: Silber et al., github.com/maurosilber/smo
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger("kintsugi.signal.smo")


def estimate_background_smo(
    image: np.ndarray,
    kernel_size: int = 7,
    sigma: float = 1.0,
    return_background: bool = False,
) -> dict:
    """
    Estimate and subtract background using the Silver Mountain Operator.

    SMO recovers an unbiased estimation of the background distribution from a
    single image. Unlike blank subtraction, it requires no reference channel.

    Parameters
    ----------
    image : np.ndarray
        2D single-channel image (uint16 or float).
    kernel_size : int
        SMO averaging window size (default 7). Robust across range 5-15 for
        most microscopy images. Larger values estimate smoother backgrounds.
    sigma : float
        Gaussian kernel sigma for SMO (default 1.0).
    return_background : bool
        If True, also return the estimated background image.

    Returns
    -------
    dict with:
        corrected : np.ndarray (background-subtracted image, same dtype as input)
        background_mean : float (estimated mean background level)
        background_std : float (estimated background noise level)
        kernel_size_used : int
        background : np.ndarray or None (if return_background=True)
    """
    try:
        from smo import SMO
    except ImportError:
        raise ImportError(
            "SMO package not installed. Install with: pip install smo\n"
            "Or install the optimize extra: pip install kintsugi[optimize]\n"
            "SMO provides parameter-free background estimation for fluorescence microscopy.\n"
            "Reference: Silber et al., github.com/maurosilber/smo"
        )

    original_dtype = image.dtype
    img_f = image.astype(np.float64)

    # Create SMO instance with image shape
    smo = SMO(sigma=sigma, size=kernel_size, shape=img_f.shape)

    # bg_corrected returns background-subtracted image
    corrected_ma = smo.bg_corrected(img_f)
    corrected = np.asarray(corrected_ma, dtype=np.float64)

    # Get background distribution statistics
    bg_rv = smo.bg_rv(img_f)
    bg_mean_val = float(bg_rv.mean())
    bg_std_val = float(bg_rv.std())

    # Estimate per-pixel background as original minus corrected
    bg_image = img_f - corrected if return_background else None

    # Clip negative values (background can exceed signal in some regions)
    corrected = np.clip(corrected, 0, None)

    # Restore original dtype
    if np.issubdtype(original_dtype, np.integer):
        corrected = corrected.astype(original_dtype)

    return {
        "corrected": corrected,
        "background_mean": bg_mean_val,
        "background_std": bg_std_val,
        "kernel_size_used": kernel_size,
        "background": bg_image if return_background else None,
    }


def auto_select_kernel_size(
    image: np.ndarray,
    test_sizes: list | None = None,
) -> dict:
    """
    Automatically select optimal SMO kernel size by testing several values
    and selecting the one that maximizes foreground-background contrast
    in the corrected image.

    Parameters
    ----------
    image : np.ndarray
        2D single-channel image.
    test_sizes : list[int], optional
        Kernel sizes to test. Default: [3, 5, 7, 9, 11, 13].

    Returns
    -------
    dict with:
        optimal_kernel_size : int
        scores : dict (kernel_size -> score)
    """
    if test_sizes is None:
        test_sizes = [3, 5, 7, 9, 11, 13]

    from skimage.filters import threshold_otsu

    try:
        thresh = threshold_otsu(image)
    except ValueError:
        thresh = np.mean(image.astype(np.float64))

    bg_mask = image <= thresh

    scores = {}
    for ks in test_sizes:
        result = estimate_background_smo(image, kernel_size=ks)
        corrected = result["corrected"].astype(np.float64)

        # Score: high foreground-background contrast relative to background noise
        bg_var = np.var(corrected[bg_mask]) if bg_mask.any() else 1e10
        fg_mean = np.mean(corrected[~bg_mask]) if (~bg_mask).any() else 0
        bg_mean = np.mean(corrected[bg_mask]) if bg_mask.any() else 0
        contrast = (fg_mean - bg_mean) / (bg_var**0.5 + 1e-10)
        scores[ks] = float(contrast)

    optimal = max(scores, key=scores.get)
    return {"optimal_kernel_size": optimal, "scores": scores}
