"""
Autofluorescence Subtraction Core Algorithm

This module implements the core autofluorescence subtraction algorithm from
KINTSUGI's signal isolation workflow. The algorithm removes autofluorescent
signal (typically from collagen, elastin, lipofuscin) by subtracting a scaled
blank channel from the signal channel.

Algorithm Overview:
1. Clip blank channel to remove noise below threshold
2. Scale blank and subtract minimum of (signal, scaled_blank)
3. Optionally smooth low/high intensity regions
4. Optionally erode the result to clean up edges

The key insight is using min(signal, blank*scale) which prevents
over-subtraction in regions where signal > autofluorescence.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
from scipy import ndimage
from skimage import morphology
from skimage.filters import threshold_otsu
from skimage.measure import label, regionprops

logger = logging.getLogger("kintsugi.signal.autofluorescence")


def subtract_autofluorescence(
    signal: np.ndarray,
    blank: np.ndarray,
    blank_clip_factor: int = 0,
    blank_scale_factor: float = 1.0,
    smooth_low: bool = False,
    low_size: int = 2,
    low_percentile: int = 60,
    smooth_high: bool = False,
    high_size: int = 2,
    high_percentile: int = 90,
    erosion: int = 0,
) -> np.ndarray:
    """
    Subtract autofluorescence from signal using blank channel.

    This is a pure NumPy implementation of the KINTSUGI autofluorescence
    subtraction algorithm, optimized for speed while maintaining the exact
    behavior of the original.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image (2D array, typically uint16)
    blank : np.ndarray
        Blank/autofluorescence channel image (same shape as signal)
    blank_clip_factor : int
        Pixel values in blank below this are set to 0. Higher values
        remove more background noise from the blank. Default: 0
    blank_scale_factor : float
        Scale factor for blank channel before subtraction. Values > 1
        subtract more aggressively. Typical range: 0.5-2.0. Default: 1.0
    smooth_low : bool
        Apply uniform filter to low-intensity regions. Default: False
    low_size : int
        Size of uniform filter for low regions. Default: 2
    low_percentile : int
        Percentile threshold for low regions (0-100). Default: 60
    smooth_high : bool
        Apply uniform filter to high-intensity regions. Default: False
    high_size : int
        Size of uniform filter for high regions. Default: 2
    high_percentile : int
        Percentile threshold for high regions (0-100). Default: 90
    erosion : int
        Erosion disk radius for final cleanup. 0 = no erosion. Default: 0

    Returns
    -------
    np.ndarray
        Subtracted signal image (uint16)

    Notes
    -----
    The algorithm works as follows:

    1. **Blank Clipping**: Values below blank_clip_factor are zeroed:
       ```
       blank_clipped = clip(blank, blank_clip_factor, max)
       blank_clipped[blank_clipped <= blank_clip_factor] = 0
       ```

    2. **Scaled Subtraction**: Key step that prevents over-subtraction:
       ```
       result = signal - min(signal, blank_clipped * scale_factor)
       ```
       Using min() ensures we never subtract more than the signal value.

    3. **Percentile Smoothing** (optional): Smooths noisy regions:
       - Low regions: pixels < low_percentile get uniform filtered
       - High regions: pixels > high_percentile get uniform filtered

    4. **Erosion** (optional): Shrinks signal regions to remove edge artifacts.

    Examples
    --------
    >>> # Basic subtraction
    >>> result = subtract_autofluorescence(cd3_image, blank_image,
    ...     blank_clip_factor=5000, blank_scale_factor=1.2)

    >>> # With smoothing for noisy data
    >>> result = subtract_autofluorescence(cd3_image, blank_image,
    ...     blank_clip_factor=5000, blank_scale_factor=1.2,
    ...     smooth_low=True, low_size=3, low_percentile=50)
    """
    # Validate inputs
    if signal.shape != blank.shape:
        raise ValueError(f"Signal shape {signal.shape} != blank shape {blank.shape}")

    # Ensure consistent dtype for computation
    signal_float = signal.astype(np.float32)
    blank_float = blank.astype(np.float32)

    # Step 1: Clip blank channel
    blank_clipped = np.clip(blank_float, blank_clip_factor, blank_float.max())
    blank_clipped[blank_clipped <= blank_clip_factor] = 0

    # Step 2: Scaled subtraction using min() to prevent over-subtraction
    scaled_blank = blank_clipped * blank_scale_factor
    subtracted = signal_float - np.minimum(signal_float, scaled_blank)

    # Step 3a: Smooth low-intensity regions
    if smooth_low and low_percentile > 0:
        low_threshold = np.percentile(signal_float.ravel(), low_percentile)
        low_mask = signal_float < low_threshold
        # Dilate mask slightly for better coverage
        low_mask = morphology.binary_dilation(low_mask, morphology.disk(1))
        smoothed_low = ndimage.uniform_filter(subtracted, size=low_size)
        subtracted = np.where(low_mask, smoothed_low, subtracted)

    # Step 3b: Smooth high-intensity regions
    if smooth_high and high_percentile < 100:
        high_threshold = np.percentile(signal_float.ravel(), high_percentile)
        high_mask = signal_float > high_threshold
        high_mask = morphology.binary_dilation(high_mask, morphology.disk(1))
        smoothed_high = ndimage.uniform_filter(subtracted, size=high_size)
        subtracted = np.where(high_mask, smoothed_high, subtracted)

    # Step 4: Erosion cleanup
    if erosion > 0:
        signal_mask = subtracted > 0
        eroded_mask = morphology.binary_erosion(signal_mask, morphology.disk(erosion))
        subtracted = np.where(eroded_mask, subtracted, 0)

    # Ensure non-negative and convert to uint16
    subtracted = np.clip(subtracted, 0, 65535)
    return subtracted.astype(np.uint16)


def subtract_autofluorescence_dask(
    signal,  # dask array
    blank,   # dask array
    blank_clip_factor: int = 0,
    blank_scale_factor: float = 1.0,
    smooth_low: bool = False,
    low_size: int = 2,
    low_percentile: int = 60,
    smooth_high: bool = False,
    high_size: int = 2,
    high_percentile: int = 90,
    erosion: int = 0,
):
    """
    Dask-compatible version of autofluorescence subtraction.

    This version works with dask arrays for memory-efficient processing
    of large images. Mirrors the original Kutils.ini_params behavior.

    Parameters
    ----------
    signal : dask.array.Array
        Signal channel as dask array
    blank : dask.array.Array
        Blank channel as dask array
    (other parameters same as subtract_autofluorescence)

    Returns
    -------
    dask.array.Array
        Subtracted signal as dask array (uint16)
    """
    import dask.array as da

    try:
        import dask_image.ndfilters
        import dask_image.ndmorph
        HAS_DASK_IMAGE = True
    except ImportError:
        HAS_DASK_IMAGE = False
        logger.warning("dask-image not available, using numpy fallback")

    # Ensure consistent dtype
    signal_float = signal.astype(np.float32)
    blank_float = blank.astype(np.float32)

    # Step 1: Clip blank channel
    blank_clipped = da.clip(blank_float, blank_clip_factor, blank_float.max())
    blank_clipped = da.where(blank_clipped <= blank_clip_factor, 0, blank_clipped)

    # Step 2: Scaled subtraction
    scaled_blank = blank_clipped * blank_scale_factor
    subtracted = signal_float - da.minimum(signal_float, scaled_blank)

    # Step 3a: Smooth low-intensity regions
    if smooth_low and low_percentile > 0 and HAS_DASK_IMAGE:
        low_threshold = da.percentile(signal_float.ravel(), low_percentile)
        low_mask = da.where(signal_float < low_threshold, 1, 0)
        low_mask = dask_image.ndmorph.binary_dilation(low_mask, morphology.disk(1))
        smoothed_low = dask_image.ndfilters.uniform_filter(subtracted, size=low_size)
        subtracted = da.where(low_mask, smoothed_low, subtracted)

    # Step 3b: Smooth high-intensity regions
    if smooth_high and high_percentile < 100 and HAS_DASK_IMAGE:
        high_threshold = da.percentile(signal_float.ravel(), high_percentile)
        high_mask = da.where(signal_float > high_threshold, 1, 0)
        high_mask = dask_image.ndmorph.binary_dilation(high_mask, morphology.disk(1))
        smoothed_high = dask_image.ndfilters.uniform_filter(subtracted, size=high_size)
        subtracted = da.where(high_mask, smoothed_high, subtracted)

    # Step 4: Erosion cleanup
    if erosion > 0 and HAS_DASK_IMAGE:
        signal_mask = da.where(subtracted > 0, 1, 0)
        eroded_mask = dask_image.ndmorph.binary_erosion(signal_mask, morphology.disk(erosion))
        subtracted = da.where(eroded_mask, subtracted, 0)

    # Ensure non-negative and convert to uint16
    subtracted = da.clip(subtracted, 0, 65535)
    return subtracted.astype(np.uint16)


def analyze_for_subtraction(
    signal: np.ndarray,
    blank: np.ndarray,
    tissue_type: Optional[str] = None,
    marker_name: Optional[str] = None,
) -> Dict:
    """
    Analyze signal and blank images to suggest optimal subtraction parameters.

    This function examines the statistical properties of both images to
    recommend parameters that will effectively remove autofluorescence
    while preserving true signal.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image
    blank : np.ndarray
        Blank/autofluorescence channel image
    tissue_type : str, optional
        Tissue type for tissue-specific recommendations
    marker_name : str, optional
        Marker name for marker-specific recommendations

    Returns
    -------
    dict
        Suggested parameters with confidence scores:
        - blank_clip_factor: Suggested clip threshold
        - blank_scale_factor: Suggested scale factor
        - smooth_low: Whether to enable low smoothing
        - low_size: Suggested low filter size
        - low_percentile: Suggested low percentile
        - smooth_high: Whether to enable high smoothing
        - high_size: Suggested high filter size
        - high_percentile: Suggested high percentile
        - erosion: Suggested erosion radius
        - confidence: Overall confidence in suggestions (0-1)
        - analysis: Detailed analysis results

    Examples
    --------
    >>> params = analyze_for_subtraction(cd3_image, blank_image,
    ...     tissue_type='tonsil', marker_name='CD3')
    >>> print(f"Suggested scale factor: {params['blank_scale_factor']}")
    >>> print(f"Confidence: {params['confidence']:.2f}")
    """
    # Compute statistics
    signal_stats = _compute_image_stats(signal)
    blank_stats = _compute_image_stats(blank)

    # Analyze correlation between signal and blank
    correlation = _compute_correlation(signal, blank)

    # Estimate autofluorescence contribution
    af_contribution = _estimate_af_contribution(signal, blank)

    # Determine blank_clip_factor
    # Use percentile-based approach: clip at noise floor of blank
    blank_noise_floor = np.percentile(blank, 5)
    blank_p50 = np.percentile(blank, 50)

    # Suggested clip: between noise floor and median, biased toward removing noise
    suggested_clip = int(blank_noise_floor + (blank_p50 - blank_noise_floor) * 0.3)

    # Determine blank_scale_factor
    # Based on ratio of signal/blank intensities in overlapping regions
    # and the correlation between them
    if correlation > 0.7:
        # High correlation: autofluorescence is significant
        # Use higher scale factor
        suggested_scale = 1.0 + (correlation - 0.7) * 1.5  # 1.0 - 1.45
    elif correlation > 0.4:
        # Moderate correlation
        suggested_scale = 0.8 + (correlation - 0.4) * 0.67  # 0.8 - 1.0
    else:
        # Low correlation: less autofluorescence
        suggested_scale = 0.5 + correlation * 0.75  # 0.5 - 0.8

    # Analyze noise characteristics for smoothing recommendations
    signal_noise = _estimate_noise_level(signal)
    blank_noise = _estimate_noise_level(blank)

    # Recommend smoothing if noise is high
    noise_ratio = signal_noise / max(signal_stats['mean'], 1)
    suggest_smooth_low = noise_ratio > 0.15
    suggest_smooth_high = blank_stats['max'] > signal_stats['p95']

    # Determine erosion based on edge characteristics
    suggest_erosion = 1 if correlation > 0.6 else 0

    # Calculate confidence based on analysis quality
    confidence = _calculate_suggestion_confidence(
        signal_stats, blank_stats, correlation, af_contribution
    )

    # Tissue/marker specific adjustments
    if tissue_type and marker_name:
        suggested_clip, suggested_scale = _apply_tissue_marker_adjustments(
            suggested_clip, suggested_scale, tissue_type, marker_name
        )

    return {
        'blank_clip_factor': max(0, suggested_clip),
        'blank_scale_factor': round(max(0.1, min(3.0, suggested_scale)), 2),
        'smooth_low': suggest_smooth_low,
        'low_size': 2 if suggest_smooth_low else 1,
        'low_percentile': 60,
        'smooth_high': suggest_smooth_high,
        'high_size': 2 if suggest_smooth_high else 1,
        'high_percentile': 90,
        'erosion': suggest_erosion,
        'confidence': round(confidence, 3),
        'analysis': {
            'signal_stats': signal_stats,
            'blank_stats': blank_stats,
            'correlation': round(correlation, 4),
            'af_contribution': round(af_contribution, 4),
            'signal_noise_ratio': round(noise_ratio, 4),
        }
    }


def suggest_blank_channel(
    signal: np.ndarray,
    blank_channels: Dict[str, np.ndarray],
) -> Tuple[str, float]:
    """
    Suggest the best blank channel for a given signal channel.

    Analyzes correlation and overlap between signal and multiple blank
    channels to determine which blank best captures the autofluorescence.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image
    blank_channels : dict
        Dictionary of {channel_name: image_array} for available blanks

    Returns
    -------
    tuple
        (best_blank_name, correlation_score)

    Examples
    --------
    >>> blanks = {'Blank1a': blank1a, 'Blank1b': blank1b, 'Blank13c': blank13c}
    >>> best_blank, score = suggest_blank_channel(cd3_image, blanks)
    >>> print(f"Use {best_blank} (correlation: {score:.3f})")
    """
    best_blank = None
    best_score = -1

    for name, blank in blank_channels.items():
        if blank.shape != signal.shape:
            continue

        # Compute correlation
        correlation = _compute_correlation(signal, blank)

        # Higher correlation = better match for autofluorescence
        if correlation > best_score:
            best_score = correlation
            best_blank = name

    return best_blank, best_score


def compute_subtraction_quality(
    original_signal: np.ndarray,
    subtracted_signal: np.ndarray,
    blank: np.ndarray,
) -> Dict:
    """
    Compute quality metrics for a subtraction result.

    Evaluates how well the subtraction removed autofluorescence while
    preserving true signal.

    Parameters
    ----------
    original_signal : np.ndarray
        Original signal image before subtraction
    subtracted_signal : np.ndarray
        Signal image after subtraction
    blank : np.ndarray
        Blank channel used for subtraction

    Returns
    -------
    dict
        Quality metrics:
        - snr_improvement: Change in signal-to-noise ratio
        - af_removal: Estimated autofluorescence removal (0-1)
        - signal_preservation: Estimated true signal preserved (0-1)
        - residual_correlation: Correlation with blank after subtraction
        - quality_score: Overall quality score (0-1)
    """
    # SNR improvement
    orig_snr = _compute_snr(original_signal)
    sub_snr = _compute_snr(subtracted_signal)
    snr_improvement = sub_snr / max(orig_snr, 0.001) - 1

    # Residual correlation with blank (lower is better)
    residual_corr = _compute_correlation(subtracted_signal, blank)

    # Autofluorescence removal estimate
    orig_corr = _compute_correlation(original_signal, blank)
    af_removal = 1 - (residual_corr / max(orig_corr, 0.001))
    af_removal = max(0, min(1, af_removal))

    # Signal preservation (how much non-zero signal remains)
    orig_nonzero = np.sum(original_signal > 0)
    sub_nonzero = np.sum(subtracted_signal > 0)
    signal_preservation = sub_nonzero / max(orig_nonzero, 1)

    # Overall quality score
    quality_score = (
        0.3 * min(1, max(0, snr_improvement + 1) / 2) +  # SNR improvement
        0.3 * af_removal +  # AF removal
        0.2 * signal_preservation +  # Signal preservation
        0.2 * (1 - residual_corr)  # Low residual correlation
    )

    return {
        'snr_improvement': round(snr_improvement, 4),
        'af_removal': round(af_removal, 4),
        'signal_preservation': round(signal_preservation, 4),
        'residual_correlation': round(residual_corr, 4),
        'quality_score': round(quality_score, 4),
    }


# ============================================================================
# Internal helper functions
# ============================================================================

def _compute_image_stats(image: np.ndarray) -> Dict:
    """Compute basic statistics for an image."""
    flat = image.ravel()
    return {
        'min': float(np.min(flat)),
        'max': float(np.max(flat)),
        'mean': float(np.mean(flat)),
        'std': float(np.std(flat)),
        'p5': float(np.percentile(flat, 5)),
        'p25': float(np.percentile(flat, 25)),
        'p50': float(np.percentile(flat, 50)),
        'p75': float(np.percentile(flat, 75)),
        'p95': float(np.percentile(flat, 95)),
        'nonzero_fraction': float(np.sum(flat > 0) / len(flat)),
    }


def _compute_correlation(img1: np.ndarray, img2: np.ndarray) -> float:
    """Compute Pearson correlation between two images."""
    flat1 = img1.ravel().astype(np.float64)
    flat2 = img2.ravel().astype(np.float64)

    # Remove zero regions from both
    mask = (flat1 > 0) | (flat2 > 0)
    flat1 = flat1[mask]
    flat2 = flat2[mask]

    if len(flat1) < 100:
        return 0.0

    # Pearson correlation
    mean1, mean2 = np.mean(flat1), np.mean(flat2)
    std1, std2 = np.std(flat1), np.std(flat2)

    if std1 < 1e-6 or std2 < 1e-6:
        return 0.0

    correlation = np.mean((flat1 - mean1) * (flat2 - mean2)) / (std1 * std2)
    return max(-1, min(1, correlation))


def _estimate_af_contribution(signal: np.ndarray, blank: np.ndarray) -> float:
    """Estimate fraction of signal that is autofluorescence."""
    # Look at regions where blank is high
    blank_threshold = np.percentile(blank, 75)
    high_af_mask = blank > blank_threshold

    if np.sum(high_af_mask) < 100:
        return 0.0

    # Compare signal in high-AF regions vs low-AF regions
    low_af_mask = blank < np.percentile(blank, 25)

    signal_in_high_af = np.mean(signal[high_af_mask])
    signal_in_low_af = np.mean(signal[low_af_mask])

    if signal_in_low_af < 1:
        return 0.0

    # Excess signal in high-AF regions suggests AF contribution
    af_contribution = (signal_in_high_af - signal_in_low_af) / signal_in_high_af
    return max(0, min(1, af_contribution))


def _estimate_noise_level(image: np.ndarray) -> float:
    """Estimate noise level using MAD (median absolute deviation)."""
    from scipy.ndimage import laplacian

    # Use Laplacian to isolate noise
    lap = laplacian(image.astype(np.float64))
    mad = np.median(np.abs(lap - np.median(lap)))
    return mad * 1.4826  # Scale to std estimate


def _compute_snr(image: np.ndarray) -> float:
    """Compute signal-to-noise ratio."""
    # Signal: mean of top 10% pixels
    signal = np.mean(np.sort(image.ravel())[-int(image.size * 0.1):])

    # Noise: std of bottom 50% pixels
    noise = np.std(np.sort(image.ravel())[:int(image.size * 0.5)])

    if noise < 1:
        return signal
    return signal / noise


def _calculate_suggestion_confidence(
    signal_stats: Dict,
    blank_stats: Dict,
    correlation: float,
    af_contribution: float,
) -> float:
    """Calculate confidence in parameter suggestions."""
    confidence = 0.5  # Base confidence

    # Higher confidence if clear autofluorescence pattern
    if correlation > 0.5:
        confidence += 0.2
    if af_contribution > 0.2:
        confidence += 0.15

    # Higher confidence if good signal range
    if signal_stats['nonzero_fraction'] > 0.1:
        confidence += 0.1

    # Lower confidence for edge cases
    if signal_stats['max'] < 1000:
        confidence -= 0.1
    if blank_stats['max'] < 500:
        confidence -= 0.1

    return max(0.1, min(1.0, confidence))


def _apply_tissue_marker_adjustments(
    clip: int,
    scale: float,
    tissue_type: str,
    marker_name: str,
) -> Tuple[int, float]:
    """Apply tissue and marker-specific adjustments."""
    tissue_type = tissue_type.lower()
    marker_name = marker_name.upper()

    # Tissue-specific adjustments
    if 'tonsil' in tissue_type or 'lymph' in tissue_type:
        # Lymphoid tissue often has high autofluorescence
        scale *= 1.1
    elif 'skin' in tissue_type:
        # Skin has collagen autofluorescence
        scale *= 1.2
    elif 'liver' in tissue_type or 'kidney' in tissue_type:
        # These tissues have lipofuscin
        clip = int(clip * 1.2)

    # Marker-specific adjustments
    # Dim markers need gentler subtraction
    dim_markers = {'FOXP3', 'CD163', 'CD11C', 'CD1C'}
    if marker_name in dim_markers:
        scale *= 0.9

    # Bright markers can tolerate more aggressive subtraction
    bright_markers = {'CD3E', 'CD3', 'CD20', 'DAPI', 'PANCK'}
    if marker_name in bright_markers:
        scale *= 1.05

    return clip, scale
