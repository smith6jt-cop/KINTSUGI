"""
Signal Isolation Utility Functions

Low-level utility functions for the signal isolation workflow.
"""

from __future__ import annotations

import numpy as np
from scipy import ndimage
from skimage import morphology


def clip_and_scale_blank(
    blank: np.ndarray,
    clip_factor: int,
    scale_factor: float,
) -> np.ndarray:
    """
    Clip and scale blank channel for subtraction.

    Parameters
    ----------
    blank : np.ndarray
        Blank channel image
    clip_factor : int
        Values below this are zeroed
    scale_factor : float
        Multiply blank by this factor

    Returns
    -------
    np.ndarray
        Processed blank ready for subtraction
    """
    result = blank.astype(np.float32)
    result = np.clip(result, clip_factor, result.max())
    result[result <= clip_factor] = 0
    result *= scale_factor
    return result


def apply_percentile_smoothing(
    image: np.ndarray,
    reference: np.ndarray,
    percentile: int,
    filter_size: int,
    mode: str = "low",
) -> np.ndarray:
    """
    Apply smoothing to pixels below/above a percentile threshold.

    Parameters
    ----------
    image : np.ndarray
        Image to smooth
    reference : np.ndarray
        Reference image for threshold calculation
    percentile : int
        Percentile threshold (0-100)
    filter_size : int
        Size of uniform filter
    mode : str
        'low' to smooth below percentile, 'high' to smooth above

    Returns
    -------
    np.ndarray
        Smoothed image
    """
    threshold = np.percentile(reference.ravel(), percentile)

    if mode == "low":
        mask = reference < threshold
    else:
        mask = reference > threshold

    # Dilate mask for better coverage
    mask = morphology.binary_dilation(mask, morphology.disk(1))

    # Apply smoothing
    smoothed = ndimage.uniform_filter(image.astype(np.float32), size=filter_size)

    # Blend
    result = np.where(mask, smoothed, image)
    return result


def apply_erosion_mask(
    image: np.ndarray,
    erosion_radius: int,
) -> np.ndarray:
    """
    Apply erosion to non-zero regions of image.

    Parameters
    ----------
    image : np.ndarray
        Image to erode
    erosion_radius : int
        Radius of erosion disk

    Returns
    -------
    np.ndarray
        Eroded image
    """
    if erosion_radius <= 0:
        return image

    mask = image > 0
    eroded_mask = morphology.binary_erosion(mask, morphology.disk(erosion_radius))
    return np.where(eroded_mask, image, 0)


def estimate_autofluorescence_level(
    signal: np.ndarray,
    blank: np.ndarray,
    method: str = "correlation",
) -> tuple[float, dict]:
    """
    Estimate the level of autofluorescence in signal channel.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image
    blank : np.ndarray
        Blank channel image
    method : str
        Estimation method: 'correlation', 'ratio', or 'regression'

    Returns
    -------
    tuple
        (af_level, details) where af_level is 0-1 estimate of AF contribution

    Examples
    --------
    >>> af_level, details = estimate_autofluorescence_level(cd3, blank)
    >>> print(f"Autofluorescence: {af_level*100:.1f}%")
    """
    signal_flat = signal.ravel().astype(np.float64)
    blank_flat = blank.ravel().astype(np.float64)

    # Mask out zero regions
    mask = (signal_flat > 0) | (blank_flat > 0)
    signal_masked = signal_flat[mask]
    blank_masked = blank_flat[mask]

    if len(signal_masked) < 100:
        return 0.0, {"method": method, "error": "insufficient_data"}

    details = {"method": method}

    if method == "correlation":
        # Pearson correlation
        mean_s, mean_b = np.mean(signal_masked), np.mean(blank_masked)
        std_s, std_b = np.std(signal_masked), np.std(blank_masked)

        if std_s < 1e-6 or std_b < 1e-6:
            return 0.0, details

        corr = np.mean((signal_masked - mean_s) * (blank_masked - mean_b)) / (std_s * std_b)
        af_level = max(0, corr)
        details["correlation"] = float(corr)

    elif method == "ratio":
        # Ratio of blank contribution
        high_blank_mask = blank_flat[mask] > np.percentile(blank_masked, 75)

        if np.sum(high_blank_mask) < 10:
            return 0.0, details

        signal_in_high_blank = np.mean(signal_masked[high_blank_mask])
        signal_overall = np.mean(signal_masked)

        af_level = (signal_in_high_blank - signal_overall) / signal_in_high_blank
        af_level = max(0, min(1, af_level))
        details["signal_in_high_blank"] = float(signal_in_high_blank)
        details["signal_overall"] = float(signal_overall)

    elif method == "regression":
        # Linear regression slope
        try:
            from scipy.stats import linregress

            slope, intercept, r_value, p_value, std_err = linregress(blank_masked, signal_masked)
            af_level = max(0, min(1, slope * np.mean(blank_masked) / np.mean(signal_masked)))
            details["slope"] = float(slope)
            details["r_squared"] = float(r_value**2)
        except Exception:
            af_level = 0.0
            details["error"] = "regression_failed"

    else:
        raise ValueError(f"Unknown method: {method}")

    return af_level, details


def compute_optimal_clip_factor(
    blank: np.ndarray,
    noise_percentile: int = 10,
    min_clip: int = 0,
    max_clip: int | None = None,
) -> int:
    """
    Compute optimal clip factor for blank channel.

    The clip factor should remove noise while preserving autofluorescence signal.

    Parameters
    ----------
    blank : np.ndarray
        Blank channel image
    noise_percentile : int
        Percentile to consider as noise floor (default: 10)
    min_clip : int
        Minimum clip value
    max_clip : int, optional
        Maximum clip value

    Returns
    -------
    int
        Suggested clip factor
    """
    if max_clip is None:
        max_clip = int(np.percentile(blank, 95))

    # Noise floor estimate
    noise_floor = np.percentile(blank, noise_percentile)

    # Signal floor estimate (where actual AF starts)
    signal_floor = np.percentile(blank, 25)

    # Clip between noise floor and signal floor
    suggested = int(noise_floor + (signal_floor - noise_floor) * 0.5)

    return max(min_clip, min(max_clip, suggested))


def compute_optimal_scale_factor(
    signal: np.ndarray,
    blank: np.ndarray,
    target_residual: float = 0.1,
) -> float:
    """
    Compute optimal scale factor for blank subtraction.

    Finds scale factor that removes autofluorescence while minimizing
    true signal loss.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel
    blank : np.ndarray
        Blank channel
    target_residual : float
        Target residual correlation (0-1), lower = more aggressive

    Returns
    -------
    float
        Suggested scale factor
    """
    # Test a range of scale factors
    test_scales = np.arange(0.5, 2.1, 0.1)
    best_scale = 1.0
    best_score = float("inf")

    signal_float = signal.astype(np.float32)
    blank_float = blank.astype(np.float32)

    for scale in test_scales:
        # Simple subtraction
        subtracted = signal_float - np.minimum(signal_float, blank_float * scale)
        subtracted = np.clip(subtracted, 0, 65535)

        # Compute residual correlation
        mask = (subtracted > 0) | (blank_float > 0)
        if np.sum(mask) < 100:
            continue

        sub_flat = subtracted.ravel()[mask.ravel()]
        blank_flat = blank_float.ravel()[mask.ravel()]

        # Correlation
        corr = np.corrcoef(sub_flat, blank_flat)[0, 1]
        corr = 0 if np.isnan(corr) else corr

        # Score: distance from target residual
        score = abs(corr - target_residual)

        if score < best_score:
            best_score = score
            best_scale = scale

    return round(best_scale, 2)


def adaptive_range_boundaries(
    signal: np.ndarray,
    n_ranges: int = 5,
    method: str = "percentile",
) -> list[tuple[float, float, str]]:
    """
    Compute intensity range boundaries for weighted subtraction.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image
    n_ranges : int
        Number of intensity ranges (default: 5)
    method : str
        Boundary computation method: "percentile" or "otsu"

    Returns
    -------
    list of (lower, upper, label)
        Non-overlapping intensity ranges covering the full signal range.
    """
    flat = signal.ravel().astype(np.float64)
    vmin, vmax = float(np.min(flat)), float(np.max(flat))

    if vmax - vmin < 1:
        # Near-uniform image: single range
        return [(vmin, vmax + 1, "uniform")]

    labels = ["background", "very_dim", "dim", "medium", "bright"]
    # Extend labels if n_ranges differs
    if n_ranges <= 3:
        labels = ["background", "dim", "bright"]
    elif n_ranges == 4:
        labels = ["background", "dim", "medium", "bright"]
    elif n_ranges > 5:
        labels = labels[:5] + [f"range_{i}" for i in range(5, n_ranges)]
    labels = labels[:n_ranges]

    if method == "percentile":
        percentiles = np.linspace(0, 100, n_ranges + 1)
        boundaries = [float(np.percentile(flat, p)) for p in percentiles]
    elif method == "otsu":
        from skimage.filters import threshold_multiotsu

        try:
            thresholds = threshold_multiotsu(signal, classes=min(n_ranges, 5))
            boundaries = [vmin] + [float(t) for t in thresholds] + [vmax + 1]
            # Pad or trim to match n_ranges
            while len(boundaries) < n_ranges + 1:
                # Subdivide largest range
                widths = [boundaries[i + 1] - boundaries[i] for i in range(len(boundaries) - 1)]
                largest = int(np.argmax(widths))
                mid = (boundaries[largest] + boundaries[largest + 1]) / 2
                boundaries.insert(largest + 1, mid)
            boundaries = boundaries[: n_ranges + 1]
        except (ValueError, RuntimeError):
            # Fallback to percentile if Otsu fails
            percentiles = np.linspace(0, 100, n_ranges + 1)
            boundaries = [float(np.percentile(flat, p)) for p in percentiles]
    else:
        raise ValueError(f"Unknown method: {method}. Use 'percentile' or 'otsu'.")

    # Ensure strictly increasing boundaries
    for i in range(1, len(boundaries)):
        if boundaries[i] <= boundaries[i - 1]:
            boundaries[i] = boundaries[i - 1] + 1

    ranges = []
    for i in range(n_ranges):
        ranges.append((boundaries[i], boundaries[i + 1], labels[i]))

    return ranges


def smooth_membership(
    value: np.ndarray,
    lower: float,
    upper: float,
    transition_width: float,
) -> np.ndarray:
    """
    Compute cosine membership for soft range transitions.

    Returns a smooth weight between 0 and 1 indicating how much `value`
    belongs to the range [lower, upper], with cosine transitions at edges.

    Parameters
    ----------
    value : np.ndarray
        Pixel intensity values
    lower : float
        Lower boundary of the range
    upper : float
        Upper boundary of the range (exclusive)
    transition_width : float
        Fraction of the range width used for blending at edges (0-0.5)

    Returns
    -------
    np.ndarray
        Membership weights in [0, 1]
    """
    range_width = upper - lower
    if range_width <= 0:
        return np.zeros_like(value, dtype=np.float32)

    tw = max(0.0, min(0.5, transition_width)) * range_width

    membership = np.ones_like(value, dtype=np.float32)

    if tw > 0:
        # Rising edge at lower boundary
        rising = (value >= lower - tw) & (value < lower + tw)
        if np.any(rising):
            t = (value[rising] - (lower - tw)) / (2 * tw)
            membership[rising] = (1 - np.cos(np.pi * t)) / 2

        # Falling edge at upper boundary
        falling = (value >= upper - tw) & (value < upper + tw)
        if np.any(falling):
            t = (value[falling] - (upper - tw)) / (2 * tw)
            membership[falling] = (1 + np.cos(np.pi * t)) / 2

    # Outside range
    membership[value < lower - tw] = 0
    membership[value >= upper + tw] = 0

    return membership


def validate_subtraction_params(
    blank_clip_factor: int,
    blank_scale_factor: float,
    low_percentile: int,
    high_percentile: int,
    erosion: int,
) -> tuple[bool, str]:
    """
    Validate subtraction parameters are within reasonable ranges.

    Returns
    -------
    tuple
        (is_valid, error_message)
    """
    errors = []

    if blank_clip_factor < 0:
        errors.append("blank_clip_factor must be >= 0")
    if blank_clip_factor > 65535:
        errors.append("blank_clip_factor must be <= 65535")

    if blank_scale_factor < 0:
        errors.append("blank_scale_factor must be >= 0")
    if blank_scale_factor > 10:
        errors.append("blank_scale_factor > 10 is unusually high")

    if not 0 <= low_percentile <= 100:
        errors.append("low_percentile must be 0-100")
    if not 0 <= high_percentile <= 100:
        errors.append("high_percentile must be 0-100")

    if erosion < 0:
        errors.append("erosion must be >= 0")
    if erosion > 10:
        errors.append("erosion > 10 is unusually high")

    if errors:
        return False, "; ".join(errors)
    return True, ""
