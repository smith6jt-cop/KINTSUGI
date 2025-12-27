"""
Stripe artifact detection for microscopy images.

This module provides the core detection algorithm used by ArtifactScanner.
For project-level scanning, use the ArtifactScanner class instead.

The detection uses row-based statistical analysis to identify horizontal
stripes caused by microscope vibration or other acquisition issues.
"""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np


@dataclass
class StripeArtifactResult:
    """Result of stripe artifact detection for a single image."""

    has_artifact: bool
    severity: Literal["none", "mild", "moderate", "severe"]
    score: float  # 0-1, higher = worse
    metrics: dict = field(default_factory=dict)
    frequency_peaks: list = field(default_factory=list)  # For backwards compatibility

    @property
    def stripe_score(self) -> float:
        """Backwards compatibility alias for score."""
        return self.score

    def __str__(self) -> str:
        if not self.has_artifact:
            return "No stripe artifact detected"
        return f"Stripe artifact: {self.severity} (score={self.score:.3f})"


def _compute_highfreq_stripe_metric(img: np.ndarray) -> float:
    """
    Compute high-frequency stripe metric for deconvolution ringing artifacts.

    Downsamples image to 1/4 dimensions for speed, then analyzes patches
    in low-intensity areas where artifacts are most visible.

    Returns a normalized score where higher values indicate stronger
    high-frequency horizontal banding.
    """
    from scipy.fft import fft2, fftshift

    # Downsample to 1/4 dimensions (1/16 pixels) for speed
    img_ds = img[::4, ::4].astype(np.float32)

    # Get non-zero pixels to determine intensity thresholds
    nonzero = img_ds[img_ds > 0]
    if len(nonzero) < 1000:
        return 0.0

    # Focus on lower 25% of non-zero values where artifacts are most visible
    low_thresh = np.percentile(nonzero, 25)
    low_mask = (img_ds > 0) & (img_ds < low_thresh)

    # Scan patches on downsampled image
    patch_size = 128  # Smaller patches for downsampled image
    step_size = patch_size
    best_metric = 0.0

    y_steps = range(0, img_ds.shape[0] - patch_size, step_size)
    x_steps = range(0, img_ds.shape[1] - patch_size, step_size)

    for y in y_steps:
        for x in x_steps:
            mask_patch = low_mask[y:y + patch_size, x:x + patch_size]

            # Only analyze if >50% of patch is low-intensity
            if np.mean(mask_patch) < 0.5:
                continue

            patch = img_ds[y:y + patch_size, x:x + patch_size]

            # 2D FFT of this patch
            fft_result = fftshift(fft2(patch))
            power = np.abs(fft_result) ** 2

            ph, pw = power.shape
            cy, cx = ph // 2, pw // 2

            # Horizontal stripe energy (vertical axis in FFT, excludes DC)
            fx_band = 3
            fy_high_start = ph // 4
            h_stripe = power[cy + fy_high_start:, cx - fx_band:cx + fx_band].sum()
            h_stripe += power[:cy - fy_high_start, cx - fx_band:cx + fx_band].sum()

            # Vertical stripe energy (horizontal axis in FFT)
            fy_band = 3
            fx_high_start = pw // 4
            v_stripe = power[cy - fy_band:cy + fy_band, cx + fx_high_start:].sum()
            v_stripe += power[cy - fy_band:cy + fy_band, :cx - fx_high_start].sum()

            # Ratio of horizontal to vertical stripe energy
            if v_stripe > 0:
                ratio = h_stripe / v_stripe
            else:
                ratio = 1.0

            best_metric = max(best_metric, ratio)

    # Normalize: ratio of 1.0 = no bias, higher = horizontal stripes
    return float(max(0.0, best_metric - 1.0))


def detect_stripe_artifact(
    image: np.ndarray,
    baseline_row_diff: float | None = None,
    baseline_highfreq: float | None = None,
    threshold_sigma: float = 2.0,
) -> StripeArtifactResult:
    """
    Detect horizontal stripe artifacts in a single image.

    Uses both low-frequency (row-mean based) and high-frequency
    (pixel-level alternating rows) metrics to detect stripe artifacts
    from various sources including vibration and deconvolution ringing.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    baseline_row_diff : float, optional
        Expected baseline row-diff metric. If provided, uses comparative
        detection. If None, uses absolute thresholds.
    baseline_highfreq : float, optional
        Expected baseline high-frequency metric for comparative detection.
    threshold_sigma : float
        For comparative detection, number of std devs above baseline.

    Returns
    -------
    StripeArtifactResult
        Detection result with severity and metrics
    """
    if image.ndim != 2:
        raise ValueError(f"Expected 2D image, got {image.ndim}D")

    img = image.astype(np.float32)

    # Compute normalized row-difference metric (low-frequency stripes)
    row_means = np.mean(img, axis=1)
    row_diff_std = np.std(np.diff(row_means))
    signal_std = np.std(img)

    if signal_std > 0:
        normalized_row_diff = row_diff_std / signal_std
    else:
        normalized_row_diff = 0.0

    # Compute high-frequency stripe metric (deconvolution ringing)
    highfreq_metric = _compute_highfreq_stripe_metric(img)

    # Compute gradient ratio (row vs column variance)
    col_means = np.mean(img, axis=0)
    v_gradient_var = np.var(np.diff(row_means))
    h_gradient_var = np.var(np.diff(col_means))

    if h_gradient_var > 0:
        gradient_ratio = v_gradient_var / h_gradient_var
    else:
        gradient_ratio = v_gradient_var if v_gradient_var > 0 else 1.0

    metrics = {
        "normalized_row_diff": float(normalized_row_diff),
        "highfreq_stripe": float(highfreq_metric),
        "gradient_ratio": float(gradient_ratio),
        "row_diff_std": float(row_diff_std),
        "signal_std": float(signal_std),
    }

    # Determine if artifact is present using both metrics
    lowfreq_artifact = False
    highfreq_artifact = False
    lowfreq_score = 0.0
    highfreq_score = 0.0

    if baseline_row_diff is not None:
        # Comparative detection for low-frequency
        deviation = normalized_row_diff - baseline_row_diff
        lowfreq_score = min(1.0, max(0.0, deviation * 50))
        lowfreq_artifact = normalized_row_diff > baseline_row_diff * (1 + threshold_sigma * 0.1)
    else:
        # Absolute detection for low-frequency
        lowfreq_score = min(1.0, max(0.0, (normalized_row_diff - 0.005) * 50))
        lowfreq_artifact = normalized_row_diff > 0.012

    # High-frequency detection (deconvolution ringing)
    # Use absolute thresholds calibrated on typical data:
    # - Clean images: highfreq < 10
    # - Mild artifacts: 10-20
    # - Moderate artifacts: 20-30
    # - Severe artifacts: > 30
    # Comparative detection can contaminate baseline when multiple planes are affected
    highfreq_artifact = highfreq_metric > 15.0
    if highfreq_metric > 30:
        highfreq_score = min(1.0, 0.7 + (highfreq_metric - 30) / 50.0)  # Severe
    elif highfreq_metric > 20:
        highfreq_score = 0.4 + (highfreq_metric - 20) / 33.0  # Moderate
    elif highfreq_metric > 10:
        highfreq_score = 0.1 + (highfreq_metric - 10) / 33.0  # Mild
    else:
        highfreq_score = 0.0

    # Combine scores - take the maximum as either type indicates artifact
    has_artifact = lowfreq_artifact or highfreq_artifact
    score = max(lowfreq_score, highfreq_score)

    # Classify severity
    if score > 0.6:
        severity = "severe"
    elif score > 0.3:
        severity = "moderate"
    elif has_artifact:
        severity = "mild"
    else:
        severity = "none"

    return StripeArtifactResult(
        has_artifact=has_artifact,
        severity=severity,
        score=score,
        metrics=metrics,
    )


def compute_zstack_baseline(z_stack: np.ndarray) -> dict:
    """
    Compute baseline statistics for a z-stack.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)

    Returns
    -------
    dict
        Dictionary with baseline statistics:
        - 'lowfreq_mean', 'lowfreq_std': normalized_row_diff stats
        - 'highfreq_mean', 'highfreq_std': high-frequency stripe stats
    """
    if z_stack.ndim != 3:
        raise ValueError(f"Expected 3D z-stack, got {z_stack.ndim}D")

    nz = z_stack.shape[0]
    stack = z_stack.astype(np.float32)

    # Compute low-frequency metric for each plane
    row_means = np.mean(stack, axis=2)
    row_diffs = np.diff(row_means, axis=1)
    row_diff_stds = np.std(row_diffs, axis=1)
    signal_stds = np.std(stack.reshape(nz, -1), axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        lowfreq_metrics = row_diff_stds / (signal_stds + 1e-10)

    # Compute high-frequency metric for each plane
    highfreq_metrics = np.array([
        _compute_highfreq_stripe_metric(stack[z]) for z in range(nz)
    ])

    return {
        'lowfreq_mean': float(np.mean(lowfreq_metrics)),
        'lowfreq_std': float(np.std(lowfreq_metrics)),
        'highfreq_mean': float(np.mean(highfreq_metrics)),
        'highfreq_std': float(np.std(highfreq_metrics)),
    }


def scan_zstack(
    z_stack: np.ndarray,
    threshold_sigma: float = 2.0,
) -> tuple[list[int], dict[int, StripeArtifactResult]]:
    """
    Scan a z-stack for stripe artifacts using comparative analysis.

    This is the recommended method for artifact detection as it
    automatically determines thresholds from the data. Detects both
    low-frequency (vibration) and high-frequency (deconvolution) artifacts.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)
    threshold_sigma : float
        Number of standard deviations above baseline to flag

    Returns
    -------
    tuple
        (affected_planes, results_dict)
        - affected_planes: List of z-indices with artifacts
        - results_dict: Dict mapping z-index to StripeArtifactResult
    """
    if z_stack.ndim != 3:
        raise ValueError(f"Expected 3D z-stack, got {z_stack.ndim}D")

    nz = z_stack.shape[0]

    # Compute baseline for low frequency metric (comparative detection)
    baseline = compute_zstack_baseline(z_stack)
    # Use both comparative threshold AND minimum absolute threshold
    # to avoid false positives when baseline variance is very tight
    lowfreq_comparative = baseline['lowfreq_mean'] + threshold_sigma * baseline['lowfreq_std']
    lowfreq_absolute_min = 0.012  # Minimum to flag as artifact
    lowfreq_threshold = max(lowfreq_comparative, lowfreq_absolute_min)

    # High-frequency uses absolute thresholds (comparative doesn't work well
    # when multiple planes have artifacts, which contaminates the baseline)
    highfreq_threshold = 15.0

    # Scan each plane
    results = {}
    affected = []

    for z in range(nz):
        result = detect_stripe_artifact(
            z_stack[z],
            baseline_row_diff=baseline['lowfreq_mean'],
            baseline_highfreq=None,  # Use absolute thresholds
            threshold_sigma=threshold_sigma,
        )

        # Override low-frequency detection using comparative threshold
        lowfreq_metric = result.metrics["normalized_row_diff"]
        highfreq_metric = result.metrics["highfreq_stripe"]

        lowfreq_exceeds = lowfreq_metric > lowfreq_threshold
        highfreq_exceeds = highfreq_metric > highfreq_threshold

        if lowfreq_exceeds or highfreq_exceeds:
            # Determine severity based on metrics
            # Low-frequency: use comparative deviation
            lowfreq_dev = (lowfreq_metric - baseline['lowfreq_mean']) / (baseline['lowfreq_std'] + 1e-10)

            # High-frequency: use absolute thresholds
            if highfreq_metric > 30:
                highfreq_severity_score = 0.7 + (highfreq_metric - 30) / 50.0
            elif highfreq_metric > 20:
                highfreq_severity_score = 0.4 + (highfreq_metric - 20) / 33.0
            elif highfreq_metric > 10:
                highfreq_severity_score = 0.1 + (highfreq_metric - 10) / 33.0
            else:
                highfreq_severity_score = 0.0

            lowfreq_severity_score = min(1.0, lowfreq_dev / 5.0)
            max_score = max(lowfreq_severity_score, highfreq_severity_score)

            result.has_artifact = True
            result.score = min(1.0, max_score)

            if result.score > 0.6:
                result.severity = "severe"
            elif result.score > 0.3:
                result.severity = "moderate"
            else:
                result.severity = "mild"

            affected.append(z)
        else:
            result.has_artifact = False
            result.severity = "none"
            result.score = 0.0

        results[z] = result

    return affected, results
