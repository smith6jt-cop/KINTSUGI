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

    def __str__(self) -> str:
        if not self.has_artifact:
            return "No stripe artifact detected"
        return f"Stripe artifact: {self.severity} (score={self.score:.3f})"


def detect_stripe_artifact(
    image: np.ndarray,
    baseline_row_diff: float | None = None,
    threshold_sigma: float = 2.0,
) -> StripeArtifactResult:
    """
    Detect horizontal stripe artifacts in a single image.

    Uses normalized row-difference metric to detect anomalous
    horizontal banding patterns.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    baseline_row_diff : float, optional
        Expected baseline row-diff metric. If provided, uses comparative
        detection. If None, uses absolute thresholds.
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

    # Compute normalized row-difference metric
    row_means = np.mean(img, axis=1)
    row_diff_std = np.std(np.diff(row_means))
    signal_std = np.std(img)

    if signal_std > 0:
        normalized_row_diff = row_diff_std / signal_std
    else:
        normalized_row_diff = 0.0

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
        "gradient_ratio": float(gradient_ratio),
        "row_diff_std": float(row_diff_std),
        "signal_std": float(signal_std),
    }

    # Determine if artifact is present
    if baseline_row_diff is not None:
        # Comparative detection (preferred)
        # Score based on deviation from baseline
        deviation = normalized_row_diff - baseline_row_diff
        score = min(1.0, max(0.0, deviation * 50))  # Scale to 0-1

        # Threshold based on provided baseline
        has_artifact = normalized_row_diff > baseline_row_diff * (1 + threshold_sigma * 0.1)
    else:
        # Absolute detection (fallback)
        # Calibrated for typical microscopy data
        # normalized_row_diff typically 0.005-0.015
        score = min(1.0, max(0.0, (normalized_row_diff - 0.005) * 50))
        has_artifact = normalized_row_diff > 0.012  # Empirical threshold

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


def compute_zstack_baseline(z_stack: np.ndarray) -> tuple[float, float]:
    """
    Compute baseline statistics for a z-stack.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)

    Returns
    -------
    tuple[float, float]
        (mean, std) of normalized_row_diff across z-planes
    """
    if z_stack.ndim != 3:
        raise ValueError(f"Expected 3D z-stack, got {z_stack.ndim}D")

    nz = z_stack.shape[0]
    stack = z_stack.astype(np.float32)

    # Compute metric for each plane
    row_means = np.mean(stack, axis=2)
    row_diffs = np.diff(row_means, axis=1)
    row_diff_stds = np.std(row_diffs, axis=1)
    signal_stds = np.std(stack.reshape(nz, -1), axis=1)

    with np.errstate(divide="ignore", invalid="ignore"):
        metrics = row_diff_stds / (signal_stds + 1e-10)

    return float(np.mean(metrics)), float(np.std(metrics))


def scan_zstack(
    z_stack: np.ndarray,
    threshold_sigma: float = 2.0,
) -> tuple[list[int], dict[int, StripeArtifactResult]]:
    """
    Scan a z-stack for stripe artifacts using comparative analysis.

    This is the recommended method for artifact detection as it
    automatically determines thresholds from the data.

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

    # Compute baseline
    baseline_mean, baseline_std = compute_zstack_baseline(z_stack)
    threshold = baseline_mean + threshold_sigma * baseline_std

    # Scan each plane
    results = {}
    affected = []

    for z in range(nz):
        result = detect_stripe_artifact(
            z_stack[z],
            baseline_row_diff=baseline_mean,
            threshold_sigma=threshold_sigma,
        )

        # Override detection using comparative threshold
        metric = result.metrics["normalized_row_diff"]
        if metric > threshold:
            deviation = (metric - baseline_mean) / (baseline_std + 1e-10)
            result.has_artifact = True
            result.score = min(1.0, deviation / 5.0)

            if deviation > 4:
                result.severity = "severe"
            elif deviation > 3:
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
