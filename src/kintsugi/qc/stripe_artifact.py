"""
Stripe artifact detection for microscopy images.

Detects periodic horizontal stripes caused by microscope vibration or other
mechanical issues during acquisition. Uses fast row-based analysis with
optional FFT confirmation.

Performance: ~1-5ms per image (vs 100ms+ for full 2D FFT approach)
"""

from dataclasses import dataclass, field
from typing import Literal

import numpy as np
from scipy.signal import find_peaks


@dataclass
class StripeArtifactResult:
    """Result of stripe artifact detection for a single image."""

    has_artifact: bool
    severity: Literal["none", "mild", "moderate", "severe"]
    stripe_score: float  # 0-1, higher = worse
    frequency_peaks: list[float] = field(default_factory=list)  # Normalized frequencies
    peak_indices: list[int] = field(default_factory=list)  # Raw peak positions in FFT
    metrics: dict = field(default_factory=dict)

    def __str__(self) -> str:
        if not self.has_artifact:
            return "No stripe artifact detected"
        return (
            f"Stripe artifact: {self.severity} (score={self.stripe_score:.3f}, "
            f"{len(self.frequency_peaks)} peaks)"
        )


@dataclass
class ZStackArtifactReport:
    """Artifact detection results for an entire z-stack."""

    total_planes: int
    affected_planes: list[int] = field(default_factory=list)
    results: dict[int, StripeArtifactResult] = field(default_factory=dict)
    worst_severity: Literal["none", "mild", "moderate", "severe"] = "none"
    summary: str = ""

    @property
    def has_artifacts(self) -> bool:
        return len(self.affected_planes) > 0

    def __str__(self) -> str:
        if not self.has_artifacts:
            return f"No artifacts detected in {self.total_planes} z-planes"
        return (
            f"Artifacts detected in {len(self.affected_planes)}/{self.total_planes} "
            f"z-planes: {self.affected_planes} (worst: {self.worst_severity})"
        )


def detect_stripe_artifacts(
    image: np.ndarray,
    sensitivity: float = 0.5,
    min_frequency: float = 0.01,
    max_frequency: float = 0.5,
) -> StripeArtifactResult:
    """
    Detect horizontal stripe artifacts using fast row-based analysis.

    Uses multiple complementary metrics:
    1. Row difference variance ratio (horizontal vs vertical gradients)
    2. High-frequency power in row means
    3. Alternating row pattern detection

    Parameters
    ----------
    image : np.ndarray
        2D image array (single z-plane)
    sensitivity : float
        Detection sensitivity (0-1). Higher = more sensitive.
        Default 0.5 is balanced for typical microscopy data.
    min_frequency : float
        Minimum normalized frequency (kept for API compatibility)
    max_frequency : float
        Maximum normalized frequency (kept for API compatibility)

    Returns
    -------
    StripeArtifactResult
        Detection result with severity, score, and metrics
    """
    if image.ndim != 2:
        raise ValueError(f"Expected 2D image, got {image.ndim}D")

    img = image.astype(np.float32)
    height, width = img.shape

    # === Metric 1: Directional gradient ratio ===
    # Horizontal stripes have strong vertical gradients (row-to-row changes)
    # compared to horizontal gradients (within-row changes)
    row_means = np.mean(img, axis=1)
    col_means = np.mean(img, axis=0)

    # Variance of row-to-row differences vs column-to-column differences
    v_gradient_var = np.var(np.diff(row_means))
    h_gradient_var = np.var(np.diff(col_means))

    # Ratio > 1 means more vertical structure (horizontal stripes)
    if h_gradient_var > 0:
        gradient_ratio = v_gradient_var / h_gradient_var
    else:
        gradient_ratio = v_gradient_var if v_gradient_var > 0 else 1.0

    # === Metric 2: High-frequency power ratio ===
    # Stripes show up as high-frequency content in row means
    row_centered = row_means - row_means.mean()
    fft = np.fft.rfft(row_centered)
    power = np.abs(fft) ** 2

    # Split into low and high frequency bands
    n_freqs = len(power)
    low_freq_end = max(1, n_freqs // 8)  # Bottom 12.5%
    mid_freq_end = n_freqs // 2  # Bottom 50%

    low_power = power[1:low_freq_end].sum() if low_freq_end > 1 else 0
    mid_power = power[low_freq_end:mid_freq_end].sum()
    high_power = power[mid_freq_end:].sum()
    total_power = power[1:].sum()

    if total_power > 0:
        hf_ratio = high_power / total_power
        mf_ratio = mid_power / total_power
    else:
        hf_ratio = 0.0
        mf_ratio = 0.0

    # === Metric 3: Row difference normalized by signal ===
    # How much does each row differ from its neighbors relative to signal level?
    signal_std = np.std(img)
    row_diff_std = np.std(np.diff(row_means))

    if signal_std > 0:
        normalized_row_diff = row_diff_std / signal_std
    else:
        normalized_row_diff = 0.0

    # === Metric 4: Alternating row pattern (every-other-row stripes) ===
    even_mean = np.mean(row_means[::2])
    odd_mean = np.mean(row_means[1::2])
    alt_diff = abs(even_mean - odd_mean)
    signal_range = np.ptp(row_means)
    if signal_range > 0:
        alt_pattern_score = alt_diff / signal_range
    else:
        alt_pattern_score = 0.0

    # === Combine metrics into stripe score ===
    # Adjust thresholds based on sensitivity
    # sensitivity 0 = strict (fewer detections), 1 = permissive (more detections)
    threshold_scale = 1.0 + (0.5 - sensitivity)  # 0.5 to 1.5

    # Score components (calibrated for microscopy artifacts)
    # gradient_ratio > 1 means vertical gradient dominates (horizontal stripes)
    # Values are dataset-dependent, so we use relative scaling
    #
    # For absolute detection (without baseline comparison):
    # - normalized_row_diff typically 0.005-0.015 for microscopy
    # - Higher values (>0.012) suggest potential artifacts
    # - Use scan_zstack_for_artifacts_fast() for relative comparison

    gradient_score = min(1.0, max(0, gradient_ratio - 0.8) / (2.0 * threshold_scale))
    hf_score = min(1.0, hf_ratio * 200 / threshold_scale)  # Very sensitive to HF
    mf_score = min(1.0, mf_ratio * 20 / threshold_scale)
    # Row diff: 0.010 -> 0.50, 0.015 -> 0.75, 0.020 -> 1.0
    diff_score = min(1.0, (normalized_row_diff - 0.005) * 50 / threshold_scale)
    diff_score = max(0.0, diff_score)
    alt_score = min(1.0, alt_pattern_score * 100 / threshold_scale)

    # Weighted combination - row-diff is primary indicator
    stripe_score = (
        0.15 * gradient_score
        + 0.20 * hf_score
        + 0.10 * mf_score
        + 0.45 * diff_score
        + 0.10 * alt_score
    )

    # Clamp to [0, 1]
    stripe_score = max(0.0, min(1.0, stripe_score))

    # === Find frequency peaks for detailed analysis ===
    frequency_peaks = []
    peak_indices = []

    if stripe_score > 0.05:
        # Only do peak finding if there's potential signal
        power_norm = power / (power.max() + 1e-10)
        peaks, _ = find_peaks(
            power_norm[1:],
            height=0.1 / threshold_scale,
            prominence=0.05 / threshold_scale,
        )
        peak_indices = (peaks + 1).tolist()
        freqs = np.fft.rfftfreq(height)
        frequency_peaks = [float(freqs[i]) for i in peak_indices if i < len(freqs)]

    # === Classify severity ===
    if stripe_score > 0.4:
        severity = "severe"
    elif stripe_score > 0.2:
        severity = "moderate"
    elif stripe_score > 0.08:
        severity = "mild"
    else:
        severity = "none"

    has_artifact = stripe_score > 0.08

    return StripeArtifactResult(
        has_artifact=has_artifact,
        severity=severity,
        stripe_score=stripe_score,
        frequency_peaks=frequency_peaks,
        peak_indices=peak_indices,
        metrics={
            "gradient_ratio": float(gradient_ratio),
            "hf_ratio": float(hf_ratio),
            "mf_ratio": float(mf_ratio),
            "normalized_row_diff": float(normalized_row_diff),
            "alt_pattern_score": float(alt_pattern_score),
            "sensitivity": sensitivity,
        },
    )


def detect_stripe_artifacts_batch(
    images: list[np.ndarray],
    sensitivity: float = 0.5,
) -> list[StripeArtifactResult]:
    """
    Detect stripe artifacts in a batch of images efficiently.

    Parameters
    ----------
    images : list[np.ndarray]
        List of 2D image arrays
    sensitivity : float
        Detection sensitivity (0-1)

    Returns
    -------
    list[StripeArtifactResult]
        Detection results for each image
    """
    return [detect_stripe_artifacts(img, sensitivity=sensitivity) for img in images]


def scan_zstack_for_artifacts(
    z_stack: np.ndarray,
    sensitivity: float = 0.5,
    verbose: bool = True,
) -> ZStackArtifactReport:
    """
    Scan all z-planes in a stack for stripe artifacts.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)
    sensitivity : float
        Detection sensitivity (0-1)
    verbose : bool
        Print progress information

    Returns
    -------
    ZStackArtifactReport
        Comprehensive report of artifacts in the z-stack
    """
    if z_stack.ndim != 3:
        raise ValueError(f"Expected 3D z-stack, got {z_stack.ndim}D")

    nz = z_stack.shape[0]
    results = {}
    affected_planes = []
    worst_severity = "none"
    severity_order = {"none": 0, "mild": 1, "moderate": 2, "severe": 3}

    for z in range(nz):
        result = detect_stripe_artifacts(z_stack[z], sensitivity=sensitivity)
        results[z] = result

        if result.has_artifact:
            affected_planes.append(z)
            if severity_order[result.severity] > severity_order[worst_severity]:
                worst_severity = result.severity

        if verbose and result.has_artifact:
            print(f"  Z-plane {z}: {result}")

    # Generate summary
    if affected_planes:
        summary = (
            f"Detected stripe artifacts in {len(affected_planes)}/{nz} z-planes. "
            f"Affected: {affected_planes}. Worst severity: {worst_severity}."
        )
    else:
        summary = f"No stripe artifacts detected in {nz} z-planes."

    return ZStackArtifactReport(
        total_planes=nz,
        affected_planes=affected_planes,
        results=results,
        worst_severity=worst_severity,
        summary=summary,
    )


def scan_zstack_for_artifacts_fast(
    z_stack: np.ndarray,
    baseline_planes: list[int] | None = None,
    threshold_sigma: float = 2.0,
) -> ZStackArtifactReport:
    """
    Fast z-stack scanning using statistical comparison across planes.

    Compares each plane's stripe metrics against the stack baseline,
    flagging planes that deviate significantly.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)
    baseline_planes : list[int], optional
        Planes to use as baseline (default: use all planes)
    threshold_sigma : float
        Number of standard deviations above baseline to flag as artifact

    Returns
    -------
    ZStackArtifactReport
        Report with affected planes identified by comparison
    """
    if z_stack.ndim != 3:
        raise ValueError(f"Expected 3D z-stack, got {z_stack.ndim}D")

    nz, height, width = z_stack.shape
    stack = z_stack.astype(np.float32)

    # Compute row means for each z-plane
    row_means = np.mean(stack, axis=2)  # Shape: (nz, height)

    # Compute row-diff std for each plane
    row_diffs = np.diff(row_means, axis=1)  # Shape: (nz, height-1)
    row_diff_stds = np.std(row_diffs, axis=1)  # Shape: (nz,)

    # Compute signal std for normalization
    signal_stds = np.std(stack.reshape(nz, -1), axis=1)  # Shape: (nz,)

    # Normalized metric
    with np.errstate(divide="ignore", invalid="ignore"):
        metrics = row_diff_stds / (signal_stds + 1e-10)

    # Compute baseline statistics
    if baseline_planes is None:
        baseline_values = metrics
    else:
        baseline_values = metrics[baseline_planes]

    baseline_mean = np.mean(baseline_values)
    baseline_std = np.std(baseline_values)

    # Flag outliers
    threshold = baseline_mean + threshold_sigma * baseline_std
    affected_planes = np.where(metrics > threshold)[0].tolist()

    # Classify severity for each affected plane
    results = {}
    worst_severity = "none"
    severity_order = {"none": 0, "mild": 1, "moderate": 2, "severe": 3}

    for z in range(nz):
        if z in affected_planes:
            # Convert deviation to severity
            deviation = (metrics[z] - baseline_mean) / (baseline_std + 1e-10)
            if deviation > 4:
                severity = "severe"
            elif deviation > 3:
                severity = "moderate"
            else:
                severity = "mild"

            score = min(1.0, deviation / 5.0)
        else:
            severity = "none"
            score = 0.0

        results[z] = StripeArtifactResult(
            has_artifact=z in affected_planes,
            severity=severity,
            stripe_score=score,
            metrics={"normalized_row_diff": float(metrics[z])},
        )

        if severity_order[severity] > severity_order[worst_severity]:
            worst_severity = severity

    summary = (
        f"Detected {len(affected_planes)}/{nz} affected planes "
        f"(threshold: {threshold:.4f}, worst: {worst_severity})"
    )

    return ZStackArtifactReport(
        total_planes=nz,
        affected_planes=affected_planes,
        results=results,
        worst_severity=worst_severity,
        summary=summary,
    )


def get_fft_visualization(
    image: np.ndarray,
    result: StripeArtifactResult | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Get FFT visualization data for displaying artifact analysis.

    Parameters
    ----------
    image : np.ndarray
        2D image
    result : StripeArtifactResult, optional
        Detection result (will compute if not provided)

    Returns
    -------
    tuple
        (log_power_spectrum, vertical_profile, peak_indices)
    """
    if result is None:
        result = detect_stripe_artifacts(image)

    height, width = image.shape

    # Compute 2D FFT for visualization only
    fft = np.fft.fft2(image.astype(np.float64))
    fft_shifted = np.fft.fftshift(fft)
    power_spectrum = np.abs(fft_shifted) ** 2

    # Log scale for visualization
    log_power = np.log1p(power_spectrum)

    # Vertical profile from row means FFT
    row_means = np.mean(image, axis=1)
    fft_1d = np.fft.rfft(row_means - row_means.mean())
    power_1d = np.abs(fft_1d) ** 2
    if power_1d.max() > 0:
        power_1d = power_1d / power_1d.max()

    return log_power, power_1d, result.peak_indices


def compute_stripe_direction(
    image: np.ndarray,
) -> tuple[str, float]:
    """
    Determine if stripes are primarily horizontal or vertical.

    Uses fast gradient-based analysis instead of full FFT.

    Parameters
    ----------
    image : np.ndarray
        2D image

    Returns
    -------
    tuple
        (direction, ratio)
        - direction: "horizontal", "vertical", or "none"
        - ratio: Strength of directional bias (>1 = horizontal, <1 = vertical)
    """
    img = image.astype(np.float32)

    # Compare row vs column gradient variance
    row_means = np.mean(img, axis=1)
    col_means = np.mean(img, axis=0)

    v_var = np.var(np.diff(row_means))
    h_var = np.var(np.diff(col_means))

    if h_var == 0:
        ratio = float("inf") if v_var > 0 else 1.0
    else:
        ratio = v_var / h_var

    if ratio > 2.0:
        direction = "horizontal"
    elif ratio < 0.5:
        direction = "vertical"
    else:
        direction = "none"

    return direction, float(ratio)
