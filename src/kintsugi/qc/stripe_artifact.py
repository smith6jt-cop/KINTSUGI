"""
Stripe artifact detection for microscopy images.

Detects periodic horizontal stripes caused by microscope vibration or other
mechanical issues during acquisition. Uses FFT analysis to identify anomalous
frequency peaks in the vertical direction.
"""

from dataclasses import dataclass, field
from typing import Literal, Optional
import numpy as np
from scipy.signal import find_peaks
from scipy.ndimage import gaussian_filter1d


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
    Detect periodic horizontal stripe artifacts using FFT analysis.

    Horizontal stripes appear as peaks in the vertical frequency spectrum
    (perpendicular to the stripe direction).

    Parameters
    ----------
    image : np.ndarray
        2D image array (single z-plane)
    sensitivity : float
        Detection sensitivity (0-1). Higher = more sensitive, more false positives.
        Default 0.5 is balanced for typical microscopy data.
    min_frequency : float
        Minimum normalized frequency to consider (excludes low-frequency variations)
    max_frequency : float
        Maximum normalized frequency to consider (excludes high-frequency noise)

    Returns
    -------
    StripeArtifactResult
        Detection result with severity, score, and frequency information
    """
    if image.ndim != 2:
        raise ValueError(f"Expected 2D image, got {image.ndim}D")

    height, width = image.shape

    # 1. Compute 2D FFT and shift DC to center
    fft = np.fft.fft2(image.astype(np.float64))
    fft_shifted = np.fft.fftshift(fft)
    power_spectrum = np.abs(fft_shifted) ** 2

    # 2. Extract vertical frequency profile
    # Horizontal stripes → energy along vertical axis in FFT
    # Sum along horizontal axis to get 1D vertical profile
    vertical_profile = np.sum(power_spectrum, axis=1)

    center_y = height // 2

    # 3. Remove DC component and low frequencies
    min_idx = int(min_frequency * height)
    max_idx = int(max_frequency * height)

    # Work with positive frequencies only (profile is symmetric)
    positive_profile = vertical_profile[center_y:].copy()
    positive_profile[0] = 0  # Remove DC

    # Apply smoothing to reduce noise
    smoothed = gaussian_filter1d(positive_profile, sigma=2)

    # Normalize
    if smoothed.max() > 0:
        profile_norm = smoothed / smoothed.max()
    else:
        profile_norm = smoothed

    # 4. Find peaks in the frequency profile
    # Adjust thresholds based on sensitivity
    height_threshold = (1.0 - sensitivity) * 0.05 + 0.01  # 0.01 to 0.06
    prominence_threshold = (1.0 - sensitivity) * 0.03 + 0.01  # 0.01 to 0.04

    peaks, properties = find_peaks(
        profile_norm[min_idx:max_idx],
        height=height_threshold,
        prominence=prominence_threshold,
        distance=3,  # Minimum 3 bins between peaks
    )

    # Adjust peak indices to account for offset
    peaks = peaks + min_idx

    # 5. Calculate stripe score from peak prominences
    if len(peaks) > 0:
        prominences = properties.get('prominences', np.array([]))
        if len(prominences) > 0:
            # Weight score by both prominence and number of peaks
            max_prominence = np.max(prominences)
            mean_prominence = np.mean(prominences)
            peak_count_factor = min(1.0, len(peaks) / 5)  # Saturate at 5 peaks

            stripe_score = 0.5 * max_prominence + 0.3 * mean_prominence + 0.2 * peak_count_factor
            stripe_score = min(1.0, stripe_score * 2)  # Scale up for visibility
        else:
            stripe_score = 0.0
    else:
        stripe_score = 0.0

    # 6. Convert peak indices to normalized frequencies
    frequencies = [p / height for p in peaks]

    # 7. Classify severity based on score
    if stripe_score > 0.4:
        severity = "severe"
    elif stripe_score > 0.2:
        severity = "moderate"
    elif stripe_score > 0.08:
        severity = "mild"
    else:
        severity = "none"

    has_artifact = stripe_score > 0.08

    # Safely extract prominence metrics (handle empty arrays)
    prominences = properties.get('prominences', np.array([]))
    if len(prominences) > 0:
        max_prominence = float(np.max(prominences))
        mean_prominence = float(np.mean(prominences))
    else:
        max_prominence = 0.0
        mean_prominence = 0.0

    return StripeArtifactResult(
        has_artifact=has_artifact,
        severity=severity,
        stripe_score=stripe_score,
        frequency_peaks=frequencies,
        peak_indices=peaks.tolist(),
        metrics={
            'peak_count': len(peaks),
            'max_prominence': max_prominence,
            'mean_prominence': mean_prominence,
            'profile_energy': float(np.sum(profile_norm[min_idx:max_idx])),
            'sensitivity': sensitivity,
        }
    )


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


def get_fft_visualization(
    image: np.ndarray,
    result: Optional[StripeArtifactResult] = None,
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
        - log_power_spectrum: Log-scaled FFT magnitude for display
        - vertical_profile: 1D profile along vertical axis
        - peak_indices: Indices of detected peaks in profile
    """
    if result is None:
        result = detect_stripe_artifacts(image)

    height, width = image.shape

    # Compute FFT
    fft = np.fft.fft2(image.astype(np.float64))
    fft_shifted = np.fft.fftshift(fft)
    power_spectrum = np.abs(fft_shifted) ** 2

    # Log scale for visualization
    log_power = np.log1p(power_spectrum)

    # Vertical profile
    vertical_profile = np.sum(power_spectrum, axis=1)
    center_y = height // 2
    positive_profile = vertical_profile[center_y:]
    positive_profile[0] = 0  # Remove DC

    # Normalize for display
    if positive_profile.max() > 0:
        positive_profile = positive_profile / positive_profile.max()

    return log_power, positive_profile, result.peak_indices


def compute_stripe_direction(
    image: np.ndarray,
) -> tuple[str, float]:
    """
    Determine if stripes are primarily horizontal or vertical.

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
    fft = np.fft.fft2(image.astype(np.float64))
    fft_shifted = np.fft.fftshift(fft)
    power = np.abs(fft_shifted) ** 2

    height, width = image.shape
    center_y, center_x = height // 2, width // 2

    # Exclude DC and very low frequencies
    exclude_radius = 5

    # Sum energy along vertical axis (horizontal stripes)
    vertical_energy = 0
    for y in range(exclude_radius, height // 4):
        vertical_energy += power[center_y - y, center_x] + power[center_y + y, center_x]

    # Sum energy along horizontal axis (vertical stripes)
    horizontal_energy = 0
    for x in range(exclude_radius, width // 4):
        horizontal_energy += power[center_y, center_x - x] + power[center_y, center_x + x]

    if horizontal_energy == 0:
        ratio = float('inf') if vertical_energy > 0 else 1.0
    else:
        ratio = vertical_energy / horizontal_energy

    if ratio > 2.0:
        direction = "horizontal"
    elif ratio < 0.5:
        direction = "vertical"
    else:
        direction = "none"

    return direction, ratio
