"""
Marker-Level Quality Control

Quality assessment for individual markers in multiplex imaging.
Inspired by CyLinter's marker validation approaches.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

import numpy as np
from scipy import stats

logger = logging.getLogger("kintsugi.qc.marker_qc")


@dataclass
class MarkerQCResult:
    """Results from marker-level quality control."""

    marker_name: str
    passed: bool
    quality_score: float

    # Metrics
    metrics: dict[str, float] = field(default_factory=dict)

    # Issues
    issues: list[str] = field(default_factory=list)

    # Expected vs observed
    expected_positive_rate: float | None = None
    observed_positive_rate: float | None = None


# Known marker characteristics for reference
MARKER_REFERENCE = {
    # Nuclear markers
    "dapi": {"location": "nuclear", "expected_positive": 0.95, "ubiquitous": True},
    "hoechst": {"location": "nuclear", "expected_positive": 0.95, "ubiquitous": True},
    # T-cell markers
    "cd3": {"location": "membrane", "expected_positive": (0.05, 0.40), "ubiquitous": False},
    "cd4": {"location": "membrane", "expected_positive": (0.02, 0.25), "ubiquitous": False},
    "cd8": {"location": "membrane", "expected_positive": (0.02, 0.20), "ubiquitous": False},
    # B-cell markers
    "cd20": {"location": "membrane", "expected_positive": (0.02, 0.30), "ubiquitous": False},
    "cd19": {"location": "membrane", "expected_positive": (0.02, 0.30), "ubiquitous": False},
    # Myeloid markers
    "cd68": {"location": "cytoplasm", "expected_positive": (0.01, 0.15), "ubiquitous": False},
    "cd163": {"location": "membrane", "expected_positive": (0.01, 0.10), "ubiquitous": False},
    "cd11c": {"location": "membrane", "expected_positive": (0.01, 0.10), "ubiquitous": False},
    # Proliferation
    "ki67": {"location": "nuclear", "expected_positive": (0.01, 0.30), "ubiquitous": False},
    # Structural
    "panck": {"location": "cytoplasm", "expected_positive": (0.1, 0.8), "ubiquitous": False},
    "cytokeratin": {"location": "cytoplasm", "expected_positive": (0.1, 0.8), "ubiquitous": False},
    "sma": {"location": "cytoplasm", "expected_positive": (0.01, 0.20), "ubiquitous": False},
    "vimentin": {"location": "cytoplasm", "expected_positive": (0.05, 0.50), "ubiquitous": False},
}


class MarkerQC:
    """
    Quality control for individual markers in multiplex imaging.

    Validates:
    - Expression patterns
    - Background levels
    - Signal-to-noise
    - Expected vs observed positive rates
    - Crosstalk with other markers

    Example
    -------
    >>> qc = MarkerQC()
    >>> result = qc.assess(
    ...     intensities=cd3_intensities,
    ...     marker_name="CD3",
    ...     tissue_type="tonsil",
    ... )
    >>> if not result.passed:
    ...     print(result.issues)
    """

    def __init__(
        self,
        min_snr: float = 3.0,
        max_background_pct: float = 0.3,
        positive_rate_tolerance: float = 0.5,
    ):
        """
        Initialize Marker QC.

        Parameters
        ----------
        min_snr : float
            Minimum signal-to-noise ratio
        max_background_pct : float
            Maximum fraction of dynamic range as background
        positive_rate_tolerance : float
            Tolerance for expected vs observed positive rate
        """
        self.min_snr = min_snr
        self.max_background_pct = max_background_pct
        self.positive_rate_tolerance = positive_rate_tolerance

    def assess(
        self,
        intensities: np.ndarray,
        marker_name: str,
        tissue_type: str | None = None,
        background_intensities: np.ndarray | None = None,
    ) -> MarkerQCResult:
        """
        Assess marker quality.

        Parameters
        ----------
        intensities : np.ndarray
            Per-cell or per-pixel intensities
        marker_name : str
            Marker name
        tissue_type : str, optional
            Tissue type for expected expression patterns
        background_intensities : np.ndarray, optional
            Background/control channel intensities for comparison

        Returns
        -------
        MarkerQCResult
            Quality assessment results
        """
        intensities = np.asarray(intensities, dtype=float)
        issues = []
        metrics = {}

        # Normalize marker name
        marker_key = marker_name.lower().replace("-", "").replace("_", "")

        # Get reference if available
        reference = MARKER_REFERENCE.get(marker_key, {})

        # 1. Basic statistics
        p5, p25, p50, p75, p95 = np.percentile(intensities, [5, 25, 50, 75, 95])
        metrics["median"] = float(p50)
        metrics["iqr"] = float(p75 - p25)
        metrics["p5"] = float(p5)
        metrics["p95"] = float(p95)

        # 2. SNR estimation
        background_region = intensities[intensities < p25]
        signal_region = intensities[intensities > p75]

        if len(background_region) > 10:
            noise_std = np.std(background_region)
            background_mean = np.mean(background_region)
        else:
            noise_std = np.std(intensities) * 0.1
            background_mean = p5

        signal_mean = np.mean(signal_region) if len(signal_region) > 10 else p95

        snr = (signal_mean - background_mean) / (noise_std + 1e-8)
        metrics["snr"] = float(snr)
        metrics["background_level"] = float(background_mean)
        metrics["signal_level"] = float(signal_mean)

        if snr < self.min_snr:
            issues.append(f"Low SNR ({snr:.1f} < {self.min_snr})")

        # 3. Background ratio
        dynamic_range = p95 - p5
        if dynamic_range > 0:
            background_ratio = (background_mean - p5) / dynamic_range
            metrics["background_ratio"] = float(background_ratio)

            if background_ratio > self.max_background_pct:
                issues.append(f"High background ({background_ratio * 100:.0f}%)")

        # 4. Distribution quality
        distribution_quality = self._assess_distribution(intensities)
        metrics.update(distribution_quality)

        if distribution_quality.get("bimodality_score", 0) < 0.1 and not reference.get(
            "ubiquitous", False
        ):
            issues.append("Poor separation of positive/negative populations")

        # 5. Expected positive rate
        if reference and "expected_positive" in reference:
            expected = reference["expected_positive"]
            observed = self._estimate_positive_rate(intensities)
            metrics["observed_positive_rate"] = float(observed)

            if isinstance(expected, tuple):
                expected_low, expected_high = expected
                if observed < expected_low * (1 - self.positive_rate_tolerance):
                    issues.append(
                        f"Lower than expected positive rate ({observed * 100:.1f}% < {expected_low * 100:.0f}%)"
                    )
                elif observed > expected_high * (1 + self.positive_rate_tolerance):
                    issues.append(
                        f"Higher than expected positive rate ({observed * 100:.1f}% > {expected_high * 100:.0f}%)"
                    )
            else:
                metrics["expected_positive_rate"] = expected

        # 6. Check for artifacts
        artifact_score = self._check_for_artifacts(intensities)
        metrics["artifact_score"] = float(artifact_score)

        if artifact_score > 0.5:
            issues.append("Potential artifacts detected in intensity distribution")

        # Overall score
        quality_score = self._compute_quality_score(metrics, issues)
        metrics["quality_score"] = quality_score

        return MarkerQCResult(
            marker_name=marker_name,
            passed=len(issues) == 0 and quality_score > 0.5,
            quality_score=quality_score,
            metrics=metrics,
            issues=issues,
            expected_positive_rate=metrics.get("expected_positive_rate"),
            observed_positive_rate=metrics.get("observed_positive_rate"),
        )

    def _assess_distribution(self, intensities: np.ndarray) -> dict[str, float]:
        """Assess quality of intensity distribution."""
        metrics = {}

        # Coefficient of variation
        cv = np.std(intensities) / (np.mean(intensities) + 1e-8)
        metrics["coefficient_of_variation"] = float(cv)

        # Bimodality using Hartigan's dip statistic approximation
        # Use histogram-based approach for efficiency
        hist, edges = np.histogram(intensities, bins=50)
        hist / (hist.sum() + 1e-8)

        # Find local minima as potential separation points
        peaks = []
        for i in range(1, len(hist) - 1):
            if hist[i] > hist[i - 1] and hist[i] > hist[i + 1]:
                peaks.append(i)

        # Bimodality score based on peak separation
        if len(peaks) >= 2:
            # Check if there's a clear valley between peaks
            valley_min = min(hist[peaks[0] : peaks[1] + 1])
            peak_avg = (hist[peaks[0]] + hist[peaks[1]]) / 2
            bimodality = 1 - (valley_min / (peak_avg + 1e-8))
        else:
            bimodality = 0.0

        metrics["bimodality_score"] = float(np.clip(bimodality, 0, 1))
        metrics["n_peaks_detected"] = len(peaks)

        # Skewness
        skewness = stats.skew(intensities)
        metrics["skewness"] = float(skewness)

        return metrics

    def _estimate_positive_rate(self, intensities: np.ndarray) -> float:
        """Estimate fraction of positive cells using Otsu's method."""
        # Simplified Otsu threshold
        hist, edges = np.histogram(intensities, bins=256)
        hist = hist.astype(float)

        total = hist.sum()
        current_max = 0
        threshold_idx = 0

        sum_total = np.sum(np.arange(256) * hist)
        sum_bg = 0
        weight_bg = 0

        for i in range(256):
            weight_bg += hist[i]
            if weight_bg == 0:
                continue

            weight_fg = total - weight_bg
            if weight_fg == 0:
                break

            sum_bg += i * hist[i]
            mean_bg = sum_bg / weight_bg
            mean_fg = (sum_total - sum_bg) / weight_fg

            variance_between = weight_bg * weight_fg * (mean_bg - mean_fg) ** 2

            if variance_between > current_max:
                current_max = variance_between
                threshold_idx = i

        threshold = edges[threshold_idx]
        positive_rate = np.sum(intensities > threshold) / len(intensities)

        return float(positive_rate)

    def _check_for_artifacts(self, intensities: np.ndarray) -> float:
        """Check for intensity distribution artifacts."""
        # Look for suspicious patterns
        artifact_score = 0.0

        # 1. Check for spike at zero
        zero_fraction = np.sum(intensities == 0) / len(intensities)
        if zero_fraction > 0.1:
            artifact_score += 0.3

        # 2. Check for saturation spike
        p99 = np.percentile(intensities, 99)
        if np.sum(intensities >= p99) / len(intensities) > 0.05:
            artifact_score += 0.3

        # 3. Check for unusual histogram patterns (multiple sharp spikes)
        hist, _ = np.histogram(intensities, bins=256)
        hist_smooth = np.convolve(hist, np.ones(5) / 5, mode="same")
        diff = np.abs(hist - hist_smooth)
        spike_fraction = np.sum(diff > hist.max() * 0.3) / len(hist)

        if spike_fraction > 0.1:
            artifact_score += 0.4

        return float(np.clip(artifact_score, 0, 1))

    def _compute_quality_score(
        self,
        metrics: dict[str, float],
        issues: list[str],
    ) -> float:
        """Compute overall marker quality score."""
        score = 1.0

        # SNR contribution
        snr = metrics.get("snr", 0)
        snr_score = min(1.0, snr / 10.0)
        score *= 0.7 + 0.3 * snr_score

        # Bimodality contribution (for non-ubiquitous markers)
        bimodality = metrics.get("bimodality_score", 0)
        score *= 0.7 + 0.3 * bimodality

        # Artifact penalty
        artifact_score = metrics.get("artifact_score", 0)
        score *= 1 - artifact_score * 0.5

        # Issue penalties
        score *= 1 - len(issues) * 0.15

        return float(np.clip(score, 0, 1))


def validate_marker_expression(
    intensities: np.ndarray,
    marker_name: str,
    tissue_type: str | None = None,
) -> MarkerQCResult:
    """
    Convenience function for marker validation.

    Parameters
    ----------
    intensities : np.ndarray
        Marker intensities
    marker_name : str
        Marker name
    tissue_type : str, optional
        Tissue type

    Returns
    -------
    MarkerQCResult
        Validation results
    """
    qc = MarkerQC()
    return qc.assess(intensities, marker_name, tissue_type)


def detect_crosstalk(
    intensities_a: np.ndarray,
    intensities_b: np.ndarray,
    marker_a: str,
    marker_b: str,
    correlation_threshold: float = 0.7,
) -> dict[str, Any]:
    """
    Detect potential crosstalk between two markers.

    High correlation between markers that shouldn't co-express
    may indicate spectral crosstalk.

    Parameters
    ----------
    intensities_a, intensities_b : np.ndarray
        Intensity values for each marker
    marker_a, marker_b : str
        Marker names
    correlation_threshold : float
        Correlation above which crosstalk is flagged

    Returns
    -------
    dict
        Crosstalk analysis results
    """
    # Pearson correlation
    correlation = np.corrcoef(intensities_a, intensities_b)[0, 1]

    # Check if markers should be mutually exclusive
    marker_a_key = marker_a.lower().replace("-", "").replace("_", "")
    marker_b_key = marker_b.lower().replace("-", "").replace("_", "")

    # Known exclusive pairs
    exclusive_pairs = [
        {"cd3", "cd20"},  # T-cells vs B-cells
        {"cd4", "cd8"},  # Helper vs cytotoxic T-cells
        {"panck", "cd45"},  # Epithelial vs immune
    ]

    should_be_exclusive = any({marker_a_key, marker_b_key} == pair for pair in exclusive_pairs)

    crosstalk_detected = correlation > correlation_threshold

    return {
        "marker_a": marker_a,
        "marker_b": marker_b,
        "correlation": float(correlation),
        "crosstalk_detected": crosstalk_detected,
        "should_be_exclusive": should_be_exclusive,
        "severity": (
            "high"
            if crosstalk_detected and should_be_exclusive
            else ("medium" if crosstalk_detected else "low")
        ),
    }


def assess_marker_specificity(
    marker_intensities: np.ndarray,
    nuclear_intensities: np.ndarray,
    marker_name: str,
) -> dict[str, float]:
    """
    Assess whether marker staining is specific or shows non-specific binding.

    Parameters
    ----------
    marker_intensities : np.ndarray
        Marker intensities per cell
    nuclear_intensities : np.ndarray
        Nuclear marker (e.g., DAPI) intensities
    marker_name : str
        Marker name

    Returns
    -------
    dict
        Specificity metrics
    """
    # Normalize marker name
    marker_key = marker_name.lower().replace("-", "").replace("_", "")
    reference = MARKER_REFERENCE.get(marker_key, {})

    # Expected location
    expected_location = reference.get("location", "unknown")

    # Correlation with nuclear marker
    nuclear_correlation = np.corrcoef(marker_intensities, nuclear_intensities)[0, 1]

    # If marker should be membrane/cytoplasm but correlates highly with nuclear...
    location_mismatch = False
    if expected_location in ["membrane", "cytoplasm"] and nuclear_correlation > 0.8:
        location_mismatch = True

    # Estimate positive fraction variance
    # Specific markers should have distinct positive populations
    p50 = np.median(marker_intensities)
    above_median = marker_intensities > p50
    below_median = marker_intensities <= p50

    separation = (
        np.mean(marker_intensities[above_median]) - np.mean(marker_intensities[below_median])
    ) / (np.std(marker_intensities) + 1e-8)

    specificity_score = np.clip(separation / 3.0, 0, 1)

    if location_mismatch:
        specificity_score *= 0.5

    return {
        "specificity_score": float(specificity_score),
        "nuclear_correlation": float(nuclear_correlation),
        "expected_location": expected_location,
        "location_mismatch": location_mismatch,
        "population_separation": float(separation),
    }
