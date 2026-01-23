"""
Quality Gate for Pre-Processing Validation

This module provides a quality gate that can be used before deconvolution
and EDF processing to detect and handle problematic images, preventing
the amplification of artifacts in downstream processing.

Usage:
    from kintsugi.qc import QualityGate

    # Create a gate with default thresholds
    gate = QualityGate()

    # Check a z-stack before processing
    result = gate.check_zstack(z_stack, channel="CD3", cycle="cyc01")

    if result.passed:
        # Safe to process
        deconvolved = decon.process_array(z_stack)
    else:
        print(f"Quality gate failed: {result.issues}")
        # Handle the failure (skip, flag, or attempt mitigation)
        if result.can_mitigate:
            mitigated_stack = gate.mitigate(z_stack, result)
            deconvolved = decon.process_array(mitigated_stack)

The gate checks for:
1. Stripe artifacts (horizontal banding from acquisition issues)
2. Corruption artifacts (organic distortion patterns)
3. Saturation (>50% pixels at max value)
4. Extremely low signal (potential blank/failed acquisition)
5. Tile grid patterns (stitching failures)
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Literal

import numpy as np

if TYPE_CHECKING:
    pass

logger = logging.getLogger(__name__)


@dataclass
class QualityGateResult:
    """Result of quality gate check."""

    passed: bool
    quality_score: float  # 0-1, higher is better
    issues: list[str] = field(default_factory=list)
    affected_planes: list[int] = field(default_factory=list)
    worst_severity: Literal["none", "mild", "moderate", "severe"] = "none"
    metrics: dict[str, float] = field(default_factory=dict)
    plane_results: dict[int, dict] = field(default_factory=dict)
    recommendations: list[str] = field(default_factory=list)

    @property
    def can_mitigate(self) -> bool:
        """Whether the issues can potentially be mitigated."""
        # Stripe artifacts can sometimes be mitigated, corruption cannot
        unmitigatable = ["corruption", "extreme_corruption", "blank"]
        return self.worst_severity != "severe" or not any(
            issue in unmitigatable for issue in self.issues
        )

    def __str__(self) -> str:
        if self.passed:
            return f"Quality gate PASSED (score={self.quality_score:.2f})"
        return f"Quality gate FAILED (score={self.quality_score:.2f}): " f"{'; '.join(self.issues)}"


class QualityGate:
    """
    Pre-processing quality gate for detecting problematic images.

    This gate should be used before deconvolution and EDF processing to
    detect images that would produce artifacts or failures. It integrates
    with the existing KINTSUGI QC infrastructure but provides a simple
    pass/fail interface suitable for batch processing.

    Parameters
    ----------
    max_saturation_pct : float
        Maximum percentage of pixels allowed at saturation (default 50%)
    min_signal_mean : float
        Minimum mean intensity relative to dtype max (default 0.001)
    stripe_threshold : float
        Threshold for stripe artifact detection (default 15.0)
    corruption_threshold : float
        Threshold for corruption detection (default 0.3)
    min_planes_required : int
        Minimum non-affected planes required to pass (default 3)

    Example
    -------
    >>> from kintsugi.qc import QualityGate
    >>> gate = QualityGate()
    >>> result = gate.check_zstack(stack, channel="CH3", cycle="cyc01")
    >>> if not result.passed:
    ...     print(f"Issues: {result.issues}")
    ...     print(f"Affected planes: {result.affected_planes}")
    """

    def __init__(
        self,
        max_saturation_pct: float = 50.0,
        min_signal_mean: float = 0.001,
        stripe_threshold: float = 15.0,
        corruption_threshold: float = 0.3,
        min_planes_required: int = 3,
    ):
        self.max_saturation_pct = max_saturation_pct
        self.min_signal_mean = min_signal_mean
        self.stripe_threshold = stripe_threshold
        self.corruption_threshold = corruption_threshold
        self.min_planes_required = min_planes_required

    def check_zstack(
        self,
        z_stack: np.ndarray,
        channel: str | None = None,
        cycle: str | None = None,
        verbose: bool = True,
    ) -> QualityGateResult:
        """
        Check a z-stack for quality issues before processing.

        Parameters
        ----------
        z_stack : np.ndarray
            3D image stack with shape (z, y, x)
        channel : str, optional
            Channel name for logging
        cycle : str, optional
            Cycle name for logging
        verbose : bool
            Whether to log detailed results

        Returns
        -------
        QualityGateResult
            Result with pass/fail status and detailed metrics
        """
        if z_stack.ndim != 3:
            raise ValueError(f"Expected 3D z-stack, got shape {z_stack.shape}")

        nz, height, width = z_stack.shape
        location = f"{cycle or ''}/{channel or ''}" if cycle or channel else "stack"

        issues = []
        affected_planes = []
        plane_results = {}
        metrics = {}
        recommendations = []

        # Get dtype info for saturation check
        if z_stack.dtype == np.uint16:
            max_val = 65535
        elif z_stack.dtype == np.uint8:
            max_val = 255
        else:
            max_val = z_stack.max()

        # === Check 1: Global saturation ===
        saturation_pct = 100.0 * np.sum(z_stack >= max_val * 0.99) / z_stack.size
        metrics["saturation_pct"] = saturation_pct
        if saturation_pct > self.max_saturation_pct:
            issues.append(f"saturation ({saturation_pct:.1f}%)")
            recommendations.append(
                "Check BaSiC flatfield minimum (BASIC_FLATFIELD_MIN) - may need to increase"
            )

        # === Check 2: Low signal / blank ===
        mean_intensity = z_stack.mean() / max_val
        metrics["mean_intensity_normalized"] = mean_intensity
        if mean_intensity < self.min_signal_mean:
            issues.append(f"blank/low_signal (mean={mean_intensity:.4f})")
            recommendations.append("Check raw data - may be failed acquisition")

        # === Check 3: Per-plane stripe artifact detection ===
        # Uses HP-FFT detection for improved sensitivity to faint patterns
        from kintsugi.qc.stripe_artifact import scan_zstack

        stripe_affected, stripe_results = scan_zstack(z_stack, use_hp_fft=True)
        for z_idx in stripe_affected:
            result = stripe_results[z_idx]
            plane_results[z_idx] = {
                "stripe_score": result.score,
                "stripe_severity": result.severity,
                "metrics": result.metrics,
            }
            if z_idx not in affected_planes:
                affected_planes.append(z_idx)

        # Track HP-FFT specific detections
        hp_fft_detected = []
        for z_idx, result in stripe_results.items():
            hp_row_sev = result.metrics.get("hp_row_severity", "clean")
            hp_col_sev = result.metrics.get("hp_col_severity", "clean")
            if hp_row_sev != "clean" or hp_col_sev != "clean":
                hp_fft_detected.append(z_idx)
                # Log HP-FFT specific metrics
                hp_row_peak = result.metrics.get("hp_row_peak", 0)
                hp_col_peak = result.metrics.get("hp_col_peak", 0)
                hp_row_z = result.metrics.get("hp_row_zscore", 0)
                hp_col_z = result.metrics.get("hp_col_zscore", 0)
                if z_idx not in plane_results:
                    plane_results[z_idx] = {}
                plane_results[z_idx]["hp_row_peak"] = hp_row_peak
                plane_results[z_idx]["hp_col_peak"] = hp_col_peak
                plane_results[z_idx]["hp_row_zscore"] = hp_row_z
                plane_results[z_idx]["hp_col_zscore"] = hp_col_z

        if hp_fft_detected:
            metrics["hp_fft_detected_planes"] = len(hp_fft_detected)

        if stripe_affected:
            worst_stripe = max(stripe_results[z].severity for z in stripe_affected)
            metrics["stripe_affected_planes"] = len(stripe_affected)
            metrics["stripe_worst_severity"] = worst_stripe

            # Check if any detections were HP-FFT only (faint patterns)
            hp_fft_only = set(hp_fft_detected) - set(stripe_affected)
            if hp_fft_only:
                issues.append(
                    f"stripe_artifact ({len(stripe_affected)} planes, {worst_stripe}; "
                    f"+{len(hp_fft_only)} faint via HP-FFT)"
                )
            else:
                issues.append(f"stripe_artifact ({len(stripe_affected)} planes, {worst_stripe})")

            if worst_stripe == "severe":
                recommendations.append(
                    "Severe stripe artifacts detected - consider excluding this stack"
                )
            else:
                recommendations.append(
                    "Stripe artifacts can be mitigated using interpolate_zplane() or "
                    "apply_directional_filter()"
                )

        # === Check 4: Corruption detection (organic distortion) ===
        corruption_scores = []
        for z_idx in range(nz):
            plane = z_stack[z_idx].astype(np.float32)
            corruption_score = self._detect_corruption(plane)
            corruption_scores.append(corruption_score)

            if corruption_score > self.corruption_threshold:
                if z_idx not in plane_results:
                    plane_results[z_idx] = {}
                plane_results[z_idx]["corruption_score"] = corruption_score
                if z_idx not in affected_planes:
                    affected_planes.append(z_idx)

        max_corruption = max(corruption_scores)
        metrics["max_corruption_score"] = max_corruption
        if max_corruption > self.corruption_threshold:
            n_corrupt = sum(1 for s in corruption_scores if s > self.corruption_threshold)
            if max_corruption > 0.6:
                issues.append(f"extreme_corruption ({n_corrupt} planes)")
                recommendations.append(
                    "Extreme corruption detected - this appears to be a data/acquisition "
                    "issue. Consider excluding this stack or reacquiring."
                )
            else:
                issues.append(f"corruption ({n_corrupt} planes)")
                recommendations.append(
                    "Moderate corruption detected - may be optical or sample issue"
                )

        # === Check 5: Tile grid pattern (stitching issues) ===
        grid_score = self._detect_tile_grid(z_stack.max(axis=0))  # Use MIP
        metrics["tile_grid_score"] = grid_score
        if grid_score > 0.8:
            issues.append(f"tile_grid_pattern (score={grid_score:.2f})")
            recommendations.append(
                "Tile grid pattern visible - check BaSiC flatfield and stitch blend sigma"
            )

        # === Determine overall result ===
        affected_planes = sorted(set(affected_planes))
        n_clean_planes = nz - len(affected_planes)

        # Calculate quality score
        quality_score = self._compute_quality_score(metrics, n_clean_planes, nz)
        metrics["quality_score"] = quality_score

        # Determine worst severity
        severity_order = {"none": 0, "mild": 1, "moderate": 2, "severe": 3}
        worst_severity = "none"
        for pr in plane_results.values():
            if "stripe_severity" in pr:
                if severity_order.get(pr["stripe_severity"], 0) > severity_order.get(
                    worst_severity, 0
                ):
                    worst_severity = pr["stripe_severity"]
            if pr.get("corruption_score", 0) > 0.6:
                worst_severity = "severe"
            elif pr.get("corruption_score", 0) > self.corruption_threshold:
                if severity_order.get("moderate", 0) > severity_order.get(worst_severity, 0):
                    worst_severity = "moderate"

        # Pass/fail decision
        passed = len(issues) == 0 or (
            quality_score > 0.7 and n_clean_planes >= self.min_planes_required
        )

        # Adjust pass for severe issues
        if "extreme_corruption" in str(issues) or "blank" in str(issues):
            passed = False

        result = QualityGateResult(
            passed=passed,
            quality_score=quality_score,
            issues=issues,
            affected_planes=affected_planes,
            worst_severity=worst_severity,
            metrics=metrics,
            plane_results=plane_results,
            recommendations=recommendations,
        )

        if verbose:
            logger.info(f"{location}: {result}")
            if not passed:
                for rec in recommendations:
                    logger.info(f"  → {rec}")

        return result

    def check_single_plane(
        self,
        image: np.ndarray,
        channel: str | None = None,
        cycle: str | None = None,
    ) -> QualityGateResult:
        """
        Check a single 2D image for quality issues.

        Parameters
        ----------
        image : np.ndarray
            2D image
        channel : str, optional
            Channel name for logging
        cycle : str, optional
            Cycle name for logging

        Returns
        -------
        QualityGateResult
            Result with pass/fail status
        """
        # Wrap in 3D and check
        stack_3d = image[np.newaxis, :, :]
        return self.check_zstack(stack_3d, channel, cycle, verbose=False)

    def _detect_corruption(self, image: np.ndarray) -> float:
        """
        Detect organic/liquid-like corruption patterns.

        Uses gradient coherence analysis - corrupted images show unusual
        gradient patterns that look organic/flowing rather than cellular.

        Returns score 0-1, higher = more likely corrupted.
        """
        # Downsample for speed
        step = max(1, min(image.shape) // 500)
        small = image[::step, ::step]

        # Compute gradients
        gy = np.diff(small, axis=0)
        gx = np.diff(small, axis=1)

        # Gradient magnitude and direction
        # Crop to same size
        gy = gy[:, :-1]
        gx = gx[:-1, :]
        mag = np.sqrt(gx**2 + gy**2)

        if mag.mean() < 1e-6:
            return 0.0

        # Compute gradient direction coherence
        # Corrupted images often show smooth, flowing gradient patterns
        with np.errstate(divide="ignore", invalid="ignore"):
            direction = np.arctan2(gy, gx)

        # Local direction variance (low = suspicious coherent flow)
        kernel_size = 5
        from scipy.ndimage import uniform_filter

        dir_local_mean = uniform_filter(direction, size=kernel_size)
        dir_local_var = uniform_filter((direction - dir_local_mean) ** 2, size=kernel_size)

        # Corrupted: low variance (coherent) but high magnitude
        coherence = 1.0 - np.clip(dir_local_var.mean() / (np.pi**2 / 3), 0, 1)

        # Also check for unusual spatial frequency distribution
        # Corrupted images often have unusual low-frequency dominance
        from scipy.fft import fft2, fftshift

        fft = fftshift(fft2(small))
        power = np.abs(fft) ** 2

        h, w = power.shape
        cy, cx = h // 2, w // 2

        # Ratio of energy in low vs high frequencies
        r_low = min(h, w) // 8
        r_high = min(h, w) // 4

        y, x = np.ogrid[:h, :w]
        low_mask = ((y - cy) ** 2 + (x - cx) ** 2) < r_low**2
        high_mask = ((y - cy) ** 2 + (x - cx) ** 2) > r_high**2

        low_energy = power[low_mask].sum()
        high_energy = power[high_mask].sum()

        if high_energy > 0:
            freq_ratio = low_energy / (high_energy + low_energy)
        else:
            freq_ratio = 1.0

        # Combine metrics
        # Normal images: coherence < 0.3, freq_ratio < 0.7
        # Corrupted: coherence > 0.4 or freq_ratio > 0.85
        corruption_score = coherence * 0.5 + max(0, freq_ratio - 0.5) * 1.0
        return float(np.clip(corruption_score, 0, 1))

    def _detect_tile_grid(self, image: np.ndarray) -> float:
        """
        Detect visible tile grid patterns from stitching.

        Returns score 0-1, higher = more visible grid.
        """
        # Compute row and column intensity profiles
        row_means = np.mean(image, axis=1)
        col_means = np.mean(image, axis=0)

        # Look for periodic patterns
        from scipy.fft import fft

        row_fft = np.abs(fft(row_means - row_means.mean()))
        col_fft = np.abs(fft(col_means - col_means.mean()))

        # Find peaks (excluding DC)
        n = len(row_fft) // 2
        row_peaks = row_fft[1:n]
        col_peaks = col_fft[1:n]

        # Check for prominent periodic components
        row_periodicity = row_peaks.max() / (row_peaks.mean() + 1e-10)
        col_periodicity = col_peaks.max() / (col_peaks.mean() + 1e-10)

        # Also check std/mean ratio (from stitched-image-qc skill)
        std_mean_ratio = image.std() / (image.mean() + 1e-10)

        # Combine: high periodicity OR high std/mean ratio
        grid_score = max(
            (row_periodicity - 3) / 10,
            (col_periodicity - 3) / 10,
            (std_mean_ratio - 0.5) / 0.5,
        )
        return float(np.clip(grid_score, 0, 1))

    def _compute_quality_score(
        self,
        metrics: dict[str, float],
        n_clean_planes: int,
        total_planes: int,
    ) -> float:
        """Compute overall quality score 0-1."""
        scores = []

        # Saturation penalty
        sat_pct = metrics.get("saturation_pct", 0)
        scores.append(max(0, 1 - sat_pct / self.max_saturation_pct))

        # Signal level
        mean_norm = metrics.get("mean_intensity_normalized", 0)
        scores.append(min(1, mean_norm / (self.min_signal_mean * 10)))

        # Corruption penalty
        corruption = metrics.get("max_corruption_score", 0)
        scores.append(max(0, 1 - corruption / 0.5))

        # Clean plane ratio
        if total_planes > 0:
            scores.append(n_clean_planes / total_planes)

        # Tile grid penalty
        grid = metrics.get("tile_grid_score", 0)
        scores.append(max(0, 1 - grid))

        return float(np.mean(scores))

    def mitigate(
        self,
        z_stack: np.ndarray,
        result: QualityGateResult,
        method: str = "auto",
    ) -> np.ndarray:
        """
        Attempt to mitigate detected issues in a z-stack.

        Parameters
        ----------
        z_stack : np.ndarray
            Original z-stack
        result : QualityGateResult
            Result from check_zstack()
        method : str
            Mitigation method: 'auto', 'interpolate', 'filter', or 'exclude'

        Returns
        -------
        np.ndarray
            Mitigated z-stack (may have fewer planes if exclusion was used)
        """
        if result.passed:
            return z_stack

        if not result.can_mitigate:
            logger.warning(
                "Issues cannot be mitigated - returning original stack. "
                "Consider excluding this data."
            )
            return z_stack

        from kintsugi.qc.stripe_mitigation import (
            apply_directional_filter,
            interpolate_zplane,
        )

        mitigated = z_stack.copy()

        for z_idx in result.affected_planes:
            plane_info = result.plane_results.get(z_idx, {})

            # Choose method based on issue type
            if method == "auto":
                if plane_info.get("corruption_score", 0) > self.corruption_threshold:
                    # Corruption: use z-interpolation if possible
                    chosen_method = "interpolate"
                elif plane_info.get("stripe_severity") in ("moderate", "severe"):
                    # Severe stripe: use directional filter
                    chosen_method = "filter"
                else:
                    # Mild stripe: try interpolation first
                    chosen_method = "interpolate"
            else:
                chosen_method = method

            # Apply mitigation
            if chosen_method == "interpolate":
                try:
                    mitigated_result = interpolate_zplane(mitigated, z_idx)
                    mitigated[z_idx] = mitigated_result.result
                    logger.debug(f"Z-plane {z_idx}: interpolated")
                except Exception as e:
                    logger.warning(f"Z-plane {z_idx}: interpolation failed ({e}), trying filter")
                    chosen_method = "filter"

            if chosen_method == "filter":
                try:
                    mitigated_result = apply_directional_filter(mitigated[z_idx])
                    mitigated[z_idx] = mitigated_result.result
                    logger.debug(f"Z-plane {z_idx}: directional filter applied")
                except Exception as e:
                    logger.warning(f"Z-plane {z_idx}: filter failed ({e})")

        return mitigated


def check_before_processing(
    z_stack: np.ndarray,
    channel: str | None = None,
    cycle: str | None = None,
    fail_action: Literal["raise", "warn", "skip", "mitigate"] = "warn",
) -> tuple[np.ndarray | None, QualityGateResult]:
    """
    Convenience function to check quality before processing.

    Parameters
    ----------
    z_stack : np.ndarray
        3D image stack
    channel : str, optional
        Channel name
    cycle : str, optional
        Cycle name
    fail_action : str
        What to do on failure:
        - 'raise': Raise ValueError
        - 'warn': Log warning and return original
        - 'skip': Return None
        - 'mitigate': Attempt mitigation and return result

    Returns
    -------
    tuple
        (processed_stack, result) where processed_stack may be None if skipped,
        the original if passed/warned, or mitigated if that action was chosen
    """
    gate = QualityGate()
    result = gate.check_zstack(z_stack, channel, cycle)

    if result.passed:
        return z_stack, result

    if fail_action == "raise":
        raise ValueError(
            f"Quality gate failed for {cycle or ''}/{channel or ''}: " f"{'; '.join(result.issues)}"
        )
    elif fail_action == "skip":
        logger.warning(f"Skipping {cycle or ''}/{channel or ''}: {result}")
        return None, result
    elif fail_action == "mitigate":
        mitigated = gate.mitigate(z_stack, result)
        return mitigated, result
    else:  # warn
        logger.warning(f"{result}")
        return z_stack, result
