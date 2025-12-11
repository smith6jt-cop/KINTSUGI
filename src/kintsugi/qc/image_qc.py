"""
Image-Level Quality Control

Quality assessment at the whole-image level for fluorescence microscopy.
Inspired by CyLinter's image QC approaches.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Union

import numpy as np
from scipy import ndimage

if TYPE_CHECKING:
    import dask.array

logger = logging.getLogger("kintsugi.qc.image_qc")

# Type alias
ArrayLike = Union[np.ndarray, "dask.array.Array"]


@dataclass
class ImageQCResult:
    """Results from image-level quality control."""

    # Overall assessment
    passed: bool
    quality_score: float
    issues: list[str] = field(default_factory=list)

    # Detailed metrics
    metrics: dict[str, float] = field(default_factory=dict)

    # Artifact regions (x, y, w, h)
    artifact_regions: list[tuple[int, int, int, int]] = field(default_factory=list)

    # Recommendations
    recommendations: list[str] = field(default_factory=list)


def _to_numpy(arr: ArrayLike) -> np.ndarray:
    """Convert to numpy, computing if dask."""
    if hasattr(arr, "compute"):
        return arr.compute()
    return np.asarray(arr)


class ImageQC:
    """
    Comprehensive image-level quality control for fluorescence microscopy.

    Assesses:
    - Focus quality (sharpness)
    - Tissue coverage
    - Intensity distribution
    - Artifacts (hot pixels, bands, etc.)
    - Illumination uniformity
    - Dynamic range utilization

    Example
    -------
    >>> qc = ImageQC()
    >>> result = qc.assess(image)
    >>> if not result.passed:
    ...     print(f"Issues: {result.issues}")
    ...     print(f"Recommendations: {result.recommendations}")
    """

    def __init__(
        self,
        min_focus_score: float = 0.3,
        min_tissue_coverage: float = 0.1,
        max_saturation_pct: float = 0.01,
        min_dynamic_range_pct: float = 0.1,
        max_artifact_area_pct: float = 0.05,
    ):
        """
        Initialize Image QC.

        Parameters
        ----------
        min_focus_score : float
            Minimum acceptable focus score (0-1)
        min_tissue_coverage : float
            Minimum tissue coverage fraction
        max_saturation_pct : float
            Maximum acceptable saturated pixel percentage
        min_dynamic_range_pct : float
            Minimum dynamic range utilization
        max_artifact_area_pct : float
            Maximum acceptable artifact area percentage
        """
        self.min_focus_score = min_focus_score
        self.min_tissue_coverage = min_tissue_coverage
        self.max_saturation_pct = max_saturation_pct
        self.min_dynamic_range_pct = min_dynamic_range_pct
        self.max_artifact_area_pct = max_artifact_area_pct

    def assess(
        self,
        image: ArrayLike,
        marker: str | None = None,
        tissue: str | None = None,
    ) -> ImageQCResult:
        """
        Perform comprehensive quality assessment on an image.

        Parameters
        ----------
        image : array-like
            2D fluorescence image
        marker : str, optional
            Marker name for marker-specific thresholds
        tissue : str, optional
            Tissue type for tissue-specific thresholds

        Returns
        -------
        ImageQCResult
            Detailed QC results
        """
        img = _to_numpy(image)
        issues = []
        recommendations = []
        metrics = {}

        # 1. Focus quality
        focus_score = compute_focus_score(img)
        metrics["focus_score"] = focus_score
        if focus_score < self.min_focus_score:
            issues.append(f"Poor focus (score={focus_score:.3f})")
            recommendations.append("Consider refocusing or adjusting Z-position")

        # 2. Tissue coverage
        coverage = compute_tissue_coverage(img)
        metrics["tissue_coverage"] = coverage
        if coverage < self.min_tissue_coverage:
            issues.append(f"Low tissue coverage ({coverage*100:.1f}%)")
            recommendations.append("Check tissue mounting or field of view selection")

        # 3. Saturation
        saturation_metrics = self._check_saturation(img)
        metrics.update(saturation_metrics)
        if saturation_metrics["saturated_pct"] > self.max_saturation_pct:
            issues.append(f"High saturation ({saturation_metrics['saturated_pct']*100:.2f}%)")
            recommendations.append("Reduce exposure time or laser power")

        # 4. Dynamic range
        dr_metrics = self._check_dynamic_range(img)
        metrics.update(dr_metrics)
        if dr_metrics["dynamic_range_utilization"] < self.min_dynamic_range_pct:
            issues.append(
                f"Poor dynamic range utilization ({dr_metrics['dynamic_range_utilization']*100:.1f}%)"
            )
            recommendations.append("Increase exposure or check staining quality")

        # 5. Artifacts
        artifact_regions, artifact_metrics = detect_artifacts(img)
        metrics.update(artifact_metrics)
        artifact_area = sum(w * h for _, _, w, h in artifact_regions)
        artifact_pct = artifact_area / (img.shape[0] * img.shape[1])
        if artifact_pct > self.max_artifact_area_pct:
            issues.append(f"Significant artifacts detected ({artifact_pct*100:.1f}%)")
            recommendations.append("Review acquisition settings or clean optical path")

        # 6. Illumination uniformity
        uniformity = self._check_illumination_uniformity(img)
        metrics["illumination_uniformity"] = uniformity
        if uniformity < 0.7:
            issues.append(f"Non-uniform illumination (uniformity={uniformity:.2f})")
            recommendations.append("Apply illumination correction or check light source")

        # 7. SNR
        snr = self._compute_snr(img)
        metrics["snr"] = snr
        if snr < 3.0:
            issues.append(f"Low SNR ({snr:.1f})")
            recommendations.append("Increase signal or reduce noise (longer exposure, averaging)")

        # Overall quality score
        quality_score = self._compute_overall_score(metrics)
        metrics["overall_quality_score"] = quality_score

        # Pass/fail decision
        passed = len(issues) == 0 and quality_score > 0.5

        return ImageQCResult(
            passed=passed,
            quality_score=quality_score,
            issues=issues,
            metrics=metrics,
            artifact_regions=artifact_regions,
            recommendations=recommendations,
        )

    def _check_saturation(self, image: np.ndarray) -> dict[str, float]:
        """Check for saturated pixels."""
        if image.dtype == np.uint16:
            max_val = 65535
        elif image.dtype == np.uint8:
            max_val = 255
        else:
            max_val = image.max()

        saturated = np.sum(image >= max_val * 0.99)
        zero = np.sum(image == 0)
        total = image.size

        return {
            "saturated_pct": saturated / total,
            "zero_pct": zero / total,
            "saturated_count": int(saturated),
        }

    def _check_dynamic_range(self, image: np.ndarray) -> dict[str, float]:
        """Check dynamic range utilization."""
        if image.dtype == np.uint16:
            max_possible = 65535
        elif image.dtype == np.uint8:
            max_possible = 255
        else:
            max_possible = 1.0

        p1, p99 = np.percentile(image, [1, 99])
        used_range = p99 - p1
        utilization = used_range / max_possible

        return {
            "intensity_p1": float(p1),
            "intensity_p99": float(p99),
            "used_dynamic_range": float(used_range),
            "dynamic_range_utilization": float(utilization),
        }

    def _check_illumination_uniformity(self, image: np.ndarray) -> float:
        """
        Check illumination uniformity across the image.

        Returns value between 0 (non-uniform) and 1 (perfectly uniform).
        """
        # Downsample for efficiency
        step = max(1, min(image.shape) // 100)
        small = image[::step, ::step].astype(np.float64)

        # Compute local means in a grid
        grid_size = min(10, small.shape[0] // 10, small.shape[1] // 10)
        if grid_size < 3:
            return 1.0

        h, w = small.shape
        tile_h, tile_w = h // grid_size, w // grid_size

        means = []
        for i in range(grid_size):
            for j in range(grid_size):
                tile = small[i * tile_h : (i + 1) * tile_h, j * tile_w : (j + 1) * tile_w]
                means.append(np.mean(tile))

        means = np.array(means)
        if means.max() == 0:
            return 1.0

        # Uniformity as ratio of min to max mean
        uniformity = means.min() / (means.max() + 1e-8)
        return float(np.clip(uniformity, 0, 1))

    def _compute_snr(self, image: np.ndarray) -> float:
        """Compute Signal-to-Noise Ratio."""
        p10 = np.percentile(image, 10)
        p90 = np.percentile(image, 90)

        # Background is low-intensity region
        background = image[image < p10]
        if len(background) < 100:
            background = image.ravel()[:1000]

        noise_std = np.std(background)
        signal = p90 - p10

        if noise_std > 0:
            return float(signal / noise_std)
        return 0.0

    def _compute_overall_score(self, metrics: dict[str, float]) -> float:
        """Compute weighted overall quality score."""
        weights = {
            "focus_score": 0.25,
            "tissue_coverage": 0.15,
            "dynamic_range_utilization": 0.15,
            "illumination_uniformity": 0.15,
            "snr": 0.20,
        }

        # Normalize SNR to 0-1 (assuming 10 is excellent)
        snr_norm = min(1.0, metrics.get("snr", 0) / 10.0)

        scores = {
            "focus_score": metrics.get("focus_score", 0),
            "tissue_coverage": min(1.0, metrics.get("tissue_coverage", 0) / 0.5),
            "dynamic_range_utilization": min(
                1.0, metrics.get("dynamic_range_utilization", 0) / 0.3
            ),
            "illumination_uniformity": metrics.get("illumination_uniformity", 0),
            "snr": snr_norm,
        }

        # Penalty for saturation
        saturation_penalty = max(0, metrics.get("saturated_pct", 0) - self.max_saturation_pct) * 10

        weighted_sum = sum(scores.get(k, 0) * v for k, v in weights.items())
        total_weight = sum(weights.values())

        return float(np.clip(weighted_sum / total_weight - saturation_penalty, 0, 1))


def assess_image_quality(
    image: ArrayLike,
    marker: str | None = None,
    tissue: str | None = None,
) -> ImageQCResult:
    """
    Convenience function for quick image quality assessment.

    Parameters
    ----------
    image : array-like
        2D fluorescence image
    marker : str, optional
        Marker name
    tissue : str, optional
        Tissue type

    Returns
    -------
    ImageQCResult
        Quality assessment results
    """
    qc = ImageQC()
    return qc.assess(image, marker=marker, tissue=tissue)


def compute_focus_score(image: ArrayLike) -> float:
    """
    Compute focus quality score using Laplacian variance.

    Higher values indicate sharper focus.
    Normalized to [0, 1] range.

    Parameters
    ----------
    image : array-like
        2D image

    Returns
    -------
    float
        Focus score (0-1)
    """
    img = _to_numpy(image).astype(np.float64)

    # Downsample large images
    if img.shape[0] > 2000 or img.shape[1] > 2000:
        step = max(img.shape) // 2000
        img = img[::step, ::step]

    # Laplacian variance
    laplacian = ndimage.laplace(img)
    variance = np.var(laplacian)

    # Normalize (empirically, variance > 1000 is sharp for 16-bit)
    if img.max() > 255:
        norm_factor = 1000
    else:
        norm_factor = 100

    score = np.tanh(variance / norm_factor)
    return float(score)


def compute_tissue_coverage(image: ArrayLike, threshold_percentile: float = 10) -> float:
    """
    Estimate tissue coverage as fraction of non-background pixels.

    Parameters
    ----------
    image : array-like
        2D image
    threshold_percentile : float
        Percentile above which pixels are considered tissue

    Returns
    -------
    float
        Fraction of image covered by tissue (0-1)
    """
    img = _to_numpy(image)

    threshold = np.percentile(img, threshold_percentile)
    tissue_pixels = np.sum(img > threshold)
    total_pixels = img.size

    return float(tissue_pixels / total_pixels)


def detect_artifacts(
    image: ArrayLike,
    hot_pixel_threshold: float = 5.0,
    band_detection: bool = True,
) -> tuple[list[tuple[int, int, int, int]], dict[str, float]]:
    """
    Detect various artifacts in the image.

    Detects:
    - Hot pixels (extremely bright single pixels)
    - Dark pixels (extremely dark single pixels)
    - Banding artifacts
    - Dust/debris

    Parameters
    ----------
    image : array-like
        2D image
    hot_pixel_threshold : float
        Standard deviations above median to flag as hot pixel
    band_detection : bool
        Whether to detect horizontal/vertical banding

    Returns
    -------
    tuple
        List of artifact regions (x, y, w, h) and metrics dict
    """
    img = _to_numpy(image).astype(np.float64)
    artifact_regions = []
    metrics = {}

    # Hot pixel detection
    median = np.median(img)
    mad = np.median(np.abs(img - median)) * 1.4826  # MAD to std estimate
    hot_threshold = median + hot_pixel_threshold * mad

    hot_pixels = img > hot_threshold
    hot_count = np.sum(hot_pixels)
    metrics["hot_pixel_count"] = int(hot_count)
    metrics["hot_pixel_pct"] = float(hot_count / img.size)

    # Find connected hot pixel regions
    if hot_count > 0:
        labeled, num_features = ndimage.label(hot_pixels)
        for i in range(1, num_features + 1):
            coords = np.where(labeled == i)
            if len(coords[0]) > 0:
                y_min, y_max = coords[0].min(), coords[0].max()
                x_min, x_max = coords[1].min(), coords[1].max()
                artifact_regions.append((x_min, y_min, x_max - x_min + 1, y_max - y_min + 1))

    # Band detection
    if band_detection:
        h_variance, v_variance = _detect_banding(img)
        metrics["horizontal_banding_score"] = float(h_variance)
        metrics["vertical_banding_score"] = float(v_variance)

        if h_variance > 0.1:
            # Find rows with significant banding
            row_means = np.mean(img, axis=1)
            row_median = np.median(row_means)
            row_mad = np.median(np.abs(row_means - row_median)) * 1.4826
            bad_rows = np.where(np.abs(row_means - row_median) > 3 * row_mad)[0]
            for row in bad_rows:
                artifact_regions.append((0, int(row), img.shape[1], 1))

        if v_variance > 0.1:
            col_means = np.mean(img, axis=0)
            col_median = np.median(col_means)
            col_mad = np.median(np.abs(col_means - col_median)) * 1.4826
            bad_cols = np.where(np.abs(col_means - col_median) > 3 * col_mad)[0]
            for col in bad_cols:
                artifact_regions.append((int(col), 0, 1, img.shape[0]))

    return artifact_regions, metrics


def _detect_banding(image: np.ndarray) -> tuple[float, float]:
    """Detect horizontal and vertical banding artifacts."""
    # Row means variation
    row_means = np.mean(image, axis=1)
    row_diff = np.diff(row_means)
    h_variance = np.std(row_diff) / (np.mean(row_means) + 1e-8)

    # Column means variation
    col_means = np.mean(image, axis=0)
    col_diff = np.diff(col_means)
    v_variance = np.std(col_diff) / (np.mean(col_means) + 1e-8)

    return float(h_variance), float(v_variance)
