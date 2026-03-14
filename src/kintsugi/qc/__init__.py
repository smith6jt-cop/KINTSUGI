"""
KINTSUGI Quality Control Module

Pure Python implementations of quality control algorithms for
multiplex immunofluorescence imaging.

Features:
- Image-level quality assessment
- Cell-level quality control
- Artifact detection (stripe/banding)
- Batch effect detection
- Marker quality validation

Artifact Detection:
    # Project-level scanning (recommended)
    from kintsugi.project import KintsugiProject
    from kintsugi.qc import ArtifactScanner

    project = KintsugiProject.load("/path/to/project")
    scanner = ArtifactScanner(project)
    report = scanner.scan_raw_data()
    print(report.summary())

    # Low-level z-stack scanning
    from kintsugi.qc import scan_zstack
    affected, results = scan_zstack(z_stack)

Image QC:
    from kintsugi.qc import ImageQC, assess_image_quality
    qc = ImageQC()
    result = qc.assess(image, marker="CD3", tissue="tonsil")

Cell QC:
    from kintsugi.qc import CellQC
    cell_qc = CellQC()
    filtered_cells = cell_qc.filter_outliers(cell_data)
"""

from kintsugi.qc.artifact_scanner import (
    ArtifactItem,
    ArtifactReport,
    ArtifactScanner,
    ChannelScanResult,
    detect_general_artifacts,
    scan_project_artifacts,
)
from kintsugi.qc.batch_qc import (
    BatchQC,
    compute_batch_statistics,
    detect_batch_effects,
)
from kintsugi.qc.cell_qc import (
    CellQC,
    detect_cycle_dropout,
    detect_doublets,
    detect_outliers,
    filter_by_intensity,
    filter_by_morphology,
    filter_cells_pipeline,
    prune_marker_outliers,
)
from kintsugi.qc.image_qc import (
    ImageQC,
    assess_image_quality,
    compute_focus_score,
    compute_tissue_coverage,
    detect_artifacts,
)
from kintsugi.qc.marker_qc import (
    MarkerQC,
    assess_marker_specificity,
    detect_crosstalk,
    validate_marker_expression,
)
from kintsugi.qc.quality_gate import (
    QualityGate,
    QualityGateResult,
    check_before_processing,
)
from kintsugi.qc.quick_metrics import (
    QuickQCResult,
    average_with_neighbors,
    compute_quick_qc_metrics,
    fix_bad_zplanes_in_tiff,
    fix_bad_zplanes_in_zarr,
)
from kintsugi.qc.signal_quality import compute_signal_isolation_quality
from kintsugi.qc.stripe_artifact import (
    StripeArtifactResult,
    classify_hp_fft_severity,
    compute_hp_fft_peak_strength,
    compute_zstack_baseline,
    detect_stripe_artifact,
    scan_zstack,
)
from kintsugi.qc.stripe_mitigation import (
    MitigationResult,
    apply_directional_filter,
    apply_notch_filter,
    get_recommended_method,
    interpolate_zplane,
    mitigate_artifact,
)


# Backwards compatibility aliases (deprecated)
def detect_stripe_artifacts(image, sensitivity=0.5, **kwargs):
    """Deprecated: Use detect_stripe_artifact() instead."""
    return detect_stripe_artifact(image, **kwargs)


def scan_zstack_for_artifacts(z_stack, sensitivity=0.5, verbose=True, **kwargs):
    """Deprecated: Use scan_zstack() instead."""
    affected, results = scan_zstack(z_stack, **kwargs)
    # Return a compatible report object
    from dataclasses import dataclass, field
    from typing import Literal

    @dataclass
    class _LegacyReport:
        total_planes: int
        affected_planes: list
        results: dict
        worst_severity: str
        summary: str

        @property
        def has_artifacts(self):
            return len(self.affected_planes) > 0

        def __str__(self):
            if not self.has_artifacts:
                return f"No artifacts detected in {self.total_planes} z-planes"
            return f"Artifacts in {len(self.affected_planes)}/{self.total_planes} planes"

    worst = "none"
    severity_order = {"none": 0, "mild": 1, "moderate": 2, "severe": 3}
    for r in results.values():
        if severity_order.get(r.severity, 0) > severity_order.get(worst, 0):
            worst = r.severity

    if verbose and affected:
        for z in affected:
            print(f"  Z-plane {z}: {results[z]}")

    return _LegacyReport(
        total_planes=len(results),
        affected_planes=affected,
        results=results,
        worst_severity=worst,
        summary=f"Detected {len(affected)}/{len(results)} affected planes",
    )


__all__ = [
    # Quick QC Metrics (inline stitching pipeline)
    "QuickQCResult",
    "compute_quick_qc_metrics",
    "average_with_neighbors",
    "fix_bad_zplanes_in_zarr",
    "fix_bad_zplanes_in_tiff",
    # Quality Gate (pre-processing validation)
    "QualityGate",
    "QualityGateResult",
    "check_before_processing",
    # Artifact Scanner (primary interface)
    "ArtifactScanner",
    "ArtifactReport",
    "ArtifactItem",
    "ChannelScanResult",
    "scan_project_artifacts",
    "detect_general_artifacts",
    # Stripe Detection (low-level)
    "StripeArtifactResult",
    "detect_stripe_artifact",
    "scan_zstack",
    "compute_zstack_baseline",
    "compute_hp_fft_peak_strength",
    "classify_hp_fft_severity",
    # Deprecated aliases (backwards compatibility)
    "detect_stripe_artifacts",
    "scan_zstack_for_artifacts",
    # Stripe Mitigation
    "MitigationResult",
    "apply_notch_filter",
    "apply_directional_filter",
    "interpolate_zplane",
    "mitigate_artifact",
    "get_recommended_method",
    # Image QC
    "ImageQC",
    "assess_image_quality",
    "detect_artifacts",
    "compute_focus_score",
    "compute_tissue_coverage",
    # Cell QC
    "CellQC",
    "detect_outliers",
    "filter_by_intensity",
    "filter_by_morphology",
    "detect_doublets",
    "detect_cycle_dropout",
    "prune_marker_outliers",
    "filter_cells_pipeline",
    # Marker QC
    "MarkerQC",
    "validate_marker_expression",
    "detect_crosstalk",
    "assess_marker_specificity",
    # Batch QC
    "BatchQC",
    "detect_batch_effects",
    "compute_batch_statistics",
    # Signal isolation quality
    "compute_signal_isolation_quality",
]
