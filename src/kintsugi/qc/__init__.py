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
    scan_project_artifacts,
)
from kintsugi.qc.batch_qc import (
    BatchQC,
    compute_batch_statistics,
    detect_batch_effects,
)
from kintsugi.qc.cell_qc import (
    CellQC,
    detect_doublets,
    detect_outliers,
    filter_by_intensity,
    filter_by_morphology,
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
from kintsugi.qc.stripe_artifact import (
    StripeArtifactResult,
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

__all__ = [
    # Artifact Scanner (primary interface)
    "ArtifactScanner",
    "ArtifactReport",
    "ArtifactItem",
    "ChannelScanResult",
    "scan_project_artifacts",
    # Stripe Detection (low-level)
    "StripeArtifactResult",
    "detect_stripe_artifact",
    "scan_zstack",
    "compute_zstack_baseline",
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
    # Marker QC
    "MarkerQC",
    "validate_marker_expression",
    "detect_crosstalk",
    "assess_marker_specificity",
    # Batch QC
    "BatchQC",
    "detect_batch_effects",
    "compute_batch_statistics",
]
