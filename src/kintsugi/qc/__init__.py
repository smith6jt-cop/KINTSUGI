"""
KINTSUGI Quality Control Module

Pure Python implementations of quality control algorithms for
multiplex immunofluorescence imaging, inspired by CyLinter.

Features:
- Image-level quality assessment
- Cell-level quality control
- Outlier detection
- Artifact identification
- Batch effect detection
- Marker quality validation

Usage:
    from kintsugi.qc import (
        ImageQC,
        CellQC,
        MarkerQC,
        detect_outliers,
        assess_image_quality,
    )

    # Image-level QC
    qc = ImageQC()
    result = qc.assess(image, marker="CD3", tissue="tonsil")

    # Cell-level QC
    cell_qc = CellQC()
    filtered_cells = cell_qc.filter_outliers(cell_data)
"""

from kintsugi.qc.image_qc import (
    ImageQC,
    assess_image_quality,
    detect_artifacts,
    compute_focus_score,
    compute_tissue_coverage,
)

from kintsugi.qc.cell_qc import (
    CellQC,
    detect_outliers,
    filter_by_intensity,
    filter_by_morphology,
    detect_doublets,
)

from kintsugi.qc.marker_qc import (
    MarkerQC,
    validate_marker_expression,
    detect_crosstalk,
    assess_marker_specificity,
)

from kintsugi.qc.batch_qc import (
    BatchQC,
    detect_batch_effects,
    compute_batch_statistics,
)

__all__ = [
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
