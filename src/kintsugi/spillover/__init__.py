"""
KINTSUGI Spillover Correction Module

Post-segmentation lateral signal spillover correction using REDSEA
(REinforcement Dynamic Spillover EliminAtion). Corrects for signal from
neighboring cells bleeding into a cell's segmentation region due to membrane
interleaving, cytoplasm overlap in z, and optical PSF spreading.

REDSEA operates after segmentation and produces corrected per-cell
quantification tables that replace standard feature extraction intensity
values for downstream clustering and spatial analysis.

Main Functions:
- run_redsea_correction: Core REDSEA spillover correction
- estimate_element_size: Estimate optimal element size from cell mask
- identify_surface_markers: Identify surface markers for correction
- compute_spillover_metrics: QC metrics for correction evaluation

Example:
    >>> from kintsugi.spillover import run_redsea_correction, estimate_element_size
    >>>
    >>> # Estimate parameters
    >>> element_size = estimate_element_size(cell_mask)
    >>>
    >>> # Run correction
    >>> result = run_redsea_correction(
    ...     multichannel_image="image.ome.tif",
    ...     segmentation_mask=cell_mask,
    ...     marker_names=["CD3", "CD20", "CD8", "DAPI"],
    ...     element_size=element_size,
    ... )
    >>> corrected_df = result.corrected

Reference:
    Bai Y et al. Adjacent Cell Marker Lateral Spillover Compensation and
    Reinforcement for Multiplexed Images. Front. Immunol. 12:652631 (2021).
    PMID: 34295327
"""

from .redsea_wrapper import SpilloverResult, run_redsea_correction
from .spillover_qc import (
    DEFAULT_EXCLUSIVE_PAIRS,
    compute_spillover_metrics,
    plot_spillover_comparison,
    plot_spillover_summary_heatmap,
)
from .utils import (
    KNOWN_NUCLEAR_MARKERS,
    KNOWN_SURFACE_MARKERS,
    estimate_element_size,
    identify_surface_markers,
    load_markers_csv,
    validate_segmentation_mask,
    write_markers_csv,
)

__all__ = [
    # Core correction
    "run_redsea_correction",
    "SpilloverResult",
    # Parameter estimation
    "estimate_element_size",
    "identify_surface_markers",
    # QC and visualization
    "compute_spillover_metrics",
    "plot_spillover_comparison",
    "plot_spillover_summary_heatmap",
    "DEFAULT_EXCLUSIVE_PAIRS",
    # Utilities
    "validate_segmentation_mask",
    "write_markers_csv",
    "load_markers_csv",
    # Constants
    "KNOWN_SURFACE_MARKERS",
    "KNOWN_NUCLEAR_MARKERS",
]
