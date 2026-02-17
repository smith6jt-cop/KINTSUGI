"""
Utility functions for REDSEA spillover correction.

Helper functions for marker list I/O, mask validation, and surface marker
identification used by the REDSEA wrapper and QC modules.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Default known surface/membrane markers that are candidates for spillover
# correction. Nuclear and intracellular markers are excluded because they
# do not spill laterally between cells.
KNOWN_SURFACE_MARKERS = {
    # T cell markers
    "CD3", "CD3e", "CD4", "CD8", "CD8a",
    # B cell markers
    "CD19", "CD20", "CD21",
    # Myeloid / macrophage markers
    "CD11b", "CD11c", "CD14", "CD15", "CD16", "CD33", "CD68", "CD163",
    # General leukocyte markers
    "CD45", "CD45RA", "CD45RO",
    # NK cell markers
    "CD56",
    # Other membrane markers
    "CD31", "CD34", "CD35", "CD44",
    "PanCK", "Pan-CK", "PanKeratin", "E-cadherin", "EpCAM",
    "HLA-DR", "PD-L1", "PD-1", "CTLA-4", "LAG-3", "TIM-3",
    "CD1c", "CD25", "CD27", "CD38", "CD57", "CD103", "CD138",
}

# Markers known to be nuclear/intracellular (should NOT be corrected)
KNOWN_NUCLEAR_MARKERS = {
    "DAPI", "Hoechst", "FoxP3", "FOXP3", "Ki67", "Ki-67",
    "Granzyme B", "GranzymeB", "GrB",
    "pSTAT3", "pSTAT5", "pERK", "pS6",
    "Blank", "Empty", "AF", "Autofluorescence",
}


def identify_surface_markers(
    marker_names: list[str],
    known_surface_markers: Optional[set[str]] = None,
) -> list[str]:
    """Identify which markers are surface/membrane markers for spillover correction.

    Surface markers are candidates for REDSEA correction because they reside on
    the cell membrane and can spill laterally into neighboring cells. Nuclear and
    intracellular markers generally do not exhibit lateral spillover.

    Parameters
    ----------
    marker_names : list of str
        Full list of marker names from the panel.
    known_surface_markers : set of str, optional
        Custom set of surface marker names. If None, uses the built-in set.

    Returns
    -------
    list of str
        Subset of marker_names that are surface/membrane markers.
    """
    if known_surface_markers is None:
        known_surface_markers = KNOWN_SURFACE_MARKERS

    surface = []
    for name in marker_names:
        # Skip blanks and known nuclear markers
        name_clean = name.strip()
        if name_clean in KNOWN_NUCLEAR_MARKERS:
            continue
        if name_clean.lower() in {"blank", "empty", "af", "autofluorescence"}:
            continue
        # Check if it matches a known surface marker (case-insensitive prefix match)
        if name_clean in known_surface_markers:
            surface.append(name_clean)
        elif any(name_clean.startswith(s) for s in known_surface_markers):
            surface.append(name_clean)
    return surface


def estimate_element_size(segmentation_mask: np.ndarray) -> int:
    """Estimate optimal REDSEA element size from mean cell area.

    Based on the empirical guidelines from Bai et al. (2021):
      - mean area ~107 pixels -> element_size = 2
      - mean area ~325 pixels -> element_size = 3-4

    Parameters
    ----------
    segmentation_mask : np.ndarray
        2D integer array where each cell has a unique label (0 = background).

    Returns
    -------
    int
        Recommended element_size parameter for REDSEA (typically 2-4).
    """
    from skimage.measure import regionprops

    # Get unique cell labels, excluding background (0)
    props = regionprops(segmentation_mask)
    if not props:
        logger.warning("No cells found in mask; defaulting to element_size=2")
        return 2

    areas = np.array([p.area for p in props])
    mean_area = np.mean(areas)
    median_area = np.median(areas)

    logger.info(
        f"Cell statistics: n={len(areas)}, mean_area={mean_area:.0f}, "
        f"median_area={median_area:.0f}"
    )

    # Empirical mapping from Bai et al.
    if median_area < 150:
        size = 2
    elif median_area < 400:
        size = 3
    else:
        size = 4

    logger.info(f"Recommended element_size={size} (median_area={median_area:.0f})")
    return size


def validate_segmentation_mask(mask: np.ndarray) -> dict:
    """Validate a segmentation mask for REDSEA compatibility.

    Parameters
    ----------
    mask : np.ndarray
        2D integer segmentation mask (0 = background, positive = cell labels).

    Returns
    -------
    dict
        Validation result with keys: 'valid' (bool), 'n_cells' (int),
        'issues' (list of str).
    """
    issues = []

    if mask.ndim != 2:
        issues.append(f"Expected 2D mask, got {mask.ndim}D")

    if not np.issubdtype(mask.dtype, np.integer):
        issues.append(f"Expected integer dtype, got {mask.dtype}")

    labels = np.unique(mask)
    n_cells = len(labels[labels > 0])

    if n_cells == 0:
        issues.append("No cells found (no labels > 0)")

    if np.any(mask < 0):
        issues.append("Mask contains negative labels")

    return {
        "valid": len(issues) == 0,
        "n_cells": n_cells,
        "issues": issues,
    }


def write_markers_csv(
    marker_names: list[str],
    output_path: Union[str, Path],
) -> Path:
    """Write marker names to a CSV file in the format expected by redseapy.

    Parameters
    ----------
    marker_names : list of str
        List of marker names.
    output_path : str or Path
        Path to write the CSV file.

    Returns
    -------
    Path
        Path to the written CSV file.
    """
    output_path = Path(output_path)
    df = pd.DataFrame({"marker_name": marker_names})
    df.to_csv(output_path, index=False)
    return output_path


def load_markers_csv(csv_path: Union[str, Path]) -> list[str]:
    """Load marker names from a CSV file.

    Parameters
    ----------
    csv_path : str or Path
        Path to CSV file with 'marker_name' column.

    Returns
    -------
    list of str
        List of marker names.
    """
    df = pd.read_csv(csv_path)
    if "marker_name" not in df.columns:
        raise ValueError(
            f"CSV must have 'marker_name' column. Found: {list(df.columns)}"
        )
    return df["marker_name"].tolist()
