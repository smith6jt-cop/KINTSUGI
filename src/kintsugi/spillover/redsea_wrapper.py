"""
REDSEA spillover correction wrapper for KINTSUGI.

Wraps the redseapy package (Bai et al., 2021) to provide spillover correction
for multiplexed immunofluorescence data. Supports both file-based and in-memory
array inputs.

Reference:
    Bai Y et al. Adjacent Cell Marker Lateral Spillover Compensation and
    Reinforcement for Multiplexed Images. Front. Immunol. 12:652631 (2021).
    PMID: 34295327
"""

from __future__ import annotations

import logging
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


@dataclass
class SpilloverResult:
    """Result of REDSEA spillover correction.

    Attributes
    ----------
    corrected : pd.DataFrame
        Per-cell quantification after spillover correction (cells x markers).
    uncorrected : pd.DataFrame
        Per-cell quantification before correction (cells x markers).
    parameters : dict
        Parameters used for correction.
    n_cells : int
        Number of cells processed.
    n_markers : int
        Number of markers corrected.
    """

    corrected: pd.DataFrame
    uncorrected: pd.DataFrame
    parameters: dict = field(default_factory=dict)
    n_cells: int = 0
    n_markers: int = 0


def _check_redseapy():
    """Check that redseapy is installed and importable."""
    try:
        from redseapy.redsea import run_redsea  # noqa: F401

        return True
    except ImportError:
        raise ImportError(
            "redseapy is required for spillover correction. "
            "Install with: pip install redseapy  "
            "or: pip install git+https://github.com/labsyspharm/redseapy.git"
        )


def run_redsea_correction(
    multichannel_image: Union[np.ndarray, str, Path],
    segmentation_mask: Union[np.ndarray, str, Path],
    marker_names: list[str],
    markers_of_interest: Optional[list[str]] = None,
    element_shape: str = "star",
    element_size: int = 2,
    output_dir: Optional[Union[str, Path]] = None,
) -> SpilloverResult:
    """Run REDSEA spillover correction on multiplexed image data.

    Wraps the redseapy ``run_redsea`` function to provide a more Pythonic
    interface with in-memory array support and structured results.

    Parameters
    ----------
    multichannel_image : np.ndarray or str or Path
        Multichannel image as a 3D array (Y, X, C) or path to a multichannel
        TIFF file.
    segmentation_mask : np.ndarray or str or Path
        2D integer segmentation mask (0 = background) or path to a TIFF file.
    marker_names : list of str
        Marker names corresponding to each channel.
    markers_of_interest : list of str, optional
        Subset of markers to correct. If None, all markers are corrected.
    element_shape : str
        Morphological element shape: "star" (diamond, recommended) or "square".
        Default: "star".
    element_size : int
        Extension size for the morphological element in pixels. Scale with cell
        size; 2-4 recommended. Default: 2.
    output_dir : str or Path, optional
        Directory to save output CSV files. If None, uses a temporary directory.

    Returns
    -------
    SpilloverResult
        Structured result with corrected/uncorrected DataFrames and metadata.

    Raises
    ------
    ImportError
        If redseapy is not installed.
    ValueError
        If inputs are invalid (shape mismatch, invalid parameters).
    """
    _check_redseapy()
    from redseapy.redsea import run_redsea

    # --- Validate element_shape ---
    shape_map = {"star": 2, "diamond": 2, "square": 1}
    if element_shape not in shape_map:
        raise ValueError(
            f"element_shape must be 'star' or 'square', got '{element_shape}'"
        )
    element_shape_int = shape_map[element_shape]

    # --- Handle in-memory arrays vs file paths ---
    use_temp = False
    temp_dir_obj = None

    try:
        if isinstance(multichannel_image, np.ndarray) or isinstance(
            segmentation_mask, np.ndarray
        ):
            # Need temp files for redseapy (it expects file paths)
            import tifffile

            temp_dir_obj = tempfile.TemporaryDirectory()
            temp_path = Path(temp_dir_obj.name)
            use_temp = True

            if isinstance(multichannel_image, np.ndarray):
                img_path = temp_path / "image.tif"
                # redseapy expects (Y, X, C) or OME-TIFF
                tifffile.imwrite(str(img_path), multichannel_image)
            else:
                img_path = Path(multichannel_image)

            if isinstance(segmentation_mask, np.ndarray):
                mask_path = temp_path / "mask.tif"
                tifffile.imwrite(str(mask_path), segmentation_mask.astype(np.int32))
            else:
                mask_path = Path(segmentation_mask)
        else:
            img_path = Path(multichannel_image)
            mask_path = Path(segmentation_mask)

        # --- Write markers CSV ---
        from .utils import write_markers_csv

        if use_temp:
            markers_csv_path = Path(temp_dir_obj.name) / "markers.csv"
        else:
            markers_csv_path = (
                Path(output_dir) / "markers.csv" if output_dir else Path(tempfile.mktemp(suffix=".csv"))
            )
        write_markers_csv(marker_names, markers_csv_path)

        # --- Set up output directory ---
        if output_dir is not None:
            out_dir = Path(output_dir)
            out_dir.mkdir(parents=True, exist_ok=True)
        elif use_temp:
            out_dir = Path(temp_dir_obj.name) / "output"
            out_dir.mkdir(parents=True, exist_ok=True)
        else:
            temp_dir_obj = tempfile.TemporaryDirectory()
            use_temp = True
            out_dir = Path(temp_dir_obj.name) / "output"
            out_dir.mkdir(parents=True, exist_ok=True)

        # --- Run REDSEA ---
        logger.info(
            f"Running REDSEA: element_shape={element_shape} ({element_shape_int}), "
            f"element_size={element_size}, markers_of_interest={markers_of_interest}"
        )

        run_redsea(
            tiff=str(img_path),
            seg_mask=str(mask_path),
            markers_csv=str(markers_csv_path),
            output_dir=str(out_dir),
            markers_of_interest=markers_of_interest,
            element_shape=element_shape_int,
            element_size=element_size,
        )

        # --- Load results ---
        corrected_path = out_dir / "single_cell_after_redsea.csv"
        uncorrected_path = out_dir / "single_cell_before_redsea.csv"

        if not corrected_path.exists():
            raise RuntimeError(
                f"REDSEA did not produce expected output: {corrected_path}"
            )

        corrected_df = pd.read_csv(corrected_path)
        uncorrected_df = pd.read_csv(uncorrected_path)

        # Build parameter record
        parameters = {
            "element_shape": element_shape,
            "element_size": element_size,
            "markers_of_interest": markers_of_interest,
            "n_markers_total": len(marker_names),
            "marker_names": marker_names,
        }

        result = SpilloverResult(
            corrected=corrected_df,
            uncorrected=uncorrected_df,
            parameters=parameters,
            n_cells=len(corrected_df),
            n_markers=len(corrected_df.columns),
        )

        logger.info(
            f"REDSEA complete: {result.n_cells} cells x {result.n_markers} markers"
        )
        return result

    finally:
        # Clean up temp directory if we created one (and user didn't specify output_dir)
        if use_temp and temp_dir_obj is not None and output_dir is None:
            try:
                temp_dir_obj.cleanup()
            except Exception:
                pass
