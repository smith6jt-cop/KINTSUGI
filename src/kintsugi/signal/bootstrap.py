"""
Bootstrap Learning Engine for Weighted Autofluorescence Subtraction

Pre-populates the parameter learning database from existing before/after pairs
on disk. This enables the learning system to make recommendations for new
channels without requiring interactive MCP sessions.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

logger = logging.getLogger("kintsugi.signal.bootstrap")


def bootstrap_from_project(
    project_dir: str | Path,
    tissue_type: str = "unknown",
    dry_run: bool = False,
    verbose: bool = False,
) -> dict[str, Any]:
    """
    Scan a KINTSUGI project for processable signal/blank pairs and bootstrap
    the learning database with weighted subtraction parameters.

    Looks for:
    - Signal: data/processed/edf/cyc##/{marker}.tif
    - Blank: any channel matching "Blank*" in the same cycle
    - Previously approved results: subtraction_summary.json

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory
    tissue_type : str
        Tissue type to record. Default: "unknown"
    dry_run : bool
        If True, only report what would be done. Default: False
    verbose : bool
        If True, print detailed progress. Default: False

    Returns
    -------
    dict
        Summary with counts of pairs found, processed, and recorded.
    """
    import tifffile

    project_dir = Path(project_dir)

    # Find EDF outputs
    edf_dir = project_dir / "data" / "processed" / "edf"
    if not edf_dir.exists():
        return {
            "status": "no_data",
            "message": f"No EDF directory found at {edf_dir}",
            "pairs_found": 0,
        }

    # Find all cycle directories
    cycle_dirs = sorted([d for d in edf_dir.iterdir() if d.is_dir()])
    if not cycle_dirs:
        return {
            "status": "no_data",
            "message": "No cycle directories found in EDF output",
            "pairs_found": 0,
        }

    # Scan for signal/blank pairs
    pairs: list[dict[str, Any]] = []

    for cycle_dir in cycle_dirs:
        cycle_name = cycle_dir.name
        tiff_files = list(cycle_dir.glob("*.tif")) + list(cycle_dir.glob("*.tiff"))

        # Identify blank channels
        blank_files = [f for f in tiff_files if "blank" in f.stem.lower()]
        signal_files = [f for f in tiff_files if "blank" not in f.stem.lower()]

        if not blank_files:
            if verbose:
                logger.info(f"  {cycle_name}: no blank channel found, skipping")
            continue

        # Use first blank as reference
        blank_path = blank_files[0]

        for signal_path in signal_files:
            marker_name = signal_path.stem
            # Skip DAPI and other nuclear stains
            if marker_name.upper().startswith("DAPI"):
                continue

            pairs.append(
                {
                    "cycle": cycle_name,
                    "marker": marker_name,
                    "signal_path": signal_path,
                    "blank_path": blank_path,
                }
            )

    if verbose:
        logger.info(f"Found {len(pairs)} signal/blank pairs across {len(cycle_dirs)} cycles")

    if dry_run:
        return {
            "status": "dry_run",
            "pairs_found": len(pairs),
            "cycles": len(cycle_dirs),
            "pairs": [{"cycle": p["cycle"], "marker": p["marker"]} for p in pairs],
        }

    # Load pairs and compute weighted parameters
    results = []
    for pair in pairs:
        try:
            signal = tifffile.imread(str(pair["signal_path"]))
            blank = tifffile.imread(str(pair["blank_path"]))

            if signal.shape != blank.shape:
                logger.warning(
                    f"Shape mismatch for {pair['marker']}: "
                    f"signal={signal.shape}, blank={blank.shape}"
                )
                continue

            result = _analyze_and_record_pair(
                signal=signal,
                blank=blank,
                marker_name=pair["marker"],
                tissue_type=tissue_type,
                project_dir=project_dir,
            )
            results.append(
                {
                    "cycle": pair["cycle"],
                    "marker": pair["marker"],
                    "recorded": result.get("recorded", False),
                    "quality": result.get("quality_score", 0),
                }
            )

            if verbose:
                logger.info(
                    f"  {pair['cycle']}/{pair['marker']}: "
                    f"quality={result.get('quality_score', 0):.3f}"
                )

        except Exception as e:
            logger.warning(f"Failed to process {pair['cycle']}/{pair['marker']}: {e}")
            results.append(
                {
                    "cycle": pair["cycle"],
                    "marker": pair["marker"],
                    "recorded": False,
                    "error": str(e),
                }
            )

    recorded = sum(1 for r in results if r.get("recorded"))
    return {
        "status": "success",
        "pairs_found": len(pairs),
        "pairs_processed": len(results),
        "pairs_recorded": recorded,
        "results": results,
    }


def bootstrap_from_pairs(
    pairs: list[tuple[np.ndarray, np.ndarray, str]],
    tissue_type: str = "unknown",
    project_dir: str | Path | None = None,
) -> list[dict[str, Any]]:
    """
    Compute weighted parameters for explicit (signal, blank, marker) tuples
    and record to the learning database.

    Parameters
    ----------
    pairs : list of (signal, blank, marker_name) tuples
    tissue_type : str
        Tissue type to record
    project_dir : str or Path, optional
        Project directory for learning database

    Returns
    -------
    list of dict
        Results for each pair with quality scores and recorded status
    """
    results = []
    for signal, blank, marker_name in pairs:
        try:
            result = _analyze_and_record_pair(
                signal=signal,
                blank=blank,
                marker_name=marker_name,
                tissue_type=tissue_type,
                project_dir=Path(project_dir) if project_dir else None,
            )
            results.append(
                {
                    "marker": marker_name,
                    "recorded": result.get("recorded", False),
                    "quality": result.get("quality_score", 0),
                }
            )
        except Exception as e:
            results.append(
                {
                    "marker": marker_name,
                    "recorded": False,
                    "error": str(e),
                }
            )
    return results


def _analyze_and_record_pair(
    signal: np.ndarray,
    blank: np.ndarray,
    marker_name: str,
    tissue_type: str,
    project_dir: Path | None = None,
) -> dict[str, Any]:
    """Analyze a signal/blank pair, compute weighted params, record to DB."""
    from .autofluorescence import (
        analyze_for_weighted_subtraction,
        compute_weighted_subtraction_quality,
        subtract_autofluorescence_weighted,
    )

    # Analyze
    analysis = analyze_for_weighted_subtraction(
        signal, blank, tissue_type=tissue_type, marker_name=marker_name
    )

    # Perform subtraction
    subtracted = subtract_autofluorescence_weighted(
        signal,
        blank,
        blank_clip_factor=analysis["blank_clip_factor"],
        base_scale_factor=analysis["base_scale_factor"],
        n_ranges=analysis["n_ranges"],
        range_method=analysis["range_method"],
        transition_width=analysis.get("transition_width", 0.1),
    )

    # Compute quality
    quality = compute_weighted_subtraction_quality(signal, subtracted, blank, analysis["ranges"])
    quality_score = quality.get("global", {}).get("quality_score", 0)

    # Build parameter dict for recording
    params_dict = {
        "algorithm_version": "weighted_v1",
        "blank_clip_factor": analysis["blank_clip_factor"],
        "base_scale_factor": analysis["base_scale_factor"],
        "n_ranges": analysis["n_ranges"],
        "range_method": analysis["range_method"],
        "transition_width": analysis.get("transition_width", 0.1),
        "ranges": analysis["ranges"],
    }

    # Record to learning database
    recorded = False
    try:
        from kintsugi.claude.parameter_learning import ParameterLearningEngine

        engine = ParameterLearningEngine(str(project_dir) if project_dir else None)
        engine.record_parameters(
            tissue_type=tissue_type,
            marker_name=marker_name,
            operation="blank_subtraction",
            parameters=params_dict,
            quality_after=quality_score,
            user_approved=False,  # Bootstrap — not human-verified
            user_notes="bootstrap: auto-computed from existing data",
            algorithm_version="weighted_v1",
        )
        recorded = True
    except Exception as e:
        logger.warning(f"Failed to record to learning DB: {e}")

    return {
        "recorded": recorded,
        "quality_score": quality_score,
        "parameters": params_dict,
    }
