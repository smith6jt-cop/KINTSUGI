"""
Quality control metrics and visualization for REDSEA spillover correction.

Provides functions to compare pre/post correction quantification, compute
double-positive reduction metrics for mutually exclusive marker pairs,
and generate biaxial scatter plots and summary heatmaps.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

# Common mutually exclusive marker pairs for QC validation.
# These pairs should NOT co-express on the same cell; double-positives
# indicate spillover contamination.
DEFAULT_EXCLUSIVE_PAIRS = [
    ("CD3", "CD20"),      # T cell / B cell
    ("CD3e", "CD20"),     # T cell / B cell (alternative name)
    ("CD4", "CD8"),       # Helper T / Cytotoxic T
    ("CD4", "CD8a"),      # Helper T / Cytotoxic T (alternative name)
    ("CD3", "PanCK"),     # T cell / Epithelial
    ("CD3e", "PanCK"),    # T cell / Epithelial
    ("CD20", "CD68"),     # B cell / Macrophage
    ("CD3", "CD68"),      # T cell / Macrophage
    ("CD3e", "CD68"),     # T cell / Macrophage
]


def _find_column(df: pd.DataFrame, marker: str) -> Optional[str]:
    """Find column matching a marker name (handles suffixes like _mean)."""
    if marker in df.columns:
        return marker
    # Try common suffixes
    for suffix in ["_mean", "_median", "_sum"]:
        candidate = marker + suffix
        if candidate in df.columns:
            return candidate
    return None


def _get_threshold(series: pd.Series, method: str = "otsu") -> float:
    """Compute a positivity threshold for a marker.

    Parameters
    ----------
    series : pd.Series
        Marker intensity values.
    method : str
        Thresholding method: "otsu" or "percentile".

    Returns
    -------
    float
        Threshold value above which cells are considered positive.
    """
    values = series.dropna().values
    if len(values) == 0:
        return 0.0

    if method == "otsu":
        try:
            from skimage.filters import threshold_otsu

            return float(threshold_otsu(values[values > 0]) if np.any(values > 0) else 0.0)
        except (ValueError, ImportError):
            # Fallback to percentile if Otsu fails
            return float(np.percentile(values[values > 0], 75)) if np.any(values > 0) else 0.0
    else:
        return float(np.percentile(values[values > 0], 75)) if np.any(values > 0) else 0.0


def compute_spillover_metrics(
    corrected: pd.DataFrame,
    uncorrected: pd.DataFrame,
    mutually_exclusive_pairs: Optional[list[tuple[str, str]]] = None,
    threshold_method: str = "otsu",
) -> dict:
    """Compute QC metrics comparing pre/post spillover correction.

    For each pair of mutually exclusive markers, computes the double-positive
    rate before and after correction. A successful correction should reduce
    double-positive rates while preserving single-positive populations.

    Parameters
    ----------
    corrected : pd.DataFrame
        Per-cell quantification after REDSEA correction.
    uncorrected : pd.DataFrame
        Per-cell quantification before correction.
    mutually_exclusive_pairs : list of (str, str), optional
        Marker pairs that should not co-express. If None, uses defaults
        filtered to markers present in the data.
    threshold_method : str
        Method for computing positivity thresholds: "otsu" or "percentile".

    Returns
    -------
    dict
        Dictionary with keys:
        - 'pair_metrics': dict mapping (marker_a, marker_b) to per-pair metrics
        - 'summary': dict with overall statistics
        - 'per_marker_change': dict mapping marker name to mean intensity change
    """
    if mutually_exclusive_pairs is None:
        mutually_exclusive_pairs = DEFAULT_EXCLUSIVE_PAIRS

    pair_metrics = {}
    pairs_evaluated = 0

    for marker_a, marker_b in mutually_exclusive_pairs:
        col_a_corr = _find_column(corrected, marker_a)
        col_b_corr = _find_column(corrected, marker_b)
        col_a_uncorr = _find_column(uncorrected, marker_a)
        col_b_uncorr = _find_column(uncorrected, marker_b)

        if not all([col_a_corr, col_b_corr, col_a_uncorr, col_b_uncorr]):
            continue

        # Compute thresholds from uncorrected data
        thresh_a = _get_threshold(uncorrected[col_a_uncorr], threshold_method)
        thresh_b = _get_threshold(uncorrected[col_b_uncorr], threshold_method)

        # Double-positive rates
        uncorr_dp = (
            (uncorrected[col_a_uncorr] > thresh_a)
            & (uncorrected[col_b_uncorr] > thresh_b)
        ).sum()
        corr_dp = (
            (corrected[col_a_corr] > thresh_a)
            & (corrected[col_b_corr] > thresh_b)
        ).sum()

        n_cells = len(uncorrected)
        uncorr_dp_pct = 100 * uncorr_dp / n_cells if n_cells > 0 else 0
        corr_dp_pct = 100 * corr_dp / n_cells if n_cells > 0 else 0

        # Single-positive preservation
        uncorr_sp_a = (
            (uncorrected[col_a_uncorr] > thresh_a)
            & (uncorrected[col_b_uncorr] <= thresh_b)
        ).sum()
        corr_sp_a = (
            (corrected[col_a_corr] > thresh_a)
            & (corrected[col_b_corr] <= thresh_b)
        ).sum()

        uncorr_sp_b = (
            (uncorrected[col_a_uncorr] <= thresh_a)
            & (uncorrected[col_b_uncorr] > thresh_b)
        ).sum()
        corr_sp_b = (
            (corrected[col_a_corr] <= thresh_a)
            & (corrected[col_b_corr] > thresh_b)
        ).sum()

        pair_metrics[(marker_a, marker_b)] = {
            "uncorrected_dp_count": int(uncorr_dp),
            "corrected_dp_count": int(corr_dp),
            "uncorrected_dp_pct": float(uncorr_dp_pct),
            "corrected_dp_pct": float(corr_dp_pct),
            "dp_reduction_pct": float(uncorr_dp_pct - corr_dp_pct),
            "threshold_a": float(thresh_a),
            "threshold_b": float(thresh_b),
            "uncorrected_sp_a": int(uncorr_sp_a),
            "corrected_sp_a": int(corr_sp_a),
            "uncorrected_sp_b": int(uncorr_sp_b),
            "corrected_sp_b": int(corr_sp_b),
        }
        pairs_evaluated += 1

    # Per-marker mean intensity change
    per_marker_change = {}
    common_cols = set(corrected.columns) & set(uncorrected.columns)
    for col in common_cols:
        if corrected[col].dtype.kind in "iufb":
            mean_before = uncorrected[col].mean()
            mean_after = corrected[col].mean()
            if mean_before > 0:
                per_marker_change[col] = {
                    "mean_before": float(mean_before),
                    "mean_after": float(mean_after),
                    "fold_change": float(mean_after / mean_before),
                    "pct_change": float(100 * (mean_after - mean_before) / mean_before),
                }

    # Summary
    total_dp_reduction = sum(
        m["dp_reduction_pct"] for m in pair_metrics.values()
    )
    avg_dp_reduction = total_dp_reduction / pairs_evaluated if pairs_evaluated > 0 else 0

    return {
        "pair_metrics": pair_metrics,
        "per_marker_change": per_marker_change,
        "summary": {
            "pairs_evaluated": pairs_evaluated,
            "avg_dp_reduction_pct": float(avg_dp_reduction),
            "n_cells": len(corrected),
        },
    }


def plot_spillover_comparison(
    corrected: pd.DataFrame,
    uncorrected: pd.DataFrame,
    marker_pair: tuple[str, str],
    threshold_method: str = "otsu",
    save_path: Optional[Union[str, Path]] = None,
):
    """Generate biaxial scatter plots comparing pre/post correction for a marker pair.

    Creates a side-by-side figure with scatter plots showing marker co-expression
    before (left) and after (right) REDSEA correction. Quadrant gates show
    double-positive populations.

    Parameters
    ----------
    corrected : pd.DataFrame
        Per-cell quantification after correction.
    uncorrected : pd.DataFrame
        Per-cell quantification before correction.
    marker_pair : tuple of (str, str)
        Pair of marker names to plot.
    threshold_method : str
        Thresholding method for quadrant gates.
    save_path : str or Path, optional
        Path to save the figure. If None, returns figure without saving.

    Returns
    -------
    matplotlib.figure.Figure
        The generated figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    marker_a, marker_b = marker_pair

    col_a_corr = _find_column(corrected, marker_a)
    col_b_corr = _find_column(corrected, marker_b)
    col_a_uncorr = _find_column(uncorrected, marker_a)
    col_b_uncorr = _find_column(uncorrected, marker_b)

    if not all([col_a_corr, col_b_corr, col_a_uncorr, col_b_uncorr]):
        raise ValueError(
            f"Markers {marker_a}/{marker_b} not found in both DataFrames"
        )

    # Thresholds from uncorrected data
    thresh_a = _get_threshold(uncorrected[col_a_uncorr], threshold_method)
    thresh_b = _get_threshold(uncorrected[col_b_uncorr], threshold_method)

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    for ax, df, col_a, col_b, title in [
        (axes[0], uncorrected, col_a_uncorr, col_b_uncorr, "Before REDSEA"),
        (axes[1], corrected, col_a_corr, col_b_corr, "After REDSEA"),
    ]:
        x = df[col_a].values
        y = df[col_b].values

        # Scatter
        ax.scatter(x, y, s=1, alpha=0.3, c="steelblue", rasterized=True)

        # Quadrant lines
        ax.axvline(thresh_a, color="red", linestyle="--", linewidth=0.8, alpha=0.7)
        ax.axhline(thresh_b, color="red", linestyle="--", linewidth=0.8, alpha=0.7)

        # Count quadrants
        n_total = len(df)
        dp = ((x > thresh_a) & (y > thresh_b)).sum()
        dp_pct = 100 * dp / n_total if n_total > 0 else 0

        ax.set_xlabel(marker_a)
        ax.set_ylabel(marker_b)
        ax.set_title(f"{title}\nDP: {dp} ({dp_pct:.1f}%)")

        # Annotate quadrant
        ax.text(
            0.95, 0.95, f"{dp_pct:.1f}%",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=12, fontweight="bold", color="red",
        )

    fig.suptitle(f"Spillover Correction: {marker_a} vs {marker_b}", fontsize=14)
    plt.tight_layout()

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved spillover QC plot: {save_path}")

    return fig


def plot_spillover_summary_heatmap(
    corrected: pd.DataFrame,
    uncorrected: pd.DataFrame,
    marker_names: list[str],
    save_path: Optional[Union[str, Path]] = None,
):
    """Generate a heatmap of fold-change in co-expression for all marker pairs.

    Parameters
    ----------
    corrected : pd.DataFrame
        Per-cell quantification after correction.
    uncorrected : pd.DataFrame
        Per-cell quantification before correction.
    marker_names : list of str
        Marker names to include in the heatmap.
    save_path : str or Path, optional
        Path to save the figure.

    Returns
    -------
    matplotlib.figure.Figure
        The generated heatmap figure.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Filter to markers present in both dataframes
    available = []
    col_map_corr = {}
    col_map_uncorr = {}
    for m in marker_names:
        c_corr = _find_column(corrected, m)
        c_uncorr = _find_column(uncorrected, m)
        if c_corr and c_uncorr:
            available.append(m)
            col_map_corr[m] = c_corr
            col_map_uncorr[m] = c_uncorr

    if len(available) < 2:
        logger.warning("Need at least 2 markers for heatmap")
        fig, ax = plt.subplots()
        ax.text(0.5, 0.5, "Insufficient markers", ha="center", va="center")
        return fig

    n = len(available)
    fc_matrix = np.ones((n, n))

    for i, mi in enumerate(available):
        for j, mj in enumerate(available):
            if i == j:
                continue

            # Mean co-expression (product of means as proxy)
            uncorr_mean_i = uncorrected[col_map_uncorr[mi]].mean()
            uncorr_mean_j = uncorrected[col_map_uncorr[mj]].mean()
            corr_mean_i = corrected[col_map_corr[mi]].mean()
            corr_mean_j = corrected[col_map_corr[mj]].mean()

            if uncorr_mean_i > 0 and uncorr_mean_j > 0:
                uncorr_coexp = uncorr_mean_i * uncorr_mean_j
                corr_coexp = corr_mean_i * corr_mean_j
                fc_matrix[i, j] = corr_coexp / uncorr_coexp if uncorr_coexp > 0 else 1.0

    fig, ax = plt.subplots(figsize=(max(8, n * 0.6), max(6, n * 0.5)))
    im = ax.imshow(fc_matrix, cmap="RdBu_r", vmin=0.5, vmax=1.5, aspect="auto")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(available, rotation=45, ha="right", fontsize=8)
    ax.set_yticklabels(available, fontsize=8)

    plt.colorbar(im, ax=ax, label="Fold change (corrected / uncorrected)")
    ax.set_title("REDSEA Correction: Co-expression Fold Change")
    plt.tight_layout()

    if save_path:
        save_path = Path(save_path)
        save_path.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info(f"Saved spillover summary heatmap: {save_path}")

    return fig
