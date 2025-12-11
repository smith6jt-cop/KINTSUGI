"""
Batch-Level Quality Control

Detection and assessment of batch effects in multiplex imaging experiments.
Inspired by CyLinter's batch effect detection approaches.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Any

import numpy as np
from scipy import stats

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger("kintsugi.qc.batch_qc")


@dataclass
class BatchQCResult:
    """Results from batch-level quality control."""

    # Overall assessment
    batch_effects_detected: bool
    severity: str  # "none", "mild", "moderate", "severe"

    # Batch statistics
    batch_stats: dict[str, dict[str, float]] = field(default_factory=dict)

    # Statistical tests
    statistical_tests: dict[str, dict[str, float]] = field(default_factory=dict)

    # Affected markers
    affected_markers: list[str] = field(default_factory=list)

    # Recommendations
    recommendations: list[str] = field(default_factory=list)


class BatchQC:
    """
    Batch effect detection and assessment for multiplex imaging.

    Detects:
    - Intensity shifts between batches
    - Distribution differences
    - Cycle-to-cycle variations
    - Technical artifacts affecting batches

    Example
    -------
    >>> qc = BatchQC()
    >>> result = qc.assess(
    ...     data=cell_data,
    ...     batch_column="batch_id",
    ...     marker_columns=["CD3", "CD20", "DAPI"],
    ... )
    >>> if result.batch_effects_detected:
    ...     print(f"Affected markers: {result.affected_markers}")
    """

    def __init__(
        self,
        p_value_threshold: float = 0.01,
        effect_size_threshold: float = 0.3,
        max_cv_difference: float = 0.5,
    ):
        """
        Initialize Batch QC.

        Parameters
        ----------
        p_value_threshold : float
            P-value threshold for statistical significance
        effect_size_threshold : float
            Effect size (Cohen's d) threshold for practical significance
        max_cv_difference : float
            Maximum coefficient of variation difference between batches
        """
        self.p_value_threshold = p_value_threshold
        self.effect_size_threshold = effect_size_threshold
        self.max_cv_difference = max_cv_difference

    def assess(
        self,
        data: np.ndarray | pd.DataFrame,
        batch_column: str | int,
        marker_columns: list[str | int],
    ) -> BatchQCResult:
        """
        Assess batch effects in the data.

        Parameters
        ----------
        data : array-like or DataFrame
            Cell/observation data
        batch_column : str or int
            Column name/index for batch identifier
        marker_columns : list
            Column names/indices for markers to assess

        Returns
        -------
        BatchQCResult
            Batch effect assessment results
        """
        # Convert to numpy arrays
        if hasattr(data, "values"):
            df = data
            batch_col_idx = (
                df.columns.get_loc(batch_column) if isinstance(batch_column, str) else batch_column
            )
            marker_col_indices = [
                df.columns.get_loc(c) if isinstance(c, str) else c for c in marker_columns
            ]
            marker_names = [str(c) for c in marker_columns]
            values = df.values
        else:
            df = None
            batch_col_idx = batch_column
            marker_col_indices = marker_columns
            marker_names = [f"marker_{i}" for i in marker_columns]
            values = np.asarray(data)

        # Extract batch labels and marker data
        batch_labels = values[:, batch_col_idx]
        unique_batches = np.unique(batch_labels)

        if len(unique_batches) < 2:
            return BatchQCResult(
                batch_effects_detected=False,
                severity="none",
                recommendations=["Only one batch detected, cannot assess batch effects"],
            )

        batch_stats = {}
        statistical_tests = {}
        affected_markers = []
        recommendations = []

        # Analyze each marker
        for marker_idx, marker_name in zip(marker_col_indices, marker_names, strict=False):
            marker_values = values[:, marker_idx].astype(float)

            # Compute per-batch statistics
            batch_data = {}
            for batch in unique_batches:
                batch_mask = batch_labels == batch
                batch_vals = marker_values[batch_mask]
                batch_data[str(batch)] = {
                    "n": len(batch_vals),
                    "mean": float(np.mean(batch_vals)),
                    "std": float(np.std(batch_vals)),
                    "median": float(np.median(batch_vals)),
                    "cv": float(np.std(batch_vals) / (np.mean(batch_vals) + 1e-8)),
                    "p5": float(np.percentile(batch_vals, 5)),
                    "p95": float(np.percentile(batch_vals, 95)),
                }

            batch_stats[marker_name] = batch_data

            # Statistical tests for batch effects
            marker_tests = self._test_batch_effects(marker_values, batch_labels, unique_batches)
            statistical_tests[marker_name] = marker_tests

            # Check if marker is affected
            if marker_tests["batch_effect_detected"]:
                affected_markers.append(marker_name)

        # Overall assessment
        n_affected = len(affected_markers)
        n_markers = len(marker_names)
        affected_ratio = n_affected / n_markers

        if affected_ratio == 0:
            severity = "none"
            batch_effects_detected = False
        elif affected_ratio < 0.25:
            severity = "mild"
            batch_effects_detected = True
            recommendations.append(
                f"Mild batch effects in {n_affected}/{n_markers} markers. "
                "Consider marker-specific normalization."
            )
        elif affected_ratio < 0.5:
            severity = "moderate"
            batch_effects_detected = True
            recommendations.append(
                f"Moderate batch effects in {n_affected}/{n_markers} markers. "
                "Recommend batch correction before downstream analysis."
            )
        else:
            severity = "severe"
            batch_effects_detected = True
            recommendations.append(
                f"Severe batch effects in {n_affected}/{n_markers} markers. "
                "Critical: Apply batch correction or investigate technical issues."
            )

        return BatchQCResult(
            batch_effects_detected=batch_effects_detected,
            severity=severity,
            batch_stats=batch_stats,
            statistical_tests=statistical_tests,
            affected_markers=affected_markers,
            recommendations=recommendations,
        )

    def _test_batch_effects(
        self,
        values: np.ndarray,
        batch_labels: np.ndarray,
        unique_batches: np.ndarray,
    ) -> dict[str, Any]:
        """Run statistical tests for batch effects on a single marker."""
        results = {}

        # Group data by batch
        batch_groups = [values[batch_labels == b] for b in unique_batches]

        # 1. Kruskal-Wallis test (non-parametric ANOVA)
        try:
            h_stat, p_value = stats.kruskal(*batch_groups)
            results["kruskal_wallis_h"] = float(h_stat)
            results["kruskal_wallis_p"] = float(p_value)
        except Exception:
            results["kruskal_wallis_h"] = 0.0
            results["kruskal_wallis_p"] = 1.0

        # 2. Effect size (eta-squared for Kruskal-Wallis)
        n_total = len(values)
        k = len(unique_batches)
        if n_total > k:
            eta_squared = (h_stat - k + 1) / (n_total - k)
            results["eta_squared"] = float(np.clip(eta_squared, 0, 1))
        else:
            results["eta_squared"] = 0.0

        # 3. Pairwise Cohen's d (effect sizes)
        max_cohens_d = 0.0
        for i in range(len(batch_groups)):
            for j in range(i + 1, len(batch_groups)):
                d = self._cohens_d(batch_groups[i], batch_groups[j])
                max_cohens_d = max(max_cohens_d, abs(d))

        results["max_cohens_d"] = float(max_cohens_d)

        # 4. CV difference
        cvs = [np.std(g) / (np.mean(g) + 1e-8) for g in batch_groups]
        cv_range = max(cvs) - min(cvs)
        results["cv_range"] = float(cv_range)

        # 5. Levene's test for variance homogeneity
        try:
            levene_stat, levene_p = stats.levene(*batch_groups)
            results["levene_stat"] = float(levene_stat)
            results["levene_p"] = float(levene_p)
        except Exception:
            results["levene_stat"] = 0.0
            results["levene_p"] = 1.0

        # Overall batch effect detection
        significant_p = results["kruskal_wallis_p"] < self.p_value_threshold
        large_effect = results["max_cohens_d"] > self.effect_size_threshold
        variance_diff = results["levene_p"] < self.p_value_threshold

        results["batch_effect_detected"] = significant_p and (large_effect or variance_diff)

        return results

    def _cohens_d(self, group1: np.ndarray, group2: np.ndarray) -> float:
        """Calculate Cohen's d effect size."""
        n1, n2 = len(group1), len(group2)
        var1, var2 = np.var(group1, ddof=1), np.var(group2, ddof=1)

        # Pooled standard deviation
        pooled_std = np.sqrt(((n1 - 1) * var1 + (n2 - 1) * var2) / (n1 + n2 - 2))

        if pooled_std < 1e-8:
            return 0.0

        return float((np.mean(group1) - np.mean(group2)) / pooled_std)


def detect_batch_effects(
    data: np.ndarray | pd.DataFrame,
    batch_column: str | int,
    marker_columns: list[str | int],
) -> BatchQCResult:
    """
    Convenience function for batch effect detection.

    Parameters
    ----------
    data : array-like or DataFrame
        Cell/observation data
    batch_column : str or int
        Column name/index for batch identifier
    marker_columns : list
        Column names/indices for markers to assess

    Returns
    -------
    BatchQCResult
        Batch effect assessment results
    """
    qc = BatchQC()
    return qc.assess(data, batch_column, marker_columns)


def compute_batch_statistics(
    data: np.ndarray | pd.DataFrame,
    batch_column: str | int,
    value_column: str | int,
) -> dict[str, dict[str, float]]:
    """
    Compute per-batch statistics for a single variable.

    Parameters
    ----------
    data : array-like or DataFrame
        Data
    batch_column : str or int
        Batch identifier column
    value_column : str or int
        Value column to summarize

    Returns
    -------
    dict
        Statistics per batch
    """
    if hasattr(data, "values"):
        df = data
        batch_col_idx = (
            df.columns.get_loc(batch_column) if isinstance(batch_column, str) else batch_column
        )
        value_col_idx = (
            df.columns.get_loc(value_column) if isinstance(value_column, str) else value_column
        )
        values = df.values
    else:
        batch_col_idx = batch_column
        value_col_idx = value_column
        values = np.asarray(data)

    batch_labels = values[:, batch_col_idx]
    data_values = values[:, value_col_idx].astype(float)

    unique_batches = np.unique(batch_labels)
    stats_dict = {}

    for batch in unique_batches:
        batch_mask = batch_labels == batch
        batch_vals = data_values[batch_mask]

        stats_dict[str(batch)] = {
            "n": len(batch_vals),
            "mean": float(np.mean(batch_vals)),
            "std": float(np.std(batch_vals)),
            "median": float(np.median(batch_vals)),
            "min": float(np.min(batch_vals)),
            "max": float(np.max(batch_vals)),
            "p25": float(np.percentile(batch_vals, 25)),
            "p75": float(np.percentile(batch_vals, 75)),
            "cv": float(np.std(batch_vals) / (np.mean(batch_vals) + 1e-8)),
        }

    return stats_dict


def normalize_batches(
    data: np.ndarray,
    batch_labels: np.ndarray,
    method: str = "quantile",
    reference_batch: str | None = None,
) -> np.ndarray:
    """
    Normalize data to reduce batch effects.

    Parameters
    ----------
    data : np.ndarray
        Data to normalize (N x M features)
    batch_labels : np.ndarray
        Batch labels (N,)
    method : str
        Normalization method:
        - "zscore": Z-score normalization per batch
        - "quantile": Quantile normalization to reference
        - "median": Median centering
    reference_batch : str, optional
        Batch to use as reference (for quantile method)

    Returns
    -------
    np.ndarray
        Normalized data
    """
    unique_batches = np.unique(batch_labels)
    normalized = np.zeros_like(data, dtype=float)

    if method == "zscore":
        # Z-score normalize each batch independently
        for batch in unique_batches:
            mask = batch_labels == batch
            batch_data = data[mask].astype(float)

            mean = np.mean(batch_data, axis=0)
            std = np.std(batch_data, axis=0)
            std[std < 1e-8] = 1.0

            normalized[mask] = (batch_data - mean) / std

    elif method == "median":
        # Center each batch to global median
        global_median = np.median(data, axis=0)

        for batch in unique_batches:
            mask = batch_labels == batch
            batch_data = data[mask].astype(float)
            batch_median = np.median(batch_data, axis=0)

            normalized[mask] = batch_data - batch_median + global_median

    elif method == "quantile":
        # Quantile normalization to reference batch
        if reference_batch is None:
            reference_batch = unique_batches[0]

        ref_mask = batch_labels == reference_batch
        ref_data = data[ref_mask].astype(float)

        # Compute reference quantiles
        ref_sorted = np.sort(ref_data, axis=0)

        for batch in unique_batches:
            mask = batch_labels == batch
            batch_data = data[mask].astype(float)
            n = batch_data.shape[0]

            # Rank transform to reference distribution
            for col in range(batch_data.shape[1]):
                ranks = stats.rankdata(batch_data[:, col], method="ordinal") - 1

                # Map ranks to reference quantiles
                ref_n = len(ref_sorted[:, col])
                ref_indices = (ranks * (ref_n - 1) / (n - 1)).astype(int)
                ref_indices = np.clip(ref_indices, 0, ref_n - 1)

                normalized[mask, col] = ref_sorted[ref_indices, col]

    else:
        raise ValueError(f"Unknown method: {method}")

    return normalized
