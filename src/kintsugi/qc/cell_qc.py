"""
Cell-Level Quality Control

Quality assessment at the single-cell level for multiplex imaging.
Inspired by CyLinter's cell filtering approaches.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Union

import numpy as np
from scipy.spatial.distance import cdist

if TYPE_CHECKING:
    import pandas as pd

logger = logging.getLogger("kintsugi.qc.cell_qc")

# Type alias
ArrayLike = Union[np.ndarray, "pd.DataFrame"]


@dataclass
class CellQCResult:
    """Results from cell-level quality control."""

    # Counts
    total_cells: int
    passed_cells: int
    filtered_cells: int

    # Filter breakdown
    filter_counts: dict[str, int] = field(default_factory=dict)

    # Indices
    passed_indices: np.ndarray = field(default_factory=lambda: np.array([]))
    filtered_indices: np.ndarray = field(default_factory=lambda: np.array([]))

    # Reasons for filtering
    filter_reasons: dict[int, list[str]] = field(default_factory=dict)

    @property
    def pass_rate(self) -> float:
        """Fraction of cells that passed QC."""
        return self.passed_cells / max(1, self.total_cells)


class CellQC:
    """
    Comprehensive cell-level quality control for multiplex imaging.

    Filters cells based on:
    - Intensity outliers
    - Morphological features (area, eccentricity)
    - Doublet detection
    - Edge artifacts
    - Technical artifacts

    Example
    -------
    >>> qc = CellQC()
    >>> result = qc.filter(
    ...     cell_data,
    ...     intensity_columns=["DAPI", "CD3", "CD20"],
    ...     area_column="area",
    ... )
    >>> clean_data = cell_data[result.passed_indices]
    """

    def __init__(
        self,
        intensity_mad_threshold: float = 5.0,
        min_area: float = 50,
        max_area: float = 5000,
        max_eccentricity: float = 0.95,
        edge_margin: int = 10,
        doublet_area_factor: float = 1.8,
    ):
        """
        Initialize Cell QC.

        Parameters
        ----------
        intensity_mad_threshold : float
            MAD threshold for intensity outlier detection
        min_area : float
            Minimum cell area in pixels
        max_area : float
            Maximum cell area in pixels
        max_eccentricity : float
            Maximum eccentricity (0=circle, 1=line)
        edge_margin : int
            Margin from image edge to filter cells
        doublet_area_factor : float
            Area factor above median to flag as potential doublet
        """
        self.intensity_mad_threshold = intensity_mad_threshold
        self.min_area = min_area
        self.max_area = max_area
        self.max_eccentricity = max_eccentricity
        self.edge_margin = edge_margin
        self.doublet_area_factor = doublet_area_factor

    def filter(
        self,
        cell_data: ArrayLike,
        intensity_columns: list[str] | None = None,
        area_column: str | None = None,
        eccentricity_column: str | None = None,
        x_column: str | None = None,
        y_column: str | None = None,
        image_shape: tuple[int, int] | None = None,
    ) -> CellQCResult:
        """
        Apply all QC filters to cell data.

        Parameters
        ----------
        cell_data : array-like or DataFrame
            Cell measurements with one row per cell
        intensity_columns : list of str, optional
            Column names for intensity measurements
        area_column : str, optional
            Column name for cell area
        eccentricity_column : str, optional
            Column name for cell eccentricity
        x_column, y_column : str, optional
            Column names for cell centroid coordinates
        image_shape : tuple, optional
            Image dimensions for edge filtering

        Returns
        -------
        CellQCResult
            Filtering results
        """
        # Convert to numpy if DataFrame
        if hasattr(cell_data, "values"):
            df = cell_data
            data = df.values
        else:
            df = None
            data = np.asarray(cell_data)

        n_cells = len(data)
        passed = np.ones(n_cells, dtype=bool)
        filter_counts = {}
        filter_reasons = {i: [] for i in range(n_cells)}

        # Get column indices
        def get_col_idx(name):
            if name is None:
                return None
            if df is not None:
                return df.columns.get_loc(name) if name in df.columns else None
            return None

        # 1. Intensity outlier filtering
        if intensity_columns:
            for col in intensity_columns:
                col_idx = get_col_idx(col)
                if col_idx is not None:
                    intensities = data[:, col_idx].astype(float)
                    outliers = self._detect_intensity_outliers(intensities)
                    n_outliers = np.sum(outliers)
                    filter_counts[f"intensity_outlier_{col}"] = n_outliers

                    for idx in np.where(outliers)[0]:
                        filter_reasons[idx].append(f"intensity_outlier_{col}")

                    passed &= ~outliers

        # 2. Area filtering
        if area_column:
            area_idx = get_col_idx(area_column)
            if area_idx is not None:
                areas = data[:, area_idx].astype(float)

                too_small = areas < self.min_area
                too_large = areas > self.max_area
                doublets = self._detect_doublets_by_area(areas)

                filter_counts["too_small"] = np.sum(too_small)
                filter_counts["too_large"] = np.sum(too_large)
                filter_counts["doublets"] = np.sum(doublets & ~too_large)

                for idx in np.where(too_small)[0]:
                    filter_reasons[idx].append("too_small")
                for idx in np.where(too_large)[0]:
                    filter_reasons[idx].append("too_large")
                for idx in np.where(doublets & ~too_large)[0]:
                    filter_reasons[idx].append("doublet")

                passed &= ~too_small & ~too_large & ~doublets

        # 3. Eccentricity filtering
        if eccentricity_column:
            ecc_idx = get_col_idx(eccentricity_column)
            if ecc_idx is not None:
                eccentricities = data[:, ecc_idx].astype(float)
                too_eccentric = eccentricities > self.max_eccentricity
                filter_counts["too_eccentric"] = np.sum(too_eccentric)

                for idx in np.where(too_eccentric)[0]:
                    filter_reasons[idx].append("too_eccentric")

                passed &= ~too_eccentric

        # 4. Edge filtering
        if x_column and y_column and image_shape:
            x_idx = get_col_idx(x_column)
            y_idx = get_col_idx(y_column)
            if x_idx is not None and y_idx is not None:
                x_coords = data[:, x_idx].astype(float)
                y_coords = data[:, y_idx].astype(float)
                h, w = image_shape

                on_edge = (
                    (x_coords < self.edge_margin)
                    | (x_coords > w - self.edge_margin)
                    | (y_coords < self.edge_margin)
                    | (y_coords > h - self.edge_margin)
                )
                filter_counts["on_edge"] = np.sum(on_edge)

                for idx in np.where(on_edge)[0]:
                    filter_reasons[idx].append("on_edge")

                passed &= ~on_edge

        # Build result
        passed_indices = np.where(passed)[0]
        filtered_indices = np.where(~passed)[0]

        # Only include cells that were actually filtered
        filter_reasons = {
            idx: reasons
            for idx, reasons in filter_reasons.items()
            if reasons and idx in filtered_indices
        }

        return CellQCResult(
            total_cells=n_cells,
            passed_cells=int(np.sum(passed)),
            filtered_cells=int(np.sum(~passed)),
            filter_counts=filter_counts,
            passed_indices=passed_indices,
            filtered_indices=filtered_indices,
            filter_reasons=filter_reasons,
        )

    def _detect_intensity_outliers(self, intensities: np.ndarray) -> np.ndarray:
        """Detect intensity outliers using MAD."""
        median = np.median(intensities)
        mad = np.median(np.abs(intensities - median)) * 1.4826

        if mad < 1e-8:
            return np.zeros(len(intensities), dtype=bool)

        z_scores = np.abs(intensities - median) / mad
        return z_scores > self.intensity_mad_threshold

    def _detect_doublets_by_area(self, areas: np.ndarray) -> np.ndarray:
        """Detect potential doublets based on unusually large area."""
        median_area = np.median(areas)
        return areas > median_area * self.doublet_area_factor


def detect_outliers(
    data: np.ndarray,
    method: str = "mad",
    threshold: float = 5.0,
    per_column: bool = True,
) -> np.ndarray:
    """
    Detect outliers in numeric data.

    Parameters
    ----------
    data : np.ndarray
        Data array (1D or 2D)
    method : str
        Detection method: "mad" (median absolute deviation),
        "iqr" (interquartile range), "zscore"
    threshold : float
        Number of MADs/IQRs/standard deviations for outlier threshold
    per_column : bool
        For 2D data, detect outliers per column

    Returns
    -------
    np.ndarray
        Boolean array marking outliers
    """
    data = np.asarray(data, dtype=float)

    if data.ndim == 1:
        return _detect_outliers_1d(data, method, threshold)

    if per_column:
        outliers = np.zeros(data.shape, dtype=bool)
        for i in range(data.shape[1]):
            outliers[:, i] = _detect_outliers_1d(data[:, i], method, threshold)
        # A row is an outlier if ANY column is an outlier
        return np.any(outliers, axis=1)
    else:
        # Flatten and detect
        flat_outliers = _detect_outliers_1d(data.ravel(), method, threshold)
        outlier_rows = np.any(flat_outliers.reshape(data.shape), axis=1)
        return outlier_rows


def _detect_outliers_1d(
    data: np.ndarray,
    method: str,
    threshold: float,
) -> np.ndarray:
    """Detect outliers in 1D data."""
    if method == "mad":
        median = np.median(data)
        mad = np.median(np.abs(data - median)) * 1.4826
        if mad < 1e-8:
            return np.zeros(len(data), dtype=bool)
        z = np.abs(data - median) / mad
        return z > threshold

    elif method == "iqr":
        q1, q3 = np.percentile(data, [25, 75])
        iqr = q3 - q1
        lower = q1 - threshold * iqr
        upper = q3 + threshold * iqr
        return (data < lower) | (data > upper)

    elif method == "zscore":
        mean = np.mean(data)
        std = np.std(data)
        if std < 1e-8:
            return np.zeros(len(data), dtype=bool)
        z = np.abs(data - mean) / std
        return z > threshold

    else:
        raise ValueError(f"Unknown method: {method}")


def filter_by_intensity(
    intensities: np.ndarray,
    min_intensity: float | None = None,
    max_intensity: float | None = None,
    min_percentile: float | None = None,
    max_percentile: float | None = None,
) -> np.ndarray:
    """
    Filter cells by intensity thresholds.

    Parameters
    ----------
    intensities : np.ndarray
        Intensity values
    min_intensity, max_intensity : float, optional
        Absolute intensity thresholds
    min_percentile, max_percentile : float, optional
        Percentile thresholds (0-100)

    Returns
    -------
    np.ndarray
        Boolean array of cells passing filter
    """
    passed = np.ones(len(intensities), dtype=bool)

    if min_percentile is not None:
        min_val = np.percentile(intensities, min_percentile)
        passed &= intensities >= min_val

    if max_percentile is not None:
        max_val = np.percentile(intensities, max_percentile)
        passed &= intensities <= max_val

    if min_intensity is not None:
        passed &= intensities >= min_intensity

    if max_intensity is not None:
        passed &= intensities <= max_intensity

    return passed


def filter_by_morphology(
    areas: np.ndarray | None = None,
    eccentricities: np.ndarray | None = None,
    min_area: float = 50,
    max_area: float = 5000,
    max_eccentricity: float = 0.95,
) -> np.ndarray:
    """
    Filter cells by morphological features.

    Parameters
    ----------
    areas : np.ndarray, optional
        Cell areas
    eccentricities : np.ndarray, optional
        Cell eccentricities
    min_area, max_area : float
        Area bounds
    max_eccentricity : float
        Maximum eccentricity

    Returns
    -------
    np.ndarray
        Boolean array of cells passing filter
    """
    n = len(areas) if areas is not None else len(eccentricities)
    passed = np.ones(n, dtype=bool)

    if areas is not None:
        passed &= (areas >= min_area) & (areas <= max_area)

    if eccentricities is not None:
        passed &= eccentricities <= max_eccentricity

    return passed


def detect_doublets(
    areas: np.ndarray,
    intensities: np.ndarray | None = None,
    area_factor: float = 1.8,
    intensity_factor: float = 1.5,
) -> np.ndarray:
    """
    Detect potential doublets (two cells merged as one).

    Parameters
    ----------
    areas : np.ndarray
        Cell areas
    intensities : np.ndarray, optional
        Total intensities (e.g., DAPI)
    area_factor : float
        Factor above median area to flag as doublet
    intensity_factor : float
        Factor above median intensity to flag as doublet

    Returns
    -------
    np.ndarray
        Boolean array marking potential doublets
    """
    median_area = np.median(areas)
    doublets = areas > median_area * area_factor

    if intensities is not None:
        median_intensity = np.median(intensities)
        intensity_high = intensities > median_intensity * intensity_factor
        doublets &= intensity_high  # Both criteria must be met

    return doublets


def detect_spatial_outliers(
    coordinates: np.ndarray,
    values: np.ndarray,
    n_neighbors: int = 10,
    threshold: float = 3.0,
) -> np.ndarray:
    """
    Detect cells with values that differ significantly from their neighbors.

    Useful for detecting technical artifacts that affect isolated cells.

    Parameters
    ----------
    coordinates : np.ndarray
        Cell coordinates (N x 2)
    values : np.ndarray
        Values to compare (N,)
    n_neighbors : int
        Number of neighbors to consider
    threshold : float
        Number of MADs to flag as outlier

    Returns
    -------
    np.ndarray
        Boolean array marking spatial outliers
    """
    n_cells = len(coordinates)
    if n_cells < n_neighbors + 1:
        return np.zeros(n_cells, dtype=bool)

    # Compute pairwise distances
    distances = cdist(coordinates, coordinates)

    outliers = np.zeros(n_cells, dtype=bool)

    for i in range(n_cells):
        # Get nearest neighbors
        neighbor_indices = np.argsort(distances[i])[1 : n_neighbors + 1]

        # Compare value to neighbors
        neighbor_values = values[neighbor_indices]
        neighbor_median = np.median(neighbor_values)
        neighbor_mad = np.median(np.abs(neighbor_values - neighbor_median)) * 1.4826

        if neighbor_mad > 1e-8:
            z = np.abs(values[i] - neighbor_median) / neighbor_mad
            outliers[i] = z > threshold

    return outliers
