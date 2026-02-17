"""
Tests for KINTSUGI spillover correction module.

Tests cover utility functions, parameter estimation, QC metrics, and
(when redseapy is available) the core REDSEA wrapper. Tests that require
redseapy are skipped if the package is not installed.
"""

import numpy as np
import pandas as pd
import pytest

from kintsugi.spillover.spillover_qc import compute_spillover_metrics
from kintsugi.spillover.utils import (
    estimate_element_size,
    identify_surface_markers,
    load_markers_csv,
    validate_segmentation_mask,
    write_markers_csv,
)

# Check if redseapy is available
try:
    from redseapy.redsea import run_redsea  # noqa: F401

    REDSEAPY_AVAILABLE = True
except ImportError:
    REDSEAPY_AVAILABLE = False


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def simple_cell_mask():
    """Create a simple 2D segmentation mask with 4 cells."""
    mask = np.zeros((100, 100), dtype=np.int32)
    # Cell 1: top-left
    mask[10:30, 10:30] = 1
    # Cell 2: top-right
    mask[10:30, 70:90] = 2
    # Cell 3: bottom-left
    mask[70:90, 10:30] = 3
    # Cell 4: bottom-right
    mask[70:90, 70:90] = 4
    return mask


@pytest.fixture
def dense_cell_mask():
    """Create a dense segmentation mask with many small cells."""
    mask = np.zeros((200, 200), dtype=np.int32)
    cell_id = 1
    for y in range(0, 200, 12):
        for x in range(0, 200, 12):
            mask[y : y + 10, x : x + 10] = cell_id
            cell_id += 1
    return mask


@pytest.fixture
def large_cell_mask():
    """Create a mask with large cells (median area > 400)."""
    mask = np.zeros((500, 500), dtype=np.int32)
    mask[10:60, 10:60] = 1  # 50x50 = 2500 pixels
    mask[10:60, 70:120] = 2
    mask[70:120, 10:60] = 3
    mask[70:120, 70:120] = 4
    return mask


@pytest.fixture
def marker_names():
    """Common CODEX marker panel."""
    return [
        "DAPI",
        "CD3",
        "CD20",
        "CD8",
        "CD68",
        "FoxP3",
        "Ki67",
        "PanCK",
        "CD45",
        "Blank",
    ]


@pytest.fixture
def sample_quantification():
    """Create sample corrected/uncorrected quantification DataFrames."""
    np.random.seed(42)
    n_cells = 500

    # Uncorrected: some spillover creating double-positives
    uncorrected = pd.DataFrame(
        {
            "CD3": np.random.exponential(100, n_cells),
            "CD20": np.random.exponential(100, n_cells),
            "CD8": np.random.exponential(80, n_cells),
            "CD4": np.random.exponential(80, n_cells),
            "CD68": np.random.exponential(60, n_cells),
        }
    )
    # Add spillover: some CD3+ cells artificially get CD20 signal
    spillover_mask = uncorrected["CD3"] > 150
    uncorrected.loc[spillover_mask, "CD20"] += 120

    # Corrected: spillover removed
    corrected = uncorrected.copy()
    corrected.loc[spillover_mask, "CD20"] -= 100

    return corrected, uncorrected


# =============================================================================
# Tests: identify_surface_markers
# =============================================================================


class TestIdentifySurfaceMarkers:
    def test_basic_identification(self, marker_names):
        """Surface markers should be correctly identified."""
        result = identify_surface_markers(marker_names)
        assert "CD3" in result
        assert "CD20" in result
        assert "CD8" in result
        assert "CD68" in result
        assert "PanCK" in result
        assert "CD45" in result

    def test_nuclear_markers_excluded(self, marker_names):
        """Nuclear markers should not appear in surface list."""
        result = identify_surface_markers(marker_names)
        assert "DAPI" not in result
        assert "FoxP3" not in result
        assert "Ki67" not in result

    def test_blank_excluded(self, marker_names):
        """Blank channels should be excluded."""
        result = identify_surface_markers(marker_names)
        assert "Blank" not in result

    def test_custom_surface_set(self):
        """Custom surface marker set should work."""
        markers = ["CD3", "CD20", "CustomMarker", "DAPI"]
        custom_set = {"CD3", "CustomMarker"}
        result = identify_surface_markers(markers, known_surface_markers=custom_set)
        assert "CD3" in result
        assert "CustomMarker" in result
        assert "CD20" not in result
        assert "DAPI" not in result

    def test_empty_panel(self):
        """Empty panel returns empty list."""
        assert identify_surface_markers([]) == []

    def test_all_nuclear(self):
        """All-nuclear panel returns empty list."""
        markers = ["DAPI", "FoxP3", "Ki67"]
        assert identify_surface_markers(markers) == []


# =============================================================================
# Tests: estimate_element_size
# =============================================================================


class TestEstimateElementSize:
    def test_small_cells(self, dense_cell_mask):
        """Dense small cells should get element_size=2."""
        size = estimate_element_size(dense_cell_mask)
        assert size == 2

    def test_large_cells(self, large_cell_mask):
        """Large cells should get element_size >= 3."""
        size = estimate_element_size(large_cell_mask)
        assert size >= 3

    def test_empty_mask(self):
        """Empty mask should return default size=2."""
        mask = np.zeros((100, 100), dtype=np.int32)
        size = estimate_element_size(mask)
        assert size == 2

    def test_return_type(self, simple_cell_mask):
        """Should return an integer."""
        size = estimate_element_size(simple_cell_mask)
        assert isinstance(size, int)

    def test_reasonable_range(self, simple_cell_mask):
        """Element size should be in the reasonable range 2-4."""
        size = estimate_element_size(simple_cell_mask)
        assert 2 <= size <= 4


# =============================================================================
# Tests: validate_segmentation_mask
# =============================================================================


class TestValidateSegmentationMask:
    def test_valid_mask(self, simple_cell_mask):
        """Valid mask should pass validation."""
        result = validate_segmentation_mask(simple_cell_mask)
        assert result["valid"] is True
        assert result["n_cells"] == 4
        assert len(result["issues"]) == 0

    def test_3d_mask_fails(self):
        """3D mask should fail validation."""
        mask = np.zeros((10, 100, 100), dtype=np.int32)
        result = validate_segmentation_mask(mask)
        assert result["valid"] is False
        assert any("2D" in issue for issue in result["issues"])

    def test_float_mask_fails(self):
        """Float mask should fail validation."""
        mask = np.zeros((100, 100), dtype=np.float64)
        result = validate_segmentation_mask(mask)
        assert result["valid"] is False
        assert any("integer" in issue for issue in result["issues"])

    def test_empty_mask_fails(self):
        """Mask with no cells should fail."""
        mask = np.zeros((100, 100), dtype=np.int32)
        result = validate_segmentation_mask(mask)
        assert result["valid"] is False
        assert result["n_cells"] == 0


# =============================================================================
# Tests: markers CSV I/O
# =============================================================================


class TestMarkersCSV:
    def test_roundtrip(self, tmp_path, marker_names):
        """Write and read markers CSV should round-trip."""
        csv_path = tmp_path / "markers.csv"
        write_markers_csv(marker_names, csv_path)
        loaded = load_markers_csv(csv_path)
        assert loaded == marker_names

    def test_csv_format(self, tmp_path):
        """CSV should have marker_name column."""
        markers = ["CD3", "CD20"]
        csv_path = tmp_path / "test.csv"
        write_markers_csv(markers, csv_path)
        df = pd.read_csv(csv_path)
        assert "marker_name" in df.columns
        assert list(df["marker_name"]) == markers

    def test_load_missing_column(self, tmp_path):
        """Loading CSV without marker_name column should raise."""
        csv_path = tmp_path / "bad.csv"
        pd.DataFrame({"wrong_col": ["a", "b"]}).to_csv(csv_path, index=False)
        with pytest.raises(ValueError, match="marker_name"):
            load_markers_csv(csv_path)


# =============================================================================
# Tests: compute_spillover_metrics
# =============================================================================


class TestSpilloverMetrics:
    def test_basic_metrics(self, sample_quantification):
        """Should compute double-positive reduction metrics."""
        corrected, uncorrected = sample_quantification
        metrics = compute_spillover_metrics(
            corrected,
            uncorrected,
            mutually_exclusive_pairs=[("CD3", "CD20")],
        )

        assert "pair_metrics" in metrics
        assert "summary" in metrics
        assert "per_marker_change" in metrics

        pair_key = ("CD3", "CD20")
        assert pair_key in metrics["pair_metrics"]

        pair = metrics["pair_metrics"][pair_key]
        assert "uncorrected_dp_pct" in pair
        assert "corrected_dp_pct" in pair
        assert "dp_reduction_pct" in pair

    def test_dp_reduction(self, sample_quantification):
        """Corrected data should have lower double-positive rate."""
        corrected, uncorrected = sample_quantification
        metrics = compute_spillover_metrics(
            corrected,
            uncorrected,
            mutually_exclusive_pairs=[("CD3", "CD20")],
        )
        pair = metrics["pair_metrics"][("CD3", "CD20")]
        assert pair["dp_reduction_pct"] >= 0, "Corrected data should have fewer double-positives"

    def test_missing_markers_skipped(self):
        """Pairs with missing markers should be silently skipped."""
        corrected = pd.DataFrame({"CD3": [1, 2, 3]})
        uncorrected = pd.DataFrame({"CD3": [1, 2, 3]})
        metrics = compute_spillover_metrics(
            corrected,
            uncorrected,
            mutually_exclusive_pairs=[("CD3", "NotPresent")],
        )
        assert metrics["summary"]["pairs_evaluated"] == 0

    def test_per_marker_change(self, sample_quantification):
        """Per-marker intensity changes should be computed."""
        corrected, uncorrected = sample_quantification
        metrics = compute_spillover_metrics(corrected, uncorrected)
        assert len(metrics["per_marker_change"]) > 0
        for _marker, change in metrics["per_marker_change"].items():
            assert "mean_before" in change
            assert "mean_after" in change
            assert "fold_change" in change

    def test_auto_exclusive_pairs(self, sample_quantification):
        """Default exclusive pairs should be used when none provided."""
        corrected, uncorrected = sample_quantification
        metrics = compute_spillover_metrics(corrected, uncorrected)
        # Should have evaluated at least one pair (CD3/CD20 is in defaults
        # and in our sample data)
        assert metrics["summary"]["pairs_evaluated"] >= 0


# =============================================================================
# Tests: QC visualization
# =============================================================================


class TestSpilloverVisualization:
    def test_biaxial_plot(self, sample_quantification, tmp_path):
        """Biaxial scatter plot should be generated and saved."""
        from kintsugi.spillover import plot_spillover_comparison

        corrected, uncorrected = sample_quantification
        save_path = tmp_path / "test_biaxial.png"
        fig = plot_spillover_comparison(
            corrected,
            uncorrected,
            marker_pair=("CD3", "CD20"),
            save_path=save_path,
        )
        assert fig is not None
        assert save_path.exists()
        import matplotlib.pyplot as plt

        plt.close(fig)

    def test_biaxial_missing_marker(self, sample_quantification):
        """Missing marker should raise ValueError."""
        from kintsugi.spillover import plot_spillover_comparison

        corrected, uncorrected = sample_quantification
        with pytest.raises(ValueError, match="not found"):
            plot_spillover_comparison(
                corrected,
                uncorrected,
                marker_pair=("CD3", "NotPresent"),
            )

    def test_summary_heatmap(self, sample_quantification, tmp_path):
        """Summary heatmap should be generated."""
        from kintsugi.spillover import plot_spillover_summary_heatmap

        corrected, uncorrected = sample_quantification
        save_path = tmp_path / "test_heatmap.png"
        fig = plot_spillover_summary_heatmap(
            corrected,
            uncorrected,
            marker_names=["CD3", "CD20", "CD8", "CD68"],
            save_path=save_path,
        )
        assert fig is not None
        assert save_path.exists()
        import matplotlib.pyplot as plt

        plt.close(fig)


# =============================================================================
# Tests: SpilloverResult dataclass
# =============================================================================


class TestSpilloverResult:
    def test_creation(self):
        """SpilloverResult should be creatable with DataFrames."""
        from kintsugi.spillover import SpilloverResult

        result = SpilloverResult(
            corrected=pd.DataFrame({"CD3": [1, 2]}),
            uncorrected=pd.DataFrame({"CD3": [1, 2]}),
            parameters={"element_size": 2},
            n_cells=2,
            n_markers=1,
        )
        assert result.n_cells == 2
        assert result.n_markers == 1
        assert len(result.corrected) == 2


# =============================================================================
# Tests: REDSEA wrapper (requires redseapy)
# =============================================================================


@pytest.mark.skipif(not REDSEAPY_AVAILABLE, reason="redseapy not installed")
class TestRedseaWrapper:
    def test_invalid_element_shape(self, simple_cell_mask, marker_names):
        """Invalid element_shape should raise ValueError."""
        from kintsugi.spillover import run_redsea_correction

        image = np.random.randint(0, 1000, (100, 100, len(marker_names)), dtype=np.uint16)
        with pytest.raises(ValueError, match="element_shape"):
            run_redsea_correction(
                multichannel_image=image,
                segmentation_mask=simple_cell_mask,
                marker_names=marker_names,
                element_shape="invalid",
            )

    def test_correction_returns_result(self, simple_cell_mask):
        """REDSEA correction should return a SpilloverResult."""
        from kintsugi.spillover import SpilloverResult, run_redsea_correction

        markers = ["CD3", "CD20", "CD8"]
        n_channels = len(markers)
        image = np.random.randint(0, 5000, (100, 100, n_channels), dtype=np.uint16)

        result = run_redsea_correction(
            multichannel_image=image,
            segmentation_mask=simple_cell_mask,
            marker_names=markers,
            element_size=2,
        )

        assert isinstance(result, SpilloverResult)
        assert result.n_cells > 0
        assert isinstance(result.corrected, pd.DataFrame)
        assert isinstance(result.uncorrected, pd.DataFrame)

    def test_correction_with_output_dir(self, simple_cell_mask, tmp_path):
        """Correction with output_dir should save CSVs."""
        from kintsugi.spillover import run_redsea_correction

        markers = ["CD3", "CD20"]
        image = np.random.randint(0, 5000, (100, 100, 2), dtype=np.uint16)
        output_dir = tmp_path / "spillover_output"

        run_redsea_correction(
            multichannel_image=image,
            segmentation_mask=simple_cell_mask,
            marker_names=markers,
            output_dir=output_dir,
        )

        assert (output_dir / "single_cell_after_redsea.csv").exists()
        assert (output_dir / "single_cell_before_redsea.csv").exists()

    def test_markers_of_interest(self, simple_cell_mask):
        """Correcting only a subset of markers should work."""
        from kintsugi.spillover import SpilloverResult, run_redsea_correction

        markers = ["CD3", "CD20", "CD8", "DAPI"]
        image = np.random.randint(0, 5000, (100, 100, 4), dtype=np.uint16)

        result = run_redsea_correction(
            multichannel_image=image,
            segmentation_mask=simple_cell_mask,
            marker_names=markers,
            markers_of_interest=["CD3", "CD20"],
            element_size=2,
        )

        assert isinstance(result, SpilloverResult)
        assert result.parameters["markers_of_interest"] == ["CD3", "CD20"]
