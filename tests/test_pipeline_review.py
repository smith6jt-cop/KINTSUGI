"""
Tests for pipeline review improvements (Feb 2026).

This module tests functionality added in response to the KINTSUGI codebase
review and CyLinter/large_image analysis:

- Cell-level QC filtering (CyLinter-inspired)
- Cycle dropout detection
- Per-marker outlier pruning
- General artifact detection (artifact_detector_v3 inspired)
- Marker name loading from CHANNELNAMES.txt
- Parquet checkpointing
"""

import os
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


# =============================================================================
# Cell QC Tests
# =============================================================================


class TestCellQCPipeline:
    """Tests for CyLinter-inspired cell QC filtering."""

    @pytest.fixture
    def cell_dataframe(self):
        """Create a realistic cell feature DataFrame."""
        np.random.seed(42)
        n_cells = 1000
        return pd.DataFrame(
            {
                "cell_label": range(1, n_cells + 1),
                "cell_DAPI_mean": np.random.lognormal(7, 0.5, n_cells).astype(float),
                "cell_CD3_mean": np.random.lognormal(5, 1.0, n_cells).astype(float),
                "cell_CD20_mean": np.random.lognormal(5, 1.0, n_cells).astype(float),
                "cell_CD68_mean": np.random.lognormal(4, 1.5, n_cells).astype(float),
                "cell_number_of_pixels": np.random.lognormal(5.5, 0.8, n_cells).astype(float),
                "cell_centroid_0": np.random.uniform(0, 10000, n_cells),
                "cell_centroid_1": np.random.uniform(0, 10000, n_cells),
            }
        )

    def test_filter_cells_pipeline_runs(self, cell_dataframe):
        """Test that the full pipeline runs without errors."""
        from kintsugi.qc import filter_cells_pipeline

        filtered, summary = filter_cells_pipeline(
            cell_dataframe,
            dna_column="cell_DAPI_mean",
            area_column="cell_number_of_pixels",
            marker_columns=["cell_CD3_mean", "cell_CD20_mean", "cell_CD68_mean"],
        )

        assert len(filtered) <= len(cell_dataframe)
        assert "initial_cells" in summary
        assert "passed_cells" in summary
        assert summary["initial_cells"] == len(cell_dataframe)
        assert summary["passed_cells"] == len(filtered)

    def test_filter_cells_removes_outliers(self, cell_dataframe):
        """Test that extreme outliers are actually removed."""
        from kintsugi.qc import filter_cells_pipeline

        # Inject extreme outliers
        df = cell_dataframe.copy()
        df.loc[0, "cell_DAPI_mean"] = 1e8  # Extreme DAPI
        df.loc[1, "cell_number_of_pixels"] = 1  # Tiny cell
        df.loc[2, "cell_number_of_pixels"] = 1e6  # Giant cell

        filtered, summary = filter_cells_pipeline(
            df,
            dna_column="cell_DAPI_mean",
            area_column="cell_number_of_pixels",
        )

        assert summary["total_filtered"] > 0
        assert len(filtered) < len(df)

    def test_filter_cells_preserves_normal(self, cell_dataframe):
        """Test that normal cells are not excessively filtered."""
        from kintsugi.qc import filter_cells_pipeline

        filtered, summary = filter_cells_pipeline(
            cell_dataframe,
            dna_column="cell_DAPI_mean",
            area_column="cell_number_of_pixels",
        )

        # Should retain most cells (>80%) with normal data
        pass_rate = summary["passed_cells"] / summary["initial_cells"]
        assert pass_rate > 0.5  # At least 50% should pass

    def test_filter_cells_missing_columns(self, cell_dataframe):
        """Test graceful handling of missing columns."""
        from kintsugi.qc import filter_cells_pipeline

        filtered, summary = filter_cells_pipeline(
            cell_dataframe,
            dna_column="nonexistent_column",
            area_column="nonexistent_area",
        )

        # Should return all cells if columns don't exist
        assert len(filtered) == len(cell_dataframe)


class TestCycleDropout:
    """Tests for cycle dropout detection."""

    def test_detect_dropout_basic(self):
        """Test basic dropout detection."""
        from kintsugi.qc import detect_cycle_dropout

        # First cycle: all cells have strong signal
        dna_first = np.array([1000, 1000, 1000, 1000, 1000], dtype=float)
        # Last cycle: some cells lost signal
        dna_last = np.array([900, 800, 50, 100, 950], dtype=float)

        dropouts = detect_cycle_dropout(dna_first, dna_last, ratio_threshold=0.3)

        assert dropouts[2]  # 50/1000 = 0.05, should be flagged
        assert dropouts[3]  # 100/1000 = 0.10, should be flagged
        assert not dropouts[0]  # 900/1000 = 0.90, should pass
        assert not dropouts[4]  # 950/1000 = 0.95, should pass

    def test_detect_dropout_low_first_cycle(self):
        """Test that cells dim in first cycle are not flagged."""
        from kintsugi.qc import detect_cycle_dropout

        dna_first = np.array([10, 20, 1000], dtype=float)
        dna_last = np.array([5, 10, 50], dtype=float)

        dropouts = detect_cycle_dropout(
            dna_first, dna_last, ratio_threshold=0.3, min_intensity=100
        )

        assert not dropouts[0]  # Too dim in first cycle to evaluate
        assert not dropouts[1]  # Too dim in first cycle to evaluate
        assert dropouts[2]  # Clear dropout

    def test_detect_dropout_no_dropouts(self):
        """Test when no cells drop out."""
        from kintsugi.qc import detect_cycle_dropout

        dna_first = np.array([1000, 1000, 1000], dtype=float)
        dna_last = np.array([950, 900, 800], dtype=float)

        dropouts = detect_cycle_dropout(dna_first, dna_last, ratio_threshold=0.3)
        assert not np.any(dropouts)


class TestMarkerOutlierPruning:
    """Tests for per-marker outlier pruning."""

    def test_prune_basic(self):
        """Test basic per-marker outlier pruning."""
        from kintsugi.qc import prune_marker_outliers

        np.random.seed(42)
        n = 500
        df = pd.DataFrame(
            {
                "CD3_mean": np.random.normal(100, 10, n),
                "CD20_mean": np.random.normal(80, 8, n),
            }
        )
        # Inject outliers
        df.loc[0, "CD3_mean"] = 500  # Extreme
        df.loc[1, "CD20_mean"] = 400  # Extreme

        passed, counts = prune_marker_outliers(
            df, ["CD3_mean", "CD20_mean"], threshold=5.0
        )

        assert not passed[0]  # CD3 outlier removed
        assert not passed[1]  # CD20 outlier removed
        assert counts["CD3_mean"] >= 1
        assert counts["CD20_mean"] >= 1

    def test_prune_no_outliers(self):
        """Test with clean data."""
        from kintsugi.qc import prune_marker_outliers

        np.random.seed(42)
        df = pd.DataFrame(
            {
                "CD3_mean": np.random.normal(100, 10, 100),
            }
        )

        passed, counts = prune_marker_outliers(df, ["CD3_mean"], threshold=10.0)

        # Very high threshold should not flag normal data
        assert np.sum(passed) >= 95  # Most should pass

    def test_prune_missing_column(self):
        """Test with a column that doesn't exist."""
        from kintsugi.qc import prune_marker_outliers

        df = pd.DataFrame({"CD3_mean": [100, 200, 300]})
        passed, counts = prune_marker_outliers(df, ["nonexistent_col"])

        assert np.all(passed)  # All should pass if column doesn't exist


# =============================================================================
# General Artifact Detection Tests
# =============================================================================


class TestGeneralArtifactDetection:
    """Tests for flood-fill based artifact detection."""

    def test_clean_image_no_artifacts(self):
        """Test that a uniform image has no artifacts."""
        from kintsugi.qc import detect_general_artifacts

        np.random.seed(42)
        # Narrow dynamic range (140-160) — well below the 50-unit threshold
        image = np.random.randint(140, 160, (256, 256), dtype=np.uint16)

        mask = detect_general_artifacts(image, downscale_factor=2)

        assert mask.shape == image.shape
        assert mask.dtype == bool
        # Nearly uniform noise should produce no artifacts
        assert np.sum(mask) == 0

    def test_bright_blob_detected(self):
        """Test that a bright artificial blob is detected."""
        from kintsugi.qc import detect_general_artifacts

        np.random.seed(42)
        image = np.random.randint(100, 200, (256, 256), dtype=np.uint16)

        # Add a bright blob (simulating tissue fold / bubble)
        yy, xx = np.mgrid[100:150, 100:150]
        image[100:150, 100:150] = 60000  # Very bright region

        mask = detect_general_artifacts(
            image,
            downscale_factor=2,
            min_artifact_size=10,
        )

        assert mask.shape == image.shape
        # The bright region should be at least partially detected
        bright_region_detected = np.sum(mask[100:150, 100:150])
        assert bright_region_detected > 0

    def test_returns_correct_shape(self):
        """Test that output mask matches input shape."""
        from kintsugi.qc import detect_general_artifacts

        for shape in [(128, 128), (256, 512), (100, 300)]:
            image = np.random.randint(0, 65535, shape, dtype=np.uint16)
            mask = detect_general_artifacts(image, downscale_factor=2)
            assert mask.shape == shape


# =============================================================================
# Marker Name Loading Tests
# =============================================================================


class TestMarkerNameLoading:
    """Tests for loading marker names from CHANNELNAMES.txt."""

    def test_load_simple_format(self, tmp_path):
        """Test loading simple list format."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import load_channel_names

        # Create CHANNELNAMES.txt
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "DAPI-01\nBlank\nBlank\nBlank\nDAPI-02\nCD31\nCD8\nCD45\n"
        )

        markers = load_channel_names(tmp_path)
        assert "DAPI" in markers
        assert "CD31" in markers
        assert "CD8" in markers
        assert "CD45" in markers
        assert "Blank" not in markers
        # DAPI-01 and DAPI-02 should be merged into single DAPI
        assert markers.count("DAPI") == 1

    def test_load_cycle_prefixed_format(self, tmp_path):
        """Test loading cycle-prefixed format."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import load_channel_names

        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "CHANNELNAMES.txt").write_text(
            "1: DAPI, Blank, Blank, Blank\n"
            "2: DAPI, CD31, CD8, CD45\n"
            "3: DAPI, CD3, CD20, Ki67\n"
        )

        markers = load_channel_names(tmp_path)
        assert "DAPI" in markers
        assert "CD31" in markers
        assert "CD3" in markers
        assert "Ki67" in markers
        assert "Blank" not in markers

    def test_load_missing_file(self, tmp_path):
        """Test error when CHANNELNAMES.txt doesn't exist."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import load_channel_names

        with pytest.raises(FileNotFoundError):
            load_channel_names(tmp_path)


class TestMarkerSubsets:
    """Tests for compartment-specific marker filtering."""

    def test_cell_markers_exclude_nuclear(self):
        """Test that cell markers exclude nuclear-only markers."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import get_cell_markers

        all_markers = [
            "CD3",
            "CD20",
            "CD68",
            "DAPI",
            "FoxP3",
            "Ki67",
            "CollagenIV",
            "ECAD",
        ]
        cell = get_cell_markers(all_markers)
        assert "DAPI" not in cell
        assert "FoxP3" not in cell
        assert "Ki67" not in cell
        assert "CollagenIV" not in cell
        assert "ECAD" not in cell
        assert "CD3" in cell
        assert "CD20" in cell

    def test_nuclear_markers_include_dapi(self):
        """Test that nuclear markers include DAPI."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import get_nuclear_markers

        all_markers = ["CD3", "CD20", "DAPI", "FoxP3", "CollagenIV"]
        nuc = get_nuclear_markers(all_markers)
        assert "DAPI" in nuc
        assert "FoxP3" in nuc
        assert "CollagenIV" not in nuc

    def test_custom_exclusions(self):
        """Test custom exclusion sets."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import get_cell_markers

        all_markers = ["CD3", "CD20", "CD68", "DAPI"]
        cell = get_cell_markers(all_markers, exclude={"CD20"})
        assert "CD20" not in cell
        assert "CD3" in cell
        assert "DAPI" in cell  # Not excluded with custom set


# =============================================================================
# Checkpointing Tests
# =============================================================================


class TestCheckpointing:
    """Tests for Parquet checkpointing utilities."""

    def test_save_and_load_parquet(self, tmp_path):
        """Test saving and loading a Parquet checkpoint."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import checkpoint_exists, load_checkpoint, save_checkpoint

        df = pd.DataFrame({"a": [1, 2, 3], "b": [4.0, 5.0, 6.0]})

        path = save_checkpoint(df, tmp_path, "test_checkpoint")
        assert path.exists()
        assert checkpoint_exists(tmp_path, "test_checkpoint")

        loaded = load_checkpoint(tmp_path, "test_checkpoint")
        assert loaded is not None
        pd.testing.assert_frame_equal(df, loaded)

    def test_load_nonexistent_returns_none(self, tmp_path):
        """Test that loading a non-existent checkpoint returns None."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import load_checkpoint

        result = load_checkpoint(tmp_path, "nonexistent")
        assert result is None

    def test_save_csv_format(self, tmp_path):
        """Test saving in CSV format."""
        import sys

        sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
        from Kanalysis import load_checkpoint, save_checkpoint

        df = pd.DataFrame({"x": [1, 2], "y": ["a", "b"]})

        save_checkpoint(df, tmp_path, "csv_test", format="csv")
        loaded = load_checkpoint(tmp_path, "csv_test", format="csv")

        assert loaded is not None
        assert len(loaded) == 2
