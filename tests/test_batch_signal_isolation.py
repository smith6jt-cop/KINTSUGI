"""
Tests for batch signal isolation processing.

Tests cover: channel discovery, method selection, per-channel processing,
tile smoothing, batch processing, QC visualization, and manifest I/O.
"""

import json

import numpy as np
import pytest
import tifffile
import yaml

from kintsugi.signal.batch import (
    BatchResult,
    ChannelResult,
    ChannelSpec,
    discover_channels,
    process_batch,
    process_channel,
    select_method,
    smooth_blank_for_subtraction,
)
from kintsugi.signal.isolation_qc import _self_normalize, generate_qc_pages


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def synthetic_project(tmp_path):
    """Create a minimal project with registered images and config.yaml.

    Structure:
    - cyc01: DAPI, Blank_1a, Blank_1b, Blank_1c
    - cyc02: DAPI_cyc02, CD31 (bright signal), FoxP3 (dim, blank >> signal), Empty_2c
    """
    project = tmp_path / "test_project"
    project.mkdir()

    registered = project / "data" / "processed" / "registered"
    cyc01 = registered / "cyc01"
    cyc02 = registered / "cyc02"
    cyc01.mkdir(parents=True)
    cyc02.mkdir(parents=True)

    rng = np.random.RandomState(42)
    shape = (256, 256)

    # Cycle 1: blanks
    blank_1a = rng.randint(2000, 8000, shape, dtype=np.uint16)
    blank_1b = rng.randint(8000, 20000, shape, dtype=np.uint16)  # High AF
    blank_1c = rng.randint(3000, 10000, shape, dtype=np.uint16)
    dapi_01 = rng.randint(5000, 30000, shape, dtype=np.uint16)

    tifffile.imwrite(str(cyc01 / "DAPI-01.tif"), dapi_01)
    tifffile.imwrite(str(cyc01 / "Blank_1a.tif"), blank_1a)
    tifffile.imwrite(str(cyc01 / "Blank_1b.tif"), blank_1b)
    tifffile.imwrite(str(cyc01 / "Blank_1c.tif"), blank_1c)

    # Cycle 2: CD31 (bright), FoxP3 (dim, blank >> signal)
    cd31 = rng.randint(5000, 40000, shape, dtype=np.uint16)
    cd31[100:150, 100:150] = 50000  # Bright signal region

    foxp3 = rng.randint(2000, 8000, shape, dtype=np.uint16)  # Dim: blank_1b dominates
    foxp3[120:130, 120:130] = 15000  # Small bright region

    dapi_cyc02 = rng.randint(5000, 30000, shape, dtype=np.uint16)
    empty_2c = np.zeros(shape, dtype=np.uint16)

    tifffile.imwrite(str(cyc02 / "DAPI_cyc02.tif"), dapi_cyc02)
    tifffile.imwrite(str(cyc02 / "CD31.tif"), cd31)
    tifffile.imwrite(str(cyc02 / "FoxP3.tif"), foxp3)
    tifffile.imwrite(str(cyc02 / "Empty_2c.tif"), empty_2c)

    # Config
    workflow_dir = project / "workflow"
    workflow_dir.mkdir(parents=True)
    config = {
        "project_dir": str(project),
        "channel_names": {
            1: ["DAPI-01", "Blank_1a", "Blank_1b", "Blank_1c"],
            2: ["DAPI_cyc02", "CD31", "FoxP3", "Empty_2c"],
        },
    }
    with open(workflow_dir / "config.yaml", "w") as f:
        yaml.dump(config, f)

    return project


@pytest.fixture
def bright_signal():
    """Signal image where signal clearly dominates blank."""
    rng = np.random.RandomState(10)
    signal = rng.randint(10000, 50000, (256, 256), dtype=np.uint16)
    signal[100:150, 100:150] = 60000
    return signal


@pytest.fixture
def dim_signal():
    """Dim signal image where blank dominates."""
    rng = np.random.RandomState(11)
    return rng.randint(2000, 8000, (256, 256), dtype=np.uint16)


@pytest.fixture
def high_blank():
    """High AF blank image."""
    rng = np.random.RandomState(12)
    return rng.randint(8000, 20000, (256, 256), dtype=np.uint16)


@pytest.fixture
def low_blank():
    """Low AF blank image."""
    rng = np.random.RandomState(13)
    return rng.randint(2000, 6000, (256, 256), dtype=np.uint16)


# =============================================================================
# TestDiscoverChannels
# =============================================================================


class TestDiscoverChannels:
    def test_finds_signal_markers(self, synthetic_project):
        specs = discover_channels(synthetic_project)
        marker_names = {s.marker_name for s in specs}
        assert "CD31" in marker_names
        assert "FoxP3" in marker_names

    def test_skips_dapi_blank_empty(self, synthetic_project):
        specs = discover_channels(synthetic_project)
        marker_names = {s.marker_name for s in specs}
        assert "DAPI-01" not in marker_names
        assert "DAPI_cyc02" not in marker_names
        assert "Blank_1a" not in marker_names
        assert "Empty_2c" not in marker_names

    def test_correct_blank_mapping(self, synthetic_project):
        specs = discover_channels(synthetic_project)
        by_name = {s.marker_name: s for s in specs}

        # CD31 is CH2 (idx 1) → Blank_1a
        assert by_name["CD31"].blank_name == "Blank_1a"
        # FoxP3 is CH3 (idx 2) → Blank_1b
        assert by_name["FoxP3"].blank_name == "Blank_1b"

    def test_paths_exist(self, synthetic_project):
        specs = discover_channels(synthetic_project)
        for spec in specs:
            assert spec.signal_path.exists()
            assert spec.blank_path.exists()

    def test_missing_config_raises(self, tmp_path):
        with pytest.raises(FileNotFoundError):
            discover_channels(tmp_path)

    def test_missing_channel_names_raises(self, tmp_path):
        workflow = tmp_path / "workflow"
        workflow.mkdir()
        with open(workflow / "config.yaml", "w") as f:
            yaml.dump({"project_dir": str(tmp_path)}, f)
        with pytest.raises(ValueError, match="No channel_names"):
            discover_channels(tmp_path)


# =============================================================================
# TestSelectMethod
# =============================================================================


class TestSelectMethod:
    def test_weighted_when_blank_dominates(self, dim_signal, high_blank):
        from kintsugi.signal.autofluorescence import analyze_for_subtraction

        analysis = analyze_for_subtraction(dim_signal, high_blank)
        method = select_method(dim_signal, high_blank, analysis)
        assert method == "weighted"

    def test_global_when_signal_dominates(self, bright_signal, low_blank):
        from kintsugi.signal.autofluorescence import analyze_for_subtraction

        analysis = analyze_for_subtraction(bright_signal, low_blank)
        method = select_method(bright_signal, low_blank, analysis)
        assert method == "global"


# =============================================================================
# TestProcessChannel
# =============================================================================


class TestProcessChannel:
    def test_produces_valid_uint16(self, synthetic_project, tmp_path):
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]
        output_dir = tmp_path / "output"
        result = process_channel(spec, method="global", output_dir=output_dir)

        assert result.status == "success"
        assert result.method == "global"
        # Verify output file is uint16
        out_path = output_dir / "CD31.tif"
        assert out_path.exists()
        img = tifffile.imread(str(out_path))
        assert img.dtype == np.uint16

    def test_quality_score_in_range(self, synthetic_project):
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]
        result = process_channel(spec, method="global")

        assert 0 <= result.quality_metrics["quality_score"] <= 1

    def test_weighted_fewer_zeros_for_dim(self, synthetic_project, tmp_path):
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "FoxP3"][0]

        # Process with global
        out_global = tmp_path / "global"
        result_global = process_channel(
            spec, method="global", output_dir=out_global, tile_smooth_sigma=0
        )

        # Process with weighted
        out_weighted = tmp_path / "weighted"
        result_weighted = process_channel(
            spec, method="weighted", output_dir=out_weighted, tile_smooth_sigma=0
        )

        # Weighted should preserve more signal (fewer zeros)
        assert result_weighted.zero_percent <= result_global.zero_percent

    def test_dry_run_no_output(self, synthetic_project, tmp_path):
        specs = discover_channels(synthetic_project)
        spec = specs[0]
        output_dir = tmp_path / "output"
        result = process_channel(spec, dry_run=True, output_dir=output_dir)

        assert result.status == "dry_run"
        assert not output_dir.exists()


# =============================================================================
# TestTileSmoothing
# =============================================================================


class TestTileSmoothing:
    def test_removes_tile_pattern(self):
        """Synthetic 5x5 tile grid with per-tile offsets → variation reduced significantly."""
        rng = np.random.RandomState(50)
        tile_h, tile_w = 200, 200
        grid = np.zeros((1000, 1000), dtype=np.uint16)

        # Create tiles with different base intensities
        offsets = rng.randint(8000, 15000, (5, 5))
        for r in range(5):
            for c in range(5):
                y0, x0 = r * tile_h, c * tile_w
                grid[y0:y0 + tile_h, x0:x0 + tile_w] = offsets[r, c]

        # Add noise
        grid = (grid.astype(np.float64) + rng.normal(0, 500, grid.shape)).clip(0, 65535).astype(
            np.uint16
        )

        # sigma ~ tile_pitch / 2 = 100
        smoothed = smooth_blank_for_subtraction(grid, sigma=100)

        # Measure tile-to-tile variation after smoothing (inner 3x3 tiles to avoid edges)
        tile_means = []
        for r in range(1, 4):
            for c in range(1, 4):
                y0, x0 = r * tile_h, c * tile_w
                margin = 40
                tile_center = smoothed[y0 + margin:y0 + tile_h - margin,
                                       x0 + margin:x0 + tile_w - margin]
                tile_means.append(np.mean(tile_center))

        tile_means = np.array(tile_means)
        cv = np.std(tile_means) / np.mean(tile_means)
        # Original CV is ~0.20; smoothing should reduce to < 0.10
        assert cv < 0.10, f"Tile-to-tile variation too high: CV={cv:.3f}"

    def test_preserves_gradient(self):
        """Linear AF gradient + tile noise → gradient preserved after smoothing."""
        # Create linear gradient
        gradient = np.linspace(5000, 25000, 300).astype(np.float64)
        image = np.tile(gradient, (300, 1))

        # Add tile noise
        rng = np.random.RandomState(51)
        tile_noise = np.zeros_like(image)
        for r in range(3):
            for c in range(3):
                offset = rng.randint(-2000, 2000)
                tile_noise[r * 100:(r + 1) * 100, c * 100:(c + 1) * 100] = offset

        image = np.clip(image + tile_noise, 0, 65535).astype(np.uint16)

        smoothed = smooth_blank_for_subtraction(image, sigma=50)

        # Check gradient preservation via R² of column means vs linear fit
        col_means = np.mean(smoothed.astype(np.float64), axis=0)
        x = np.arange(len(col_means))
        coeffs = np.polyfit(x, col_means, 1)
        predicted = np.polyval(coeffs, x)
        ss_res = np.sum((col_means - predicted) ** 2)
        ss_tot = np.sum((col_means - np.mean(col_means)) ** 2)
        r_squared = 1 - ss_res / ss_tot

        assert r_squared > 0.95, f"Gradient not preserved: R²={r_squared:.4f}"

    def test_preserves_mean(self):
        """Output mean within 1% of input mean when preserve_mean=True."""
        rng = np.random.RandomState(52)
        image = rng.randint(5000, 20000, (256, 256), dtype=np.uint16)

        smoothed = smooth_blank_for_subtraction(image, sigma=50, preserve_mean=True)

        orig_mean = np.mean(image[image > 0].astype(np.float64))
        smooth_mean = np.mean(smoothed[smoothed > 0].astype(np.float64))

        rel_diff = abs(smooth_mean - orig_mean) / orig_mean
        assert rel_diff < 0.01, f"Mean not preserved: rel_diff={rel_diff:.4f}"

    def test_sigma_zero_is_noop(self):
        """sigma=0 returns original image unchanged."""
        rng = np.random.RandomState(53)
        image = rng.randint(5000, 20000, (256, 256), dtype=np.uint16)

        result = smooth_blank_for_subtraction(image, sigma=0)
        np.testing.assert_array_equal(result, image)


# =============================================================================
# TestProcessBatch
# =============================================================================


class TestProcessBatch:
    def test_dry_run_writes_no_files(self, synthetic_project):
        output_dir = synthetic_project / "data" / "processed" / "signal_isolated"
        result = process_batch(synthetic_project, dry_run=True)

        assert not output_dir.exists() or not list(output_dir.glob("*.tif"))
        assert result.summary["total"] > 0

    def test_manifest_written(self, synthetic_project):
        result = process_batch(synthetic_project, force=True)
        manifest_path = (
            synthetic_project / "data" / "processed" / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        assert manifest_path.exists()

        with open(manifest_path) as f:
            data = json.load(f)
        assert "channels" in data
        assert "summary" in data
        assert len(data["channels"]) == result.summary["total"]

    def test_skip_existing(self, synthetic_project):
        # First run
        process_batch(synthetic_project, force=True)
        # Second run without force
        result = process_batch(synthetic_project, force=False)
        assert result.summary["skipped"] > 0

    def test_channel_filter(self, synthetic_project):
        result = process_batch(
            synthetic_project, channels=["CD31"], force=True
        )
        assert result.summary["total"] == 1
        assert "CD31" in result.channels


# =============================================================================
# TestQCVisualization
# =============================================================================


class TestQCVisualization:
    def test_self_normalize_range(self):
        image = np.random.randint(1000, 50000, (100, 100), dtype=np.uint16)
        normalized = _self_normalize(image)
        assert normalized.min() >= 0
        assert normalized.max() <= 1

    def test_self_normalize_visible(self):
        """Normalized image has visible content (mean > 0.1)."""
        image = np.random.randint(5000, 30000, (100, 100), dtype=np.uint16)
        normalized = _self_normalize(image)
        assert normalized.mean() > 0.1

    def test_generates_pages(self, synthetic_project):
        """Process then generate QC pages."""
        process_batch(synthetic_project, force=True)
        pages = generate_qc_pages(synthetic_project, page_size=6, downsample=1)
        assert len(pages) >= 1
        for page in pages:
            assert page.exists()
            assert page.suffix == ".png"


# =============================================================================
# TestManifest
# =============================================================================


class TestManifest:
    def test_manifest_valid_json(self, synthetic_project):
        process_batch(synthetic_project, force=True)
        manifest_path = (
            synthetic_project / "data" / "processed" / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        with open(manifest_path) as f:
            data = json.load(f)

        # Top-level fields
        assert "timestamp" in data
        assert "method_override" in data
        assert "tile_smooth_sigma" in data
        assert "summary" in data
        assert "channels" in data

    def test_quality_metrics_present(self, synthetic_project):
        process_batch(synthetic_project, force=True)
        manifest_path = (
            synthetic_project / "data" / "processed" / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        with open(manifest_path) as f:
            data = json.load(f)

        for marker, info in data["channels"].items():
            assert "quality_metrics" in info
            assert "quality_score" in info["quality_metrics"]
            assert "parameters" in info
            assert "zero_percent" in info

    def test_batch_result_roundtrip(self, synthetic_project):
        result = process_batch(synthetic_project, force=True)
        d = result.to_dict()
        # Verify serializable
        json_str = json.dumps(d, default=str)
        loaded = json.loads(json_str)
        assert loaded["summary"]["total"] == result.summary["total"]


# =============================================================================
# TestDataclasses
# =============================================================================


class TestDataclasses:
    def test_channel_spec_to_dict(self, tmp_path):
        spec = ChannelSpec(
            marker_name="CD31",
            cycle=2,
            channel_position=2,
            signal_path=tmp_path / "signal.tif",
            blank_name="Blank_1a",
            blank_path=tmp_path / "blank.tif",
        )
        d = spec.to_dict()
        assert d["marker_name"] == "CD31"
        assert isinstance(d["signal_path"], str)

    def test_channel_result_to_dict(self):
        result = ChannelResult(
            marker_name="CD31",
            cycle=2,
            blank_used="Blank_1a",
            method="global",
            parameters={"blank_scale_factor": 0.85},
            quality_metrics={"quality_score": 0.72},
            zero_percent=12.3,
        )
        d = result.to_dict()
        assert d["method"] == "global"
        assert d["zero_percent"] == 12.3

    def test_batch_result_to_dict(self):
        batch = BatchResult(
            timestamp="2026-02-20T12:00:00",
            method_override="auto",
            summary={"total": 5},
        )
        d = batch.to_dict()
        assert d["timestamp"] == "2026-02-20T12:00:00"
        assert d["summary"]["total"] == 5
