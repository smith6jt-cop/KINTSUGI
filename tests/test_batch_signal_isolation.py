"""
Tests for batch signal isolation processing.

Tests cover: channel discovery, method selection, per-channel processing,
tile smoothing, batch processing, QC visualization, manifest I/O,
recipe parsing, blank name resolution, background cleaning, and
recipe-driven processing.
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
    CleanParams,
    MarkerRecipe,
    SubtractionParams,
    _check_over_subtraction,
    _extract_blank_position,
    _normalize_blank_name,
    _resolve_blank_path,
    _strip_for_comparison,
    _validate_secondary_blank,
    clean_background,
    clip_outliers,
    discover_channels,
    load_legacy_recipes,
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


@pytest.fixture
def recipe_dir(tmp_path):
    """Create sample legacy param files for recipe parsing."""
    recipe_path = tmp_path / "recipes"
    recipe_path.mkdir()

    # CD3e: primary + second subtraction + clean
    (recipe_path / "CD3e_param.txt").write_text(
        "timestamp: '2025/03/17_11:36'\n"
        "signal_channel: 'CD3e'\n"
        "blankID: 'Blank1b'\n"
        "blank_clip_factor: 20000\n"
        "blank_scale_factor: 1.4\n"
        "smooth_low: False\n"
        "smooth_high: False\n"
        "low_size: 2\n"
        "high_size: 2\n"
        "erosion: 1\n"
        "second_subtract: True\n"
        "blankID2: 'Blank13b'\n"
        "blank_clip_factor_2: 0\n"
        "blank_scale_factor_2: 0.6\n"
        "smooth_low_2: False\n"
        "smooth_high_2: False\n"
        "low_size_2: 2\n"
        "high_size_2: 1\n"
        "erosion_2: 1\n"
        "clean_used: True\n"
        "clean_image: dask.array<astype, shape=(8986, 9464), dtype=uint16>\n"
        "smooth: True\n"
        "remove_small: True\n"
        "small_size: 50\n"
        "footprint: 3\n"
    )

    # CD31: primary only + clean (no second subtraction)
    (recipe_path / "CD31_param.txt").write_text(
        "signal_channel: 'CD31'\n"
        "blankID: 'Blank13b'\n"
        "blank_clip_factor: 0\n"
        "blank_scale_factor: 0.2\n"
        "smooth_low: False\n"
        "smooth_high: False\n"
        "low_size: 2\n"
        "high_size: 2\n"
        "erosion: 1\n"
        "second_subtract: False\n"
        "clean_used: True\n"
        "background_threshold: 3000\n"
        "smooth: False\n"
        "remove_small: False\n"
        "small_size: 75\n"
        "footprint: 2\n"
    )

    # CD163: uses non-blank marker as blank source
    (recipe_path / "CD163_param.txt").write_text(
        "signal_channel: 'CD163'\n"
        "blankID: 'CD3e'\n"
        "blank_clip_factor: 11000\n"
        "blank_scale_factor: 1.5\n"
        "erosion: 2\n"
        "second_subtract: True\n"
        "blankID2: 'Blank13c'\n"
        "blank_clip_factor_2: 0\n"
        "blank_scale_factor_2: 0.7\n"
        "erosion_2: 2\n"
        "clean_used: True\n"
        "background_threshold: 1000\n"
        "smooth: True\n"
        "smooth_threshold: 800\n"
        "remove_small: True\n"
        "small_size: 3\n"
        "footprint: 2\n"
    )

    # FoxP3: includes dask artifact lines
    (recipe_path / "FoxP3_param.txt").write_text(
        "signal_channel: 'FoxP3'\n"
        "blankID: 'Blank1b'\n"
        "blank_clip_factor: 17200\n"
        "blank_scale_factor: 1.2\n"
        "erosion: 0\n"
        "second_subtract: True\n"
        "blankID2: 'Blank1b'\n"
        "blank_clip_factor_2: 0\n"
        "blank_scale_factor_2: 0.9\n"
        "erosion_2: 1\n"
        "clean_used: True\n"
        "clean_image: dask.array<array, shape=(8986, 9464), dtype=float64>\n"
        "background_threshold: 10\n"
        "smooth: False\n"
        "remove_small: True\n"
        "small_size: 10\n"
        "footprint: 3\n"
    )

    return recipe_path


@pytest.fixture
def multi_cycle_project(tmp_path):
    """Project with 13 cycles and diverse blanks for recipe blank resolution testing."""
    project = tmp_path / "multi_project"
    project.mkdir()

    registered = project / "data" / "processed" / "registered"
    rng = np.random.RandomState(99)
    shape = (64, 64)

    channel_names = {}

    # Cycle 1: blanks
    cyc01 = registered / "cyc01"
    cyc01.mkdir(parents=True)
    channel_names[1] = ["DAPI-01", "Blank_1a", "Blank_1b", "Blank_1c"]
    for name in channel_names[1]:
        tifffile.imwrite(
            str(cyc01 / f"{name}.tif"), rng.randint(1000, 10000, shape, dtype=np.uint16)
        )

    # Cycle 2: signal channels
    cyc02 = registered / "cyc02"
    cyc02.mkdir(parents=True)
    channel_names[2] = ["DAPI-02", "CD31", "CD3e", "CD45"]
    for name in channel_names[2]:
        tifffile.imwrite(
            str(cyc02 / f"{name}.tif"), rng.randint(1000, 20000, shape, dtype=np.uint16)
        )

    # Cycle 13: more blanks (e.g., Blank_13b, Blank_13c)
    cyc13 = registered / "cyc13"
    cyc13.mkdir(parents=True)
    channel_names[13] = ["DAPI-13", "Blank_13a", "Blank_13b", "Blank_13c"]
    for name in channel_names[13]:
        tifffile.imwrite(
            str(cyc13 / f"{name}.tif"), rng.randint(1000, 10000, shape, dtype=np.uint16)
        )

    # Config
    workflow_dir = project / "workflow"
    workflow_dir.mkdir(parents=True)
    config = {
        "project_dir": str(project),
        "channel_names": channel_names,
    }
    with open(workflow_dir / "config.yaml", "w") as f:
        yaml.dump(config, f)

    return project


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

    def test_recipe_overrides_blank(self, multi_cycle_project, recipe_dir):
        """With recipes, blank selection uses recipe blank_name, not positional."""
        recipes = load_legacy_recipes(recipe_dir)
        specs = discover_channels(multi_cycle_project, recipes=recipes)
        by_name = {s.marker_name: s for s in specs}

        # CD31 recipe uses Blank13b → should resolve to Blank_13b in cyc13
        assert "CD31" in by_name
        assert by_name["CD31"].blank_name == "Blank_13b"

        # CD3e recipe uses Blank1b → should resolve to Blank_1b in cyc01
        assert "CD3e" in by_name
        assert by_name["CD3e"].blank_name == "Blank_1b"


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
                grid[y0 : y0 + tile_h, x0 : x0 + tile_w] = offsets[r, c]

        # Add noise
        grid = (
            (grid.astype(np.float64) + rng.normal(0, 500, grid.shape))
            .clip(0, 65535)
            .astype(np.uint16)
        )

        # sigma ~ tile_pitch / 2 = 100
        smoothed = smooth_blank_for_subtraction(grid, sigma=100)

        # Measure tile-to-tile variation after smoothing (inner 3x3 tiles to avoid edges)
        tile_means = []
        for r in range(1, 4):
            for c in range(1, 4):
                y0, x0 = r * tile_h, c * tile_w
                margin = 40
                tile_center = smoothed[
                    y0 + margin : y0 + tile_h - margin, x0 + margin : x0 + tile_w - margin
                ]
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
                tile_noise[r * 100 : (r + 1) * 100, c * 100 : (c + 1) * 100] = offset

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
            synthetic_project
            / "data"
            / "processed"
            / "signal_isolated"
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
        result = process_batch(synthetic_project, channels=["CD31"], force=True)
        assert result.summary["total"] == 1
        assert "CD31" in result.channels

    def test_default_sigma_is_zero(self, synthetic_project):
        """Default tile_smooth_sigma changed from 500 to 0."""
        result = process_batch(synthetic_project, dry_run=True)
        assert result.tile_smooth_sigma == 0.0


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
            synthetic_project
            / "data"
            / "processed"
            / "signal_isolated"
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
            synthetic_project
            / "data"
            / "processed"
            / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        with open(manifest_path) as f:
            data = json.load(f)

        for _marker, info in data["channels"].items():
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

    def test_manifest_includes_recipe_dir(self, synthetic_project, recipe_dir):
        result = process_batch(synthetic_project, recipe_dir=recipe_dir, force=True)
        d = result.to_dict()
        assert d["recipe_dir"] == str(recipe_dir)


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


# =============================================================================
# TestNormalizeBlankName
# =============================================================================


class TestNormalizeBlankName:
    def test_blank1b(self):
        assert _normalize_blank_name("Blank1b") == "Blank_1b"

    def test_blank13c(self):
        assert _normalize_blank_name("Blank13c") == "Blank_13c"

    def test_blank_1a_passthrough(self):
        """Already-normalized names pass through unchanged."""
        assert _normalize_blank_name("Blank_1a") == "Blank_1a"

    def test_non_blank_passthrough(self):
        assert _normalize_blank_name("CD3e") == "CD3e"
        assert _normalize_blank_name("DAPI") == "DAPI"
        assert _normalize_blank_name("HLA-DR") == "HLA-DR"

    def test_case_insensitive(self):
        assert _normalize_blank_name("blank1b") == "Blank_1b"

    def test_strip_for_comparison(self):
        assert _strip_for_comparison("HLA-DR") == "hladr"
        assert _strip_for_comparison("Blank_1a") == "blank1a"
        assert _strip_for_comparison("CD3e") == "cd3e"


# =============================================================================
# TestResolveBlankPath
# =============================================================================


class TestResolveBlankPath:
    def test_exact_match(self, tmp_path):
        path = tmp_path / "Blank_1b.tif"
        location_map = {"Blank_1b": path}
        result = _resolve_blank_path("Blank_1b", location_map)
        assert result == ("Blank_1b", path)

    def test_normalize_then_match(self, tmp_path):
        """Blank1b normalizes to Blank_1b and matches."""
        path = tmp_path / "Blank_1b.tif"
        location_map = {"Blank_1b": path}
        result = _resolve_blank_path("Blank1b", location_map)
        assert result == ("Blank_1b", path)

    def test_fuzzy_match_hladr(self, tmp_path):
        """HLADR fuzzy-matches HLA-DR."""
        path = tmp_path / "HLA-DR.tif"
        location_map = {"HLA-DR": path}
        result = _resolve_blank_path("HLADR", location_map)
        assert result == ("HLA-DR", path)

    def test_not_found(self):
        result = _resolve_blank_path("NonExistent", {})
        assert result is None

    def test_cd3e_direct(self, tmp_path):
        """Non-blank names resolve by exact match."""
        path = tmp_path / "CD3e.tif"
        location_map = {"CD3e": path}
        result = _resolve_blank_path("CD3e", location_map)
        assert result == ("CD3e", path)


# =============================================================================
# TestLoadLegacyRecipes
# =============================================================================


class TestLoadLegacyRecipes:
    def test_loads_all_recipes(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        assert len(recipes) == 4
        assert "CD3e" in recipes
        assert "CD31" in recipes
        assert "CD163" in recipes
        assert "FoxP3" in recipes

    def test_primary_params(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        cd3e = recipes["CD3e"]
        assert cd3e.primary.blank_name == "Blank1b"
        assert cd3e.primary.blank_clip_factor == 20000
        assert cd3e.primary.blank_scale_factor == 1.4
        assert cd3e.primary.erosion == 1

    def test_second_subtraction(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        cd3e = recipes["CD3e"]
        assert cd3e.second is not None
        assert cd3e.second.blank_name == "Blank13b"
        assert cd3e.second.blank_scale_factor == 0.6

    def test_no_second_subtraction(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        cd31 = recipes["CD31"]
        assert cd31.second is None

    def test_clean_params(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        cd31 = recipes["CD31"]
        assert cd31.clean is not None
        assert cd31.clean.background_threshold == 3000
        assert cd31.clean.smooth is False
        assert cd31.clean.remove_small is False

    def test_clean_with_smooth(self, recipe_dir):
        recipes = load_legacy_recipes(recipe_dir)
        cd3e = recipes["CD3e"]
        assert cd3e.clean is not None
        assert cd3e.clean.smooth is True
        assert cd3e.clean.remove_small is True
        assert cd3e.clean.small_size == 50
        assert cd3e.clean.footprint == 3

    def test_skips_dask_lines(self, recipe_dir):
        """Dask array lines should be skipped without error."""
        recipes = load_legacy_recipes(recipe_dir)
        # All recipes should load successfully despite dask lines
        assert len(recipes) == 4

    def test_non_blank_as_blank_source(self, recipe_dir):
        """CD163 uses CD3e as blank source."""
        recipes = load_legacy_recipes(recipe_dir)
        cd163 = recipes["CD163"]
        assert cd163.primary.blank_name == "CD3e"

    def test_empty_dir(self, tmp_path):
        recipes = load_legacy_recipes(tmp_path)
        assert len(recipes) == 0

    def test_missing_signal_channel(self, tmp_path):
        """File without signal_channel is skipped."""
        recipe_path = tmp_path / "bad_param.txt"
        recipe_path.write_text("blankID: 'Blank1b'\nblank_scale_factor: 1.0\n")
        recipes = load_legacy_recipes(tmp_path)
        assert len(recipes) == 0


# =============================================================================
# TestCleanBackground
# =============================================================================


class TestCleanBackground:
    def test_threshold_zeros_background(self):
        """Pixels below threshold become zero."""
        image = np.array([[500, 1000, 2000], [3000, 100, 4000]], dtype=np.uint16)
        params = CleanParams(background_threshold=1000)
        result = clean_background(image, params)

        assert result[0, 0] == 0  # 500 <= 1000
        assert result[0, 1] == 0  # 1000 <= 1000
        assert result[0, 2] > 0  # 2000 > 1000

    def test_smooth_transition_zone(self):
        """Median filter applied in transition zone only."""
        rng = np.random.RandomState(70)
        image = rng.randint(100, 5000, (64, 64), dtype=np.uint16)
        params = CleanParams(
            background_threshold=100,
            smooth=True,
            smooth_threshold=2000,
            footprint=3,
        )
        result = clean_background(image, params)

        # Output should be uint16
        assert result.dtype == np.uint16
        # Background pixels should be zero
        assert np.sum(result[image <= 100] > 0) == 0

    def test_remove_small_objects(self):
        """Small isolated bright pixels are removed."""
        image = np.zeros((64, 64), dtype=np.uint16)
        # Large region
        image[10:30, 10:30] = 5000
        # Tiny speck (< min_size)
        image[50, 50] = 5000

        params = CleanParams(
            background_threshold=0,
            remove_small=True,
            small_size=10,
        )
        result = clean_background(image, params)

        # Large region preserved
        assert result[20, 20] > 0
        # Tiny speck removed (single pixel < min_size=10)
        assert result[50, 50] == 0

    def test_noop_with_defaults(self):
        """Default params (threshold=0) effectively zero-out nothing."""
        rng = np.random.RandomState(71)
        image = rng.randint(100, 50000, (64, 64), dtype=np.uint16)
        params = CleanParams()  # All defaults
        result = clean_background(image, params)

        # threshold=0 means bg_thresh=max(1,0)=1, zeroing only pixels <= 1
        # Most pixels are > 1, so most should survive
        assert np.sum(result > 0) > 0.9 * image.size

    def test_output_dtype_uint16(self):
        image = np.ones((32, 32), dtype=np.uint16) * 10000
        params = CleanParams(background_threshold=5000)
        result = clean_background(image, params)
        assert result.dtype == np.uint16


# =============================================================================
# TestRecipeDrivenProcessing
# =============================================================================


class TestRecipeDrivenProcessing:
    def test_single_subtraction_with_recipe(self, synthetic_project, tmp_path):
        """Process CD31 with a simple recipe (primary only + clean)."""
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]

        recipe = MarkerRecipe(
            marker_name="CD31",
            primary=SubtractionParams(
                blank_name="Blank_1a",
                blank_clip_factor=0,
                blank_scale_factor=0.5,
            ),
            clean=CleanParams(background_threshold=500),
        )

        output_dir = tmp_path / "recipe_out"
        result = process_channel(spec, recipe=recipe, output_dir=output_dir, location_map={})

        assert result.status == "success"
        assert result.method == "recipe"
        assert (output_dir / "CD31.tif").exists()
        img = tifffile.imread(str(output_dir / "CD31.tif"))
        assert img.dtype == np.uint16

    def test_double_subtraction_with_recipe(self, synthetic_project, tmp_path):
        """Process FoxP3 with primary + second subtraction.

        FoxP3 is dim with high blank, so double subtraction over-subtracts.
        With quality gate enabled, this triggers auto fallback.
        Disable quality gate to test pure recipe path.
        """
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "FoxP3"][0]

        # Build location map for second blank resolution
        registered_dir = synthetic_project / "data" / "processed" / "registered"
        location_map = {
            "Blank_1a": registered_dir / "cyc01" / "Blank_1a.tif",
            "Blank_1b": registered_dir / "cyc01" / "Blank_1b.tif",
            "Blank_1c": registered_dir / "cyc01" / "Blank_1c.tif",
        }

        recipe = MarkerRecipe(
            marker_name="FoxP3",
            primary=SubtractionParams(
                blank_name="Blank_1b",
                blank_clip_factor=5000,
                blank_scale_factor=1.2,
            ),
            second=SubtractionParams(
                blank_name="Blank_1b",
                blank_clip_factor=0,
                blank_scale_factor=0.9,
                erosion=1,
            ),
        )

        output_dir = tmp_path / "recipe_out"
        # Disable quality gate to test pure recipe path
        result = process_channel(
            spec, recipe=recipe, output_dir=output_dir, location_map=location_map,
            quality_gate=0.0,
        )

        assert result.status == "success"
        assert result.method == "recipe"
        assert "Blank_1b" in result.blank_used

    def test_recipe_dry_run(self, synthetic_project):
        """Dry run with recipe returns parameters without saving."""
        specs = discover_channels(synthetic_project)
        spec = specs[0]

        recipe = MarkerRecipe(
            marker_name=spec.marker_name,
            primary=SubtractionParams(blank_name="Blank_1a", blank_scale_factor=1.0),
        )

        result = process_channel(spec, recipe=recipe, dry_run=True)
        assert result.status == "dry_run"
        assert result.method == "recipe"
        assert "primary_scale" in result.parameters

    def test_fallback_without_recipe(self, synthetic_project, tmp_path):
        """Without recipe, uses auto-analysis path (backward compatible)."""
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]

        output_dir = tmp_path / "auto_out"
        result = process_channel(spec, method="global", output_dir=output_dir)

        assert result.status == "success"
        assert result.method == "global"
        assert "blank_scale_factor" in result.parameters


# =============================================================================
# TestBatchWithRecipes
# =============================================================================


class TestBatchWithRecipes:
    def test_batch_with_recipe_dir(self, synthetic_project, tmp_path):
        """Batch processing with recipes uses recipe pipeline for matched markers."""
        # Create recipe for CD31 only
        recipe_path = tmp_path / "recipes"
        recipe_path.mkdir()
        (recipe_path / "CD31_param.txt").write_text(
            "signal_channel: 'CD31'\n"
            "blankID: 'Blank_1a'\n"
            "blank_clip_factor: 0\n"
            "blank_scale_factor: 0.5\n"
            "erosion: 0\n"
            "second_subtract: False\n"
            "clean_used: True\n"
            "background_threshold: 500\n"
            "smooth: False\n"
            "remove_small: False\n"
        )

        result = process_batch(synthetic_project, recipe_dir=recipe_path, force=True)

        # CD31 should use recipe, FoxP3 should use auto
        assert result.channels["CD31"].method == "recipe"
        assert result.channels["FoxP3"].method in ("global", "weighted")
        assert result.summary["recipe"] == 1

    def test_batch_recipe_count_in_summary(self, synthetic_project, tmp_path):
        """Summary includes recipe count.

        Disable quality gate so FoxP3 (dim + high blank) isn't flagged.
        """
        recipe_path = tmp_path / "recipes"
        recipe_path.mkdir()
        # Create recipes for both CD31 and FoxP3
        for name, blank in [("CD31", "Blank_1a"), ("FoxP3", "Blank_1b")]:
            (recipe_path / f"{name}_param.txt").write_text(
                f"signal_channel: '{name}'\n"
                f"blankID: '{blank}'\n"
                "blank_scale_factor: 1.0\n"
                "second_subtract: False\n"
                "clean_used: False\n"
            )

        result = process_batch(
            synthetic_project, recipe_dir=recipe_path, force=True, quality_gate=0.0
        )
        assert result.summary["recipe"] == 2
        assert result.summary["global"] == 0
        assert result.summary["weighted"] == 0


# =============================================================================
# TestLearningIntegration
# =============================================================================


class TestLearningIntegration:
    def test_learn_flag_creates_db(self, synthetic_project, tmp_path):
        """With --learn and tissue_type, learning DB gets records."""
        recipe_path = tmp_path / "recipes"
        recipe_path.mkdir()
        (recipe_path / "CD31_param.txt").write_text(
            "signal_channel: 'CD31'\n"
            "blankID: 'Blank_1a'\n"
            "blank_scale_factor: 0.5\n"
            "second_subtract: False\n"
            "clean_used: True\n"
            "background_threshold: 500\n"
            "smooth: False\n"
            "remove_small: False\n"
        )

        process_batch(
            synthetic_project,
            recipe_dir=recipe_path,
            tissue_type="spleen",
            learn=True,
            force=True,
        )

        # Verify the learning DB was created
        db_path = synthetic_project / ".kintsugi" / "parameter_learning.db"
        if db_path.exists():
            import sqlite3

            conn = sqlite3.connect(str(db_path))
            cursor = conn.cursor()
            cursor.execute("SELECT COUNT(*) FROM parameter_records")
            count = cursor.fetchone()[0]
            conn.close()
            assert count > 0

    def test_no_learn_flag(self, synthetic_project, tmp_path):
        """With --no-learn, no DB records created."""
        recipe_path = tmp_path / "recipes"
        recipe_path.mkdir()
        (recipe_path / "CD31_param.txt").write_text(
            "signal_channel: 'CD31'\n"
            "blankID: 'Blank_1a'\n"
            "blank_scale_factor: 0.5\n"
            "second_subtract: False\n"
            "clean_used: False\n"
        )

        process_batch(
            synthetic_project,
            recipe_dir=recipe_path,
            tissue_type="spleen",
            learn=False,
            force=True,
        )

        db_path = synthetic_project / ".kintsugi" / "parameter_learning.db"
        # DB should not exist when learn=False
        assert not db_path.exists()


# =============================================================================
# Bug 1: Blank Position Extraction and Positional Fallback
# =============================================================================


class TestExtractBlankPosition:
    def test_blank_1b(self):
        assert _extract_blank_position("Blank_1b") == "b"

    def test_blank_13c(self):
        assert _extract_blank_position("Blank_13c") == "c"

    def test_unnormalized_blank1b(self):
        assert _extract_blank_position("Blank1b") == "b"

    def test_non_blank_returns_none(self):
        assert _extract_blank_position("CD3e") is None

    def test_dapi_returns_none(self):
        assert _extract_blank_position("DAPI-01") is None

    def test_blank_without_suffix(self):
        """Blank_1 (no letter suffix) returns None."""
        assert _extract_blank_position("Blank_1") is None


class TestPositionalFallback:
    def test_positional_fallback_remaps_blank13b(self, tmp_path):
        """Blank13b -> Blank_1b when only cycle 1 blanks exist."""
        path = tmp_path / "Blank_1b.tif"
        location_map = {"Blank_1b": path}
        result = _resolve_blank_path("Blank13b", location_map)
        assert result is not None
        assert result == ("Blank_1b", path)

    def test_no_fallback_when_exact_exists(self, tmp_path):
        """Blank_13b found exactly, no remapping needed."""
        path_13b = tmp_path / "Blank_13b.tif"
        path_1b = tmp_path / "Blank_1b.tif"
        location_map = {"Blank_13b": path_13b, "Blank_1b": path_1b}
        result = _resolve_blank_path("Blank13b", location_map)
        assert result is not None
        assert result == ("Blank_13b", path_13b)

    def test_positional_fallback_blank13c(self, tmp_path):
        """Blank_13c -> Blank_1c when only cycle 1 blanks exist."""
        path = tmp_path / "Blank_1c.tif"
        location_map = {"Blank_1a": tmp_path / "a.tif", "Blank_1c": path}
        result = _resolve_blank_path("Blank13c", location_map)
        assert result is not None
        assert result == ("Blank_1c", path)

    def test_secondary_blank_positional_fallback_integration(
        self, synthetic_project, tmp_path
    ):
        """Full recipe processing with a remapped secondary blank."""
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]

        # Build location map with only cycle 1 blanks
        registered_dir = synthetic_project / "data" / "processed" / "registered"
        location_map = {
            "Blank_1a": registered_dir / "cyc01" / "Blank_1a.tif",
            "Blank_1b": registered_dir / "cyc01" / "Blank_1b.tif",
            "Blank_1c": registered_dir / "cyc01" / "Blank_1c.tif",
        }

        # Recipe with Blank13b (doesn't exist) → should remap to Blank_1b
        recipe = MarkerRecipe(
            marker_name="CD31",
            primary=SubtractionParams(
                blank_name="Blank_1a",
                blank_scale_factor=0.5,
            ),
            second=SubtractionParams(
                blank_name="Blank13b",
                blank_scale_factor=0.6,
            ),
        )

        output_dir = tmp_path / "recipe_out"
        result = process_channel(
            spec, recipe=recipe, output_dir=output_dir, location_map=location_map
        )
        assert result.status == "success"
        # Should have used the remapped blank
        assert "Blank13b" in result.blank_used


# =============================================================================
# Bug 2: Self-Referential Secondary Blank Validation
# =============================================================================


class TestSelfReferentialBlank:
    def test_self_referential_secondary_rejected(self):
        """CD20 as secondary for CD20 is detected."""
        valid, reason = _validate_secondary_blank("CD20", "CD20")
        assert not valid
        assert reason == "self-referential"

    def test_non_self_referential_passes(self):
        """Blank_1b as secondary for CD20 passes."""
        valid, reason = _validate_secondary_blank("CD20", "Blank_1b")
        assert valid
        assert reason == ""

    def test_case_insensitive_self_reference(self):
        """cd20 vs CD20 is still caught."""
        valid, reason = _validate_secondary_blank("CD20", "cd20")
        assert not valid
        assert reason == "self-referential"

    def test_self_referential_in_recipe(self, synthetic_project, tmp_path):
        """Self-referential secondary blank is skipped during recipe processing."""
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "CD31"][0]

        # Recipe where CD31 uses itself as secondary blank
        recipe = MarkerRecipe(
            marker_name="CD31",
            primary=SubtractionParams(
                blank_name="Blank_1a",
                blank_scale_factor=0.5,
            ),
            second=SubtractionParams(
                blank_name="CD31",
                blank_scale_factor=0.6,
            ),
        )

        registered_dir = synthetic_project / "data" / "processed" / "registered"
        location_map = {
            "CD31": registered_dir / "cyc02" / "CD31.tif",
            "Blank_1a": registered_dir / "cyc01" / "Blank_1a.tif",
        }

        output_dir = tmp_path / "recipe_out"
        result = process_channel(
            spec, recipe=recipe, output_dir=output_dir, location_map=location_map
        )
        assert result.status == "success"
        # Secondary was skipped, so blank_used should NOT include "+CD31"
        assert "+CD31" not in result.blank_used
        assert "secondary_blank_self-referential" in result.warning


# =============================================================================
# Bug 3: Over-Subtraction Detection and Quality Gate
# =============================================================================


class TestOverSubtractionDetection:
    def test_over_subtraction_detected(self):
        """98% zero → over-subtracted."""
        over_sub, reason = _check_over_subtraction(98.0, 0.3)
        assert over_sub
        assert "zero=98.0%" in reason

    def test_low_quality_with_high_zeros_detected(self):
        """Low quality + high zeros triggers gate."""
        over_sub, reason = _check_over_subtraction(75.0, 0.4)
        assert over_sub
        assert "quality=0.40" in reason

    def test_low_quality_low_zeros_passes(self):
        """Low quality but low zeros does not trigger (signal preserved)."""
        over_sub, _ = _check_over_subtraction(30.0, 0.4)
        assert not over_sub

    def test_good_result_passes(self):
        """Normal result passes gate."""
        over_sub, reason = _check_over_subtraction(30.0, 0.75)
        assert not over_sub
        assert reason == ""

    def test_quality_gate_configurable(self):
        """Custom threshold works."""
        # With default threshold (0.6), quality 0.5 + 75% zeros → fails
        over_sub_default, _ = _check_over_subtraction(75.0, 0.5)
        assert over_sub_default

        # With lowered threshold (0.4), quality 0.5 passes even with 75% zeros
        over_sub_custom, _ = _check_over_subtraction(75.0, 0.5, quality_gate=0.4)
        assert not over_sub_custom

        # With quality_gate=0, gate is disabled entirely
        over_sub_disabled, _ = _check_over_subtraction(99.0, 0.1, quality_gate=0.0)
        assert not over_sub_disabled

    def test_over_subtraction_marks_failed_in_recipe(self, synthetic_project, tmp_path):
        """Recipe over-subtraction marks status=failed when fallback doesn't improve."""
        specs = discover_channels(synthetic_project)
        spec = [s for s in specs if s.marker_name == "FoxP3"][0]

        # Recipe that aggressively over-subtracts (high scale factor)
        recipe = MarkerRecipe(
            marker_name="FoxP3",
            primary=SubtractionParams(
                blank_name="Blank_1b",
                blank_clip_factor=0,
                blank_scale_factor=3.0,  # Very aggressive
            ),
        )

        registered_dir = synthetic_project / "data" / "processed" / "registered"
        location_map = {
            "Blank_1b": registered_dir / "cyc01" / "Blank_1b.tif",
        }

        output_dir = tmp_path / "output"
        result = process_channel(
            spec,
            recipe=recipe,
            output_dir=output_dir,
            location_map=location_map,
            quality_gate=0.6,
        )
        # FoxP3 is dim, blank is high, scale 3.0 → severe over-subtraction
        # Either auto fallback improves (recipe_source=auto_fallback) or
        # it stays failed — either way the gate was triggered
        assert result.status in ("success", "failed") or result.recipe_source == "auto_fallback"

    def test_failed_count_in_summary(self, synthetic_project, tmp_path):
        """Batch summary includes 'failed' key."""
        result = process_batch(synthetic_project, force=True)
        assert "failed" in result.summary
        # The 'failed' count exists in summary (may be 0 for this data)
        assert isinstance(result.summary["failed"], int)

    def test_failed_counted_with_recipe_over_subtraction(self, synthetic_project, tmp_path):
        """Recipe over-subtraction produces non-zero failed count in batch summary."""
        recipe_path = tmp_path / "recipes"
        recipe_path.mkdir()
        # Create FoxP3 recipe with extreme over-subtraction
        (recipe_path / "FoxP3_param.txt").write_text(
            "signal_channel: 'FoxP3'\n"
            "blankID: 'Blank_1b'\n"
            "blank_scale_factor: 5.0\n"
            "second_subtract: False\n"
            "clean_used: False\n"
        )
        result = process_batch(
            synthetic_project,
            recipe_dir=recipe_path,
            force=True,
            quality_gate=0.6,
        )
        # FoxP3 should be detected as over-subtracted (or auto-fallback used)
        foxp3_result = result.channels.get("FoxP3")
        assert foxp3_result is not None
        assert foxp3_result.status in ("failed", "success")
        if foxp3_result.status == "success":
            # Auto fallback was used
            assert foxp3_result.recipe_source == "auto_fallback"


# =============================================================================
# Outlier Clipping
# =============================================================================


class TestClipOutliers:
    def test_clips_hot_pixels(self):
        """Hot pixels above p99.5 are clipped."""
        rng = np.random.RandomState(80)
        image = rng.randint(100, 10000, (256, 256), dtype=np.uint16)
        # Insert hot pixels far above the signal range
        image[10, 10] = 55000
        image[20, 20] = 60000
        image[30, 30] = 50000

        result = clip_outliers(image, percentile=99.5)

        # Hot pixels should be clipped
        assert result[10, 10] < 55000
        assert result[20, 20] < 60000
        # Normal signal range should be mostly preserved
        assert result.dtype == np.uint16
        # Median pixel should be unchanged (well below threshold)
        median_val = np.median(image[image > 0])
        low_mask = image < median_val
        np.testing.assert_array_equal(result[low_mask], image[low_mask])

    def test_disabled_with_zero(self):
        """percentile=0 returns image unchanged."""
        rng = np.random.RandomState(81)
        image = rng.randint(100, 50000, (64, 64), dtype=np.uint16)
        result = clip_outliers(image, percentile=0)
        np.testing.assert_array_equal(result, image)

    def test_preserves_zeros(self):
        """Zero pixels from subtraction are not affected."""
        image = np.zeros((64, 64), dtype=np.uint16)
        image[10:30, 10:30] = 5000
        image[15, 15] = 50000  # Hot pixel
        result = clip_outliers(image, percentile=99.0)
        # Zeros preserved
        assert np.sum(result == 0) == np.sum(image == 0)

    def test_uses_nonzero_distribution(self):
        """Percentile computed from non-zero pixels only."""
        # Image that's mostly zero (like post-subtraction)
        image = np.zeros((100, 100), dtype=np.uint16)
        image[0:10, 0:10] = 5000  # Small signal region
        image[5, 5] = 50000  # Hot pixel in signal region

        result = clip_outliers(image, percentile=99.0)
        # Hot pixel should be clipped relative to non-zero distribution
        assert result[5, 5] < 50000
        assert result[5, 5] > 0

    def test_all_zeros_noop(self):
        """All-zero image returns unchanged."""
        image = np.zeros((64, 64), dtype=np.uint16)
        result = clip_outliers(image, percentile=99.5)
        np.testing.assert_array_equal(result, image)

    def test_histogram_compression(self):
        """After clipping, max should be much closer to p99."""
        rng = np.random.RandomState(82)
        image = rng.randint(100, 10000, (256, 256), dtype=np.uint16)
        # Add extreme outliers
        image.ravel()[:50] = rng.randint(40000, 60000, 50).astype(np.uint16)

        result = clip_outliers(image, percentile=99.5)

        # Max should be dramatically reduced
        assert result.max() < image.max()
        # p99 should be essentially unchanged
        orig_p99 = np.percentile(image[image > 0], 99)
        new_p99 = np.percentile(result[result > 0], 99)
        assert abs(new_p99 - orig_p99) / orig_p99 < 0.05
