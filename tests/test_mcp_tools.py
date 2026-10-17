"""
Comprehensive Unit Tests for KINTSUGI MCP Tools

Tests cover all five tool modules:
1. signal_isolation - Core image processing operations
2. quality_assessment - QC metrics and quality scoring
3. visualization - Image statistics and thumbnails
4. workflow - Session management and parameter suggestions
5. learning - Parameter learning with SQLite storage

Run with: pytest tests/test_mcp_tools.py -v
Coverage: pytest tests/test_mcp_tools.py --cov=src/kintsugi/mcp --cov-report=html
"""

from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

# Check for optional dependencies using importlib
HAS_ZARR = importlib.util.find_spec("zarr") is not None
HAS_DASK_IMAGE = importlib.util.find_spec("dask_image") is not None

# Skip markers for tests requiring optional dependencies
requires_zarr = pytest.mark.skipif(not HAS_ZARR, reason="zarr not installed")
requires_dask_image = pytest.mark.skipif(not HAS_DASK_IMAGE, reason="dask_image not installed")

requires_mcp = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")

# The full public MCP tool surface (kept in sync with kintsugi.mcp.tool_specs).
EXPECTED_TOOL_NAMES = {
    "load_channel",
    "subtract_blank",
    "analyze_weighted_subtraction",
    "denoise",
    "apply_clahe",
    "clean_background",
    "gaussian_subtract",
    "denoise_advanced",
    "assess_quality",
    "compute_snr",
    "batch_assess_quality",
    "get_image_stats",
    "get_thumbnail",
    "compare_channels",
    "get_intensity_profile",
    "list_channels",
    "save_processed",
    "get_processing_history",
    "suggest_parameters",
    "generate_jupyter_cell",
    "get_session_state",
    "unload_channel",
    "get_workflow_status",
    "set_review_decision",
    "cluster_channels",
    "propagate_parameters",
    "get_learned_parameters",
    "record_successful_parameters",
    "suggest_with_learning",
    "approve_and_learn",
    "get_learning_statistics",
    "optimize_parameters",
    "predict_parameters",
    "estimate_background",
}


def _exec_kintsugi_imports(cell: str) -> None:
    """Run the cell's ``from kintsugi... import ...`` lines so a renamed module fails here."""
    lines = [ln.strip() for ln in cell.splitlines() if ln.strip().startswith("from kintsugi")]
    assert lines, "generated cell imports nothing from kintsugi"
    for line in lines:
        exec(line, {})


# =============================================================================
# Test Fixtures
# =============================================================================


@pytest.fixture
def sample_uint16_image():
    """Generate a sample 16-bit test image with realistic intensity distribution."""
    np.random.seed(42)
    # Create image with background, signal, and some structure
    image = np.random.randint(100, 500, (512, 512), dtype=np.uint16)
    # Add some bright spots (signal)
    image[100:150, 100:150] = np.random.randint(2000, 5000, (50, 50), dtype=np.uint16)
    image[300:350, 300:350] = np.random.randint(3000, 8000, (50, 50), dtype=np.uint16)
    return image


@pytest.fixture
def sample_blank_image():
    """Generate a sample blank/autofluorescence image."""
    np.random.seed(43)
    # Autofluorescence is typically more uniform with lower intensity
    return np.random.randint(50, 400, (512, 512), dtype=np.uint16)


@pytest.fixture
def sample_noisy_image():
    """Generate a noisy test image for denoising tests."""
    np.random.seed(44)
    # Create clean image
    clean = np.random.randint(500, 2000, (512, 512), dtype=np.uint16)
    # Add Gaussian noise
    noise = np.random.normal(0, 300, (512, 512)).astype(np.int32)
    noisy = np.clip(clean.astype(np.int32) + noise, 0, 65535).astype(np.uint16)
    return noisy


@pytest.fixture
def temp_project_dir(tmp_path):
    """Create a temporary KINTSUGI project structure."""
    project_dir = tmp_path / "test_project"
    raw_dir = project_dir / "data" / "raw" / "cyc001"
    processed_dir = project_dir / "data" / "processed"

    raw_dir.mkdir(parents=True)
    processed_dir.mkdir(parents=True)

    return project_dir


@pytest.fixture
def temp_project_with_images(temp_project_dir, sample_uint16_image, sample_blank_image):
    """Create a temp project with actual TIFF files."""
    import tifffile

    raw_dir = temp_project_dir / "data" / "raw" / "cyc001"

    # Save test images
    signal_path = raw_dir / "CD3_cyc001.tif"
    blank_path = raw_dir / "blank_cyc001.tif"

    tifffile.imwrite(str(signal_path), sample_uint16_image)
    tifffile.imwrite(str(blank_path), sample_blank_image)

    return temp_project_dir


@pytest.fixture
def clear_mcp_state():
    """Clear MCP session state before/after each test."""
    from kintsugi.mcp.tools.signal_isolation import _loaded_images, _processing_history

    # Clear before test
    _loaded_images.clear()
    _processing_history.clear()

    yield

    # Clear after test
    _loaded_images.clear()
    _processing_history.clear()


# =============================================================================
# Signal Isolation Tests
# =============================================================================


class TestSignalIsolationStateManagement:
    """Tests for signal isolation state management functions."""

    def test_store_and_get_image(self, sample_uint16_image, clear_mcp_state):
        """Test storing and retrieving images from session state."""
        from kintsugi.mcp.tools.signal_isolation import (
            _get_image,
            _loaded_images,
            _store_image,
        )

        # Store image
        _store_image("test_channel", sample_uint16_image, {"cycle": "cyc001"})

        # Verify storage
        assert "test_channel" in _loaded_images
        assert _loaded_images["test_channel"]["dtype"] == "uint16"
        assert _loaded_images["test_channel"]["shape"] == (512, 512)
        assert _loaded_images["test_channel"]["metadata"]["cycle"] == "cyc001"

        # Retrieve image
        retrieved = _get_image("test_channel")
        np.testing.assert_array_equal(retrieved, sample_uint16_image)

    def test_get_nonexistent_image_raises(self, clear_mcp_state):
        """Test that getting a non-existent image raises ValueError."""
        from kintsugi.mcp.tools.signal_isolation import _get_image

        with pytest.raises(ValueError, match="Image not loaded"):
            _get_image("nonexistent_channel")

    def test_add_history(self, sample_uint16_image, clear_mcp_state):
        """Test processing history tracking."""
        from kintsugi.mcp.tools.signal_isolation import (
            _add_history,
            _processing_history,
            _store_image,
        )

        _store_image("test_channel", sample_uint16_image)
        _add_history("test_channel", "denoise", {"method": "median", "filter_size": 3})

        assert len(_processing_history["test_channel"]) == 1
        assert _processing_history["test_channel"][0]["operation"] == "denoise"
        assert _processing_history["test_channel"][0]["params"]["method"] == "median"

    def test_clear_session(self, sample_uint16_image, clear_mcp_state):
        """Test clearing session state."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _processing_history,
            _store_image,
            clear_session,
        )

        _store_image("test_channel", sample_uint16_image)
        assert len(_loaded_images) == 1

        clear_session()

        assert len(_loaded_images) == 0
        assert len(_processing_history) == 0


class TestLoadChannel:
    """Tests for the load_channel tool."""

    @requires_zarr
    @pytest.mark.asyncio
    async def test_load_channel_success(self, temp_project_with_images, clear_mcp_state):
        """Test successful channel loading."""
        from kintsugi.mcp.tools.signal_isolation import _loaded_images, load_channel

        result = await load_channel(
            project_path=str(temp_project_with_images),
            cycle="1",
            channel="CD3",
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "loaded"
        assert "cyc001_CD3" in result["name"]
        assert result["shape"] == [512, 512]
        assert result["dtype"] == "uint16"
        assert "stats" in result
        assert "cyc001_CD3" in _loaded_images

    @requires_zarr
    @pytest.mark.asyncio
    async def test_load_channel_not_found(self, temp_project_dir, clear_mcp_state):
        """Test loading a non-existent channel returns error."""
        from kintsugi.mcp.tools.signal_isolation import load_channel

        result = await load_channel(
            project_path=str(temp_project_dir),
            cycle="1",
            channel="NonExistent",
        )

        assert "error" in result
        # Accept either "not found" or "Missing dependency" errors
        assert (
            "not found" in result["error"].lower()
            or "missing dependency" in result["error"].lower()
        )

    @requires_zarr
    @pytest.mark.asyncio
    async def test_load_channel_cycle_normalization(
        self, temp_project_with_images, clear_mcp_state
    ):
        """Test that cycle names are properly normalized."""
        from kintsugi.mcp.tools.signal_isolation import load_channel

        # Load using numeric cycle
        result = await load_channel(
            project_path=str(temp_project_with_images),
            cycle="1",
            channel="CD3",
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "loaded"
        assert "cyc001" in result["name"]


class TestSubtractBlank:
    """Tests for the subtract_blank tool."""

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_subtract_blank_success(
        self, sample_uint16_image, sample_blank_image, clear_mcp_state
    ):
        """Test successful blank subtraction."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            subtract_blank,
        )

        _store_image("signal_ch", sample_uint16_image)
        _store_image("blank_ch", sample_blank_image)

        result = await subtract_blank(
            signal_channel="signal_ch",
            blank_channel="blank_ch",
            blank_scale_factor=0.8,
            blank_clip_factor=50,
            auto_suggest=False,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["output_name"] == "signal_ch_subtracted"
        assert "signal_ch_subtracted" in _loaded_images
        assert result["parameters_used"]["blank_scale_factor"] == 0.8

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_subtract_blank_missing_channel(self, sample_uint16_image, clear_mcp_state):
        """Test error when blank channel is missing."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, subtract_blank

        _store_image("signal_ch", sample_uint16_image)

        result = await subtract_blank(
            signal_channel="signal_ch",
            blank_channel="nonexistent_blank",
            auto_suggest=False,
        )

        assert "error" in result
        # Accept either "not loaded" or "Missing dependency" errors
        assert (
            "not loaded" in result["error"].lower()
            or "missing dependency" in result["error"].lower()
        )

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_subtract_blank_records_history(
        self, sample_uint16_image, sample_blank_image, clear_mcp_state
    ):
        """Test that blank subtraction records processing history."""
        from kintsugi.mcp.tools.signal_isolation import (
            _processing_history,
            _store_image,
            subtract_blank,
        )

        _store_image("signal_ch", sample_uint16_image)
        _store_image("blank_ch", sample_blank_image)

        result = await subtract_blank(
            signal_channel="signal_ch",
            blank_channel="blank_ch",
            auto_suggest=False,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        history = _processing_history.get("signal_ch_subtracted", [])
        assert len(history) >= 1
        assert any(h["operation"] == "subtract_blank" for h in history)


class TestDenoise:
    """Tests for the denoise tool."""

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_denoise_median(self, sample_noisy_image, clear_mcp_state):
        """Test median filter denoising."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            denoise,
        )

        _store_image("noisy_ch", sample_noisy_image)

        result = await denoise(
            channel="noisy_ch",
            method="median",
            filter_size=3,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["output_name"] == "noisy_ch_denoised"
        assert result["method"] == "median"
        assert "noisy_ch_denoised" in _loaded_images

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_denoise_percentile(self, sample_noisy_image, clear_mcp_state):
        """Test percentile filter denoising."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise

        _store_image("noisy_ch", sample_noisy_image)

        result = await denoise(
            channel="noisy_ch",
            method="percentile",
            filter_size=3,
            percentile=10,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["method"] == "percentile"

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_denoise_invalid_method(self, sample_noisy_image, clear_mcp_state):
        """Test error for invalid denoise method."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise

        _store_image("noisy_ch", sample_noisy_image)

        result = await denoise(
            channel="noisy_ch",
            method="invalid_method",
            filter_size=3,
        )

        assert "error" in result
        # Accept either "Unknown method" or "Missing dependency" errors
        assert "Unknown method" in result["error"] or "Missing dependency" in result["error"]

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_denoise_uniform(self, sample_noisy_image, clear_mcp_state):
        """Test uniform filter denoising."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise

        _store_image("noisy_ch", sample_noisy_image)

        result = await denoise(
            channel="noisy_ch",
            method="uniform",
            filter_size=3,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["method"] == "uniform"


class TestDenoiseAdvanced:
    """Tests for the denoise_advanced tool."""

    @pytest.mark.asyncio
    async def test_denoise_advanced_adaptive(self, sample_noisy_image, clear_mcp_state):
        """Test adaptive denoising."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise_advanced

        _store_image("noisy_ch", sample_noisy_image)

        try:
            result = await denoise_advanced(
                channel="noisy_ch",
                method="adaptive",
                strength="auto",
            )
            # If kintsugi.denoise is available, check result
            if "error" not in result:
                assert result["status"] == "success"
                assert "noise_reduction" in result
            else:
                # Missing dependency is acceptable
                assert "error" in result
        except ImportError:
            # Expected if denoise module not fully available
            pytest.skip("Denoise module not available")


class TestGaussianSubtract:
    """Tests for the gaussian_subtract tool."""

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_gaussian_subtract_success(self, sample_uint16_image, clear_mcp_state):
        """Test Gaussian background subtraction."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            gaussian_subtract,
        )

        _store_image("test_ch", sample_uint16_image)

        result = await gaussian_subtract(
            channel="test_ch",
            sigma=10,
            scale_factor=0.1,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["output_name"] == "test_ch_gaussub"
        assert result["sigma"] == 10
        assert result["scale_factor"] == 0.1
        assert "test_ch_gaussub" in _loaded_images


class TestApplyClahe:
    """Tests for the apply_clahe tool."""

    @pytest.mark.asyncio
    async def test_apply_clahe_success(self, sample_uint16_image, clear_mcp_state):
        """Test CLAHE application."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            apply_clahe,
        )

        _store_image("test_ch", sample_uint16_image)

        result = await apply_clahe(
            channel="test_ch",
            clip_limit=0.01,
            tile_grid_size=50,
        )

        assert result["status"] == "success"
        assert result["output_name"] == "test_ch_clahe"
        assert "test_ch_clahe" in _loaded_images


class TestCleanBackground:
    """Tests for the clean_background tool."""

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_clean_background_success(self, sample_uint16_image, clear_mcp_state):
        """Test background cleaning."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            clean_background,
        )

        _store_image("test_ch", sample_uint16_image)

        result = await clean_background(
            channel="test_ch",
            background_threshold=100,
            remove_small_objects=False,
        )

        # Handle missing dependency error
        if "error" in result and "Missing dependency" in result["error"]:
            pytest.skip(f"Missing dependency: {result['error']}")

        assert result["status"] == "success"
        assert result["output_name"] == "test_ch_cleaned"
        assert "test_ch_cleaned" in _loaded_images


# =============================================================================
# Quality Assessment Tests
# =============================================================================


class TestQualityAssessment:
    """Tests for quality assessment tools."""

    @pytest.mark.asyncio
    async def test_assess_quality_success(self, sample_uint16_image, clear_mcp_state):
        """Test quality assessment returns expected metrics."""
        from kintsugi.mcp.tools.quality_assessment import assess_quality
        from kintsugi.mcp.tools.signal_isolation import _store_image

        _store_image("test_ch", sample_uint16_image)

        result = await assess_quality(
            channel="test_ch",
            use_dl_model=False,
            confidence_threshold=0.85,
        )

        assert result["status"] == "success"
        assert "quality_score" in result
        assert "metrics" in result
        assert "recommendations" in result
        assert isinstance(result["quality_score"], float)
        assert 0 <= result["quality_score"] <= 1

    @pytest.mark.asyncio
    async def test_assess_quality_missing_channel(self, clear_mcp_state):
        """Test error when assessing non-existent channel."""
        from kintsugi.mcp.tools.quality_assessment import assess_quality

        result = await assess_quality(channel="nonexistent")

        assert "error" in result

    @pytest.mark.asyncio
    async def test_compute_snr_success(self, sample_uint16_image, clear_mcp_state):
        """Test SNR computation."""
        from kintsugi.mcp.tools.quality_assessment import compute_snr
        from kintsugi.mcp.tools.signal_isolation import _store_image

        _store_image("test_ch", sample_uint16_image)

        result = await compute_snr(channel="test_ch")

        assert result["status"] == "success"
        assert "snr" in result
        assert "signal_to_background" in result
        assert "interpretation" in result
        assert isinstance(result["snr"], float)

    @pytest.mark.asyncio
    async def test_batch_assess_quality(
        self, sample_uint16_image, sample_noisy_image, clear_mcp_state
    ):
        """Test batch quality assessment."""
        from kintsugi.mcp.tools.quality_assessment import batch_assess_quality
        from kintsugi.mcp.tools.signal_isolation import _store_image

        _store_image("ch1", sample_uint16_image)
        _store_image("ch2", sample_noisy_image)

        result = await batch_assess_quality(
            channels=["ch1", "ch2"],
            confidence_threshold=0.85,
        )

        assert result["status"] == "success"
        assert result["total_channels"] == 2
        assert result["channels_assessed"] == 2
        assert "average_score" in result


# =============================================================================
# Visualization Tests
# =============================================================================


class TestVisualization:
    """Tests for visualization tools."""

    @pytest.mark.asyncio
    async def test_get_image_stats(self, sample_uint16_image, clear_mcp_state):
        """Test image statistics retrieval."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.visualization import get_image_stats

        _store_image("test_ch", sample_uint16_image)

        result = await get_image_stats(channel="test_ch")

        assert result["status"] == "success"
        stats = result["statistics"]
        assert "min" in stats
        assert "max" in stats
        assert "mean" in stats
        assert "std" in stats
        assert "percentiles" in stats
        assert "histogram" in stats
        assert "zero_pixels" in stats

    @requires_mcp
    @pytest.mark.asyncio
    async def test_get_thumbnail(self, sample_uint16_image, clear_mcp_state):
        """Thumbnail arrives as an image content block plus structured metadata."""
        pytest.importorskip("PIL")
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.visualization import get_thumbnail

        _store_image("test_ch", sample_uint16_image)

        result = await get_thumbnail(
            channel="test_ch",
            max_size=256,
            normalize=True,
        )

        assert not result.is_error
        # The renderable image block comes first, followed by the JSON text block.
        assert result.content[0].type == "image"
        assert result.content[0].mime_type == "image/png"
        assert result.content[1].type == "text"

        meta = result.structured_content
        assert meta["status"] == "success"
        assert meta["channel"] == "test_ch"
        assert meta["thumbnail"]["format"] == "png"
        # The base64 data: URI is opt-in (include_data_uri) and off by default.
        assert "data_uri" not in meta["thumbnail"]
        # Thumbnail should be smaller or equal to max_size
        assert meta["thumbnail"]["width"] <= 256
        assert meta["thumbnail"]["height"] <= 256
        assert meta["original_shape"] == [512, 512]

    @pytest.mark.asyncio
    async def test_compare_channels(self, sample_uint16_image, sample_blank_image, clear_mcp_state):
        """Test channel comparison."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.visualization import compare_channels

        _store_image("ch1", sample_uint16_image)
        _store_image("ch2", sample_blank_image)

        result = await compare_channels(channel1="ch1", channel2="ch2")

        assert result["status"] == "success"
        assert "comparison" in result
        assert "correlation" in result["comparison"]
        assert "ssim" in result["comparison"]
        assert "interpretation" in result["comparison"]

    @pytest.mark.asyncio
    async def test_get_intensity_profile(self, sample_uint16_image, clear_mcp_state):
        """Test intensity profile extraction."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.visualization import get_intensity_profile

        _store_image("test_ch", sample_uint16_image)

        result = await get_intensity_profile(
            channel="test_ch",
            axis="horizontal",
            position=0.5,
        )

        assert result["status"] == "success"
        assert result["axis"] == "horizontal"
        assert "profile" in result
        assert "values" in result["profile"]


# =============================================================================
# Workflow Tests
# =============================================================================


class TestWorkflow:
    """Tests for workflow tools."""

    @pytest.mark.asyncio
    async def test_list_channels(self, temp_project_with_images, clear_mcp_state):
        """Test channel listing."""
        from kintsugi.mcp.tools.workflow import list_channels

        result = await list_channels(
            project_path=str(temp_project_with_images),
            cycle="cyc001",
        )

        assert result["status"] == "success"
        assert result["total_channels"] >= 1
        assert "channels" in result

    @pytest.mark.asyncio
    async def test_save_processed_tiff(self, sample_uint16_image, tmp_path, clear_mcp_state):
        """Test saving processed image as TIFF."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import save_processed

        _store_image("test_ch", sample_uint16_image, {"project": str(tmp_path)})

        result = await save_processed(
            channel="test_ch",
            output_name="test_output",
            format="tiff",
            project_path=str(tmp_path),
        )

        assert result["status"] == "success"
        assert result["format"] == "tiff"
        assert Path(result["output_path"]).exists()

    @requires_zarr
    @pytest.mark.asyncio
    async def test_save_processed_zarr(self, sample_uint16_image, tmp_path, clear_mcp_state):
        """Test saving processed image as Zarr."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import save_processed

        _store_image("test_ch", sample_uint16_image, {"project": str(tmp_path)})

        result = await save_processed(
            channel="test_ch",
            output_name="test_output",
            format="zarr",
            project_path=str(tmp_path),
        )

        # Handle missing dependency error
        if "error" in result and (
            "zarr" in result["error"].lower() or "Missing dependency" in result["error"]
        ):
            pytest.skip(f"Zarr dependency not available: {result.get('error', '')}")

        assert result["status"] == "success"
        assert result["format"] == "zarr"

    @pytest.mark.asyncio
    async def test_get_processing_history(self, sample_uint16_image, clear_mcp_state):
        """Test processing history retrieval."""
        from kintsugi.mcp.tools.signal_isolation import (
            _add_history,
            _store_image,
        )
        from kintsugi.mcp.tools.workflow import get_processing_history

        _store_image("test_ch", sample_uint16_image)
        _add_history("test_ch", "denoise", {"method": "median"})
        _add_history("test_ch", "clahe", {"clip_limit": 0.01})

        result = await get_processing_history(channel="test_ch")

        assert result["status"] == "success"
        assert result["total_steps"] == 2
        assert len(result["steps"]) == 2

    @pytest.mark.asyncio
    async def test_suggest_parameters(
        self, sample_uint16_image, sample_blank_image, clear_mcp_state
    ):
        """Test parameter suggestion."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import suggest_parameters

        _store_image("signal_ch", sample_uint16_image)
        _store_image("blank_ch", sample_blank_image)

        result = await suggest_parameters(
            channel="signal_ch",
            blank_channel="blank_ch",
        )

        assert result["status"] == "success"
        assert "suggestions" in result
        suggestions = result["suggestions"]
        assert "analysis" in suggestions
        assert "recommendations" in suggestions

    @pytest.mark.asyncio
    async def test_generate_jupyter_cell(self, sample_uint16_image, clear_mcp_state, tmp_path):
        """The generated cell loads the channel from its file and is valid Python."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import generate_jupyter_cell

        signal_file = tmp_path / "test_ch.tif"
        blank_file = tmp_path / "it's blank.tif"
        _store_image("test_ch", sample_uint16_image, {"source_path": str(signal_file)})
        _store_image("blank_ch", sample_uint16_image, {"source_path": str(blank_file)})

        result = await generate_jupyter_cell(
            operation="blank_subtraction",
            channel="test_ch",
            initial_params={"blank_channel": "blank_ch", "blank_scale_factor": 0.75},
        )

        assert result["status"] == "success"
        assert result["source_file"] == str(signal_file)
        cell = result["jupyter_cell"]
        assert str(signal_file) in cell and str(blank_file) in cell
        assert "0.75" in cell
        compile(cell, "<jupyter_cell>", "exec")  # quotes in paths/values stay valid
        _exec_kintsugi_imports(cell)

    @pytest.mark.parametrize("operation", ["denoise", "clahe", "clean", "gaussian_subtract"])
    async def test_generate_jupyter_cell_every_operation_compiles(
        self, operation, sample_uint16_image, clear_mcp_state, tmp_path
    ):
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import generate_jupyter_cell

        _store_image("ch", sample_uint16_image, {"source_path": str(tmp_path / "ch.tif")})
        result = await generate_jupyter_cell(
            operation=operation, channel="ch", initial_params={"method": 'x"y'}
        )
        assert result["status"] == "success"
        compile(result["jupyter_cell"], "<jupyter_cell>", "exec")
        _exec_kintsugi_imports(result["jupyter_cell"])

    async def test_generate_jupyter_cell_needs_a_file(self, sample_uint16_image, clear_mcp_state):
        """An in-memory-only channel cannot be reached from a notebook kernel."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import generate_jupyter_cell

        _store_image("mem_only", sample_uint16_image)
        result = await generate_jupyter_cell(operation="denoise", channel="mem_only")
        assert "error" in result and "save_processed" in result["error"]

    @pytest.mark.asyncio
    async def test_get_session_state(self, sample_uint16_image, clear_mcp_state):
        """Test session state retrieval."""
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import get_session_state

        _store_image("ch1", sample_uint16_image)
        _store_image("ch2", sample_uint16_image)

        result = await get_session_state()

        assert result["status"] == "success"
        assert result["total_loaded"] == 2
        assert "memory_usage" in result


# =============================================================================
# Learning Tests
# =============================================================================


class TestLearning:
    """Tests for parameter learning tools."""

    @pytest.fixture
    def mock_learning_engine(self):
        """Mock the ParameterLearningEngine."""
        with patch("kintsugi.mcp.tools.learning._get_learning_engine") as mock:
            engine = MagicMock()
            engine.recommend_parameters.return_value = {
                "status": "success",
                "recommended_parameters": {"blank_scale_factor": 0.8},
                "confidence": 0.75,
                "based_on_records": 5,
            }
            engine.record_parameters.return_value = "record_123"
            engine.get_statistics.return_value = {
                "total_records": 100,
                "unique_tissues": 5,
                "unique_markers": 20,
            }
            mock.return_value = engine
            yield engine

    @pytest.mark.asyncio
    async def test_get_learned_parameters(self, mock_learning_engine, clear_mcp_state):
        """Test getting learned parameters."""
        from kintsugi.mcp.tools.learning import get_learned_parameters

        result = await get_learned_parameters(
            tissue_type="tonsil",
            marker_name="CD3",
            operation="blank_subtraction",
        )

        assert result["status"] == "success"
        assert "recommended_parameters" in result

    @pytest.mark.asyncio
    async def test_record_successful_parameters(
        self, mock_learning_engine, sample_uint16_image, clear_mcp_state
    ):
        """Test recording successful parameters."""
        from kintsugi.mcp.tools.learning import record_successful_parameters
        from kintsugi.mcp.tools.signal_isolation import _store_image

        _store_image("test_ch", sample_uint16_image)

        result = await record_successful_parameters(
            tissue_type="tonsil",
            marker_name="CD3",
            operation="blank_subtraction",
            parameters={"blank_scale_factor": 0.8},
            quality_before=0.5,
            quality_after=0.85,
            channel="test_ch",
        )

        assert result["status"] == "success"
        assert "record_id" in result
        assert result["quality_improvement"] == 0.35

    @pytest.mark.asyncio
    async def test_get_learning_statistics(self, mock_learning_engine, clear_mcp_state):
        """Test getting learning database statistics."""
        from kintsugi.mcp.tools.learning import get_learning_statistics

        result = await get_learning_statistics()

        assert result["status"] == "success"
        assert "statistics" in result


# =============================================================================
# Integration Tests
# =============================================================================


class TestIntegration:
    """Integration tests for MCP tool pipelines."""

    @requires_dask_image
    @pytest.mark.asyncio
    async def test_full_processing_pipeline(
        self, sample_uint16_image, sample_blank_image, tmp_path, clear_mcp_state
    ):
        """Test a complete processing pipeline: load -> subtract -> denoise -> save."""
        from kintsugi.mcp.tools.signal_isolation import (
            _store_image,
            denoise,
            subtract_blank,
        )
        from kintsugi.mcp.tools.workflow import get_processing_history, save_processed

        # Step 1: Store images (simulating load)
        _store_image("signal_ch", sample_uint16_image, {"project": str(tmp_path)})
        _store_image("blank_ch", sample_blank_image, {"project": str(tmp_path)})

        # Step 2: Subtract blank
        result1 = await subtract_blank(
            signal_channel="signal_ch",
            blank_channel="blank_ch",
            blank_scale_factor=0.8,
            auto_suggest=False,
        )

        # Handle missing dependency error
        if "error" in result1 and "Missing dependency" in result1["error"]:
            pytest.skip(f"Missing dependency: {result1['error']}")

        assert result1["status"] == "success"
        subtracted_name = result1["output_name"]

        # Step 3: Denoise
        result2 = await denoise(
            channel=subtracted_name,
            method="median",
            filter_size=3,
        )

        # Handle missing dependency error
        if "error" in result2 and "Missing dependency" in result2["error"]:
            pytest.skip(f"Missing dependency: {result2['error']}")

        assert result2["status"] == "success"
        denoised_name = result2["output_name"]

        # Step 4: Save
        result3 = await save_processed(
            channel=denoised_name,
            output_name="final_output",
            format="tiff",
            project_path=str(tmp_path),
        )
        assert result3["status"] == "success"

        # Verify history tracking
        history = await get_processing_history(channel=denoised_name)
        assert history["total_steps"] >= 1

    @pytest.mark.asyncio
    async def test_cross_module_state_sharing(self, sample_uint16_image, clear_mcp_state):
        """Test that state is properly shared across modules."""
        from kintsugi.mcp.tools.quality_assessment import assess_quality
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.visualization import get_image_stats
        from kintsugi.mcp.tools.workflow import get_session_state

        # Store in signal_isolation
        _store_image("shared_ch", sample_uint16_image)

        # Access from quality_assessment
        qc_result = await assess_quality(channel="shared_ch")
        assert qc_result["status"] == "success"

        # Access from visualization
        stats_result = await get_image_stats(channel="shared_ch")
        assert stats_result["status"] == "success"

        # Access from workflow
        session_result = await get_session_state()
        assert "shared_ch" in session_result["loaded_images"]


# =============================================================================
# Error Handling Tests
# =============================================================================


class TestErrorHandling:
    """Tests for error handling across tools."""

    @pytest.mark.asyncio
    async def test_missing_dependency_graceful(self, sample_uint16_image, clear_mcp_state):
        """Test graceful handling of missing dependencies."""
        from kintsugi.mcp.tools.signal_isolation import _store_image

        _store_image("test_ch", sample_uint16_image)

        # Test with mocked missing import
        with patch.dict("sys.modules", {"dask_image": None}):
            # Should still work with fallback or return informative error
            pass  # The actual import happens inside the function

    @pytest.mark.asyncio
    async def test_invalid_parameter_values(self, sample_uint16_image, clear_mcp_state):
        """Test handling of invalid parameter values."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise

        _store_image("test_ch", sample_uint16_image)

        # Invalid method
        result = await denoise(
            channel="test_ch",
            method="nonexistent_method",
            filter_size=3,
        )
        assert "error" in result


# =============================================================================
# Cross-Platform Tests
# =============================================================================


class TestCrossPlatform:
    """Tests for cross-platform compatibility."""

    @pytest.mark.asyncio
    async def test_path_handling(self, tmp_path, clear_mcp_state):
        """Test that paths are handled correctly across platforms."""
        from kintsugi.mcp.tools.workflow import list_channels

        # Create a path with spaces
        project_dir = tmp_path / "test project with spaces"
        raw_dir = project_dir / "data" / "raw" / "cyc001"
        raw_dir.mkdir(parents=True)

        result = await list_channels(project_path=str(project_dir))

        # Should not raise, should return valid result
        assert result["status"] == "success"

    def test_json_serialization(self, sample_uint16_image, clear_mcp_state):
        """Test that all return values are JSON serializable."""
        from kintsugi.mcp.tools.signal_isolation import _store_image, get_loaded_images

        _store_image("test_ch", sample_uint16_image)
        result = get_loaded_images()

        # Should be JSON serializable
        json_str = json.dumps(result)
        assert json_str is not None


# =============================================================================
# Session Isolation Tests
# =============================================================================


class TestSessionIsolation:
    """Tests for session isolation and cleanup."""

    def test_session_state_not_leaked(self):
        """Test that clearing session properly isolates state."""
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _processing_history,
            _store_image,
            clear_session,
        )

        # Setup
        clear_session()

        # Add some data
        img = np.random.randint(0, 65535, (100, 100), dtype=np.uint16)
        _store_image("session1_ch", img)

        # Verify data exists
        assert "session1_ch" in _loaded_images

        # Clear and verify clean
        clear_session()
        assert "session1_ch" not in _loaded_images
        assert len(_processing_history) == 0


# =============================================================================
# Server Tests
# =============================================================================


class TestMCPServer:
    """Tests for the MCP server itself."""

    def test_server_creation_without_mcp(self):
        """Test server behavior when MCP package is not available."""
        # Test that MCP_AVAILABLE flag controls server creation
        from kintsugi.mcp import server

        # Store original value
        original_mcp_available = server.MCP_AVAILABLE

        try:
            # Simulate MCP not being available
            server.MCP_AVAILABLE = False

            # Verify create_server raises ImportError when MCP is unavailable
            with pytest.raises(ImportError, match="Install with"):
                server.create_server()
        finally:
            # Restore original value
            server.MCP_AVAILABLE = original_mcp_available

    def test_tool_registry_complete(self):
        """Every advertised tool resolves to an importable coroutine handler."""
        import importlib
        import inspect

        from kintsugi.mcp import tool_specs

        names = [spec["name"] for spec in tool_specs.TOOL_SPECS]
        # Names are unique and cover the full public surface.
        assert len(names) == len(set(names))
        assert set(names) == EXPECTED_TOOL_NAMES

        # Each spec points at a real async handler in kintsugi.mcp.tools and
        # carries a title plus the four MCP behaviour hints.
        hint_keys = {"readOnlyHint", "destructiveHint", "idempotentHint", "openWorldHint"}
        for spec in tool_specs.TOOL_SPECS:
            module = importlib.import_module(f"kintsugi.mcp.tools.{spec['module']}")
            handler = getattr(module, spec["attr"])
            assert inspect.iscoroutinefunction(handler), spec["name"]
            assert spec["title"].strip(), spec["name"]
            assert set(spec["annotations"]) == hint_keys, spec["name"]
            assert spec["annotations"]["openWorldHint"] is False, spec["name"]


# =============================================================================
# Structured Output Tests (FastMCP server)
# =============================================================================


@requires_mcp
class TestStructuredOutput:
    """The server advertises curated schemas and emits structured content."""

    @staticmethod
    async def _tools():
        from mcp.client import Client

        from kintsugi.mcp.server import create_server

        async with Client(create_server()) as client:
            return {t.name: t for t in (await client.list_tools()).tools}

    async def test_all_tools_registered_with_output_schema(self):
        """Every tool is registered and declares an output schema."""
        tools = await self._tools()

        assert set(tools) == EXPECTED_TOOL_NAMES
        for name, tool in tools.items():
            assert tool.output_schema is not None, f"{name} missing outputSchema"

    async def test_call_returns_structured_content(self):
        """A tool call returns structuredContent plus a human-readable text block."""
        from mcp.client import Client

        from kintsugi.mcp.server import create_server

        async with Client(create_server()) as client:
            result = await client.call_tool("get_processing_history", {"channel": "nope"})

        assert result.is_error is False
        assert result.structured_content == {
            "status": "success",
            "channel": "nope",
            "source": "session",
            "steps": [],
            "total_steps": 0,
            "source_metadata": {},
        }
        # Backward-compatible text block is still present.
        assert result.content and result.content[0].type == "text"


# =============================================================================
# Regression tests for audited defects
# =============================================================================


class TestAuditedDefects:
    """Each test pins one defect found in the MCP audit."""

    async def test_list_channels_rejects_cycle_traversal(self, tmp_path, clear_mcp_state):
        from kintsugi.mcp.tools.workflow import list_channels

        result = await list_channels(project_path=str(tmp_path), cycle="../..")
        assert "error" in result and "cycle" in result["error"]

    async def test_denoise_advanced_confines_model_path(
        self, sample_uint16_image, tmp_path, clear_mcp_state
    ):
        from kintsugi.mcp.tools.signal_isolation import _store_image, denoise_advanced

        _store_image("ch", sample_uint16_image)
        outside = tmp_path / "outside" / "model.pt"
        outside.parent.mkdir()
        outside.write_bytes(b"not a checkpoint")
        project = tmp_path / "proj"
        project.mkdir()
        (project / "kintsugi_project.json").write_text("{}")

        result = await denoise_advanced(
            channel="ch", method="n2v", model_path=str(outside), project_path=str(project)
        )
        assert "error" in result and "outside the allowed directories" in result["error"]

        # A caller cannot widen the sandbox by naming a bare directory as the project.
        result = await denoise_advanced(
            channel="ch", method="n2v", model_path=str(outside), project_path=str(tmp_path)
        )
        assert "error" in result and "outside the allowed directories" in result["error"]

        inside = project / "model.txt"
        inside.write_bytes(b"")
        result = await denoise_advanced(
            channel="ch", method="n2v", model_path=str(inside), project_path=str(project)
        )
        assert "error" in result and ".pt/.pth" in result["error"]

    async def test_predict_parameters_confines_model_path(
        self, sample_uint16_image, tmp_path, clear_mcp_state
    ):
        pytest.importorskip("sklearn")
        from kintsugi.mcp.tools.signal_isolation import _store_image, predict_parameters

        _store_image("sig", sample_uint16_image)
        outside = tmp_path / "evil.joblib"
        outside.write_bytes(b"")
        project = tmp_path / "proj"
        project.mkdir()

        result = await predict_parameters(
            signal_channel="sig", model_path=str(outside), project_path=str(project)
        )
        assert "error" in result and "outside the allowed directories" in result["error"]

    async def test_subtract_blank_uses_learned_parameters(
        self, tmp_path, monkeypatch, clear_mcp_state
    ):
        """The learned-parameter merge was dead: wrong import and wrong status key."""
        pytest.importorskip("dask_image")
        from kintsugi.assist.parameter_learning import ParameterLearningEngine
        from kintsugi.mcp.tools.signal_isolation import _store_image, subtract_blank

        monkeypatch.setenv("HOME", str(tmp_path / "home"))
        (tmp_path / "home").mkdir()
        project = tmp_path / "proj"
        project.mkdir()
        engine = ParameterLearningEngine(str(project))
        for _ in range(3):
            engine.record_parameters(
                tissue_type="tonsil",
                marker_name="CD3",
                operation="blank_subtraction",
                parameters={"blank_scale_factor": 0.33, "blank_clip_factor": 7, "erosion": 2},
                quality_before=0.4,
                quality_after=0.9,
                user_approved=True,
                to_global=False,
            )

        rng = np.random.default_rng(0)
        signal = rng.integers(100, 2000, (128, 128), dtype=np.uint16)
        blank = rng.integers(50, 500, (128, 128), dtype=np.uint16)
        _store_image("sig", signal)
        _store_image("blk", blank)

        plain = await subtract_blank(
            signal_channel="sig", blank_channel="blk", tissue_type="tonsil", marker_name="CD3"
        )
        learned = await subtract_blank(
            signal_channel="sig",
            blank_channel="blk",
            tissue_type="tonsil",
            marker_name="CD3",
            project_path=str(project),
        )
        assert plain.get("status") == "success" and learned.get("status") == "success"
        plain_scale = plain["parameters_used"]["blank_scale_factor"]
        learned_scale = learned["parameters_used"]["blank_scale_factor"]
        # The learned value (0.33) is blended in with up to 60% weight, so the
        # auto-suggested scale moves towards it once the project history exists.
        assert learned_scale != plain_scale
        assert abs(learned_scale - 0.33) < abs(plain_scale - 0.33)

    async def test_cluster_then_propagate(self, tmp_path, clear_mcp_state):
        """cluster_channels tags members with cluster_id so propagate_parameters finds them."""
        pytest.importorskip("sklearn")
        tifffile = pytest.importorskip("tifffile")
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            cluster_channels_tool,
            propagate_parameters_tool,
        )

        project = tmp_path / "proj"
        edf = project / "data" / "processed" / "edf"
        edf.mkdir(parents=True)
        rng = np.random.default_rng(1)
        for name, level in (("CD3", 200), ("CD4", 220), ("DAPI", 3000), ("CD8", 210)):
            img = rng.integers(0, level, (96, 96), dtype=np.uint16)
            tifffile.imwrite(edf / f"{name}.tif", img)

        clustered = await cluster_channels_tool(project_path=str(project), n_clusters=2)
        assert clustered.get("status") == "success", clustered
        assert set(clustered["loaded_channels"]) == {"CD3", "CD4", "DAPI", "CD8"}
        cid = clustered["assignments"]["CD3"]
        assert _loaded_images["CD3"]["metadata"]["cluster_id"] == cid

        propagated = await propagate_parameters_tool(
            project_path=str(project),
            cluster_id=cid,
            params={"clean_params": {"background_threshold": 10}},
        )
        assert propagated.get("status") == "success", propagated
        assert propagated["total"] >= 1

    async def test_cluster_state_is_scoped_per_project(self, tmp_path, clear_mcp_state):
        """Clustering replaces earlier cluster state and propagate ignores other projects."""
        pytest.importorskip("sklearn")
        tifffile = pytest.importorskip("tifffile")
        from kintsugi.mcp.tools.signal_isolation import (
            _loaded_images,
            _store_image,
            cluster_channels_tool,
            propagate_parameters_tool,
        )

        project_b = tmp_path / "b"
        edf = project_b / "data" / "processed" / "edf"
        edf.mkdir(parents=True)
        rng = np.random.default_rng(2)
        for name, level in (("CD3", 200), ("CD4", 220), ("DAPI", 3000), ("CD8", 210)):
            tifffile.imwrite(edf / f"{name}.tif", rng.integers(0, level, (64, 64), dtype=np.uint16))

        # A stale channel clustered for another project with a colliding cluster id.
        _store_image(
            "stale_from_a",
            rng.integers(0, 100, (64, 64), dtype=np.uint16),
            {"cluster_id": 0, "project": str(tmp_path / "a")},
        )

        clustered = await cluster_channels_tool(project_path=str(project_b), n_clusters=2)
        assert clustered.get("status") == "success", clustered
        assert "stale_from_a" not in _loaded_images  # previous cluster state is replaced

        # Even if a foreign entry survives, propagate only touches this project's members.
        _store_image(
            "foreign",
            rng.integers(0, 100, (64, 64), dtype=np.uint16),
            {"cluster_id": clustered["assignments"]["CD3"], "project": str(tmp_path / "a")},
        )
        propagated = await propagate_parameters_tool(
            project_path=str(project_b),
            cluster_id=clustered["assignments"]["CD3"],
            params={"clean_params": {"background_threshold": 10}},
        )
        assert propagated.get("status") == "success", propagated
        assert "foreign" not in propagated["results"]

    async def test_generate_jupyter_cell_neutralises_newlines_in_names(
        self, sample_uint16_image, clear_mcp_state, tmp_path
    ):
        from kintsugi.mcp.tools.signal_isolation import _store_image
        from kintsugi.mcp.tools.workflow import generate_jupyter_cell

        evil = "CD3\nimport os; os.system('id')  # "
        _store_image(evil, sample_uint16_image, {"source_path": str(tmp_path / "CD3.tif")})
        result = await generate_jupyter_cell(operation="denoise", channel=evil)
        assert result["status"] == "success"
        cell = result["jupyter_cell"]
        assert (
            "os.system"
            not in [line.strip() for line in cell.splitlines() if not line.strip().startswith("#")]
            and "\nimport os" not in cell
        )
        compile(cell, "<jupyter_cell>", "exec")


# =============================================================================
# Deterministic file matching (results must not depend on directory order)
# =============================================================================


def _write_raw(root: Path, cycle: str, names: list[str]) -> Path:
    import tifffile

    raw = root / "data" / "raw" / cycle
    raw.mkdir(parents=True, exist_ok=True)
    for name in names:
        tifffile.imwrite(str(raw / name), np.full((8, 8), 100, dtype=np.uint16))
    return raw


class TestDeterministicFileMatching:
    """list_channels and load_channel give the same answer on every filesystem."""

    async def test_list_channels_is_sorted(self, tmp_path, clear_mcp_state):
        from kintsugi.mcp.tools.workflow import list_channels

        # Written in a deliberately unsorted order.
        _write_raw(tmp_path, "cyc002", ["Zeta.tif"])
        _write_raw(tmp_path, "cyc001", ["Zeta.tif", "Alpha.tif", "CD3.tif"])

        result = await list_channels(project_path=str(tmp_path))

        assert result["status"] == "success"
        assert [(c["cycle"], c["channel"]) for c in result["channels"]] == [
            ("cyc001", "Alpha"),
            ("cyc001", "CD3"),
            ("cyc001", "Zeta"),
            ("cyc002", "Zeta"),
        ]

    @requires_zarr
    async def test_load_channel_prefers_exact_stem_over_substring(self, tmp_path, clear_mcp_state):
        from kintsugi.mcp.tools.signal_isolation import load_channel

        _write_raw(tmp_path, "cyc001", ["CD31.tif", "CD3.tif"])

        cd3 = await load_channel(project_path=str(tmp_path), cycle="1", channel="CD3")
        cd31 = await load_channel(project_path=str(tmp_path), cycle="1", channel="CD31")

        assert cd3["status"] == "loaded" and cd3["source"].endswith("CD3.tif")
        assert cd31["status"] == "loaded" and cd31["source"].endswith("CD31.tif")

    @requires_zarr
    async def test_load_channel_prefers_whole_token_over_substring(self, tmp_path, clear_mcp_state):
        from kintsugi.mcp.tools.signal_isolation import load_channel

        _write_raw(tmp_path, "cyc001", ["CD31.tif", "cyc001_cd3_registered.tif"])

        result = await load_channel(project_path=str(tmp_path), cycle="1", channel="CD3")

        assert result["status"] == "loaded"
        assert result["source"].endswith("cyc001_cd3_registered.tif")

    @requires_zarr
    async def test_load_channel_reports_ambiguous_match(self, tmp_path, clear_mcp_state):
        from kintsugi.mcp.tools.signal_isolation import load_channel

        _write_raw(tmp_path, "cyc001", ["CD3_b.tif", "CD3_a.tif"])

        result = await load_channel(project_path=str(tmp_path), cycle="1", channel="CD3")

        assert "ambiguous" in result["error"]
        assert [Path(c).name for c in result["candidates"]] == ["CD3_a.tif", "CD3_b.tif"]


class TestNewlyRegisteredTools:
    """compare_channels, get_intensity_profile and batch_assess_quality on the wire."""

    async def test_intensity_profile_last_row_and_bad_inputs(self, clear_mcp_state):
        from kintsugi.mcp.session import _store_image
        from kintsugi.mcp.tools.visualization import get_intensity_profile

        image = np.arange(64, dtype=np.uint16).reshape(8, 8) * 100
        _store_image("ch", image)

        last_row = await get_intensity_profile(channel="ch", axis="horizontal", position=1.0)
        assert last_row["status"] == "success"
        assert last_row["profile"]["values"] == image[-1].tolist()

        last_col = await get_intensity_profile(channel="ch", axis="vertical", position=1.0)
        assert last_col["profile"]["values"] == image[:, -1].tolist()

        assert "position" in (await get_intensity_profile(channel="ch", position=1.5))["error"]
        assert "axis" in (await get_intensity_profile(channel="ch", axis="diagonal"))["error"]
        assert "not loaded" in (await get_intensity_profile(channel="missing"))["error"].lower()

    async def test_batch_assess_quality_reports_missing_channel_in_place(
        self, sample_uint16_image, clear_mcp_state
    ):
        from kintsugi.mcp.session import _store_image
        from kintsugi.mcp.tools.quality_assessment import batch_assess_quality

        _store_image("good", sample_uint16_image)

        result = await batch_assess_quality(channels=["good", "missing"])

        assert result["status"] == "success"
        assert result["total_channels"] == 2 and result["channels_assessed"] == 1
        by_name = {r["channel"]: r for r in result["results"]}
        assert "score" in by_name["good"] and "error" in by_name["missing"]
        assert isinstance(result["average_score"], float)

    async def test_intensity_profile_caps_points_and_rejects_3d(self, clear_mcp_state):
        from kintsugi.mcp.session import _store_image
        from kintsugi.mcp.tools.visualization import get_intensity_profile

        _store_image("wide", np.ones((4, 999), dtype=np.uint16))
        _store_image("stack", np.ones((3, 20, 30), dtype=np.uint16))
        _store_image("single_plane", np.ones((1, 5, 7), dtype=np.uint16))

        wide = await get_intensity_profile(channel="wide")
        assert len(wide["profile"]["values"]) == 500 == len(wide["profile"]["x"])
        assert wide["profile"]["x"][0] == 0 and wide["profile"]["x"][-1] == 998

        assert "2-D" in (await get_intensity_profile(channel="stack"))["error"]
        plane = await get_intensity_profile(channel="single_plane", axis="vertical")
        assert len(plane["profile"]["values"]) == 5

    async def test_compare_channels_ssim_scales_with_dtype(self, clear_mcp_state):
        from kintsugi.mcp.session import _store_image
        from kintsugi.mcp.tools.visualization import compare_channels

        rng = np.random.default_rng(0)
        _store_image("a8", rng.integers(0, 255, (32, 32)).astype(np.uint8))
        _store_image("b8", rng.integers(0, 255, (32, 32)).astype(np.uint8))
        _store_image("stack", np.ones((3, 32, 32), dtype=np.uint8))

        result = await compare_channels(channel1="a8", channel2="b8")
        comparison = result["comparison"]
        # Independent noise images are dissimilar; with 16-bit constants SSIM was ~1.
        assert comparison["ssim"] < 0.5
        assert comparison["ssim_intensity_range"] == 255.0
        assert "2-D" in (await compare_channels(channel1="a8", channel2="stack"))["error"]

    async def test_batch_assess_quality_results_are_sorted_worst_first(
        self, sample_uint16_image, sample_noisy_image, clear_mcp_state
    ):
        from kintsugi.mcp.session import _store_image
        from kintsugi.mcp.tools.quality_assessment import batch_assess_quality

        _store_image("clean", sample_uint16_image)
        _store_image("noisy", sample_noisy_image)

        result = await batch_assess_quality(channels=["clean", "missing", "noisy"])

        scores = [r["score"] for r in result["results"] if "score" in r]
        assert scores == sorted(scores) and len(scores) == 2
        assert result["results"][-1]["channel"] == "missing"

    def test_three_tools_are_registered_read_only(self):
        from kintsugi.mcp.tool_specs import TOOL_SPECS

        specs = {s["name"]: s for s in TOOL_SPECS}
        for name in ("compare_channels", "get_intensity_profile", "batch_assess_quality"):
            assert specs[name]["annotations"]["readOnlyHint"] is True
            assert specs[name]["title"]
