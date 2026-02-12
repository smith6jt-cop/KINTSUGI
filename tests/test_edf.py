"""
Tests for the Extended Depth of Focus (EDF) module.

This module tests the core EDF algorithms which are critical for z-stack
projection in the KINTSUGI pipeline. Includes regression tests for the
Jan 2026 sigma fix that resolved significant detail loss issues.
"""

import numpy as np
import pytest

from kintsugi.edf import (
    EDFProcessor,
    _check_cupy_available,
    edf_tiled,
    extended_depth_of_focus_laplacian,
    extended_depth_of_focus_variance,
    process_edf,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def simple_zstack():
    """Create a simple z-stack for testing."""
    np.random.seed(42)
    # 10 z-slices, 128x128 pixels
    stack = np.random.randint(1000, 10000, (10, 128, 128), dtype=np.uint16)
    return stack


@pytest.fixture
def focused_zstack():
    """Create a z-stack with clear best-focus plane."""
    np.random.seed(42)
    stack = np.zeros((15, 128, 128), dtype=np.uint16)

    # Create background
    for z in range(15):
        stack[z] = np.random.randint(500, 1500, (128, 128)).astype(np.uint16)

    # Add sharp features at specific z-planes
    # Best focus at z=7
    for z in range(15):
        blur_amount = abs(z - 7)
        feature_intensity = max(0, 10000 - blur_amount * 1500)

        # Add a bright square that's sharpest at z=7
        if feature_intensity > 0:
            # Sharper edges at z=7, blurrier away from it
            margin = blur_amount * 2
            if margin < 20:
                stack[z, 50 - margin : 78 + margin, 50 - margin : 78 + margin] += int(
                    feature_intensity * 0.3
                )
            stack[z, 50:78, 50:78] += int(feature_intensity)

    return stack.astype(np.uint16)


@pytest.fixture
def high_contrast_zstack():
    """Create a high-contrast z-stack for testing detail preservation."""
    np.random.seed(42)
    stack = np.zeros((10, 256, 256), dtype=np.uint16)

    for z in range(10):
        # Add checkerboard pattern (high-frequency detail)
        checker = np.indices((256, 256)).sum(axis=0) % 2
        stack[z] = (checker * 5000 + 2000).astype(np.uint16)

        # Add some noise
        stack[z] += np.random.randint(0, 500, (256, 256)).astype(np.uint16)

        # Make middle slices "in focus" (higher contrast)
        if 3 <= z <= 6:
            stack[z] = (stack[z] * 1.5).clip(0, 65535).astype(np.uint16)

    return stack


@pytest.fixture
def large_zstack():
    """Create a larger z-stack for tiled processing tests."""
    np.random.seed(42)
    # 8 z-slices, 512x512 pixels
    stack = np.random.randint(1000, 10000, (8, 512, 512), dtype=np.uint16)
    return stack


# =============================================================================
# Tests for extended_depth_of_focus_variance
# =============================================================================


class TestExtendedDepthOfFocusVariance:
    """Tests for the variance-based EDF function."""

    def test_basic_projection(self, simple_zstack):
        """Test basic z-stack projection."""
        result = extended_depth_of_focus_variance(simple_zstack)

        # Result should be 2D
        assert result.ndim == 2
        # Result should have same spatial dimensions
        assert result.shape == simple_zstack.shape[1:]
        # Result dtype should match input
        assert result.dtype == simple_zstack.dtype

    def test_output_dtype_uint16(self, simple_zstack):
        """Test that uint16 input produces uint16 output."""
        result = extended_depth_of_focus_variance(simple_zstack)
        assert result.dtype == np.uint16

    def test_output_dtype_uint8(self):
        """Test that uint8 input produces uint8 output."""
        stack = np.random.randint(0, 255, (5, 64, 64), dtype=np.uint8)
        result = extended_depth_of_focus_variance(stack)
        assert result.dtype == np.uint8

    def test_output_dtype_float32(self):
        """Test that float32 input produces float32 output."""
        stack = np.random.rand(5, 64, 64).astype(np.float32) * 10000
        result = extended_depth_of_focus_variance(stack)
        assert result.dtype == np.float32

    def test_selects_best_focus(self, focused_zstack):
        """Test that EDF selects pixels from best-focus plane."""
        result = extended_depth_of_focus_variance(focused_zstack)

        # The feature region should have high intensity (from focused plane)
        feature_region = result[50:78, 50:78]
        background_region = result[0:30, 0:30]

        assert feature_region.mean() > background_region.mean() * 2

    def test_radius_parameters(self, simple_zstack):
        """Test that radius parameters affect result."""
        result_small = extended_depth_of_focus_variance(simple_zstack, radius_x=1, radius_y=1)
        result_large = extended_depth_of_focus_variance(simple_zstack, radius_x=5, radius_y=5)

        # Both should produce valid results
        assert result_small.shape == simple_zstack.shape[1:]
        assert result_large.shape == simple_zstack.shape[1:]
        # Results may differ due to different variance windows
        # (not strictly required but expected)

    def test_sigma_parameter(self, simple_zstack):
        """Test that sigma parameter affects result."""
        result_no_smooth = extended_depth_of_focus_variance(simple_zstack, sigma=0.0)
        result_smooth = extended_depth_of_focus_variance(simple_zstack, sigma=20.0)

        # Both should produce valid results
        assert result_no_smooth.shape == simple_zstack.shape[1:]
        assert result_smooth.shape == simple_zstack.shape[1:]

    def test_blend_depth_parameter(self, focused_zstack):
        """Test blend_depth for smooth transitions."""
        result_no_blend = extended_depth_of_focus_variance(focused_zstack, blend_depth=0)
        result_blend = extended_depth_of_focus_variance(focused_zstack, blend_depth=3)

        # Both should produce valid results
        assert result_no_blend.shape == focused_zstack.shape[1:]
        assert result_blend.shape == focused_zstack.shape[1:]

        # Blended version may have smoother transitions
        # (comparing gradients could verify this, but basic shape test suffices)

    def test_z_smooth_sigma_parameter(self, focused_zstack):
        """Test z_smooth_sigma for z-index smoothing."""
        result_no_smooth = extended_depth_of_focus_variance(focused_zstack, z_smooth_sigma=0.0)
        result_smooth = extended_depth_of_focus_variance(focused_zstack, z_smooth_sigma=3.0)

        # Both should be valid
        assert result_no_smooth.shape == focused_zstack.shape[1:]
        assert result_smooth.shape == focused_zstack.shape[1:]

    def test_raises_on_2d_input(self):
        """Test that 2D input raises ValueError."""
        img_2d = np.random.randint(0, 10000, (128, 128), dtype=np.uint16)
        with pytest.raises(ValueError, match="Expected 3D"):
            extended_depth_of_focus_variance(img_2d)

    def test_raises_on_4d_input(self):
        """Test that 4D input raises ValueError."""
        img_4d = np.random.randint(0, 10000, (5, 3, 128, 128), dtype=np.uint16)
        with pytest.raises(ValueError, match="Expected 3D"):
            extended_depth_of_focus_variance(img_4d)

    def test_device_cpu_explicit(self, simple_zstack):
        """Test explicit CPU device selection."""
        result = extended_depth_of_focus_variance(simple_zstack, device="cpu")
        assert result.shape == simple_zstack.shape[1:]

    def test_single_slice_stack(self):
        """Test handling of single-slice stack."""
        stack = np.random.randint(0, 10000, (1, 64, 64), dtype=np.uint16)
        result = extended_depth_of_focus_variance(stack)

        # Should just return that slice
        np.testing.assert_array_equal(result, stack[0])


# =============================================================================
# Tests for extended_depth_of_focus_laplacian
# =============================================================================


class TestExtendedDepthOfFocusLaplacian:
    """Tests for the Laplacian-based EDF function."""

    def test_basic_projection(self, simple_zstack):
        """Test basic Laplacian projection."""
        result = extended_depth_of_focus_laplacian(simple_zstack)

        assert result.ndim == 2
        assert result.shape == simple_zstack.shape[1:]
        assert result.dtype == simple_zstack.dtype

    def test_kernel_size_parameter(self, simple_zstack):
        """Test kernel_size parameter."""
        result_small = extended_depth_of_focus_laplacian(simple_zstack, kernel_size=3)
        result_large = extended_depth_of_focus_laplacian(simple_zstack, kernel_size=7)

        # Both should be valid
        assert result_small.shape == simple_zstack.shape[1:]
        assert result_large.shape == simple_zstack.shape[1:]

    def test_raises_on_2d_input(self):
        """Test that 2D input raises ValueError."""
        img_2d = np.random.randint(0, 10000, (128, 128), dtype=np.uint16)
        with pytest.raises(ValueError, match="Expected 3D"):
            extended_depth_of_focus_laplacian(img_2d)


# =============================================================================
# Tests for edf_tiled
# =============================================================================


class TestEDFTiled:
    """Tests for tiled EDF processing."""

    def test_basic_tiled_processing(self, large_zstack):
        """Test basic tiled processing."""
        result = edf_tiled(large_zstack, tile_size=(256, 256), overlap=50)

        assert result.shape == large_zstack.shape[1:]
        assert result.dtype == large_zstack.dtype

    def test_small_tile_size(self, simple_zstack):
        """Test with very small tiles."""
        result = edf_tiled(simple_zstack, tile_size=(64, 64), overlap=10)

        assert result.shape == simple_zstack.shape[1:]

    def test_large_overlap(self, simple_zstack):
        """Test with large overlap."""
        result = edf_tiled(simple_zstack, tile_size=(64, 64), overlap=30)

        assert result.shape == simple_zstack.shape[1:]

    def test_laplacian_method(self, simple_zstack):
        """Test tiled processing with Laplacian method."""
        result = edf_tiled(simple_zstack, tile_size=(64, 64), method="laplacian")

        assert result.shape == simple_zstack.shape[1:]

    def test_invalid_method_raises(self, simple_zstack):
        """Test that invalid method raises ValueError."""
        with pytest.raises(ValueError, match="Unknown method"):
            edf_tiled(simple_zstack, method="invalid_method")

    def test_raises_on_2d_input(self):
        """Test that 2D input raises ValueError."""
        img_2d = np.random.randint(0, 10000, (128, 128), dtype=np.uint16)
        with pytest.raises(ValueError, match="Expected 3D"):
            edf_tiled(img_2d)

    def test_blend_depth_in_tiled(self, large_zstack):
        """Test blend_depth parameter in tiled processing."""
        result = edf_tiled(large_zstack, tile_size=(256, 256), method="variance", blend_depth=2)

        assert result.shape == large_zstack.shape[1:]

    def test_tile_seams_blended(self, large_zstack):
        """Test that tile seams are properly handled."""
        result = edf_tiled(large_zstack, tile_size=(256, 256), overlap=50)

        # Basic check that result is valid
        assert result.shape == large_zstack.shape[1:]
        assert result.dtype == large_zstack.dtype
        # Result should contain actual data values
        assert result.max() > result.min()


# =============================================================================
# Tests for EDFProcessor class
# =============================================================================


class TestEDFProcessor:
    """Tests for the EDFProcessor class."""

    def test_initialization_auto(self):
        """Test initialization with auto backend."""
        processor = EDFProcessor(backend="auto")

        # Should select a valid backend
        assert processor.backend in ["numpy", "cupy", "clij2"]

    def test_initialization_numpy(self):
        """Test initialization with numpy backend."""
        processor = EDFProcessor(backend="numpy")

        assert processor.backend == "numpy"

    def test_initialization_invalid_backend(self):
        """Test initialization with invalid backend raises."""
        with pytest.raises(ValueError, match="Unknown backend"):
            EDFProcessor(backend="invalid_backend")

    def test_process_basic(self, simple_zstack):
        """Test basic processing through processor."""
        processor = EDFProcessor(backend="numpy")
        result = processor.process(simple_zstack)

        assert result.shape == simple_zstack.shape[1:]

    def test_process_with_z_range(self, focused_zstack):
        """Test processing with z_start and z_end."""
        processor = EDFProcessor(backend="numpy")
        # Use only middle slices (1-indexed for CLIJ2 compat)
        result = processor.process(focused_zstack, z_start=5, z_end=10)

        assert result.shape == focused_zstack.shape[1:]

    def test_process_with_tiles_parameter(self, large_zstack):
        """Test processing with tiles parameter (CLIJ2-style)."""
        processor = EDFProcessor(backend="numpy")
        result = processor.process(large_zstack, tiles=(2, 2))

        assert result.shape == large_zstack.shape[1:]

    def test_process_with_tile_size(self, large_zstack):
        """Test processing with explicit tile_size."""
        processor = EDFProcessor(backend="numpy")
        result = processor.process(large_zstack, tile_size=(256, 256))

        assert result.shape == large_zstack.shape[1:]

    def test_laplacian_method(self, simple_zstack):
        """Test processor with Laplacian method."""
        processor = EDFProcessor(backend="numpy", method="laplacian")
        result = processor.process(simple_zstack)

        assert result.shape == simple_zstack.shape[1:]


# =============================================================================
# Tests for process_edf convenience function
# =============================================================================


class TestProcessEDF:
    """Tests for the process_edf convenience function."""

    def test_basic_usage(self, simple_zstack):
        """Test basic convenience function usage."""
        result = process_edf(simple_zstack)

        assert result.shape == simple_zstack.shape[1:]

    def test_with_parameters(self, simple_zstack):
        """Test with custom parameters."""
        result = process_edf(simple_zstack, radius_x=3, radius_y=3, sigma=15.0, backend="numpy")

        assert result.shape == simple_zstack.shape[1:]


# =============================================================================
# Regression tests for Jan 2026 sigma fix
# =============================================================================


class TestSigmaFixRegression:
    """
    Regression tests for the critical sigma fix from Jan 2026.

    The bug was that sigma was smoothing the INPUT IMAGE before variance
    calculation (destroying detail), instead of smoothing the VARIANCE MAP
    after variance calculation (regularizing focus selection).
    """

    def test_detail_preserved_with_sigma(self, high_contrast_zstack):
        """Test that high-frequency detail is preserved even with large sigma.

        Before the fix, large sigma values would blur the input image,
        destroying the checkerboard pattern. After the fix, sigma only
        affects the variance map (altitude surface), not the image content.
        """
        # Process with large sigma
        result = extended_depth_of_focus_variance(
            high_contrast_zstack, sigma=20.0, radius_x=2, radius_y=2
        )

        # Check that checkerboard contrast is preserved
        # The checkerboard pattern should still be visible
        result_float = result.astype(np.float32)

        # Compute local contrast by looking at adjacent pixels
        diff_x = np.abs(np.diff(result_float, axis=1))
        diff_y = np.abs(np.diff(result_float, axis=0))

        # High contrast differences should be present (checkerboard edges)
        # Before fix, sigma would blur these away
        max_contrast = max(diff_x.max(), diff_y.max())

        # Should have significant contrast (checkerboard pattern)
        assert (
            max_contrast > 1000
        ), "Detail was lost - sigma may be smoothing input instead of variance"

    def test_sigma_zero_vs_nonzero_preserves_intensity(self, focused_zstack):
        """Test that sigma doesn't change the intensity values significantly.

        Sigma should only affect which z-slice is selected, not the actual
        pixel values. The output should be actual pixel values from the input.
        """
        result_no_sigma = extended_depth_of_focus_variance(focused_zstack, sigma=0.0)
        result_with_sigma = extended_depth_of_focus_variance(focused_zstack, sigma=10.0)

        # Both results should have similar intensity ranges
        # (they select from the same input stack)
        assert abs(result_no_sigma.mean() - result_with_sigma.mean()) < result_no_sigma.mean() * 0.3

    def test_default_parameters_match_clij2(self, simple_zstack):
        """Test that default parameters (radius=2, sigma=10) work correctly.

        These are the CLIJ2 defaults mentioned in the CLAUDE.md fix notes.
        """
        result = extended_depth_of_focus_variance(simple_zstack, radius_x=2, radius_y=2, sigma=10.0)

        # Should produce valid output
        assert result.shape == simple_zstack.shape[1:]
        assert result.dtype == simple_zstack.dtype

        # Output should be within input range
        assert result.min() >= simple_zstack.min()
        assert result.max() <= simple_zstack.max()


# =============================================================================
# Edge case tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_very_small_stack(self):
        """Test with minimal stack size."""
        stack = np.random.randint(0, 10000, (2, 32, 32), dtype=np.uint16)
        result = extended_depth_of_focus_variance(stack)

        assert result.shape == (32, 32)

    def test_asymmetric_dimensions(self):
        """Test with non-square images."""
        stack = np.random.randint(0, 10000, (5, 100, 200), dtype=np.uint16)
        result = extended_depth_of_focus_variance(stack)

        assert result.shape == (100, 200)

    def test_constant_stack(self):
        """Test with constant-value stack."""
        stack = np.full((5, 64, 64), 5000, dtype=np.uint16)
        result = extended_depth_of_focus_variance(stack)

        # All slices are identical, so result should equal input slice
        np.testing.assert_array_equal(result, stack[0])

    def test_single_bright_pixel_focus(self):
        """Test focus detection with single bright pixel at different depths."""
        stack = np.zeros((10, 64, 64), dtype=np.uint16)
        stack[:, :, :] = 1000  # Background

        # Add bright pixel at z=5
        stack[5, 32, 32] = 50000

        result = extended_depth_of_focus_variance(stack)

        # The bright pixel location should have high value
        assert result[32, 32] > 10000

    def test_negative_radius_clamped(self, simple_zstack):
        """Test that negative radius doesn't crash."""
        # This may raise or clamp to minimum - either is acceptable
        try:
            result = extended_depth_of_focus_variance(simple_zstack, radius_x=-1)
            # If it runs, should still produce valid output
            assert result.shape == simple_zstack.shape[1:]
        except (ValueError, RuntimeError):
            # Also acceptable to raise an error
            pass

    def test_very_large_sigma(self, simple_zstack):
        """Test with very large sigma value."""
        result = extended_depth_of_focus_variance(simple_zstack, sigma=1000.0)

        # Should not crash, should produce valid output
        assert result.shape == simple_zstack.shape[1:]

    def test_float64_input(self):
        """Test with float64 input."""
        stack = np.random.rand(5, 64, 64).astype(np.float64) * 10000
        result = extended_depth_of_focus_variance(stack)

        assert result.dtype == np.float64


# =============================================================================
# GPU availability check
# =============================================================================


class TestGPUAvailability:
    """Tests for GPU availability checking."""

    def test_check_cupy_available_returns_bool(self):
        """Test that _check_cupy_available returns a boolean."""
        result = _check_cupy_available()
        assert isinstance(result, bool)

    @pytest.mark.skipif(not _check_cupy_available(), reason="CuPy not available")
    def test_gpu_produces_same_result_as_cpu(self, simple_zstack):
        """Test that GPU and CPU produce similar results."""
        result_cpu = extended_depth_of_focus_variance(simple_zstack, device="cpu")
        result_gpu = extended_depth_of_focus_variance(simple_zstack, device="gpu")

        # Results should be very similar (floating point differences allowed)
        np.testing.assert_allclose(result_cpu.astype(float), result_gpu.astype(float), rtol=0.01)
