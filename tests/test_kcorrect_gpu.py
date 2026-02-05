"""
Tests for the GPU-accelerated BaSiC illumination correction module.

This module tests the KCorrectGPU class which implements the BaSiC algorithm
for flatfield and darkfield correction. Includes regression tests for the
Jan 2026 fix for negative flatfield values.
"""

import numpy as np
import pytest

from kintsugi.kcorrect_gpu import (
    CUPY_AVAILABLE,
    KCorrectGPU,
    KCorrectGPUFunc,
    check_gpu,
    get_available_gpus,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def tile_images():
    """Create a stack of tile images with illumination bias."""
    np.random.seed(42)
    n_images = 25
    height, width = 128, 128

    # Create illumination pattern (vignetting - darker at edges)
    y, x = np.ogrid[:height, :width]
    center_y, center_x = height / 2, width / 2
    distance = np.sqrt((y - center_y) ** 2 + (x - center_x) ** 2)
    max_dist = np.sqrt(center_y**2 + center_x**2)

    # Flatfield: center is brighter than edges
    flatfield_pattern = 1.0 - 0.3 * (distance / max_dist)

    # Create images with this illumination pattern
    images = np.zeros((n_images, height, width), dtype=np.float64)
    for i in range(n_images):
        # Random tissue signal
        signal = np.random.randint(2000, 8000, (height, width)).astype(np.float64)
        # Apply illumination pattern
        images[i] = signal * flatfield_pattern
        # Add some noise
        images[i] += np.random.randn(height, width) * 100

    return images.astype(np.float64)


@pytest.fixture
def uniform_images():
    """Create images with uniform content (for flatfield estimation)."""
    np.random.seed(42)
    n_images = 20
    height, width = 64, 64

    # Create images with slight illumination variation but uniform content
    images = np.zeros((n_images, height, width), dtype=np.float64)

    # Create a simple illumination gradient (top-left bright, bottom-right dim)
    y, x = np.ogrid[:height, :width]
    gradient = 1.0 - 0.2 * (y / height + x / width)

    for i in range(n_images):
        # Uniform base intensity with noise
        base = 5000 + np.random.randn(height, width) * 500
        images[i] = base * gradient

    return images.astype(np.float64)


@pytest.fixture
def low_signal_images():
    """Create low-signal images (sparse markers or blank channels)."""
    np.random.seed(42)
    n_images = 15
    height, width = 64, 64

    # Very low signal - mostly background
    images = np.random.randint(50, 200, (n_images, height, width)).astype(np.float64)

    # Add a few bright spots (sparse signal)
    for i in range(n_images):
        # Random bright spots
        n_spots = np.random.randint(1, 5)
        for _ in range(n_spots):
            cy = np.random.randint(10, height - 10)
            cx = np.random.randint(10, width - 10)
            images[i, cy - 2 : cy + 2, cx - 2 : cx + 2] = 2000

    return images


@pytest.fixture
def high_intensity_images():
    """Create high-intensity images."""
    np.random.seed(42)
    n_images = 15
    height, width = 64, 64

    # High intensity throughout
    images = np.random.randint(30000, 50000, (n_images, height, width)).astype(np.float64)

    return images


# =============================================================================
# Tests for KCorrectGPU class
# =============================================================================


class TestKCorrectGPU:
    """Tests for the KCorrectGPU class."""

    def test_initialization_cpu(self):
        """Test initialization with CPU backend."""
        corrector = KCorrectGPU(use_gpu=False)

        assert corrector.device == "cpu"
        assert corrector.xp == np

    def test_initialization_auto(self):
        """Test initialization with auto backend selection."""
        corrector = KCorrectGPU(use_gpu="auto")

        # Should select GPU if available, else CPU
        assert corrector.device in ["gpu", "cpu"]

    def test_working_size_parameter(self):
        """Test working_size parameter is stored."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)

        assert corrector.working_size == 64

    def test_fit_returns_flatfield_darkfield(self, tile_images):
        """Test that fit returns flatfield and darkfield arrays."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, if_darkfield=True)

        # Should return arrays of correct shape
        assert flatfield.shape == tile_images.shape[1:]
        assert darkfield.shape == tile_images.shape[1:]

    def test_fit_without_darkfield(self, tile_images):
        """Test fit with darkfield disabled."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, if_darkfield=False)

        # Darkfield should be zeros
        np.testing.assert_array_equal(darkfield, np.zeros_like(darkfield))

    def test_flatfield_normalized(self, tile_images):
        """Test that flatfield is normalized (mean ~1.0)."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, _ = corrector.fit(tile_images, max_iterations=50)

        # Flatfield should have mean close to 1.0
        assert 0.9 < flatfield.mean() < 1.1

    def test_flatfield_positive(self, tile_images):
        """Test that flatfield values are positive."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, _ = corrector.fit(tile_images, max_iterations=50)

        assert flatfield.min() > 0

    def test_transform_applies_correction(self, tile_images):
        """Test that transform applies flatfield/darkfield correction."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, max_iterations=50)

        corrected = corrector.transform(tile_images, flatfield, darkfield)

        # Corrected images should have same shape
        assert corrected.shape == tile_images.shape
        # Should be non-negative
        assert corrected.min() >= 0

    def test_transform_single_image(self, tile_images):
        """Test transform on a single 2D image."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, max_iterations=50)

        # Transform single image
        single_image = tile_images[0]
        corrected = corrector.transform(single_image, flatfield, darkfield)

        # Should return 2D array
        assert corrected.ndim == 2
        assert corrected.shape == single_image.shape

    def test_transform_flatfield_min_parameter(self, tile_images):
        """Test that flatfield_min prevents extreme amplification."""
        corrector = KCorrectGPU(use_gpu=False, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, max_iterations=50)

        # Create a flatfield with very low values
        flatfield_low = flatfield.copy()
        flatfield_low[flatfield_low < 0.5] = 0.01  # Very low values

        # With flatfield_min=0.1, amplification limited to 10x
        corrected = corrector.transform(
            tile_images, flatfield_low, darkfield, flatfield_min=0.1
        )

        # Should not have extreme amplification
        assert corrected.max() < tile_images.max() * 15

    def test_reduced_iterations(self, uniform_images):
        """Test with reduced iterations for speed."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32)
        flatfield, darkfield = corrector.fit(
            uniform_images, max_iterations=10, max_reweight_iterations=2
        )

        # Should still produce valid output
        assert flatfield.shape == uniform_images.shape[1:]
        assert flatfield.min() > 0


# =============================================================================
# Tests for KCorrectGPUFunc functional interface
# =============================================================================


class TestKCorrectGPUFunc:
    """Tests for the functional interface."""

    def test_basic_usage(self, tile_images):
        """Test basic functional interface usage."""
        flatfield, darkfield = KCorrectGPUFunc(
            tile_images, use_gpu=False, max_iterations=20
        )

        assert flatfield.shape == tile_images.shape[1:]
        assert darkfield.shape == tile_images.shape[1:]

    def test_list_input(self, tile_images):
        """Test with list of images as input."""
        images_list = [tile_images[i] for i in range(tile_images.shape[0])]

        flatfield, darkfield = KCorrectGPUFunc(
            images_list, use_gpu=False, max_iterations=20
        )

        assert flatfield.shape == tile_images.shape[1:]

    def test_transposed_input(self):
        """Test handling of transposed input (height, width, n_images)."""
        np.random.seed(42)
        # Shape (height, width, n_images) - wrong order
        images = np.random.randint(1000, 5000, (64, 64, 15)).astype(np.float64)

        flatfield, darkfield = KCorrectGPUFunc(
            images, use_gpu=False, max_iterations=20
        )

        # Should handle the transpose and return correct shape
        # After transpose, images become (15, 64, 64)
        assert flatfield.shape == (64, 64)


# =============================================================================
# Regression tests for negative flatfield fix (Jan 2026)
# =============================================================================


class TestNegativeFlatfieldFix:
    """
    Regression tests for the negative flatfield fix.

    The BaSiC algorithm can produce negative flatfield values for low-signal
    images (sparse markers, blank channels). The fix detects this and falls
    back to a uniform flatfield (no correction).
    """

    def test_low_signal_images_no_crash(self, low_signal_images):
        """Test that low-signal images don't crash the algorithm."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32, verbose=False)

        # Should not crash
        flatfield, darkfield = corrector.fit(
            low_signal_images, max_iterations=20, max_reweight_iterations=3
        )

        assert flatfield.shape == low_signal_images.shape[1:]

    def test_low_signal_flatfield_positive(self, low_signal_images):
        """Test that flatfield is always positive even for low-signal data."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32, verbose=False)

        flatfield, _ = corrector.fit(
            low_signal_images, max_iterations=20, max_reweight_iterations=3
        )

        # Flatfield should never be negative (the fix ensures this)
        assert flatfield.min() > 0, "Flatfield has negative values - fix may have regressed"

    def test_low_signal_uses_uniform_flatfield(self, low_signal_images):
        """Test that very low signal images may use uniform flatfield."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32, verbose=False)

        flatfield, _ = corrector.fit(
            low_signal_images, max_iterations=20, max_reweight_iterations=3
        )

        # For very low signal, the fix may produce uniform flatfield
        # Check that it's either uniform (all 1s) or has reasonable variation
        flatfield_std = flatfield.std()

        # Either uniform (std near 0) or reasonable variation (std < 0.5)
        assert flatfield_std < 0.5, f"Flatfield has unusual variation: std={flatfield_std}"

    def test_transform_with_clamped_flatfield(self, low_signal_images):
        """Test that transform works with potentially clamped flatfield."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32, verbose=False)

        flatfield, darkfield = corrector.fit(
            low_signal_images, max_iterations=20, max_reweight_iterations=3
        )

        # Transform should work without issues
        corrected = corrector.transform(low_signal_images, flatfield, darkfield)

        assert corrected.shape == low_signal_images.shape
        assert corrected.min() >= 0


# =============================================================================
# Tests for utility functions
# =============================================================================


class TestUtilityFunctions:
    """Tests for utility functions."""

    def test_check_gpu_returns_tuple(self):
        """Test that check_gpu returns (bool, str) tuple."""
        available, message = check_gpu()

        assert isinstance(available, bool)
        assert isinstance(message, str)

    def test_get_available_gpus_returns_list(self):
        """Test that get_available_gpus returns a list."""
        gpus = get_available_gpus()

        assert isinstance(gpus, list)
        # All elements should be integers
        for gpu_id in gpus:
            assert isinstance(gpu_id, int)

    def test_cupy_available_flag(self):
        """Test that CUPY_AVAILABLE is a boolean."""
        assert isinstance(CUPY_AVAILABLE, bool)


# =============================================================================
# Tests for GPU backend (conditional)
# =============================================================================


@pytest.mark.skipif(not CUPY_AVAILABLE, reason="CuPy not available")
class TestGPUBackend:
    """Tests for GPU backend (only run if CuPy available)."""

    def test_gpu_initialization(self):
        """Test GPU initialization."""
        corrector = KCorrectGPU(use_gpu=True)

        assert corrector.device == "gpu"

    def test_gpu_fit(self, tile_images):
        """Test fit on GPU."""
        corrector = KCorrectGPU(use_gpu=True, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, max_iterations=20)

        assert flatfield.shape == tile_images.shape[1:]
        # Result should be numpy array (transferred back from GPU)
        assert isinstance(flatfield, np.ndarray)

    def test_gpu_transform(self, tile_images):
        """Test transform on GPU."""
        corrector = KCorrectGPU(use_gpu=True, working_size=64)
        flatfield, darkfield = corrector.fit(tile_images, max_iterations=20)
        corrected = corrector.transform(tile_images, flatfield, darkfield)

        assert corrected.shape == tile_images.shape
        assert isinstance(corrected, np.ndarray)

    def test_gpu_cpu_consistency(self, uniform_images):
        """Test that GPU and CPU produce similar results."""
        # CPU
        corrector_cpu = KCorrectGPU(use_gpu=False, working_size=32)
        flatfield_cpu, darkfield_cpu = corrector_cpu.fit(
            uniform_images, max_iterations=50, max_reweight_iterations=5
        )

        # GPU
        corrector_gpu = KCorrectGPU(use_gpu=True, working_size=32)
        flatfield_gpu, darkfield_gpu = corrector_gpu.fit(
            uniform_images, max_iterations=50, max_reweight_iterations=5
        )

        # Results should be similar (allow for numerical differences)
        np.testing.assert_allclose(flatfield_cpu, flatfield_gpu, rtol=0.1, atol=0.05)


# =============================================================================
# Edge case tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_single_image(self):
        """Test with single image (minimum input)."""
        np.random.seed(42)
        images = np.random.randint(1000, 5000, (1, 64, 64)).astype(np.float64)

        corrector = KCorrectGPU(use_gpu=False, working_size=32)

        # Single image is unusual but shouldn't crash
        flatfield, darkfield = corrector.fit(images, max_iterations=10)

        assert flatfield.shape == (64, 64)

    def test_very_small_images(self):
        """Test with very small images."""
        np.random.seed(42)
        images = np.random.randint(1000, 5000, (10, 16, 16)).astype(np.float64)

        corrector = KCorrectGPU(use_gpu=False, working_size=16)
        flatfield, darkfield = corrector.fit(images, max_iterations=10)

        assert flatfield.shape == (16, 16)

    def test_constant_images(self):
        """Test with constant-value images."""
        images = np.full((10, 64, 64), 5000.0, dtype=np.float64)
        # Add tiny noise to avoid division issues
        images += np.random.randn(10, 64, 64) * 0.1

        corrector = KCorrectGPU(use_gpu=False, working_size=32)
        flatfield, darkfield = corrector.fit(images, max_iterations=10)

        # Should produce roughly uniform flatfield
        assert flatfield.std() < 0.5

    def test_high_contrast_images(self, high_intensity_images):
        """Test with high-intensity images."""
        corrector = KCorrectGPU(use_gpu=False, working_size=32)
        flatfield, darkfield = corrector.fit(high_intensity_images, max_iterations=20)

        assert flatfield.shape == high_intensity_images.shape[1:]
        assert flatfield.min() > 0

    def test_different_dtypes(self):
        """Test with different input dtypes."""
        np.random.seed(42)

        for dtype in [np.float32, np.float64, np.uint16]:
            images = np.random.randint(1000, 5000, (10, 32, 32)).astype(dtype)

            corrector = KCorrectGPU(use_gpu=False, working_size=32)
            flatfield, darkfield = corrector.fit(images, max_iterations=10)

            # Should work regardless of input dtype
            assert flatfield.dtype == np.float64
            assert flatfield.shape == (32, 32)


# =============================================================================
# Tests for DCT operations
# =============================================================================


class TestDCTOperations:
    """Tests for internal DCT operations."""

    def test_dct_idct_roundtrip(self):
        """Test that DCT followed by IDCT recovers original."""
        corrector = KCorrectGPU(use_gpu=False)

        # Create test matrix
        np.random.seed(42)
        original = np.random.randn(32, 32).astype(np.float64)

        # DCT then IDCT
        dct = corrector._dct2d(original)
        recovered = corrector._idct2d(dct)

        # Should recover original
        np.testing.assert_allclose(original, recovered, rtol=1e-10)

    def test_dct_shape_preserved(self):
        """Test that DCT preserves shape."""
        corrector = KCorrectGPU(use_gpu=False)

        matrix = np.random.randn(64, 64).astype(np.float64)
        dct = corrector._dct2d(matrix)

        assert dct.shape == matrix.shape

    def test_dct_rejects_3d_input(self):
        """Test that DCT raises on 3D input."""
        corrector = KCorrectGPU(use_gpu=False)

        matrix_3d = np.random.randn(10, 32, 32).astype(np.float64)

        with pytest.raises(ValueError, match="2D"):
            corrector._dct2d(matrix_3d)


# =============================================================================
# Tests for shrinkage operator
# =============================================================================


class TestShrinkageOperator:
    """Tests for the soft thresholding shrinkage operator."""

    def test_shrinkage_zero_epsilon(self):
        """Test shrinkage with zero threshold."""
        corrector = KCorrectGPU(use_gpu=False)

        matrix = np.array([[-2, -1, 0, 1, 2]], dtype=np.float64)
        result = corrector._shrinkage_operator(matrix, 0)

        # With zero threshold, should return original
        np.testing.assert_array_equal(result, matrix)

    def test_shrinkage_positive_values(self):
        """Test shrinkage on positive values."""
        corrector = KCorrectGPU(use_gpu=False)

        matrix = np.array([[0, 1, 2, 3, 4]], dtype=np.float64)
        result = corrector._shrinkage_operator(matrix, 1.5)

        # Values <= threshold become 0, others reduced by threshold
        expected = np.array([[0, 0, 0.5, 1.5, 2.5]], dtype=np.float64)
        np.testing.assert_array_almost_equal(result, expected)

    def test_shrinkage_negative_values(self):
        """Test shrinkage on negative values."""
        corrector = KCorrectGPU(use_gpu=False)

        matrix = np.array([[-4, -3, -2, -1, 0]], dtype=np.float64)
        result = corrector._shrinkage_operator(matrix, 1.5)

        # Negative values: |value| <= threshold become 0, others increased by threshold
        expected = np.array([[-2.5, -1.5, -0.5, 0, 0]], dtype=np.float64)
        np.testing.assert_array_almost_equal(result, expected)
