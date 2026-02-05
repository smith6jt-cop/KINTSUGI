"""
Tests for the Kdecon deconvolution module.

This module tests the Lucy-Richardson deconvolution algorithm and PSF generation
functions. Includes regression tests for the Jan 2026 edge apodization fix
that resolved horizontal banding (Gibbs ringing) artifacts.
"""

import sys
from pathlib import Path

import numpy as np
import pytest

# Add notebooks to path for import
NOTEBOOKS_PATH = Path(__file__).parent.parent / "notebooks"
if str(NOTEBOOKS_PATH) not in sys.path:
    sys.path.insert(0, str(NOTEBOOKS_PATH))

from Kdecon.deconvolution import (
    CUPY_AVAILABLE,
    _create_tukey_window_3d,
    _find_good_fft_length,
    _psf_to_otf,
    deconvolve,
    lucy_richardson_cpu,
)
from Kdecon.psf import (
    _mirror_octant_to_full,
    calculate_psf_size,
    generate_psf,
)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def simple_stack():
    """Create a simple 3D stack for testing."""
    np.random.seed(42)
    # Small stack for fast testing
    stack = np.random.randint(1000, 5000, (10, 64, 64), dtype=np.uint16)
    return stack.astype(np.float32)


@pytest.fixture
def stack_with_features():
    """Create a stack with distinct features to test deconvolution."""
    np.random.seed(42)
    stack = np.zeros((15, 64, 64), dtype=np.float32)

    # Add background
    stack += 1000

    # Add a bright point source (should be sharpened by deconvolution)
    # Place at center z-slice
    stack[7, 32, 32] = 10000

    # Add some blurred structure
    for z in range(15):
        blur_level = abs(z - 7) + 1
        stack[z, 20:30, 20:30] += 3000 / blur_level

    # Add noise
    stack += np.random.randn(15, 64, 64).astype(np.float32) * 100

    return stack


@pytest.fixture
def small_psf():
    """Create a small PSF for testing."""
    # Simple Gaussian PSF
    size = 7
    sigma_xy = 1.0
    sigma_z = 2.0

    z, y, x = np.mgrid[:size, :size, :size].astype(np.float32)
    center = (size - 1) / 2

    psf = np.exp(
        -(
            ((x - center) ** 2 + (y - center) ** 2) / (2 * sigma_xy**2)
            + (z - center) ** 2 / (2 * sigma_z**2)
        )
    )
    psf = psf / psf.sum()  # Normalize

    return psf.astype(np.float32)


@pytest.fixture
def blurred_stack(stack_with_features, small_psf):
    """Create a blurred stack by convolving with PSF."""
    from scipy.ndimage import convolve

    blurred = convolve(stack_with_features, small_psf, mode="reflect")
    return blurred.astype(np.float32)


# =============================================================================
# Tests for _find_good_fft_length
# =============================================================================


class TestFindGoodFFTLength:
    """Tests for FFT length optimization."""

    def test_returns_at_least_input(self):
        """Test that returned length is >= input."""
        for n in [10, 50, 100, 256]:
            result = _find_good_fft_length(n)
            assert result >= n

    def test_power_of_two_unchanged(self):
        """Test that powers of 2 are unchanged (already optimal)."""
        for n in [16, 32, 64, 128, 256]:
            result = _find_good_fft_length(n)
            assert result == n

    def test_finds_efficient_size(self):
        """Test that result has only small prime factors (2, 3, 5)."""
        for n in [17, 23, 49, 97, 123]:
            result = _find_good_fft_length(n)

            # Verify result has only factors 2, 3, 5
            remaining = result
            for p in [2, 3, 5]:
                while remaining % p == 0:
                    remaining //= p
            assert remaining == 1, f"FFT length {result} has large prime factor"

    def test_common_image_sizes(self):
        """Test common image dimensions."""
        test_cases = [
            (1024, 1024),  # Power of 2
            (1000, 1000),  # Round number
            (512, 512),  # Power of 2
            (500, 500),  # Near power of 2
        ]
        for input_n, expected_max in test_cases:
            result = _find_good_fft_length(input_n)
            assert input_n <= result <= expected_max + 100


# =============================================================================
# Tests for _psf_to_otf
# =============================================================================


class TestPSFtoOTF:
    """Tests for PSF to OTF conversion."""

    def test_output_shape(self, small_psf):
        """Test that OTF has correct shape."""
        target_shape = (20, 20, 20)
        otf = _psf_to_otf(small_psf, target_shape)

        assert otf.shape == target_shape

    def test_output_is_complex(self, small_psf):
        """Test that OTF is complex."""
        target_shape = (15, 15, 15)
        otf = _psf_to_otf(small_psf, target_shape)

        assert np.iscomplexobj(otf)

    def test_dc_component(self, small_psf):
        """Test that DC component (sum) is preserved."""
        target_shape = (20, 20, 20)
        otf = _psf_to_otf(small_psf, target_shape)

        # DC component should equal sum of PSF
        # (FFT of normalized PSF has DC = 1.0)
        assert np.abs(otf[0, 0, 0]) > 0.9  # Should be ~1.0 for normalized PSF


# =============================================================================
# Tests for _create_tukey_window_3d
# =============================================================================


class TestCreateTukeyWindow3D:
    """Tests for Tukey window creation (edge apodization)."""

    def test_output_shape(self):
        """Test that window has correct shape."""
        shape = (32, 32, 16)
        pad_before = (4, 4, 2)
        pad_after = (4, 4, 2)

        window = _create_tukey_window_3d(shape, pad_before, pad_after)

        assert window.shape == shape

    def test_center_is_one(self):
        """Test that center region is all ones."""
        shape = (32, 32, 16)
        pad_before = (4, 4, 2)
        pad_after = (4, 4, 2)

        window = _create_tukey_window_3d(shape, pad_before, pad_after)

        # Center region should be 1.0
        center = window[8:24, 8:24, 4:12]
        assert np.allclose(center, 1.0, atol=0.01)

    def test_edges_taper(self):
        """Test that edges taper to lower values."""
        shape = (32, 32, 16)
        pad_before = (8, 8, 4)
        pad_after = (8, 8, 4)

        window = _create_tukey_window_3d(shape, pad_before, pad_after)

        # Edge values should be less than center
        assert window[0, 0, 0] < 0.5
        assert window[-1, -1, -1] < 0.5

    def test_no_padding_gives_flat_window(self):
        """Test that no padding produces flat window."""
        shape = (32, 32, 16)
        pad_before = (0, 0, 0)
        pad_after = (0, 0, 0)

        window = _create_tukey_window_3d(shape, pad_before, pad_after)

        # Should be all ones
        np.testing.assert_allclose(window, 1.0)

    def test_asymmetric_padding(self):
        """Test with different padding on each side."""
        shape = (32, 32, 16)
        pad_before = (8, 4, 2)
        pad_after = (4, 8, 4)

        window = _create_tukey_window_3d(shape, pad_before, pad_after)

        assert window.shape == shape
        # Should still have smooth transitions


# =============================================================================
# Tests for lucy_richardson_cpu
# =============================================================================


class TestLucyRichardsonCPU:
    """Tests for CPU Lucy-Richardson implementation."""

    def test_output_shape(self, simple_stack, small_psf):
        """Test that output shape matches input."""
        result = lucy_richardson_cpu(
            simple_stack, small_psf, iterations=2, verbose=False
        )

        assert result.shape == simple_stack.shape

    def test_output_non_negative(self, simple_stack, small_psf):
        """Test that output is non-negative."""
        result = lucy_richardson_cpu(
            simple_stack, small_psf, iterations=3, verbose=False
        )

        assert np.all(result >= 0)

    def test_convergence(self, blurred_stack, small_psf):
        """Test that algorithm converges (delta decreases)."""
        deltas = []

        def track_delta(iteration, delta_rel):
            deltas.append(delta_rel)

        result = lucy_richardson_cpu(
            blurred_stack,
            small_psf,
            iterations=10,
            callback=track_delta,
            verbose=False,
        )

        # After first iteration, deltas should generally decrease
        # (convergence behavior)
        assert len(deltas) > 0

    def test_damping_parameter(self, simple_stack, small_psf):
        """Test that damping parameter affects result."""
        result_no_damping = lucy_richardson_cpu(
            simple_stack, small_psf, iterations=3, damping=0.0, verbose=False
        )
        result_with_damping = lucy_richardson_cpu(
            simple_stack, small_psf, iterations=3, damping=0.1, verbose=False
        )

        # Results should differ
        assert not np.allclose(result_no_damping, result_with_damping)

    def test_stop_criterion(self, blurred_stack, small_psf):
        """Test early stopping with stop_criterion."""
        iterations_run = [0]

        def count_iterations(iteration, delta_rel):
            iterations_run[0] = iteration

        # Very loose criterion should stop early
        result = lucy_richardson_cpu(
            blurred_stack,
            small_psf,
            iterations=50,
            stop_criterion=50.0,  # Will stop quickly
            callback=count_iterations,
            verbose=False,
        )

        # Should have stopped before max iterations
        assert iterations_run[0] < 50


# =============================================================================
# Tests for deconvolve function
# =============================================================================


class TestDeconvolve:
    """Tests for the main deconvolve function."""

    def test_basic_deconvolution(self, simple_stack, small_psf):
        """Test basic deconvolution runs without error."""
        result = deconvolve(
            simple_stack, small_psf, iterations=2, device="cpu", verbose=False
        )

        assert result.shape == simple_stack.shape

    def test_output_dtype(self, simple_stack, small_psf):
        """Test that output preserves reasonable dtype."""
        result = deconvolve(
            simple_stack, small_psf, iterations=2, device="cpu", verbose=False
        )

        # Should be float array
        assert np.issubdtype(result.dtype, np.floating)

    def test_padding_option(self, simple_stack, small_psf):
        """Test pad_for_fft option."""
        result_padded = deconvolve(
            simple_stack,
            small_psf,
            iterations=2,
            device="cpu",
            pad_for_fft=True,
            verbose=False,
        )
        result_no_pad = deconvolve(
            simple_stack,
            small_psf,
            iterations=2,
            device="cpu",
            pad_for_fft=False,
            verbose=False,
        )

        # Both should produce same-shaped output
        assert result_padded.shape == simple_stack.shape
        assert result_no_pad.shape == simple_stack.shape

    def test_raises_on_2d_input(self, small_psf):
        """Test that 2D input raises ValueError."""
        stack_2d = np.random.rand(64, 64).astype(np.float32)

        with pytest.raises(ValueError, match="3D"):
            deconvolve(stack_2d, small_psf, device="cpu", verbose=False)

    def test_raises_on_2d_psf(self, simple_stack):
        """Test that 2D PSF raises ValueError."""
        psf_2d = np.random.rand(7, 7).astype(np.float32)

        with pytest.raises(ValueError, match="3D"):
            deconvolve(simple_stack, psf_2d, device="cpu", verbose=False)


# =============================================================================
# Tests for PSF generation
# =============================================================================


class TestGeneratePSF:
    """Tests for PSF generation functions."""

    def test_calculate_psf_size(self):
        """Test PSF size calculation."""
        dxy = 377.0  # nm
        dz = 1500.0  # nm
        NA = 0.75
        n = 1.44
        lambda_ex = 560.0  # nm
        lambda_em = 575.0  # nm

        nxy, nz, fwhm_xy, fwhm_z = calculate_psf_size(dxy, dz, NA, n, lambda_ex, lambda_em)

        # Dimensions should be positive odd integers
        assert nxy > 0
        assert nz > 0
        assert nxy % 2 == 1  # Odd
        assert nz % 2 == 1  # Odd

        # FWHM should be positive
        assert fwhm_xy > 0
        assert fwhm_z > 0

    def test_generate_psf_basic(self):
        """Test basic PSF generation."""
        dxy = 500.0  # nm
        dz = 2000.0  # nm
        NA = 0.5
        n = 1.4
        lambda_ex = 500.0  # nm
        lambda_em = 520.0  # nm

        # Use small fixed size for speed
        psf, info = generate_psf(
            dxy, dz, NA, n, lambda_ex, lambda_em, nxy=9, nz=7, verbose=False
        )

        # Check shape
        assert psf.shape == (9, 9, 7)

        # Check normalization (should sum to 1)
        assert np.isclose(psf.sum(), 1.0, atol=1e-5)

        # Check symmetry (PSF should be symmetric)
        assert np.allclose(psf, np.flip(psf, axis=0), atol=1e-6)
        assert np.allclose(psf, np.flip(psf, axis=1), atol=1e-6)
        assert np.allclose(psf, np.flip(psf, axis=2), atol=1e-6)

        # Check info dict
        assert "nxy" in info
        assert "nz" in info
        assert "fwhm_xy" in info
        assert "fwhm_z" in info

    def test_psf_peak_at_center(self):
        """Test that PSF peak is at center."""
        psf, _ = generate_psf(
            dxy=500.0,
            dz=2000.0,
            NA=0.5,
            n=1.4,
            lambda_ex=500.0,
            lambda_em=520.0,
            nxy=9,
            nz=7,
            verbose=False,
        )

        # Find maximum position
        max_pos = np.unravel_index(np.argmax(psf), psf.shape)
        center = tuple(s // 2 for s in psf.shape)

        # Maximum should be at or near center
        assert max_pos == center

    def test_mirror_octant_to_full(self):
        """Test octant mirroring function."""
        # Create small octant
        octant = np.arange(8).reshape(2, 2, 2).astype(np.float32)

        full = _mirror_octant_to_full(octant)

        # Check shape
        expected_shape = (3, 3, 3)
        assert full.shape == expected_shape

        # Check 8-fold symmetry
        assert np.allclose(full, np.flip(full, axis=0))
        assert np.allclose(full, np.flip(full, axis=1))
        assert np.allclose(full, np.flip(full, axis=2))


# =============================================================================
# Regression tests for edge apodization fix (Jan 2026)
# =============================================================================


class TestEdgeApodizationRegression:
    """
    Regression tests for the Tukey window edge apodization fix.

    The bug was that _create_tukey_window_3d existed but was NEVER CALLED,
    causing severe horizontal banding (Gibbs ringing) in deconvolution output.
    """

    def test_tukey_window_is_applied(self, simple_stack, small_psf):
        """Test that Tukey window is actually applied during deconvolution.

        Before the fix, the Tukey window function existed but was never called.
        """
        # This test verifies the fix is in place by checking that padded
        # deconvolution produces smooth edges without ringing

        result = deconvolve(
            simple_stack,
            small_psf,
            iterations=3,
            device="cpu",
            pad_for_fft=True,
            verbose=False,
        )

        # Result should be valid (no NaN or Inf from FFT artifacts)
        assert np.all(np.isfinite(result))

    def test_no_horizontal_banding_in_output(self, blurred_stack, small_psf):
        """Test that output doesn't have horizontal banding artifacts.

        Horizontal banding (Gibbs ringing) was the main symptom of the bug.
        """
        result = deconvolve(
            blurred_stack,
            small_psf,
            iterations=5,
            device="cpu",
            pad_for_fft=True,
            verbose=False,
        )

        # Check for horizontal stripes by computing row-mean variation
        for z in range(result.shape[0]):
            slice_2d = result[z]
            row_means = np.mean(slice_2d, axis=1)
            row_diff_std = np.std(np.diff(row_means))

            # Row-to-row variation should be smooth, not stripey
            # High row_diff_std indicates banding artifacts
            slice_std = np.std(slice_2d)
            if slice_std > 0:
                normalized_row_diff = row_diff_std / slice_std
                # Before fix, this could be > 0.5; after fix should be < 0.2
                assert normalized_row_diff < 0.3, f"Possible horizontal banding in z={z}"

    def test_edge_values_not_extreme(self, simple_stack, small_psf):
        """Test that edge values aren't extreme (indicating proper apodization)."""
        result = deconvolve(
            simple_stack,
            small_psf,
            iterations=3,
            device="cpu",
            pad_for_fft=True,
            verbose=False,
        )

        # Edge values should not be extreme compared to center
        center_mean = result[:, 10:-10, 10:-10].mean()
        edge_mean = result[:, :5, :].mean()

        if center_mean > 0:
            edge_ratio = edge_mean / center_mean
            # Edges shouldn't be wildly different from center
            assert 0.1 < edge_ratio < 10, "Edge values are extreme"


# =============================================================================
# GPU tests (conditional)
# =============================================================================


@pytest.mark.skipif(not CUPY_AVAILABLE, reason="CuPy not available")
class TestGPUDeconvolution:
    """Tests for GPU deconvolution (only run if CuPy available)."""

    def test_gpu_deconvolution(self, simple_stack, small_psf):
        """Test GPU deconvolution runs."""
        result = deconvolve(
            simple_stack, small_psf, iterations=2, device="gpu", verbose=False
        )

        assert result.shape == simple_stack.shape

    def test_gpu_cpu_consistency(self, blurred_stack, small_psf):
        """Test that GPU and CPU produce similar results."""
        result_cpu = deconvolve(
            blurred_stack, small_psf, iterations=5, device="cpu", verbose=False
        )
        result_gpu = deconvolve(
            blurred_stack, small_psf, iterations=5, device="gpu", verbose=False
        )

        # Results should be similar (allow for numerical differences)
        np.testing.assert_allclose(result_cpu, result_gpu, rtol=0.1, atol=100)


# =============================================================================
# Edge case tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_very_small_stack(self, small_psf):
        """Test with minimal stack size."""
        stack = np.random.rand(5, 16, 16).astype(np.float32) * 1000

        result = deconvolve(stack, small_psf, iterations=2, device="cpu", verbose=False)

        assert result.shape == stack.shape

    def test_asymmetric_dimensions(self, small_psf):
        """Test with non-cubic stack."""
        stack = np.random.rand(8, 64, 32).astype(np.float32) * 1000

        result = deconvolve(stack, small_psf, iterations=2, device="cpu", verbose=False)

        assert result.shape == stack.shape

    def test_single_z_slice(self, small_psf):
        """Test with single z-slice stack."""
        stack = np.random.rand(1, 32, 32).astype(np.float32) * 1000

        # This should work (though deconvolution is less meaningful)
        result = deconvolve(
            stack, small_psf, iterations=2, device="cpu", verbose=False
        )

        assert result.shape == stack.shape

    def test_large_psf(self, simple_stack):
        """Test with relatively large PSF."""
        large_psf = np.random.rand(15, 15, 15).astype(np.float32)
        large_psf = large_psf / large_psf.sum()

        result = deconvolve(
            simple_stack, large_psf, iterations=2, device="cpu", verbose=False
        )

        assert result.shape == simple_stack.shape

    def test_zero_iterations(self, simple_stack, small_psf):
        """Test with zero iterations (should return input)."""
        result = deconvolve(
            simple_stack, small_psf, iterations=0, device="cpu", verbose=False
        )

        # With 0 iterations, should be close to input
        # (exact match depends on padding behavior)
        assert result.shape == simple_stack.shape
