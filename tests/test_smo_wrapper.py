"""
Unit tests for SMO (Silver Mountain Operator) background estimation wrapper.
"""

from __future__ import annotations

from unittest.mock import MagicMock, patch

import numpy as np
import pytest


# --- Fixtures ---


@pytest.fixture
def synthetic_image_with_background():
    """Create a synthetic image with known background structure."""
    np.random.seed(42)
    # Known background: smooth gradient
    yy, xx = np.mgrid[0:128, 0:128]
    background = (xx + yy).astype(np.float64) / 256 * 2000 + 500

    # Known foreground: bright spots
    signal = np.zeros((128, 128), dtype=np.float64)
    for cx, cy, r, intensity in [
        (32, 32, 10, 8000),
        (96, 96, 15, 12000),
        (64, 64, 8, 6000),
    ]:
        mask = (xx - cx) ** 2 + (yy - cy) ** 2 < r**2
        signal[mask] = intensity

    # Combined: background + signal + noise
    combined = background + signal + np.random.poisson(100, (128, 128))
    return combined.clip(0, 65535).astype(np.uint16), background


@pytest.fixture
def uniform_background_image():
    """Image with uniform background (easy case)."""
    np.random.seed(42)
    img = np.random.poisson(500, (64, 64)).astype(np.uint16)
    # Add a bright spot
    yy, xx = np.mgrid[0:64, 0:64]
    mask = (xx - 32) ** 2 + (yy - 32) ** 2 < 100
    img[mask] = 10000
    return img


# --- Tests ---


class TestEstimateBackgroundSMO:
    @pytest.fixture(autouse=True)
    def _check_smo(self):
        """Skip tests if smo package is not installed."""
        pytest.importorskip("smo", reason="SMO package not installed")

    def test_basic_output_structure(self, synthetic_image_with_background):
        """Output should contain expected keys."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img)

        assert "corrected" in result
        assert "background_mean" in result
        assert "background_std" in result
        assert "kernel_size_used" in result
        assert result["background"] is None  # return_background=False default

    def test_corrected_dtype_preserved(self, synthetic_image_with_background):
        """Corrected image should match input dtype."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img)
        assert result["corrected"].dtype == img.dtype

    def test_corrected_non_negative(self, synthetic_image_with_background):
        """Corrected image should have no negative values."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img)
        assert np.all(result["corrected"] >= 0)

    def test_background_lower_variance(self, synthetic_image_with_background):
        """Background region of corrected image should have lower variance."""
        from kintsugi.signal.smo import estimate_background_smo

        img, true_bg = synthetic_image_with_background

        # Identify background pixels (below median)
        bg_mask = img < np.median(img)

        result = estimate_background_smo(img)
        corrected = result["corrected"].astype(np.float64)

        # Background variance should decrease after correction
        orig_bg_var = np.var(img.astype(np.float64)[bg_mask])
        corr_bg_var = np.var(corrected[bg_mask])
        assert corr_bg_var <= orig_bg_var

    def test_return_background_image(self, synthetic_image_with_background):
        """Should return background image when requested."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img, return_background=True)
        assert result["background"] is not None
        assert result["background"].shape == img.shape

    def test_background_spatially_smooth(self, synthetic_image_with_background):
        """Estimated background should be smoother than input."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img, return_background=True)
        bg = result["background"]

        # Check smoothness via gradient magnitude
        grad_img = np.gradient(img.astype(np.float64))
        grad_bg = np.gradient(bg)
        img_grad_mag = np.sqrt(grad_img[0] ** 2 + grad_img[1] ** 2).mean()
        bg_grad_mag = np.sqrt(grad_bg[0] ** 2 + grad_bg[1] ** 2).mean()
        assert bg_grad_mag < img_grad_mag

    def test_kernel_size_parameter(self, uniform_background_image):
        """Different kernel sizes should produce different results."""
        from kintsugi.signal.smo import estimate_background_smo

        img = uniform_background_image
        r1 = estimate_background_smo(img, kernel_size=3)
        r2 = estimate_background_smo(img, kernel_size=11)
        assert r1["kernel_size_used"] == 3
        assert r2["kernel_size_used"] == 11
        # Results should differ
        assert not np.array_equal(r1["corrected"], r2["corrected"])

    def test_float_input(self, synthetic_image_with_background):
        """Should handle float input images."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        img_f = img.astype(np.float64)
        result = estimate_background_smo(img_f)
        assert result["corrected"].dtype == np.float64
        assert np.all(result["corrected"] >= 0)

    def test_background_mean_positive(self, synthetic_image_with_background):
        """Estimated background mean should be positive."""
        from kintsugi.signal.smo import estimate_background_smo

        img, _ = synthetic_image_with_background
        result = estimate_background_smo(img)
        assert result["background_mean"] > 0
        assert result["background_std"] >= 0


class TestAutoSelectKernelSize:
    @pytest.fixture(autouse=True)
    def _check_smo(self):
        """Skip tests if smo package is not installed."""
        pytest.importorskip("smo", reason="SMO package not installed")

    def test_basic_output(self, uniform_background_image):
        """Should return optimal kernel size and scores dict."""
        from kintsugi.signal.smo import auto_select_kernel_size

        result = auto_select_kernel_size(uniform_background_image)
        assert "optimal_kernel_size" in result
        assert "scores" in result
        assert result["optimal_kernel_size"] in result["scores"]

    def test_custom_test_sizes(self, uniform_background_image):
        """Should test only specified kernel sizes."""
        from kintsugi.signal.smo import auto_select_kernel_size

        result = auto_select_kernel_size(uniform_background_image, test_sizes=[5, 7, 9])
        assert set(result["scores"].keys()) == {5, 7, 9}
        assert result["optimal_kernel_size"] in [5, 7, 9]

    def test_optimal_within_range(self, synthetic_image_with_background):
        """Optimal kernel should be within tested range."""
        from kintsugi.signal.smo import auto_select_kernel_size

        img, _ = synthetic_image_with_background
        test_sizes = [3, 5, 7, 9, 11]
        result = auto_select_kernel_size(img, test_sizes=test_sizes)
        assert result["optimal_kernel_size"] in test_sizes


class TestSMOImportError:
    def test_missing_smo_raises_import_error(self):
        """Should raise ImportError with instructions if smo not installed."""
        with patch.dict("sys.modules", {"smo": None}):
            # Force reimport
            import importlib

            from kintsugi.signal import smo

            importlib.reload(smo)
            img = np.zeros((64, 64), dtype=np.uint16)
            with pytest.raises(ImportError, match="SMO package not installed"):
                smo.estimate_background_smo(img)
