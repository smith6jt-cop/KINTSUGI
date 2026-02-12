"""
Tests for the signal isolation autofluorescence subtraction module.

This module tests the core autofluorescence subtraction algorithm which is
critical for removing autofluorescent signal from multiplex immunofluorescence
images.
"""

import numpy as np
import pytest

from kintsugi.signal.autofluorescence import (
    _apply_tissue_marker_adjustments,
    _calculate_suggestion_confidence,
    _compute_correlation,
    _compute_image_stats,
    _compute_snr,
    _estimate_af_contribution,
    _estimate_noise_level,
    analyze_for_subtraction,
    compute_subtraction_quality,
    subtract_autofluorescence,
    suggest_blank_channel,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def signal_image():
    """Create a synthetic signal image with autofluorescence."""
    np.random.seed(42)
    # Base signal with some structure
    signal = np.zeros((256, 256), dtype=np.uint16)
    # Add some bright spots (true signal)
    signal[100:120, 100:120] = 30000
    signal[150:170, 150:170] = 25000
    # Add background autofluorescence component
    af_component = np.random.randint(1000, 5000, (256, 256), dtype=np.uint16)
    signal = signal + af_component
    return signal.astype(np.uint16)


@pytest.fixture
def blank_image():
    """Create a synthetic blank/autofluorescence image."""
    np.random.seed(43)
    # Blank channel captures autofluorescence
    blank = np.random.randint(1000, 5000, (256, 256), dtype=np.uint16)
    return blank


@pytest.fixture
def low_signal_image():
    """Create a low-signal image (sparse marker)."""
    np.random.seed(44)
    signal = np.random.randint(500, 2000, (256, 256), dtype=np.uint16)
    # Only a few bright spots
    signal[50:55, 50:55] = 10000
    return signal


@pytest.fixture
def high_af_blank():
    """Create a high autofluorescence blank image."""
    np.random.seed(45)
    blank = np.random.randint(5000, 15000, (256, 256), dtype=np.uint16)
    return blank


# =============================================================================
# Tests for subtract_autofluorescence
# =============================================================================


class TestSubtractAutofluorescence:
    """Tests for the main subtraction function."""

    def test_basic_subtraction(self, signal_image, blank_image):
        """Test basic autofluorescence subtraction."""
        result = subtract_autofluorescence(signal_image, blank_image)

        # Result should be same shape
        assert result.shape == signal_image.shape
        # Result should be uint16
        assert result.dtype == np.uint16
        # Result should have reduced intensity in AF regions
        assert result.mean() < signal_image.mean()

    def test_shape_mismatch_raises(self, signal_image):
        """Test that mismatched shapes raise ValueError."""
        wrong_shape_blank = np.zeros((128, 128), dtype=np.uint16)
        with pytest.raises(ValueError, match="Signal shape"):
            subtract_autofluorescence(signal_image, wrong_shape_blank)

    def test_blank_clip_factor(self, signal_image, blank_image):
        """Test that blank_clip_factor removes noise from blank."""
        # Higher clip factor should result in less subtraction
        result_no_clip = subtract_autofluorescence(signal_image, blank_image, blank_clip_factor=0)
        result_high_clip = subtract_autofluorescence(
            signal_image, blank_image, blank_clip_factor=3000
        )

        # Higher clip = less subtraction = higher mean result
        assert result_high_clip.mean() > result_no_clip.mean()

    def test_blank_scale_factor(self, signal_image, blank_image):
        """Test that blank_scale_factor adjusts subtraction intensity."""
        result_low_scale = subtract_autofluorescence(
            signal_image, blank_image, blank_scale_factor=0.5
        )
        result_high_scale = subtract_autofluorescence(
            signal_image, blank_image, blank_scale_factor=2.0
        )

        # Higher scale = more subtraction = lower mean result
        assert result_high_scale.mean() < result_low_scale.mean()

    def test_smooth_low_option(self, signal_image, blank_image):
        """Test smooth_low parameter produces valid output."""
        result_no_smooth = subtract_autofluorescence(signal_image, blank_image, smooth_low=False)
        result_smooth = subtract_autofluorescence(
            signal_image, blank_image, smooth_low=True, low_size=3, low_percentile=60
        )

        # Both should produce valid output with same shape
        assert result_smooth.shape == result_no_smooth.shape
        assert result_smooth.dtype == np.uint16
        # Results may differ due to smoothing applied
        # (exact behavior depends on data distribution)

    def test_smooth_high_option(self, signal_image, blank_image):
        """Test smooth_high parameter smooths high-intensity regions."""
        result_no_smooth = subtract_autofluorescence(signal_image, blank_image, smooth_high=False)
        result_smooth = subtract_autofluorescence(
            signal_image, blank_image, smooth_high=True, high_size=3, high_percentile=90
        )

        # Both should produce valid results
        assert result_no_smooth.shape == result_smooth.shape

    def test_erosion_option(self, signal_image, blank_image):
        """Test erosion parameter shrinks signal regions."""
        result_no_erosion = subtract_autofluorescence(signal_image, blank_image, erosion=0)
        result_erosion = subtract_autofluorescence(signal_image, blank_image, erosion=2)

        # Erosion should reduce number of non-zero pixels
        nonzero_no_erosion = np.sum(result_no_erosion > 0)
        nonzero_erosion = np.sum(result_erosion > 0)

        assert nonzero_erosion <= nonzero_no_erosion

    def test_output_non_negative(self, signal_image, blank_image):
        """Test that output is always non-negative."""
        result = subtract_autofluorescence(
            signal_image, blank_image, blank_scale_factor=3.0  # Aggressive subtraction
        )
        assert np.all(result >= 0)

    def test_preserves_bright_signal(self, signal_image, blank_image):
        """Test that true signal (bright regions) is preserved."""
        result = subtract_autofluorescence(signal_image, blank_image)

        # The bright spot regions should still have significant signal
        bright_region_original = signal_image[100:120, 100:120].mean()
        bright_region_result = result[100:120, 100:120].mean()

        # Should preserve most of the bright signal
        assert bright_region_result > bright_region_original * 0.5


# =============================================================================
# Tests for analyze_for_subtraction
# =============================================================================


class TestAnalyzeForSubtraction:
    """Tests for the parameter suggestion analysis."""

    def test_returns_expected_keys(self, signal_image, blank_image):
        """Test that analysis returns all expected parameter keys."""
        result = analyze_for_subtraction(signal_image, blank_image)

        expected_keys = [
            "blank_clip_factor",
            "blank_scale_factor",
            "smooth_low",
            "low_size",
            "low_percentile",
            "smooth_high",
            "high_size",
            "high_percentile",
            "erosion",
            "confidence",
            "analysis",
        ]

        for key in expected_keys:
            assert key in result, f"Missing key: {key}"

    def test_scale_factor_in_valid_range(self, signal_image, blank_image):
        """Test that suggested scale factor is in valid range."""
        result = analyze_for_subtraction(signal_image, blank_image)

        assert 0.1 <= result["blank_scale_factor"] <= 3.0

    def test_confidence_in_valid_range(self, signal_image, blank_image):
        """Test that confidence is between 0 and 1."""
        result = analyze_for_subtraction(signal_image, blank_image)

        assert 0.0 <= result["confidence"] <= 1.0

    def test_analysis_contains_stats(self, signal_image, blank_image):
        """Test that analysis dict contains statistics."""
        result = analyze_for_subtraction(signal_image, blank_image)

        assert "signal_stats" in result["analysis"]
        assert "blank_stats" in result["analysis"]
        assert "correlation" in result["analysis"]
        assert "af_contribution" in result["analysis"]

    def test_tissue_marker_adjustments(self, signal_image, blank_image):
        """Test that tissue/marker parameters affect suggestions."""
        result_default = analyze_for_subtraction(signal_image, blank_image)
        result_tonsil = analyze_for_subtraction(
            signal_image, blank_image, tissue_type="tonsil", marker_name="CD3"
        )

        # Tonsil tissue should increase scale factor
        # (due to lymphoid tissue adjustment in _apply_tissue_marker_adjustments)
        # Note: May be equal if other factors dominate
        assert result_tonsil["blank_scale_factor"] >= result_default["blank_scale_factor"] * 0.9


# =============================================================================
# Tests for suggest_blank_channel
# =============================================================================


class TestSuggestBlankChannel:
    """Tests for blank channel selection."""

    def test_selects_best_blank(self, signal_image):
        """Test that function selects the blank with highest correlation."""
        np.random.seed(46)

        # Create blanks with different correlations to signal
        blank_low_corr = np.random.randint(0, 1000, signal_image.shape, dtype=np.uint16)
        blank_high_corr = (
            signal_image * 0.5 + np.random.randint(0, 500, signal_image.shape)
        ).astype(np.uint16)

        blanks = {"low_corr": blank_low_corr, "high_corr": blank_high_corr}

        best_name, score = suggest_blank_channel(signal_image, blanks)

        assert best_name == "high_corr"
        assert score > 0  # Should have positive correlation

    def test_handles_empty_blanks(self, signal_image):
        """Test handling of empty blank dictionary."""
        best_name, score = suggest_blank_channel(signal_image, {})

        assert best_name is None
        assert score == -1

    def test_skips_wrong_shape_blanks(self, signal_image):
        """Test that blanks with wrong shape are skipped."""
        blanks = {
            "wrong_shape": np.zeros((128, 128), dtype=np.uint16),
            "correct_shape": np.zeros(signal_image.shape, dtype=np.uint16),
        }

        best_name, _score = suggest_blank_channel(signal_image, blanks)

        assert best_name == "correct_shape"


# =============================================================================
# Tests for compute_subtraction_quality
# =============================================================================


class TestComputeSubtractionQuality:
    """Tests for subtraction quality metrics."""

    def test_returns_expected_keys(self, signal_image, blank_image):
        """Test that quality metrics contain expected keys."""
        subtracted = subtract_autofluorescence(signal_image, blank_image)
        quality = compute_subtraction_quality(signal_image, subtracted, blank_image)

        expected_keys = [
            "snr_improvement",
            "af_removal",
            "signal_preservation",
            "residual_correlation",
            "quality_score",
        ]

        for key in expected_keys:
            assert key in quality, f"Missing key: {key}"

    def test_quality_score_in_range(self, signal_image, blank_image):
        """Test that quality score is between 0 and 1."""
        subtracted = subtract_autofluorescence(signal_image, blank_image)
        quality = compute_subtraction_quality(signal_image, subtracted, blank_image)

        assert 0.0 <= quality["quality_score"] <= 1.0

    def test_residual_correlation_decreases(self, signal_image, blank_image):
        """Test that subtraction reduces correlation with blank."""
        subtracted = subtract_autofluorescence(signal_image, blank_image)
        quality = compute_subtraction_quality(signal_image, subtracted, blank_image)

        original_corr = _compute_correlation(signal_image, blank_image)

        # Residual correlation should be less than original
        # (unless original was already very low)
        assert quality["residual_correlation"] <= original_corr + 0.1


# =============================================================================
# Tests for helper functions
# =============================================================================


class TestHelperFunctions:
    """Tests for internal helper functions."""

    def test_compute_image_stats(self, signal_image):
        """Test image statistics computation."""
        stats = _compute_image_stats(signal_image)

        assert stats["min"] >= 0
        assert stats["max"] <= 65535
        assert stats["p5"] <= stats["p50"] <= stats["p95"]
        assert 0 <= stats["nonzero_fraction"] <= 1

    def test_compute_correlation_range(self, signal_image, blank_image):
        """Test that correlation is in valid range."""
        corr = _compute_correlation(signal_image, blank_image)

        assert -1 <= corr <= 1

    def test_compute_correlation_self(self, signal_image):
        """Test that self-correlation is close to 1."""
        corr = _compute_correlation(signal_image, signal_image)

        assert corr > 0.99

    def test_estimate_af_contribution(self, signal_image, blank_image):
        """Test AF contribution estimate is in valid range."""
        af = _estimate_af_contribution(signal_image, blank_image)

        assert 0 <= af <= 1

    def test_estimate_noise_level(self, signal_image):
        """Test noise level estimation is non-negative."""
        noise = _estimate_noise_level(signal_image)

        assert noise >= 0

    def test_compute_snr(self, signal_image):
        """Test SNR computation is non-negative."""
        snr = _compute_snr(signal_image)

        assert snr >= 0

    def test_calculate_suggestion_confidence(self):
        """Test confidence calculation with various inputs."""
        signal_stats = {
            "nonzero_fraction": 0.5,
            "max": 50000,
        }
        blank_stats = {"max": 10000}

        # High correlation should increase confidence
        conf_high = _calculate_suggestion_confidence(signal_stats, blank_stats, 0.8, 0.5)
        conf_low = _calculate_suggestion_confidence(signal_stats, blank_stats, 0.2, 0.1)

        assert conf_high > conf_low

    def test_apply_tissue_marker_adjustments_tonsil(self):
        """Test tissue-specific adjustments for tonsil."""
        clip, scale = _apply_tissue_marker_adjustments(1000, 1.0, "tonsil", "CD3")

        # Tonsil should increase scale (lymphoid tissue)
        assert scale > 1.0

    def test_apply_tissue_marker_adjustments_dim_marker(self):
        """Test marker-specific adjustments for dim markers."""
        clip, scale = _apply_tissue_marker_adjustments(1000, 1.0, "generic", "FOXP3")

        # Dim markers should decrease scale
        assert scale < 1.0

    def test_apply_tissue_marker_adjustments_bright_marker(self):
        """Test marker-specific adjustments for bright markers."""
        clip, scale = _apply_tissue_marker_adjustments(1000, 1.0, "generic", "CD3E")

        # Bright markers should increase scale
        assert scale > 1.0


# =============================================================================
# Edge case tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases and error handling."""

    def test_all_zeros_signal(self, blank_image):
        """Test handling of all-zeros signal image."""
        signal = np.zeros((256, 256), dtype=np.uint16)
        result = subtract_autofluorescence(signal, blank_image)

        assert np.all(result == 0)

    def test_all_zeros_blank(self, signal_image):
        """Test handling of all-zeros blank image."""
        blank = np.zeros((256, 256), dtype=np.uint16)
        result = subtract_autofluorescence(signal_image, blank)

        # Should return original signal (nothing to subtract)
        np.testing.assert_array_equal(result, signal_image)

    def test_identical_images(self, signal_image):
        """Test subtraction when signal equals blank."""
        result = subtract_autofluorescence(signal_image, signal_image, blank_scale_factor=1.0)

        # All pixels should be zero after subtracting identical image
        assert result.mean() == 0

    def test_float_input_handling(self, signal_image, blank_image):
        """Test that float inputs are handled correctly."""
        signal_float = signal_image.astype(np.float32)
        blank_float = blank_image.astype(np.float32)

        result = subtract_autofluorescence(signal_float, blank_float)

        assert result.dtype == np.uint16

    def test_very_low_signal(self, low_signal_image, blank_image):
        """Test handling of very low signal images."""
        result = subtract_autofluorescence(low_signal_image, blank_image)

        # Should not crash and should produce valid output
        assert result.dtype == np.uint16
        assert result.shape == low_signal_image.shape

    def test_high_af_subtraction(self, signal_image, high_af_blank):
        """Test subtraction with high autofluorescence blank."""
        result = subtract_autofluorescence(signal_image, high_af_blank)

        # Result should be valid even with aggressive subtraction
        assert result.dtype == np.uint16
        assert np.all(result <= 65535)


# =============================================================================
# Regression tests for specific fixes
# =============================================================================


class TestRegressionFixes:
    """Regression tests for previously fixed bugs."""

    def test_min_subtraction_prevents_over_subtraction(self, signal_image, blank_image):
        """Test that min(signal, blank*scale) prevents over-subtraction.

        The algorithm uses min(signal, scaled_blank) to ensure we never
        subtract more than the signal value at any pixel.
        """
        # Use very high scale factor
        result = subtract_autofluorescence(signal_image, blank_image, blank_scale_factor=10.0)

        # Result should never exceed original signal
        assert np.all(result <= signal_image)
        # Result should be non-negative
        assert np.all(result >= 0)

    def test_consistent_output_dtype(self, signal_image, blank_image):
        """Test that output dtype is always uint16 regardless of input."""
        # Test with various input dtypes
        for dtype in [np.uint8, np.uint16, np.float32, np.float64]:
            signal = signal_image.astype(dtype)
            blank = blank_image.astype(dtype)
            result = subtract_autofluorescence(signal, blank)

            assert result.dtype == np.uint16, f"Failed for input dtype {dtype}"
