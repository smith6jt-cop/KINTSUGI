"""Tests for signal isolation quality scoring."""

import numpy as np

from kintsugi.qc.signal_quality import compute_signal_isolation_quality


def _make_signal_image(size=256, spot_intensity=10000, bg_level=500, seed=42):
    """Create a synthetic signal image with a gaussian spot."""
    rng = np.random.RandomState(seed)
    img = np.full((size, size), bg_level, dtype=np.float64)
    # Add gaussian spot in center
    y, x = np.mgrid[:size, :size]
    center = size // 2
    sigma = size // 8
    spot = spot_intensity * np.exp(-((x - center) ** 2 + (y - center) ** 2) / (2 * sigma**2))
    img += spot
    # Add noise
    img += rng.normal(0, 50, img.shape)
    img = np.clip(img, 0, 65535)
    return img.astype(np.uint16)


def _make_blank_image(size=256, level=800, seed=43):
    """Create a synthetic blank/autofluorescence image."""
    rng = np.random.RandomState(seed)
    img = np.full((size, size), level, dtype=np.float64)
    img += rng.normal(0, 30, img.shape)
    img = np.clip(img, 0, 65535)
    return img.astype(np.uint16)


class TestComputeSignalIsolationQuality:
    def test_good_result_has_positive_score(self):
        """Well-processed image should have positive quality score with good structure."""
        original = _make_signal_image(spot_intensity=30000, bg_level=2000)
        blank = _make_blank_image(level=2000)
        # Good processing: subtract blank, preserve signal structure
        processed = np.clip(original.astype(np.float64) - blank * 0.8, 0, 65535).astype(np.uint16)
        result = compute_signal_isolation_quality(original, processed, blank=blank)

        assert result["quality_score"] > 0.0
        assert result["structure_preservation"] > 0.5
        assert result["artifact_flag"] is False
        assert "components" in result
        assert "snr" in result

    def test_over_subtracted_scores_low(self):
        """Over-subtracted image (mostly zeros) should score low and flag artifact."""
        original = _make_signal_image(spot_intensity=5000, bg_level=500)
        blank = _make_blank_image(level=800)
        # Over-subtraction: result is mostly zeros
        processed = np.zeros_like(original)
        result = compute_signal_isolation_quality(original, processed, blank=blank)

        assert result["quality_score"] < 0.4
        assert result["artifact_flag"] is True

    def test_identical_images_high_structure(self):
        """Identical original and processed should have high structure preservation."""
        original = _make_signal_image()
        result = compute_signal_isolation_quality(original, original.copy())

        assert result["structure_preservation"] > 0.9

    def test_custom_weights(self):
        """Custom weights should be used in composite calculation."""
        original = _make_signal_image()
        processed = original.copy()
        weights = {
            "snr": 0.5,
            "dynamic_range": 0.1,
            "bg_uniformity": 0.1,
            "structure": 0.2,
            "residual_af": 0.1,
        }
        result = compute_signal_isolation_quality(original, processed, weights=weights)
        assert "quality_score" in result
        assert isinstance(result["quality_score"], float)

    def test_with_mask(self):
        """Providing an explicit mask should work."""
        original = _make_signal_image()
        processed = original.copy()
        mask = original > np.median(original)
        result = compute_signal_isolation_quality(original, processed, mask=mask)
        assert "quality_score" in result

    def test_artifact_flag_negative_values(self):
        """Negative values in processed image should trigger artifact flag."""
        original = _make_signal_image()
        processed = original.astype(np.float64) - 1000
        # Keep as float to allow negatives
        result = compute_signal_isolation_quality(original, processed.astype(np.float64))
        assert result["artifact_flag"] is True

    def test_uniform_image(self):
        """Uniform image should not crash."""
        uniform = np.full((64, 64), 1000, dtype=np.uint16)
        result = compute_signal_isolation_quality(uniform, uniform.copy())
        assert "quality_score" in result

    def test_result_keys(self):
        """All expected keys should be present in result."""
        original = _make_signal_image(size=64)
        processed = original.copy()
        result = compute_signal_isolation_quality(original, processed)

        expected_keys = {
            "quality_score",
            "snr",
            "dynamic_range_utilization",
            "background_uniformity",
            "structure_preservation",
            "residual_af_score",
            "artifact_flag",
            "components",
            "pass_threshold",
            "passed",
        }
        assert expected_keys.issubset(set(result.keys()))

    def test_quality_score_bounded(self):
        """Quality score should be between 0 and 1."""
        original = _make_signal_image(size=64)
        blank = _make_blank_image(size=64)
        processed = np.clip(original.astype(np.float64) - blank * 0.5, 0, 65535).astype(np.uint16)
        result = compute_signal_isolation_quality(original, processed, blank=blank)
        assert 0.0 <= result["quality_score"] <= 1.0

    def test_without_blank(self):
        """Should work without blank channel (residual_af neutral)."""
        original = _make_signal_image(size=64)
        processed = original.copy()
        result = compute_signal_isolation_quality(original, processed)
        assert result["residual_af_score"] == 0.5  # neutral

    def test_saturated_image_flags_artifact(self):
        """Image with >1% saturated pixels should flag artifact."""
        original = _make_signal_image(size=64)
        processed = np.full_like(original, 65535)
        # Set a few pixels low so it's not 100% saturated
        processed[:5, :5] = 0
        result = compute_signal_isolation_quality(original, processed)
        assert result["artifact_flag"] is True
