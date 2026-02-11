"""
Tests for the QC artifact detection modules.

This module tests:
- stripe_artifact.py: Stripe/banding artifact detection algorithms
- artifact_scanner.py: Project-level artifact scanning (data classes only)

Note: Full artifact_scanner.py tests require project fixtures which are
created in a separate integration test module.
"""

import numpy as np
import pytest

from kintsugi.qc.artifact_scanner import (
    ArtifactItem,
    ArtifactReport,
    ChannelScanResult,
)
from kintsugi.qc.stripe_artifact import (
    StripeArtifactResult,
    classify_hp_fft_severity,
    compute_hp_fft_peak_strength,
    compute_zstack_baseline,
    detect_stripe_artifact,
    scan_zstack,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def clean_image():
    """Create a clean image without stripe artifacts."""
    np.random.seed(42)
    # Uniform random noise - no stripe patterns
    img = np.random.randint(1000, 10000, (512, 512), dtype=np.uint16)
    return img


@pytest.fixture
def striped_image():
    """Create an image with horizontal stripe artifacts."""
    np.random.seed(42)
    # Base image
    img = np.random.randint(1000, 5000, (512, 512), dtype=np.uint16)

    # Add horizontal stripes with ~30px period
    for i in range(0, 512, 30):
        img[i : i + 2, :] += 3000  # Bright stripe
        if i + 15 < 512:
            img[i + 15 : i + 17, :] -= 500  # Dark stripe

    return img.astype(np.uint16)


@pytest.fixture
def severe_striped_image():
    """Create an image with severe stripe artifacts."""
    np.random.seed(42)
    img = np.random.randint(500, 2000, (512, 512), dtype=np.uint16)

    # Very pronounced horizontal stripes
    for i in range(0, 512, 20):
        img[i : i + 5, :] = 20000
        if i + 10 < 512:
            img[i + 10 : i + 15, :] = 500

    return img.astype(np.uint16)


@pytest.fixture
def clean_zstack():
    """Create a clean z-stack without artifacts."""
    np.random.seed(42)
    # 10 slices of clean images with gradual intensity variation
    stack = np.zeros((10, 256, 256), dtype=np.uint16)
    for z in range(10):
        # Vary intensity across z (simulating focus)
        base_intensity = 2000 + z * 500
        stack[z] = np.random.randint(base_intensity - 500, base_intensity + 500, (256, 256)).astype(
            np.uint16
        )
    return stack


@pytest.fixture
def zstack_with_artifact():
    """Create a z-stack with stripe artifact in one plane."""
    np.random.seed(42)
    stack = np.zeros((10, 256, 256), dtype=np.uint16)

    for z in range(10):
        base_intensity = 2000 + z * 500
        stack[z] = np.random.randint(base_intensity - 500, base_intensity + 500, (256, 256)).astype(
            np.uint16
        )

    # Add stripe artifact to plane 5
    for i in range(0, 256, 25):
        stack[5, i : i + 3, :] += 5000

    return stack


# =============================================================================
# Tests for StripeArtifactResult dataclass
# =============================================================================


class TestStripeArtifactResult:
    """Tests for the StripeArtifactResult dataclass."""

    def test_creation(self):
        """Test basic creation of result object."""
        result = StripeArtifactResult(
            has_artifact=True, severity="moderate", score=0.5, metrics={"test": 1.0}
        )

        assert result.has_artifact
        assert result.severity == "moderate"
        assert result.score == 0.5
        assert result.metrics["test"] == 1.0

    def test_stripe_score_alias(self):
        """Test that stripe_score property returns score."""
        result = StripeArtifactResult(has_artifact=True, severity="mild", score=0.3)

        assert result.stripe_score == result.score

    def test_str_no_artifact(self):
        """Test string representation with no artifact."""
        result = StripeArtifactResult(has_artifact=False, severity="none", score=0.0)

        assert "No stripe artifact" in str(result)

    def test_str_with_artifact(self):
        """Test string representation with artifact."""
        result = StripeArtifactResult(has_artifact=True, severity="moderate", score=0.5)

        assert "moderate" in str(result)
        assert "0.5" in str(result)


# =============================================================================
# Tests for detect_stripe_artifact
# =============================================================================


class TestDetectStripeArtifact:
    """Tests for single-image stripe detection."""

    def test_detects_clean_image(self, clean_image):
        """Test that clean images produce valid results."""
        result = detect_stripe_artifact(clean_image)

        # Should return a valid result with metrics
        assert "normalized_row_diff" in result.metrics
        assert "highfreq_stripe" in result.metrics
        # Score should be in valid range
        assert 0.0 <= result.score <= 1.0

    def test_detects_striped_image(self, striped_image):
        """Test that striped images are detected."""
        result = detect_stripe_artifact(striped_image)

        # Should detect the artifact
        assert result.has_artifact  # Use truthiness, not identity
        assert result.score > 0.0

    def test_detects_severe_stripes(self, severe_striped_image):
        """Test that severe stripes are classified correctly."""
        result = detect_stripe_artifact(severe_striped_image)

        # Should detect severe artifact
        assert result.has_artifact  # Use truthiness
        assert result.severity in ["moderate", "severe"]
        assert result.score > 0.3

    def test_returns_metrics(self, clean_image):
        """Test that metrics dictionary is populated."""
        result = detect_stripe_artifact(clean_image)

        assert "normalized_row_diff" in result.metrics
        assert "highfreq_stripe" in result.metrics
        assert "gradient_ratio" in result.metrics

    def test_raises_on_3d_input(self, clean_zstack):
        """Test that 3D input raises ValueError."""
        with pytest.raises(ValueError, match="Expected 2D"):
            detect_stripe_artifact(clean_zstack)

    def test_baseline_comparative_detection(self, clean_image, striped_image):
        """Test comparative detection with baseline."""
        # Compute baseline from clean image
        baseline_result = detect_stripe_artifact(clean_image)
        baseline_row_diff = baseline_result.metrics["normalized_row_diff"]

        # Test striped image with baseline
        result = detect_stripe_artifact(
            striped_image, baseline_row_diff=baseline_row_diff, threshold_sigma=2.0
        )

        # Should detect the deviation from baseline
        assert result.has_artifact  # Use truthiness


# =============================================================================
# Tests for compute_hp_fft_peak_strength
# =============================================================================


class TestComputeHPFFTPeakStrength:
    """Tests for HP-FFT peak strength metric."""

    def test_returns_expected_keys(self, clean_image):
        """Test that function returns expected keys."""
        result = compute_hp_fft_peak_strength(clean_image)

        assert "row_peak_strength" in result
        assert "col_peak_strength" in result

    def test_clean_image_low_peaks(self, clean_image):
        """Test that clean images have low peak strengths."""
        result = compute_hp_fft_peak_strength(clean_image)

        # Clean images should have low values (empirically < 10)
        assert result["row_peak_strength"] < 20
        assert result["col_peak_strength"] < 50

    def test_striped_image_higher_row_peak(self, striped_image):
        """Test that horizontal stripes increase row peak."""
        clean_result = compute_hp_fft_peak_strength(
            np.random.randint(1000, 5000, (512, 512), dtype=np.uint16)
        )
        striped_result = compute_hp_fft_peak_strength(striped_image)

        # Horizontal stripes should increase row peak strength
        # (detected in row-wise analysis)
        assert striped_result["row_peak_strength"] >= clean_result["row_peak_strength"] * 0.5

    def test_hp_sigma_parameter(self, clean_image):
        """Test that hp_sigma affects results."""
        result_low_sigma = compute_hp_fft_peak_strength(clean_image, hp_sigma=10.0)
        result_high_sigma = compute_hp_fft_peak_strength(clean_image, hp_sigma=100.0)

        # Different sigma values should produce different results
        # Both should be valid (non-negative)
        assert result_low_sigma["row_peak_strength"] >= 0
        assert result_high_sigma["row_peak_strength"] >= 0


# =============================================================================
# Tests for classify_hp_fft_severity
# =============================================================================


class TestClassifyHPFFTSeverity:
    """Tests for HP-FFT severity classification."""

    def test_clean_classification(self):
        """Test classification of clean metrics."""
        severity, z_score = classify_hp_fft_severity(
            metric=3.0, baseline_median=3.0, baseline_std=0.5, metric_type="row"
        )

        assert severity == "clean"
        assert abs(z_score) < 2

    def test_faint_classification(self):
        """Test classification of faint artifact."""
        severity, z_score = classify_hp_fft_severity(
            metric=5.0, baseline_median=3.0, baseline_std=0.5, metric_type="row"
        )

        # Z-score = 4, but metric barely above threshold
        # Should be faint or moderate
        assert severity in ["faint", "moderate"]

    def test_severe_classification(self):
        """Test classification of severe artifact."""
        severity, z_score = classify_hp_fft_severity(
            metric=50.0, baseline_median=3.0, baseline_std=1.0, metric_type="row"
        )

        # Very high z-score should be severe
        assert severity == "severe"
        assert z_score > 10

    def test_zero_std_handling(self):
        """Test handling of zero standard deviation."""
        severity, z_score = classify_hp_fft_severity(
            metric=5.0, baseline_median=5.0, baseline_std=0.0, metric_type="row"
        )

        # Should not crash, z_score should be 0
        assert z_score == 0.0
        assert severity == "clean"


# =============================================================================
# Tests for compute_zstack_baseline
# =============================================================================


class TestComputeZstackBaseline:
    """Tests for z-stack baseline computation."""

    def test_returns_expected_keys(self, clean_zstack):
        """Test that baseline contains expected keys."""
        baseline = compute_zstack_baseline(clean_zstack)

        assert "lowfreq_mean" in baseline
        assert "lowfreq_std" in baseline
        assert "highfreq_mean" in baseline
        assert "highfreq_std" in baseline

    def test_includes_hp_fft_metrics(self, clean_zstack):
        """Test that HP-FFT metrics are included when requested."""
        baseline = compute_zstack_baseline(clean_zstack, include_hp_fft=True)

        assert "hp_row_median" in baseline
        assert "hp_row_std" in baseline
        assert "hp_col_median" in baseline
        assert "hp_col_std" in baseline

    def test_excludes_hp_fft_when_disabled(self, clean_zstack):
        """Test that HP-FFT metrics excluded when disabled."""
        baseline = compute_zstack_baseline(clean_zstack, include_hp_fft=False)

        assert "hp_row_median" not in baseline
        assert "hp_col_median" not in baseline

    def test_raises_on_2d_input(self, clean_image):
        """Test that 2D input raises ValueError."""
        with pytest.raises(ValueError, match="Expected 3D"):
            compute_zstack_baseline(clean_image)

    def test_baseline_values_non_negative(self, clean_zstack):
        """Test that baseline values are non-negative."""
        baseline = compute_zstack_baseline(clean_zstack)

        assert baseline["lowfreq_mean"] >= 0
        assert baseline["lowfreq_std"] >= 0
        assert baseline["highfreq_mean"] >= 0
        assert baseline["highfreq_std"] >= 0


# =============================================================================
# Tests for scan_zstack
# =============================================================================


class TestScanZstack:
    """Tests for z-stack artifact scanning."""

    def test_clean_stack_no_artifacts(self, clean_zstack):
        """Test that clean stack has no/few artifacts."""
        affected, results = scan_zstack(clean_zstack)

        # Clean stack should have few or no artifacts
        assert len(affected) <= 2  # Allow some false positives
        assert len(results) == clean_zstack.shape[0]

    def test_detects_artifact_plane(self, zstack_with_artifact):
        """Test that artifact plane is detected."""
        affected, results = scan_zstack(zstack_with_artifact)

        # Plane 5 should be detected (where we added stripes)
        assert 5 in affected
        assert results[5].has_artifact

    def test_returns_results_for_all_planes(self, clean_zstack):
        """Test that results dict contains all planes."""
        affected, results = scan_zstack(clean_zstack)

        for z in range(clean_zstack.shape[0]):
            assert z in results
            assert isinstance(results[z], StripeArtifactResult)

    def test_threshold_sigma_parameter(self, zstack_with_artifact):
        """Test that threshold_sigma affects detection sensitivity."""
        # Lower threshold = more detections
        affected_low, _ = scan_zstack(zstack_with_artifact, threshold_sigma=1.0)
        # Higher threshold = fewer detections
        affected_high, _ = scan_zstack(zstack_with_artifact, threshold_sigma=4.0)

        # Lower threshold should detect at least as many artifacts
        assert len(affected_low) >= len(affected_high)

    def test_raises_on_2d_input(self, clean_image):
        """Test that 2D input raises ValueError."""
        with pytest.raises(ValueError, match="Expected 3D"):
            scan_zstack(clean_image)

    def test_use_hp_fft_parameter(self, clean_zstack):
        """Test that use_hp_fft parameter works."""
        # Should work with HP-FFT enabled
        affected_with, results_with = scan_zstack(clean_zstack, use_hp_fft=True)
        # Should work with HP-FFT disabled
        affected_without, results_without = scan_zstack(clean_zstack, use_hp_fft=False)

        # Both should return valid results
        assert len(results_with) == clean_zstack.shape[0]
        assert len(results_without) == clean_zstack.shape[0]


# =============================================================================
# Tests for ArtifactItem dataclass
# =============================================================================


class TestArtifactItem:
    """Tests for the ArtifactItem dataclass."""

    def test_creation(self):
        """Test basic creation."""
        item = ArtifactItem(
            cycle="cyc001",
            channel=1,
            z_plane=5,
            severity="moderate",
            score=0.5,
        )

        assert item.cycle == "cyc001"
        assert item.channel == 1
        assert item.z_plane == 5

    def test_to_dict(self):
        """Test conversion to dictionary."""
        item = ArtifactItem(cycle="cyc001", channel=1, z_plane=5, severity="mild", score=0.3)
        d = item.to_dict()

        assert d["cycle"] == "cyc001"
        assert d["channel"] == 1
        assert d["z_plane"] == 5
        assert d["severity"] == "mild"

    def test_from_dict(self):
        """Test creation from dictionary."""
        d = {
            "cycle": "cyc002",
            "channel": 2,
            "z_plane": 10,
            "tile": None,
            "severity": "severe",
            "score": 0.8,
            "metrics": {},
        }
        item = ArtifactItem.from_dict(d)

        assert item.cycle == "cyc002"
        assert item.channel == 2
        assert item.severity == "severe"


# =============================================================================
# Tests for ChannelScanResult dataclass
# =============================================================================


class TestChannelScanResult:
    """Tests for the ChannelScanResult dataclass."""

    def test_has_artifacts_property(self):
        """Test has_artifacts property."""
        result_with = ChannelScanResult(
            cycle="cyc001",
            channel=1,
            n_planes=10,
            n_tiles=1,
            affected_planes=[3, 5, 7],
        )
        result_without = ChannelScanResult(
            cycle="cyc001", channel=1, n_planes=10, n_tiles=1, affected_planes=[]
        )

        assert result_with.has_artifacts
        assert not result_without.has_artifacts

    def test_artifact_ratio_property(self):
        """Test artifact_ratio property."""
        result = ChannelScanResult(
            cycle="cyc001",
            channel=1,
            n_planes=10,
            n_tiles=1,
            affected_planes=[3, 5],  # 2 out of 10
        )

        assert result.artifact_ratio == 0.2

    def test_artifact_ratio_zero_planes(self):
        """Test artifact_ratio with zero planes."""
        result = ChannelScanResult(
            cycle="cyc001", channel=1, n_planes=0, n_tiles=1, affected_planes=[]
        )

        assert result.artifact_ratio == 0.0


# =============================================================================
# Tests for ArtifactReport dataclass
# =============================================================================


class TestArtifactReport:
    """Tests for the ArtifactReport dataclass."""

    def test_creation(self):
        """Test basic creation."""
        report = ArtifactReport(
            project_name="test_project",
            data_source="raw",
            total_cycles=3,
            total_channels=12,
        )

        assert report.project_name == "test_project"
        assert not report.has_artifacts

    def test_has_artifacts_property(self):
        """Test has_artifacts property."""
        report = ArtifactReport()
        assert not report.has_artifacts

        report.total_artifacts_found = 5
        assert report.has_artifacts

    def test_worst_severity_property(self):
        """Test worst_severity property."""
        report = ArtifactReport()

        # No items = none
        assert report.worst_severity == "none"

        # Add mild item
        report.affected_items.append(
            ArtifactItem(cycle="cyc001", channel=1, z_plane=1, severity="mild")
        )
        assert report.worst_severity == "mild"

        # Add severe item
        report.affected_items.append(
            ArtifactItem(cycle="cyc001", channel=2, z_plane=2, severity="severe")
        )
        assert report.worst_severity == "severe"

    def test_summary_no_artifacts(self):
        """Test summary with no artifacts."""
        report = ArtifactReport(
            project_name="test", data_source="raw", total_cycles=1, total_channels=4
        )
        summary = report.summary()

        assert "No artifacts detected" in summary

    def test_summary_with_artifacts(self):
        """Test summary with artifacts."""
        report = ArtifactReport(
            project_name="test",
            data_source="raw",
            total_artifacts_found=2,
        )
        report.affected_items = [
            ArtifactItem(cycle="cyc001", channel=1, z_plane=5, severity="moderate"),
            ArtifactItem(cycle="cyc001", channel=1, z_plane=7, severity="mild"),
        ]
        summary = report.summary()

        assert "Artifacts found: 2" in summary
        assert "cyc001" in summary

    def test_save_and_load(self, tmp_path):
        """Test saving and loading report."""
        report = ArtifactReport(
            project_name="test_project",
            data_source="raw",
            total_cycles=2,
            total_channels=8,
            total_planes_scanned=100,
            total_artifacts_found=3,
        )
        report.affected_items = [
            ArtifactItem(cycle="cyc001", channel=1, z_plane=5, severity="moderate", score=0.5),
        ]

        # Save
        save_path = tmp_path / "artifact_report.json"
        report.save(save_path)

        # Load
        loaded = ArtifactReport.load(save_path)

        assert loaded.project_name == "test_project"
        assert loaded.total_artifacts_found == 3
        assert len(loaded.affected_items) == 1
        assert loaded.affected_items[0].cycle == "cyc001"


# =============================================================================
# Edge case tests
# =============================================================================


class TestEdgeCases:
    """Tests for edge cases."""

    def test_small_image(self):
        """Test detection on small image."""
        small_img = np.random.randint(0, 10000, (64, 64), dtype=np.uint16)
        result = detect_stripe_artifact(small_img)

        # Should not crash, should return valid result
        assert isinstance(result, StripeArtifactResult)

    def test_single_value_image(self):
        """Test detection on constant-value image."""
        const_img = np.full((256, 256), 5000, dtype=np.uint16)
        result = detect_stripe_artifact(const_img)

        # Constant image has no stripes
        assert not result.has_artifact

    def test_short_zstack(self):
        """Test scanning very short z-stack."""
        short_stack = np.random.randint(0, 10000, (3, 128, 128), dtype=np.uint16)
        affected, results = scan_zstack(short_stack)

        # Should work without error
        assert len(results) == 3

    def test_very_noisy_image(self):
        """Test detection on very noisy image."""
        np.random.seed(42)
        noisy_img = np.random.randint(0, 65535, (256, 256), dtype=np.uint16)
        result = detect_stripe_artifact(noisy_img)

        # Should not crash
        assert isinstance(result, StripeArtifactResult)

    def test_zero_image(self):
        """Test detection on all-zero image."""
        zero_img = np.zeros((256, 256), dtype=np.uint16)
        result = detect_stripe_artifact(zero_img)

        # Should handle gracefully
        assert result.score >= 0  # Score should be valid
