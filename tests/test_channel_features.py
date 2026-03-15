"""Tests for channel feature extraction."""

import numpy as np
import pytest

from kintsugi.signal.features import (
    batch_extract_features,
    extract_channel_features,
    features_to_vector,
    infer_wavelength_group,
)


class TestInferWavelengthGroup:
    def test_dapi(self):
        assert infer_wavelength_group("DAPI-01") == "405"
        assert infer_wavelength_group("Hoechst") == "405"

    def test_488(self):
        assert infer_wavelength_group("CD3-AF488") == "488"
        assert infer_wavelength_group("GFP-channel") == "488"
        assert infer_wavelength_group("FITC") == "488"

    def test_550(self):
        assert infer_wavelength_group("PE-marker") == "550"
        assert infer_wavelength_group("Cy3-labeled") == "550"

    def test_594(self):
        assert infer_wavelength_group("AF594-CD20") == "594"
        assert infer_wavelength_group("TexasRed") == "594"

    def test_650(self):
        assert infer_wavelength_group("Cy5-probe") == "650"
        assert infer_wavelength_group("AF647") == "650"
        assert infer_wavelength_group("APC-CD45") == "650"

    def test_750(self):
        assert infer_wavelength_group("Cy7-marker") == "750"
        assert infer_wavelength_group("IR800") == "750"

    def test_unknown(self):
        assert infer_wavelength_group("CD3") == "unknown"
        assert infer_wavelength_group("SomeMarker") == "unknown"
        assert infer_wavelength_group("Blank") == "unknown"

    def test_case_insensitive(self):
        assert infer_wavelength_group("DAPI") == "405"
        assert infer_wavelength_group("dapi") == "405"
        assert infer_wavelength_group("Dapi") == "405"


class TestExtractChannelFeatures:
    @pytest.fixture
    def signal_image(self):
        rng = np.random.RandomState(42)
        img = np.zeros((128, 128), dtype=np.uint16)
        # Add some structure
        img[30:70, 30:70] = 5000
        img += rng.randint(0, 200, img.shape, dtype=np.uint16)
        return img

    @pytest.fixture
    def blank_image(self):
        rng = np.random.RandomState(43)
        return rng.randint(100, 800, (128, 128), dtype=np.uint16)

    def test_correct_keys(self, signal_image):
        features = extract_channel_features(signal_image)
        expected_keys = {
            "intensity_percentiles",
            "mean",
            "std",
            "skewness",
            "kurtosis",
            "snr_estimate",
            "background_level",
            "foreground_fraction",
            "texture_entropy",
            "gradient_magnitude_mean",
            "blank_correlation",
            "blank_ratio_mean",
            "wavelength_group",
            "dynamic_range",
        }
        assert expected_keys == set(features.keys())

    def test_percentiles_length(self, signal_image):
        features = extract_channel_features(signal_image)
        assert len(features["intensity_percentiles"]) == 9

    def test_reproducibility(self, signal_image):
        f1 = extract_channel_features(signal_image)
        f2 = extract_channel_features(signal_image)
        for key in ["mean", "std", "skewness", "snr_estimate"]:
            assert abs(f1[key] - f2[key]) < f1[key] * 0.01 + 1e-10

    def test_uint8_input(self):
        img = np.random.RandomState(42).randint(0, 255, (64, 64), dtype=np.uint8)
        features = extract_channel_features(img)
        assert features["mean"] > 0

    def test_float_input(self):
        img = np.random.RandomState(42).random((64, 64)).astype(np.float32)
        features = extract_channel_features(img)
        assert isinstance(features["mean"], float)

    def test_with_blank(self, signal_image, blank_image):
        features = extract_channel_features(signal_image, blank=blank_image)
        assert features["blank_correlation"] is not None
        assert features["blank_ratio_mean"] is not None
        assert -1.0 <= features["blank_correlation"] <= 1.0

    def test_without_blank(self, signal_image):
        features = extract_channel_features(signal_image)
        assert features["blank_correlation"] is None
        assert features["blank_ratio_mean"] is None

    def test_channel_name_wavelength(self, signal_image):
        features = extract_channel_features(signal_image, channel_name="CD3-AF488")
        assert features["wavelength_group"] == "488"

    def test_no_channel_name(self, signal_image):
        features = extract_channel_features(signal_image)
        assert features["wavelength_group"] == "unknown"

    def test_foreground_fraction_range(self, signal_image):
        features = extract_channel_features(signal_image)
        assert 0.0 <= features["foreground_fraction"] <= 1.0

    def test_snr_positive_for_structured_image(self, signal_image):
        features = extract_channel_features(signal_image)
        assert features["snr_estimate"] > 0

    def test_downsample_parameter(self, signal_image):
        f1 = extract_channel_features(signal_image, downsample=1)
        f8 = extract_channel_features(signal_image, downsample=8)
        # Gradient should differ with different downsample, but other stats same
        assert abs(f1["mean"] - f8["mean"]) < 1e-10
        # Gradient may differ slightly
        assert isinstance(f8["gradient_magnitude_mean"], float)


class TestFeaturesToVector:
    def test_output_shape(self):
        img = np.random.RandomState(42).randint(0, 5000, (64, 64), dtype=np.uint16)
        features = extract_channel_features(img, channel_name="DAPI-01")
        vec = features_to_vector(features)
        assert vec.shape == (22,)
        assert vec.dtype == np.float64

    def test_none_handling(self):
        img = np.random.RandomState(42).randint(0, 5000, (64, 64), dtype=np.uint16)
        features = extract_channel_features(img)  # no blank -> None values
        vec = features_to_vector(features)
        assert not np.any(np.isnan(vec))

    def test_wavelength_one_hot(self):
        img = np.random.RandomState(42).randint(0, 5000, (64, 64), dtype=np.uint16)
        features = extract_channel_features(img, channel_name="DAPI")
        vec = features_to_vector(features)
        # Last 7 elements are one-hot wavelength
        wl_part = vec[-7:]
        assert wl_part.sum() == 1.0
        # DAPI -> 405 -> index 0
        assert wl_part[0] == 1.0

    def test_unknown_wavelength(self):
        img = np.random.RandomState(42).randint(0, 5000, (64, 64), dtype=np.uint16)
        features = extract_channel_features(img, channel_name="CD3")
        vec = features_to_vector(features)
        wl_part = vec[-7:]
        assert wl_part[-1] == 1.0  # unknown is last

    def test_numeric_values_correct(self):
        img = np.random.RandomState(42).randint(0, 5000, (64, 64), dtype=np.uint16)
        features = extract_channel_features(img)
        vec = features_to_vector(features)
        assert vec[0] == features["mean"]
        assert vec[1] == features["std"]
        assert vec[4] == features["snr_estimate"]


class TestBatchExtractFeatures:
    def test_multi_channel(self):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(0, 5000, (64, 64), dtype=np.uint16),
            "CD20": rng.randint(0, 3000, (64, 64), dtype=np.uint16),
            "DAPI": rng.randint(0, 10000, (64, 64), dtype=np.uint16),
        }
        features = batch_extract_features(marker_dict, progress=False)
        assert set(features.keys()) == {"CD3", "CD20", "DAPI"}
        for _name, f in features.items():
            assert "mean" in f
            assert "snr_estimate" in f

    def test_with_blank_map(self):
        rng = np.random.RandomState(42)
        marker_dict = {
            "CD3": rng.randint(0, 5000, (64, 64), dtype=np.uint16),
            "Blank1": rng.randint(0, 800, (64, 64), dtype=np.uint16),
        }
        blank_map = {"CD3": "Blank1"}
        features = batch_extract_features(marker_dict, blank_map=blank_map, progress=False)
        assert features["CD3"]["blank_correlation"] is not None

    def test_empty_dict(self):
        features = batch_extract_features({}, progress=False)
        assert features == {}
