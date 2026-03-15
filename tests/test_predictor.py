"""
Unit tests for parameter prediction from image features.
"""

from __future__ import annotations

import json
import sqlite3
import tempfile
from pathlib import Path

import numpy as np
import pytest

from kintsugi.signal.predictor import (
    ParameterPredictor,
    train_from_sqlite,
)

# --- Fixtures ---


def _make_feature_vector(seed=0):
    """Create a synthetic 22-element feature vector."""
    np.random.seed(seed)
    # 15 numeric features + 7 one-hot wavelength
    numeric = np.random.rand(15) * np.array(
        [65535, 65535, 65535, 100, 0.5, 5, 65535, 100, 1.0, 1.0, 65535, 1.0, 65535, 65535, 65535]
    )
    wavelength = np.zeros(7)
    wavelength[seed % 7] = 1.0
    return np.concatenate([numeric, wavelength])


def _make_params(seed=0):
    """Create synthetic parameter dict."""
    np.random.seed(seed)
    return {
        "blank_params": {
            "blank_clip_factor": int(np.random.randint(0, 30000)),
            "blank_scale_factor": float(np.round(np.random.uniform(0.5, 3.0), 2)),
            "smooth_low": bool(np.random.choice([True, False])),
            "low_size": 2,
            "low_percentile": int(np.random.choice([20, 40, 60, 80])),
            "smooth_high": bool(np.random.choice([True, False])),
            "high_size": 2,
            "high_percentile": int(np.random.choice([50, 70, 90])),
            "erosion": int(np.random.randint(0, 5)),
        },
        "clean_params": {
            "backgrnd_thresh": int(np.random.randint(0, 2000)),
            "smooth_thresh": 0,
            "smooth": False,
            "remove_small": True,
            "small_size": 64,
            "footprint": 3,
        },
    }


@pytest.fixture
def training_data():
    """Generate training data with 30 examples."""
    n = 30
    features = [_make_feature_vector(seed=i) for i in range(n)]
    params = [_make_params(seed=i) for i in range(n)]
    quality_scores = [0.5 + np.random.rand() * 0.4 for _ in range(n)]
    return features, params, quality_scores


@pytest.fixture
def small_training_data():
    """Generate insufficient training data (< MIN_TRAINING_EXAMPLES)."""
    n = 5
    features = [_make_feature_vector(seed=i) for i in range(n)]
    params = [_make_params(seed=i) for i in range(n)]
    return features, params


# --- Tests for ParameterPredictor ---


class TestParameterPredictor:
    def test_untrained_predictor_not_trained(self):
        """Untrained predictor should report is_trained=False."""
        predictor = ParameterPredictor()
        assert not predictor.is_trained
        assert predictor.n_training_examples == 0

    def test_insufficient_data_returns_status(self, small_training_data):
        """Training with too few examples should report insufficient_data."""
        features, params = small_training_data
        predictor = ParameterPredictor()
        result = predictor.train(features, params)
        assert result["status"] == "insufficient_data"
        assert result["n_examples"] == len(features)
        assert not predictor.is_trained

    def test_training_succeeds(self, training_data):
        """Training with sufficient data should succeed."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        result = predictor.train(features, params, quality_scores)

        assert result["status"] == "success"
        assert result["n_training"] > 0
        assert result["n_validation"] > 0
        assert "train_r2" in result
        assert "validation_r2" in result
        assert "per_param_r2" in result
        assert predictor.is_trained
        assert predictor.n_training_examples == 30

    def test_training_with_quality_weights(self, training_data):
        """Quality scores should be used as sample weights."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        result = predictor.train(features, params, quality_scores)
        assert result["status"] == "success"

    def test_predict_untrained_returns_error(self):
        """Predicting before training should return error."""
        predictor = ParameterPredictor()
        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv)
        assert result["status"] == "model_not_trained"
        assert result["predicted_params"] is None
        assert result["confidence"] == 0.0

    def test_predict_returns_valid_params(self, training_data):
        """Predictions should contain valid parameter structure."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        predictor.train(features, params, quality_scores)

        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv)

        assert result["status"] == "success"
        predicted = result["predicted_params"]
        assert "blank_params" in predicted
        assert "clean_params" in predicted
        assert "blank_clip_factor" in predicted["blank_params"]
        assert "blank_scale_factor" in predicted["blank_params"]
        assert "backgrnd_thresh" in predicted["clean_params"]

    def test_predict_params_in_valid_range(self, training_data):
        """Predicted parameters should be within valid ranges."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        predictor.train(features, params, quality_scores)

        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv)
        predicted = result["predicted_params"]
        bp = predicted["blank_params"]

        assert bp["blank_clip_factor"] >= 0
        assert isinstance(bp["smooth_low"], bool)
        assert isinstance(bp["smooth_high"], bool)
        assert bp["erosion"] >= 0

    def test_predict_returns_uncertainty(self, training_data):
        """Prediction should include per-parameter uncertainty."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        predictor.train(features, params, quality_scores)

        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv, return_uncertainty=True)

        assert "uncertainty" in result
        assert len(result["uncertainty"]) > 0
        for _param_name, std_val in result["uncertainty"].items():
            assert std_val >= 0.0

    def test_predict_confidence_between_0_and_1(self, training_data):
        """Confidence should be between 0 and 1."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        predictor.train(features, params, quality_scores)

        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv)

        assert 0.0 <= result["confidence"] <= 1.0

    def test_save_and_load_roundtrip(self, training_data):
        """Model should produce identical predictions after save/load."""
        features, params, quality_scores = training_data

        with tempfile.TemporaryDirectory() as tmpdir:
            model_path = Path(tmpdir) / "model.joblib"

            # Train and save
            predictor1 = ParameterPredictor(model_path=model_path)
            predictor1.train(features, params, quality_scores)

            fv = _make_feature_vector(seed=99)
            pred1 = predictor1.predict(fv, return_uncertainty=False)

            # Load and predict
            predictor2 = ParameterPredictor(model_path=model_path)
            assert predictor2.is_trained
            pred2 = predictor2.predict(fv, return_uncertainty=False)

            # Predictions should match
            bp1 = pred1["predicted_params"]["blank_params"]
            bp2 = pred2["predicted_params"]["blank_params"]
            assert bp1["blank_clip_factor"] == bp2["blank_clip_factor"]
            assert bp1["blank_scale_factor"] == bp2["blank_scale_factor"]

    def test_default_params_filled(self, training_data):
        """Non-predicted parameters should have sensible defaults."""
        features, params, quality_scores = training_data
        predictor = ParameterPredictor()
        predictor.train(features, params, quality_scores)

        fv = _make_feature_vector(seed=99)
        result = predictor.predict(fv)
        predicted = result["predicted_params"]

        # These should be filled with defaults
        assert "low_size" in predicted["blank_params"]
        assert "high_size" in predicted["blank_params"]
        assert "footprint" in predicted["clean_params"]


# --- Tests for train_from_sqlite ---


class TestTrainFromSQLite:
    def test_missing_table_returns_no_table(self):
        """Should handle missing tuning_examples table gracefully."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "params.db"
            model_path = Path(tmpdir) / "model.joblib"

            # Create empty database
            conn = sqlite3.connect(str(db_path))
            conn.close()

            result = train_from_sqlite(db_path, model_path)
            assert result["status"] == "no_table"

    def test_empty_table_returns_no_data(self):
        """Should handle empty tuning_examples table."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "params.db"
            model_path = Path(tmpdir) / "model.joblib"

            conn = sqlite3.connect(str(db_path))
            conn.execute("""
                CREATE TABLE tuning_examples (
                    features TEXT,
                    optimal_params TEXT,
                    quality_score REAL,
                    tissue_type TEXT
                )
            """)
            conn.commit()
            conn.close()

            result = train_from_sqlite(db_path, model_path)
            assert result["status"] == "no_training_data"

    def test_train_from_populated_db(self):
        """Should train successfully from populated database."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "params.db"
            model_path = Path(tmpdir) / "model.joblib"

            conn = sqlite3.connect(str(db_path))
            conn.execute("""
                CREATE TABLE tuning_examples (
                    features TEXT,
                    optimal_params TEXT,
                    quality_score REAL,
                    tissue_type TEXT
                )
            """)

            # Insert 20 examples

            for i in range(20):
                fv = _make_feature_vector(seed=i)
                # Build feature dict that features_to_vector can consume
                feat_dict = {
                    "mean": float(fv[0]),
                    "std": float(fv[1]),
                    "skewness": float(fv[2]),
                    "kurtosis": float(fv[3]),
                    "snr_estimate": float(fv[4]),
                    "background_level": float(fv[5]),
                    "foreground_fraction": float(fv[6]),
                    "texture_entropy": float(fv[7]),
                    "gradient_magnitude_mean": float(fv[8]),
                    "blank_correlation": float(fv[9]),
                    "blank_ratio_mean": float(fv[10]),
                    "dynamic_range": float(fv[11]),
                    "intensity_percentiles": [
                        float(fv[12]),
                        0,
                        0,
                        0,
                        float(fv[13]),
                        0,
                        0,
                        0,
                        float(fv[14]),
                    ],
                    "wavelength_group": ["405", "488", "550", "594", "650", "750", "unknown"][
                        i % 7
                    ],
                }
                params = _make_params(seed=i)
                conn.execute(
                    "INSERT INTO tuning_examples VALUES (?, ?, ?, ?)",
                    (json.dumps(feat_dict), json.dumps(params), 0.7, "tonsil"),
                )

            conn.commit()
            conn.close()

            result = train_from_sqlite(db_path, model_path, tissue_type="tonsil")
            assert result["status"] == "success"
            assert model_path.exists()

    def test_tissue_type_filter(self):
        """Should filter by tissue type when specified."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "params.db"
            model_path = Path(tmpdir) / "model.joblib"

            conn = sqlite3.connect(str(db_path))
            conn.execute("""
                CREATE TABLE tuning_examples (
                    features TEXT,
                    optimal_params TEXT,
                    quality_score REAL,
                    tissue_type TEXT
                )
            """)

            # Insert examples for different tissue types
            feat_dict = {
                "mean": 1000.0,
                "std": 200.0,
                "skewness": 0.5,
                "kurtosis": 0.3,
                "snr_estimate": 5.0,
                "background_level": 500.0,
                "foreground_fraction": 0.3,
                "texture_entropy": 3.0,
                "gradient_magnitude_mean": 100.0,
                "blank_correlation": 0.2,
                "blank_ratio_mean": 0.1,
                "dynamic_range": 0.06,
                "intensity_percentiles": [100.0, 200, 300, 500, 800.0, 1200, 2000, 3000, 4500.0],
                "wavelength_group": "488",
            }
            params = _make_params(seed=0)

            for _i in range(5):
                conn.execute(
                    "INSERT INTO tuning_examples VALUES (?, ?, ?, ?)",
                    (json.dumps(feat_dict), json.dumps(params), 0.7, "tonsil"),
                )
            for _i in range(3):
                conn.execute(
                    "INSERT INTO tuning_examples VALUES (?, ?, ?, ?)",
                    (json.dumps(feat_dict), json.dumps(params), 0.7, "spleen"),
                )

            conn.commit()
            conn.close()

            # Filter for spleen: only 3 examples, insufficient
            result = train_from_sqlite(db_path, model_path, tissue_type="spleen")
            assert result["status"] == "insufficient_data"
            assert result["n_examples"] == 3
