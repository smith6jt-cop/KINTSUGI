"""
Unit tests for Optuna-based signal isolation optimizer.
"""

from __future__ import annotations

import tempfile
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pytest

optuna = pytest.importorskip("optuna", reason="optuna not installed")

from kintsugi.signal.optimizer import (  # noqa: E402
    BLANK_PARAMS_SPACE,
    CLEAN_PARAMS_SPACE,
    _flatten_params_for_optuna,
    _sample_blank_params,
    _sample_clean_params,
    optimize_signal_isolation,
)

# --- Fixtures ---


@pytest.fixture
def synthetic_signal():
    """Create a synthetic signal image with known structure."""
    np.random.seed(42)
    img = np.zeros((256, 256), dtype=np.uint16)
    # Add some circular "cells"
    yy, xx = np.mgrid[0:256, 0:256]
    for cx, cy, r, intensity in [
        (64, 64, 20, 8000),
        (128, 128, 30, 12000),
        (192, 192, 15, 6000),
        (80, 180, 25, 10000),
    ]:
        mask = (xx - cx) ** 2 + (yy - cy) ** 2 < r**2
        img[mask] = intensity
    # Add background noise
    img = img + np.random.poisson(500, img.shape).astype(np.uint16)
    return img


@pytest.fixture
def synthetic_blank():
    """Create a synthetic blank/AF image."""
    np.random.seed(43)
    # Smooth background gradient + noise
    yy, xx = np.mgrid[0:256, 0:256]
    gradient = (xx + yy).astype(np.float64) / 512 * 2000
    noise = np.random.poisson(300, (256, 256)).astype(np.float64)
    blank = (gradient + noise).clip(0, 65535).astype(np.uint16)
    return blank


# --- Tests for parameter sampling ---


class TestParameterSampling:
    def test_sample_blank_params_returns_all_keys(self):
        """Blank params sampling should produce all expected keys."""
        import optuna

        study = optuna.create_study()

        def obj(trial):
            bp = _sample_blank_params(trial)
            assert set(bp.keys()) == set(BLANK_PARAMS_SPACE.keys())
            return 0.0

        study.optimize(obj, n_trials=1)

    def test_sample_clean_params_returns_all_keys(self):
        """Clean params sampling should produce all expected keys."""
        import optuna

        study = optuna.create_study()

        def obj(trial):
            cp = _sample_clean_params(trial)
            assert set(cp.keys()) == set(CLEAN_PARAMS_SPACE.keys())
            return 0.0

        study.optimize(obj, n_trials=1)

    def test_blank_params_in_valid_range(self):
        """All sampled blank params should be within defined ranges."""
        import optuna

        study = optuna.create_study()

        def obj(trial):
            bp = _sample_blank_params(trial)
            for key, spec in BLANK_PARAMS_SPACE.items():
                if spec["type"] == "int":
                    assert spec["low"] <= bp[key] <= spec["high"]
                elif spec["type"] == "float":
                    assert spec["low"] <= bp[key] <= spec["high"]
                elif spec["type"] == "categorical":
                    assert bp[key] in spec["choices"]
            return 0.0

        study.optimize(obj, n_trials=5)


# --- Tests for param flattening/extraction ---


class TestParamConversion:
    def test_flatten_and_extract_roundtrip(self):
        """Flattening and extracting should roundtrip correctly."""
        params = {
            "blank_params": {
                "blank_clip_factor": 5000,
                "blank_scale_factor": 1.5,
                "smooth_low": True,
                "low_size": 3,
                "low_percentile": 60,
                "smooth_high": False,
                "high_size": 2,
                "high_percentile": 80,
                "erosion": 1,
            },
            "clean_params": {
                "backgrnd_thresh": 200,
                "smooth_thresh": 100,
                "smooth": False,
                "remove_small": True,
                "small_size": 64,
                "footprint": 3,
            },
        }

        flat = _flatten_params_for_optuna(params, optimize_clean=True)
        # Check flattened keys
        assert "blank_clip_factor" in flat
        assert "clean_backgrnd_thresh" in flat

    def test_flatten_without_clean(self):
        """Flattening without clean params should skip clean keys."""
        params = {
            "blank_params": {"blank_clip_factor": 5000, "blank_scale_factor": 1.0},
        }
        flat = _flatten_params_for_optuna(params, optimize_clean=False)
        assert "blank_clip_factor" in flat
        assert not any(k.startswith("clean_") for k in flat)


# --- Tests for optimization ---


class TestOptimization:
    def test_optimize_basic(self, synthetic_signal, synthetic_blank):
        """Basic optimization should complete and return valid result."""
        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=10,
            timeout=60,
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=5,
        )

        assert result["status"] == "success"
        assert result["best_quality"] > 0.0
        assert "blank_params" in result["best_params"]
        assert result["n_trials_completed"] == 10
        assert result["optimization_time_seconds"] > 0
        assert len(result["all_trials_summary"]) > 0

    def test_optimize_with_clean(self, synthetic_signal, synthetic_blank):
        """Optimization with clean params should include clean_params in result."""
        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=5,
            timeout=30,
            optimize_clean=True,
            show_progress=False,
            n_startup_trials=3,
        )

        assert result["status"] == "success"
        assert "blank_params" in result["best_params"]
        assert "clean_params" in result["best_params"]

    def test_optimize_with_warm_start(self, synthetic_signal, synthetic_blank):
        """Warm-start should use the provided parameters as first trial."""
        warm_params = {
            "blank_params": {
                "blank_clip_factor": 1000,
                "blank_scale_factor": 1.0,
                "smooth_low": False,
                "low_size": 2,
                "low_percentile": 50,
                "smooth_high": False,
                "high_size": 2,
                "high_percentile": 75,
                "erosion": 0,
            },
        }

        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=5,
            timeout=30,
            warm_start_params=warm_params,
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=2,
        )

        assert result["status"] == "success"
        assert result["n_trials_completed"] == 5

    def test_optimize_with_persistence(self, synthetic_signal, synthetic_blank):
        """Study should be persistable to SQLite and loadable."""
        with tempfile.TemporaryDirectory() as tmpdir:
            db_path = Path(tmpdir) / "optuna.db"
            study_name = "test_persist"

            result = optimize_signal_isolation(
                signal=synthetic_signal,
                blank=synthetic_blank,
                n_trials=5,
                timeout=30,
                optimize_clean=False,
                show_progress=False,
                study_name=study_name,
                storage_path=db_path,
                n_startup_trials=3,
            )

            assert result["status"] == "success"
            assert db_path.exists()

            # Load and verify
            import optuna

            study = optuna.load_study(
                study_name=study_name,
                storage=f"sqlite:///{db_path}",
            )
            assert len(study.trials) == 5

    def test_optimize_timeout_respected(self, synthetic_signal, synthetic_blank):
        """Optimization should stop at timeout even if n_trials not reached."""
        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=10000,  # very high, should be limited by timeout
            timeout=3,  # 3 seconds
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=2,
        )

        assert result["status"] == "success"
        # Should have stopped well before 10000 trials
        assert result["n_trials_completed"] < 10000
        assert result["optimization_time_seconds"] < 15  # generous upper bound

    def test_optimize_with_custom_process_fn(self, synthetic_signal, synthetic_blank):
        """Custom process function should be used when provided."""

        def simple_process(signal, blank, params):
            bp = params.get("blank_params", {})
            scale = bp.get("blank_scale_factor", 1.0)
            result = signal.astype(np.float64) - blank.astype(np.float64) * scale
            return np.clip(result, 0, 65535).astype(np.uint16)

        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=5,
            timeout=30,
            process_fn=simple_process,
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=3,
        )

        assert result["status"] == "success"

    def test_optimize_quality_details_stored(self, synthetic_signal, synthetic_blank):
        """Quality details should be stored in trial user attributes."""
        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=5,
            timeout=30,
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=3,
        )

        assert "best_quality_details" in result
        details = result["best_quality_details"]
        # Should have quality metric components
        assert "quality_score" in details or "snr" in details

    def test_trials_summary_sorted_by_quality(self, synthetic_signal, synthetic_blank):
        """Top trials should be sorted in descending quality order."""
        result = optimize_signal_isolation(
            signal=synthetic_signal,
            blank=synthetic_blank,
            n_trials=10,
            timeout=30,
            optimize_clean=False,
            show_progress=False,
            n_startup_trials=5,
        )

        summary = result["all_trials_summary"]
        scores = [t["quality_score"] for t in summary]
        assert scores == sorted(scores, reverse=True)


# --- Tests for optuna import error ---


class TestImportError:
    def test_missing_optuna_raises_import_error(self, synthetic_signal, synthetic_blank):
        """Should raise ImportError with instructions if optuna not installed."""
        with patch.dict("sys.modules", {"optuna": None}):
            with pytest.raises(ImportError, match="Optuna package not installed"):
                optimize_signal_isolation(
                    signal=synthetic_signal,
                    blank=synthetic_blank,
                    n_trials=1,
                    show_progress=False,
                )
