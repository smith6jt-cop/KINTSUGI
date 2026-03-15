"""
Parameter Prediction from Image Features

Train a lightweight regression model (Random Forest) on accumulated
manually-tuned examples to predict signal isolation parameters from image
features. Predictions serve as warm-start values for Optuna or as direct
parameter suggestions for low-risk channels.

Design Principles:
    1. No deep learning required — scikit-learn Random Forest
    2. Trains incrementally as new manual tuning results are recorded
    3. Returns confidence/uncertainty estimates alongside predictions
    4. Falls back to tissue-type-specific defaults when insufficient data
    5. Model stored as joblib alongside the SQLite parameter learning database
"""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np

logger = logging.getLogger("kintsugi.signal.predictor")

# Minimum training examples before prediction is considered reliable
MIN_TRAINING_EXAMPLES = 15

# Minimum R² on held-out validation before predictions are offered
MIN_R2_THRESHOLD = 0.3


class ParameterPredictor:
    """
    Predict signal isolation parameters from channel image features.

    Uses a Random Forest regressor trained on accumulated (feature_vector,
    optimal_params) pairs from manual tuning sessions. Supports multi-output
    regression (one model predicts all parameters simultaneously).

    Attributes
    ----------
    model : sklearn.ensemble.RandomForestRegressor or None
    n_training_examples : int
    validation_r2 : float or None
    is_trained : bool
    """

    # Parameters to predict (subset of blank_params + clean_params)
    TARGET_PARAMS = [
        "blank_clip_factor",
        "blank_scale_factor",
        "smooth_low",  # predicted as 0/1, rounded to bool
        "low_percentile",
        "smooth_high",  # predicted as 0/1, rounded to bool
        "high_percentile",
        "erosion",
        "backgrnd_thresh",
    ]

    def __init__(self, model_path: Optional[Path] = None):
        """
        Parameters
        ----------
        model_path : Path, optional
            Path to saved model (joblib). If exists, loads it.
            If None, model must be trained before prediction.
        """
        self.model = None
        self.scaler = None
        self.n_training_examples = 0
        self.validation_r2 = None
        self._model_path = model_path
        if model_path and Path(model_path).exists():
            self._load(model_path)

    @property
    def is_trained(self) -> bool:
        return self.model is not None and self.n_training_examples >= MIN_TRAINING_EXAMPLES

    def train(
        self,
        features: list[np.ndarray],
        params: list[dict],
        quality_scores: Optional[list[float]] = None,
        validation_split: float = 0.2,
    ) -> dict:
        """
        Train the prediction model on accumulated tuning examples.

        Parameters
        ----------
        features : list[np.ndarray]
            Feature vectors from extract_channel_features -> features_to_vector.
        params : list[dict]
            Corresponding optimal parameters (blank_params + clean_params dicts).
        quality_scores : list[float], optional
            Quality scores for each example. If provided, examples are weighted
            by quality (higher quality = more influence on model).
        validation_split : float
            Fraction of data for validation (default 0.2).

        Returns
        -------
        dict with:
            status : str
            n_training : int
            n_validation : int
            train_r2 : float
            validation_r2 : float
            per_param_r2 : dict (param_name -> validation R²)
            model_path : str (if saved)
        """
        from sklearn.ensemble import RandomForestRegressor
        from sklearn.metrics import r2_score
        from sklearn.model_selection import train_test_split
        from sklearn.preprocessing import StandardScaler

        X = np.array(features)
        y = self._params_to_matrix(params)

        if len(X) < MIN_TRAINING_EXAMPLES:
            return {
                "status": "insufficient_data",
                "n_examples": len(X),
                "min_required": MIN_TRAINING_EXAMPLES,
            }

        # Sample weights from quality scores
        sample_weight = None
        if quality_scores:
            sample_weight = np.array(quality_scores)
            sample_weight = sample_weight / sample_weight.sum() * len(sample_weight)

        # Split
        X_train, X_val, y_train, y_val = train_test_split(
            X, y, test_size=validation_split, random_state=42
        )
        w_train = None
        if sample_weight is not None:
            w_train, _ = train_test_split(
                sample_weight, test_size=validation_split, random_state=42
            )

        # Scale features
        self.scaler = StandardScaler()
        X_train_s = self.scaler.fit_transform(X_train)
        X_val_s = self.scaler.transform(X_val)

        # Train multi-output Random Forest
        self.model = RandomForestRegressor(
            n_estimators=100,
            max_depth=10,
            min_samples_leaf=3,
            random_state=42,
            n_jobs=-1,
        )
        self.model.fit(X_train_s, y_train, sample_weight=w_train)

        # Evaluate
        train_pred = self.model.predict(X_train_s)
        val_pred = self.model.predict(X_val_s)
        train_r2 = r2_score(y_train, train_pred, multioutput="uniform_average")
        val_r2 = r2_score(y_val, val_pred, multioutput="uniform_average")

        # Per-parameter R²
        per_param = {}
        for i, name in enumerate(self.TARGET_PARAMS):
            per_param[name] = float(r2_score(y_val[:, i], val_pred[:, i]))

        self.n_training_examples = len(X)
        self.validation_r2 = val_r2

        # Save
        if self._model_path:
            self._save(self._model_path)

        return {
            "status": "success",
            "n_training": len(X_train),
            "n_validation": len(X_val),
            "train_r2": float(train_r2),
            "validation_r2": float(val_r2),
            "per_param_r2": per_param,
            "model_path": str(self._model_path) if self._model_path else None,
        }

    def predict(
        self,
        features: np.ndarray,
        return_uncertainty: bool = True,
    ) -> dict:
        """
        Predict parameters from a channel feature vector.

        Parameters
        ----------
        features : np.ndarray
            Feature vector from features_to_vector().
        return_uncertainty : bool
            If True, compute prediction uncertainty from Random Forest
            tree variance (requires individual tree predictions).

        Returns
        -------
        dict with:
            status : str
            predicted_params : dict (blank_params + clean_params)
            confidence : float 0-1 (based on model R² and tree agreement)
            uncertainty : dict (per-parameter std from tree predictions)
            warning : str or None
        """
        if not self.is_trained:
            return {
                "status": "model_not_trained",
                "predicted_params": None,
                "confidence": 0.0,
                "warning": (
                    f"Need {MIN_TRAINING_EXAMPLES} examples, "
                    f"have {self.n_training_examples}"
                ),
            }

        X = self.scaler.transform(features.reshape(1, -1))
        y_pred = self.model.predict(X)[0]

        # Convert prediction back to parameter dict
        predicted = self._vector_to_params(y_pred)

        # Uncertainty from tree variance
        uncertainty = {}
        confidence = min(1.0, max(0.0, self.validation_r2 or 0.0))

        if return_uncertainty:
            tree_preds = np.array([t.predict(X)[0] for t in self.model.estimators_])
            tree_std = tree_preds.std(axis=0)
            for i, name in enumerate(self.TARGET_PARAMS):
                uncertainty[name] = float(tree_std[i])

            # Reduce confidence if trees disagree strongly
            mean_cv = np.mean(tree_std / (np.abs(y_pred) + 1e-10))
            confidence *= max(0.2, 1.0 - mean_cv)

        warning = None
        if confidence < 0.5:
            warning = "Low prediction confidence; manual verification strongly recommended."

        return {
            "status": "success",
            "predicted_params": predicted,
            "confidence": float(confidence),
            "uncertainty": uncertainty,
            "warning": warning,
        }

    def _params_to_matrix(self, params_list: list[dict]) -> np.ndarray:
        """Convert list of param dicts to numpy matrix for training."""
        rows = []
        for p in params_list:
            bp = p.get("blank_params", p)
            cp = p.get("clean_params", {})
            merged = {**bp, **cp}
            row = []
            for name in self.TARGET_PARAMS:
                val = merged.get(name, 0)
                if isinstance(val, bool):
                    val = float(val)
                row.append(float(val))
            rows.append(row)
        return np.array(rows)

    def _vector_to_params(self, y: np.ndarray) -> dict:
        """Convert prediction vector back to structured param dict."""
        blank_keys = {
            "blank_clip_factor",
            "blank_scale_factor",
            "smooth_low",
            "low_percentile",
            "smooth_high",
            "high_percentile",
            "erosion",
        }
        bp = {}
        cp = {}

        for i, name in enumerate(self.TARGET_PARAMS):
            val = y[i]
            # Type enforcement
            if name in ("smooth_low", "smooth_high"):
                val = bool(round(val))
            elif name in (
                "blank_clip_factor",
                "low_percentile",
                "high_percentile",
                "erosion",
                "backgrnd_thresh",
            ):
                val = max(0, int(round(val)))
            else:
                val = float(round(val, 2))

            if name in blank_keys:
                bp[name] = val
            else:
                cp[name] = val

        # Fill defaults for parameters not predicted
        bp.setdefault("low_size", 2)
        bp.setdefault("high_size", 2)
        cp.setdefault("smooth_thresh", 0)
        cp.setdefault("smooth", False)
        cp.setdefault("remove_small", True)
        cp.setdefault("small_size", 64)
        cp.setdefault("footprint", 3)

        return {"blank_params": bp, "clean_params": cp}

    def _save(self, path: Path):
        """Save model, scaler, and metadata to disk."""
        import joblib

        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        joblib.dump(
            {
                "model": self.model,
                "scaler": self.scaler,
                "n_training_examples": self.n_training_examples,
                "validation_r2": self.validation_r2,
                "target_params": self.TARGET_PARAMS,
            },
            path,
        )

    def _load(self, path: Path):
        """Load saved model from disk."""
        import joblib

        data = joblib.load(path)
        self.model = data["model"]
        self.scaler = data["scaler"]
        self.n_training_examples = data["n_training_examples"]
        self.validation_r2 = data["validation_r2"]


def train_from_sqlite(
    db_path: Path,
    model_path: Path,
    tissue_type: Optional[str] = None,
) -> dict:
    """
    Train a ParameterPredictor from the KINTSUGI parameter learning SQLite database.

    Queries the tuning_examples table for feature vectors and verified optimal
    parameters.

    Parameters
    ----------
    db_path : Path
        Path to KINTSUGI parameter learning SQLite database.
    model_path : Path
        Path to save trained model.
    tissue_type : str, optional
        Filter training data by tissue type (e.g., 'tonsil').

    Returns
    -------
    dict : training result from ParameterPredictor.train()
    """
    import sqlite3

    from kintsugi.signal.features import features_to_vector

    conn = sqlite3.connect(str(db_path))
    cursor = conn.cursor()

    # Check if table exists
    cursor.execute(
        "SELECT name FROM sqlite_master WHERE type='table' AND name='tuning_examples'"
    )
    if not cursor.fetchone():
        conn.close()
        return {"status": "no_table", "n_examples": 0}

    query = "SELECT features, optimal_params, quality_score FROM tuning_examples"
    params_q = []
    if tissue_type:
        query += " WHERE tissue_type = ?"
        params_q.append(tissue_type)

    cursor.execute(query, params_q)
    rows = cursor.fetchall()
    conn.close()

    if not rows:
        return {"status": "no_training_data", "n_examples": 0}

    features = []
    params = []
    quality_scores = []

    for feat_json, param_json, score in rows:
        feat_dict = json.loads(feat_json)
        features.append(features_to_vector(feat_dict))
        params.append(json.loads(param_json))
        quality_scores.append(score)

    predictor = ParameterPredictor(model_path=model_path)
    return predictor.train(features, params, quality_scores)
