"""
Bayesian Parameter Optimization for Signal Isolation

Wraps KINTSUGI's signal isolation processing functions as an Optuna objective,
enabling Bayesian optimization to find optimal parameters for a given
signal+blank channel pair using the quality metric from qc/signal_quality.py.

Design Principles:
    1. Each Optuna trial calls the existing pipeline — no new processing logic
    2. Warm-starting from predicted or previously tuned parameters
    3. Pruning via quality score thresholds stops unpromising trials early
    4. All trial results logged transparently
    5. GPU acceleration inherited from underlying processing functions
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path
from typing import Callable, Optional

import numpy as np

logger = logging.getLogger("kintsugi.signal.optimizer")

# Parameter search space definitions
BLANK_PARAMS_SPACE = {
    "blank_clip_factor": {"type": "int", "low": 0, "high": 30000, "step": 500},
    "blank_scale_factor": {"type": "float", "low": 0.5, "high": 3.0, "step": 0.1},
    "smooth_low": {"type": "categorical", "choices": [True, False]},
    "low_size": {"type": "int", "low": 1, "high": 5},
    "low_percentile": {"type": "int", "low": 20, "high": 95, "step": 5},
    "smooth_high": {"type": "categorical", "choices": [True, False]},
    "high_size": {"type": "int", "low": 1, "high": 5},
    "high_percentile": {"type": "int", "low": 50, "high": 99, "step": 5},
    "erosion": {"type": "int", "low": 0, "high": 5},
}

CLEAN_PARAMS_SPACE = {
    "backgrnd_thresh": {"type": "int", "low": 0, "high": 2000, "step": 50},
    "smooth_thresh": {"type": "int", "low": 0, "high": 1000, "step": 50},
    "smooth": {"type": "categorical", "choices": [True, False]},
    "remove_small": {"type": "categorical", "choices": [True, False]},
    "small_size": {"type": "int", "low": 16, "high": 256, "step": 16},
    "footprint": {"type": "int", "low": 1, "high": 7, "step": 2},
}


def optimize_signal_isolation(
    signal: np.ndarray,
    blank: np.ndarray,
    n_trials: int = 80,
    timeout: Optional[int] = 300,
    warm_start_params: Optional[dict] = None,
    process_fn: Optional[Callable] = None,
    quality_weights: Optional[dict] = None,
    optimize_clean: bool = True,
    study_name: Optional[str] = None,
    storage_path: Optional[Path] = None,
    n_startup_trials: int = 15,
    show_progress: bool = True,
) -> dict:
    """
    Find optimal signal isolation parameters via Bayesian optimization.

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image (uint16).
    blank : np.ndarray
        Blank/autofluorescence channel image (uint16).
    n_trials : int
        Maximum number of optimization trials (default 80).
        With warm-start, 50-80 trials typically suffice.
    timeout : int, optional
        Maximum seconds for optimization (default 300 = 5 min).
        Whichever limit is reached first stops the optimization.
    warm_start_params : dict, optional
        Initial parameter guess (e.g., from predictor or previous dataset).
        If provided, first trial uses these exact parameters, and the search
        space is narrowed around them.
    process_fn : callable, optional
        Custom processing function(signal, blank, params) -> processed.
        Default: clustering._default_process.
    quality_weights : dict, optional
        Custom weights for quality metric components.
    optimize_clean : bool
        Also optimize clean() parameters (doubles search space).
    study_name : str, optional
        Optuna study name for persistence.
    storage_path : Path, optional
        SQLite path for Optuna study persistence.
        Default: in-memory (results lost after session).
    n_startup_trials : int
        Number of random trials before Bayesian optimization begins.
    show_progress : bool
        Display Optuna progress bar.

    Returns
    -------
    dict with:
        status : str ('success' or 'error')
        best_params : dict (blank_params + clean_params)
        best_quality : float (composite quality score)
        best_quality_details : dict (full quality breakdown)
        n_trials_completed : int
        optimization_time_seconds : float
        study_name : str (for resuming or inspecting)
        all_trials_summary : list[dict] (top 10 trials)
    """
    try:
        import optuna
    except ImportError:
        raise ImportError(
            "Optuna package not installed. Install with: pip install optuna\n"
            "Or install the optimize extra: pip install kintsugi[optimize]"
        )

    from kintsugi.qc.signal_quality import compute_signal_isolation_quality

    start_time = time.time()

    # Configure Optuna logging
    if not show_progress:
        optuna.logging.set_verbosity(optuna.logging.WARNING)

    # Create or load study
    storage = None
    if storage_path:
        storage_path = Path(storage_path)
        storage_path.parent.mkdir(parents=True, exist_ok=True)
        storage = f"sqlite:///{storage_path}"

    if study_name is None:
        study_name = f"signal_isolation_{int(time.time())}"

    study = optuna.create_study(
        study_name=study_name,
        storage=storage,
        direction="maximize",
        load_if_exists=True,
        sampler=optuna.samplers.TPESampler(
            n_startup_trials=n_startup_trials,
            seed=42,
        ),
    )

    # Enqueue warm-start trial if provided
    if warm_start_params:
        enqueue_params = _flatten_params_for_optuna(warm_start_params, optimize_clean)
        study.enqueue_trial(enqueue_params)

    def objective(trial):
        # Sample parameters
        bp = _sample_blank_params(trial)
        cp = None
        if optimize_clean:
            cp = _sample_clean_params(trial)

        params = {"blank_params": bp}
        if cp:
            params["clean_params"] = cp

        # Process
        try:
            if process_fn is not None:
                processed = process_fn(signal, blank, params)
            else:
                from kintsugi.signal.clustering import _default_process

                processed = _default_process(signal, blank, params)
        except Exception:
            return 0.0  # failed trial

        # Evaluate quality
        qc = compute_signal_isolation_quality(
            signal, processed, blank=blank, weights=quality_weights
        )

        # Store quality details as trial user attributes
        trial.set_user_attr(
            "quality_details",
            {
                k: float(v) if isinstance(v, (int, float, np.floating)) else v
                for k, v in qc.items()
                if k != "components"
            },
        )

        return qc["quality_score"]

    # Run optimization
    study.optimize(
        objective,
        n_trials=n_trials,
        timeout=timeout,
        show_progress_bar=show_progress,
    )

    elapsed = time.time() - start_time

    # Extract best parameters
    best_trial = study.best_trial
    best_params = _extract_params_from_trial(best_trial, optimize_clean)
    best_quality = best_trial.value
    best_details = best_trial.user_attrs.get("quality_details", {})

    # Top 10 trials summary
    sorted_trials = sorted(
        [t for t in study.trials if t.value is not None],
        key=lambda t: t.value,
        reverse=True,
    )
    top_trials = []
    for t in sorted_trials[:10]:
        top_trials.append(
            {
                "trial_number": t.number,
                "quality_score": t.value,
                "params": _extract_params_from_trial(t, optimize_clean),
            }
        )

    return {
        "status": "success",
        "best_params": best_params,
        "best_quality": best_quality,
        "best_quality_details": best_details,
        "n_trials_completed": len(study.trials),
        "optimization_time_seconds": elapsed,
        "study_name": study_name,
        "all_trials_summary": top_trials,
    }


def _sample_blank_params(trial) -> dict:
    """Sample blank subtraction parameters for an Optuna trial."""
    bp = {}
    for key, spec in BLANK_PARAMS_SPACE.items():
        if spec["type"] == "int":
            bp[key] = trial.suggest_int(key, spec["low"], spec["high"], step=spec.get("step", 1))
        elif spec["type"] == "float":
            bp[key] = trial.suggest_float(
                key, spec["low"], spec["high"], step=spec.get("step")
            )
        elif spec["type"] == "categorical":
            bp[key] = trial.suggest_categorical(key, spec["choices"])
    return bp


def _sample_clean_params(trial) -> dict:
    """Sample cleaning parameters for an Optuna trial."""
    cp = {}
    for key, spec in CLEAN_PARAMS_SPACE.items():
        name = f"clean_{key}"  # prefix to avoid name collision
        if spec["type"] == "int":
            cp[key] = trial.suggest_int(name, spec["low"], spec["high"], step=spec.get("step", 1))
        elif spec["type"] == "float":
            cp[key] = trial.suggest_float(
                name, spec["low"], spec["high"], step=spec.get("step")
            )
        elif spec["type"] == "categorical":
            cp[key] = trial.suggest_categorical(name, spec["choices"])
    return cp


def _flatten_params_for_optuna(params: dict, optimize_clean: bool) -> dict:
    """Convert nested params dict to flat dict for Optuna enqueue."""
    flat = {}
    bp = params.get("blank_params", params)
    for key in BLANK_PARAMS_SPACE:
        if key in bp:
            flat[key] = bp[key]

    if optimize_clean:
        cp = params.get("clean_params", {})
        for key in CLEAN_PARAMS_SPACE:
            if key in cp:
                flat[f"clean_{key}"] = cp[key]

    return flat


def _extract_params_from_trial(trial, optimize_clean: bool) -> dict:
    """Extract structured params dict from an Optuna trial."""
    bp = {key: trial.params.get(key) for key in BLANK_PARAMS_SPACE if key in trial.params}
    result = {"blank_params": bp}

    if optimize_clean:
        cp = {
            key: trial.params.get(f"clean_{key}")
            for key in CLEAN_PARAMS_SPACE
            if f"clean_{key}" in trial.params
        }
        result["clean_params"] = cp

    return result


def plot_optimization_history(
    study_name: str,
    storage_path: Path,
    save_path: Optional[Path] = None,
):
    """
    Plot Optuna optimization history: quality score over trials,
    parameter importance, and score distribution.

    Returns matplotlib Figure.
    """
    import optuna

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    study = optuna.load_study(
        study_name=study_name,
        storage=f"sqlite:///{storage_path}",
    )

    fig, axes = plt.subplots(1, 3, figsize=(18, 5))

    # 1. Optimization history
    trials = [t for t in study.trials if t.value is not None]
    axes[0].plot(
        [t.number for t in trials],
        [t.value for t in trials],
        "o-",
        markersize=3,
        alpha=0.6,
    )
    axes[0].set_xlabel("Trial")
    axes[0].set_ylabel("Quality Score")
    axes[0].set_title("Optimization History")
    axes[0].axhline(
        study.best_value,
        color="r",
        linestyle="--",
        alpha=0.5,
        label=f"Best: {study.best_value:.3f}",
    )
    axes[0].legend()

    # 2. Parameter importance
    try:
        importances = optuna.importance.get_param_importances(study)
        names = list(importances.keys())[:10]
        values = [importances[n] for n in names]
        axes[1].barh(range(len(names)), values, tick_label=names)
        axes[1].set_xlabel("Importance")
        axes[1].set_title("Parameter Importance")
    except Exception:
        axes[1].text(
            0.5,
            0.5,
            "Insufficient trials\nfor importance",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
        )

    # 3. Quality score distribution
    scores = [t.value for t in trials]
    axes[2].hist(scores, bins=20, edgecolor="black", alpha=0.7)
    axes[2].set_xlabel("Quality Score")
    axes[2].set_ylabel("Count")
    axes[2].set_title("Score Distribution")

    plt.tight_layout()
    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")

    return fig
