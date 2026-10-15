"""
KINTSUGI assisted tuning and parameter learning.

Shared by the MCP server (any MCP client), the Jupyter notebooks and the batch
pipeline:

* ``param_tuner`` / ``artifact_tuner`` - interactive Jupyter widgets for tuning
  signal-isolation parameters by eye;
* ``parameter_learning`` - the SQLite-backed engine that records approved
  parameters per tissue and marker and recommends them for new channels;
* ``workflow_state`` - per-project processing status and audit trail;
* ``state_dir`` - where a project's state lives (``<project>/.kintsugi/``).

Names are resolved lazily, so ``from kintsugi.assist import
ParameterLearningEngine`` does not import ipywidgets, IPython or stackview.

Usage in Jupyter::

    from kintsugi.assist import blank_subtraction_tuner

    result = blank_subtraction_tuner(
        signal_channel="cyc001_CD3",
        blank_channel="cyc001_blank",
    )
    params = result.get_params()

Parameter learning::

    from kintsugi.assist import ParameterLearningEngine

    engine = ParameterLearningEngine("/path/to/project")
    rec = engine.recommend_parameters(
        tissue_type="tonsil",
        marker_name="CD3",
        operation="blank_subtraction",
    )
    print(rec["recommended_parameters"], rec["confidence"])
"""

from __future__ import annotations

import importlib
from typing import Any

# Public name -> submodule that defines it.
_EXPORTS = {
    # Parameter tuners (Jupyter widgets)
    "blank_subtraction_tuner": "param_tuner",
    "denoise_tuner": "param_tuner",
    "clahe_tuner": "param_tuner",
    "clean_tuner": "param_tuner",
    "gaussian_subtract_tuner": "param_tuner",
    # Workflow state
    "WorkflowState": "workflow_state",
    # Parameter learning
    "ParameterLearningEngine": "parameter_learning",
    "ImageCharacteristics": "parameter_learning",
    "ParameterRecord": "parameter_learning",
    "normalize_tissue_type": "parameter_learning",
    "normalize_marker_name": "parameter_learning",
    # Project state directory
    "project_state_dir": "state_dir",
    "STATE_DIR_NAME": "state_dir",
    "LEGACY_STATE_DIR_NAME": "state_dir",
}

__all__ = list(_EXPORTS)


def __getattr__(name: str) -> Any:
    module = _EXPORTS.get(name)
    if module is None:
        # AttributeError (not ImportError) lets ``from kintsugi.assist import
        # param_tuner`` fall back to importing the submodule.
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(importlib.import_module(f"{__name__}.{module}"), name)
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
