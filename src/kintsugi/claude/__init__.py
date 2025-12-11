"""
KINTSUGI Claude Integration

Jupyter bridge for Claude-guided image processing workflows.
Provides interactive parameter tuning widgets that integrate with
the MCP server tools.

Usage in Jupyter:
    from kintsugi.claude import param_tuner

    # Interactive blank subtraction tuning
    result = param_tuner.blank_subtraction_tuner(
        signal_channel="cyc001_CD3",
        blank_channel="cyc001_blank",
    )

    # Get optimized parameters
    params = result.get_params()

Parameter Learning:
    from kintsugi.claude import ParameterLearningEngine

    # Initialize learning engine
    engine = ParameterLearningEngine("/path/to/project")

    # Get recommendations based on tissue/marker
    rec = engine.recommend_parameters(
        tissue_type="tonsil",
        marker_name="CD3",
        operation="blank_subtraction",
    )
    print(rec["recommended_parameters"])
    print(f"Confidence: {rec['confidence']}")
"""

from kintsugi.claude.param_tuner import (
    blank_subtraction_tuner,
    denoise_tuner,
    clahe_tuner,
    clean_tuner,
    gaussian_subtract_tuner,
)

from kintsugi.claude.workflow_state import WorkflowState

from kintsugi.claude.parameter_learning import (
    ParameterLearningEngine,
    ImageCharacteristics,
    ParameterRecord,
    normalize_tissue_type,
    normalize_marker_name,
)

__all__ = [
    # Parameter tuners
    "blank_subtraction_tuner",
    "denoise_tuner",
    "clahe_tuner",
    "clean_tuner",
    "gaussian_subtract_tuner",
    # Workflow state
    "WorkflowState",
    # Parameter learning
    "ParameterLearningEngine",
    "ImageCharacteristics",
    "ParameterRecord",
    "normalize_tissue_type",
    "normalize_marker_name",
]
