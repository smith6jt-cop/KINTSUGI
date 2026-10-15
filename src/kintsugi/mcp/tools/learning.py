"""
Parameter Learning Tools for KINTSUGI MCP Server

Provides tools for recording and recommending processing parameters
based on tissue type and marker name.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Annotated, Any, Literal

import numpy as np
from pydantic import Field

from kintsugi.mcp.heuristics import (
    heuristic_suggestion,
    sample_image,
    sample_step,
    summarize_intensity,
)
from kintsugi.mcp.schemas import LearningSuggestion
from kintsugi.mcp.session import _loaded_images

if TYPE_CHECKING:
    from kintsugi.claude.parameter_learning import ImageCharacteristics

logger = logging.getLogger("kintsugi.mcp.learning")


def _plain(value: Any) -> Any:
    """Recursively turn numpy scalars/arrays into JSON-serialisable Python values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


# Learning engines, one per project (key "" = the user-wide database only). A single
# cached instance used to be reused across projects when project_path was omitted,
# which silently wrote records into the previous project's database.
_learning_engines: dict[str, Any] = {}


def _engine_key(project_path: str | None) -> str:
    from pathlib import Path

    return str(Path(project_path).resolve()) if project_path else ""


def _get_learning_engine(project_path: str | None = None):
    """Get or create the learning engine for ``project_path`` (or the global one)."""
    key = _engine_key(project_path)
    engine = _learning_engines.get(key)
    if engine is None:
        from kintsugi.claude.parameter_learning import ParameterLearningEngine

        engine = _learning_engines[key] = ParameterLearningEngine(project_path)
    return engine


def _channel_characteristics(channel: str | None) -> ImageCharacteristics | None:
    """``ImageCharacteristics`` of a loaded channel (sub-sampled), or None."""
    if not channel or channel not in _loaded_images:
        return None
    try:
        from kintsugi.claude.parameter_learning import ImageCharacteristics

        data = _loaded_images[channel]["data"]
        if hasattr(data, "compute"):
            step = max(1, data.shape[0] // 1000)
            sample = data[::step, ::step].compute()
        else:
            sample = data[::10, ::10]
        return ImageCharacteristics.from_image(sample)
    except Exception as e:
        logger.warning(f"Failed to compute characteristics: {e}")
        return None


async def get_learned_parameters(
    tissue_type: Annotated[
        str, Field(description="Type of tissue (e.g., 'tonsil', 'lymph node', 'pancreas')")
    ],
    marker_name: Annotated[str, Field(description="Marker name (e.g., 'CD3', 'CD20', 'DAPI')")],
    operation: Annotated[
        Literal["blank_subtraction", "denoise", "clahe", "clean_background", "gaussian_subtract"],
        Field(description="Processing operation whose learned parameters to look up"),
    ],
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory holding a project-specific learning "
            "database; None consults the user-wide global database only"
        ),
    ] = None,
    channel: Annotated[
        str | None,
        Field(
            description="Name of a loaded channel whose image characteristics (SNR, dynamic "
            "range, sparsity) refine the match; None matches on tissue and marker only"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Get recommended parameters based on learned history.

    Queries the learning database for similar tissue/marker combinations
    and returns weighted recommendations with confidence scores.

    Args:
        tissue_type: Type of tissue (e.g., "tonsil", "lymph node")
        marker_name: Marker name (e.g., "CD3", "DAPI")
        operation: Processing operation (e.g., "blank_subtraction")
        project_path: Project path for project-specific learning
        channel: Channel name to get image characteristics from

    Returns:
        Recommendation with parameters and confidence
    """
    try:
        from kintsugi.claude.parameter_learning import ImageCharacteristics
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    engine = _get_learning_engine(project_path)

    # Get image characteristics if channel is loaded
    characteristics = None
    if channel and channel in _loaded_images:
        try:
            data = _loaded_images[channel]["data"]
            if hasattr(data, "compute"):
                # Sample for speed
                step = max(1, data.shape[0] // 1000)
                sample = data[::step, ::step].compute()
            else:
                sample = data[::10, ::10]
            characteristics = ImageCharacteristics.from_image(sample)
        except Exception as e:
            logger.warning(f"Failed to compute characteristics: {e}")

    # Get recommendations
    recommendations = engine.recommend_parameters(
        tissue_type=tissue_type,
        marker_name=marker_name,
        operation=operation,
        characteristics=characteristics,
    )

    return recommendations


async def record_successful_parameters(
    tissue_type: Annotated[str, Field(description="Type of tissue processed (e.g., 'tonsil')")],
    marker_name: Annotated[str, Field(description="Marker name (e.g., 'CD3', 'DAPI')")],
    operation: Annotated[
        Literal["blank_subtraction", "denoise", "clahe", "clean_background", "gaussian_subtract"],
        Field(description="Processing operation the parameters belong to"),
    ],
    parameters: Annotated[
        dict[str, Any], Field(description="Parameter dictionary that was used for the operation")
    ],
    quality_before: Annotated[
        float,
        Field(
            description="Quality score before processing on a 0.0-1.0 scale",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ] = 0.0,
    quality_after: Annotated[
        float,
        Field(
            description="Quality score after processing on a 0.0-1.0 scale",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ] = 0.0,
    user_approved: Annotated[
        bool,
        Field(
            description="Whether the user approved the result; approved records carry full "
            "weight in recommendations, unapproved (bootstrap) records are included at "
            "reduced weight"
        ),
    ] = True,
    user_notes: Annotated[str, Field(description="Optional free-text notes from the user")] = "",
    channel: Annotated[
        str | None,
        Field(
            description="Name of a loaded channel whose image characteristics are stored with "
            "the record; None stores the record without characteristics"
        ),
    ] = None,
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory holding a project-specific learning "
            "database; None records to the user-wide global database only"
        ),
    ] = None,
    algorithm_version: Annotated[
        str | None,
        Field(
            description="Free-form tag of the algorithm variant, e.g. 'weighted_v1'; None "
            "falls back to parameters['algorithm_version'] when present"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Record a successful parameter set for future learning.

    Should be called after user approves processing results.

    Args:
        tissue_type: Type of tissue processed
        marker_name: Marker name
        operation: Processing operation
        parameters: Parameter dictionary that was used
        quality_before: Quality score before processing
        quality_after: Quality score after processing
        user_approved: Whether user approved the result
        user_notes: Optional notes from user
        channel: Channel name (for characteristics)
        project_path: Project path
        algorithm_version: Free-form tag of the algorithm variant (e.g. "weighted_v1")

    Returns:
        Confirmation with record ID
    """
    try:
        from kintsugi.claude.parameter_learning import ImageCharacteristics  # noqa: F401
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    engine = _get_learning_engine(project_path)

    characteristics = _channel_characteristics(channel)

    # Record parameters
    record_id = engine.record_parameters(
        tissue_type=tissue_type,
        marker_name=marker_name,
        operation=operation,
        parameters=parameters,
        quality_before=quality_before,
        quality_after=quality_after,
        user_approved=user_approved,
        user_notes=user_notes,
        characteristics=characteristics,
        channel_name=channel or "",
        algorithm_version=algorithm_version or parameters.get("algorithm_version"),
    )

    return {
        "status": "success",
        "record_id": record_id,
        "message": f"Recorded parameters for {tissue_type}/{marker_name}/{operation}",
        "quality_improvement": quality_after - quality_before,
        "will_influence_future": user_approved,
    }


async def get_learning_statistics(
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory holding a project-specific learning "
            "database to summarize; None summarizes the user-wide global database only"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Get statistics about the parameter learning database.

    Returns counts of records, unique tissues/markers, and operation stats.
    """
    try:
        from kintsugi.claude.parameter_learning import ParameterLearningEngine  # noqa: F401
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    engine = _get_learning_engine(project_path)
    stats = engine.get_statistics()

    return {
        "status": "success",
        "statistics": stats,
    }


async def suggest_with_learning(
    channel: Annotated[str, Field(description="Name of the loaded channel to analyze")],
    tissue_type: Annotated[str, Field(description="Tissue type (e.g., 'tonsil', 'lymph node')")],
    marker_name: Annotated[str, Field(description="Marker name (e.g., 'CD3', 'DAPI')")],
    blank_channel: Annotated[
        str | None,
        Field(
            description="Name of a loaded blank/autofluorescence channel used to estimate "
            "blank-subtraction parameters; None uses default blank-subtraction heuristics"
        ),
    ] = None,
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory holding a project-specific learning "
            "database; None consults the user-wide global database only"
        ),
    ] = None,
) -> LearningSuggestion:
    """
    Suggest processing parameters combining image analysis and learned history.

    This is the main recommendation function that:
    1. Analyzes the current image
    2. Queries learned parameters for tissue/marker
    3. Combines both sources with confidence weighting

    The return annotation is the advertised ``outputSchema``: keys outside it
    are dropped on the wire and non-conforming values turn the whole result
    into an opaque error, so every value is coerced to a plain Python type.

    Args:
        channel: Loaded channel name
        tissue_type: Tissue type
        marker_name: Marker name
        blank_channel: Blank channel for AF estimation
        project_path: Project path

    Returns:
        Comprehensive recommendations for all operations
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

    try:
        from kintsugi.claude.parameter_learning import (
            LEARNING_OPERATIONS,
            ImageCharacteristics,
            normalize_marker_name,
            normalize_tissue_type,
        )
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    engine = _get_learning_engine(project_path)

    data = _loaded_images[channel]["data"]
    step = sample_step(data)
    sample = sample_image(data, step=step)
    characteristics = ImageCharacteristics.from_image(sample)
    summary = summarize_intensity(sample)
    blank_sample = None
    if blank_channel:
        if blank_channel not in _loaded_images:
            return {"error": f"Blank channel not loaded: {blank_channel}"}
        blank_sample = sample_image(_loaded_images[blank_channel]["data"], step=step)

    # Get learned parameters for each operation
    operations = list(LEARNING_OPERATIONS)

    recommendations: dict[str, Any] = {
        "status": "success",
        "channel": channel,
        "tissue_type": tissue_type,
        "tissue_type_normalized": normalize_tissue_type(tissue_type),
        "marker_name": marker_name,
        "marker_name_normalized": normalize_marker_name(marker_name),
        "image_characteristics": {
            "snr": round(float(characteristics.snr), 2),
            "dynamic_range": round(float(characteristics.dynamic_range), 0),
            "sparsity": round(float(characteristics.sparsity), 3),
            "background_level": round(float(characteristics.background_level), 0),
        },
        "operations": {},
    }

    # One query per database for all operations instead of one per operation.
    learned_by_operation = engine.recommend_parameters_many(
        tissue_type=tissue_type,
        marker_name=marker_name,
        operations=operations,
        characteristics=characteristics,
    )

    for operation in operations:
        learned = learned_by_operation[operation]
        heuristic = heuristic_suggestion(operation, sample, summary, blank_sample)

        # Combine learned and heuristic
        if learned["status"] == "success" and learned["confidence"] > 0.5:
            # High confidence from learning - use learned with heuristic fallback
            combined_params = learned["recommended_parameters"].copy()
            source = "learned"
            confidence = learned["confidence"]

            # Fill in any missing params from heuristic
            for key, value in heuristic.get("parameters", {}).items():
                if key not in combined_params:
                    combined_params[key] = value

        elif learned["status"] == "success" and learned["confidence"] > 0.3:
            # Medium confidence - blend learned and heuristic
            combined_params = {}
            learned_params = learned["recommended_parameters"] or {}
            heuristic_params = heuristic.get("parameters", {})

            # Weighted average for numeric params
            for key in set(learned_params.keys()) | set(heuristic_params.keys()):
                l_val = learned_params.get(key)
                h_val = heuristic_params.get(key)

                if l_val is not None and h_val is not None:
                    if isinstance(l_val, (int, float)) and isinstance(h_val, (int, float)):
                        # Blend with learned weight based on confidence
                        w = learned["confidence"]
                        combined_params[key] = l_val * w + h_val * (1 - w)
                        if isinstance(l_val, int) and isinstance(h_val, int):
                            combined_params[key] = int(round(combined_params[key]))
                    else:
                        # Non-numeric: prefer learned if confident enough
                        combined_params[key] = l_val if learned["confidence"] > 0.4 else h_val
                elif l_val is not None:
                    combined_params[key] = l_val
                else:
                    combined_params[key] = h_val

            source = "blended"
            confidence = (learned["confidence"] + 0.5) / 2  # Boost slightly

        else:
            # No good learned data - use heuristic only
            combined_params = heuristic.get("parameters", {})
            source = "heuristic"
            confidence = 0.3

        recommendations["operations"][operation] = {
            "recommended_parameters": _plain(combined_params),
            "confidence": round(float(confidence), 3),
            "source": source,
            "rationale": str(heuristic.get("rationale", "")),
            "learned_records": int(learned.get("based_on_records", 0) or 0),
        }

    # Suggest processing order based on image
    recommendations["suggested_order"] = [
        str(op) for op in _suggest_processing_order(characteristics, recommendations["operations"])
    ]

    return recommendations


def _suggest_processing_order(
    characteristics: ImageCharacteristics,
    operations: dict[str, Any],
) -> list[str]:
    """Suggest optimal processing order based on image characteristics."""
    order = []

    # Always start with blank subtraction if available and confident
    if "blank_subtraction" in operations:
        if operations["blank_subtraction"]["confidence"] > 0.3:
            order.append("blank_subtraction")

    # Denoise early if noisy
    if characteristics.snr < 5:
        order.append("denoise")

    # Gaussian subtract for structured background
    if characteristics.background_level > 500:
        order.append("gaussian_subtract")

    # CLAHE for contrast
    if characteristics.dynamic_range < 30000:
        order.append("clahe")

    # Clean at the end
    order.append("clean_background")

    # Add any operations not yet included
    for op in operations:
        if op not in order:
            order.append(op)

    return order


async def approve_and_learn(
    channel: Annotated[str, Field(description="Name of the loaded channel that was processed")],
    tissue_type: Annotated[str, Field(description="Tissue type (e.g., 'tonsil', 'lymph node')")],
    marker_name: Annotated[str, Field(description="Marker name (e.g., 'CD3', 'DAPI')")],
    operations_params: Annotated[
        dict[str, dict[str, Any]],
        Field(
            description="Mapping of operation name to the parameter dict that was used; keys may "
            "be learning operations (blank_subtraction, denoise, clahe, clean_background, "
            "gaussian_subtract) or the tool names recorded in the session history "
            "(subtract_blank, subtract_blank_weighted, denoise_advanced, ...)"
        ),
    ],
    quality_before: Annotated[
        float,
        Field(
            description="Quality score before all processing on a 0.0-1.0 scale",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ],
    quality_after: Annotated[
        float,
        Field(
            description="Quality score after all processing on a 0.0-1.0 scale",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ],
    notes: Annotated[str, Field(description="Free-text user notes stored with every record")] = "",
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory holding a project-specific learning "
            "database; None records to the user-wide global database only"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Approve processing results and record all parameters for learning.

    Call this after user approves the final processed result.

    Args:
        channel: Channel name
        tissue_type: Tissue type
        marker_name: Marker name
        operations_params: Dict mapping operation names to their parameters. Keys may
            be learning operations or session-history tool names (e.g. "subtract_blank"),
            which are mapped through HISTORY_TO_LEARNING_OP before recording.
        quality_before: Quality score before all processing
        quality_after: Quality score after all processing
        notes: User notes
        project_path: Project path

    Returns:
        Confirmation with learning status
    """
    try:
        from kintsugi.claude.parameter_learning import (
            HISTORY_TO_LEARNING_OP,
            LEARNING_OPERATIONS,
        )
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    # The session history records tool names (e.g. "subtract_blank") while the
    # learning database is keyed by operation names; map every key before any
    # record is written so an unknown key cannot leave a partial set behind.
    mapped_operations: list[tuple[str, dict[str, Any]]] = []
    for operation, parameters in operations_params.items():
        learning_op = HISTORY_TO_LEARNING_OP.get(operation, operation)
        if learning_op not in LEARNING_OPERATIONS:
            return {
                "error": (
                    f"Unknown operation '{operation}'. Allowed operations: "
                    f"{', '.join(LEARNING_OPERATIONS)}; history tool names "
                    f"{', '.join(sorted(HISTORY_TO_LEARNING_OP))} are also accepted"
                )
            }
        mapped_operations.append((learning_op, parameters))

    # All operations are written in one transaction per database, with the image
    # characteristics computed once rather than once per operation.
    engine = _get_learning_engine(project_path)
    record_ids = engine.record_parameters_many(
        [{"operation": op, "parameters": params} for op, params in mapped_operations],
        tissue_type=tissue_type,
        marker_name=marker_name,
        quality_before=quality_before,
        quality_after=quality_after,
        user_approved=True,
        user_notes=notes,
        characteristics=_channel_characteristics(channel),
        channel_name=channel or "",
    )
    results = [
        {"operation": op, "recorded": True, "record_id": record_id}
        for (op, _params), record_id in zip(mapped_operations, record_ids, strict=True)
    ]

    return {
        "status": "success",
        "message": "Processing approved and parameters recorded for learning",
        "tissue_type": tissue_type,
        "marker_name": marker_name,
        "quality_improvement": quality_after - quality_before,
        "operations_recorded": results,
    }
