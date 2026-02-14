"""
Parameter Learning Tools for KINTSUGI MCP Server

Provides tools for recording and recommending processing parameters
based on tissue type and marker name.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Any

import numpy as np

# Import session state from signal_isolation
from kintsugi.mcp.tools.signal_isolation import _loaded_images

if TYPE_CHECKING:
    from kintsugi.claude.parameter_learning import ImageCharacteristics

logger = logging.getLogger("kintsugi.mcp.learning")

# Global learning engine instance
_learning_engine = None


def _get_learning_engine(project_path: str | None = None):
    """Get or create the learning engine."""
    global _learning_engine

    if _learning_engine is None or (
        project_path and str(_learning_engine.project_path) != project_path
    ):
        from kintsugi.claude.parameter_learning import ParameterLearningEngine

        _learning_engine = ParameterLearningEngine(project_path)

    return _learning_engine


async def get_learned_parameters(
    tissue_type: str,
    marker_name: str,
    operation: str,
    project_path: str | None = None,
    channel: str | None = None,
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
    tissue_type: str,
    marker_name: str,
    operation: str,
    parameters: dict[str, Any],
    quality_before: float = 0.0,
    quality_after: float = 0.0,
    user_approved: bool = True,
    user_notes: str = "",
    channel: str | None = None,
    project_path: str | None = None,
    algorithm_version: str | None = None,
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

    Returns:
        Confirmation with record ID
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
                step = max(1, data.shape[0] // 1000)
                sample = data[::step, ::step].compute()
            else:
                sample = data[::10, ::10]
            characteristics = ImageCharacteristics.from_image(sample)
        except Exception as e:
            logger.warning(f"Failed to compute characteristics: {e}")

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
    project_path: str | None = None,
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
    channel: str,
    tissue_type: str,
    marker_name: str,
    blank_channel: str | None = None,
    project_path: str | None = None,
) -> dict[str, Any]:
    """
    Suggest processing parameters combining image analysis and learned history.

    This is the main recommendation function that:
    1. Analyzes the current image
    2. Queries learned parameters for tissue/marker
    3. Combines both sources with confidence weighting

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
            ImageCharacteristics,
            normalize_marker_name,
            normalize_tissue_type,
        )
    except ImportError as e:
        return {"error": f"Learning module not available: {e}"}

    engine = _get_learning_engine(project_path)

    # Get image data
    data = _loaded_images[channel]["data"]
    if hasattr(data, "compute"):
        step = max(1, data.shape[0] // 2000)
        sample = data[::step, ::step].compute()
    else:
        sample = data

    # Compute characteristics
    characteristics = ImageCharacteristics.from_image(sample)

    # Analyze image for heuristic suggestions
    p1, p10, p50, p90, p99 = np.percentile(sample, [1, 10, 50, 90, 99])
    dynamic_range = p99 - p1
    noise_estimate = float(np.std(sample[sample < p10]))

    # Get learned parameters for each operation
    operations = ["blank_subtraction", "denoise", "clahe", "clean_background", "gaussian_subtract"]

    recommendations = {
        "channel": channel,
        "tissue_type": tissue_type,
        "tissue_type_normalized": normalize_tissue_type(tissue_type),
        "marker_name": marker_name,
        "marker_name_normalized": normalize_marker_name(marker_name),
        "image_characteristics": {
            "snr": round(characteristics.snr, 2),
            "dynamic_range": round(characteristics.dynamic_range, 0),
            "sparsity": round(characteristics.sparsity, 3),
            "background_level": round(characteristics.background_level, 0),
        },
        "operations": {},
    }

    for operation in operations:
        # Get learned parameters
        learned = engine.recommend_parameters(
            tissue_type=tissue_type,
            marker_name=marker_name,
            operation=operation,
            characteristics=characteristics,
        )

        # Get heuristic suggestion
        heuristic = _get_heuristic_suggestion(
            operation,
            sample,
            p10,
            p90,
            dynamic_range,
            noise_estimate,
            blank_channel,
        )

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
            "recommended_parameters": combined_params,
            "confidence": round(confidence, 3),
            "source": source,
            "rationale": heuristic.get("rationale", ""),
            "learned_records": learned.get("based_on_records", 0),
        }

    # Suggest processing order based on image
    recommendations["suggested_order"] = _suggest_processing_order(
        characteristics,
        recommendations["operations"],
    )

    return recommendations


def _get_heuristic_suggestion(
    operation: str,
    sample: np.ndarray,
    p10: float,
    p90: float,
    dynamic_range: float,
    noise_estimate: float,
    blank_channel: str | None,
    method: str = "global",
) -> dict[str, Any]:
    """Get heuristic suggestion for an operation based on image analysis."""

    if operation == "blank_subtraction" and method == "weighted":
        # Weighted subtraction heuristic
        if blank_channel and blank_channel in _loaded_images:
            blank_data = _loaded_images[blank_channel]["data"]
            if hasattr(blank_data, "compute"):
                blank_sample = blank_data[::10, ::10].compute()
            else:
                blank_sample = blank_data[::10, ::10]

            try:
                from kintsugi.signal.autofluorescence import (
                    analyze_for_weighted_subtraction,
                )

                analysis = analyze_for_weighted_subtraction(sample, blank_sample)
                return {
                    "parameters": {
                        "method": "weighted",
                        "base_scale_factor": analysis["base_scale_factor"],
                        "blank_clip_factor": analysis["blank_clip_factor"],
                        "n_ranges": analysis["n_ranges"],
                        "range_method": analysis["range_method"],
                        "ranges": analysis["ranges"],
                        "transition_width": analysis["transition_width"],
                        "smooth_low": analysis["smooth_low"],
                        "smooth_high": analysis["smooth_high"],
                        "erosion": analysis["erosion"],
                    },
                    "rationale": (
                        f"Weighted: {analysis['n_ranges']} ranges, "
                        f"dim weight={analysis['ranges'][1]['weight']:.2f} "
                        f"bright weight={analysis['ranges'][-1]['weight']:.2f}"
                        if len(analysis["ranges"]) > 1
                        else "Weighted subtraction"
                    ),
                }
            except Exception as e:
                logger.warning(f"Weighted analysis failed, falling back: {e}")

        return {
            "parameters": {
                "method": "weighted",
                "base_scale_factor": 1.0,
                "blank_clip_factor": 0,
                "n_ranges": 5,
                "range_method": "percentile",
                "transition_width": 0.1,
            },
            "rationale": "Default weighted parameters (no blank loaded)",
        }

    if operation == "blank_subtraction":
        if blank_channel and blank_channel in _loaded_images:
            blank_data = _loaded_images[blank_channel]["data"]
            if hasattr(blank_data, "compute"):
                blank_sample = blank_data[::10, ::10].compute()
            else:
                blank_sample = blank_data[::10, ::10]
            blank_mean = float(np.mean(blank_sample))
            signal_mean = float(np.mean(sample))

            if blank_mean > 0:
                suggested_scale = min(1.0, max(0.5, signal_mean / blank_mean * 0.8))
            else:
                suggested_scale = 0.8

            return {
                "parameters": {
                    "blank_scale_factor": round(suggested_scale, 2),
                    "blank_clip_factor": int(np.percentile(blank_sample, 5)),
                    "smooth_low": noise_estimate > 300,
                    "smooth_high": False,
                    "erosion": 1 if suggested_scale > 0.7 else 0,
                },
                "rationale": f"Blank mean={blank_mean:.0f}, Signal mean={signal_mean:.0f}",
            }
        else:
            return {
                "parameters": {
                    "blank_scale_factor": 0.8,
                    "blank_clip_factor": 0,
                    "smooth_low": False,
                    "smooth_high": False,
                    "erosion": 0,
                },
                "rationale": "No blank channel loaded - using defaults",
            }

    elif operation == "denoise":
        if noise_estimate > 500:
            return {
                "parameters": {"method": "median", "filter_size": 3},
                "rationale": f"High noise (std={noise_estimate:.0f})",
            }
        elif noise_estimate > 200:
            return {
                "parameters": {"method": "median", "filter_size": 2},
                "rationale": f"Moderate noise (std={noise_estimate:.0f})",
            }
        else:
            return {
                "parameters": {"method": "median", "filter_size": 1},
                "rationale": f"Low noise (std={noise_estimate:.0f}) - minimal filtering",
            }

    elif operation == "clahe":
        if dynamic_range < 10000:
            return {
                "parameters": {"clip_limit": 0.02, "tile_grid_size": 50, "nbins": 128},
                "rationale": f"Low dynamic range ({dynamic_range:.0f}) - aggressive CLAHE",
            }
        elif dynamic_range < 30000:
            return {
                "parameters": {"clip_limit": 0.01, "tile_grid_size": 70, "nbins": 128},
                "rationale": f"Moderate dynamic range ({dynamic_range:.0f})",
            }
        else:
            return {
                "parameters": {"clip_limit": 0.005, "tile_grid_size": 100, "nbins": 256},
                "rationale": f"Good dynamic range ({dynamic_range:.0f}) - gentle CLAHE",
            }

    elif operation == "clean_background":
        return {
            "parameters": {
                "background_threshold": int(p10 * 0.5),
                "remove_small_objects": True,
                "min_object_size": 30,
            },
            "rationale": f"Background level at p10={p10:.0f}",
        }

    elif operation == "gaussian_subtract":
        return {
            "parameters": {
                "sigma": 10,
                "scale_factor": 0.1 if dynamic_range > 20000 else 0.2,
            },
            "rationale": "Standard Gaussian subtraction",
        }

    return {"parameters": {}, "rationale": "Unknown operation"}


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
    channel: str,
    tissue_type: str,
    marker_name: str,
    operations_params: dict[str, dict[str, Any]],
    quality_before: float,
    quality_after: float,
    notes: str = "",
    project_path: str | None = None,
) -> dict[str, Any]:
    """
    Approve processing results and record all parameters for learning.

    Call this after user approves the final processed result.

    Args:
        channel: Channel name
        tissue_type: Tissue type
        marker_name: Marker name
        operations_params: Dict mapping operation names to their parameters
        quality_before: Quality score before all processing
        quality_after: Quality score after all processing
        notes: User notes
        project_path: Project path

    Returns:
        Confirmation with learning status
    """
    results = []

    for operation, parameters in operations_params.items():
        result = await record_successful_parameters(
            tissue_type=tissue_type,
            marker_name=marker_name,
            operation=operation,
            parameters=parameters,
            quality_before=quality_before,
            quality_after=quality_after,
            user_approved=True,
            user_notes=notes,
            channel=channel,
            project_path=project_path,
        )
        results.append(
            {
                "operation": operation,
                "recorded": result.get("status") == "success",
            }
        )

    return {
        "status": "success",
        "message": "Processing approved and parameters recorded for learning",
        "tissue_type": tissue_type,
        "marker_name": marker_name,
        "quality_improvement": quality_after - quality_before,
        "operations_recorded": results,
    }
