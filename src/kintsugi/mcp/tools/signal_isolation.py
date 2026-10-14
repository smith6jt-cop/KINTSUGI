"""
Signal Isolation Tools for KINTSUGI MCP Server

Wraps Kutils functions for blank subtraction, denoising, and signal enhancement.
"""

from __future__ import annotations

import asyncio
import logging
from collections.abc import Callable
from pathlib import Path
from typing import Annotated, Any, Literal

import numpy as np
from pydantic import Field

from kintsugi.mcp._sdk import Context
from kintsugi.mcp.path_safety import (
    PathSafetyError,
    allowed_model_roots,
    clamp_float,
    clamp_int,
    ensure_within,
    ensure_within_any,
    safe_filename,
    validate_choice,
)

logger = logging.getLogger("kintsugi.mcp.signal_isolation")

# Session state - shared with server
_loaded_images: dict[str, Any] = {}
_processing_history: dict[str, list] = {}


def _get_image(name: str) -> Any:
    """Get a loaded image by name."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


def _store_image(name: str, data: Any, metadata: dict[str, Any] = None):
    """Store an image with metadata."""
    _loaded_images[name] = {
        "data": data,
        "metadata": metadata or {},
        "dtype": str(data.dtype) if hasattr(data, "dtype") else "unknown",
        "shape": data.shape if hasattr(data, "shape") else None,
    }
    if name not in _processing_history:
        _processing_history[name] = []


def _add_history(name: str, operation: str, params: dict[str, Any]):
    """Add an operation to processing history."""
    if name not in _processing_history:
        _processing_history[name] = []
    _processing_history[name].append(
        {
            "operation": operation,
            "params": params,
        }
    )


def _as_numpy(data: Any) -> np.ndarray:
    """Materialize a loaded image (dask or array-like) as a NumPy array.

    Runs the (possibly heavy) ``compute()`` so callers can push it onto a
    worker thread with :func:`asyncio.to_thread`.
    """
    if hasattr(data, "compute"):
        return data.compute()
    return np.array(data)


async def _report_progress(
    ctx: Context | None, progress: float, total: float, message: str
) -> None:
    """Send a progress notification when the call carries a request context.

    Direct in-process calls pass no ``ctx`` and get no notifications. The SDK
    itself drops the notification when the client sent no progress token, so
    this is safe to call unconditionally inside a handler.
    """
    if ctx is None:
        return
    try:
        await ctx.report_progress(progress, total, message)
    except Exception:  # progress is advisory; never fail the tool over it
        logger.debug("progress notification failed", exc_info=True)


def _log_progress_failure(future: Any) -> None:
    if not future.cancelled() and future.exception() is not None:
        logger.debug("progress notification failed: %s", future.exception())


def _thread_progress_reporter(
    ctx: Context | None, loop: asyncio.AbstractEventLoop
) -> Callable[[float, float, str], None]:
    """Build a plain callable that a worker thread can use to report progress.

    ``ctx.report_progress`` is a coroutine that must run on the event loop, so
    the callable hands it back to ``loop`` with
    :func:`asyncio.run_coroutine_threadsafe` instead of awaiting in the thread.
    """

    def report(progress: float, total: float, message: str) -> None:
        if ctx is None:
            return
        try:
            future = asyncio.run_coroutine_threadsafe(
                ctx.report_progress(progress, total, message), loop
            )
        except RuntimeError:  # loop already closed (request cancelled/finished)
            return
        future.add_done_callback(_log_progress_failure)

    return report


async def load_channel(
    cycle: Annotated[str, Field(description="Cycle name or number (e.g., 'cyc001' or '1')")],
    channel: Annotated[str, Field(description="Channel name (e.g., 'CD3', 'DAPI')")],
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project directory; None = the Claude Code project "
            "(CLAUDE_PROJECT_DIR) or the current directory"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Load a channel image from a KINTSUGI project.

    Returns image metadata and prepares it for processing. When ``project_path``
    is omitted the server's own project (``CLAUDE_PROJECT_DIR``, else the
    working directory) is used.
    """
    # Imported lazily so this module imports without the session helper.
    from kintsugi.mcp.session import resolve_project_path

    project_path = resolve_project_path(project_path)

    # cycle and channel become path/glob components below, so reject any value
    # that could traverse out of the project tree or broaden the glob search.
    # Validate cheap inputs before importing heavy optional dependencies.
    try:
        cycle = safe_filename(str(cycle), "cycle")
        channel = safe_filename(channel, "channel")
    except (PathSafetyError, ValueError) as e:
        return {"error": str(e)}
    for _value, _name in ((cycle, "cycle"), (channel, "channel")):
        if any(ch in _value for ch in "*?[]"):
            return {"error": f"{_name} must not contain glob wildcards (*?[]): {_value!r}"}

    try:
        import dask.array as da  # noqa: F401
        import tifffile  # noqa: F401
        import zarr  # noqa: F401
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    project_path = Path(project_path).resolve()

    # Try to load from KINTSUGI project structure
    try:
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_path)
        raw_dir = project.paths.raw
        processed_dir = project.paths.processed
    except Exception:
        # Fallback to direct path
        raw_dir = project_path / "data" / "raw"
        processed_dir = project_path / "data" / "processed"

    # Normalize cycle name
    if cycle.isdigit():
        cycle_name = f"cyc{int(cycle):03d}"
    else:
        cycle_name = cycle

    # Search for the channel file
    search_paths = [
        raw_dir / cycle_name,
        processed_dir / "registered" / cycle_name,
        processed_dir / "stitched" / cycle_name,
        project_path / cycle_name,
    ]

    image_path = None
    for search_dir in search_paths:
        if not search_dir.exists():
            continue

        # Look for files matching the channel name
        patterns = [
            f"*{channel}*.tif*",
            f"*{channel}*.zarr",
            f"*{channel.upper()}*.tif*",
            f"*{channel.lower()}*.tif*",
        ]

        for pattern in patterns:
            matches = list(search_dir.glob(pattern))
            if matches:
                image_path = matches[0]
                break

        if image_path:
            break

    if image_path is None:
        return {
            "error": f"Channel '{channel}' not found in cycle '{cycle_name}'",
            "searched_paths": [str(p) for p in search_paths if p.exists()],
        }

    # Guard against a symlinked match resolving outside the project tree.
    try:
        ensure_within(image_path, project_path, "channel file")
    except PathSafetyError as e:
        return {"error": str(e)}

    # Load the image
    try:
        if str(image_path).endswith(".zarr"):
            # Load as zarr
            z = zarr.open(str(image_path), mode="r")
            data = da.from_zarr(z)
        else:
            # Load as TIFF - use dask for lazy loading of large files
            with tifffile.TiffFile(str(image_path)) as tif:
                if tif.is_ome:
                    # OME-TIFF - may be multichannel
                    data = tifffile.imread(str(image_path))
                    if data.ndim > 2:
                        # Need to extract the correct channel
                        logger.info(f"Multichannel image detected: {data.shape}")
                else:
                    data = tifffile.imread(str(image_path))

            # Convert to dask array for efficient processing
            data = da.from_array(data, chunks=(1000, 1000))

    except Exception as e:
        return {"error": f"Failed to load image: {e}"}

    # Create unique name for this loaded image
    image_name = f"{cycle_name}_{channel}"

    # Store the image
    _store_image(
        image_name,
        data,
        {
            "source_path": str(image_path),
            "cycle": cycle_name,
            "channel": channel,
            "project": str(project_path),
        },
    )

    # Compute basic stats (sample for large images)
    try:
        if hasattr(data, "compute"):
            # Dask array - compute stats on a sample
            sample = data[::10, ::10].compute()
        else:
            sample = data[::10, ::10]

        stats = {
            "min": float(np.min(sample)),
            "max": float(np.max(sample)),
            "mean": float(np.mean(sample)),
            "std": float(np.std(sample)),
        }
    except Exception:
        stats = {}

    return {
        "status": "loaded",
        "name": image_name,
        "shape": list(data.shape),
        "dtype": str(data.dtype),
        "source": str(image_path),
        "stats": stats,
    }


async def suggest_subtraction_parameters(
    signal_channel: str,
    blank_channel: str,
    tissue_type: str = "unknown",
    marker_name: str | None = None,
    use_learning: bool = True,
    project_path: str | None = None,
) -> dict[str, Any]:
    """
    Analyze images and suggest optimal autofluorescence subtraction parameters.

    This uses intelligent analysis of the signal and blank images combined
    with learned parameters from previous successful runs (if available).

    Parameters:
    - signal_channel: Name of loaded signal channel
    - blank_channel: Name of loaded blank channel
    - tissue_type: Type of tissue (e.g., 'tonsil', 'skin', 'lymph_node')
    - marker_name: Name of the marker (e.g., 'CD3', 'CD20')
    - use_learning: Whether to use learned parameters from history
    - project_path: Path to project for parameter learning

    Returns:
    - Suggested parameters with confidence scores and analysis details
    """
    try:
        signal_data = _get_image(signal_channel)
        blank_data = _get_image(blank_channel)
    except ValueError as e:
        return {"error": str(e)}

    # Convert to numpy for analysis
    if hasattr(signal_data, "compute"):
        # Sample for large images
        if signal_data.shape[0] > 2000 or signal_data.shape[1] > 2000:
            signal_np = signal_data[::2, ::2].compute()
            blank_np = blank_data[::2, ::2].compute()
        else:
            signal_np = signal_data.compute()
            blank_np = blank_data.compute()
    else:
        signal_np = np.array(signal_data)
        blank_np = np.array(blank_data)

    # Use the intelligent analysis
    try:
        from kintsugi.signal import analyze_for_subtraction
    except ImportError:
        # Fallback to basic analysis
        logger.warning("kintsugi.signal module not available, using basic analysis")
        return _basic_parameter_suggestion(signal_np, blank_np)

    # Get analysis-based suggestions
    analysis = analyze_for_subtraction(
        signal_np,
        blank_np,
        tissue_type=tissue_type,
        marker_name=marker_name,
    )

    # Try to get learned parameters
    learned_params = None
    learned_confidence = 0
    if use_learning and project_path and marker_name:
        try:
            from kintsugi.claude.parameter_learning import ParameterLearningEngine

            engine = ParameterLearningEngine(project_path)
            recommendation = engine.recommend_parameters(
                operation="blank_subtraction",
                tissue_type=tissue_type,
                marker_name=marker_name,
            )
            if recommendation.get("status") == "success":
                learned_params = recommendation.get("recommended_parameters", {})
                learned_confidence = recommendation.get("confidence", 0)
        except Exception as e:
            logger.warning(f"Failed to get learned parameters: {e}")

    # Combine analysis and learned parameters
    if learned_params and learned_confidence > 0.5:
        # Weighted merge
        weight = learned_confidence * 0.6  # learned params get up to 60% weight
        suggested = {
            "blank_clip_factor": int(
                analysis["blank_clip_factor"] * (1 - weight)
                + learned_params.get("blank_clip_factor", analysis["blank_clip_factor"]) * weight
            ),
            "blank_scale_factor": round(
                analysis["blank_scale_factor"] * (1 - weight)
                + learned_params.get("blank_scale_factor", analysis["blank_scale_factor"]) * weight,
                2,
            ),
            "smooth_low": learned_params.get("smooth_low", analysis["smooth_low"]),
            "low_size": learned_params.get("low_size", analysis["low_size"]),
            "low_percentile": learned_params.get("low_percentile", analysis["low_percentile"]),
            "smooth_high": learned_params.get("smooth_high", analysis["smooth_high"]),
            "high_size": learned_params.get("high_size", analysis["high_size"]),
            "high_percentile": learned_params.get("high_percentile", analysis["high_percentile"]),
            "erosion": learned_params.get("erosion", analysis["erosion"]),
        }
        source = "analysis+learned"
    else:
        suggested = {
            "blank_clip_factor": analysis["blank_clip_factor"],
            "blank_scale_factor": analysis["blank_scale_factor"],
            "smooth_low": analysis["smooth_low"],
            "low_size": analysis["low_size"],
            "low_percentile": analysis["low_percentile"],
            "smooth_high": analysis["smooth_high"],
            "high_size": analysis["high_size"],
            "high_percentile": analysis["high_percentile"],
            "erosion": analysis["erosion"],
        }
        source = "analysis"

    return {
        "status": "success",
        "suggested_parameters": suggested,
        "confidence": analysis["confidence"],
        "source": source,
        "analysis": {
            "signal_blank_correlation": analysis["analysis"]["correlation"],
            "autofluorescence_contribution": analysis["analysis"]["af_contribution"],
            "signal_noise_ratio": analysis["analysis"]["signal_noise_ratio"],
        },
        "learned_parameters": (
            {
                "available": learned_params is not None,
                "confidence": learned_confidence,
            }
            if use_learning
            else None
        ),
        "tissue_type": tissue_type,
        "marker_name": marker_name,
    }


def _basic_parameter_suggestion(signal: np.ndarray, blank: np.ndarray) -> dict:
    """Fallback basic parameter suggestion without the signal module."""
    # Simple correlation
    flat_s = signal.ravel().astype(np.float64)
    flat_b = blank.ravel().astype(np.float64)
    mask = (flat_s > 0) | (flat_b > 0)
    corr = np.corrcoef(flat_s[mask], flat_b[mask])[0, 1] if np.sum(mask) > 100 else 0

    # Basic suggestions
    clip = int(np.percentile(blank, 10))
    scale = 1.0 + max(0, corr - 0.5) * 0.8

    return {
        "status": "success",
        "suggested_parameters": {
            "blank_clip_factor": clip,
            "blank_scale_factor": round(scale, 2),
            "smooth_low": False,
            "low_size": 2,
            "low_percentile": 60,
            "smooth_high": False,
            "high_size": 2,
            "high_percentile": 90,
            "erosion": 1 if corr > 0.6 else 0,
        },
        "confidence": 0.4,
        "source": "basic_analysis",
        "analysis": {"correlation": round(corr, 4)},
    }


async def analyze_weighted_subtraction(
    signal_channel: Annotated[str, Field(description="Name of loaded signal channel")],
    blank_channel: Annotated[str, Field(description="Name of loaded blank channel")],
    n_ranges: Annotated[
        int, Field(description="Number of intensity ranges (default: 5)", ge=3, le=10)
    ] = 5,
    range_method: Annotated[
        Literal["percentile", "otsu"],
        Field(description="Range boundary method (default: percentile)"),
    ] = "percentile",
    tissue_type: Annotated[
        str, Field(description="Tissue type for context (optional)")
    ] = "unknown",
    marker_name: Annotated[
        str | None, Field(description="Marker name for context (optional)")
    ] = None,
) -> dict[str, Any]:
    """
    Preview range assignments and suggested weights without performing subtraction.

    Shows per-range pixel percentages, AF ratios, and weights so the user
    can understand the plan before applying.

    Parameters:
    - signal_channel: Name of loaded signal channel
    - blank_channel: Name of loaded blank channel
    - n_ranges: Number of intensity ranges (default: 5)
    - range_method: Boundary method: "percentile" or "otsu"
    - tissue_type: Tissue type for context
    - marker_name: Marker name for context
    """
    try:
        clamp_int(n_ranges, "n_ranges", minimum=3, maximum=10)
    except ValueError as e:
        return {"error": str(e)}

    try:
        signal_data = _get_image(signal_channel)
        blank_data = _get_image(blank_channel)
    except ValueError as e:
        return {"error": str(e)}

    # Convert to numpy
    if hasattr(signal_data, "compute"):
        if signal_data.shape[0] > 2000 or signal_data.shape[1] > 2000:
            signal_np = signal_data[::2, ::2].compute()
            blank_np = blank_data[::2, ::2].compute()
        else:
            signal_np = signal_data.compute()
            blank_np = blank_data.compute()
    else:
        signal_np = np.array(signal_data)
        blank_np = np.array(blank_data)

    try:
        from kintsugi.signal.autofluorescence import analyze_for_weighted_subtraction

        analysis = analyze_for_weighted_subtraction(
            signal_np,
            blank_np,
            tissue_type=tissue_type,
            marker_name=marker_name,
            n_ranges=n_ranges,
            range_method=range_method,
        )

        return {
            "status": "success",
            "method": "weighted",
            "n_ranges": analysis["n_ranges"],
            "range_method": analysis["range_method"],
            "base_scale_factor": analysis["base_scale_factor"],
            "blank_clip_factor": analysis["blank_clip_factor"],
            "transition_width": analysis["transition_width"],
            "ranges": analysis["ranges"],
            "confidence": analysis["confidence"],
            "analysis": {
                "signal_blank_correlation": analysis["analysis"]["correlation"],
                "af_contribution": analysis["analysis"]["af_contribution"],
            },
            "tissue_type": tissue_type,
            "marker_name": marker_name,
        }
    except Exception as e:
        return {"error": f"Weighted analysis failed: {e}"}


async def subtract_blank(
    signal_channel: Annotated[str, Field(description="Name of loaded signal channel")],
    blank_channel: Annotated[str, Field(description="Name of loaded blank channel")],
    blank_clip_factor: Annotated[
        int | None,
        Field(
            description="Clip blank values below this threshold to 0; None = auto-suggest "
            "(default: 0)",
            ge=0,
        ),
    ] = None,
    blank_scale_factor: Annotated[
        float | None,
        Field(
            description="Scale factor applied to the blank before subtraction; None = "
            "auto-suggest (default: 1.0)",
            ge=0.0,
            allow_inf_nan=False,
        ),
    ] = None,
    smooth_low: Annotated[
        bool | None,
        Field(description="Apply smoothing to low-intensity regions; None = auto-suggest"),
    ] = None,
    low_size: Annotated[
        int | None,
        Field(
            description="Uniform filter size in pixels for low-intensity smoothing; None = "
            "auto-suggest (default: 2)",
            ge=1,
        ),
    ] = None,
    low_percentile: Annotated[
        int | None,
        Field(
            description="Intensity percentile (0-100) below which pixels count as low-intensity; "
            "None = auto-suggest (default: 60)",
            ge=0,
            le=100,
        ),
    ] = None,
    smooth_high: Annotated[
        bool | None,
        Field(description="Apply smoothing to high-intensity regions; None = auto-suggest"),
    ] = None,
    high_size: Annotated[
        int | None,
        Field(
            description="Uniform filter size in pixels for high-intensity smoothing; None = "
            "auto-suggest (default: 2)",
            ge=1,
        ),
    ] = None,
    high_percentile: Annotated[
        int | None,
        Field(
            description="Intensity percentile (0-100) above which pixels count as "
            "high-intensity; None = auto-suggest (default: 90)",
            ge=0,
            le=100,
        ),
    ] = None,
    erosion: Annotated[
        int | None,
        Field(
            description="Erosion disk size in pixels for void removal; None = auto-suggest "
            "(default: 0)",
            ge=0,
            le=50,
        ),
    ] = None,
    auto_suggest: Annotated[
        bool,
        Field(
            description="Auto-suggest every parameter left as None from image analysis and "
            "learned history (default: true)"
        ),
    ] = True,
    tissue_type: Annotated[
        str,
        Field(
            description="Tissue type (e.g., 'tonsil') used to look up learned parameters "
            "(default: 'unknown')"
        ),
    ] = "unknown",
    marker_name: Annotated[
        str | None,
        Field(description="Marker name (e.g., 'CD3') used to look up learned parameters"),
    ] = None,
    project_path: Annotated[
        str | None,
        Field(description="KINTSUGI project directory used for parameter learning"),
    ] = None,
    method: Annotated[
        Literal["global", "weighted"],
        Field(
            description="'global' applies one scale factor everywhere; 'weighted' applies "
            "per-intensity-range weights, which analyze_weighted_subtraction previews "
            "(default: global)"
        ),
    ] = "global",
    n_ranges: Annotated[
        int,
        Field(
            description="Number of intensity ranges for method='weighted' (default: 5)",
            ge=3,
            le=10,
        ),
    ] = 5,
    range_method: Annotated[
        Literal["percentile", "otsu"],
        Field(description="Range boundary method for method='weighted' (default: percentile)"),
    ] = "percentile",
    transition_width: Annotated[
        float,
        Field(
            description="Fractional width of the blend between adjacent intensity ranges for "
            "method='weighted' (default: 0.1)",
            ge=0.0,
            le=0.5,
            allow_inf_nan=False,
        ),
    ] = 0.1,
    ctx: Context | None = None,
) -> dict[str, Any]:
    """
    Subtract autofluorescence/blank channel from signal.

    Uses the KINTSUGI autofluorescence subtraction algorithm (ini_params).

    If auto_suggest=True or clip/scale factors are None, parameters will be
    automatically suggested based on image analysis.

    Parameters:
    - signal_channel: Name of loaded signal channel
    - blank_channel: Name of loaded blank channel
    - blank_clip_factor: Clip blank values below this (None = auto-suggest)
    - blank_scale_factor: Scale blank before subtraction (None = auto-suggest)
    - smooth_low: Smooth low-intensity regions (None = auto-suggest)
    - low_size: Filter size for low smoothing
    - low_percentile: Percentile threshold for low regions
    - smooth_high: Smooth high-intensity regions (None = auto-suggest)
    - high_size: Filter size for high smoothing
    - high_percentile: Percentile threshold for high regions
    - erosion: Erosion radius for edge cleanup (None = auto-suggest)
    - auto_suggest: Use intelligent parameter suggestion for None params
    - tissue_type: Tissue type for parameter suggestion
    - marker_name: Marker name for parameter suggestion
    - project_path: Project path for parameter learning
    - method: "global" (single scale factor) or "weighted" (per-intensity-range weights)
    - n_ranges: Number of intensity ranges for the weighted method
    - range_method: Range boundary method for the weighted method
    - transition_width: Blend width between adjacent ranges for the weighted method
    """
    # Validate untrusted arguments (the low-level server does not enforce schema).
    try:
        validate_choice(method, {"global", "weighted"}, "method")
        validate_choice(range_method, {"percentile", "otsu"}, "range_method")
        clamp_int(n_ranges, "n_ranges", minimum=3, maximum=10)
        clamp_float(transition_width, "transition_width", minimum=0.0, maximum=0.5)
        if blank_scale_factor is not None:
            clamp_float(blank_scale_factor, "blank_scale_factor", minimum=0.0)
        if erosion is not None:
            clamp_int(erosion, "erosion", minimum=0, maximum=50)
    except ValueError as e:
        return {"error": str(e)}

    # Get loaded images
    try:
        signal_data = _get_image(signal_channel)
        blank_data = _get_image(blank_channel)
    except ValueError as e:
        return {"error": str(e)}

    # ---- Weighted method path ----
    if method == "weighted":
        await _report_progress(ctx, 0, 3, "loading images")
        # Materialize on a worker thread so the event loop keeps serving pings.
        signal_np = await asyncio.to_thread(_as_numpy, signal_data)
        blank_np = await asyncio.to_thread(_as_numpy, blank_data)

        try:
            from kintsugi.signal.autofluorescence import (
                compute_intensity_ranges,
                compute_weighted_subtraction_quality,
                subtract_autofluorescence_weighted,
            )

            # Auto-suggest clip factor if needed
            if blank_clip_factor is None:
                suggestion_result = await suggest_subtraction_parameters(
                    signal_channel=signal_channel,
                    blank_channel=blank_channel,
                    tissue_type=tissue_type,
                    marker_name=marker_name,
                    project_path=project_path,
                )
                if suggestion_result.get("status") == "success":
                    sp = suggestion_result.get("suggested_parameters", {})
                    blank_clip_factor = sp.get("blank_clip_factor", 0)

            blank_clip_factor = blank_clip_factor or 0
            base_scale = blank_scale_factor if blank_scale_factor is not None else 1.0

            await _report_progress(ctx, 1, 3, "weighted subtraction")
            result_np = await asyncio.to_thread(
                subtract_autofluorescence_weighted,
                signal_np,
                blank_np,
                blank_clip_factor=blank_clip_factor,
                base_scale_factor=base_scale,
                n_ranges=n_ranges,
                range_method=range_method,
                transition_width=transition_width,
                smooth_low=smooth_low or False,
                low_size=low_size or 2,
                low_percentile=low_percentile or 60,
                smooth_high=smooth_high or False,
                high_size=high_size or 2,
                high_percentile=high_percentile or 90,
                erosion=erosion or 0,
            )
            await _report_progress(ctx, 2, 3, "scoring result")

            # Compute quality with ranges
            signal_f = signal_np.astype(np.float32)
            blank_f = blank_np.astype(np.float32)
            blank_clipped = np.clip(blank_f, blank_clip_factor, blank_f.max())
            blank_clipped[blank_clipped <= blank_clip_factor] = 0
            ranges = compute_intensity_ranges(
                signal_f, blank_clipped, n_ranges=n_ranges, method=range_method
            )
            quality = compute_weighted_subtraction_quality(signal_np, result_np, blank_np, ranges)

            # Store as dask for consistency
            try:
                import dask.array as da

                result = da.from_array(result_np, chunks=(1000, 1000))
            except ImportError:
                result = result_np

            output_name = f"{signal_channel}_subtracted"
            _store_image(
                output_name,
                result,
                {
                    "source": signal_channel,
                    "blank": blank_channel,
                    "operation": "subtract_blank_weighted",
                },
            )
            _add_history(
                output_name,
                "subtract_blank",
                {
                    "signal_channel": signal_channel,
                    "blank_channel": blank_channel,
                    "method": "weighted",
                    "blank_clip_factor": blank_clip_factor,
                    "base_scale_factor": base_scale,
                    "n_ranges": n_ranges,
                    "range_method": range_method,
                },
            )

            stats = {
                "min": float(np.min(result_np)),
                "max": float(np.max(result_np)),
                "mean": float(np.mean(result_np)),
                "std": float(np.std(result_np)),
            }

            await _report_progress(ctx, 3, 3, "done")
            return {
                "status": "success",
                "output_name": output_name,
                "method": "weighted",
                "shape": list(result_np.shape),
                "stats": stats,
                "parameters_used": {
                    "method": "weighted",
                    "blank_clip_factor": blank_clip_factor,
                    "base_scale_factor": base_scale,
                    "n_ranges": n_ranges,
                    "range_method": range_method,
                    "transition_width": transition_width,
                },
                "range_metrics": quality.get("per_range", []),
                "quality": quality.get("global", {}),
            }

        except ImportError as e:
            return {"error": f"Weighted subtraction requires kintsugi.signal: {e}"}
        except Exception as e:
            return {"error": f"Weighted subtraction failed: {e}"}

    # ---- Global method path (original) ----
    try:
        import dask.array as da
        import dask_image.ndfilters
        import dask_image.ndmorph
        from skimage import morphology
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    await _report_progress(ctx, 0, 2, "preparing subtraction")

    # Auto-suggest parameters if needed
    suggested_params = None
    if auto_suggest and (blank_clip_factor is None or blank_scale_factor is None):
        suggestion_result = await suggest_subtraction_parameters(
            signal_channel=signal_channel,
            blank_channel=blank_channel,
            tissue_type=tissue_type,
            marker_name=marker_name,
            use_learning=True,
            project_path=project_path,
        )
        if suggestion_result.get("status") == "success":
            suggested_params = suggestion_result.get("suggested_parameters", {})
            logger.info(f"Using suggested parameters: {suggested_params}")

    # Apply suggested or default parameters
    if suggested_params:
        if blank_clip_factor is None:
            blank_clip_factor = suggested_params.get("blank_clip_factor", 0)
        if blank_scale_factor is None:
            blank_scale_factor = suggested_params.get("blank_scale_factor", 1.0)
        if smooth_low is None:
            smooth_low = suggested_params.get("smooth_low", False)
        if low_size is None:
            low_size = suggested_params.get("low_size", 2)
        if low_percentile is None:
            low_percentile = suggested_params.get("low_percentile", 60)
        if smooth_high is None:
            smooth_high = suggested_params.get("smooth_high", False)
        if high_size is None:
            high_size = suggested_params.get("high_size", 2)
        if high_percentile is None:
            high_percentile = suggested_params.get("high_percentile", 90)
        if erosion is None:
            erosion = suggested_params.get("erosion", 0)
    else:
        # Use defaults if no suggestion available
        blank_clip_factor = blank_clip_factor if blank_clip_factor is not None else 0
        blank_scale_factor = blank_scale_factor if blank_scale_factor is not None else 1.0
        smooth_low = smooth_low if smooth_low is not None else False
        low_size = low_size if low_size is not None else 2
        low_percentile = low_percentile if low_percentile is not None else 60
        smooth_high = smooth_high if smooth_high is not None else False
        high_size = high_size if high_size is not None else 2
        high_percentile = high_percentile if high_percentile is not None else 90
        erosion = erosion if erosion is not None else 0

    # Ensure data are dask arrays for downstream dask operations
    if not isinstance(blank_data, da.Array):
        blank_data = da.from_array(np.asarray(blank_data), chunks=(1000, 1000))
    if not isinstance(signal_data, da.Array):
        signal_data = da.from_array(np.asarray(signal_data), chunks=(1000, 1000))

    # Ensure compatible types
    blank_copy = blank_data.copy()
    signal_copy = signal_data.astype(blank_data.dtype)

    # Apply blank subtraction (from Kutils.ini_params)
    blank_copy = da.clip(blank_copy, blank_clip_factor, blank_copy.max())
    blank_copy = da.where(blank_copy <= blank_clip_factor, 0, blank_copy)
    result = signal_copy - da.minimum(signal_copy, blank_copy * blank_scale_factor)

    # Low intensity smoothing
    if smooth_low:
        low_mask = da.where(signal_copy < da.percentile(signal_copy.ravel(), low_percentile), 1, 0)
        low_mask = dask_image.ndmorph.binary_dilation(low_mask, morphology.disk(1))
        result = da.where(
            low_mask,
            dask_image.ndfilters.uniform_filter(result, size=low_size),
            result,
        )

    # High intensity smoothing
    if smooth_high:
        high_mask = da.where(
            signal_copy > da.percentile(signal_copy.ravel(), high_percentile), 1, 0
        )
        high_mask = dask_image.ndmorph.binary_dilation(high_mask, morphology.disk(1))
        result = da.where(
            high_mask,
            dask_image.ndfilters.uniform_filter(result, size=high_size),
            result,
        )

    # Erosion for void removal
    if erosion > 0:
        erod_mask = da.where(result > 0, 1, 0)
        signal_mask = dask_image.ndmorph.binary_erosion(erod_mask, morphology.disk(erosion))
        result = da.where(signal_mask, result, 0)

    result = result.astype(np.uint16)

    # Store result with new name
    output_name = f"{signal_channel}_subtracted"
    _store_image(
        output_name,
        result,
        {
            "source": signal_channel,
            "blank": blank_channel,
            "operation": "subtract_blank",
        },
    )

    # Record history
    _add_history(
        output_name,
        "subtract_blank",
        {
            "signal_channel": signal_channel,
            "blank_channel": blank_channel,
            "method": "global",
            "blank_clip_factor": blank_clip_factor,
            "blank_scale_factor": blank_scale_factor,
            "smooth_low": smooth_low,
            "smooth_high": smooth_high,
            "erosion": erosion,
        },
    )

    # Compute stats for result (this evaluates the dask graph, so keep it off
    # the event loop).
    try:
        sample = await asyncio.to_thread(result[::10, ::10].compute)
        stats = {
            "min": float(np.min(sample)),
            "max": float(np.max(sample)),
            "mean": float(np.mean(sample)),
            "std": float(np.std(sample)),
        }
    except Exception:
        stats = {}
    await _report_progress(ctx, 1, 2, "subtraction computed")

    await _report_progress(ctx, 2, 2, "done")
    return {
        "status": "success",
        "output_name": output_name,
        "method": "global",
        "shape": list(result.shape),
        "stats": stats,
        "parameters_used": {
            "method": "global",
            "blank_clip_factor": blank_clip_factor,
            "blank_scale_factor": blank_scale_factor,
            "smooth_low": smooth_low,
            "smooth_high": smooth_high,
            "erosion": erosion,
        },
    }


async def denoise(
    channel: Annotated[str, Field(description="Name of loaded channel to denoise")],
    method: Annotated[
        Literal["percentile", "uniform", "median"],
        Field(description="Denoising method (default: median)"),
    ] = "median",
    filter_size: Annotated[
        int, Field(description="Filter kernel size in pixels (default: 3)", ge=1, le=99)
    ] = 3,
    percentile: Annotated[
        int,
        Field(description="Percentile value for percentile filter (default: 10)", ge=0, le=100),
    ] = 10,
    upper: Annotated[
        bool,
        Field(
            description="Use the upper percentile band instead of the lower one; only affects "
            "method='percentile' (default: false)"
        ),
    ] = False,
) -> dict[str, Any]:
    """
    Apply denoising filters to reduce noise.

    Methods:
    - percentile: Percentile filter (good for salt-and-pepper noise)
    - uniform: Uniform/box filter (smoothing)
    - median: Median filter (edge-preserving)
    """
    try:
        clamp_int(filter_size, "filter_size", minimum=1, maximum=99)
        clamp_int(percentile, "percentile", minimum=0, maximum=100)
    except ValueError as e:
        return {"error": str(e)}

    try:
        import dask.array as da
        import dask_image.ndfilters
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Ensure data is a dask array for downstream dask operations
    if not isinstance(data, da.Array):
        data = da.from_array(np.asarray(data), chunks=(1000, 1000))

    data_copy = data.copy()

    if method == "percentile":
        p = -percentile if upper else percentile
        result = dask_image.ndfilters.percentile_filter(data_copy, p, filter_size)
    elif method == "uniform":
        result = dask_image.ndfilters.uniform_filter(data_copy, filter_size)
    elif method == "median":
        result = dask_image.ndfilters.median_filter(data_copy, filter_size)
    else:
        return {"error": f"Unknown method: {method}"}

    # Store result
    output_name = f"{channel}_denoised"
    _store_image(
        output_name,
        result,
        {
            "source": channel,
            "operation": "denoise",
            "method": method,
        },
    )

    _add_history(
        output_name,
        "denoise",
        {
            "method": method,
            "filter_size": filter_size,
            "percentile": percentile if method == "percentile" else None,
        },
    )

    return {
        "status": "success",
        "output_name": output_name,
        "method": method,
        "filter_size": filter_size,
    }


async def apply_clahe(
    channel: Annotated[str, Field(description="Name of loaded channel")],
    clip_limit: Annotated[
        float,
        Field(
            description="Clip limit for contrast limiting, as a fraction of the histogram "
            "(default: 0.01)",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ] = 0.01,
    tile_grid_size: Annotated[
        int,
        Field(
            description="Size of grid for histogram equalization; the kernel is image height "
            "divided by this value (default: 70)",
            ge=1,
            le=10000,
        ),
    ] = 70,
    nbins: Annotated[
        int, Field(description="Number of histogram bins (default: 128)", ge=2, le=65536)
    ] = 128,
) -> dict[str, Any]:
    """
    Apply Contrast Limited Adaptive Histogram Equalization.

    CLAHE enhances local contrast while limiting noise amplification.
    """
    try:
        clamp_float(clip_limit, "clip_limit", minimum=0.0, maximum=1.0)
        clamp_int(tile_grid_size, "tile_grid_size", minimum=1, maximum=10000)
        clamp_int(nbins, "nbins", minimum=2, maximum=65536)
    except ValueError as e:
        return {"error": str(e)}

    try:
        from skimage.exposure import equalize_adapthist, rescale_intensity
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute if dask array
    if hasattr(data, "compute"):
        data = data.compute()

    data_uint16 = data.astype(np.uint16)
    kernel_size = data_uint16.shape[0] // tile_grid_size

    # Apply CLAHE
    result = equalize_adapthist(
        data_uint16,
        clip_limit=clip_limit,
        kernel_size=kernel_size,
        nbins=nbins,
    )

    # Rescale back to original range
    result = rescale_intensity(
        result,
        out_range=(float(data_uint16.min()), float(data_uint16.max())),
    ).astype(np.float64)

    # Store result
    output_name = f"{channel}_clahe"
    _store_image(
        output_name,
        result,
        {
            "source": channel,
            "operation": "clahe",
        },
    )

    _add_history(
        output_name,
        "clahe",
        {
            "clip_limit": clip_limit,
            "tile_grid_size": tile_grid_size,
            "nbins": nbins,
            "actual_kernel_size": kernel_size,
        },
    )

    return {
        "status": "success",
        "output_name": output_name,
        "tile_size": f"{kernel_size}x{kernel_size}",
        "clip_limit": clip_limit,
    }


async def clean_background(
    channel: Annotated[str, Field(description="Name of loaded channel")],
    background_threshold: Annotated[
        int,
        Field(
            description="Intensity threshold below which pixels are set to 0 (default: 100)",
            ge=0,
        ),
    ] = 100,
    smooth: Annotated[
        bool,
        Field(
            description="Apply a 3x3 median filter to pixels below smooth_threshold to soften "
            "transition zones after thresholding (default: false)"
        ),
    ] = False,
    smooth_threshold: Annotated[
        int,
        Field(
            description="Intensity below which pixels are median-smoothed when smooth=True "
            "(default: 100)",
            ge=0,
        ),
    ] = 100,
    remove_small_objects: Annotated[
        bool, Field(description="Whether to remove small connected components")
    ] = False,
    min_object_size: Annotated[
        int, Field(description="Minimum object size in pixels to keep (default: 30)", ge=1)
    ] = 30,
) -> dict[str, Any]:
    """
    Remove background and small objects from the image.
    """
    try:
        import dask.array as da
        import dask_image.ndfilters
        from skimage import morphology
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Ensure data is a dask array for downstream dask operations
    if not isinstance(data, da.Array):
        data = da.from_array(np.asarray(data), chunks=(1000, 1000))

    result = data.copy()

    # Background thresholding
    result = da.where(result <= background_threshold, 0, result)

    # Smooth transition zones
    if smooth:
        transition_mask = da.where((result < smooth_threshold), result, 0)
        smoothed = dask_image.ndfilters.median_filter(result, size=3)
        result = da.where(transition_mask, smoothed, result)

    # Remove small objects
    if remove_small_objects:
        binary_mask = result > 0

        def remove_small_chunk(chunk):
            return morphology.remove_small_objects(chunk, min_size=min_object_size, connectivity=2)

        filtered_mask = binary_mask.map_blocks(remove_small_chunk, dtype=bool)
        result = da.where(filtered_mask, result, 0)

        def apply_closing(chunk):
            return morphology.closing(chunk, morphology.disk(1))

        result = result.map_blocks(apply_closing, dtype=result.dtype)

    # Store result
    output_name = f"{channel}_cleaned"
    _store_image(
        output_name,
        result,
        {
            "source": channel,
            "operation": "clean_background",
        },
    )

    _add_history(
        output_name,
        "clean_background",
        {
            "background_threshold": background_threshold,
            "smooth": smooth,
            "remove_small_objects": remove_small_objects,
            "min_object_size": min_object_size,
        },
    )

    return {
        "status": "success",
        "output_name": output_name,
        "background_threshold": background_threshold,
        "small_objects_removed": remove_small_objects,
    }


async def gaussian_subtract(
    channel: Annotated[str, Field(description="Name of loaded channel")],
    sigma: Annotated[
        int, Field(description="Gaussian sigma value in pixels (default: 10)", ge=1, le=1000)
    ] = 10,
    scale_factor: Annotated[
        float,
        Field(
            description="Scale factor for subtraction (default: 0.1)",
            ge=0.0,
            allow_inf_nan=False,
        ),
    ] = 0.1,
    low_clip: Annotated[
        int,
        Field(
            description="Lower intensity clip bound on the Gaussian-blurred background: blurred "
            "values at or below it are zeroed before subtraction (16-bit range, default: 0)",
            ge=0,
        ),
    ] = 0,
    high_clip: Annotated[
        int,
        Field(
            description="Upper intensity clip bound on the Gaussian-blurred background: blurred "
            "values at or above it are zeroed before subtraction (16-bit range, default: 65535)",
            ge=0,
        ),
    ] = 65535,
) -> dict[str, Any]:
    """
    Subtract Gaussian-blurred version to remove structured background.
    """
    try:
        clamp_int(sigma, "sigma", minimum=1, maximum=1000)
    except ValueError as e:
        return {"error": str(e)}

    try:
        import dask.array as da
        import dask_image.ndfilters
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Ensure data is a dask array for downstream dask operations
    if not isinstance(data, da.Array):
        data = da.from_array(np.asarray(data), chunks=(1000, 1000))

    # Create Gaussian-blurred version
    blurred = dask_image.ndfilters.gaussian_filter(data, sigma=sigma)

    # Clip extreme values
    blurred = da.where(low_clip >= blurred, 0, blurred)
    blurred = da.where(blurred >= high_clip, 0, blurred)

    # Subtract scaled blurred from original
    result = data - da.minimum(data, blurred * scale_factor)

    # Store result
    output_name = f"{channel}_gaussub"
    _store_image(
        output_name,
        result,
        {
            "source": channel,
            "operation": "gaussian_subtract",
        },
    )

    _add_history(
        output_name,
        "gaussian_subtract",
        {
            "sigma": sigma,
            "scale_factor": scale_factor,
            "low_clip": low_clip,
            "high_clip": high_clip,
        },
    )

    return {
        "status": "success",
        "output_name": output_name,
        "sigma": sigma,
        "scale_factor": scale_factor,
    }


async def denoise_advanced(
    channel: Annotated[str, Field(description="Name of loaded channel to denoise")],
    method: Annotated[
        Literal["adaptive", "bilateral", "nlm", "bm3d", "n2v", "patch_svd"],
        Field(description="Advanced denoising method (default: adaptive)"),
    ] = "adaptive",
    strength: Annotated[
        Literal["light", "medium", "strong", "auto"],
        Field(description="Denoising strength (default: auto)"),
    ] = "auto",
    preserve_edges: Annotated[
        bool, Field(description="Prioritize edge preservation (default: true)")
    ] = True,
    n2v_epochs: Annotated[
        int, Field(description="Training epochs for N2V method (default: 50)", ge=1, le=10000)
    ] = 50,
    model_path: Annotated[
        str | None,
        Field(
            description="Path to a pre-trained N2V checkpoint (.pt/.pth); skips training. Must "
            "lie inside project_path, the Claude Code project, or ~/.kintsugi/models"
        ),
    ] = None,
    project_path: Annotated[
        str | None,
        Field(description="KINTSUGI project directory that model_path may be loaded from"),
    ] = None,
    ctx: Context | None = None,
) -> dict[str, Any]:
    """
    Apply advanced denoising using N2V, CARE, or adaptive methods.

    Methods:
    - adaptive: Automatically select best traditional filter based on noise level
    - bilateral: Edge-preserving bilateral filter
    - nlm: Non-local means denoising
    - bm3d: BM3D-inspired block matching denoising
    - n2v: Self-supervised Noise2Void (requires PyTorch, slower but very effective)
    - patch_svd: SVD-based patch denoising

    Parameters:
    - channel: Name of loaded channel to denoise
    - method: Denoising algorithm to use
    - strength: "light", "medium", "strong", or "auto"
    - preserve_edges: Prioritize edge preservation (for bilateral/adaptive)
    - n2v_epochs: Training epochs for N2V (only used with n2v method)
    - model_path: Path to pre-trained N2V model (skip training if provided)

    Returns:
    - status, output name, noise estimates, and method details
    """
    try:
        clamp_int(n2v_epochs, "n2v_epochs", minimum=1, maximum=10000)
        if model_path is not None:
            # Loading a checkpoint deserializes a file the model chose: confine it
            # to directories the user controls before torch sees it.
            resolved_model = ensure_within_any(
                model_path, allowed_model_roots(project_path), "model_path"
            )
            if resolved_model.suffix.lower() not in {".pt", ".pth"}:
                return {"error": f"model_path must be a .pt/.pth checkpoint: {model_path}"}
            if not resolved_model.is_file():
                return {"error": f"model_path does not exist or is not a file: {model_path}"}
            model_path = str(resolved_model)
    except ValueError as e:  # PathSafetyError is a ValueError
        return {"error": str(e)}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    await _report_progress(ctx, 0, 3, "loading image")

    # Convert to numpy if needed (off the event loop: compute() can be heavy)
    if hasattr(data, "compute") and (data.shape[0] > 4000 or data.shape[1] > 4000):
        # For large images, process a reasonable subset or use tiled approach
        logger.warning("Large image detected, computing full array...")
    img = await asyncio.to_thread(_as_numpy, data)

    # Estimate noise level before denoising
    try:
        from kintsugi.denoise.filters import estimate_noise_level

        noise_before = estimate_noise_level(img, method="mad")
    except ImportError:
        noise_before = float(np.std(img[img < np.percentile(img, 20)]))

    # Pick the denoiser, then run it on a worker thread so the server keeps
    # answering pings/cancellation while it works.
    await _report_progress(ctx, 1, 3, f"denoising ({method})")
    try:
        denoiser: Callable[..., np.ndarray]
        kwargs: dict[str, Any]
        if method == "adaptive":
            from kintsugi.denoise.filters import adaptive_denoise

            denoiser = adaptive_denoise
            kwargs = {"strength": strength, "preserve_edges": preserve_edges}

        elif method == "bilateral":
            from kintsugi.denoise.filters import denoise_bilateral

            sigma_spatial = 1.0 if strength == "light" else (2.0 if strength == "medium" else 3.0)
            denoiser = denoise_bilateral
            kwargs = {"sigma_spatial": sigma_spatial}

        elif method == "nlm":
            from kintsugi.denoise.filters import denoise_nlm

            patch_size = 5 if strength in ("light", "auto") else 7
            denoiser = denoise_nlm
            kwargs = {"patch_size": patch_size}

        elif method == "bm3d":
            from kintsugi.denoise.patch_based import denoise_bm3d_lite

            denoiser = denoise_bm3d_lite
            kwargs = {"sigma": noise_before if strength == "auto" else None}

        elif method == "n2v":
            from kintsugi.denoise.n2v import denoise_n2v

            denoiser = denoise_n2v
            kwargs = {"n_epochs": n2v_epochs, "model_path": model_path}

        elif method == "patch_svd":
            from kintsugi.denoise.patch_based import denoise_svd_patch

            rank = None if strength == "auto" else (5 if strength == "light" else 10)
            denoiser = denoise_svd_patch
            kwargs = {"rank": rank}

        else:
            return {"error": f"Unknown advanced denoise method: {method}"}

        result = await asyncio.to_thread(denoiser, img, **kwargs)

    except ImportError as e:
        return {
            "error": f"Required module not available: {e}",
            "hint": "Install with: pip install kintsugi[denoise] for N2V/CARE support",
        }

    await _report_progress(ctx, 2, 3, "estimating residual noise")

    # Estimate noise after
    try:
        noise_after = estimate_noise_level(result, method="mad")
    except Exception:
        noise_after = float(np.std(result[result < np.percentile(result, 20)]))

    # Store result
    output_name = f"{channel}_{method}_denoised"
    import dask.array as da

    result_dask = da.from_array(result, chunks=(1000, 1000))
    _store_image(
        output_name,
        result_dask,
        {
            "source": channel,
            "operation": f"denoise_advanced_{method}",
        },
    )

    _add_history(
        output_name,
        "denoise_advanced",
        {
            "method": method,
            "strength": strength,
            "preserve_edges": preserve_edges,
            "noise_before": noise_before,
            "noise_after": noise_after,
        },
    )

    await _report_progress(ctx, 3, 3, "done")
    return {
        "status": "success",
        "output_name": output_name,
        "method": method,
        "noise_reduction": {
            "before": round(noise_before, 2),
            "after": round(noise_after, 2),
            "reduction_pct": (
                round((1 - noise_after / noise_before) * 100, 1) if noise_before > 0 else 0
            ),
        },
        "parameters": {
            "strength": strength,
            "preserve_edges": preserve_edges,
        },
    }


# Utility functions for state management
def get_loaded_images() -> dict[str, Any]:
    """Get all loaded images and their metadata."""
    return {
        name: {
            "shape": info["shape"],
            "dtype": info["dtype"],
            "metadata": info["metadata"],
        }
        for name, info in _loaded_images.items()
    }


def get_processing_history(channel: str) -> list:
    """Get processing history for a channel."""
    return _processing_history.get(channel, [])


def clear_session():
    """Clear all loaded images and history."""
    _loaded_images.clear()
    _processing_history.clear()


async def cluster_channels_tool(
    project_path: Annotated[str, Field(description="Path to KINTSUGI project directory")],
    n_clusters: Annotated[
        int | None, Field(description="Number of clusters; None = auto-select", ge=2)
    ] = None,
    wavelength_aware: Annotated[
        bool, Field(description="Enforce same-wavelength clustering (default: true)")
    ] = True,
    ctx: Context | None = None,
) -> dict[str, Any]:
    """
    Cluster project channels by image feature similarity.

    Loads registered/EDF images, extracts features, clusters, and returns
    cluster assignments with representative channels.
    """
    try:
        from kintsugi.signal.clustering import cluster_channels, get_cluster_representatives
        from kintsugi.signal.features import batch_extract_features
    except ImportError as e:
        return {"error": f"Missing module: {e}"}

    project_path_obj = Path(project_path).resolve()

    await _report_progress(ctx, 0, 3, "discovering channel images")

    # Discover channel images
    search_dirs = [
        project_path_obj / "data" / "processed" / "registered",
        project_path_obj / "data" / "processed" / "edf",
        project_path_obj / "data" / "processed" / "stitched",
    ]

    marker_dict: dict[str, Any] = {}
    path_by_name: dict[str, Path] = {}
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        for tif in sorted(search_dir.rglob("*.tif")):
            name = tif.stem
            if name in marker_dict:
                continue
            try:
                # rglob follows symlinks; keep every file inside the project.
                tif = ensure_within(tif, project_path_obj, "channel file")
            except PathSafetyError as e:
                logger.warning(str(e))
                continue
            try:
                import tifffile

                marker_dict[name] = tifffile.imread(str(tif))
                path_by_name[name] = tif
            except Exception as e:
                logger.warning(f"Failed to load {tif}: {e}")
        if marker_dict:
            break

    if not marker_dict:
        return {"error": "No channel images found in project"}

    # Extract features (heavy, CPU-bound: run on a worker thread)
    await _report_progress(ctx, 1, 3, f"extracting features from {len(marker_dict)} channels")
    try:
        features = await asyncio.to_thread(batch_extract_features, marker_dict, progress=False)
    except Exception as e:
        return {"error": f"Feature extraction failed: {e}"}

    # Cluster
    await _report_progress(ctx, 2, 3, "clustering channels")
    try:
        assignments = cluster_channels(
            features,
            n_clusters=n_clusters,
            wavelength_aware=wavelength_aware,
        )
        representatives = get_cluster_representatives(features, assignments)
    except Exception as e:
        return {"error": f"Clustering failed: {e}"}

    # Keep every clustered channel in the session, tagged with its cluster, so
    # propagate_parameters can find the members (it filters on metadata.cluster_id).
    for name, cid in assignments.items():
        _store_image(
            name,
            marker_dict[name],
            {
                "cluster_id": int(cid),
                "source_path": str(path_by_name[name]),
                "project": str(project_path_obj),
                "wavelength_group": features[name].get("wavelength_group"),
            },
        )

    # Format result
    cluster_info = {}
    for cid, (rep_name, count) in representatives.items():
        members = [name for name, c in assignments.items() if c == cid]
        cluster_info[str(cid)] = {
            "representative": rep_name,
            "member_count": count,
            "members": members,
            "wavelength_group": features[rep_name]["wavelength_group"],
        }

    await _report_progress(ctx, 3, 3, "done")
    return {
        "status": "success",
        "n_channels": len(assignments),
        "n_clusters": len(representatives),
        "assignments": {name: int(cid) for name, cid in assignments.items()},
        "clusters": cluster_info,
        "loaded_channels": sorted(assignments),
    }


async def propagate_parameters_tool(
    project_path: Annotated[str, Field(description="Path to KINTSUGI project directory")],
    cluster_id: Annotated[int, Field(description="Cluster ID to propagate parameters to", ge=0)],
    params: Annotated[
        dict[str, Any],
        Field(description="Processing parameters (blank_params and/or clean_params dicts)"),
    ],
    output_dir: Annotated[
        str | None,
        Field(
            description="Output directory, which must lie inside the project; None = the "
            "project's data/processed/signal_isolated directory"
        ),
    ] = None,
    ctx: Context | None = None,
) -> dict[str, Any]:
    """
    Propagate tuned parameters to all members of a cluster.
    """
    try:
        from kintsugi.signal.clustering import propagate_cluster_parameters
    except ImportError as e:
        return {"error": f"Missing module: {e}"}

    project_path_obj = Path(project_path).resolve()

    try:
        if output_dir:
            out = ensure_within(output_dir, project_path_obj, "output_dir")
        else:
            out = project_path_obj / "data" / "processed" / "signal_isolated"
    except PathSafetyError as e:
        return {"error": str(e)}

    # Get member channels from loaded images
    loaded = get_loaded_images()
    if not loaded:
        return {"error": "No images loaded. Use cluster_channels first to load images."}

    # Get member channels for this cluster
    members = [
        name
        for name, info in loaded.items()
        if info.get("metadata", {}).get("cluster_id") == cluster_id
    ]

    if not members:
        return {
            "error": f"No channels found for cluster {cluster_id}. Load and cluster channels first."
        }

    marker_dict_local: dict[str, Any] = {name: _loaded_images[name]["data"] for name in members}
    blank_map = {
        name: info["metadata"].get("blank_name", "")
        for name, info in _loaded_images.items()
        if name in members
    }

    await _report_progress(ctx, 0, 2, f"propagating parameters to {len(members)} channels")
    try:
        results = await asyncio.to_thread(
            propagate_cluster_parameters, marker_dict_local, members, params, blank_map, out
        )
    except Exception as e:
        return {"error": f"Propagation failed: {e}"}
    await _report_progress(ctx, 1, 2, "propagation complete")

    summary = {
        name: {
            "status": r["status"],
            "quality_score": r.get("quality_score", 0),
            "output_path": str(r["output_path"]) if r.get("output_path") else None,
        }
        for name, r in results.items()
    }

    successes = sum(1 for r in results.values() if r["status"] == "success")
    await _report_progress(ctx, 2, 2, "done")
    return {
        "status": "success",
        "cluster_id": cluster_id,
        "processed": successes,
        "total": len(members),
        "results": summary,
    }


async def optimize_parameters(
    signal_channel: Annotated[str, Field(description="Name of loaded signal channel")],
    blank_channel: Annotated[str, Field(description="Name of loaded blank channel")],
    n_trials: Annotated[
        int, Field(description="Maximum optimization trials (default: 80)", ge=1, le=1000)
    ] = 80,
    timeout: Annotated[
        int,
        Field(description="Maximum seconds for optimization (default: 300)", ge=1, le=86400),
    ] = 300,
    optimize_clean: Annotated[
        bool,
        Field(description="Also optimize background cleaning parameters (default: true)"),
    ] = True,
    warm_start_params: Annotated[
        dict[str, Any] | None,
        Field(description="Initial parameter guess for faster convergence (optional)"),
    ] = None,
    project_path: Annotated[
        str | None,
        Field(description="Project path for storing optimization study (optional)"),
    ] = None,
    ctx: Context | None = None,
) -> dict[str, Any]:
    """
    Find optimal signal isolation parameters via Bayesian optimization (Optuna).

    Automatically searches the parameter space using a quality metric to
    evaluate each trial. Much faster than manual slider exploration.

    Parameters:
    - signal_channel: Name of loaded signal channel
    - blank_channel: Name of loaded blank channel
    - n_trials: Maximum optimization trials (default 80)
    - timeout: Maximum seconds (default 300)
    - optimize_clean: Also optimize background cleaning parameters
    - warm_start_params: Initial parameter guess for faster convergence
    - project_path: Path to project for storing optimization study

    Returns:
    - Best parameters, quality score, and optimization summary
    """
    try:
        clamp_int(n_trials, "n_trials", minimum=1, maximum=1000)
        clamp_int(timeout, "timeout", minimum=1, maximum=86400)
    except ValueError as e:
        return {"error": str(e)}

    try:
        from kintsugi.signal.optimizer import optimize_signal_isolation
    except ImportError as e:
        return {"error": f"Missing dependency: {e}. Install with: pip install kintsugi[optimize]"}

    try:
        signal_data = _get_image(signal_channel)
        blank_data = _get_image(blank_channel)
    except ValueError as e:
        return {"error": str(e)}

    await _report_progress(ctx, 0, n_trials, f"trial 0/{n_trials}")

    # Convert to numpy (off the event loop)
    signal_np = await asyncio.to_thread(_as_numpy, signal_data)
    blank_np = await asyncio.to_thread(_as_numpy, blank_data)

    # Storage path
    storage_path = None
    if project_path:
        storage_dir = Path(project_path) / "configs" / "optuna"
        storage_dir.mkdir(parents=True, exist_ok=True)
        storage_path = storage_dir / f"{signal_channel}_optimization.db"

    # Optuna runs on a worker thread; each finished trial schedules a progress
    # notification back onto the event loop (never awaited from the thread).
    report = _thread_progress_reporter(ctx, asyncio.get_running_loop())

    def _on_trial(completed: int, total: int) -> None:
        report(completed, total, f"trial {completed}/{total}")

    try:
        result = await asyncio.to_thread(
            optimize_signal_isolation,
            signal=signal_np,
            blank=blank_np,
            n_trials=n_trials,
            timeout=timeout,
            warm_start_params=warm_start_params,
            optimize_clean=optimize_clean,
            storage_path=storage_path,
            show_progress=False,
            trial_callback=_on_trial if ctx is not None else None,
        )
    except Exception as e:
        return {"error": f"Optimization failed: {e}"}

    await _report_progress(
        ctx,
        n_trials,
        n_trials,
        f"optimization complete ({result['n_trials_completed']} trials)",
    )

    _add_history(
        signal_channel,
        "optimize_parameters",
        {
            "n_trials": result["n_trials_completed"],
            "best_quality": result["best_quality"],
            "time_seconds": result["optimization_time_seconds"],
        },
    )

    return result


async def predict_parameters(
    signal_channel: Annotated[str, Field(description="Name of loaded signal channel")],
    blank_channel: Annotated[
        str | None, Field(description="Name of loaded blank channel (optional)")
    ] = None,
    project_path: Annotated[
        str | None, Field(description="Project path containing predictor model (optional)")
    ] = None,
    model_path: Annotated[
        str | None, Field(description="Path to a trained predictor .joblib (optional)")
    ] = None,
) -> dict[str, Any]:
    """
    Predict signal isolation parameters from image features using a trained model.

    Uses a Random Forest model trained on previous tuning examples to predict
    optimal parameters without manual exploration. Includes confidence scores
    and per-parameter uncertainty estimates.

    Parameters:
    - signal_channel: Name of loaded signal channel
    - blank_channel: Name of loaded blank channel (optional, improves prediction)
    - project_path: Path to project containing predictor model
    - model_path: Explicit path to predictor model (.joblib)

    Returns:
    - Predicted parameters with confidence and uncertainty
    """
    try:
        from kintsugi.signal.features import extract_channel_features, features_to_vector
        from kintsugi.signal.predictor import ParameterPredictor
    except ImportError as e:
        return {"error": f"Missing module: {e}"}

    try:
        signal_data = _get_image(signal_channel)
    except ValueError as e:
        return {"error": str(e)}

    # Convert to numpy
    if hasattr(signal_data, "compute"):
        signal_np = signal_data.compute()
    else:
        signal_np = np.array(signal_data)

    blank_np = None
    if blank_channel:
        try:
            blank_data = _get_image(blank_channel)
            blank_np = (
                blank_data.compute() if hasattr(blank_data, "compute") else np.array(blank_data)
            )
        except ValueError:
            pass

    # Extract features
    features = extract_channel_features(
        signal_np,
        blank=blank_np,
        channel_name=signal_channel,
    )
    fv = features_to_vector(features)

    # Find model path. joblib.load unpickles arbitrary objects, so an explicit
    # path must lie in a directory the user controls.
    predictor_path = None
    if model_path:
        try:
            predictor_path = ensure_within_any(
                model_path, allowed_model_roots(project_path), "model_path"
            )
        except PathSafetyError as e:
            return {"error": str(e)}
        if predictor_path.suffix.lower() != ".joblib":
            return {"error": f"model_path must be a .joblib predictor model: {model_path}"}
    elif project_path:
        predictor_path = Path(project_path) / "configs" / "parameter_predictor.joblib"

    if predictor_path is None or not predictor_path.exists():
        return {
            "status": "no_model",
            "message": "No trained predictor model found. Train one with train_from_sqlite().",
            "features": {
                k: float(v) if isinstance(v, (int, float, np.floating)) else v
                for k, v in features.items()
            },
        }

    predictor = ParameterPredictor(model_path=predictor_path)
    result = predictor.predict(fv)

    _add_history(
        signal_channel,
        "predict_parameters",
        {
            "confidence": result.get("confidence", 0),
            "status": result["status"],
        },
    )

    return result


async def estimate_background(
    channel: Annotated[str, Field(description="Name of loaded channel image")],
    kernel_size: Annotated[
        int,
        Field(
            description="SMO kernel size in pixels (default: 7, robust across 5-15)",
            ge=1,
            le=999,
        ),
    ] = 7,
    return_background: Annotated[
        bool,
        Field(description="Store estimated background as a loaded image (default: false)"),
    ] = False,
) -> dict[str, Any]:
    """
    Estimate and subtract background using SMO (Silver Mountain Operator).

    Parameter-free background estimation — no blank channel needed.
    Useful when blank channels are unavailable or unreliable.

    Parameters:
    - channel: Name of loaded channel image
    - kernel_size: SMO kernel size (default 7, robust across 5-15)
    - return_background: If True, store estimated background as a loaded image

    Returns:
    - Background statistics, corrected image stored as '{channel}_smo_corrected'
    """
    try:
        clamp_int(kernel_size, "kernel_size", minimum=1, maximum=999)
    except ValueError as e:
        return {"error": str(e)}

    try:
        from kintsugi.signal.smo import estimate_background_smo
    except ImportError as e:
        return {"error": f"Missing dependency: {e}. Install with: pip install kintsugi[optimize]"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    if hasattr(data, "compute"):
        img_np = data.compute()
    else:
        img_np = np.array(data)

    # Handle 3D stacks by using max projection
    if img_np.ndim > 2:
        img_np = img_np.max(axis=0) if img_np.ndim == 3 else img_np[0]

    try:
        result = estimate_background_smo(
            img_np,
            kernel_size=kernel_size,
            return_background=return_background,
        )
    except Exception as e:
        return {"error": f"SMO estimation failed: {e}"}

    # Store corrected image
    corrected_name = f"{channel}_smo_corrected"
    _store_image(
        corrected_name,
        result["corrected"],
        {"source_channel": channel, "method": "smo", "kernel_size": kernel_size},
    )

    _add_history(
        channel,
        "estimate_background_smo",
        {"kernel_size": kernel_size, "background_mean": result["background_mean"]},
    )

    # Store background if requested
    if return_background and result["background"] is not None:
        bg_name = f"{channel}_smo_background"
        _store_image(
            bg_name, result["background"], {"source_channel": channel, "type": "background"}
        )

    return {
        "status": "success",
        "corrected_name": corrected_name,
        "background_mean": result["background_mean"],
        "background_std": result["background_std"],
        "kernel_size_used": result["kernel_size_used"],
        "corrected_shape": list(result["corrected"].shape),
        "corrected_dtype": str(result["corrected"].dtype),
    }
