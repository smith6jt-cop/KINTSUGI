"""
Workflow Tools for KINTSUGI MCP Server

State management, channel discovery, saving results, and Jupyter integration.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Annotated, Any, Literal

import numpy as np
from pydantic import Field

from kintsugi.mcp.heuristics import (
    HEURISTIC_OPERATIONS,
    heuristic_suggestion,
    sample_image,
    sample_step,
    summarize_intensity,
)
from kintsugi.mcp.path_safety import (
    PathSafetyError,
    ensure_within,
    safe_filename,
    validate_choice,
)
from kintsugi.mcp.schemas import SessionStateResult
from kintsugi.mcp.session import (
    NO_PROJECT_MESSAGE,
    PROJECT_PATH_HELP,
    ProjectNotConfigured,
    _loaded_images,
    _processing_history,
    clear_session,
    configured_project_root,
    estimate_memory_usage,
    resolve_project_path,
)

logger = logging.getLogger("kintsugi.mcp.workflow")


async def list_channels(
    project_path: Annotated[str | None, Field(description=PROJECT_PATH_HELP)] = None,
    cycle: Annotated[
        str | None,
        Field(
            description="Cycle directory name such as 'cyc001', or a bare number; None lists "
            "the channels of every cycle"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    List all available channels in the current project.

    Searches raw and processed directories for image files.
    """
    try:
        project_path = Path(resolve_project_path(project_path))
    except ProjectNotConfigured as e:
        return {"error": str(e)}

    # Try to load project configuration
    try:
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_path)
        raw_dir = project.paths.raw
        processed_dir = project.paths.processed
        registered_dir = project.paths.registered
    except Exception:
        raw_dir = project_path / "data" / "raw"
        processed_dir = project_path / "data" / "processed"
        registered_dir = processed_dir / "registered"

    # A cycle is a bare directory name, never a path (rejects "../..").
    if cycle:
        try:
            safe_filename(cycle, "cycle")
        except PathSafetyError as e:
            return {"error": str(e)}

    channels = {}

    # Search directories
    search_dirs = [
        ("raw", raw_dir),
        ("registered", registered_dir),
        ("processed", processed_dir),
    ]

    for source_name, search_dir in search_dirs:
        if not search_dir.exists():
            continue

        # If cycle specified, only search that cycle
        if cycle:
            cycle_dirs = [search_dir / cycle]
            if cycle.isdigit():
                cycle_dirs.append(search_dir / f"cyc{int(cycle):03d}")
        else:
            cycle_dirs = sorted(d for d in search_dir.iterdir() if d.is_dir())

        for cycle_dir in cycle_dirs:
            if not cycle_dir.exists():
                continue

            cycle_name = cycle_dir.name

            # Find image files - single directory scan instead of multiple glob calls.
            # Sorted so the listing does not depend on filesystem iteration order.
            valid_extensions = {".tif", ".tiff", ".zarr"}
            valid_suffixes = {".ome.tif", ".ome.tiff"}
            for f in sorted(cycle_dir.iterdir()):
                # Check if file matches our patterns
                name_lower = f.name.lower()
                suffix_lower = f.suffix.lower()
                is_valid = suffix_lower in valid_extensions or any(
                    name_lower.endswith(s) for s in valid_suffixes
                )
                if not is_valid:
                    continue

                # Extract channel name from filename
                channel_name = _extract_channel_name(f.stem)

                key = f"{cycle_name}/{channel_name}"
                if key not in channels:
                    channels[key] = {
                        "cycle": cycle_name,
                        "channel": channel_name,
                        "source": source_name,
                        "path": str(f),
                        "size_mb": f.stat().st_size / (1024 * 1024) if f.is_file() else 0,
                    }

    # Also list currently loaded images
    loaded = list(_loaded_images.keys())

    return {
        "status": "success",
        "project": str(project_path),
        "total_channels": len(channels),
        "channels": list(channels.values()),
        "currently_loaded": loaded,
    }


def _extract_channel_name(filename: str) -> str:
    """Extract channel name from filename."""
    # Common patterns: CD3_cyc001, cyc001_CD3, channel_CD3, etc.
    parts = filename.replace("-", "_").split("_")

    # Look for known marker patterns
    markers = ["DAPI", "CD3", "CD4", "CD8", "CD20", "CD45", "CD68", "CK", "blank", "AF"]

    for part in parts:
        for marker in markers:
            if marker.lower() in part.lower():
                return part

    # Return last meaningful part
    for part in reversed(parts):
        if part and not part.isdigit() and len(part) > 1:
            return part

    return filename


async def save_processed(
    channel: Annotated[str, Field(description="Name of the loaded processed channel to save")],
    output_name: Annotated[
        str | None,
        Field(
            description="Bare filename without extension or directory; None auto-generates "
            "'<channel>_<YYYYmmdd_HHMMSS>'"
        ),
    ] = None,
    format: Annotated[
        Literal["tiff", "tif", "ome-tiff", "zarr"],
        Field(
            description="Output format; 'tiff' and 'tif' both write a plain 16-bit .tiff, "
            "'ome-tiff' adds OME-XML metadata, 'zarr' writes a chunked .zarr store"
        ),
    ] = "tiff",
    project_path: Annotated[
        str | None,
        Field(
            description="KINTSUGI project whose data/processed/signal_isolated directory "
            "receives the file; defaults to the project the channel was loaded from, else "
            "the project the server was started for"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Save the processed channel to the project output directory.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

    # Validate untrusted inputs up front so failures are reported in-band.
    try:
        validate_choice(format, {"tiff", "tif", "ome-tiff", "zarr"}, "format")
        if output_name is not None:
            output_name = safe_filename(output_name, "output_name")
    except (PathSafetyError, ValueError) as e:
        return {"error": str(e)}

    image_info = _loaded_images[channel]
    data = image_info["data"]
    metadata = image_info.get("metadata", {})

    # Determine output path
    if project_path:
        project_path = Path(project_path)
        try:
            from kintsugi.project import KintsugiProject

            project = KintsugiProject.load(project_path)
            output_dir = project.paths.processed / "signal_isolated"
        except Exception:
            output_dir = project_path / "data" / "processed" / "signal_isolated"
    else:
        # Use metadata to find project
        src_project = metadata.get("project") or configured_project_root()
        if not src_project:
            # Never create data/processed/... under whatever directory the client
            # happened to start the server in.
            return {"error": NO_PROJECT_MESSAGE}
        output_dir = Path(src_project) / "data" / "processed" / "signal_isolated"

    output_dir.mkdir(parents=True, exist_ok=True)
    output_dir = output_dir.resolve()

    # Generate output filename
    if output_name is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_name = f"{channel}_{timestamp}"

    # Compute if dask array
    if hasattr(data, "compute"):
        data = data.compute()

    # Save based on format
    try:
        if format == "tiff" or format == "tif":
            import tifffile

            output_path = ensure_within(output_dir / f"{output_name}.tiff", output_dir)
            tifffile.imwrite(str(output_path), data.astype(np.uint16))

        elif format == "ome-tiff":
            import tifffile

            output_path = ensure_within(output_dir / f"{output_name}.ome.tiff", output_dir)
            # Create basic OME-XML
            tifffile.imwrite(
                str(output_path),
                data.astype(np.uint16),
                ome=True,
                metadata={"axes": "YX"},
            )

        elif format == "zarr":
            import zarr

            output_path = ensure_within(output_dir / f"{output_name}.zarr", output_dir)
            z = zarr.open(
                str(output_path),
                mode="w",
                shape=data.shape,
                dtype=data.dtype,
                chunks=(1000, 1000),
            )
            z[:] = data

        else:
            return {"error": f"Unknown format: {format}"}

        # Save processing history alongside
        history_path = output_dir / f"{output_name}_history.json"
        history = {
            "channel": channel,
            "saved": datetime.now().isoformat(),
            "format": format,
            "shape": list(data.shape),
            "dtype": str(data.dtype),
            "processing_steps": _processing_history.get(channel, []),
            "source_metadata": metadata,
        }
        with open(history_path, "w") as f:
            json.dump(history, f, indent=2)

        _loaded_images[channel].setdefault("metadata", {})["saved_path"] = str(output_path)

        return {
            "status": "success",
            "output_path": str(output_path),
            "history_path": str(history_path),
            "format": format,
            "shape": list(data.shape),
        }

    except Exception as e:
        return {"error": f"Failed to save: {e}"}


async def get_processing_history(
    channel: Annotated[
        str, Field(description="Name of a loaded channel whose recorded operations to return")
    ],
) -> dict[str, Any]:
    """
    Get the processing history for a channel.

    Shows all operations applied in order.
    """
    history = _processing_history.get(channel, [])

    # Also get metadata if available
    metadata = {}
    if channel in _loaded_images:
        metadata = _loaded_images[channel].get("metadata", {})

    return {
        "status": "success",
        "channel": channel,
        "steps": history,
        "total_steps": len(history),
        "source_metadata": metadata,
    }


async def suggest_parameters(
    channel: Annotated[str, Field(description="Name of the loaded channel to analyze")],
    blank_channel: Annotated[
        str | None,
        Field(
            description="Name of a loaded blank/autofluorescence channel used to estimate "
            "blank-subtraction parameters; None omits the blank-subtraction recommendation"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Rule-based parameter suggestions from the channel's intensity statistics.

    Uses the same rules as ``suggest_with_learning`` (``kintsugi.mcp.heuristics``)
    without consulting the learning database, so it needs no tissue type or
    marker name. Returns the percentile analysis and one recommendation per
    operation, each with a one-line rationale.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}
    data = _loaded_images[channel]["data"]
    step = sample_step(data)
    blank_sample = None
    if blank_channel:
        if blank_channel not in _loaded_images:
            return {"error": f"Blank channel not loaded: {blank_channel}"}
        blank_sample = sample_image(_loaded_images[blank_channel]["data"], step=step)

    sample = sample_image(data, step=step)
    summary = summarize_intensity(sample)

    recommendations: dict[str, Any] = {}
    for operation in HEURISTIC_OPERATIONS:
        if operation == "blank_subtraction" and blank_sample is None:
            continue
        suggestion = heuristic_suggestion(operation, sample, summary, blank_sample)
        recommendations[operation] = {
            **suggestion["parameters"],
            "rationale": suggestion["rationale"],
        }

    return {
        "status": "success",
        "suggestions": {
            "channel": channel,
            "analysis": {"intensity": summary.to_dict()},
            "recommendations": recommendations,
        },
    }


async def generate_jupyter_cell(
    operation: Annotated[
        Literal["blank_subtraction", "denoise", "clahe", "clean", "gaussian_subtract"],
        Field(description="Processing operation to generate an interactive Kview2 tuner cell for"),
    ],
    channel: Annotated[str, Field(description="Name of the loaded channel to tune")],
    initial_params: Annotated[
        dict[str, Any] | None,
        Field(
            description="Initial widget values keyed by parameter name (e.g. 'blank_channel', "
            "'sigma'); None starts the tuner from its built-in defaults"
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Generate a Jupyter notebook cell for interactive parameter tuning.

    Creates code that uses Kview2 widgets for visual parameter adjustment.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

    # The notebook kernel is a different process from this server, so the cell
    # cannot see the session's in-memory images: it loads the channel from the
    # file it came from (or was saved to with save_processed).
    metadata = _loaded_images[channel].get("metadata", {})
    source = metadata.get("saved_path") or metadata.get("source_path")
    if not source:
        return {
            "error": (
                f"Channel {channel!r} exists only in server memory. Call save_processed "
                "first, then generate the cell for the saved channel."
            )
        }
    initial_params = initial_params or {}

    def value(key: str, default: Any) -> str:
        """repr() keeps the generated code valid for any parameter value."""
        return repr(initial_params.get(key, default))

    def comment(text: str) -> str:
        """Collapse whitespace so a name can never break out of a ``#`` comment."""
        return " ".join(str(text).split())

    channel_note = comment(channel)

    if operation == "blank_subtraction":
        blank_channel = initial_params.get("blank_channel", "")
        blank_source = ""
        if blank_channel in _loaded_images:
            blank_meta = _loaded_images[blank_channel].get("metadata", {})
            blank_source = blank_meta.get("saved_path") or blank_meta.get("source_path") or ""
        code = f"""# Interactive Blank Subtraction Parameter Tuning
from kintsugi.assist.param_tuner import blank_subtraction_tuner

# Image files (the tuner loads them; MCP session names are not visible here)
signal = {source!r}  # {channel_note}
blank = {blank_source!r}  # {comment(blank_channel) or "set the blank/AF channel file"}

# Interactive tuner - adjust sliders and see results in real-time
result = blank_subtraction_tuner(
    signal_channel=signal,
    blank_channel=blank,
    initial_params={{
        "blank_clip_factor": {value("blank_clip_factor", 0)},
        "blank_scale_factor": {value("blank_scale_factor", 0.8)},
        "smooth_low": {value("smooth_low", False)},
        "smooth_high": {value("smooth_high", False)},
        "erosion": {value("erosion", 0)},
    }}
)

# After tuning, print final parameters
print("Final parameters:", result.get_params())
"""

    elif operation == "denoise":
        code = f"""# Interactive Denoising Parameter Tuning
from kintsugi.assist.param_tuner import denoise_tuner

channel = {source!r}  # {channel_note}

# Interactive tuner
result = denoise_tuner(
    channel=channel,
    initial_params={{
        "method": {value("method", "median")},
        "filter_size": {value("filter_size", 3)},
        "percentile": {value("percentile", 10)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "clahe":
        code = f"""# Interactive CLAHE Parameter Tuning
from kintsugi.assist.param_tuner import clahe_tuner

channel = {source!r}  # {channel_note}

# Interactive tuner
result = clahe_tuner(
    channel=channel,
    initial_params={{
        "clip_limit": {value("clip_limit", 0.01)},
        "tile_grid_size": {value("tile_grid_size", 70)},
        "nbins": {value("nbins", 128)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "clean":
        code = f"""# Interactive Background Cleaning
from kintsugi.assist.param_tuner import clean_tuner

channel = {source!r}  # {channel_note}

# Interactive tuner
result = clean_tuner(
    channel=channel,
    initial_params={{
        "background_threshold": {value("background_threshold", 100)},
        "remove_small_objects": {value("remove_small_objects", False)},
        "min_object_size": {value("min_object_size", 30)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "gaussian_subtract":
        code = f"""# Interactive Gaussian Subtraction
from kintsugi.assist.param_tuner import gaussian_subtract_tuner

channel = {source!r}  # {channel_note}

# Interactive tuner
result = gaussian_subtract_tuner(
    channel=channel,
    initial_params={{
        "sigma": {value("sigma", 10)},
        "scale_factor": {value("scale_factor", 0.1)},
    }}
)

print("Final parameters:", result.get_params())
"""

    else:
        return {"error": f"Unknown operation: {operation}"}

    return {
        "status": "success",
        "operation": operation,
        "channel": channel,
        "source_file": source,
        "jupyter_cell": code,
        "instructions": [
            "Copy this code to a Jupyter notebook cell",
            "Run the cell to launch interactive parameter tuning",
            "Adjust sliders/widgets to find optimal parameters",
            "Use the final parameters in subsequent processing",
        ],
    }


async def get_session_state() -> SessionStateResult:
    """
    Get the current session state.

    Shows the project root the server serves, loaded images (shape, dtype,
    metadata), per-channel processing history counts, and memory usage.
    """
    root = configured_project_root()
    loaded = {
        name: {
            "shape": [int(x) for x in info["shape"]] if info.get("shape") else [],
            "dtype": str(info["dtype"]),
            "metadata": info.get("metadata", {}),
        }
        for name, info in _loaded_images.items()
    }

    history_summary = {channel: len(steps) for channel, steps in _processing_history.items()}
    memory = estimate_memory_usage()

    return {
        "status": "success",
        "project_root": str(root) if root is not None else None,
        "loaded_images": loaded,
        "total_loaded": len(loaded),
        "processing_history": history_summary,
        "memory_usage": {
            "total_mb": float(memory["total_mb"]),
            "total_gb": float(memory["total_gb"]),
        },
    }


async def unload_channel(
    channel: Annotated[
        str | None,
        Field(description="Name of the loaded channel to drop from the session"),
    ] = None,
    clear_all: Annotated[
        bool,
        Field(description="Drop every loaded channel and all processing history (default: false)"),
    ] = False,
) -> dict[str, Any]:
    """
    Free session memory by dropping loaded channels.

    Loaded images live in server memory for the whole session; call this when
    a channel (or everything) is no longer needed. Files on disk are untouched.
    """

    before = estimate_memory_usage()["total_mb"]
    if clear_all:
        dropped = sorted(_loaded_images)
        clear_session()
    elif channel is None:
        return {"error": "Pass a channel name, or clear_all=true to drop everything."}
    elif channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}
    else:
        dropped = [channel]
        del _loaded_images[channel]
        _processing_history.pop(channel, None)
    after = estimate_memory_usage()["total_mb"]

    return {
        "status": "success",
        "unloaded": dropped,
        "remaining": sorted(_loaded_images),
        "freed_mb": round(float(before - after), 2),
    }
