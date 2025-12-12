"""
Workflow Tools for KINTSUGI MCP Server

State management, channel discovery, saving results, and Jupyter integration.
"""

from __future__ import annotations

import json
import logging
from datetime import datetime
from pathlib import Path
from typing import Any

import numpy as np

# Import session state from signal_isolation
from kintsugi.mcp.tools.signal_isolation import (
    _loaded_images,
    _processing_history,
)

logger = logging.getLogger("kintsugi.mcp.workflow")


async def list_channels(
    project_path: str,
    cycle: str | None = None,
) -> dict[str, Any]:
    """
    List all available channels in the current project.

    Searches raw and processed directories for image files.
    """
    project_path = Path(project_path)

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
            cycle_dirs = [d for d in search_dir.iterdir() if d.is_dir()]

        for cycle_dir in cycle_dirs:
            if not cycle_dir.exists():
                continue

            cycle_name = cycle_dir.name

            # Find image files - single directory scan instead of multiple glob calls
            valid_extensions = {".tif", ".tiff", ".zarr"}
            valid_suffixes = {".ome.tif", ".ome.tiff"}
            for f in cycle_dir.iterdir():
                # Check if file matches our patterns
                name_lower = f.name.lower()
                suffix_lower = f.suffix.lower()
                is_valid = (
                    suffix_lower in valid_extensions
                    or any(name_lower.endswith(s) for s in valid_suffixes)
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
    channel: str,
    output_name: str | None = None,
    format: str = "tiff",
    project_path: str | None = None,
) -> dict[str, Any]:
    """
    Save the processed channel to the project output directory.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

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
        src_project = metadata.get("project")
        if src_project:
            output_dir = Path(src_project) / "data" / "processed" / "signal_isolated"
        else:
            output_dir = Path.cwd() / "output"

    output_dir.mkdir(parents=True, exist_ok=True)

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

            output_path = output_dir / f"{output_name}.tiff"
            tifffile.imwrite(str(output_path), data.astype(np.uint16))

        elif format == "ome-tiff":
            import tifffile

            output_path = output_dir / f"{output_name}.ome.tiff"
            # Create basic OME-XML
            tifffile.imwrite(
                str(output_path),
                data.astype(np.uint16),
                ome=True,
                metadata={"axes": "YX"},
            )

        elif format == "zarr":
            import zarr

            output_path = output_dir / f"{output_name}.zarr"
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

        return {
            "status": "success",
            "output_path": str(output_path),
            "history_path": str(history_path),
            "format": format,
            "shape": list(data.shape),
        }

    except Exception as e:
        return {"error": f"Failed to save: {e}"}


async def get_processing_history(channel: str) -> dict[str, Any]:
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
    channel: str,
    blank_channel: str | None = None,
) -> dict[str, Any]:
    """
    Analyze channel and suggest processing parameters.

    Uses image statistics to recommend appropriate settings.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

    data = _loaded_images[channel]["data"]

    # Compute sample for analysis
    if hasattr(data, "compute"):
        step = max(1, data.shape[0] // 2000)
        sample = data[::step, ::step].compute()
    else:
        sample = data

    suggestions = {
        "channel": channel,
        "analysis": {},
        "recommendations": {},
    }

    # Analyze intensity distribution
    p1, p10, p50, p90, p99 = np.percentile(sample, [1, 10, 50, 90, 99])
    suggestions["analysis"]["intensity"] = {
        "p1": float(p1),
        "p10": float(p10),
        "p50": float(p50),
        "p90": float(p90),
        "p99": float(p99),
        "dynamic_range": float(p99 - p1),
    }

    # Background estimation
    float(p10)
    float(p90)

    # Blank subtraction parameters
    if blank_channel and blank_channel in _loaded_images:
        blank_data = _loaded_images[blank_channel]["data"]
        if hasattr(blank_data, "compute"):
            blank_sample = blank_data[::step, ::step].compute()
        else:
            blank_sample = blank_data

        blank_mean = float(np.mean(blank_sample))
        signal_mean = float(np.mean(sample))

        if blank_mean > 0:
            suggested_scale = min(1.0, signal_mean / blank_mean * 0.8)
        else:
            suggested_scale = 0.8

        suggestions["recommendations"]["blank_subtraction"] = {
            "blank_scale_factor": round(suggested_scale, 2),
            "blank_clip_factor": int(np.percentile(blank_sample, 5)),
            "rationale": f"Blank mean={blank_mean:.0f}, Signal mean={signal_mean:.0f}",
        }

    # Denoising parameters
    noise_estimate = float(np.std(sample[sample < p10]))
    if noise_estimate > 500:
        suggestions["recommendations"]["denoise"] = {
            "method": "median",
            "filter_size": 3,
            "rationale": f"High noise detected (std={noise_estimate:.0f})",
        }
    elif noise_estimate > 200:
        suggestions["recommendations"]["denoise"] = {
            "method": "median",
            "filter_size": 2,
            "rationale": f"Moderate noise detected (std={noise_estimate:.0f})",
        }

    # CLAHE parameters
    dynamic_range = p99 - p1
    if dynamic_range < 10000:
        suggestions["recommendations"]["clahe"] = {
            "clip_limit": 0.02,
            "tile_grid_size": 50,
            "rationale": f"Low dynamic range ({dynamic_range:.0f}), aggressive CLAHE",
        }
    elif dynamic_range < 30000:
        suggestions["recommendations"]["clahe"] = {
            "clip_limit": 0.01,
            "tile_grid_size": 70,
            "rationale": f"Moderate dynamic range ({dynamic_range:.0f})",
        }

    # Background cleaning
    if p10 > 500:
        suggestions["recommendations"]["clean_background"] = {
            "background_threshold": int(p10 * 0.5),
            "remove_small_objects": True,
            "min_object_size": 50,
            "rationale": f"Elevated background (p10={p10:.0f})",
        }

    return {
        "status": "success",
        "suggestions": suggestions,
    }


async def generate_jupyter_cell(
    operation: str,
    channel: str,
    initial_params: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Generate a Jupyter notebook cell for interactive parameter tuning.

    Creates code that uses Kview2 widgets for visual parameter adjustment.
    """
    if channel not in _loaded_images:
        return {"error": f"Channel not loaded: {channel}"}

    _loaded_images[channel].get("metadata", {})
    initial_params = initial_params or {}

    # Generate code based on operation
    if operation == "blank_subtraction":
        blank_channel = initial_params.get("blank_channel", "")
        code = f"""# Interactive Blank Subtraction Parameter Tuning
from kintsugi.claude.param_tuner import blank_subtraction_tuner

# Load your channels first if not already loaded
signal = "{channel}"  # Signal channel
blank = "{blank_channel}"  # Blank/AF channel

# Interactive tuner - adjust sliders and see results in real-time
result = blank_subtraction_tuner(
    signal_channel=signal,
    blank_channel=blank,
    initial_params={{
        "blank_clip_factor": {initial_params.get("blank_clip_factor", 0)},
        "blank_scale_factor": {initial_params.get("blank_scale_factor", 0.8)},
        "smooth_low": {initial_params.get("smooth_low", False)},
        "smooth_high": {initial_params.get("smooth_high", False)},
        "erosion": {initial_params.get("erosion", 0)},
    }}
)

# After tuning, print final parameters
print("Final parameters:", result.get_params())
"""

    elif operation == "denoise":
        code = f"""# Interactive Denoising Parameter Tuning
from kintsugi.claude.param_tuner import denoise_tuner

channel = "{channel}"

# Interactive tuner
result = denoise_tuner(
    channel=channel,
    initial_params={{
        "method": "{initial_params.get("method", "median")}",
        "filter_size": {initial_params.get("filter_size", 3)},
        "percentile": {initial_params.get("percentile", 10)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "clahe":
        code = f"""# Interactive CLAHE Parameter Tuning
from kintsugi.claude.param_tuner import clahe_tuner

channel = "{channel}"

# Interactive tuner
result = clahe_tuner(
    channel=channel,
    initial_params={{
        "clip_limit": {initial_params.get("clip_limit", 0.01)},
        "tile_grid_size": {initial_params.get("tile_grid_size", 70)},
        "nbins": {initial_params.get("nbins", 128)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "clean":
        code = f"""# Interactive Background Cleaning
from kintsugi.claude.param_tuner import clean_tuner

channel = "{channel}"

# Interactive tuner
result = clean_tuner(
    channel=channel,
    initial_params={{
        "background_threshold": {initial_params.get("background_threshold", 100)},
        "remove_small_objects": {initial_params.get("remove_small_objects", False)},
        "min_object_size": {initial_params.get("min_object_size", 30)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "gaussian_subtract":
        code = f"""# Interactive Gaussian Subtraction
from kintsugi.claude.param_tuner import gaussian_subtract_tuner

channel = "{channel}"

# Interactive tuner
result = gaussian_subtract_tuner(
    channel=channel,
    initial_params={{
        "sigma": {initial_params.get("sigma", 10)},
        "scale_factor": {initial_params.get("scale_factor", 0.1)},
    }}
)

print("Final parameters:", result.get_params())
"""

    elif operation == "full_pipeline":
        code = f"""# Full Signal Isolation Pipeline with Interactive Review
from kintsugi.claude.pipeline import SignalIsolationPipeline
from kintsugi.claude.param_tuner import pipeline_tuner

channel = "{channel}"
blank_channel = "{initial_params.get("blank_channel", "")}"

# Create pipeline with suggested parameters
pipeline = SignalIsolationPipeline(
    signal_channel=channel,
    blank_channel=blank_channel,
)

# Run with interactive review at each step
result = pipeline.run_interactive()

# Final result
print("Pipeline complete!")
print("Quality score:", result.quality_score)
print("All parameters:", result.get_all_params())

# Save if satisfied
# result.save("output_name")
"""

    else:
        return {"error": f"Unknown operation: {operation}"}

    return {
        "status": "success",
        "operation": operation,
        "channel": channel,
        "jupyter_cell": code,
        "instructions": [
            "Copy this code to a Jupyter notebook cell",
            "Run the cell to launch interactive parameter tuning",
            "Adjust sliders/widgets to find optimal parameters",
            "Use the final parameters in subsequent processing",
        ],
    }


async def set_project_context(project_path: str) -> dict[str, Any]:
    """
    Set the current project context for the session.

    Loads project configuration and prepares for processing.
    """
    project_path = Path(project_path)

    try:
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_path)

        return {
            "status": "success",
            "project": {
                "name": project.config.name,
                "path": str(project_path),
                "cycles": project.list_cycles(),
                "raw_dir": str(project.paths.raw),
                "processed_dir": str(project.paths.processed),
            },
        }

    except FileNotFoundError:
        return {
            "status": "not_found",
            "message": f"No KINTSUGI project at {project_path}",
            "suggestion": "Use KintsugiProject.create() to initialize a new project",
        }

    except Exception as e:
        return {"error": str(e)}


async def get_session_state() -> dict[str, Any]:
    """
    Get the current session state.

    Shows loaded images, processing history, and session info.
    """
    loaded = {
        name: {
            "shape": info["shape"],
            "dtype": info["dtype"],
            "metadata": info.get("metadata", {}),
        }
        for name, info in _loaded_images.items()
    }

    history_summary = {channel: len(steps) for channel, steps in _processing_history.items()}

    return {
        "status": "success",
        "loaded_images": loaded,
        "total_loaded": len(loaded),
        "processing_history": history_summary,
        "memory_usage": _estimate_memory_usage(),
    }


def _estimate_memory_usage() -> dict[str, float]:
    """Estimate memory usage of loaded images."""
    total_bytes = 0

    for _name, info in _loaded_images.items():
        data = info["data"]
        if hasattr(data, "nbytes"):
            total_bytes += data.nbytes
        elif info.get("shape"):
            # Estimate from shape
            shape = info["shape"]
            dtype = info.get("dtype", "uint16")
            itemsize = {"uint8": 1, "uint16": 2, "float32": 4, "float64": 8}.get(dtype, 2)
            total_bytes += np.prod(shape) * itemsize

    return {
        "total_mb": total_bytes / (1024 * 1024),
        "total_gb": total_bytes / (1024 * 1024 * 1024),
    }
