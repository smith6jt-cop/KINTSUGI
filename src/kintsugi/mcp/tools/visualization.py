"""
Visualization Tools for KINTSUGI MCP Server

Provides image statistics and thumbnail generation for Claude to understand images.
"""

from __future__ import annotations

import base64
import io
import logging
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger("kintsugi.mcp.visualization")

# Import session state from signal_isolation
from kintsugi.mcp.tools.signal_isolation import _loaded_images


def _get_image(name: str) -> Any:
    """Get a loaded image by name."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


async def get_image_stats(channel: str) -> Dict[str, Any]:
    """
    Get comprehensive statistics for a loaded image.

    Returns intensity distribution, percentiles, and histogram data.
    """
    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute if dask array (sample for large images)
    if hasattr(data, "compute"):
        # Sample every nth pixel for efficiency
        step = max(1, data.shape[0] // 2000)
        data_sample = data[::step, ::step].compute()
    else:
        data_sample = data

    # Basic statistics
    stats = {
        "shape": list(data.shape) if hasattr(data, "shape") else None,
        "dtype": str(data.dtype) if hasattr(data, "dtype") else None,
        "min": float(np.min(data_sample)),
        "max": float(np.max(data_sample)),
        "mean": float(np.mean(data_sample)),
        "std": float(np.std(data_sample)),
        "median": float(np.median(data_sample)),
    }

    # Percentiles
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    stats["percentiles"] = {
        f"p{p}": float(np.percentile(data_sample, p)) for p in percentiles
    }

    # Histogram (for understanding intensity distribution)
    hist, bin_edges = np.histogram(data_sample.ravel(), bins=50)
    stats["histogram"] = {
        "counts": hist.tolist(),
        "bin_edges": bin_edges.tolist(),
        "bin_centers": ((bin_edges[:-1] + bin_edges[1:]) / 2).tolist(),
    }

    # Zero and saturation counts
    total_pixels = data_sample.size
    stats["zero_pixels"] = int(np.sum(data_sample == 0))
    stats["zero_ratio"] = stats["zero_pixels"] / total_pixels

    max_val = 65535 if data_sample.dtype == np.uint16 else 255
    stats["saturated_pixels"] = int(np.sum(data_sample >= max_val * 0.99))
    stats["saturation_ratio"] = stats["saturated_pixels"] / total_pixels

    # Non-zero statistics (for sparse markers)
    nonzero = data_sample[data_sample > 0]
    if len(nonzero) > 0:
        stats["nonzero_mean"] = float(np.mean(nonzero))
        stats["nonzero_std"] = float(np.std(nonzero))
        stats["nonzero_median"] = float(np.median(nonzero))
        stats["nonzero_ratio"] = len(nonzero) / total_pixels
    else:
        stats["nonzero_mean"] = 0
        stats["nonzero_std"] = 0
        stats["nonzero_median"] = 0
        stats["nonzero_ratio"] = 0

    return {
        "status": "success",
        "channel": channel,
        "statistics": stats,
    }


async def get_thumbnail(
    channel: str,
    max_size: int = 512,
    normalize: bool = True,
    format: str = "png",
) -> Dict[str, Any]:
    """
    Get a downsampled thumbnail of the image.

    Returns base64-encoded image data that Claude can view.
    """
    try:
        from PIL import Image
    except ImportError:
        return {"error": "PIL/Pillow not installed"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute if dask array
    if hasattr(data, "compute"):
        # Downsample while computing
        h, w = data.shape[:2]
        scale = max(h, w) / max_size
        if scale > 1:
            step = int(scale)
            data = data[::step, ::step].compute()
        else:
            data = data.compute()
    else:
        data = data.copy()

    # Handle 3D arrays (take middle slice or first channel)
    if data.ndim == 3:
        if data.shape[0] < data.shape[2]:
            # Likely CxHxW format
            data = data[0]
        else:
            # Likely HxWxC format - convert to grayscale
            data = data[:, :, 0]

    # Normalize to 0-255 for display
    if normalize:
        p_low, p_high = np.percentile(data, [1, 99])
        if p_high > p_low:
            data = np.clip((data - p_low) / (p_high - p_low) * 255, 0, 255)
        else:
            data = np.zeros_like(data)
        data = data.astype(np.uint8)
    else:
        # Scale to 8-bit
        if data.dtype == np.uint16:
            data = (data / 256).astype(np.uint8)
        elif data.dtype == np.float64 or data.dtype == np.float32:
            data = (np.clip(data, 0, 1) * 255).astype(np.uint8)

    # Resize if still too large
    h, w = data.shape[:2]
    if max(h, w) > max_size:
        scale = max_size / max(h, w)
        new_h, new_w = int(h * scale), int(w * scale)
        img = Image.fromarray(data)
        img = img.resize((new_w, new_h), Image.Resampling.LANCZOS)
    else:
        img = Image.fromarray(data)

    # Convert to base64
    buffer = io.BytesIO()
    img.save(buffer, format=format.upper())
    img_base64 = base64.b64encode(buffer.getvalue()).decode("utf-8")

    return {
        "status": "success",
        "channel": channel,
        "thumbnail": {
            "width": img.width,
            "height": img.height,
            "format": format,
            "data_uri": f"data:image/{format};base64,{img_base64}",
        },
        "original_shape": list(_loaded_images[channel]["data"].shape),
    }


async def compare_channels(
    channel1: str,
    channel2: str,
) -> Dict[str, Any]:
    """
    Compare two channels statistically.

    Useful for before/after comparison or cross-channel analysis.
    """
    try:
        data1 = _get_image(channel1)
        data2 = _get_image(channel2)
    except ValueError as e:
        return {"error": str(e)}

    # Compute samples
    if hasattr(data1, "compute"):
        step = max(1, data1.shape[0] // 2000)
        data1 = data1[::step, ::step].compute()
    if hasattr(data2, "compute"):
        step = max(1, data2.shape[0] // 2000)
        data2 = data2[::step, ::step].compute()

    # Ensure same shape for comparison
    min_h = min(data1.shape[0], data2.shape[0])
    min_w = min(data1.shape[1], data2.shape[1])
    data1 = data1[:min_h, :min_w]
    data2 = data2[:min_h, :min_w]

    # Calculate difference statistics
    diff = data1.astype(float) - data2.astype(float)

    comparison = {
        "channel1": {
            "name": channel1,
            "mean": float(np.mean(data1)),
            "std": float(np.std(data1)),
            "max": float(np.max(data1)),
        },
        "channel2": {
            "name": channel2,
            "mean": float(np.mean(data2)),
            "std": float(np.std(data2)),
            "max": float(np.max(data2)),
        },
        "difference": {
            "mean": float(np.mean(diff)),
            "std": float(np.std(diff)),
            "min": float(np.min(diff)),
            "max": float(np.max(diff)),
            "abs_mean": float(np.mean(np.abs(diff))),
        },
    }

    # Correlation
    if data1.std() > 0 and data2.std() > 0:
        correlation = np.corrcoef(data1.ravel(), data2.ravel())[0, 1]
        comparison["correlation"] = float(correlation)
    else:
        comparison["correlation"] = 0.0

    # Structural similarity (simplified)
    mean1, mean2 = data1.mean(), data2.mean()
    std1, std2 = data1.std(), data2.std()
    cov = np.mean((data1 - mean1) * (data2 - mean2))

    c1 = (0.01 * 65535) ** 2
    c2 = (0.03 * 65535) ** 2

    ssim = ((2 * mean1 * mean2 + c1) * (2 * cov + c2)) / (
        (mean1**2 + mean2**2 + c1) * (std1**2 + std2**2 + c2)
    )
    comparison["ssim"] = float(ssim)

    # Interpretation
    if comparison["correlation"] > 0.95:
        interpretation = "Channels are very similar"
    elif comparison["correlation"] > 0.8:
        interpretation = "Channels are moderately correlated"
    elif comparison["correlation"] > 0.5:
        interpretation = "Channels show some correlation"
    else:
        interpretation = "Channels are largely independent"

    comparison["interpretation"] = interpretation

    return {
        "status": "success",
        "comparison": comparison,
    }


async def get_intensity_profile(
    channel: str,
    axis: str = "horizontal",
    position: Optional[float] = None,
) -> Dict[str, Any]:
    """
    Get intensity profile along a line through the image.

    Useful for checking for gradients or uneven illumination.
    """
    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute if dask array
    if hasattr(data, "compute"):
        data = data.compute()

    h, w = data.shape[:2]

    # Default position is center
    if position is None:
        position = 0.5

    if axis == "horizontal":
        y = int(h * position)
        profile = data[y, :].astype(float)
        x_coords = list(range(w))
    else:  # vertical
        x = int(w * position)
        profile = data[:, x].astype(float)
        x_coords = list(range(h))

    # Downsample if too many points
    max_points = 500
    if len(profile) > max_points:
        step = len(profile) // max_points
        profile = profile[::step]
        x_coords = x_coords[::step]

    return {
        "status": "success",
        "channel": channel,
        "axis": axis,
        "position": position,
        "profile": {
            "x": x_coords,
            "values": profile.tolist(),
            "mean": float(np.mean(profile)),
            "std": float(np.std(profile)),
            "min": float(np.min(profile)),
            "max": float(np.max(profile)),
        },
    }
