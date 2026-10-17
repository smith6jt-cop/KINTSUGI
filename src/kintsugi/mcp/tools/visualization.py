"""
Visualization Tools for KINTSUGI MCP Server

Provides image statistics and thumbnails so an assistant can inspect images.
"""

from __future__ import annotations

import base64
import io
import json
import logging
from typing import Annotated, Any, Literal

import numpy as np
from pydantic import Field

from kintsugi.mcp._sdk import CallToolResult, ImageContent, TextContent
from kintsugi.mcp.path_safety import clamp_float, clamp_int, validate_choice
from kintsugi.mcp.schemas import ImageStatsResult
from kintsugi.mcp.session import _get_image

logger = logging.getLogger("kintsugi.mcp.visualization")


def _as_2d(data: Any, name: str) -> Any:
    """Return ``data`` as a 2-D image (a leading singleton axis is dropped) or raise ValueError."""
    ndim = getattr(data, "ndim", None)
    if ndim == 3 and data.shape[0] == 1:
        data = data[0]
        ndim = 2
    if ndim != 2:
        shape = tuple(getattr(data, "shape", ()))
        raise ValueError(f"{name} must be a 2-D image for this tool, got shape {shape}")
    return data


def _ssim_dynamic_range(*arrays: np.ndarray) -> float:
    """The intensity range L that stabilises SSIM: the dtype maximum for integer
    images, the observed maximum for float images."""
    largest = 0.0
    for a in arrays:
        if np.issubdtype(a.dtype, np.integer):
            largest = max(largest, float(np.iinfo(a.dtype).max))
        elif a.size:
            largest = max(largest, float(np.max(a)))
    return largest if largest > 0 else 1.0


async def get_image_stats(
    channel: Annotated[str, Field(description="Name of loaded channel")],
) -> ImageStatsResult:
    """
    Get comprehensive statistics for a loaded image.

    Returns intensity distribution, percentiles, and histogram data. The return
    annotation is the advertised ``outputSchema``; every value must be a plain
    Python scalar/list so the SDK's output validation and JSON serialisation
    accept it (numpy integers are not JSON-serialisable by pydantic).
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
    stats: dict[str, Any] = {
        "shape": [int(s) for s in data.shape] if hasattr(data, "shape") else None,
        "dtype": str(data.dtype) if hasattr(data, "dtype") else None,
        "min": float(np.min(data_sample)),
        "max": float(np.max(data_sample)),
        "mean": float(np.mean(data_sample)),
        "std": float(np.std(data_sample)),
        "median": float(np.median(data_sample)),
    }

    # Percentiles
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    stats["percentiles"] = {f"p{p}": float(np.percentile(data_sample, p)) for p in percentiles}

    # Histogram (for understanding intensity distribution)
    hist, bin_edges = np.histogram(data_sample.ravel(), bins=50)
    stats["histogram"] = {
        "counts": hist.tolist(),
        "bin_edges": bin_edges.tolist(),
        "bin_centers": ((bin_edges[:-1] + bin_edges[1:]) / 2).tolist(),
    }

    # Zero and saturation counts
    total_pixels = int(data_sample.size)
    stats["zero_pixels"] = int(np.sum(data_sample == 0))
    stats["zero_ratio"] = float(stats["zero_pixels"] / total_pixels)

    max_val = 65535 if data_sample.dtype == np.uint16 else 255
    stats["saturated_pixels"] = int(np.sum(data_sample >= max_val * 0.99))
    stats["saturation_ratio"] = float(stats["saturated_pixels"] / total_pixels)

    # Non-zero statistics (for sparse markers)
    nonzero = data_sample[data_sample > 0]
    if len(nonzero) > 0:
        stats["nonzero_mean"] = float(np.mean(nonzero))
        stats["nonzero_std"] = float(np.std(nonzero))
        stats["nonzero_median"] = float(np.median(nonzero))
        stats["nonzero_ratio"] = float(len(nonzero) / total_pixels)
    else:
        stats["nonzero_mean"] = 0.0
        stats["nonzero_std"] = 0.0
        stats["nonzero_median"] = 0.0
        stats["nonzero_ratio"] = 0.0

    return {
        "status": "success",
        "channel": channel,
        "statistics": stats,
    }


def _error_result(message: str) -> CallToolResult:
    """In-band error in the ``CallToolResult`` shape (never raised to the SDK)."""
    err = {"error": message}
    return CallToolResult(
        content=[TextContent(type="text", text=json.dumps(err))],
        structuredContent=err,
        isError=True,
    )


def _render_thumbnail(
    data: Any,
    max_size: int,
    normalize: bool,
    fmt: str,
) -> tuple[bytes, dict[str, Any]]:
    """Downsample, normalise and encode an image for display.

    Pure helper behind :func:`get_thumbnail`: it reads nothing from the session
    state and raises (``ImportError`` without Pillow, ``ValueError`` on bad
    input) instead of returning in-band errors. ``data`` may be a numpy or dask
    array; ``fmt`` is ``"png"`` or ``"jpeg"`` (``"jpg"`` is accepted as an
    alias). Returns the encoded bytes and ``{"width", "height", "format"}``.
    """
    from PIL import Image

    # Canonicalize so the accepted "jpg" alias does not fail in Pillow, which
    # only recognizes "JPEG" for that format.
    fmt = "jpeg" if fmt.lower() == "jpg" else fmt.lower()

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

    buffer = io.BytesIO()
    img.save(buffer, format=fmt.upper())
    return buffer.getvalue(), {"width": int(img.width), "height": int(img.height), "format": fmt}


async def get_thumbnail(
    channel: Annotated[str, Field(description="Name of loaded channel")],
    max_size: Annotated[
        int,
        Field(
            description="Maximum dimension of thumbnail in pixels (default: 512)", ge=16, le=4096
        ),
    ] = 512,
    normalize: Annotated[
        bool,
        Field(description="Percentile-stretch 1-99% to 8-bit before encoding (default: true)"),
    ] = True,
    format: Annotated[
        Literal["png", "jpeg", "jpg"],
        Field(
            description="Encoding of the returned image; 'jpg' is an alias of 'jpeg' (default: png)"
        ),
    ] = "png",
    include_data_uri: Annotated[
        bool,
        Field(
            description="Also embed the image as a data: URI in the structured result "
            "(large; off by default)"
        ),
    ] = False,
) -> Annotated[CallToolResult, dict[str, Any]]:
    """
    Get a downsampled thumbnail of the image.

    The thumbnail is delivered as a renderable ``image`` content block (so the
    client can display it) followed by a JSON text block; the structured
    result carries the thumbnail dimensions, format and the original shape, and
    the base64 ``data_uri`` only when ``include_data_uri`` is set. Errors are
    returned in-band as ``{"error": ...}`` with ``isError`` set.
    """
    try:
        clamp_int(max_size, "max_size", minimum=16, maximum=4096)
        validate_choice(format.lower(), {"png", "jpeg", "jpg"}, "format")
    except ValueError as e:
        return _error_result(str(e))

    fmt = "jpeg" if format.lower() == "jpg" else format.lower()

    try:
        import PIL  # noqa: F401
    except ImportError:
        return _error_result("PIL/Pillow not installed")

    try:
        data = _get_image(channel)
    except ValueError as e:
        return _error_result(str(e))

    original_shape = [int(s) for s in data.shape]
    encoded, thumbnail = _render_thumbnail(data, max_size, normalize, fmt)
    img_base64 = base64.b64encode(encoded).decode("ascii")
    if include_data_uri:
        thumbnail["data_uri"] = f"data:image/{fmt};base64,{img_base64}"

    meta: dict[str, Any] = {
        "status": "success",
        "channel": channel,
        "thumbnail": thumbnail,
        "original_shape": original_shape,
    }
    return CallToolResult(
        content=[
            ImageContent(type="image", data=img_base64, mimeType=f"image/{fmt}"),
            TextContent(type="text", text=json.dumps(meta)),
        ],
        structuredContent=meta,
    )


async def compare_channels(
    channel1: Annotated[str, Field(description="Name of the first loaded channel")],
    channel2: Annotated[
        str,
        Field(
            description="Name of the second loaded channel, e.g. the processed version of channel1"
        ),
    ],
) -> dict[str, Any]:
    """
    Compare two loaded channels statistically.

    Reports per-channel mean/std/max, the difference image statistics, Pearson
    correlation and a global SSIM, with a one-line interpretation. Useful for
    before/after comparison or cross-channel analysis.
    """
    try:
        data1 = _as_2d(_get_image(channel1), channel1)
        data2 = _as_2d(_get_image(channel2), channel2)
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

    # Stabilisers scale with the image range (uint8, uint16 or float data alike).
    intensity_range = _ssim_dynamic_range(data1, data2)
    c1 = (0.01 * intensity_range) ** 2
    c2 = (0.03 * intensity_range) ** 2

    ssim = ((2 * mean1 * mean2 + c1) * (2 * cov + c2)) / (
        (mean1**2 + mean2**2 + c1) * (std1**2 + std2**2 + c2)
    )
    comparison["ssim"] = float(ssim)
    comparison["ssim_intensity_range"] = intensity_range

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
    channel: Annotated[str, Field(description="Name of the loaded channel")],
    axis: Annotated[
        Literal["horizontal", "vertical"],
        Field(description="Line direction: a row (horizontal) or a column (vertical) of the image"),
    ] = "horizontal",
    position: Annotated[
        float | None,
        Field(
            description="Where the line sits along the other axis as a fraction 0.0-1.0; "
            "None = 0.5 (the centre)",
            ge=0.0,
            le=1.0,
            allow_inf_nan=False,
        ),
    ] = None,
) -> dict[str, Any]:
    """
    Intensity values along a horizontal or vertical line through the image.

    Useful for spotting gradients or uneven illumination. Long lines are
    down-sampled to at most 500 points.
    """
    try:
        validate_choice(axis, ("horizontal", "vertical"), "axis")
        if position is not None:
            clamp_float(position, "position", minimum=0.0, maximum=1.0)
        data = _as_2d(_get_image(channel), channel)
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
        y = min(h - 1, int(h * position))  # position=1.0 is the last row, not one past it
        profile = data[y, :].astype(float)
        x_coords = list(range(w))
    else:  # vertical
        x = min(w - 1, int(w * position))
        profile = data[:, x].astype(float)
        x_coords = list(range(h))

    # Down-sample long lines to exactly max_points evenly spaced positions.
    max_points = 500
    if len(profile) > max_points:
        index = np.linspace(0, len(profile) - 1, max_points).round().astype(int)
        profile = profile[index]
        x_coords = index.tolist()

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
