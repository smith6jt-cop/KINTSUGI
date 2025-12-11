"""
Signal Isolation Tools for KINTSUGI MCP Server

Wraps Kutils functions for blank subtraction, denoising, and signal enhancement.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Dict, Optional

import numpy as np

logger = logging.getLogger("kintsugi.mcp.signal_isolation")

# Session state - shared with server
_loaded_images: Dict[str, Any] = {}
_processing_history: Dict[str, list] = {}


def _get_image(name: str) -> Any:
    """Get a loaded image by name."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


def _store_image(name: str, data: Any, metadata: Dict[str, Any] = None):
    """Store an image with metadata."""
    _loaded_images[name] = {
        "data": data,
        "metadata": metadata or {},
        "dtype": str(data.dtype) if hasattr(data, "dtype") else "unknown",
        "shape": data.shape if hasattr(data, "shape") else None,
    }
    if name not in _processing_history:
        _processing_history[name] = []


def _add_history(name: str, operation: str, params: Dict[str, Any]):
    """Add an operation to processing history."""
    if name not in _processing_history:
        _processing_history[name] = []
    _processing_history[name].append({
        "operation": operation,
        "params": params,
    })


async def load_channel(
    project_path: str,
    cycle: str,
    channel: str,
) -> Dict[str, Any]:
    """
    Load a channel image from a KINTSUGI project.

    Returns image metadata and prepares it for processing.
    """
    from pathlib import Path

    try:
        import dask.array as da
        import zarr
        from skimage import io as skio
        import tifffile
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    project_path = Path(project_path)

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
    _store_image(image_name, data, {
        "source_path": str(image_path),
        "cycle": cycle_name,
        "channel": channel,
        "project": str(project_path),
    })

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


async def subtract_blank(
    signal_channel: str,
    blank_channel: str,
    blank_clip_factor: int = 0,
    blank_scale_factor: float = 1.0,
    smooth_low: bool = False,
    low_size: int = 1,
    low_percentile: int = 60,
    smooth_high: bool = False,
    high_size: int = 1,
    high_percentile: int = 90,
    erosion: int = 0,
) -> Dict[str, Any]:
    """
    Subtract autofluorescence/blank channel from signal.

    Uses the ini_params function from Kutils.
    """
    try:
        import dask.array as da
        import dask_image.ndfilters
        import dask_image.ndmorph
        from skimage import morphology
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    # Get loaded images
    try:
        signal_data = _get_image(signal_channel)
        blank_data = _get_image(blank_channel)
    except ValueError as e:
        return {"error": str(e)}

    # Ensure compatible types
    blank_copy = da.Array.copy(blank_data)
    signal_copy = signal_data.astype(blank_data.dtype)

    # Apply blank subtraction (from Kutils.ini_params)
    blank_copy = da.clip(blank_copy, blank_clip_factor, blank_copy.max())
    blank_copy = da.where(blank_copy <= blank_clip_factor, 0, blank_copy)
    result = signal_copy - da.minimum(signal_copy, blank_copy * blank_scale_factor)

    # Low intensity smoothing
    if smooth_low:
        low_mask = da.where(
            signal_copy < da.percentile(signal_copy.ravel(), low_percentile), 1, 0
        )
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
    _store_image(output_name, result, {
        "source": signal_channel,
        "blank": blank_channel,
        "operation": "subtract_blank",
    })

    # Record history
    _add_history(output_name, "subtract_blank", {
        "signal_channel": signal_channel,
        "blank_channel": blank_channel,
        "blank_clip_factor": blank_clip_factor,
        "blank_scale_factor": blank_scale_factor,
        "smooth_low": smooth_low,
        "smooth_high": smooth_high,
        "erosion": erosion,
    })

    # Compute stats for result
    try:
        sample = result[::10, ::10].compute()
        stats = {
            "min": float(np.min(sample)),
            "max": float(np.max(sample)),
            "mean": float(np.mean(sample)),
            "std": float(np.std(sample)),
        }
    except Exception:
        stats = {}

    return {
        "status": "success",
        "output_name": output_name,
        "shape": list(result.shape),
        "stats": stats,
        "parameters_used": {
            "blank_clip_factor": blank_clip_factor,
            "blank_scale_factor": blank_scale_factor,
            "smooth_low": smooth_low,
            "smooth_high": smooth_high,
            "erosion": erosion,
        },
    }


async def denoise(
    channel: str,
    method: str = "median",
    filter_size: int = 3,
    percentile: int = 10,
    upper: bool = False,
) -> Dict[str, Any]:
    """
    Apply denoising filters to reduce noise.

    Methods:
    - percentile: Percentile filter (good for salt-and-pepper noise)
    - uniform: Uniform/box filter (smoothing)
    - median: Median filter (edge-preserving)
    """
    try:
        import dask.array as da
        import dask_image.ndfilters
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    data_copy = da.Array.copy(data)

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
    _store_image(output_name, result, {
        "source": channel,
        "operation": "denoise",
        "method": method,
    })

    _add_history(output_name, "denoise", {
        "method": method,
        "filter_size": filter_size,
        "percentile": percentile if method == "percentile" else None,
    })

    return {
        "status": "success",
        "output_name": output_name,
        "method": method,
        "filter_size": filter_size,
    }


async def apply_clahe(
    channel: str,
    clip_limit: float = 0.01,
    tile_grid_size: int = 70,
    nbins: int = 128,
) -> Dict[str, Any]:
    """
    Apply Contrast Limited Adaptive Histogram Equalization.

    CLAHE enhances local contrast while limiting noise amplification.
    """
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
    _store_image(output_name, result, {
        "source": channel,
        "operation": "clahe",
    })

    _add_history(output_name, "clahe", {
        "clip_limit": clip_limit,
        "tile_grid_size": tile_grid_size,
        "nbins": nbins,
        "actual_kernel_size": kernel_size,
    })

    return {
        "status": "success",
        "output_name": output_name,
        "tile_size": f"{kernel_size}x{kernel_size}",
        "clip_limit": clip_limit,
    }


async def clean_background(
    channel: str,
    background_threshold: int = 100,
    smooth: bool = False,
    smooth_threshold: int = 100,
    remove_small_objects: bool = False,
    min_object_size: int = 30,
) -> Dict[str, Any]:
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

    result = da.Array.copy(data)

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
    _store_image(output_name, result, {
        "source": channel,
        "operation": "clean_background",
    })

    _add_history(output_name, "clean_background", {
        "background_threshold": background_threshold,
        "smooth": smooth,
        "remove_small_objects": remove_small_objects,
        "min_object_size": min_object_size,
    })

    return {
        "status": "success",
        "output_name": output_name,
        "background_threshold": background_threshold,
        "small_objects_removed": remove_small_objects,
    }


async def gaussian_subtract(
    channel: str,
    sigma: int = 10,
    scale_factor: float = 0.1,
    low_clip: int = 0,
    high_clip: int = 65535,
) -> Dict[str, Any]:
    """
    Subtract Gaussian-blurred version to remove structured background.
    """
    try:
        import dask.array as da
        import dask_image.ndfilters
    except ImportError as e:
        return {"error": f"Missing dependency: {e}"}

    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Create Gaussian-blurred version
    blurred = dask_image.ndfilters.gaussian_filter(data, sigma=sigma)

    # Clip extreme values
    blurred = da.where(low_clip >= blurred, 0, blurred)
    blurred = da.where(blurred >= high_clip, 0, blurred)

    # Subtract scaled blurred from original
    result = data - da.minimum(data, blurred * scale_factor)

    # Store result
    output_name = f"{channel}_gaussub"
    _store_image(output_name, result, {
        "source": channel,
        "operation": "gaussian_subtract",
    })

    _add_history(output_name, "gaussian_subtract", {
        "sigma": sigma,
        "scale_factor": scale_factor,
        "low_clip": low_clip,
        "high_clip": high_clip,
    })

    return {
        "status": "success",
        "output_name": output_name,
        "sigma": sigma,
        "scale_factor": scale_factor,
    }


async def denoise_advanced(
    channel: str,
    method: str = "adaptive",
    strength: str = "auto",
    preserve_edges: bool = True,
    n2v_epochs: int = 50,
    model_path: Optional[str] = None,
) -> Dict[str, Any]:
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
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Convert to numpy if needed
    if hasattr(data, "compute"):
        # For large images, process a reasonable subset or use tiled approach
        if data.shape[0] > 4000 or data.shape[1] > 4000:
            logger.warning("Large image detected, computing full array...")
        img = data.compute()
    else:
        img = np.array(data)

    original_dtype = img.dtype

    # Estimate noise level before denoising
    try:
        from kintsugi.denoise.filters import estimate_noise_level
        noise_before = estimate_noise_level(img, method="mad")
    except ImportError:
        noise_before = float(np.std(img[img < np.percentile(img, 20)]))

    # Apply denoising based on method
    try:
        if method == "adaptive":
            from kintsugi.denoise.filters import adaptive_denoise
            result = adaptive_denoise(img, strength=strength, preserve_edges=preserve_edges)

        elif method == "bilateral":
            from kintsugi.denoise.filters import denoise_bilateral
            sigma_spatial = 1.0 if strength == "light" else (2.0 if strength == "medium" else 3.0)
            result = denoise_bilateral(img, sigma_spatial=sigma_spatial)

        elif method == "nlm":
            from kintsugi.denoise.filters import denoise_nlm
            patch_size = 5 if strength in ("light", "auto") else 7
            result = denoise_nlm(img, patch_size=patch_size)

        elif method == "bm3d":
            from kintsugi.denoise.patch_based import denoise_bm3d_lite
            result = denoise_bm3d_lite(img, sigma=noise_before if strength == "auto" else None)

        elif method == "n2v":
            from kintsugi.denoise.n2v import denoise_n2v
            result = denoise_n2v(
                img,
                n_epochs=n2v_epochs,
                model_path=model_path,
            )

        elif method == "patch_svd":
            from kintsugi.denoise.patch_based import denoise_svd_patch
            rank = None if strength == "auto" else (5 if strength == "light" else 10)
            result = denoise_svd_patch(img, rank=rank)

        else:
            return {"error": f"Unknown advanced denoise method: {method}"}

    except ImportError as e:
        return {
            "error": f"Required module not available: {e}",
            "hint": "Install with: pip install kintsugi[denoise] for N2V/CARE support",
        }

    # Estimate noise after
    try:
        noise_after = estimate_noise_level(result, method="mad")
    except Exception:
        noise_after = float(np.std(result[result < np.percentile(result, 20)]))

    # Store result
    output_name = f"{channel}_{method}_denoised"
    import dask.array as da
    result_dask = da.from_array(result, chunks=(1000, 1000))
    _store_image(output_name, result_dask, {
        "source": channel,
        "operation": f"denoise_advanced_{method}",
    })

    _add_history(output_name, "denoise_advanced", {
        "method": method,
        "strength": strength,
        "preserve_edges": preserve_edges,
        "noise_before": noise_before,
        "noise_after": noise_after,
    })

    return {
        "status": "success",
        "output_name": output_name,
        "method": method,
        "noise_reduction": {
            "before": round(noise_before, 2),
            "after": round(noise_after, 2),
            "reduction_pct": round((1 - noise_after / noise_before) * 100, 1) if noise_before > 0 else 0,
        },
        "parameters": {
            "strength": strength,
            "preserve_edges": preserve_edges,
        },
    }


# Utility functions for state management
def get_loaded_images() -> Dict[str, Any]:
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
