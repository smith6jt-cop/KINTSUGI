"""
Image chunking utilities for processing large images that don't fit in memory.

This module provides functions to automatically split large images into
manageable chunks, process them with overlap to avoid edge artifacts,
and reassemble the results.
"""

import numpy as np
from typing import Tuple, List, Callable, Optional
import warnings

# Try to import CuPy for GPU memory detection
try:
    import cupy as cp
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None


def _find_good_fft_length(n: int) -> int:
    """
    Find the next integer >= n whose largest prime factor is <= 5.

    FFT is most efficient for sizes that are products of small primes (2, 3, 5).
    """
    def max_prime_factor(num):
        if num <= 1:
            return 1
        factor = 2
        max_factor = 1
        while factor * factor <= num:
            while num % factor == 0:
                max_factor = factor
                num //= factor
            factor += 1
        if num > 1:
            max_factor = num
        return max_factor

    while max_prime_factor(n) > 5:
        n += 1
    return n


def estimate_padded_shape(shape: Tuple[int, ...], psf_shape: Tuple[int, ...]) -> Tuple[int, ...]:
    """
    Estimate the padded shape that deconvolution will use for FFT efficiency.

    Parameters
    ----------
    shape : tuple
        Original image shape
    psf_shape : tuple
        PSF shape

    Returns
    -------
    tuple
        Padded shape for FFT-efficient processing
    """
    return tuple(
        _find_good_fft_length(s + 4 * p)
        for s, p in zip(shape, psf_shape)
    )


def detect_best_device(device: str = 'auto') -> Tuple[bool, int, str]:
    """
    Detect the best device to use, prioritizing multi-GPU setups.

    This function provides unified device detection that should be used
    by both memory estimation and actual processing to avoid mismatches.

    Parameters
    ----------
    device : str
        'GPU', 'CPU', 'multi-gpu', or 'auto'

    Returns
    -------
    tuple
        (use_gpu: bool, device_id: int, device_name: str)
        device_id is -1 for CPU, otherwise the GPU device ID
    """
    if device.lower() == 'cpu':
        return False, -1, "CPU"

    # Try to use the centralized GPU manager first
    try:
        import sys
        from pathlib import Path

        # Add src to path if needed
        repo_root = Path(__file__).resolve().parents[2]
        src_path = repo_root / 'src'
        if src_path.exists() and str(src_path) not in sys.path:
            sys.path.insert(0, str(src_path))

        from kintsugi.gpu import get_gpu_manager, has_multi_gpu

        gpu = get_gpu_manager()

        if gpu.has_gpu and gpu.cupy_available:
            # Multi-GPU takes priority when 'auto' - select GPU with most free memory
            if has_multi_gpu() and device.lower() in ('auto', 'multi-gpu'):
                best_device = gpu.primary_device
                best_free = 0

                for info in gpu.gpu_info:
                    mem_info = gpu.get_memory_info(info.device_id)
                    if mem_info['free'] > best_free:
                        best_free = mem_info['free']
                        best_device = info.device_id

                # Set CuPy to use this device
                cp.cuda.Device(best_device).use()
                device_name = f"GPU {best_device} (Multi-GPU: {gpu.device_count} available)"
                return True, best_device, device_name

            # Single GPU or explicit GPU request
            device_id = gpu.primary_device
            device_name = gpu.gpu_info[0].name if gpu.gpu_info else "GPU"
            cp.cuda.Device(device_id).use()
            return True, device_id, f"GPU ({device_name})"

    except ImportError:
        pass
    except Exception as e:
        warnings.warn(f"GPUManager not available, using direct detection: {e}")

    # Fallback to direct CuPy detection
    if device.lower() in ('gpu', 'auto', 'multi-gpu') and CUPY_AVAILABLE:
        try:
            # Test that GPU is actually accessible, not just importable
            device_count = cp.cuda.runtime.getDeviceCount()
            if device_count > 0:
                # Get the device with most free memory
                best_device = 0
                best_free = 0

                for i in range(device_count):
                    with cp.cuda.Device(i):
                        free_mem, total_mem = cp.cuda.Device(i).mem_info
                        if free_mem > best_free:
                            best_free = free_mem
                            best_device = i

                cp.cuda.Device(best_device).use()
                props = cp.cuda.Device(best_device)

                if device_count > 1:
                    device_name = f"GPU {best_device} (Multi-GPU: {device_count} available)"
                else:
                    device_name = f"GPU ({props.name if hasattr(props, 'name') else best_device})"

                return True, best_device, device_name

        except Exception as e:
            if device.lower() == 'gpu':
                warnings.warn(f"GPU requested but not accessible: {e}. Falling back to CPU.")

    # CPU fallback
    return False, -1, "CPU"


def get_available_memory(device: str = 'auto', device_id: int = None) -> Tuple[int, str]:
    """
    Get available memory on the specified device.

    Uses unified device detection to ensure consistency between memory
    estimation and actual processing.

    Parameters
    ----------
    device : str
        'GPU', 'CPU', 'multi-gpu', or 'auto'
    device_id : int, optional
        Specific GPU device ID to query. If None, auto-detects best device.

    Returns
    -------
    tuple
        (available_bytes, device_name)
    """
    # Use unified device detection
    use_gpu, detected_device_id, device_name = detect_best_device(device)

    if device_id is not None:
        detected_device_id = device_id

    if use_gpu and CUPY_AVAILABLE:
        try:
            with cp.cuda.Device(detected_device_id):
                free_mem, total_mem = cp.cuda.Device(detected_device_id).mem_info
                return free_mem, device_name
        except Exception as e:
            warnings.warn(f"Could not get GPU memory: {e}. Falling back to CPU.")
            # Update device_name to reflect fallback
            device_name = "CPU (GPU memory detection failed)"

    # CPU memory
    try:
        import psutil
        mem = psutil.virtual_memory()
        return mem.available, device_name if "CPU" in device_name else "CPU"
    except ImportError:
        # Fallback: assume 8GB available
        warnings.warn("psutil not available, assuming 8GB memory")
        return 8 * 1024**3, "CPU (estimated)"


def calculate_chunks(shape: Tuple[int, int, int],
                     dtype: np.dtype,
                     available_memory: int,
                     memory_factor: float = 6.0,
                     overlap: int = 0) -> Tuple[int, int, int]:
    """
    Calculate the number of chunks needed in each dimension.

    The deconvolution process requires approximately 6x the input data size
    in memory (input, output, OTF, intermediate buffers).

    Parameters
    ----------
    shape : tuple
        Shape of the input image (z, y, x) or (x, y, z)
    dtype : dtype
        Data type
    available_memory : int
        Available memory in bytes
    memory_factor : float
        Memory multiplier for processing (default 6.0)
    overlap : int
        Overlap between chunks in pixels

    Returns
    -------
    tuple
        (nx_chunks, ny_chunks, nz_chunks)
    """
    bytes_per_voxel = np.dtype(dtype).itemsize
    if bytes_per_voxel < 4:
        bytes_per_voxel = 4  # Deconvolution uses float32

    total_voxels = np.prod(shape)
    total_bytes = total_voxels * bytes_per_voxel * memory_factor

    # Calculate maximum chunk size
    max_chunk_bytes = int(available_memory * 0.8)  # Use 80% of available memory
    max_chunk_voxels = max_chunk_bytes / (bytes_per_voxel * memory_factor)

    if total_bytes <= max_chunk_bytes:
        return (1, 1, 1)

    # Split along largest dimensions first (x, y are typically much larger than z)
    nx, ny, nz = 1, 1, 1

    # Iteratively increase splits until chunk fits in memory
    while True:
        chunk_shape = (
            int(np.ceil(shape[0] / nx)) + overlap,
            int(np.ceil(shape[1] / ny)) + overlap,
            int(np.ceil(shape[2] / nz)) + overlap
        )
        chunk_voxels = np.prod(chunk_shape)

        if chunk_voxels <= max_chunk_voxels:
            break

        # Split along the dimension with largest chunk size
        chunk_dims = [shape[0] / nx, shape[1] / ny, shape[2] / nz]

        if chunk_dims[0] >= chunk_dims[1] and chunk_dims[0] >= chunk_dims[2]:
            nx += 1
        elif chunk_dims[1] >= chunk_dims[2]:
            ny += 1
        else:
            nz += 1

        # Safety check
        if nx * ny * nz > 1000:
            raise RuntimeError("Calculated > 1000 chunks needed. Not enough memory.")

    return (nx, ny, nz)


def get_chunk_slices(shape: Tuple[int, int, int],
                     n_chunks: Tuple[int, int, int],
                     overlap: int = 0) -> List[Tuple[slice, slice, slice]]:
    """
    Calculate slice objects for each chunk with overlap.

    Parameters
    ----------
    shape : tuple
        Shape of the full image
    n_chunks : tuple
        Number of chunks in each dimension
    overlap : int
        Overlap between chunks

    Returns
    -------
    list
        List of (slice_x, slice_y, slice_z) for each chunk
    """
    nx, ny, nz = n_chunks
    slices = []

    chunk_sizes = [
        int(np.ceil(shape[0] / nx)),
        int(np.ceil(shape[1] / ny)),
        int(np.ceil(shape[2] / nz))
    ]

    for iz in range(nz):
        z_start = iz * chunk_sizes[2]
        z_end = min((iz + 1) * chunk_sizes[2], shape[2])

        # Add overlap (but don't exceed boundaries)
        z_start_ext = max(0, z_start - overlap // 2)
        z_end_ext = min(shape[2], z_end + overlap // 2)

        for iy in range(ny):
            y_start = iy * chunk_sizes[1]
            y_end = min((iy + 1) * chunk_sizes[1], shape[1])

            y_start_ext = max(0, y_start - overlap // 2)
            y_end_ext = min(shape[1], y_end + overlap // 2)

            for ix in range(nx):
                x_start = ix * chunk_sizes[0]
                x_end = min((ix + 1) * chunk_sizes[0], shape[0])

                x_start_ext = max(0, x_start - overlap // 2)
                x_end_ext = min(shape[0], x_end + overlap // 2)

                slices.append({
                    'extended': (
                        slice(x_start_ext, x_end_ext),
                        slice(y_start_ext, y_end_ext),
                        slice(z_start_ext, z_end_ext)
                    ),
                    'core': (
                        slice(x_start, x_end),
                        slice(y_start, y_end),
                        slice(z_start, z_end)
                    ),
                    'offset': (
                        x_start - x_start_ext,
                        y_start - y_start_ext,
                        z_start - z_start_ext
                    ),
                    'core_size': (
                        x_end - x_start,
                        y_end - y_start,
                        z_end - z_start
                    )
                })

    return slices


def process_in_chunks(stack: np.ndarray,
                      psf: np.ndarray,
                      deconv_func: Callable,
                      device: str = 'auto',
                      max_memory_fraction: float = 0.8,
                      overlap: Optional[int] = None,
                      verbose: bool = True) -> np.ndarray:
    """
    Process a large image stack in chunks.

    Parameters
    ----------
    stack : ndarray
        Input 3D image stack
    psf : ndarray
        Point spread function
    deconv_func : callable
        Function that takes (stack, psf) and returns deconvolved stack
    device : str
        'GPU', 'CPU', or 'auto'
    max_memory_fraction : float
        Fraction of available memory to use
    overlap : int, optional
        Overlap between chunks (default: max PSF dimension)
    verbose : bool
        Print progress information

    Returns
    -------
    ndarray
        Deconvolved image stack
    """
    if overlap is None:
        overlap = max(psf.shape)

    # Get available memory
    available_mem, device_name = get_available_memory(device)
    available_mem = int(available_mem * max_memory_fraction)

    if verbose:
        print(f"Available memory on {device_name}: {available_mem / 1e9:.2f} GB")

    # Calculate the PADDED shape that FFT deconvolution will actually use
    # This is critical for accurate memory estimation (fixes OOM bug)
    padded_shape = estimate_padded_shape(stack.shape, psf.shape)
    padding_ratio = np.prod(padded_shape) / np.prod(stack.shape)

    if verbose:
        print(f"FFT padding estimate: {stack.shape} → {padded_shape} ({padding_ratio:.2f}x voxels)")

    # Calculate number of chunks using PADDED shape for accurate memory estimation
    n_chunks = calculate_chunks(
        padded_shape,  # Use padded shape, not original!
        stack.dtype,
        available_mem,
        memory_factor=6.0,
        overlap=overlap
    )

    total_chunks = np.prod(n_chunks)

    if verbose:
        if total_chunks > 1:
            print(f"Splitting into {n_chunks[0]} x {n_chunks[1]} x {n_chunks[2]} = {total_chunks} chunks")
        else:
            print("Processing entire image (no chunking needed)")

    if total_chunks == 1:
        # No chunking needed
        return deconv_func(stack, psf)

    # Get chunk slices
    chunk_slices = get_chunk_slices(stack.shape, n_chunks, overlap)

    # Initialize output
    result = np.zeros(stack.shape, dtype=np.float32)

    # Process each chunk
    for i, chunk_info in enumerate(chunk_slices):
        if verbose:
            print(f"\nProcessing chunk {i + 1}/{total_chunks}")

        # Extract chunk with extended boundaries
        ext_slices = chunk_info['extended']
        chunk = stack[ext_slices[0], ext_slices[1], ext_slices[2]]

        # Process chunk
        chunk_result = deconv_func(chunk, psf)

        # Extract core region (without overlap)
        offset = chunk_info['offset']
        core_size = chunk_info['core_size']

        core_result = chunk_result[
            offset[0]:offset[0] + core_size[0],
            offset[1]:offset[1] + core_size[1],
            offset[2]:offset[2] + core_size[2]
        ]

        # Place in output
        core_slices = chunk_info['core']
        result[core_slices[0], core_slices[1], core_slices[2]] = core_result

        if verbose:
            print(f"  Chunk {i + 1} complete")

    return result


def estimate_memory_usage(shape: Tuple[int, int, int],
                          dtype: np.dtype = np.float32) -> dict:
    """
    Estimate memory usage for deconvolution.

    Parameters
    ----------
    shape : tuple
        Image shape
    dtype : dtype
        Data type

    Returns
    -------
    dict
        Dictionary with memory estimates in bytes
    """
    bytes_per_voxel = np.dtype(dtype).itemsize
    total_voxels = np.prod(shape)

    estimates = {
        'input': total_voxels * bytes_per_voxel,
        'output': total_voxels * bytes_per_voxel,
        'otf': total_voxels * 8,  # Complex64
        'buffers': total_voxels * bytes_per_voxel * 2,  # Intermediate buffers
        'total': total_voxels * (bytes_per_voxel * 4 + 8),
    }

    return estimates
