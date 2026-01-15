"""
Extended Depth of Field (EDF) Processing Module

This module provides pure Python implementations of Extended Depth of Field
algorithms for multiplex immunofluorescence z-stack projection.

The primary algorithm is variance-based projection, matching CLIJ2's
`extendedDepthOfFocusVarianceProjection` function behavior.

Supports:
- GPU acceleration via CuPy (recommended, ~10x faster than CPU)
- CPU processing via NumPy/SciPy (fallback)
- Tiled processing for large images
- Unified interface with automatic backend selection

Note: The PyImageJ/CLIJ2 backend is deprecated and will be removed in a
future version. Use the CuPy or NumPy backends instead.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Literal, Union

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

if TYPE_CHECKING:
    import cupy as cp

logger = logging.getLogger(__name__)

# Type aliases
BackendType = Literal["auto", "numpy", "cupy", "clij2"]
ArrayLike = Union[np.ndarray, "cp.ndarray"]


def _check_cupy_available() -> bool:
    """Check if CuPy is available for GPU acceleration."""
    try:
        import cupy as cp

        # Test that CUDA is actually working
        _ = cp.array([1, 2, 3])
        return True
    except (ImportError, Exception):
        return False


def _detect_best_device(device: str = "auto") -> tuple[bool, int, str]:
    """
    Detect the best device to use, prioritizing multi-GPU setups.

    This provides unified device detection for consistency with deconvolution.

    Parameters
    ----------
    device : str
        'gpu', 'cpu', 'multi-gpu', or 'auto'

    Returns
    -------
    tuple
        (use_gpu: bool, device_id: int, device_name: str)
        device_id is -1 for CPU, otherwise the GPU device ID
    """
    if device.lower() == "cpu":
        return False, -1, "CPU"

    # Try to use the centralized GPU manager first
    try:
        from kintsugi.gpu import get_gpu_manager, has_multi_gpu

        gpu = get_gpu_manager()

        if gpu.has_gpu and gpu.cupy_available:
            import cupy as cp

            # Multi-GPU takes priority when 'auto' - select GPU with most free memory
            if has_multi_gpu() and device.lower() in ("auto", "multi-gpu"):
                best_device = gpu.primary_device
                best_free = 0

                for info in gpu.gpu_info:
                    mem_info = gpu.get_memory_info(info.device_id)
                    if mem_info["free"] > best_free:
                        best_free = mem_info["free"]
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
        logger.warning(f"GPUManager not available, using direct detection: {e}")

    # Fallback to direct CuPy detection
    if device.lower() in ("gpu", "auto", "multi-gpu") and _check_cupy_available():
        try:
            import cupy as cp

            # Test that GPU is actually accessible
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

                if device_count > 1:
                    device_name = f"GPU {best_device} (Multi-GPU: {device_count} available)"
                else:
                    device_name = f"GPU {best_device}"

                return True, best_device, device_name

        except Exception as e:
            if device.lower() == "gpu":
                logger.warning(f"GPU requested but not accessible: {e}. Falling back to CPU.")

    # CPU fallback
    return False, -1, "CPU"


def _check_clij2_available() -> bool:
    """Check if PyImageJ/CLIJ2 is available.

    .. deprecated::
        CLIJ2 backend is deprecated. Use CuPy or NumPy backends instead.
    """
    try:
        import warnings
        from importlib.util import find_spec

        if find_spec("imagej") is None or find_spec("scyjava") is None:
            return False
        warnings.warn(
            "CLIJ2 backend is deprecated and will be removed in a future version. "
            "Use backend='cupy' or backend='numpy' instead.",
            DeprecationWarning,
            stacklevel=2,
        )
        return True
    except ImportError:
        return False


def extended_depth_of_focus_variance(
    stack: np.ndarray,
    radius_x: int = 5,
    radius_y: int = 5,
    sigma: float = 20.0,
    use_gpu: bool = False,
    device: str = "auto",
    blend_depth: int = 0,
    z_smooth_sigma: float = 0.0,
) -> np.ndarray:
    """
    Extended depth of focus using local variance projection.

    This implementation matches CLIJ2's extendedDepthOfFocusVarianceProjection
    algorithm. For each (x,y) position, it selects the pixel from the Z-slice
    with the highest local variance (indicating best focus).

    Parameters
    ----------
    stack : np.ndarray
        3D image stack with shape (z, y, x). Can be any numeric dtype.
    radius_x : int, default=5
        Half-width of variance calculation window in X direction.
        Window width = 2 * radius_x + 1
    radius_y : int, default=5
        Half-height of variance calculation window in Y direction.
        Window height = 2 * radius_y + 1
    sigma : float, default=20.0
        Gaussian smoothing sigma applied before variance calculation.
        Larger values = more smoothing = less noise sensitivity.
    use_gpu : bool, default=False
        If True and CuPy is available, use GPU acceleration.
        Deprecated: Use `device` parameter instead.
    device : str, default='auto'
        Device selection: 'auto', 'gpu', 'cpu', or 'multi-gpu'.
        When 'auto', prioritizes multi-GPU and selects best available GPU.
    blend_depth : int, default=0
        Number of z-slices to blend for smooth transitions. When > 0, uses
        weighted blending of top N z-slices based on variance scores instead
        of hard argmax selection. This prevents abrupt transitions in areas
        of changing contrast. Recommended: 2-3 for smooth transitions.
    z_smooth_sigma : float, default=0.0
        Gaussian smoothing sigma for the z-index map. When > 0, smooths the
        selected z-indices spatially to reduce abrupt transitions between
        regions selecting different z-slices. Recommended: 2-5 pixels.

    Returns
    -------
    np.ndarray
        2D projected image with shape (y, x) and same dtype as input.

    Examples
    --------
    >>> import numpy as np
    >>> from kintsugi.edf import extended_depth_of_focus_variance
    >>> # Create a test stack (17 z-slices, 1000x1000 pixels)
    >>> stack = np.random.randint(0, 65535, (17, 1000, 1000), dtype=np.uint16)
    >>> result = extended_depth_of_focus_variance(stack, radius_x=5, radius_y=5, sigma=20.0)
    >>> result.shape
    (1000, 1000)

    Notes
    -----
    Algorithm:
    1. For each Z-slice, apply Gaussian smoothing with sigma
    2. Calculate local variance in a (2*radius_y+1) x (2*radius_x+1) window
    3. For each (x,y), find the Z-slice with maximum variance
    4. Output pixel = input pixel from that Z-slice

    The variance is calculated as: Var(X) = E[X^2] - E[X]^2

    For smooth transitions with blend_depth > 0:
    - Instead of selecting only the best z-slice, blends top N slices
    - Weights are based on variance scores (softmax normalization)
    - This eliminates abrupt boundaries in areas of varying contrast
    """
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack, got shape {stack.shape}")

    original_dtype = stack.dtype
    n_slices, height, width = stack.shape

    # Use unified device detection (prioritizes multi-GPU)
    # If use_gpu is explicitly set, convert to device string for backwards compatibility
    if use_gpu and device == "auto":
        device = "gpu"
    elif not use_gpu and device == "auto":
        # Keep as auto - let detection decide
        pass

    detected_gpu, device_id, device_name = _detect_best_device(device)

    # Select backend based on detected device
    if detected_gpu:
        import cupy as cp
        from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian_filter
        from cupyx.scipy.ndimage import uniform_filter as gpu_uniform_filter

        # Ensure we're using the selected device
        cp.cuda.Device(device_id).use()

        xp = cp
        stack_gpu = cp.asarray(stack.astype(np.float32))
        uf = gpu_uniform_filter
        gf = gpu_gaussian_filter
        logger.info(f"Using {device_name} for EDF")
    else:
        xp = np
        stack_gpu = stack.astype(np.float32)
        uf = uniform_filter
        gf = gaussian_filter
        logger.info("Using NumPy CPU backend for EDF")

    # Window size for variance calculation
    window_size = (2 * radius_y + 1, 2 * radius_x + 1)

    # Calculate variance for each slice
    variances = xp.zeros((n_slices, height, width), dtype=xp.float32)

    for z in range(n_slices):
        # Apply Gaussian smoothing
        smoothed = gf(stack_gpu[z], sigma=sigma)

        # Calculate local variance: Var(X) = E[X^2] - E[X]^2
        local_mean = uf(smoothed, size=window_size)
        local_mean_sq = uf(smoothed**2, size=window_size)
        variances[z] = local_mean_sq - local_mean**2

    # Determine blending strategy
    if blend_depth > 1 and blend_depth <= n_slices:
        # Weighted blending of top N z-slices for smooth transitions
        logger.info(f"Using weighted blending with depth={blend_depth} for smooth transitions")

        # Normalize variances to create weights (softmax-like)
        # Add small epsilon to avoid division by zero
        var_sum = xp.sum(variances, axis=0, keepdims=True)
        var_sum = xp.maximum(var_sum, 1e-10)
        weights = variances / var_sum

        # Get indices of top N z-slices by variance at each pixel
        # Use argsort and take the last N (highest variance)
        sorted_indices = xp.argsort(variances, axis=0)
        top_indices = sorted_indices[-blend_depth:]  # Shape: (blend_depth, H, W)

        # Initialize result
        result = xp.zeros((height, width), dtype=xp.float32)
        weight_sum = xp.zeros((height, width), dtype=xp.float32)

        # Create index arrays
        y_idx, x_idx = xp.meshgrid(xp.arange(height), xp.arange(width), indexing="ij")

        # Blend top N slices weighted by their variance
        for i in range(blend_depth):
            z_idx = top_indices[i]
            pixel_values = stack_gpu[z_idx, y_idx, x_idx]
            pixel_weights = weights[z_idx, y_idx, x_idx]
            result += pixel_values * pixel_weights
            weight_sum += pixel_weights

        # Normalize
        result = result / xp.maximum(weight_sum, 1e-10)

    else:
        # Standard argmax selection
        max_var_idx = xp.argmax(variances, axis=0)

        # Apply spatial smoothing to z-index map for gradual transitions
        if z_smooth_sigma > 0:
            logger.info(f"Applying z-index smoothing with sigma={z_smooth_sigma}")
            max_var_idx_float = max_var_idx.astype(xp.float32)
            max_var_idx_smooth = gf(max_var_idx_float, sigma=z_smooth_sigma)
            # Round to nearest integer and clip to valid range
            max_var_idx = xp.clip(xp.round(max_var_idx_smooth), 0, n_slices - 1).astype(xp.int32)

        # Create index arrays for advanced indexing
        y_idx, x_idx = xp.meshgrid(xp.arange(height), xp.arange(width), indexing="ij")

        # Select pixels from best-focus slices
        result = stack_gpu[max_var_idx, y_idx, x_idx]

    # Convert back to numpy if using GPU
    if detected_gpu:
        import cupy as cp

        result = cp.asnumpy(result)

    # Restore original dtype
    if np.issubdtype(original_dtype, np.integer):
        result = np.clip(result, np.iinfo(original_dtype).min, np.iinfo(original_dtype).max)
    result = result.astype(original_dtype)

    return result


def extended_depth_of_focus_laplacian(
    stack: np.ndarray,
    kernel_size: int = 5,
    use_gpu: bool = False,
    device: str = "auto",
) -> np.ndarray:
    """
    Extended depth of focus using Laplacian gradient method.

    An alternative to variance-based projection that uses the Laplacian
    (second derivative) as a focus measure. Pixels with higher absolute
    Laplacian values are considered more in-focus.

    Parameters
    ----------
    stack : np.ndarray
        3D image stack with shape (z, y, x)
    kernel_size : int, default=5
        Size of Gaussian blur kernel applied before Laplacian calculation
    use_gpu : bool, default=False
        If True and CuPy available, use GPU acceleration.
        Deprecated: Use `device` parameter instead.
    device : str, default='auto'
        Device selection: 'auto', 'gpu', 'cpu', or 'multi-gpu'.
        When 'auto', prioritizes multi-GPU and selects best available GPU.

    Returns
    -------
    np.ndarray
        2D projected image with shape (y, x)

    Notes
    -----
    This method is faster than variance-based but may be more sensitive to noise.
    """
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack, got shape {stack.shape}")

    original_dtype = stack.dtype
    n_slices, height, width = stack.shape

    # Laplacian kernel
    laplacian_kernel = np.array([[0, 1, 0], [1, -4, 1], [0, 1, 0]], dtype=np.float32)

    # Use unified device detection (prioritizes multi-GPU)
    if use_gpu and device == "auto":
        device = "gpu"

    detected_gpu, device_id, device_name = _detect_best_device(device)

    if detected_gpu:
        import cupy as cp
        from cupyx.scipy.ndimage import convolve as gpu_convolve
        from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian_filter

        # Ensure we're using the selected device
        cp.cuda.Device(device_id).use()

        xp = cp
        stack_gpu = cp.asarray(stack.astype(np.float32))
        laplacian_kernel = cp.asarray(laplacian_kernel)
        gf = gpu_gaussian_filter
        conv = gpu_convolve
        logger.info(f"Using {device_name} for EDF (Laplacian)")
    else:
        from scipy.ndimage import convolve

        xp = np
        stack_gpu = stack.astype(np.float32)
        gf = gaussian_filter
        conv = convolve
        logger.info("Using NumPy CPU backend for EDF (Laplacian)")

    # Calculate Laplacian magnitude for each slice
    focus_measures = xp.zeros((n_slices, height, width), dtype=xp.float32)

    for z in range(n_slices):
        # Apply Gaussian blur to reduce noise
        blurred = gf(stack_gpu[z], sigma=kernel_size / 3)

        # Calculate Laplacian
        laplacian = conv(blurred, laplacian_kernel)

        # Use absolute value as focus measure
        focus_measures[z] = xp.abs(laplacian)

    # Find Z-index with maximum focus measure
    max_focus_idx = xp.argmax(focus_measures, axis=0)

    # Select pixels
    y_idx, x_idx = xp.meshgrid(xp.arange(height), xp.arange(width), indexing="ij")
    result = stack_gpu[max_focus_idx, y_idx, x_idx]

    if detected_gpu:
        import cupy as cp

        result = cp.asnumpy(result)

    if np.issubdtype(original_dtype, np.integer):
        result = np.clip(result, np.iinfo(original_dtype).min, np.iinfo(original_dtype).max)

    return result.astype(original_dtype)


def edf_tiled(
    stack: np.ndarray,
    tile_size: tuple[int, int] = (3000, 3000),
    overlap: int = 50,
    method: Literal["variance", "laplacian"] = "variance",
    device: str = "auto",
    blend_depth: int = 0,
    z_smooth_sigma: float = 0.0,
    **kwargs,
) -> np.ndarray:
    """
    Process large stacks in tiles to manage memory.

    Splits the image into overlapping tiles, processes each independently,
    and blends them together using raised cosine weights for seamless output.

    Parameters
    ----------
    stack : np.ndarray
        3D image stack with shape (z, y, x)
    tile_size : tuple of int, default=(3000, 3000)
        (height, width) of each processing tile
    overlap : int, default=50
        Pixel overlap between tiles for seamless blending
    method : {'variance', 'laplacian'}, default='variance'
        EDF algorithm to use
    device : str, default='auto'
        Device selection: 'auto', 'gpu', 'cpu', or 'multi-gpu'.
        When 'auto', prioritizes multi-GPU and selects best available GPU.
    blend_depth : int, default=0
        Number of z-slices to blend for smooth transitions. When > 0, uses
        weighted blending of top N z-slices based on variance scores instead
        of hard argmax selection. This prevents abrupt transitions in areas
        of changing contrast. Recommended: 2-3 for smooth transitions.
    z_smooth_sigma : float, default=0.0
        Gaussian smoothing sigma for the z-index map. When > 0, smooths the
        selected z-indices spatially to reduce abrupt transitions between
        regions selecting different z-slices. Recommended: 2-5 pixels.
    **kwargs
        Additional arguments passed to the EDF function

    Returns
    -------
    np.ndarray
        2D EDF result with same dtype as input

    Examples
    --------
    >>> # Process a large stack in 3x3 tiles (similar to CLIJ2 xTiles=3, yTiles=3)
    >>> from kintsugi.edf import edf_tiled
    >>> result = edf_tiled(
    ...     large_stack,
    ...     tile_size=(3000, 3000),
    ...     overlap=100,
    ...     method='variance',
    ...     radius_x=5,
    ...     radius_y=5,
    ...     sigma=20.0,
    ...     blend_depth=2  # Smooth transitions
    ... )
    """
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack, got shape {stack.shape}")

    n_slices, height, width = stack.shape
    tile_h, tile_w = tile_size
    original_dtype = stack.dtype

    # Select EDF function
    if method == "variance":
        edf_func = extended_depth_of_focus_variance
    elif method == "laplacian":
        edf_func = extended_depth_of_focus_laplacian
    else:
        raise ValueError(f"Unknown method: {method}")

    # Initialize output arrays
    result = np.zeros((height, width), dtype=np.float64)
    weight = np.zeros((height, width), dtype=np.float64)

    # Create blending weights (raised cosine for smooth transitions)
    def create_blend_weights(h, w, overlap):
        blend_y = np.ones(h)
        blend_x = np.ones(w)
        if overlap > 0:
            ramp = 0.5 * (1 - np.cos(np.pi * np.arange(overlap) / overlap))
            if h > 2 * overlap:
                blend_y[:overlap] = ramp
                blend_y[-overlap:] = ramp[::-1]
            if w > 2 * overlap:
                blend_x[:overlap] = ramp
                blend_x[-overlap:] = ramp[::-1]
        return np.outer(blend_y, blend_x)

    # Calculate tile positions
    step_h = max(1, tile_h - overlap)
    step_w = max(1, tile_w - overlap)

    n_tiles_y = (height + step_h - 1) // step_h
    n_tiles_x = (width + step_w - 1) // step_w

    logger.info(f"Processing {n_tiles_y}x{n_tiles_x} tiles of size {tile_h}x{tile_w}")

    tile_count = 0
    total_tiles = n_tiles_y * n_tiles_x

    for y_start in range(0, height, step_h):
        for x_start in range(0, width, step_w):
            # Calculate tile bounds
            y_end = min(y_start + tile_h, height)
            x_end = min(x_start + tile_w, width)

            # Extract tile
            tile = stack[:, y_start:y_end, x_start:x_end]

            # Process tile (pass device and smooth transition parameters explicitly)
            if method == "variance":
                tile_result = edf_func(
                    tile, device=device, blend_depth=blend_depth,
                    z_smooth_sigma=z_smooth_sigma, **kwargs
                )
            else:
                tile_result = edf_func(tile, device=device, **kwargs)

            # Create blending weights for this tile
            tile_blend = create_blend_weights(y_end - y_start, x_end - x_start, overlap)

            # Accumulate weighted result
            result[y_start:y_end, x_start:x_end] += tile_result.astype(np.float64) * tile_blend
            weight[y_start:y_end, x_start:x_end] += tile_blend

            tile_count += 1
            if tile_count % 10 == 0 or tile_count == total_tiles:
                logger.info(f"Processed {tile_count}/{total_tiles} tiles")

    # Normalize by weights
    result = result / np.maximum(weight, 1e-10)

    # Restore dtype
    if np.issubdtype(original_dtype, np.integer):
        result = np.clip(result, np.iinfo(original_dtype).min, np.iinfo(original_dtype).max)

    return result.astype(original_dtype)


class EDFProcessor:
    """
    Unified Extended Depth of Field processor with automatic backend selection.

    This class provides a consistent interface for EDF processing, automatically
    selecting the best available backend (CLIJ2 > CuPy > NumPy) based on
    availability and user preference.

    Parameters
    ----------
    backend : {'auto', 'numpy', 'cupy', 'clij2'}, default='auto'
        Processing backend to use:
        - 'auto': Automatically select best available backend
        - 'numpy': Force NumPy/SciPy CPU processing
        - 'cupy': Force CuPy GPU processing
        - 'clij2': Force PyImageJ/CLIJ2 processing
    method : {'variance', 'laplacian'}, default='variance'
        EDF algorithm to use (only for Python backends)

    Attributes
    ----------
    backend : str
        The actual backend being used
    method : str
        The EDF algorithm being used

    Examples
    --------
    >>> from kintsugi.edf import EDFProcessor
    >>> processor = EDFProcessor(backend='auto')
    >>> result = processor.process(
    ...     stack,
    ...     radius_x=5,
    ...     radius_y=5,
    ...     sigma=20.0,
    ...     tiles=(3, 3)
    ... )
    """

    def __init__(
        self,
        backend: BackendType = "auto",
        method: Literal["variance", "laplacian"] = "variance",
    ):
        self.method = method
        self._requested_backend = backend
        self.backend = self._detect_backend(backend)
        logger.info(f"EDFProcessor initialized with backend: {self.backend}")

    def _detect_backend(self, requested: BackendType) -> str:
        """Detect and validate the processing backend."""
        if requested == "auto":
            # Priority: CuPy > NumPy (CLIJ2 deprecated)
            if _check_cupy_available():
                return "cupy"
            else:
                return "numpy"

        elif requested == "clij2":
            # CLIJ2 is deprecated but still supported for backwards compatibility
            logger.warning(
                "CLIJ2 backend is deprecated and will be removed in a future version. "
                "Consider using backend='cupy' or backend='auto' instead."
            )
            if not _check_clij2_available():
                logger.warning("CLIJ2 not available, falling back to auto-detection")
                return self._detect_backend("auto")
            return "clij2"

        elif requested == "cupy":
            if not _check_cupy_available():
                logger.warning("CuPy not available, falling back to NumPy")
                return "numpy"
            return "cupy"

        elif requested == "numpy":
            return "numpy"

        else:
            raise ValueError(f"Unknown backend: {requested}")

    def process(
        self,
        stack: np.ndarray,
        radius_x: int = 5,
        radius_y: int = 5,
        sigma: float = 20.0,
        z_start: int | None = None,
        z_end: int | None = None,
        tiles: tuple[int, int] | None = None,
        tile_size: tuple[int, int] | None = None,
        ij_instance=None,
        device: str = None,
        blend_depth: int = 0,
        z_smooth_sigma: float = 0.0,
    ) -> np.ndarray:
        """
        Process a 3D stack to produce a 2D EDF image.

        Parameters
        ----------
        stack : np.ndarray
            3D image stack with shape (z, y, x)
        radius_x : int, default=5
            Variance window half-width in X
        radius_y : int, default=5
            Variance window half-height in Y
        sigma : float, default=20.0
            Gaussian smoothing sigma
        z_start : int, optional
            Starting Z-slice (1-indexed for CLIJ2 compatibility)
        z_end : int, optional
            Ending Z-slice (1-indexed for CLIJ2 compatibility)
        tiles : tuple of int, optional
            Number of tiles (y_tiles, x_tiles) for CLIJ2-style processing.
            Automatically calculates tile_size from image dimensions.
        tile_size : tuple of int, optional
            Explicit tile size (height, width). Overrides tiles parameter.
        ij_instance : optional
            Pre-initialized ImageJ instance for CLIJ2 backend
        device : str, optional
            Device selection: 'auto', 'gpu', 'cpu', or 'multi-gpu'.
            For CLIJ2 backend, this is the OpenCL device name.
            For CuPy/NumPy backends, 'auto' prioritizes multi-GPU.
        blend_depth : int, default=0
            Number of z-slices to blend for smooth transitions. When > 0, uses
            weighted blending of top N z-slices based on variance scores.
            This prevents abrupt transitions in areas of changing contrast.
            Recommended: 2-3 for smooth transitions.
        z_smooth_sigma : float, default=0.0
            Gaussian smoothing sigma for the z-index map. When > 0, smooths the
            selected z-indices spatially to reduce abrupt transitions.
            Recommended: 2-5 pixels.

        Returns
        -------
        np.ndarray
            2D EDF result
        """
        # Slice the stack if z_start/z_end specified
        if z_start is not None or z_end is not None:
            # Convert to 0-indexed
            z_start_idx = (z_start - 1) if z_start is not None else 0
            z_end_idx = z_end if z_end is not None else stack.shape[0]
            stack = stack[z_start_idx:z_end_idx]
            logger.info(f"Using Z-slices {z_start_idx+1} to {z_end_idx} ({stack.shape[0]} slices)")

        n_slices, height, width = stack.shape

        # Determine tile size
        if tile_size is None and tiles is not None:
            y_tiles, x_tiles = tiles
            tile_size = (height // y_tiles + 1, width // x_tiles + 1)
            logger.info(f"Using {y_tiles}x{x_tiles} tiling: tile_size={tile_size}")

        # Route to appropriate backend
        if self.backend == "clij2":
            return self._process_clij2(
                stack, radius_x, radius_y, sigma, tiles or (1, 1), ij_instance, device
            )
        else:
            # Determine device for CuPy/NumPy backends
            # Use 'auto' for cupy backend, 'cpu' for numpy backend
            if device is None:
                device = "auto" if self.backend == "cupy" else "cpu"

            if tile_size is not None:
                return edf_tiled(
                    stack,
                    tile_size=tile_size,
                    method=self.method,
                    device=device,
                    blend_depth=blend_depth,
                    z_smooth_sigma=z_smooth_sigma,
                    radius_x=radius_x,
                    radius_y=radius_y,
                    sigma=sigma,
                )
            else:
                if self.method == "variance":
                    return extended_depth_of_focus_variance(
                        stack, radius_x, radius_y, sigma, device=device,
                        blend_depth=blend_depth, z_smooth_sigma=z_smooth_sigma
                    )
                else:
                    return extended_depth_of_focus_laplacian(stack, device=device)

    def _process_clij2(
        self,
        stack: np.ndarray,
        radius_x: int,
        radius_y: int,
        sigma: float,
        tiles: tuple[int, int],
        ij_instance,
        device: str,
    ) -> np.ndarray:
        """Process using PyImageJ/CLIJ2 backend."""
        try:
            import imagej  # noqa: F401
        except ImportError:
            raise RuntimeError("PyImageJ not available for CLIJ2 backend")

        # Initialize ImageJ if not provided
        if ij_instance is None:
            logger.warning(
                "No ImageJ instance provided. For batch processing, "
                "pass a pre-initialized instance to avoid startup overhead."
            )
            ij_instance = imagej.init()

        n_slices, height, width = stack.shape
        y_tiles, x_tiles = tiles

        # Build macro
        macro = f"""
        run("CLIJ2 Macro Extensions", "cl_device=[{device or ''}] automatic_output_naming=false");
        Ext.CLIJ2_clear();

        input = "input_stack";
        output = "output";
        newImage(output, "16-bit black", {width}, {height}, 1);

        numTilesZ = 1;
        numTilesX = {x_tiles};
        numTilesY = {y_tiles};

        tileDepth = round(({n_slices}/numTilesZ));
        tileWidth = round(({width}/numTilesX));
        tileHeight = round(({height}/numTilesY));

        for (x = 0; x < numTilesX; x++) {{
            for (y = 0; y < numTilesY; y++) {{
                for (z = 0; z < numTilesZ; z++) {{
                    Ext.CLIJ2_pushTile(input, x, y, z, tileWidth, tileHeight, tileDepth, 0, 0, 0);
                    Ext.CLIJ2_extendedDepthOfFocusVarianceProjection(input, output, {radius_x}, {radius_y}, {sigma});
                    Ext.CLIJ2_release(input);
                    Ext.CLIJ2_pullTile(output, x, y, z, tileWidth, tileHeight, 1, 0, 0, 0);
                    Ext.CLIJ2_release(output);
                }}
            }}
        }}

        selectImage(output);
        Ext.CLIJ2_clear();
        """

        # Convert stack to ImageJ format and run
        ij_stack = ij_instance.py.to_java(stack)
        ij_instance.ui().show("input_stack", ij_stack)

        ij_instance.py.run_macro(macro)

        # Get result
        result_ij = ij_instance.py.active_dataset()
        result = ij_instance.py.from_java(result_ij)

        return np.array(result)


# Convenience functions for common use cases


def process_edf(
    stack: np.ndarray,
    radius_x: int = 5,
    radius_y: int = 5,
    sigma: float = 20.0,
    backend: BackendType = "auto",
    **kwargs,
) -> np.ndarray:
    """
    Convenience function for one-shot EDF processing.

    Parameters
    ----------
    stack : np.ndarray
        3D image stack (z, y, x)
    radius_x, radius_y : int
        Variance window radii
    sigma : float
        Gaussian smoothing sigma
    backend : str
        Processing backend ('auto', 'numpy', 'cupy', 'clij2')
    **kwargs
        Additional arguments passed to EDFProcessor.process()

    Returns
    -------
    np.ndarray
        2D EDF result

    Examples
    --------
    >>> from kintsugi.edf import process_edf
    >>> result = process_edf(stack, radius_x=5, radius_y=5, sigma=20.0)
    """
    processor = EDFProcessor(backend=backend)
    return processor.process(stack, radius_x, radius_y, sigma, **kwargs)
