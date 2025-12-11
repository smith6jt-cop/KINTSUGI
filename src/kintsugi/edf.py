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

import logging
from typing import Literal, Optional, Tuple, Union

import numpy as np
from scipy.ndimage import gaussian_filter, uniform_filter

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


def _check_clij2_available() -> bool:
    """Check if PyImageJ/CLIJ2 is available.

    .. deprecated::
        CLIJ2 backend is deprecated. Use CuPy or NumPy backends instead.
    """
    try:
        import imagej
        import scyjava
        import warnings
        warnings.warn(
            "CLIJ2 backend is deprecated and will be removed in a future version. "
            "Use backend='cupy' or backend='numpy' instead.",
            DeprecationWarning,
            stacklevel=2
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
    """
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack, got shape {stack.shape}")

    original_dtype = stack.dtype
    n_slices, height, width = stack.shape

    # Select backend
    if use_gpu and _check_cupy_available():
        import cupy as cp
        from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian_filter
        from cupyx.scipy.ndimage import uniform_filter as gpu_uniform_filter

        xp = cp
        stack_gpu = cp.asarray(stack.astype(np.float32))
        uf = gpu_uniform_filter
        gf = gpu_gaussian_filter
        logger.info("Using CuPy GPU backend for EDF")
    else:
        xp = np
        stack_gpu = stack.astype(np.float32)
        uf = uniform_filter
        gf = gaussian_filter
        if use_gpu:
            logger.warning("GPU requested but CuPy not available, falling back to NumPy")
        else:
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

    # Find Z-index with maximum variance at each (x,y)
    max_var_idx = xp.argmax(variances, axis=0)

    # Create index arrays for advanced indexing
    y_idx, x_idx = xp.meshgrid(xp.arange(height), xp.arange(width), indexing="ij")

    # Select pixels from best-focus slices
    result = stack_gpu[max_var_idx, y_idx, x_idx]

    # Convert back to numpy if using GPU
    if use_gpu and _check_cupy_available():
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
        If True and CuPy available, use GPU acceleration

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

    if use_gpu and _check_cupy_available():
        import cupy as cp
        from cupyx.scipy.ndimage import convolve as gpu_convolve
        from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian_filter

        xp = cp
        stack_gpu = cp.asarray(stack.astype(np.float32))
        laplacian_kernel = cp.asarray(laplacian_kernel)
        gf = gpu_gaussian_filter
        conv = gpu_convolve
    else:
        from scipy.ndimage import convolve

        xp = np
        stack_gpu = stack.astype(np.float32)
        gf = gaussian_filter
        conv = convolve

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

    if use_gpu and _check_cupy_available():
        import cupy as cp

        result = cp.asnumpy(result)

    if np.issubdtype(original_dtype, np.integer):
        result = np.clip(result, np.iinfo(original_dtype).min, np.iinfo(original_dtype).max)

    return result.astype(original_dtype)


def edf_tiled(
    stack: np.ndarray,
    tile_size: Tuple[int, int] = (3000, 3000),
    overlap: int = 50,
    method: Literal["variance", "laplacian"] = "variance",
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
    ...     sigma=20.0
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

            # Process tile
            tile_result = edf_func(tile, **kwargs)

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
        z_start: Optional[int] = None,
        z_end: Optional[int] = None,
        tiles: Optional[Tuple[int, int]] = None,
        tile_size: Optional[Tuple[int, int]] = None,
        ij_instance=None,
        device: str = None,
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
            OpenCL device name for CLIJ2 backend

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
            use_gpu = self.backend == "cupy"

            if tile_size is not None:
                return edf_tiled(
                    stack,
                    tile_size=tile_size,
                    method=self.method,
                    radius_x=radius_x,
                    radius_y=radius_y,
                    sigma=sigma,
                    use_gpu=use_gpu,
                )
            else:
                if self.method == "variance":
                    return extended_depth_of_focus_variance(
                        stack, radius_x, radius_y, sigma, use_gpu
                    )
                else:
                    return extended_depth_of_focus_laplacian(stack, use_gpu=use_gpu)

    def _process_clij2(
        self,
        stack: np.ndarray,
        radius_x: int,
        radius_y: int,
        sigma: float,
        tiles: Tuple[int, int],
        ij_instance,
        device: str,
    ) -> np.ndarray:
        """Process using PyImageJ/CLIJ2 backend."""
        try:
            import imagej
            import scyjava
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
