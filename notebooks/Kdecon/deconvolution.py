"""
Lucy-Richardson deconvolution implementation with GPU and CPU support.

This module implements the Lucy-Richardson algorithm for image deconvolution
using FFT-based convolution. It supports both GPU acceleration via CuPy
and multi-threaded CPU processing via NumPy/SciPy.
"""

import numpy as np
from scipy import fft as scipy_fft
from scipy.ndimage import convolve
from typing import Tuple, Optional, Union, Callable
import warnings

# Try to import CuPy for GPU support
try:
    import cupy as cp
    from cupyx.scipy import fft as cupy_fft
    from cupyx.scipy.ndimage import convolve as cupy_convolve
    CUPY_AVAILABLE = True
except ImportError:
    CUPY_AVAILABLE = False
    cp = None


def _find_good_fft_length(n: int) -> int:
    """
    Find the next integer >= n whose largest prime factor is <= 5.

    FFT is most efficient for sizes that are products of small primes (2, 3, 5).

    Parameters
    ----------
    n : int
        Minimum size needed

    Returns
    -------
    int
        Efficient FFT size >= n
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


def _psf_to_otf(psf: np.ndarray, shape: Tuple[int, ...]) -> np.ndarray:
    """
    Convert PSF to OTF (Optical Transfer Function) via FFT.

    Parameters
    ----------
    psf : ndarray
        Point spread function (must be smaller than shape in all dimensions)
    shape : tuple
        Output shape for OTF

    Returns
    -------
    ndarray
        Complex OTF array
    """
    # Pad PSF to output shape
    pad_shape = [(0, s - p) for s, p in zip(shape, psf.shape)]
    psf_padded = np.pad(psf, pad_shape, mode='constant', constant_values=0)

    # Circularly shift so center of PSF is at origin
    for axis in range(psf.ndim):
        psf_padded = np.roll(psf_padded, -(psf.shape[axis] // 2), axis=axis)

    # FFT to get OTF
    return scipy_fft.fftn(psf_padded)


def _psf_to_otf_gpu(psf: 'cp.ndarray', shape: Tuple[int, ...]) -> 'cp.ndarray':
    """GPU version of _psf_to_otf."""
    pad_shape = [(0, s - p) for s, p in zip(shape, psf.shape)]
    psf_padded = cp.pad(psf, pad_shape, mode='constant', constant_values=0)

    for axis in range(psf.ndim):
        psf_padded = cp.roll(psf_padded, -(psf.shape[axis] // 2), axis=axis)

    return cupy_fft.fftn(psf_padded)


def _conv_fft(data: np.ndarray, otf: np.ndarray) -> np.ndarray:
    """Convolution via FFT multiplication (CPU)."""
    return np.real(scipy_fft.ifftn(otf * scipy_fft.fftn(data)))


def _conv_fft_gpu(data: 'cp.ndarray', otf: 'cp.ndarray') -> 'cp.ndarray':
    """Convolution via FFT multiplication (GPU)."""
    return cp.real(cupy_fft.ifftn(otf * cupy_fft.fftn(data)))


def _conv_same_gpu(data: 'cp.ndarray', kernel: 'cp.ndarray') -> 'cp.ndarray':
    """
    3D convolution with 'same' output size on GPU.

    Equivalent to MATLAB's convn(data, kernel, 'same').
    Uses FFT with zero-padding to avoid circular wraparound artifacts.

    Parameters
    ----------
    data : cp.ndarray
        Input 3D array
    kernel : cp.ndarray
        Convolution kernel (PSF)

    Returns
    -------
    cp.ndarray
        Convolved array with same shape as input data
    """
    # Calculate padded size to avoid circular convolution artifacts
    # Full convolution size is data.shape + kernel.shape - 1
    full_shape = tuple(d + k - 1 for d, k in zip(data.shape, kernel.shape))

    # Find FFT-efficient sizes (optional but faster)
    def next_fast_len(n):
        # Find next size with small prime factors for efficient FFT
        while True:
            m = n
            for p in (2, 3, 5):
                while m % p == 0:
                    m //= p
            if m == 1:
                return n
            n += 1

    fft_shape = tuple(next_fast_len(s) for s in full_shape)

    # FFT of zero-padded data and kernel
    data_fft = cupy_fft.fftn(data, s=fft_shape)
    kernel_fft = cupy_fft.fftn(kernel, s=fft_shape)

    # Multiply in frequency domain and inverse FFT
    result_full = cp.real(cupy_fft.ifftn(data_fft * kernel_fft))

    # Extract 'same' size output (centered on full convolution)
    # For 'same' mode: output[i] corresponds to kernel centered at data[i]
    # Start index is (kernel.shape - 1) // 2
    start = tuple((k - 1) // 2 for k in kernel.shape)
    end = tuple(s + d for s, d in zip(start, data.shape))

    return result_full[start[0]:end[0], start[1]:end[1], start[2]:end[2]]


def _create_damping_kernel() -> np.ndarray:
    """Create 3x3x3 averaging kernel for damping (excluding center)."""
    R = np.ones((3, 3, 3), dtype=np.float32) / 26
    R[1, 1, 1] = 0
    return R


def _create_tukey_window_3d(shape: Tuple[int, ...],
                            pad_before: Tuple[int, ...],
                            pad_after: Tuple[int, ...],
                            alpha: float = 0.5) -> np.ndarray:
    """
    Create a 3D Tukey (cosine-tapered) window for apodization.

    This window smoothly tapers the padded regions to prevent FFT edge
    artifacts (Gibbs ringing) that cause horizontal banding.

    Parameters
    ----------
    shape : tuple
        Shape of the padded array (padded_x, padded_y, padded_z)
    pad_before : tuple
        Padding added before data in each dimension
    pad_after : tuple
        Padding added after data in each dimension
    alpha : float
        Shape parameter (0 = rectangular, 1 = Hann window).
        Default 0.5 provides good edge smoothing.

    Returns
    -------
    ndarray
        3D window array with smooth transitions at boundaries
    """
    windows_1d = []

    for i in range(3):
        n = shape[i]
        pb = pad_before[i]
        pa = pad_after[i]

        if pb == 0 and pa == 0:
            # No padding in this dimension - use flat window
            windows_1d.append(np.ones(n, dtype=np.float32))
            continue

        w = np.ones(n, dtype=np.float32)

        # Taper width is the padding region (or at least 4 pixels)
        taper_before = max(pb, 4) if pb > 0 else 0
        taper_after = max(pa, 4) if pa > 0 else 0

        # Cosine taper at the start
        if taper_before > 0:
            t = np.linspace(0, np.pi/2, taper_before)
            w[:taper_before] = np.sin(t) ** 2

        # Cosine taper at the end
        if taper_after > 0:
            t = np.linspace(np.pi/2, 0, taper_after)
            w[-taper_after:] = np.sin(t) ** 2

        windows_1d.append(w)

    # Create 3D window as outer product of 1D windows
    window = np.outer(windows_1d[0], windows_1d[1]).reshape(shape[0], shape[1], 1)
    window = window * windows_1d[2].reshape(1, 1, shape[2])

    return window.astype(np.float32)


def lucy_richardson_cpu(stack: np.ndarray, psf: np.ndarray,
                        iterations: int = 25,
                        damping: float = 0.0,
                        stop_criterion: float = 5.0,
                        callback: Optional[Callable] = None,
                        verbose: bool = True) -> np.ndarray:
    """
    Lucy-Richardson deconvolution on CPU.

    Parameters
    ----------
    stack : ndarray
        Input 3D image stack (float32)
    psf : ndarray
        Point spread function
    iterations : int
        Maximum number of iterations
    damping : float
        Damping factor (0-1, typically 0-0.1). Increase for noisy images.
    stop_criterion : float
        Stop if relative change falls below this percentage
    callback : callable, optional
        Function called after each iteration with (iteration, delta_rel)
    verbose : bool
        Print progress information

    Returns
    -------
    ndarray
        Deconvolved image stack
    """
    # Initialize deconvolved estimate
    deconvolved = stack.astype(np.float32).copy()

    # Compute OTF
    otf = _psf_to_otf(psf.astype(np.float32), stack.shape)
    otf_conj = np.conj(otf)

    # Damping kernel
    if damping > 0:
        R = _create_damping_kernel()

    delta_last = np.inf
    # Use eps('single') as floor value - matches MATLAB exactly
    eps = np.finfo(np.float32).eps  # ~1.19e-7

    for i in range(iterations):
        # Richardson-Lucy iteration:
        # d_new = d * conj(OTF) * (stack / (OTF * d))

        # Forward model: blur estimate
        denom = _conv_fft(deconvolved, otf)
        denom = np.maximum(denom, eps)  # Protect against division by zero

        # Ratio (no clipping - matches MATLAB)
        ratio = stack / denom

        # Backward step
        if damping == 0:
            deconvolved_new = _conv_fft(ratio, otf_conj) * deconvolved
        else:
            # With damping regularization
            lr_update = _conv_fft(ratio, otf_conj) * deconvolved
            smooth_term = convolve(deconvolved, R, mode='reflect')
            deconvolved_new = (1 - damping) * lr_update + damping * smooth_term

        # Calculate convergence criterion
        delta = np.sqrt(np.sum((deconvolved - deconvolved_new)**2))

        if i == 0:
            delta_rel = 0.0
        else:
            delta_rel = (delta_last - delta) / delta_last * 100 if delta_last > 0 else 0.0

        deconvolved = deconvolved_new
        delta_last = delta

        if verbose:
            print(f"  Iteration {i+1}: delta_rel = {delta_rel:.2f}%")

        if callback:
            callback(i + 1, delta_rel)

        # Check stop criterion
        if i > 0 and delta_rel <= stop_criterion:
            if verbose:
                print(f"  Stop criterion reached at iteration {i+1}")
            break

    # Remove any imaginary artifacts
    return np.abs(deconvolved)


def lucy_richardson_gpu(stack: np.ndarray, psf: np.ndarray,
                        iterations: int = 25,
                        damping: float = 0.0,
                        stop_criterion: float = 5.0,
                        callback: Optional[Callable] = None,
                        verbose: bool = True) -> np.ndarray:
    """
    Lucy-Richardson deconvolution on GPU using CuPy.

    This implementation matches MATLAB's deconGPU exactly:
    - Uses spatial convolution (not FFT) with 'same' output size
    - Uses spatially reversed PSF for backward step
    - Uses eps('single') as floor value

    Parameters are the same as lucy_richardson_cpu().
    """
    if not CUPY_AVAILABLE:
        raise RuntimeError("CuPy is not available. Install cupy for GPU support.")

    # Transfer to GPU
    deconvolved = cp.asarray(stack.astype(np.float32))
    stack_gpu = cp.asarray(stack.astype(np.float32))
    psf_gpu = cp.asarray(psf.astype(np.float32))

    # Spatially reversed PSF for backward step (matches MATLAB: psf(end:-1:1, end:-1:1, end:-1:1))
    psf_inv = psf_gpu[::-1, ::-1, ::-1]

    # Damping kernel
    if damping > 0:
        R = cp.asarray(_create_damping_kernel())

    delta_last = np.inf
    # Use eps('single') as floor value - matches MATLAB exactly
    eps = np.finfo(np.float32).eps  # ~1.19e-7

    try:
        for i in range(iterations):
            # Forward model: convolution with PSF (matches MATLAB convGPU)
            # MATLAB: denom = convGPU(deconvolved, psf) = gather(convn(gpuArray(data), psf, 'same'))
            # Uses _conv_same_gpu which implements convn(..., 'same') via FFT with proper padding
            denom = _conv_same_gpu(deconvolved, psf_gpu)
            denom = cp.maximum(denom, eps)  # Protect against division by zero

            # Ratio (no clipping - matches MATLAB)
            ratio = stack_gpu / denom

            # Backward step: convolution with reversed PSF
            # MATLAB: convGPU(stack ./ denom, psf_inv) .* deconvolved
            if damping == 0:
                deconvolved_new = _conv_same_gpu(ratio, psf_inv) * deconvolved
            else:
                lr_update = _conv_same_gpu(ratio, psf_inv) * deconvolved
                smooth_term = cupy_convolve(deconvolved, R, mode='reflect')
                deconvolved_new = (1 - damping) * lr_update + damping * smooth_term

            # Calculate convergence
            delta = float(cp.sqrt(cp.sum((deconvolved - deconvolved_new)**2)))

            if i == 0:
                delta_rel = 0.0
            else:
                delta_rel = (delta_last - delta) / delta_last * 100 if delta_last > 0 else 0.0

            deconvolved = deconvolved_new
            delta_last = delta

            if verbose:
                print(f"  Iteration {i+1}: delta_rel = {delta_rel:.2f}%")

            if callback:
                callback(i + 1, delta_rel)

            if i > 0 and delta_rel <= stop_criterion:
                if verbose:
                    print(f"  Stop criterion reached at iteration {i+1}")
                break

        # Transfer back to CPU - get rid of imaginary artifacts (matches MATLAB: deconvolved = abs(deconvolved))
        result = cp.asnumpy(cp.abs(deconvolved))

    finally:
        # Clean up GPU memory
        del deconvolved, stack_gpu, psf_gpu, psf_inv
        if damping > 0:
            del R
        cp.get_default_memory_pool().free_all_blocks()

    return result


def deconvolve(stack: np.ndarray, psf: np.ndarray,
               iterations: int = 25,
               damping: float = 0.0,
               stop_criterion: float = 5.0,
               device: str = 'auto',
               device_id: Optional[int] = None,
               pad_for_fft: bool = True,
               callback: Optional[Callable] = None,
               verbose: bool = True) -> np.ndarray:
    """
    Deconvolve a 3D image stack using Lucy-Richardson algorithm.

    This function automatically handles FFT-efficient padding and
    selects GPU or CPU based on availability and device parameter.

    Parameters
    ----------
    stack : ndarray
        Input 3D image stack
    psf : ndarray
        Point spread function
    iterations : int
        Maximum number of iterations (default 25)
    damping : float
        Damping factor for regularization (0-1, default 0).
        Increase for noisy images.
    stop_criterion : float
        Stop if relative change falls below this percentage (default 5.0)
    device : str
        'GPU', 'CPU', or 'auto' (default 'auto')
    device_id : int, optional
        Specific GPU device ID for multi-GPU processing. If provided,
        this GPU will be used instead of auto-detection.
    pad_for_fft : bool
        Pad data for efficient FFT (default True)
    callback : callable, optional
        Called after each iteration with (iteration, delta_rel)
    verbose : bool
        Print progress (default True)

    Returns
    -------
    ndarray
        Deconvolved image stack (same shape as input)
    """
    # Validate inputs
    if stack.ndim != 3:
        raise ValueError(f"Expected 3D stack, got {stack.ndim}D")
    if psf.ndim != 3:
        raise ValueError(f"Expected 3D PSF, got {psf.ndim}D")

    # Use unified device detection from chunking module for consistency
    # This ensures memory estimation and actual processing use the same device
    from .chunking import detect_best_device

    use_gpu, detected_device_id, device_name = detect_best_device(device, explicit_device_id=device_id)

    # Set the CuPy device if using GPU
    if use_gpu and CUPY_AVAILABLE:
        try:
            cp.cuda.Device(detected_device_id).use()
        except Exception as e:
            warnings.warn(f"Failed to set GPU device {detected_device_id}: {e}. Using CPU.")
            use_gpu = False
            device_name = "CPU"

    if verbose:
        print(f"  Using {device_name} for deconvolution")

    # Convert to float32
    original_dtype = stack.dtype
    stack_float = stack.astype(np.float32)

    # Apply FFT-efficient padding
    if pad_for_fft:
        original_shape = stack_float.shape

        # Calculate padding (add 4x PSF size, then find efficient FFT length)
        padded_shape = tuple(
            _find_good_fft_length(s + 4 * p)
            for s, p in zip(stack_float.shape, psf.shape)
        )

        pad_before = tuple((ps - os) // 2 for ps, os in zip(padded_shape, original_shape))
        pad_after = tuple(ps - os - pb for ps, os, pb in zip(padded_shape, original_shape, pad_before))

        # Use symmetric padding (matches MATLAB padarray with 'symmetric' mode)
        # MATLAB: block = padarray(block, [floor(pad_x) floor(pad_y) floor(pad_z)], 'pre', 'symmetric');
        #         block = padarray(block, [ceil(pad_x) ceil(pad_y) ceil(pad_z)], 'post', 'symmetric');
        pad_width = [(pb, pa) for pb, pa in zip(pad_before, pad_after)]
        stack_padded = np.pad(stack_float, pad_width, mode='symmetric')

        if verbose:
            print(f"  Padded from {original_shape} to {padded_shape} for FFT efficiency")
    else:
        stack_padded = stack_float
        pad_before = (0, 0, 0)
        original_shape = stack_float.shape

    # Run deconvolution
    if use_gpu:
        result_padded = lucy_richardson_gpu(
            stack_padded, psf, iterations, damping, stop_criterion, callback, verbose
        )
    else:
        result_padded = lucy_richardson_cpu(
            stack_padded, psf, iterations, damping, stop_criterion, callback, verbose
        )

    # Remove padding
    if pad_for_fft:
        result = result_padded[
            pad_before[0]:pad_before[0] + original_shape[0],
            pad_before[1]:pad_before[1] + original_shape[1],
            pad_before[2]:pad_before[2] + original_shape[2]
        ]
    else:
        result = result_padded

    return result


def deconvolve_stack(stack: np.ndarray, psf: np.ndarray,
                     iterations: int = 25,
                     damping: float = 0.0,
                     stop_criterion: float = 5.0,
                     device: str = 'auto',
                     device_id: Optional[int] = None,
                     max_memory_fraction: float = 0.8,
                     blend: bool = False,
                     verbose: bool = True) -> np.ndarray:
    """
    Deconvolve a 3D stack with automatic chunking for large images.

    This function automatically splits large images into chunks that fit
    in memory, processes them sequentially, and reassembles the result.

    Parameters
    ----------
    stack : ndarray
        Input 3D image stack
    psf : ndarray
        Point spread function
    iterations : int
        Maximum iterations (default 25)
    damping : float
        Damping factor (default 0)
    stop_criterion : float
        Convergence threshold (default 5.0)
    device : str
        'GPU', 'CPU', or 'auto'
    device_id : int, optional
        Specific GPU device ID for multi-GPU processing. If provided,
        this GPU will be used instead of auto-detection.
    max_memory_fraction : float
        Fraction of available memory to use per chunk (default 0.8)
    blend : bool
        If True, use cosine-weighted blending in overlap regions to eliminate
        boundary artifacts (recommended for oblique tissue acquisitions).
        If False, use core-only extraction (original MATLAB style).
        Default False to maintain backward compatibility.
    verbose : bool
        Print progress (default True)

    Returns
    -------
    ndarray
        Deconvolved image stack
    """
    from .chunking import calculate_chunks, process_in_chunks

    return process_in_chunks(
        stack, psf,
        deconv_func=lambda s, p: deconvolve(
            s, p, iterations, damping, stop_criterion, device, device_id, True, None, verbose
        ),
        device=device,
        device_id=device_id,
        max_memory_fraction=max_memory_fraction,
        overlap=None,  # Use default (4x PSF size) for proper FFT boundary handling
        blend=blend,
        verbose=verbose
    )
