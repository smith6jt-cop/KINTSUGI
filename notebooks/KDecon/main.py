"""
Main KDecon module - high-level interface for deconvolution.

This module provides the main user-facing API that replaces the MATLAB
deconvolution implementation. It maintains compatibility with the
original parameter interface while using pure Python.
"""

import numpy as np
from pathlib import Path
from typing import Union, Optional, Tuple, Dict, Any
from datetime import datetime, timedelta
import time
import shutil
import warnings

from .psf import generate_psf, calculate_psf_size
from .deconvolution import deconvolve, deconvolve_stack, CUPY_AVAILABLE
from .chunking import get_available_memory, calculate_chunks
from .io import (
    get_stack_info, load_stack, save_stack,
    write_parameters_file
)


class KDecon:
    """
    KINTSUGI Deconvolution class for light sheet microscopy.

    This class provides GPU-accelerated or multi-threaded CPU Lucy-Richardson
    deconvolution with automatic handling of large images.

    Parameters
    ----------
    dxy : float
        XY voxel size in nanometers
    dz : float
        Z voxel size in nanometers
    NA : float
        Numerical aperture of objective
    rf : float
        Refractive index of medium/sample
    lambda_ex : float
        Excitation wavelength in nanometers
    lambda_em : float
        Emission wavelength in nanometers
    iterations : int
        Maximum Lucy-Richardson iterations (default 25)
    damping : float
        Damping factor 0-1 (default 0). Increase for noisy images.
    stop_criterion : float
        Stop if relative change < this % (default 5.0)
    device : str
        'GPU', 'CPU', or 'auto' (default 'auto')
    max_memory_fraction : float
        Max fraction of memory to use per chunk (default 0.8)
    fcyl : float, optional
        Focal length of cylinder lens (mm) - for metadata only
    slitwidth : float, optional
        Slit aperture width (mm) - for metadata only
    clipval : float
        Histogram clipping percentage 0-5 (default 0.01)
    verbose : bool
        Print progress information (default True)

    Examples
    --------
    >>> decon = KDecon(
    ...     dxy=377, dz=1500,
    ...     NA=0.75, rf=1.44,
    ...     lambda_ex=560, lambda_em=575,
    ...     device='GPU'
    ... )
    >>> decon.process('/path/to/input', '/path/to/output')
    """

    def __init__(self,
                 dxy: float,
                 dz: float,
                 NA: float,
                 rf: float,
                 lambda_ex: float,
                 lambda_em: float,
                 iterations: int = 25,
                 damping: float = 0.0,
                 stop_criterion: float = 5.0,
                 device: str = 'auto',
                 max_memory_fraction: float = 0.8,
                 fcyl: float = None,
                 slitwidth: float = None,
                 clipval: float = 0.01,
                 verbose: bool = True):

        self.dxy = dxy
        self.dz = dz
        self.NA = NA
        self.rf = rf
        self.lambda_ex = lambda_ex
        self.lambda_em = lambda_em
        self.iterations = iterations
        self.damping = damping
        self.stop_criterion = stop_criterion
        self.device = device
        self.max_memory_fraction = max_memory_fraction
        self.fcyl = fcyl
        self.slitwidth = slitwidth
        self.clipval = clipval
        self.verbose = verbose

        # Will be computed on first use
        self._psf = None
        self._psf_info = None

    @property
    def psf(self) -> np.ndarray:
        """Get or compute the PSF."""
        if self._psf is None:
            self._compute_psf()
        return self._psf

    @property
    def psf_info(self) -> dict:
        """Get PSF information."""
        if self._psf_info is None:
            self._compute_psf()
        return self._psf_info

    def _compute_psf(self):
        """Compute the PSF based on current parameters."""
        self._psf, self._psf_info = generate_psf(
            self.dxy, self.dz, self.NA, self.rf,
            self.lambda_ex, self.lambda_em,
            verbose=self.verbose
        )

    def process(self,
                input_path: Union[str, Path],
                output_path: Optional[Union[str, Path]] = None,
                bit_depth: int = 16) -> np.ndarray:
        """
        Process a stack of TIFF images.

        Parameters
        ----------
        input_path : str or Path
            Directory containing input TIFF files
        output_path : str or Path, optional
            Directory for output files. If None, creates 'deconvolved'
            subdirectory in input_path.
        bit_depth : int
            Output bit depth: 16 or 32 (default 16)

        Returns
        -------
        ndarray
            Deconvolved image stack
        """
        start_time = time.time()
        input_path = Path(input_path)

        if output_path is None:
            output_path = input_path / 'deconvolved'
        else:
            output_path = Path(output_path)

        if self.verbose:
            print("=" * 60)
            print("KDecon - KINTSUGI Python Deconvolution")
            print("=" * 60)

        # Get available memory and determine device
        available_mem, device_name = get_available_memory(self.device)

        if self.verbose:
            print(f"\nDevice: {device_name}")
            print(f"Available memory: {available_mem / 1e9:.2f} GB")

        # Load stack
        if self.verbose:
            print(f"\nInput: {input_path}")

        stack = load_stack(input_path, dtype=np.float32, verbose=self.verbose)

        # Get raw max for scaling
        raw_max = np.max(stack)

        # Determine number of chunks
        n_chunks = calculate_chunks(
            stack.shape,
            stack.dtype,
            int(available_mem * self.max_memory_fraction),
            memory_factor=6.0,
            overlap=max(self.psf.shape)
        )

        if self.verbose:
            print(f"\nImage dimensions: {stack.shape[0]} x {stack.shape[1]} x {stack.shape[2]}")
            if np.prod(n_chunks) > 1:
                print(f"Splitting into {n_chunks[0]} x {n_chunks[1]} x {n_chunks[2]} blocks")

        # Run deconvolution
        if self.verbose:
            print("\nStarting deconvolution...")

        if np.prod(n_chunks) > 1:
            result = deconvolve_stack(
                stack, self.psf,
                iterations=self.iterations,
                damping=self.damping,
                stop_criterion=self.stop_criterion,
                device=self.device,
                max_memory_fraction=self.max_memory_fraction,
                verbose=self.verbose
            )
        else:
            result = deconvolve(
                stack, self.psf,
                iterations=self.iterations,
                damping=self.damping,
                stop_criterion=self.stop_criterion,
                device=self.device,
                verbose=self.verbose
            )

        # Apply histogram clipping if requested
        deconv_max = np.max(result)
        deconv_min = np.min(result)

        if self.clipval > 0:
            if self.verbose:
                print("\nApplying histogram clipping...")

            # Calculate percentiles for clipping
            low_clip = np.percentile(result, self.clipval)
            high_clip = np.percentile(result, 100 - self.clipval)

            result = np.clip(result, low_clip, high_clip)
            result = (result - low_clip) / (high_clip - low_clip)
        else:
            # Scale to 0-1
            result = (result - deconv_min) / (deconv_max - deconv_min)

        # Scale to output range
        if raw_max <= 65535:
            scale = 65535
        else:
            scale = raw_max

        result = result * scale

        # Save results
        if self.verbose:
            print(f"\nOutput: {output_path}")

        save_stack(result, output_path, bit_depth=bit_depth, verbose=self.verbose)

        # Calculate elapsed time
        elapsed = time.time() - start_time
        elapsed_str = str(timedelta(seconds=int(elapsed)))

        # Write parameters file
        params = {
            'dxy': self.dxy,
            'dz': self.dz,
            'NA': self.NA,
            'rf': self.rf,
            'lambda_ex': self.lambda_ex,
            'lambda_em': self.lambda_em,
            'iterations': self.iterations,
            'damping': self.damping,
            'stop_criterion': self.stop_criterion,
            'fcyl': self.fcyl,
            'slitwidth': self.slitwidth,
            'clipval': self.clipval,
            'inpath': str(input_path),
        }

        write_parameters_file(
            output_path, params, self.psf_info,
            elapsed_str, device_name, n_chunks, np.max(result)
        )

        if self.verbose:
            print("\n" + "=" * 60)
            print(f"Deconvolution complete!")
            print(f"Elapsed time: {elapsed_str}")
            print("=" * 60)

        return result

    def process_array(self, stack: np.ndarray) -> np.ndarray:
        """
        Process a numpy array directly.

        Parameters
        ----------
        stack : ndarray
            3D image stack (x, y, z)

        Returns
        -------
        ndarray
            Deconvolved image stack
        """
        if stack.ndim != 3:
            raise ValueError(f"Expected 3D array, got {stack.ndim}D")

        return deconvolve_stack(
            stack, self.psf,
            iterations=self.iterations,
            damping=self.damping,
            stop_criterion=self.stop_criterion,
            device=self.device,
            max_memory_fraction=self.max_memory_fraction,
            verbose=self.verbose
        )


def decon(base_dir: Union[str, Path],
          stitch_dir: Union[str, Path],
          dec_cycle: int,
          dec_channel: int,
          xy_vox: float = 377,
          z_vox: float = 1500,
          iterations: int = 25,
          mic_NA: float = 0.75,
          tissue_RI: float = 1.44,
          slit_aper: float = 6.5,
          f_cyl: float = 1,
          damping: float = 0,
          hist_clip: float = 0.01,
          stop_criterion: float = 5.0,
          max_memory: float = 0.8,
          device: str = 'auto',
          wavelengths: Optional[Dict[int, Tuple[float, float]]] = None) -> str:
    """
    High-level deconvolution function matching the original MATLAB interface.

    This function is a drop-in replacement for the MATLAB-based decon function
    used in the KINTSUGI notebooks.

    Parameters
    ----------
    base_dir : str or Path
        Base directory for project
    stitch_dir : str or Path
        Directory containing stitched images
    dec_cycle : int
        Cycle number to process
    dec_channel : int
        Channel number to process
    xy_vox : float
        XY pixel size in nanometers (default 377)
    z_vox : float
        Z pixel size in nanometers (default 1500)
    iterations : int
        Max Lucy-Richardson iterations (default 25)
    mic_NA : float
        Microscope numerical aperture (default 0.75)
    tissue_RI : float
        Tissue refractive index (default 1.44)
    slit_aper : float
        Slit aperture width in mm (default 6.5)
    f_cyl : float
        Cylinder lens focal length in mm (default 1)
    damping : float
        Damping factor 0-10% (default 0)
    hist_clip : float
        Histogram clip percentage 0-5% (default 0.01)
    stop_criterion : float
        Stop criterion percentage (default 5.0)
    max_memory : float
        Max memory fraction to use (default 0.8)
    device : str
        'GPU', 'CPU', or 'auto' (default 'auto')
    wavelengths : dict, optional
        Channel wavelengths dict {channel: (ex, em)}. If None, uses defaults.

    Returns
    -------
    str
        Path to output directory
    """
    base_dir = Path(base_dir)
    stitch_dir = Path(stitch_dir)

    # Default wavelengths from original implementation
    if wavelengths is None:
        wavelengths = {
            1: (358, 461),
            2: (753, 775),
            3: (560, 575),
            4: (648, 668)
        }

    # Get wavelengths for this channel
    if dec_channel not in wavelengths:
        raise ValueError(f"Channel {dec_channel} not in wavelengths dict")

    lambda_ex, lambda_em = wavelengths[dec_channel]

    # Set up paths
    source = stitch_dir / f"cyc{dec_cycle:02d}" / f"CH{dec_channel}"
    decon_dir = Path(str(stitch_dir).replace('_BaSiC_Stitched', '_Decon'))
    dest = decon_dir / f"cyc{dec_cycle:02d}" / f"CH{dec_channel}"

    dest.mkdir(parents=True, exist_ok=True)

    # Convert damping from percentage to fraction
    damping_frac = damping / 100.0 if damping > 1 else damping

    # Create deconvolution object and process
    kdecon = KDecon(
        dxy=xy_vox,
        dz=z_vox,
        NA=mic_NA,
        rf=tissue_RI,
        lambda_ex=lambda_ex,
        lambda_em=lambda_em,
        iterations=iterations,
        damping=damping_frac,
        stop_criterion=stop_criterion,
        device=device,
        max_memory_fraction=max_memory,
        fcyl=f_cyl,
        slitwidth=slit_aper,
        clipval=hist_clip,
        verbose=True
    )

    # Process
    output_dir = dest / 'deconvolved'
    kdecon.process(source, output_dir)

    return str(output_dir)


def check_gpu_available() -> Tuple[bool, str]:
    """
    Check if GPU is available for deconvolution.

    Returns
    -------
    tuple
        (available, message)
    """
    if not CUPY_AVAILABLE:
        return False, "CuPy is not installed. Install with: pip install cupy-cuda11x"

    try:
        import cupy as cp
        device = cp.cuda.Device()
        mem = device.mem_info

        return True, f"GPU available: {device.name} with {mem[1] / 1e9:.1f} GB"
    except Exception as e:
        return False, f"GPU not accessible: {e}"


def print_info():
    """Print KDecon version and system information."""
    from . import __version__

    print(f"KDecon v{__version__}")
    print("-" * 40)

    gpu_available, gpu_msg = check_gpu_available()
    print(f"GPU: {gpu_msg}")

    try:
        import psutil
        mem = psutil.virtual_memory()
        print(f"CPU Memory: {mem.total / 1e9:.1f} GB total, {mem.available / 1e9:.1f} GB available")
    except ImportError:
        print("CPU Memory: psutil not installed")

    print("-" * 40)
