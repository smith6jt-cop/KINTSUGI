"""
I/O utilities for reading and writing TIFF image stacks.

This module provides functions for loading TIFF stacks from directories
and saving deconvolved results.
"""

import numpy as np
from pathlib import Path
from typing import Tuple, Optional, Union, List
import warnings

try:
    from skimage.io import imread, imsave
    SKIMAGE_AVAILABLE = True
except ImportError:
    SKIMAGE_AVAILABLE = False

try:
    import tifffile
    TIFFFILE_AVAILABLE = True
except ImportError:
    TIFFFILE_AVAILABLE = False


def get_stack_info(directory: Union[str, Path]) -> Tuple[int, int, int, List[Path]]:
    """
    Get information about TIFF files in a directory.

    Parameters
    ----------
    directory : str or Path
        Directory containing TIFF files

    Returns
    -------
    tuple
        (nx, ny, nz, file_list) - dimensions and sorted list of files
    """
    directory = Path(directory)

    if not directory.exists():
        raise FileNotFoundError(f"Directory not found: {directory}")

    # Find all TIFF files
    tiff_files = sorted(list(directory.glob("*.tif")) + list(directory.glob("*.tiff")))

    if len(tiff_files) == 0:
        return 0, 0, 0, []

    # Read first file to get dimensions
    if TIFFFILE_AVAILABLE:
        sample = tifffile.imread(tiff_files[0])
    elif SKIMAGE_AVAILABLE:
        sample = imread(tiff_files[0])
    else:
        raise ImportError("Either tifffile or scikit-image is required")

    # Dimensions: (y, x) from image, z from number of files
    ny, nx = sample.shape[:2]
    nz = len(tiff_files)

    return nx, ny, nz, tiff_files


def load_stack(directory: Union[str, Path],
               dtype: np.dtype = np.float32,
               verbose: bool = True) -> np.ndarray:
    """
    Load a TIFF stack from a directory of 2D images.

    Parameters
    ----------
    directory : str or Path
        Directory containing TIFF files
    dtype : dtype
        Output data type (default float32)
    verbose : bool
        Print progress

    Returns
    -------
    ndarray
        3D array with shape (nx, ny, nz)
    """
    nx, ny, nz, tiff_files = get_stack_info(directory)

    if nz == 0:
        raise ValueError(f"No TIFF files found in {directory}")

    if verbose:
        print(f"Loading {nz} images ({nx} x {ny})...")

    # Allocate output array (x, y, z) - MATLAB convention
    stack = np.zeros((nx, ny, nz), dtype=dtype)

    for k, fpath in enumerate(tiff_files):
        if TIFFFILE_AVAILABLE:
            img = tifffile.imread(fpath)
        else:
            img = imread(fpath)

        # Transpose from (y, x) to (x, y)
        stack[:, :, k] = img.T.astype(dtype)

        if verbose and (k + 1) % max(1, nz // 5) == 0:
            print(f"  Loaded {k + 1}/{nz} images")

    if verbose:
        print(f"  Stack loaded: shape = {stack.shape}")

    return stack


def load_block(directory: Union[str, Path],
               x_range: Tuple[int, int],
               y_range: Tuple[int, int],
               z_range: Tuple[int, int],
               dtype: np.dtype = np.float32) -> np.ndarray:
    """
    Load a block of data from a TIFF stack.

    Parameters
    ----------
    directory : str or Path
        Directory containing TIFF files
    x_range : tuple
        (start, end) in x dimension
    y_range : tuple
        (start, end) in y dimension
    z_range : tuple
        (start, end) in z dimension (file indices)
    dtype : dtype
        Output data type

    Returns
    -------
    ndarray
        3D block with shape (x_size, y_size, z_size)
    """
    directory = Path(directory)
    _, _, _, tiff_files = get_stack_info(directory)

    x1, x2 = x_range
    y1, y2 = y_range
    z1, z2 = z_range

    # Allocate output
    block = np.zeros((x2 - x1, y2 - y1, z2 - z1), dtype=dtype)

    for k in range(z1, z2):
        if TIFFFILE_AVAILABLE:
            img = tifffile.imread(tiff_files[k])
        else:
            img = imread(tiff_files[k])

        # Read region and transpose
        block[:, :, k - z1] = img[y1:y2, x1:x2].T.astype(dtype)

    return block


def save_stack(stack: np.ndarray,
               directory: Union[str, Path],
               prefix: str = "deconv_",
               bit_depth: int = 16,
               verbose: bool = True) -> None:
    """
    Save a 3D stack as a series of TIFF images.

    Parameters
    ----------
    stack : ndarray
        3D array with shape (nx, ny, nz)
    directory : str or Path
        Output directory
    prefix : str
        Filename prefix (default "deconv_")
    bit_depth : int
        Output bit depth: 16 or 32 (default 16)
    verbose : bool
        Print progress
    """
    directory = Path(directory)
    directory.mkdir(parents=True, exist_ok=True)

    nz = stack.shape[2]

    if verbose:
        print(f"Saving {nz} images to {directory}...")

    for k in range(nz):
        # Generate filename with zero-padded index
        filename = f"{prefix}{k + 1:03d}.tif"
        filepath = directory / filename

        # Extract slice and transpose back to (y, x)
        img_slice = stack[:, :, k].T

        if bit_depth == 16:
            # Scale to uint16
            img_out = np.clip(img_slice, 0, 65535).astype(np.uint16)
        elif bit_depth == 32:
            img_out = img_slice.astype(np.float32)
        else:
            raise ValueError(f"Unsupported bit depth: {bit_depth}")

        if TIFFFILE_AVAILABLE:
            if bit_depth == 32:
                tifffile.imwrite(filepath, img_out, compression='lzw')
            else:
                tifffile.imwrite(filepath, img_out)
        else:
            imsave(filepath, img_out, check_contrast=False)

        if verbose and (k + 1) % max(1, nz // 5) == 0:
            print(f"  Saved {k + 1}/{nz} images")

    if verbose:
        print(f"  Stack saved to {directory}")


def write_parameters_file(directory: Union[str, Path],
                          params: dict,
                          psf_info: dict,
                          elapsed_time: str,
                          device_used: str,
                          n_blocks: Tuple[int, int, int],
                          max_intensity: float) -> None:
    """
    Write deconvolution parameters to a text file.

    Parameters
    ----------
    directory : str or Path
        Output directory
    params : dict
        Deconvolution parameters
    psf_info : dict
        PSF information from generate_psf()
    elapsed_time : str
        Elapsed time string
    device_used : str
        Device used for processing
    n_blocks : tuple
        Number of blocks in each dimension
    max_intensity : float
        Maximum intensity in deconvolved data
    """
    from datetime import datetime

    directory = Path(directory)
    filepath = directory / "DECONV_parameters.txt"

    with open(filepath, 'w') as f:
        f.write(f"deconvolution finished at {datetime.now()}\n")
        f.write(f"data processed on: {device_used}\n")
        f.write(f"elapsed time: {elapsed_time}\n")
        f.write("\n")
        f.write(f"highest intensity value in deconvolved data: {max_intensity:.2f}\n")
        f.write(f"focal length of cylinder lens (mm): {params.get('fcyl', 'N/A')}\n")
        f.write(f"width of slit aperture (mm): {params.get('slitwidth', 'N/A')}\n")
        f.write(f"histogram clipping value (%): {params.get('clipval', 'N/A')}\n")
        f.write(f"numerical aperture: {params['NA']}\n")
        f.write(f"excitation wavelength (nm): {params['lambda_ex']}\n")
        f.write(f"emission wavelength (nm): {params['lambda_em']}\n")
        f.write(f"refractive index: {params['rf']}\n")
        f.write(f"max. iterations: {params['iterations']}\n")
        f.write(f"damping factor (%): {params['damping'] * 100}\n")
        f.write(f"stop criterion (%): {params['stop_criterion']}\n")
        f.write("\n")
        f.write(f"source data folder: {params.get('inpath', 'N/A')}\n")
        f.write(f"number of blocks (x y z): {n_blocks[0]} x {n_blocks[1]} x {n_blocks[2]}\n")
        f.write("\n")
        f.write(f"voxel size: {params['dxy']} nm x {params['dxy']} nm x {params['dz']} nm\n")
        f.write(f"size of PSF (pixel): {psf_info['nxy']} x {psf_info['nxy']} x {psf_info['nz']}\n")
        f.write(f"FWHM of PSF lateral (nm): {psf_info['fwhm_xy']:.1f}\n")
        f.write(f"FWHM of PSF axial (nm): {psf_info['fwhm_z']:.1f}\n")
        f.write(f"Rayleigh range of objective lateral (nm): {psf_info['rayleigh_xy']:.1f}\n")
        f.write(f"Rayleigh range of objective axial (nm): {psf_info['rayleigh_z']:.1f}\n")
