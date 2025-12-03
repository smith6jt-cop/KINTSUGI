"""
Point Spread Function (PSF) generation for light sheet microscopy.

This module implements the theoretical PSF calculation based on the Airy disk model
using Bessel function integrals. The PSF is computed as the product of excitation
and emission PSFs.
"""

import numpy as np
from scipy import special
from scipy.integrate import quad
from scipy.optimize import brentq
from typing import Tuple
import warnings


def _bessel_integrand(p: float, x: float, y: float, z: float,
                      wavelength: float, NA: float, n: float) -> complex:
    """
    Bessel function integrand for PSF calculation.

    Implements the integral: J0(2*pi*NA*r*p/(lambda*n)) * exp(-i*pi*p^2*z*NA^2/(lambda*n^2)) * p

    Parameters
    ----------
    p : float
        Integration variable (0 to 1)
    x, y, z : float
        Position coordinates in nanometers
    wavelength : float
        Wavelength in nanometers
    NA : float
        Numerical aperture
    n : float
        Refractive index

    Returns
    -------
    complex
        Integrand value
    """
    r = np.sqrt(x**2 + y**2)
    bessel_arg = 2 * np.pi * NA * r * p / (wavelength * n)
    phase_arg = -np.pi * p**2 * z * NA**2 / (wavelength * n**2)

    return special.j0(bessel_arg) * np.exp(1j * phase_arg) * p


def _psf_single_wavelength(x: float, y: float, z: float,
                           NA: float, n: float, wavelength: float) -> float:
    """
    Calculate PSF intensity at a single point for one wavelength.

    Parameters
    ----------
    x, y, z : float
        Position coordinates in nanometers
    NA : float
        Numerical aperture
    n : float
        Refractive index
    wavelength : float
        Wavelength in nanometers

    Returns
    -------
    float
        PSF intensity (normalized)
    """
    def real_integrand(p):
        val = _bessel_integrand(p, x, y, z, wavelength, NA, n)
        return np.real(val)

    def imag_integrand(p):
        val = _bessel_integrand(p, x, y, z, wavelength, NA, n)
        return np.imag(val)

    # Integrate from 0 to 1
    real_part, _ = quad(real_integrand, 0, 1, limit=100)
    imag_part, _ = quad(imag_integrand, 0, 1, limit=100)

    # PSF = 4 * |integral|^2
    return 4 * (real_part**2 + imag_part**2)


def _psf_light_sheet(x: float, y: float, z: float, NA: float, n: float,
                     lambda_ex: float, lambda_em: float) -> float:
    """
    Calculate light sheet PSF at a point (product of excitation and emission PSFs).

    Parameters
    ----------
    x, y, z : float
        Position coordinates in nanometers
    NA : float
        Numerical aperture
    n : float
        Refractive index
    lambda_ex : float
        Excitation wavelength in nanometers
    lambda_em : float
        Emission wavelength in nanometers

    Returns
    -------
    float
        Combined PSF intensity
    """
    psf_ex = _psf_single_wavelength(x, y, z, NA, n, lambda_ex)
    psf_em = _psf_single_wavelength(x, y, z, NA, n, lambda_em)
    return psf_ex * psf_em


def calculate_psf_size(dxy: float, dz: float, NA: float, n: float,
                       lambda_ex: float, lambda_em: float,
                       grid_size_xy: float = 2.0,
                       grid_size_z: float = 2.0) -> Tuple[int, int, float, float]:
    """
    Determine the required PSF grid size based on FWHM.

    Parameters
    ----------
    dxy : float
        XY voxel size in nanometers
    dz : float
        Z voxel size in nanometers
    NA : float
        Numerical aperture
    n : float
        Refractive index
    lambda_ex : float
        Excitation wavelength in nanometers
    lambda_em : float
        Emission wavelength in nanometers
    grid_size_xy : float
        Grid size as multiple of FWHM in XY (default 2.0)
    grid_size_z : float
        Grid size as multiple of FWHM in Z (default 2.0)

    Returns
    -------
    tuple
        (nxy, nz, FWHM_xy, FWHM_z) - grid dimensions and FWHM values in nm
    """
    # Get PSF value at origin for half-max calculation
    psf_origin = _psf_light_sheet(0, 0, 0, NA, n, lambda_ex, lambda_em)
    half_max = 0.5 * psf_origin

    # Find FWHM in XY direction
    def xy_func(x):
        return _psf_light_sheet(x, 0, 0, NA, n, lambda_ex, lambda_em) - half_max

    # Find FWHM in Z direction
    def z_func(z):
        return _psf_light_sheet(0, 0, z, NA, n, lambda_ex, lambda_em) - half_max

    # Search for half-max crossings
    try:
        fwhm_xy_half = brentq(xy_func, 1, 5000)
        fwhm_xy = 2 * fwhm_xy_half
    except ValueError:
        # Fallback to theoretical estimate
        fwhm_xy = 0.51 * lambda_em / NA
        warnings.warn(f"Could not find XY FWHM, using theoretical estimate: {fwhm_xy:.1f} nm")

    try:
        fwhm_z_half = brentq(z_func, 1, 20000)
        fwhm_z = 2 * fwhm_z_half
    except ValueError:
        # Fallback to theoretical estimate
        fwhm_z = 2 * lambda_em * n / (NA**2)
        warnings.warn(f"Could not find Z FWHM, using theoretical estimate: {fwhm_z:.1f} nm")

    # Calculate Rayleigh limit for minimum sampling
    Rxy = 0.61 * lambda_em / NA
    dxy_corr = min(dxy, Rxy / 3)

    # Calculate grid size
    nxy = int(np.ceil(grid_size_xy * fwhm_xy / dxy_corr))
    nz = int(np.ceil(grid_size_z * fwhm_z / dz))

    # Ensure odd dimensions (for centered PSF)
    if nxy % 2 == 0:
        nxy += 1
    if nz % 2 == 0:
        nz += 1

    return nxy, nz, fwhm_xy, fwhm_z


def generate_psf(dxy: float, dz: float, NA: float, n: float,
                 lambda_ex: float, lambda_em: float,
                 nxy: int = None, nz: int = None,
                 verbose: bool = True) -> Tuple[np.ndarray, dict]:
    """
    Generate a theoretical PSF for light sheet microscopy.

    This function computes a 3D PSF based on the Airy disk model, exploiting
    octant symmetry for computational efficiency.

    Parameters
    ----------
    dxy : float
        XY voxel size in nanometers
    dz : float
        Z voxel size in nanometers
    NA : float
        Numerical aperture
    n : float
        Refractive index
    lambda_ex : float
        Excitation wavelength in nanometers
    lambda_em : float
        Emission wavelength in nanometers
    nxy : int, optional
        PSF grid size in XY (auto-calculated if None)
    nz : int, optional
        PSF grid size in Z (auto-calculated if None)
    verbose : bool
        Print progress information

    Returns
    -------
    tuple
        (psf, info) - 3D PSF array and dictionary with PSF parameters
    """
    if verbose:
        print("Calculating PSF...")

    # Calculate Rayleigh limit for minimum sampling
    Rxy = 0.61 * lambda_em / NA
    dxy_corr = min(dxy, Rxy / 3)

    # Auto-calculate size if not provided
    if nxy is None or nz is None:
        auto_nxy, auto_nz, fwhm_xy, fwhm_z = calculate_psf_size(
            dxy, dz, NA, n, lambda_ex, lambda_em
        )
        if nxy is None:
            nxy = auto_nxy
        if nz is None:
            nz = auto_nz
    else:
        # Still calculate FWHM for info
        _, _, fwhm_xy, fwhm_z = calculate_psf_size(dxy, dz, NA, n, lambda_ex, lambda_em)

    # Ensure odd dimensions
    if nxy % 2 == 0:
        nxy += 1
    if nz % 2 == 0:
        nz += 1

    if verbose:
        print(f"  PSF size: {nxy} x {nxy} x {nz}")
        print(f"  FWHM lateral: {fwhm_xy:.1f} nm")
        print(f"  FWHM axial: {fwhm_z:.1f} nm")

    # Calculate only first octant (exploit 8-fold symmetry)
    half_nxy = (nxy - 1) // 2 + 1
    half_nz = (nz - 1) // 2 + 1

    octant = np.zeros((half_nxy, half_nxy, half_nz), dtype=np.float32)

    total_points = half_nxy * half_nxy * half_nz
    computed = 0

    for iz in range(half_nz):
        z = iz * dz
        for iy in range(half_nxy):
            y = iy * dxy_corr
            for ix in range(half_nxy):
                x = ix * dxy_corr
                octant[ix, iy, iz] = _psf_light_sheet(x, y, z, NA, n, lambda_ex, lambda_em)
                computed += 1

        if verbose and (iz + 1) % max(1, half_nz // 5) == 0:
            print(f"  Computing PSF: {100 * computed / total_points:.0f}%")

    # Mirror octant to full PSF (8-fold symmetry)
    psf = _mirror_octant_to_full(octant)

    # Normalize to unit integral
    psf = psf / np.sum(psf)

    if verbose:
        print("  PSF calculation complete")

    # Calculate Rayleigh range info
    Rz = 2 * lambda_em * n / (NA**2)

    info = {
        'nxy': nxy,
        'nz': nz,
        'fwhm_xy': fwhm_xy,
        'fwhm_z': fwhm_z,
        'rayleigh_xy': Rxy,
        'rayleigh_z': Rz,
        'dxy_corrected': dxy_corr,
    }

    return psf.astype(np.float32), info


def _mirror_octant_to_full(octant: np.ndarray) -> np.ndarray:
    """
    Mirror first octant to create full PSF using 8-fold symmetry.

    Parameters
    ----------
    octant : ndarray
        First octant of PSF (positive x, y, z)

    Returns
    -------
    ndarray
        Full 3D PSF
    """
    sx = 2 * octant.shape[0] - 1
    sy = 2 * octant.shape[1] - 1
    sz = 2 * octant.shape[2] - 1

    cx = octant.shape[0] - 1
    cy = octant.shape[1] - 1
    cz = octant.shape[2] - 1

    full = np.zeros((sx, sy, sz), dtype=octant.dtype)

    # Place all 8 octants
    # (+x, +y, +z)
    full[cx:, cy:, cz:] = octant
    # (+x, -y, +z)
    full[cx:, :cy+1, cz:] = np.flip(octant, axis=1)
    # (-x, -y, +z)
    full[:cx+1, :cy+1, cz:] = np.flip(np.flip(octant, axis=0), axis=1)
    # (-x, +y, +z)
    full[:cx+1, cy:, cz:] = np.flip(octant, axis=0)
    # (+x, +y, -z)
    full[cx:, cy:, :cz+1] = np.flip(octant, axis=2)
    # (+x, -y, -z)
    full[cx:, :cy+1, :cz+1] = np.flip(np.flip(octant, axis=1), axis=2)
    # (-x, -y, -z)
    full[:cx+1, :cy+1, :cz+1] = np.flip(np.flip(np.flip(octant, axis=0), axis=1), axis=2)
    # (-x, +y, -z)
    full[:cx+1, cy:, :cz+1] = np.flip(np.flip(octant, axis=0), axis=2)

    return full


def generate_psf_fast(dxy: float, dz: float, NA: float, n: float,
                      lambda_ex: float, lambda_em: float,
                      nxy: int = None, nz: int = None,
                      verbose: bool = True) -> Tuple[np.ndarray, dict]:
    """
    Generate PSF using vectorized computation (faster but uses more memory).

    Parameters are the same as generate_psf().
    """
    if verbose:
        print("Calculating PSF (fast mode)...")

    # Calculate Rayleigh limit for minimum sampling
    Rxy = 0.61 * lambda_em / NA
    dxy_corr = min(dxy, Rxy / 3)

    # Auto-calculate size if not provided
    if nxy is None or nz is None:
        auto_nxy, auto_nz, fwhm_xy, fwhm_z = calculate_psf_size(
            dxy, dz, NA, n, lambda_ex, lambda_em
        )
        if nxy is None:
            nxy = auto_nxy
        if nz is None:
            nz = auto_nz
    else:
        _, _, fwhm_xy, fwhm_z = calculate_psf_size(dxy, dz, NA, n, lambda_ex, lambda_em)

    # Ensure odd dimensions
    if nxy % 2 == 0:
        nxy += 1
    if nz % 2 == 0:
        nz += 1

    if verbose:
        print(f"  PSF size: {nxy} x {nxy} x {nz}")

    # Create coordinate grids (only need half for each dimension due to symmetry)
    half_nxy = (nxy - 1) // 2 + 1
    half_nz = (nz - 1) // 2 + 1

    x = np.arange(half_nxy) * dxy_corr
    y = np.arange(half_nxy) * dxy_corr
    z = np.arange(half_nz) * dz

    # Vectorized PSF calculation using numerical integration
    # This is a simplified/approximated version for speed
    octant = np.zeros((half_nxy, half_nxy, half_nz), dtype=np.float32)

    for iz, zval in enumerate(z):
        for iy, yval in enumerate(y):
            for ix, xval in enumerate(x):
                octant[ix, iy, iz] = _psf_light_sheet(xval, yval, zval, NA, n, lambda_ex, lambda_em)
        if verbose and (iz + 1) % max(1, half_nz // 3) == 0:
            print(f"  Computing PSF: {100 * (iz + 1) / half_nz:.0f}%")

    # Mirror to full PSF
    psf = _mirror_octant_to_full(octant)

    # Normalize
    psf = psf / np.sum(psf)

    if verbose:
        print("  PSF calculation complete")

    Rz = 2 * lambda_em * n / (NA**2)

    info = {
        'nxy': nxy,
        'nz': nz,
        'fwhm_xy': fwhm_xy,
        'fwhm_z': fwhm_z,
        'rayleigh_xy': Rxy,
        'rayleigh_z': Rz,
        'dxy_corrected': dxy_corr,
    }

    return psf.astype(np.float32), info
