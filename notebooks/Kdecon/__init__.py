"""
KDecon - Python implementation of Lucy-Richardson deconvolution for light sheet microscopy.

This module replaces the MATLAB-based deconvolution with a pure Python implementation
that supports GPU acceleration via CuPy and multi-threaded CPU processing.
"""

from .deconvolution import deconvolve, deconvolve_stack
from .psf import generate_psf, calculate_psf_size
from .main import KDecon, decon

__version__ = "1.0.0"
__all__ = [
    "KDecon",
    "decon",
    "deconvolve",
    "deconvolve_stack",
    "generate_psf",
    "calculate_psf_size",
]
