"""Kstitch - KINTSUGI Image Stitching Module.

A MIST-inspired stitching implementation in Python for microscope images,
with GPU acceleration via CuPy.

Based on m2stitch by Yohsuke T. Fukai, extended with GPU support.
"""

__author__ = """Yohsuke T. Fukai"""
__email__ = "ysk@yfukai.net"

from .stitching import stitch_images
from .stitching import check_gpu_status
from .stitching import HAS_CUPY, GPU_ERROR_MSG, GPU_DEVICE_INFO

__all__ = [
    "stitch_images",
    "check_gpu_status",
    "HAS_CUPY",
    "GPU_ERROR_MSG",
    "GPU_DEVICE_INFO",
]
