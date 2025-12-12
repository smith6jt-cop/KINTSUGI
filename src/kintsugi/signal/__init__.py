"""
KINTSUGI Signal Isolation Module

Pure Python implementation of autofluorescence subtraction and signal isolation
algorithms, mirroring the original Kutils.py functionality with added intelligence,
automation, and learning capabilities.

Main Functions:
- subtract_autofluorescence: Core blank subtraction algorithm
- analyze_for_subtraction: Intelligent parameter suggestion
- AutofluorescenceSubtractor: Class-based interface with learning

Example:
    >>> from kintsugi.signal import subtract_autofluorescence, analyze_for_subtraction
    >>>
    >>> # Get suggested parameters
    >>> params = analyze_for_subtraction(signal_image, blank_image)
    >>> print(f"Suggested: clip={params['blank_clip_factor']}, scale={params['blank_scale_factor']}")
    >>>
    >>> # Apply subtraction
    >>> result = subtract_autofluorescence(signal_image, blank_image, **params)
"""

from .autofluorescence import (
    analyze_for_subtraction,
    compute_subtraction_quality,
    subtract_autofluorescence,
    subtract_autofluorescence_dask,
    suggest_blank_channel,
)
from .subtractor import (
    AutofluorescenceSubtractor,
    SubtractionParameters,
    SubtractionResult,
)
from .utils import (
    apply_erosion_mask,
    apply_percentile_smoothing,
    clip_and_scale_blank,
    estimate_autofluorescence_level,
)

__all__ = [
    # Core functions
    "subtract_autofluorescence",
    "subtract_autofluorescence_dask",
    "analyze_for_subtraction",
    "suggest_blank_channel",
    "compute_subtraction_quality",
    # Class-based interface
    "AutofluorescenceSubtractor",
    "SubtractionResult",
    "SubtractionParameters",
    # Utilities
    "clip_and_scale_blank",
    "apply_percentile_smoothing",
    "apply_erosion_mask",
    "estimate_autofluorescence_level",
]
