"""
KINTSUGI Deprecation Module

Manages deprecation warnings and provides transition guidance for
deprecated notebooks and functionality.
"""

from __future__ import annotations

import warnings
from collections.abc import Callable
from functools import wraps


class KintsugiDeprecationWarning(UserWarning):
    """Custom warning for KINTSUGI deprecations."""

    pass


def deprecated_notebook(
    notebook_name: str,
    replacement: str,
    version: str = "2.0.0",
    details: str | None = None,
):
    """
    Decorator to mark notebook functions as deprecated.

    Parameters
    ----------
    notebook_name : str
        Name of the deprecated notebook
    replacement : str
        Name/description of the replacement
    version : str
        Version when removal is planned
    details : str, optional
        Additional details about the transition
    """

    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            message = (
                f"\n"
                f"{'=' * 60}\n"
                f"DEPRECATION NOTICE: {notebook_name}\n"
                f"{'=' * 60}\n"
                f"\n"
                f"This notebook/functionality is deprecated and will be removed\n"
                f"in KINTSUGI version {version}.\n"
                f"\n"
                f"REPLACEMENT: {replacement}\n"
            )
            if details:
                message += f"\n{details}\n"

            message += (
                f"\n"
                f"For migration guidance, see:\n"
                f"  - KINTSUGI documentation\n"
                f"  - notebooks/MIGRATION_GUIDE.md\n"
                f"{'=' * 60}\n"
            )

            warnings.warn(message, KintsugiDeprecationWarning, stacklevel=2)
            return func(*args, **kwargs)

        return wrapper

    return decorator


def emit_notebook_deprecation_warning(notebook_number: int):
    """
    Emit deprecation warning for a specific notebook.

    Parameters
    ----------
    notebook_number : int
        Notebook number (3 or 5)
    """
    if notebook_number == 3:
        message = """
================================================================================
DEPRECATION NOTICE: Notebook 3 (Signal Isolation)
================================================================================

This notebook is DEPRECATED and will be removed in KINTSUGI v2.0.0.

The Signal Isolation workflow has been replaced with:

1. PURE PYTHON MODULES:
   ```python
   from kintsugi.signal import subtract_autofluorescence

   result = subtract_autofluorescence(signal, blank)
   ```
   - kintsugi.signal: global and weighted autofluorescence subtraction, batch runs
   - kintsugi.denoise: Advanced denoising (N2V, NLM, BM3D)
   - kintsugi.qc: Quality control modules

2. JUPYTER PARAMETER TUNERS:
   - Interactive widgets for manual parameter tuning
   - Real-time preview with stackview

   ```python
   from kintsugi.assist import blank_subtraction_tuner, denoise_tuner

   # Interactive blank subtraction
   result = blank_subtraction_tuner(signal_channel, blank_channel)

   # Get optimized parameters
   params = result.get_params()
   ```

3. AI-ASSISTED WORKFLOW (optional, any MCP client):
   - The KINTSUGI MCP server lets the AI assistant you choose load channels,
     suggest parameters and record what worked for future recommendations
   - Requires the optional mcp extra: pip install -e ".[mcp]"

   ```bash
   kintsugi mcp clients                                  # the supported clients
   kintsugi mcp config /path/to/project --client NAME    # configure yours
   ```

For detailed migration instructions, see:
  notebooks/MIGRATION_GUIDE.md

================================================================================
"""
    elif notebook_number == 5:
        message = """
================================================================================
DEPRECATION NOTICE: Notebook 5 (DL Channel Refinement)
================================================================================

This notebook is DEPRECATED and will be removed in KINTSUGI v2.0.0.

The DL Channel Refinement workflow has been replaced with:

1. QC MODULE:
   ```python
   from kintsugi.qc import ImageQC, CellQC, MarkerQC

   # Image-level quality control
   qc = ImageQC()
   result = qc.assess(image)

   # Marker validation
   marker_qc = MarkerQC()
   result = marker_qc.assess(intensities, marker_name="CD3")
   ```

2. AI-ASSISTED QUALITY ASSESSMENT (optional, any MCP client):
   - The assess_quality and compute_snr tools of the KINTSUGI MCP server
   - Configure your client with: kintsugi mcp config /path/to/project --client NAME

3. PARAMETER LEARNING:
   - Parameters are automatically learned from successful runs
   - Recommendations improve over time based on tissue/marker combinations

For detailed migration instructions, see:
  notebooks/MIGRATION_GUIDE.md

================================================================================
"""
    else:
        return

    warnings.warn(message, KintsugiDeprecationWarning, stacklevel=2)


# Mapping of deprecated functions to their replacements
DEPRECATION_MAP = {
    # Notebook 3 functions
    "ini_params": {
        "replacement": "kintsugi.signal.subtract_autofluorescence",
        "message": (
            "Use kintsugi.signal.subtract_autofluorescence or "
            "kintsugi.assist.blank_subtraction_tuner"
        ),
    },
    "ini_params_alt": {
        "replacement": "kintsugi.signal.subtract_autofluorescence",
        "message": "Use kintsugi.signal.subtract_autofluorescence",
    },
    "ini_params_full": {
        "replacement": "kintsugi.signal.subtract_autofluorescence_dask",
        "message": "Use kintsugi.signal.subtract_autofluorescence_dask for Dask arrays",
    },
    "ini_params_SNR": {
        "replacement": "kintsugi.qc.image_qc.assess_image_quality",
        "message": "Use QC module for SNR computation",
    },
    # Notebook 5 functions
    "ChannelAssessor": {
        "replacement": "kintsugi.qc.ImageQC",
        "message": "Use kintsugi.qc module for quality assessment",
    },
    "BatchChannelProcessor": {
        "replacement": "kintsugi.qc.BatchQC",
        "message": "Use kintsugi.qc.BatchQC for batch processing",
    },
}


def get_replacement_info(func_name: str) -> dict | None:
    """Get replacement information for a deprecated function."""
    return DEPRECATION_MAP.get(func_name)
