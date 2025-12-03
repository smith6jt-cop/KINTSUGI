"""
Kreg module bridge.

Re-exports the Kreg registration module from notebooks/Kreg
for use with the kintsugi package interface.
"""

import sys
from pathlib import Path

# Add notebooks directory to path for imports
_notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
if str(_notebooks_path) not in sys.path:
    sys.path.insert(0, str(_notebooks_path))

# Import and re-export from notebooks/Kreg
try:
    from Kreg import (
        __version__,
        affine_optimizer,
        feature_detectors,
        feature_matcher,
        micro_rigid_registrar,
        non_rigid_registrars,
        preprocessing,
        registration,
        serial_non_rigid,
        serial_rigid,
        slide_io,
        slide_tools,
        valtils,
        viz,
        warp_tools,
    )

    # Main classes for convenience
    from Kreg.registration import Valis

    __all__ = [
        "__version__",
        "affine_optimizer",
        "feature_detectors",
        "feature_matcher",
        "non_rigid_registrars",
        "preprocessing",
        "registration",
        "serial_non_rigid",
        "serial_rigid",
        "slide_io",
        "slide_tools",
        "valtils",
        "viz",
        "warp_tools",
        "micro_rigid_registrar",
        "Valis",
    ]

except ImportError as e:
    import warnings

    warnings.warn(
        f"Failed to import Kreg module: {e}. "
        "Make sure you're running from the KINTSUGI directory or "
        "have installed the package properly.",
        stacklevel=2,
    )
    __all__ = []
