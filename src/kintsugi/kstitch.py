"""
Kstitch module bridge.

Re-exports the Kstitch stitching module from notebooks/Kstitch
for use with the kintsugi package interface.
"""

import sys
from pathlib import Path

# Add notebooks directory to path for imports
_notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
if str(_notebooks_path) not in sys.path:
    sys.path.insert(0, str(_notebooks_path))

# Import and re-export from notebooks/Kstitch
try:
    from Kstitch import stitch_images

    __all__ = ["stitch_images"]

except ImportError as e:
    import warnings
    warnings.warn(
        f"Failed to import Kstitch module: {e}. "
        "Make sure you're running from the KINTSUGI directory or "
        "have installed the package properly."
    )
    __all__ = []
