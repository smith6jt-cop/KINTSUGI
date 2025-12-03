"""
Kview2 module bridge.

Re-exports the Kview2 visualization module from notebooks/Kview2
for use with the kintsugi package interface.
"""

import sys
from pathlib import Path

# Add notebooks directory to path for imports
_notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
if str(_notebooks_path) not in sys.path:
    sys.path.insert(0, str(_notebooks_path))

# Import and re-export from notebooks/Kview2
try:
    from Kview2 import (
        __version__,
        jupyter_displayable_output,
        insight,
        merge_rgb,
        nop,
        crop,
        annotate,
        interact,
        slice,
        curtain,
        orthogonal,
        side_by_side,
        picker,
        switch,
        create_colormap,
        imshow,
        animate,
        animate_curtain,
        display_range,
        scatterplot,
        grid,
        clusterplot,
        sliceplot,
    )

    __all__ = [
        "__version__",
        "jupyter_displayable_output",
        "insight",
        "merge_rgb",
        "nop",
        "crop",
        "annotate",
        "interact",
        "slice",
        "curtain",
        "orthogonal",
        "side_by_side",
        "picker",
        "switch",
        "create_colormap",
        "imshow",
        "animate",
        "animate_curtain",
        "display_range",
        "scatterplot",
        "grid",
        "clusterplot",
        "sliceplot",
    ]

except ImportError as e:
    import warnings
    warnings.warn(
        f"Failed to import Kview2 module: {e}. "
        "Make sure you're running from the KINTSUGI directory or "
        "have installed the package properly."
    )
    __all__ = []
