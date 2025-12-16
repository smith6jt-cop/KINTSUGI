"""
Kview module bridge.

Re-exports the Kview visualization module from notebooks/Kview
for use with the kintsugi package interface.
"""

import sys
from pathlib import Path

# Add notebooks directory to path for imports
_notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
if str(_notebooks_path) not in sys.path:
    sys.path.insert(0, str(_notebooks_path))

# Import and re-export from notebooks/Kview
try:
    from Kview import (
        __version__,
        animate,
        animate_curtain,
        annotate,
        clusterplot,
        create_colormap,
        crop,
        curtain,
        display_range,
        grid,
        imshow,
        insight,
        interact,
        jupyter_displayable_output,
        merge_rgb,
        nop,
        orthogonal,
        picker,
        scatterplot,
        side_by_side,
        slice,
        sliceplot,
        switch,
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
        f"Failed to import Kview module: {e}. "
        "Make sure you're running from the KINTSUGI directory or "
        "have installed the package properly.",
        stacklevel=2,
    )
    __all__ = []
