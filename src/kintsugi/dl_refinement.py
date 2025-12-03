"""
Deep Learning Refinement module bridge.

Re-exports the dl_refinement module from notebooks/dl_refinement
for use with the kintsugi package interface.
"""

import sys
from pathlib import Path

# Add notebooks directory to path for imports
_notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
if str(_notebooks_path) not in sys.path:
    sys.path.insert(0, str(_notebooks_path))

# Import and re-export from notebooks/dl_refinement
try:
    from dl_refinement import (
        BatchChannelProcessor,
        BatchReviewInterface,
        ChannelAssessor,
        ChannelQualityResult,
        ChannelReviewInterface,
        HeuristicChannelAssessor,
        StreamingBatchProcessor,
    )

    __all__ = [
        "ChannelAssessor",
        "HeuristicChannelAssessor",
        "ChannelQualityResult",
        "BatchChannelProcessor",
        "StreamingBatchProcessor",
        "ChannelReviewInterface",
        "BatchReviewInterface",
    ]

except ImportError as e:
    import warnings

    warnings.warn(
        f"Failed to import dl_refinement module: {e}. "
        "Make sure you're running from the KINTSUGI directory or "
        "have installed the package properly.",
        stacklevel=2,
    )
    __all__ = []
