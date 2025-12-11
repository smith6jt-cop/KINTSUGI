"""
Deep Learning Refinement module (DEPRECATED).

This module has been replaced by the new `kintsugi.qc` module.
Please update your imports:

    # Old (deprecated)
    from kintsugi import dl_refinement
    assessor = dl_refinement.ChannelAssessor()

    # New
    from kintsugi.qc import ImageQC, CellQC, MarkerQC
    qc = ImageQC()
    result = qc.assess(image)

See notebooks/MIGRATION_GUIDE.md for detailed migration instructions.
"""

import warnings

warnings.warn(
    "The 'dl_refinement' module is deprecated and has been removed. "
    "Use 'kintsugi.qc' instead. See notebooks/MIGRATION_GUIDE.md for migration guidance.",
    DeprecationWarning,
    stacklevel=2,
)

# Provide helpful error messages for old imports
def __getattr__(name):
    """Provide helpful error messages for deprecated classes."""
    _migration_map = {
        "ChannelAssessor": "Use kintsugi.qc.ImageQC instead",
        "HeuristicChannelAssessor": "Use kintsugi.qc.ImageQC instead",
        "BatchChannelProcessor": "Use kintsugi.qc.BatchQC instead",
        "ChannelQualityResult": "Use kintsugi.qc.QCResult instead",
        "ChannelReviewInterface": "Use the 3_Signal_Isolation_QC.ipynb notebook",
        "BatchReviewInterface": "Use the 3_Signal_Isolation_QC.ipynb notebook",
        "StreamingBatchProcessor": "Use kintsugi.qc.BatchQC instead",
    }

    if name in _migration_map:
        raise ImportError(
            f"'{name}' has been removed from dl_refinement. "
            f"{_migration_map[name]}. "
            "See notebooks/MIGRATION_GUIDE.md for migration guidance."
        )

    raise AttributeError(f"module 'kintsugi.dl_refinement' has no attribute '{name}'")

__all__ = []
