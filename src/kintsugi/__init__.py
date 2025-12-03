"""
KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

Multi-cycle immunofluorescence image registration and analysis pipeline.
Includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution,
extended depth of focus, registration, autofluorescence removal, segmentation,
clustering, and spatial analysis.

Citation:
    Smith, J. A. et al. Protocol for processing and analyzing multiplexed images
    improves lymphatic cell identification and spatial architecture in human tissue.
    STAR Protocols 6, 103976 (2025).
"""

from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("kintsugi")
except PackageNotFoundError:
    __version__ = "1.1.0"  # Fallback for development

__author__ = "Smith JT"
__email__ = "smith6jt@cop.ufl.edu"

# Lazy imports to avoid loading heavy dependencies at startup
def __getattr__(name: str):
    """Lazy loading of submodules."""
    if name == "Kreg":
        from kintsugi import kreg
        return kreg
    elif name == "Kview2":
        from kintsugi import kview2
        return kview2
    elif name == "Kstitch":
        from kintsugi import kstitch
        return kstitch
    elif name == "dl_refinement":
        from kintsugi import dl_refinement
        return dl_refinement
    elif name == "deps":
        from kintsugi import deps
        return deps
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def check_dependencies(verbose: bool = True) -> dict:
    """
    Check all KINTSUGI dependencies and return status.

    Parameters
    ----------
    verbose : bool
        If True, print status messages.

    Returns
    -------
    dict
        Dictionary with dependency status information.
    """
    from kintsugi.deps import DependencyChecker
    checker = DependencyChecker()
    return checker.check_all(verbose=verbose)


def get_config_template() -> dict:
    """
    Get a template configuration dictionary for registration workflow.

    Returns
    -------
    dict
        Template configuration with all available options.
    """
    return {
        "src_dir": "/path/to/source/images",
        "dst_dir": "/path/to/output",
        "reference_image": "cycle1.tif",
        "image_type": "tif",
        "series": 0,
        "max_image_dim_px": 2048,
        "max_processed_image_dim_px": 2048,
        "micro_rigid_registrar_cls": "RigidRegistrar",
        "align_to_reference": True,
        "create_masks": True,
        "resolution_xyu": [0.325, 0.325, "um"],
        "channel_names": [],
        "compose_non_rigid": True,
        "crop_to_overlap": True,
    }


__all__ = [
    "__version__",
    "__author__",
    "__email__",
    "check_dependencies",
    "get_config_template",
    "Kreg",
    "Kview2",
    "Kstitch",
    "dl_refinement",
    "deps",
]
