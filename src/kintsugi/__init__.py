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

import importlib
from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("kintsugi")
except PackageNotFoundError:
    __version__ = "1.2.3"  # Fallback for development

__author__ = "Smith JT"
__email__ = "smith6jt@cop.ufl.edu"

# Mapping of public names to submodule names
_SUBMODULE_MAPPING = {
    "Kreg": "kreg",
    "Kview2": "kview2",
    "Kstitch": "kstitch",
    "deps": "deps",
    "edf": "edf",
    "zarr_io": "zarr_io",
    "kcorrect_gpu": "kcorrect_gpu",
    "parallel_io": "parallel_io",
    "project": "project",
    "signal": "signal",
    "qc": "qc",
    "denoise": "denoise",
    "segment": "segment",
    "gpu": "gpu",
    "rapids": "rapids",
}


# Lazy imports to avoid loading heavy dependencies at startup
def __getattr__(name: str):
    """Lazy loading of submodules.

    If the real submodule fails to import (optional deps missing), return a
    light proxy so attribute checks (e.g. hasattr) succeed while any actual
    use raises an informative ImportError.
    """
    if name in _SUBMODULE_MAPPING:
        submodule_name = _SUBMODULE_MAPPING[name]
        try:
            mod = importlib.import_module(f"{__name__}.{submodule_name}")
            globals()[name] = mod  # Cache the loaded module
            return mod
        except Exception as e:
            # Provide a proxy so the attribute exists (hasattr returns True)
            # but using it raises a helpful error about the missing dependency.
            class _OptionalModuleProxy:
                def __init__(self, mod_name, err):
                    self.__name__ = f"kintsugi.{mod_name}"
                    self._import_error = err

                def __getattr__(self, attr):
                    raise ImportError(
                        f"Optional submodule {submodule_name!r} could not be imported: "
                        f"{self._import_error!s}. Install its optional dependencies to "
                        "enable this functionality."
                    )

                def __repr__(self):
                    return f"<kintsugi optional submodule {submodule_name!r} (failed import)>"

            proxy = _OptionalModuleProxy(submodule_name, e)
            globals()[name] = proxy
            return proxy

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
    # Notebook module bridges
    "Kreg",
    "Kview2",
    "Kstitch",
    # Core modules
    "deps",
    "edf",
    "zarr_io",
    "kcorrect_gpu",
    "parallel_io",
    "project",
    # New processing modules
    "signal",
    "qc",
    "denoise",
    "segment",
    # GPU utilities
    "gpu",
    "rapids",
]
