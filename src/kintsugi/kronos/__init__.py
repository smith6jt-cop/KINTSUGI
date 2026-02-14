"""
KINTSUGI KRONOS Integration Module

Integration with KRONOS (Knowledge Representation for Omics via Natively
Orchestrated Segmentation), a panel-agnostic foundation model for spatial
proteomics from Harvard/Brigham and Women's Hospital.

KRONOS was trained via self-supervised learning (DINO) on 47 million
single-marker patches spanning 175 protein markers, 16 tissue types,
8 imaging platforms, and 5 institutions.

Capabilities:
- Extract rich embeddings from multiplex IF images
- Cell and tissue phenotyping with minimal labels
- Artifact detection via learned representations
- Spatial search across tissue patches
- Cross-dataset patient stratification

Prerequisites:
    pip install kronos  # or: pip install -e /path/to/KRONOS
    # Model weights require HuggingFace access: MahmoodLab/KRONOS
    # License: CC-BY-NC-ND-4.0 (academic non-commercial only)

Usage:
    from kintsugi.kronos import KronosEmbedder, MarkerMapper

    # Map KINTSUGI channel names to KRONOS marker IDs
    mapper = MarkerMapper.from_project(project_dir)

    # Extract embeddings from registered images
    embedder = KronosEmbedder()
    results = embedder.embed_project(project_dir, mapper)

    # Cluster and analyze
    from kintsugi.kronos.analysis import cluster_embeddings
    labels = cluster_embeddings(results.patch_embeddings)
"""

from __future__ import annotations

import importlib
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .embeddings import KronosEmbedder
    from .markers import MarkerMapper
    from .model import KronosModel


def _check_kronos_available() -> bool:
    """Check if the kronos package is importable."""
    try:
        import kronos  # noqa: F401

        return True
    except ImportError:
        return False


def _check_torch_available() -> bool:
    """Check if PyTorch is available."""
    try:
        import torch  # noqa: F401

        return True
    except ImportError:
        return False


# Lazy imports to avoid loading heavy deps at import time
_LAZY_IMPORTS = {
    "KronosModel": ".model",
    "KronosEmbedder": ".embeddings",
    "MarkerMapper": ".markers",
    "cluster_embeddings": ".analysis",
    "visualize_embeddings": ".analysis",
    "spatial_search": ".analysis",
}


def __getattr__(name: str):
    if name in _LAZY_IMPORTS:
        module = importlib.import_module(_LAZY_IMPORTS[name], __name__)
        obj = getattr(module, name)
        globals()[name] = obj
        return obj
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


__all__ = [
    "KronosModel",
    "KronosEmbedder",
    "MarkerMapper",
    "cluster_embeddings",
    "visualize_embeddings",
    "spatial_search",
    "_check_kronos_available",
    "_check_torch_available",
]
