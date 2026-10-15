"""Session state and project-root helpers shared by the MCP tool modules.

Two things live here so that every tool module, the server's resources and the
tests see the *same* objects without importing each other:

* the in-memory session: ``_loaded_images`` (name -> data, metadata, dtype,
  shape) and ``_processing_history`` (name -> recorded operations), with the
  accessors the tool modules use;
* the project root the server serves. Claude Code starts a project's MCP server
  in the project directory and sets ``CLAUDE_PROJECT_DIR`` in its environment,
  so tools default their ``project_path`` to that directory instead of requiring
  the model to repeat it on every call. The same root anchors the path-safety
  checks.

Nothing here imports ``mcp``; the module is importable without the ``claude``
extra.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_DIR_ENV = "CLAUDE_PROJECT_DIR"

# ---------------------------------------------------------------------------
# Session state (process-wide; one MCP server serves one Claude Code session)
# ---------------------------------------------------------------------------

_loaded_images: dict[str, dict[str, Any]] = {}
_processing_history: dict[str, list[dict[str, Any]]] = {}


def _get_image(name: str) -> Any:
    """Return the data of a loaded image or raise ``ValueError`` if it is not loaded."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


def _store_image(name: str, data: Any, metadata: dict[str, Any] | None = None) -> None:
    """Store (or replace) an image in the session with its metadata."""
    _loaded_images[name] = {
        "data": data,
        "metadata": metadata or {},
        "dtype": str(data.dtype) if hasattr(data, "dtype") else "unknown",
        "shape": data.shape if hasattr(data, "shape") else None,
    }
    if name not in _processing_history:
        _processing_history[name] = []


def _add_history(name: str, operation: str, params: dict[str, Any]) -> None:
    """Append an operation to the processing history of ``name``."""
    _processing_history.setdefault(name, []).append({"operation": operation, "params": params})


def _as_numpy(data: Any) -> np.ndarray:
    """Materialize a loaded image (dask or array-like) as a NumPy array.

    Runs the (possibly heavy) ``compute()`` so callers can push it onto a
    worker thread with :func:`asyncio.to_thread`.
    """
    if hasattr(data, "compute"):
        return data.compute()
    return np.array(data)


def clear_session() -> None:
    """Drop every loaded image and its history."""
    _loaded_images.clear()
    _processing_history.clear()


def estimate_memory_usage() -> dict[str, float]:
    """Estimate the memory held by the loaded images, in MB and GB."""
    total_bytes = 0
    for info in _loaded_images.values():
        data = info["data"]
        if hasattr(data, "nbytes"):
            total_bytes += int(data.nbytes)
        elif info.get("shape"):
            itemsize = {"uint8": 1, "uint16": 2, "float32": 4, "float64": 8}.get(
                info.get("dtype", "uint16"), 2
            )
            total_bytes += int(np.prod(info["shape"])) * itemsize
    return {
        "total_mb": total_bytes / (1024 * 1024),
        "total_gb": total_bytes / (1024 * 1024 * 1024),
    }


# ---------------------------------------------------------------------------
# Project root
# ---------------------------------------------------------------------------


def default_project_root() -> Path:
    """The project directory this server serves.

    ``CLAUDE_PROJECT_DIR`` (set by Claude Code for project-scoped servers) wins;
    otherwise the server's working directory, which is where Claude Code starts
    project servers anyway.
    """
    env = os.environ.get(PROJECT_DIR_ENV)
    return Path(env).resolve() if env else Path.cwd().resolve()


def resolve_project_path(project_path: str | os.PathLike[str] | None) -> str:
    """Return ``project_path`` as a string, defaulting to :func:`default_project_root`."""
    if project_path:
        return str(Path(project_path))
    return str(default_project_root())
