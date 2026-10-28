"""Session state and project-root helpers shared by the MCP tool modules.

Two things live here so that every tool module, the server's resources and the
tests see the *same* objects without importing each other:

* the in-memory session: ``_loaded_images`` (name -> data, metadata, dtype,
  shape) and ``_processing_history`` (name -> recorded operations), with the
  accessors the tool modules use;
* the project root the server serves, so tools can default ``project_path``
  instead of requiring the model to repeat it on every call. MCP clients differ
  in what they tell a server about the workspace (only some set a working
  directory, few export an environment variable), so the root is configured
  explicitly: ``kintsugi mcp start --project-dir PATH`` (what the generated client
  configs pass), else ``KINTSUGI_PROJECT_DIR``, else the working directory when it
  is a KINTSUGI project. Without any of these there is no default project: tools
  that need one say so instead of guessing, and nothing is written under an
  arbitrary working directory. The same root anchors the path-safety checks.

Nothing here imports ``mcp``; the module is importable without the ``mcp``
extra.
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import Any

import numpy as np

PROJECT_DIR_ENV = "KINTSUGI_PROJECT_DIR"
PROJECT_MARKER = "kintsugi_project.json"

#: Shared description of every tool's optional ``project_path`` argument.
PROJECT_PATH_HELP = (
    "KINTSUGI project directory; None = the project the server was started for "
    "(--project-dir / KINTSUGI_PROJECT_DIR, or the working directory if it is a "
    "KINTSUGI project)"
)

NO_PROJECT_MESSAGE = (
    "No project configured: pass project_path, or start the server with "
    "`kintsugi mcp start --project-dir PATH` (or set KINTSUGI_PROJECT_DIR)"
)

# ---------------------------------------------------------------------------
# Session state (process-wide; one MCP server serves one client session)
# ---------------------------------------------------------------------------

_loaded_images: dict[str, dict[str, Any]] = {}
_processing_history: dict[str, list[dict[str, Any]]] = {}


def _get_image(name: str) -> Any:
    """Return the data of a loaded image or raise ``ValueError`` if it is not loaded."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


def _store_image(name: str, data: Any, metadata: dict[str, Any] | None = None) -> None:
    """Store (or replace) an image in the session with its metadata.

    Lineage is filled in so every image knows the project it belongs to: a derived
    image (``metadata["source"]`` or ``source_channel`` naming a loaded image)
    inherits that image's ``project``, ``cycle``, ``channel`` and ``root``; an image
    read from disk (``source_path``) is its own ``root``.
    """
    metadata = dict(metadata or {})
    parent_name = metadata.get("source") or metadata.get("source_channel")
    if parent_name in _loaded_images:
        parent = _loaded_images[parent_name]["metadata"]
        for key in ("project", "cycle", "channel", "root"):
            if key in parent and key not in metadata:
                metadata[key] = parent[key]
    elif "source_path" in metadata:
        metadata.setdefault("root", name)
    _loaded_images[name] = {
        "data": data,
        "metadata": metadata,
        "dtype": str(data.dtype) if hasattr(data, "dtype") else "unknown",
        "shape": data.shape if hasattr(data, "shape") else None,
    }
    if name not in _processing_history:
        _processing_history[name] = []


def _add_history(name: str, operation: str, params: dict[str, Any]) -> None:
    """Append an operation to the processing history of ``name``.

    Operations that change an image are also recorded in the project's persistent
    workflow state (:mod:`kintsugi.mcp.workflow_store`).
    """
    _processing_history.setdefault(name, []).append({"operation": operation, "params": params})
    from kintsugi.mcp.workflow_store import TRACKED_OPERATIONS, track

    if operation in TRACKED_OPERATIONS and name in _loaded_images:
        track(name, _loaded_images[name]["metadata"], operation, params)


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


class ProjectNotConfigured(ValueError):
    """No ``project_path`` was given and the server has no configured project."""

    def __init__(self, message: str = NO_PROJECT_MESSAGE) -> None:
        super().__init__(message)


_explicit_root: Path | None = None


def set_project_root(path: str | os.PathLike[str] | None) -> None:
    """Set (or with ``None`` clear) the project given by ``mcp start --project-dir``."""
    global _explicit_root
    _explicit_root = Path(path).resolve() if path else None


def is_kintsugi_project(path: str | os.PathLike[str]) -> bool:
    """True when ``path`` is a KINTSUGI project directory (has the project marker file)."""
    try:
        return (Path(path).resolve() / PROJECT_MARKER).is_file()
    except OSError:
        return False


def configured_project_root() -> Path | None:
    """The project directory this server serves, or ``None`` if there is none.

    Precedence: the ``--project-dir`` given to ``kintsugi mcp start``, then a
    non-empty ``KINTSUGI_PROJECT_DIR``, then the working directory but only when
    it holds ``kintsugi_project.json`` (clients start servers in all kinds of
    directories, including ``/`` and the user's home).
    """
    if _explicit_root is not None:
        return _explicit_root
    env = os.environ.get(PROJECT_DIR_ENV, "").strip()
    if env:
        return Path(env).resolve()
    cwd = Path.cwd()
    return cwd.resolve() if is_kintsugi_project(cwd) else None


def resolve_project_path(project_path: str | os.PathLike[str] | None) -> str:
    """Return ``project_path``, defaulting to the configured project.

    Raises :class:`ProjectNotConfigured` when neither is available; tools turn
    that into an in-band error the model can act on.
    """
    if project_path:
        return str(Path(project_path))
    root = configured_project_root()
    if root is None:
        raise ProjectNotConfigured()
    return str(root)


def learning_project_root() -> str | None:
    """Project whose learning database the learning tools use by default.

    The configured project when it is a KINTSUGI project, else ``None`` (the
    user-wide database only).
    """
    root = configured_project_root()
    return str(root) if root is not None and is_kintsugi_project(root) else None
