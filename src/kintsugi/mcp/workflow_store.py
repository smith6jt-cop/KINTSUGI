"""Bridge between the MCP session and the persistent per-project workflow state.

Tools keep working on in-memory images (:mod:`kintsugi.mcp.session`); this module
mirrors what matters across restarts into
:class:`kintsugi.assist.workflow_state.WorkflowState` (``<project>/.kintsugi/``):
which channels were loaded, every processing step, saves, and review decisions.

Tracking is best effort: an image without a project in its metadata is not
tracked, and a failure to persist is logged, never raised into a tool. A channel
is tracked under its *root* (the loaded original a derived image came from), so
``cyc001_CD3_bg_subtracted_denoised`` records its steps on ``cyc001_CD3``.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

from kintsugi.assist.workflow_state import WorkflowState

logger = logging.getLogger("kintsugi.mcp.workflow")

# Session-history operations that change an image and so are recorded as steps.
TRACKED_OPERATIONS = frozenset(
    {
        "subtract_blank",
        "subtract_blank_weighted",
        "denoise",
        "clahe",
        "clean_background",
        "gaussian_subtract",
        "denoise_advanced",
        "optimize_parameters",
        "estimate_background_smo",
    }
)

_stores: dict[Path, WorkflowState] = {}


def workflow_store(project: str | Path | None) -> WorkflowState | None:
    """The cached :class:`WorkflowState` of ``project``.

    None without a project, and None for a directory that is not a KINTSUGI
    project (no ``kintsugi_project.json``): ``project_path`` is a model-supplied
    argument, and state must never be written into an arbitrary directory.
    """
    from kintsugi.mcp.session import is_kintsugi_project

    if not project or not is_kintsugi_project(project):
        return None
    key = Path(project).resolve()
    store = _stores.get(key)
    if store is None:
        store = _stores[key] = WorkflowState(key)
    return store


def clear_workflow_stores() -> None:
    """Forget cached stores (tests, the demo, and a simulated server restart)."""
    _stores.clear()


def image_lineage(metadata: dict[str, Any]) -> tuple[str | None, str | None, str]:
    """``(project, root, cycle)`` recorded in an image's metadata."""
    return metadata.get("project"), metadata.get("root"), metadata.get("cycle") or "unknown"


def _guard(what: str, func, *args, **kwargs) -> Any:
    try:
        return func(*args, **kwargs)
    except Exception as exc:  # never break a tool because state could not be saved
        logger.warning("Could not update workflow state (%s): %s", what, exc)
        return None


def register_channel(metadata: dict[str, Any]) -> None:
    """Register a loaded original (status ``pending`` if new)."""
    project, root, cycle = image_lineage(metadata)
    store = workflow_store(project)
    if store is not None and root:
        _guard("load", store.add_channel, root, cycle=cycle)


def track(
    image: str,
    metadata: dict[str, Any],
    operation: str,
    params: dict[str, Any],
    *,
    quality_after: float | None = None,
    status: str | None = "processing",
) -> None:
    """Record a step for the root of ``image`` in its project's workflow state."""
    project, root, cycle = image_lineage(metadata)
    store = workflow_store(project)
    if store is None or not root:
        return
    _guard(
        operation,
        store.record,
        root,
        operation,
        params,
        image=image,
        cycle=cycle,
        quality_after=quality_after,
        status=status,
    )


def resolve_channel(
    channel: str, project_path: str | Path | None = None
) -> tuple[WorkflowState | None, str]:
    """The store and tracked (root) name for ``channel``.

    A loaded image resolves through its lineage (``project`` and ``root`` in its
    metadata); any other name is taken as a root name in ``project_path`` or the
    project the server was started for. The store is None without a project.
    """
    from kintsugi.mcp.session import _loaded_images, configured_project_root

    metadata = _loaded_images.get(channel, {}).get("metadata", {})
    project = project_path or metadata.get("project") or configured_project_root()
    return workflow_store(project), metadata.get("root") or channel


def record_approval(
    channel: str,
    project_path: str | Path | None,
    *,
    final_params: dict[str, Any],
    quality_score: float | None,
    notes: str = "",
) -> str | None:
    """Mark ``channel`` approved in its project's workflow state; return the tracked name.

    A loaded channel that was never registered (for example one that came from
    another project) is registered first. Returns None when there is no project or
    the channel is unknown there; failures are logged, never raised.
    """
    from kintsugi.mcp.session import _loaded_images

    store, root = resolve_channel(channel, project_path)
    if store is None:
        return None

    def approve() -> str | None:
        if store.get_channel(root) is None:
            if channel not in _loaded_images:
                return None
            cycle = _loaded_images[channel]["metadata"].get("cycle") or "unknown"
            store.add_channel(root, cycle=cycle)
        ok = store.set_user_decision(
            root, "approved", notes, final_params=final_params, quality_score=quality_score
        )
        return root if ok else None

    return _guard("approve", approve)
