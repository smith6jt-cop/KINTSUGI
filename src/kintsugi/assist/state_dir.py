"""Location of a project's KINTSUGI state (learning database, workflow state).

Project state lives in ``<project>/.kintsugi/``. Earlier releases wrote it to
``<project>/.claude_workflow/``; that directory is still read as is, and the first
write moves it to ``.kintsugi/`` so no learned parameters are left behind.

Reads never create anything: :func:`project_state_dir` with ``create=False`` only
looks at what exists.
"""

from __future__ import annotations

import logging
from pathlib import Path

logger = logging.getLogger(__name__)

STATE_DIR_NAME = ".kintsugi"
LEGACY_STATE_DIR_NAME = ".claude_workflow"

# (reason, project) pairs already warned about (each warning once per process).
_warned: set[tuple[str, Path]] = set()


def _warn_once(reason: str, root: Path, message: str, *args: object) -> None:
    key = (reason, root.resolve())
    if key not in _warned:
        _warned.add(key)
        logger.warning(message, *args)


def project_state_dir(project: str | Path, *, create: bool = False) -> Path:
    """Return the state directory of ``project``.

    With ``create=False`` (reads) this returns ``.kintsugi/`` if it exists, else an
    existing legacy ``.claude_workflow/``, else ``.kintsugi/`` (which is not
    created). With ``create=True`` (writes) a legacy-only directory is renamed to
    ``.kintsugi/`` first, then ``.kintsugi/`` is created if needed. If the rename
    fails for a reason other than another process having done it, the legacy
    directory is returned so writes stay next to the existing data.
    """
    root = Path(project)
    current = root / STATE_DIR_NAME
    legacy = root / LEGACY_STATE_DIR_NAME

    if legacy.is_dir() and current.is_dir():
        _warn_once(
            "both",
            root,
            "%s and %s both exist; using %s. Records in %s are no longer read: merge or remove it.",
            current,
            legacy,
            current,
            legacy,
        )
        return current

    if not create:
        return legacy if legacy.is_dir() else current

    if legacy.is_dir():
        try:
            legacy.rename(current)
            logger.info("Moved project state %s to %s", legacy, current)
        except FileNotFoundError:
            pass  # another process moved it first
        except OSError as exc:
            if not current.is_dir():
                _warn_once(
                    "move",
                    root,
                    "Could not move %s to %s (%s); keeping project state in %s",
                    legacy,
                    current,
                    exc,
                    legacy,
                )
                return legacy

    current.mkdir(parents=True, exist_ok=True)
    return current
