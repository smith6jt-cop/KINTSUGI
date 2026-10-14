"""Session-level helpers shared by the MCP tool modules.

Claude Code starts a project's MCP server in the project directory and sets
``CLAUDE_PROJECT_DIR`` in its environment, so tools can default their
``project_path`` to that directory instead of requiring the model to repeat it
on every call. The same root anchors the path-safety checks.
"""

from __future__ import annotations

import os
from pathlib import Path

PROJECT_DIR_ENV = "CLAUDE_PROJECT_DIR"


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
