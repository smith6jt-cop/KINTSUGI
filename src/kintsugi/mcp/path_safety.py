"""
Path-safety and input-validation helpers for the KINTSUGI MCP server.

The SDK validates argument *types*, enums and numeric ranges from each handler's
``Annotated``/``Literal`` signature before the handler runs. Every argument is
still untrusted (the values originate from an LLM), and some checks cannot be
expressed as a schema, so these helpers guard against the two highest-risk
classes:

1. **Path traversal** — filenames/paths supplied by the model escaping the
   project tree (``../``, absolute paths, symlink escapes). This is the most
   common MCP-server vulnerability class.
2. **Out-of-range / wrong-type arguments** — numeric parameters that could
   exhaust memory or hang (huge filter sizes, unbounded trial counts). Pydantic
   runs in lax mode (``True`` -> ``1``, ``"7"`` -> ``7``), so the clamps below
   stay as the strict, in-body line of defence.

A third class, **deserialization of model files** chosen by the model
(``torch.load``/``joblib.load``), is handled by :func:`allowed_model_roots` and
:func:`ensure_within_any`: a checkpoint path must live under the project or the
user's model directory before anything is loaded.

Helpers raise ``PathSafetyError`` / ``ValueError`` on violation. Tool handlers
catch these and return the existing in-band ``{"error": ...}`` dict, so Claude
sees the failure and can self-correct.
"""

from __future__ import annotations

import math
import os
from collections.abc import Iterable
from pathlib import Path


class PathSafetyError(ValueError):
    """Raised when a path or filename fails a safety check."""


def safe_filename(name: str, param: str = "filename") -> str:
    """Validate that ``name`` is a bare filename with no directory component.

    Rejects path separators (``/`` or ``\\``), parent references (``..``),
    absolute paths, and null bytes. Returns ``name`` unchanged when valid.

    >>> safe_filename("CD3_20260101")
    'CD3_20260101'
    >>> safe_filename("../evil")            # doctest: +SKIP
    PathSafetyError
    """
    if not isinstance(name, str) or not name.strip():
        raise PathSafetyError(f"{param} must be a non-empty string")
    if "\x00" in name:
        raise PathSafetyError(f"{param} must not contain null bytes")
    if "\\" in name:
        raise PathSafetyError(f"{param} must not contain path separators: {name!r}")
    if name in (".", ".."):
        raise PathSafetyError(f"{param} is not a valid filename: {name!r}")
    # Path(...).name strips any directory component; if the result differs the
    # input contained a separator or traversal segment (e.g. "a/b", "../x",
    # "/etc/passwd" all reduce to a different ``.name``).
    if Path(name).name != name:
        raise PathSafetyError(f"{param} must not contain path separators: {name!r}")
    return name


def ensure_within(path: str | Path, root: str | Path, param: str = "path") -> Path:
    """Resolve ``path`` and require it to stay inside ``root``.

    ``Path.resolve()`` collapses ``..`` segments *and* follows symlinks, so this
    also rejects symlink escapes (the EscapeRoute / CVE-2025-53109 class). The
    resolved path is returned for downstream use.
    """
    root_resolved = Path(root).resolve()
    resolved = Path(path).resolve()
    if not resolved.is_relative_to(root_resolved):
        raise PathSafetyError(f"{param} {resolved} escapes the allowed directory {root_resolved}")
    return resolved


def ensure_within_any(path: str | Path, roots: Iterable[str | Path], param: str = "path") -> Path:
    """Resolve ``path`` and require it to stay inside at least one of ``roots``."""
    resolved = Path(path).resolve()
    root_list = [Path(r).resolve() for r in roots]
    for root in root_list:
        if resolved.is_relative_to(root):
            return resolved
    allowed = ", ".join(str(r) for r in root_list) or "(none)"
    raise PathSafetyError(f"{param} {resolved} is outside the allowed directories: {allowed}")


USER_MODEL_DIR_ENV = "KINTSUGI_MODEL_DIR"
PROJECT_MARKER = "kintsugi_project.json"


def is_kintsugi_project(path: str | Path) -> bool:
    """True when ``path`` is a KINTSUGI project directory (has the project marker file)."""
    try:
        return (Path(path).resolve() / PROJECT_MARKER).is_file()
    except OSError:
        return False


def allowed_model_roots(project_path: str | Path | None = None) -> list[Path]:
    """Directories a tool may load model checkpoints from.

    Loading a checkpoint executes a deserializer on a file the model chose, so the
    file must come from somewhere the user controls: the Claude Code project
    (``CLAUDE_PROJECT_DIR``) and the user's model directory (``~/.kintsugi/models``
    or ``$KINTSUGI_MODEL_DIR``). An explicit ``project_path`` is itself a
    model-supplied argument, so it only counts when it really is a KINTSUGI project
    (contains ``kintsugi_project.json``); ``project_path="/"`` or any other bare
    directory does not widen the sandbox.
    """
    roots: list[Path] = []
    if project_path and is_kintsugi_project(project_path):
        roots.append(Path(project_path))
    claude_project = os.environ.get("CLAUDE_PROJECT_DIR")
    if claude_project:
        roots.append(Path(claude_project))
    user_dir = os.environ.get(USER_MODEL_DIR_ENV)
    roots.append(Path(user_dir) if user_dir else Path.home() / ".kintsugi" / "models")
    return roots


def validate_choice(value: object, choices: Iterable[object], param: str) -> object:
    """Require ``value`` to be one of ``choices``; raise ``ValueError`` otherwise."""
    choice_set = set(choices)
    if value not in choice_set:
        allowed = ", ".join(repr(c) for c in sorted(choice_set, key=str))
        raise ValueError(f"{param} must be one of [{allowed}], got {value!r}")
    return value


def clamp_int(
    value: int,
    param: str,
    *,
    minimum: int | None = None,
    maximum: int | None = None,
) -> int:
    """Validate that ``value`` is an int within ``[minimum, maximum]``."""
    # bool is a subclass of int; reject it explicitly to catch type confusion.
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{param} must be an integer, got {value!r}")
    if minimum is not None and value < minimum:
        raise ValueError(f"{param} must be >= {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{param} must be <= {maximum}, got {value}")
    return value


def clamp_float(
    value: float,
    param: str,
    *,
    minimum: float | None = None,
    maximum: float | None = None,
) -> float:
    """Validate that ``value`` is a finite real number within ``[minimum, maximum]``."""
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f"{param} must be a number, got {value!r}")
    value = float(value)
    # NaN slips past range comparisons (both are False) and inf can bypass an
    # absent bound, so reject non-finite values explicitly.
    if not math.isfinite(value):
        raise ValueError(f"{param} must be a finite number, got {value!r}")
    if minimum is not None and value < minimum:
        raise ValueError(f"{param} must be >= {minimum}, got {value}")
    if maximum is not None and value > maximum:
        raise ValueError(f"{param} must be <= {maximum}, got {value}")
    return value
