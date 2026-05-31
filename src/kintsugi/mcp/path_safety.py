"""
Path-safety and input-validation helpers for the KINTSUGI MCP server.

The low-level ``mcp.server.Server`` does not validate tool arguments against
their declared ``inputSchema``, so every argument that reaches a tool handler
must be treated as untrusted (the values originate from an LLM). These helpers
guard against the two highest-risk classes:

1. **Path traversal** — filenames/paths supplied by the model escaping the
   project tree (``../``, absolute paths, symlink escapes). This is the most
   common MCP-server vulnerability class.
2. **Out-of-range / wrong-type arguments** — numeric parameters that could
   exhaust memory or hang (huge filter sizes, unbounded trial counts) and enum
   values the schema cannot enforce on its own.

Helpers raise ``PathSafetyError`` / ``ValueError`` on violation. Tool handlers
catch these and return the existing in-band ``{"error": ...}`` dict, so Claude
sees the failure and can self-correct.
"""

from __future__ import annotations

import math
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
