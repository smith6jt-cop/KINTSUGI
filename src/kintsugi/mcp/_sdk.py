"""Import gate for the MCP Python SDK.

Everything in :mod:`kintsugi.mcp` that touches the SDK imports it from here so
that the rest of the package (the tool registry, the tool handlers, and the
``kintsugi mcp tools`` listing) can be imported without the optional ``claude``
extra, and so the server can tell "``mcp`` is not installed" apart from "an
unsupported ``mcp`` 1.x is installed".

KINTSUGI targets the 2.x line of the official SDK (``mcp>=2.0,<3``), where the
high-level server is :class:`mcp.server.mcpserver.MCPServer` (renamed from
``FastMCP`` in 1.x). The 1.x line is not supported: its ``mcp.server.fastmcp``
import path no longer exists in 2.x and its in-memory test client was removed.
"""

from __future__ import annotations

from importlib import metadata

INSTALL_HINT = 'pip install -e ".[claude]"  (or: kintsugi install claude)'
REQUIRED_SPEC = "mcp>=2.0,<3"

MCP_SDK_VERSION: str | None
try:
    MCP_SDK_VERSION = metadata.version("mcp")
except metadata.PackageNotFoundError:  # pragma: no cover - depends on the environment
    MCP_SDK_VERSION = None

MCP_AVAILABLE = False
MCP_UNSUPPORTED_REASON: str | None = None

try:
    from mcp.server.mcpserver import Context, MCPServer
    from mcp.server.mcpserver.exceptions import ToolError
    from mcp.types import CallToolResult, ImageContent, TextContent, ToolAnnotations

    MCP_AVAILABLE = True
except ImportError:  # pragma: no cover - depends on the environment
    if MCP_SDK_VERSION is None:
        MCP_UNSUPPORTED_REASON = "The 'mcp' package is not installed."
    else:
        MCP_UNSUPPORTED_REASON = (
            f"mcp {MCP_SDK_VERSION} is installed but KINTSUGI requires {REQUIRED_SPEC} "
            "(the 1.x 'mcp.server.fastmcp' API is not supported)."
        )

    class Context:  # type: ignore[no-redef]
        """Placeholder so ``ctx: Context | None`` annotations resolve without the SDK."""

    class ToolError(Exception):  # type: ignore[no-redef]
        """Placeholder for :class:`mcp.server.mcpserver.exceptions.ToolError`."""

    MCPServer = None  # type: ignore[assignment,misc]
    ToolAnnotations = None  # type: ignore[assignment,misc]
    CallToolResult = None  # type: ignore[assignment,misc]
    ImageContent = None  # type: ignore[assignment,misc]
    TextContent = None  # type: ignore[assignment,misc]


def unavailable_message() -> str:
    """Human-readable explanation of why the MCP server cannot start."""
    reason = MCP_UNSUPPORTED_REASON or "The MCP SDK could not be imported."
    return f"{reason} Install with: {INSTALL_HINT}"


__all__ = [
    "INSTALL_HINT",
    "MCP_AVAILABLE",
    "MCP_SDK_VERSION",
    "MCP_UNSUPPORTED_REASON",
    "REQUIRED_SPEC",
    "CallToolResult",
    "Context",
    "ImageContent",
    "MCPServer",
    "TextContent",
    "ToolAnnotations",
    "ToolError",
    "unavailable_message",
]
