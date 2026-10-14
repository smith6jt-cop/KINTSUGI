"""
KINTSUGI MCP Server

Exposes the KINTSUGI image-processing tools to Claude Code over the Model
Context Protocol, built on :class:`mcp.server.mcpserver.MCPServer` from the
official ``mcp`` SDK (2.x).

Each tool returns a plain ``dict`` which the SDK emits as **structured content**
(``structuredContent``) alongside the human-readable JSON text block, so Claude
receives machine-readable results rather than only a serialized string. The tool
implementations live in :mod:`kintsugi.mcp.tools` and are registered from the
declarative :data:`kintsugi.mcp.tool_specs.TOOL_SPECS` registry.
"""

from __future__ import annotations

import asyncio
import logging
import sys

from kintsugi.mcp._sdk import MCP_AVAILABLE, MCPServer, unavailable_message
from kintsugi.mcp.tool_specs import TOOL_SPECS

# Configure logging. The stdio transport uses stdout for the JSON-RPC stream,
# so all logging MUST go to stderr to avoid corrupting protocol messages.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("kintsugi.mcp")

SERVER_TITLE = "KINTSUGI image processing"
SERVER_WEBSITE = "https://kintsugi.readthedocs.io/en/latest/mcp.html"

# Sent to the client at initialize time; this is the model's orientation.
SERVER_INSTRUCTIONS = """\
KINTSUGI processes multiplex immunofluorescence images. Tools act on a KINTSUGI
project directory (project_path) and keep loaded images in server memory for the
session.

Typical signal-isolation session:
1. list_channels -> load_channel for the signal channel and its blank
   (autofluorescence) channel.
2. get_image_stats / get_thumbnail to inspect them.
3. suggest_with_learning (learned + heuristic) or suggest_parameters (heuristic).
4. analyze_weighted_subtraction to preview, then subtract_blank
   (method="weighted" is preferred); optionally denoise, clean_background,
   gaussian_subtract.
5. assess_quality / compute_snr, then save_processed.
6. approve_and_learn so the next session starts from what worked.

Failures are returned as tool errors with a next-step hint.
"""


def _server_version() -> str:
    from kintsugi import __version__

    return __version__


def create_server() -> MCPServer:
    """Create and configure the KINTSUGI MCP server."""
    if not MCP_AVAILABLE:
        raise ImportError(unavailable_message())

    server = MCPServer(
        "kintsugi",
        title=SERVER_TITLE,
        instructions=SERVER_INSTRUCTIONS,
        website_url=SERVER_WEBSITE,
        version=_server_version(),
    )
    _register_tools(server)
    return server


def _register_tools(server: MCPServer) -> None:
    """Register every tool in :data:`TOOL_SPECS` with the server.

    Tool handlers stay plain ``async`` functions in :mod:`kintsugi.mcp.tools`
    (SDK-agnostic, so the unit tests can call them directly). They are
    imported lazily here so importing this module does not pull the heavy
    image-processing dependencies.
    """
    from kintsugi.mcp.tools import (
        learning,
        quality_assessment,
        signal_isolation,
        visualization,
        workflow,
    )

    modules = {
        "signal_isolation": signal_isolation,
        "quality_assessment": quality_assessment,
        "visualization": visualization,
        "workflow": workflow,
        "learning": learning,
    }

    for spec in TOOL_SPECS:
        handler = getattr(modules[spec["module"]], spec["attr"])
        # structured_output=True makes the SDK emit the returned dict as
        # ``structuredContent`` (plus the JSON text block) for every tool.
        server.add_tool(
            handler,
            name=spec["name"],
            description=spec["description"],
            structured_output=True,
        )
        _advertise_schema(server, spec["name"], spec["input_schema"])


def _advertise_schema(server: MCPServer, name: str, schema: dict) -> None:
    """Advertise the curated input schema for ``name`` verbatim.

    The SDK derives an input schema from the handler's type hints, which drops
    the per-parameter descriptions and ``enum`` constraints the KINTSUGI tools
    rely on. We overwrite the advertised schema with the curated one from
    :data:`TOOL_SPECS`; argument *validation* still runs from the handler's
    signature, so this only changes the contract clients read.

    This touches an SDK internal (the tool manager), so it is wrapped
    defensively: if a future ``mcp`` release changes the internals the server
    still starts, falling back to the auto-generated schema (and the dedicated
    test in ``tests/test_mcp_tools.py`` flags the regression).
    """
    try:
        tool = server._tool_manager.get_tool(name)
        if tool is not None:
            tool.parameters = schema
    except Exception as exc:  # pragma: no cover - depends on mcp internals
        logger.warning("Could not advertise curated input schema for tool %r: %s", name, exc)


async def run_server() -> None:
    """Run the MCP server using stdio transport.

    Nothing may be written to stdout before or during the run: the stdio
    transport owns it for JSON-RPC frames. Diagnostics go to stderr via logging.
    """
    if not MCP_AVAILABLE:
        print(f"Error: {unavailable_message()}", file=sys.stderr)
        sys.exit(1)

    server = create_server()
    logger.info("Starting KINTSUGI MCP server (%d tools)...", len(TOOL_SPECS))
    await server.run_stdio_async()


def main() -> None:
    """Entry point for the MCP server."""
    asyncio.run(run_server())


if __name__ == "__main__":
    main()
