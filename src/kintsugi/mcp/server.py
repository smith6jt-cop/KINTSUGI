"""
KINTSUGI MCP Server

Exposes image processing tools for Claude Code integration, built on the
high-level FastMCP server from the official ``mcp`` SDK.

Each tool returns a plain ``dict`` which FastMCP emits as **structured content**
(``structuredContent``) alongside the human-readable JSON text block, so Claude
receives machine-readable results rather than only a serialized string. The tool
implementations live in :mod:`kintsugi.mcp.tools` and are registered from the
declarative :data:`kintsugi.mcp.tool_specs.TOOL_SPECS` registry.
"""

from __future__ import annotations

import asyncio
import logging
import sys

from kintsugi.mcp.tool_specs import TOOL_SPECS

try:
    from mcp.server.fastmcp import FastMCP

    MCP_AVAILABLE = True
except ImportError:
    MCP_AVAILABLE = False
    FastMCP = None  # type: ignore[assignment,misc]

# Configure logging. The stdio transport uses stdout for the JSON-RPC stream,
# so all logging MUST go to stderr to avoid corrupting protocol messages.
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    stream=sys.stderr,
)
logger = logging.getLogger("kintsugi.mcp")


def create_server() -> FastMCP:
    """Create and configure the KINTSUGI FastMCP server."""
    if not MCP_AVAILABLE:
        raise ImportError("MCP package not installed. Install with: pip install kintsugi[claude]")

    server = FastMCP("kintsugi")
    _register_tools(server)
    return server


def _register_tools(server: FastMCP) -> None:
    """Register every tool in :data:`TOOL_SPECS` with the FastMCP server.

    Tool handlers stay plain ``async`` functions in :mod:`kintsugi.mcp.tools`
    (FastMCP-agnostic, so the unit tests can call them directly). They are
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
        # structured_output=True makes FastMCP emit the returned dict as
        # ``structuredContent`` (plus the JSON text block) for every tool.
        server.add_tool(
            handler,
            name=spec["name"],
            description=spec["description"],
            structured_output=True,
        )
        _advertise_schema(server, spec["name"], spec["input_schema"])


def _advertise_schema(server: FastMCP, name: str, schema: dict) -> None:
    """Advertise the curated input schema for ``name`` verbatim.

    FastMCP derives an input schema from the handler's type hints, which drops
    the per-parameter descriptions and ``enum`` constraints the KINTSUGI tools
    rely on. We overwrite the advertised schema with the curated one from
    :data:`TOOL_SPECS`; argument *validation* still runs from the handler's
    signature, so this only changes the contract clients read.

    This touches a FastMCP internal (the tool manager), so it is wrapped
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
    """Run the MCP server using stdio transport."""
    if not MCP_AVAILABLE:
        print(
            "Error: MCP package not installed. Install with: pip install kintsugi[claude]",
            file=sys.stderr,
        )
        sys.exit(1)

    server = create_server()
    logger.info("Starting KINTSUGI MCP server...")
    await server.run_stdio_async()


def main() -> None:
    """Entry point for the MCP server."""
    asyncio.run(run_server())


if __name__ == "__main__":
    main()
