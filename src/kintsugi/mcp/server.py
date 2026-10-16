"""
KINTSUGI MCP Server

Exposes the KINTSUGI image-processing tools to any MCP client over the Model
Context Protocol, built on :class:`mcp.server.mcpserver.MCPServer` from the
official ``mcp`` SDK (2.x).

Each tool returns a plain ``dict`` which the SDK emits as **structured content**
(``structuredContent``) alongside the human-readable JSON text block, so the
client receives machine-readable results rather than only a serialized string. The tool
implementations live in :mod:`kintsugi.mcp.tools`; their input schemas are
derived from the handlers' ``Annotated``/``Literal`` signatures (and enforced by
the SDK), while names, titles, descriptions and behaviour annotations come from
the declarative :data:`kintsugi.mcp.tool_specs.TOOL_SPECS` registry.
"""

from __future__ import annotations

import asyncio
import functools
import json
import logging
import sys
from collections.abc import Awaitable, Callable
from typing import Any

from kintsugi.assist._coerce import plain
from kintsugi.mcp._sdk import (
    MCP_AVAILABLE,
    MCPServer,
    ToolAnnotations,
    ToolError,
    unavailable_message,
)
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
project directory; project_path defaults to the project the server was started
for (get_session_state shows it as project_root), so pass it only to work on
another project. Loaded images stay in server memory for the session:
get_session_state shows them, unload_channel frees them. Progress is also kept
in the project (<project>/.kintsugi/) and survives restarts: get_workflow_status
shows every channel's status and the review queue, so start there to resume.

Typical signal-isolation session (also available as the
signal_isolation_workflow prompt):
1. list_channels -> load_channel for the signal channel and its blank
   (autofluorescence) channel.
2. get_image_stats / get_thumbnail to inspect them.
3. suggest_with_learning (learned + heuristic) or suggest_parameters (heuristic).
4. analyze_weighted_subtraction to preview, then subtract_blank
   (method="weighted" is preferred); optionally denoise, clean_background,
   gaussian_subtract.
5. assess_quality / compute_snr, then save_processed.
6. Show the result to the user: approve_and_learn when they approve (so the
   next session starts from what worked), set_review_decision otherwise.

Resources: kintsugi://session/images, kintsugi://session/history/{channel},
kintsugi://learning/stats, kintsugi://project/channels,
kintsugi://project/workflow.
Failures are reported as tool errors (isError) whose text is a JSON object with
an "error" message and, where possible, the next step to take.
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
    _register_resources(server)
    _register_prompts(server)
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
        # The input schema (descriptions, enums, ranges) is derived from the
        # handler's Annotated/Literal signature and enforced by the SDK.
        # structured_output=True makes the SDK emit the returned dict as
        # ``structuredContent`` (plus the JSON text block) for every tool.
        server.add_tool(
            _wrap_handler(handler),
            name=spec["name"],
            title=spec["title"],
            description=spec["description"],
            annotations=ToolAnnotations(**spec["annotations"]),
            structured_output=True,
        )


def _wrap_handler(handler: Callable[..., Awaitable[Any]]) -> Callable[..., Awaitable[Any]]:
    """Report a handler's in-band failure as an MCP tool error.

    Handlers return ``{"error": ...}`` dicts (so direct Python callers and the
    unit tests see plain data). On the wire the MCP spec reports execution
    failures with ``isError: true``; raising :class:`ToolError` here makes the SDK
    do that, with the JSON-encoded error dict as the message so the model still
    gets the next-step hint. ``functools.wraps`` preserves the handler signature
    (and, through ``__wrapped__``, its module globals) so the SDK still derives
    the schema from the handler's own annotations.

    Numpy scalars and arrays anywhere in a result dict are converted to plain
    Python values first: the SDK's serializer rejects them (a ``numpy.bool_`` in
    ``subtract_blank``'s parameters used to fail the whole call).
    """

    @functools.wraps(handler)
    async def wrapper(*args: Any, **kwargs: Any) -> Any:
        result = await handler(*args, **kwargs)
        if isinstance(result, dict):
            result = plain(result)
            if "error" in result:
                raise ToolError(json.dumps(result, default=str))
        return result

    return wrapper


def _register_resources(server: MCPServer) -> None:
    """Expose session and project state as read-only resources.

    Resources are the cheap, parameter-free way for a client to look at state
    without spending a tool call; each is backed by the same function as the
    corresponding tool.
    """
    from kintsugi.mcp.tools import learning, workflow

    @server.resource(
        "kintsugi://session/images",
        name="session_images",
        title="Loaded channels",
        description="Channels held in server memory for this session, with shape, dtype, "
        "metadata, processing step counts and memory usage.",
        mime_type="application/json",
    )
    async def session_images() -> str:
        return json.dumps(await workflow.get_session_state(), default=str)

    @server.resource(
        "kintsugi://session/history/{channel}",
        name="session_history",
        title="Processing history of a channel",
        description="Every operation (with parameters) applied to the named loaded channel.",
        mime_type="application/json",
    )
    async def session_history(channel: str) -> str:
        return json.dumps(await workflow.get_processing_history(channel), default=str)

    @server.resource(
        "kintsugi://learning/stats",
        name="learning_statistics",
        title="Parameter learning statistics",
        description="Record counts, tissues and markers in the user-wide and project "
        "parameter-learning databases.",
        mime_type="application/json",
    )
    async def learning_statistics() -> str:
        from kintsugi.mcp.session import learning_project_root

        return json.dumps(
            await learning.get_learning_statistics(project_path=learning_project_root()),
            default=str,
        )

    @server.resource(
        "kintsugi://project/channels",
        name="project_channels",
        title="Channels in the project",
        description="Image files found in the raw, registered and processed directories "
        "of the project the server was started for.",
        mime_type="application/json",
    )
    async def project_channels() -> str:
        return json.dumps(await workflow.list_channels(), default=str)

    @server.resource(
        "kintsugi://project/workflow",
        name="project_workflow",
        title="Workflow status of the project",
        description="Persistent per-channel status, quality scores and the processing and "
        "review queues of the project the server was started for (survives restarts).",
        mime_type="application/json",
    )
    async def project_workflow() -> str:
        return json.dumps(await workflow.get_workflow_status(), default=str)


def _register_prompts(server: MCPServer) -> None:
    """Register the guided signal-isolation prompt."""
    from kintsugi.mcp import prompts

    server.prompt(
        name=prompts.PROMPT_NAME,
        title=prompts.PROMPT_TITLE,
        description=prompts.PROMPT_DESCRIPTION,
    )(prompts.signal_isolation_workflow)


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
