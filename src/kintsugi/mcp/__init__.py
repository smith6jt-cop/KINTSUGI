"""
KINTSUGI MCP Server

Optional Model Context Protocol server for any MCP client (AI assistant),
providing image processing tools for signal isolation, quality assessment, and parameter
learning.

Requires the 2.x Python SDK (``mcp>=2.0,<3``); install it with
``pip install -e ".[mcp]"`` or ``kintsugi install mcp`` (never
``pip install kintsugi[mcp]``, an unrelated PyPI package).

Usage:
    # Start the server (stdio transport; stdout carries JSON-RPC only)
    kintsugi mcp start

    # Or, without the console script on PATH
    python -m kintsugi.mcp.server

Client configuration is written only for the client(s) the user chooses
(``kintsugi mcp config <project> --client NAME``, or ``kintsugi init --mcp-client
NAME``; ``kintsugi mcp clients`` lists them). Every entry has the same shape, here
in ``mcpServers`` form; ``command`` is the resolved absolute path to the
``kintsugi`` executable::

    {
        "mcpServers": {
            "kintsugi": {
                "type": "stdio",
                "command": "/path/to/env/bin/kintsugi",
                "args": ["mcp", "start", "--project-dir", "/path/to/project"]
            }
        }
    }

``--project-dir`` tells the server which project ``project_path`` defaults to,
because MCP clients differ in the working directory and environment they give a
server; ``KINTSUGI_PROJECT_DIR`` does the same where a client can set environment
variables. Without either, the working directory is used only if it is a
KINTSUGI project. Some clients have their own command for this; Claude Code's,
for example, is::

    claude mcp add --scope project --transport stdio kintsugi -- \
        kintsugi mcp start --project-dir /path/to/project

The registry is :data:`kintsugi.mcp.tool_specs.TOOL_SPECS` (``kintsugi mcp tools``
lists it); the learning engine is
:class:`kintsugi.assist.parameter_learning.ParameterLearningEngine`.
See ``docs/mcp.md`` for the full guide.
"""

from kintsugi.mcp.server import create_server, run_server

__all__ = ["create_server", "run_server"]
