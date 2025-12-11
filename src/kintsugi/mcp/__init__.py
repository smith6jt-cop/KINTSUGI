"""
KINTSUGI MCP Server

Model Context Protocol server for Claude Code integration, providing
image processing tools for signal isolation and quality assessment.

Usage:
    # Start the server
    python -m kintsugi.mcp.server

    # Or use the CLI
    kintsugi mcp start

Configuration in Claude Code settings:
    {
        "mcpServers": {
            "kintsugi": {
                "command": "python",
                "args": ["-m", "kintsugi.mcp.server"],
                "cwd": "/path/to/project"
            }
        }
    }
"""

from kintsugi.mcp.server import create_server, run_server

__all__ = ["create_server", "run_server"]
