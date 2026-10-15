#!/usr/bin/env python
"""Smoke-test an *installed* KINTSUGI over the real MCP stdio transport.

Unlike ``tests/test_mcp_stdio.py`` this has no pytest, no conftest and no
``src/`` on ``sys.path``: run it from any directory outside the checkout against
the installed wheel and it fails loudly if the SDK gate is unavailable, if the
package being exercised is not the installed one, or if the server's stdout
carries anything but JSON-RPC.

    python scripts/mcp_stdio_smoke.py
"""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
import tempfile
import threading
from pathlib import Path


def fail(message: str) -> None:
    print(f"FAIL: {message}", file=sys.stderr)
    raise SystemExit(1)


def main() -> None:
    import kintsugi
    from kintsugi.mcp._sdk import MCP_AVAILABLE, unavailable_message

    location = Path(kintsugi.__file__).resolve()
    if "site-packages" not in location.parts and "dist-packages" not in location.parts:
        fail(f"kintsugi imported from {location}, not from an installed package")
    if not MCP_AVAILABLE:
        fail(unavailable_message())

    exe = shutil.which("kintsugi")
    if not exe:
        fail("the kintsugi console script is not on PATH")

    tmp = Path(tempfile.mkdtemp(prefix="kintsugi-smoke-"))
    try:
        project = tmp / "project"
        project.mkdir()
        (project / "kintsugi_project.json").write_text(
            json.dumps({"config": {"name": "smoke"}, "kintsugi_path": ""})
        )
        bare = tmp / "bare"
        bare.mkdir()

        # 1. With --project-dir (what generated client configs pass), started from a
        #    directory that is not the project: tools default to that project.
        by_id = converse(exe, ["--project-dir", str(project)], bare)
        if by_id[1].get("result", {}).get("serverInfo", {}).get("name") != "kintsugi":
            fail(f"no valid initialize response: {by_id.get(1)}")
        tools = by_id.get(2, {}).get("result", {}).get("tools", [])
        if not tools:
            fail(f"tools/list returned no tools: {by_id.get(2)}")
        listed = by_id.get(3, {}).get("result", {})
        if listed.get("isError") or listed.get("structuredContent", {}).get("project") != str(
            project.resolve()
        ):
            fail(f"list_channels did not default to --project-dir: {listed}")

        # 2. No project at all: the server still starts and lists tools, and a tool
        #    that needs a project says so instead of writing into its cwd.
        by_id = converse(exe, [], bare)
        if not by_id.get(2, {}).get("result", {}).get("tools"):
            fail(f"tools/list failed without a project: {by_id.get(2)}")
        refused = by_id.get(3, {}).get("result", {})
        if not refused.get("isError") or "No project configured" not in str(refused):
            fail(f"list_channels without a project should be an in-band error: {refused}")
        if any(bare.iterdir()):
            fail(f"the server wrote into its working directory: {list(bare.iterdir())}")
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print(
        f"OK: kintsugi {kintsugi.__version__} from {location.parent}; {len(tools)} tools over "
        "stdio, --project-dir honoured, no-project case refused cleanly"
    )


REQUESTS = [
    {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2025-06-18",
            "capabilities": {},
            "clientInfo": {"name": "kintsugi-smoke", "version": "0"},
        },
    },
    {"jsonrpc": "2.0", "method": "notifications/initialized"},
    {"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}},
    {
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {"name": "list_channels", "arguments": {}},
    },
]


def converse(exe: str, extra_args: list[str], cwd: Path) -> dict:
    """Run ``kintsugi mcp start`` with ``extra_args``; return JSON-RPC frames by id."""
    env = {k: v for k, v in os.environ.items() if k != "KINTSUGI_PROJECT_DIR"}
    env["PYTHONUNBUFFERED"] = "1"
    proc = subprocess.Popen(
        [exe, "mcp", "start", *extra_args],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=cwd,
        env=env,
    )
    assert proc.stdin and proc.stdout and proc.stderr
    stderr_chunks: list[str] = []
    threading.Thread(target=lambda: stderr_chunks.extend(proc.stderr), daemon=True).start()  # type: ignore[arg-type]
    proc.stdin.write("".join(json.dumps(r) + "\n" for r in REQUESTS))
    proc.stdin.flush()

    lines: list[str] = []
    for line in proc.stdout:
        lines.append(line)
        try:
            if json.loads(line).get("id") == 3:
                break
        except json.JSONDecodeError:
            break
    proc.stdin.close()
    try:
        proc.wait(timeout=60)
    except subprocess.TimeoutExpired:
        proc.kill()

    frames = []
    for line in lines:
        try:
            frames.append(json.loads(line))
        except json.JSONDecodeError:
            fail(f"non-JSON line on stdout: {line!r}\n--- stderr ---\n{''.join(stderr_chunks)}")
    by_id = {f.get("id"): f for f in frames if "id" in f}
    if 1 not in by_id:
        fail(f"no initialize response; stderr:\n{''.join(stderr_chunks)}")
    return by_id


if __name__ == "__main__":
    main()
