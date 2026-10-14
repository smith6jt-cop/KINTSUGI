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
    requests = [
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
    ]
    proc = subprocess.Popen(
        [exe, "mcp", "start"],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=tmp,
        env={**os.environ, "CLAUDE_PROJECT_DIR": str(tmp), "PYTHONUNBUFFERED": "1"},
    )
    assert proc.stdin and proc.stdout and proc.stderr
    stderr_chunks: list[str] = []
    threading.Thread(target=lambda: stderr_chunks.extend(proc.stderr), daemon=True).start()  # type: ignore[arg-type]
    proc.stdin.write("".join(json.dumps(r) + "\n" for r in requests))
    proc.stdin.flush()

    lines: list[str] = []
    for line in proc.stdout:
        lines.append(line)
        try:
            if json.loads(line).get("id") == 2:
                break
        except json.JSONDecodeError:
            break
    proc.stdin.close()
    try:
        proc.wait(timeout=60)
    except subprocess.TimeoutExpired:
        proc.kill()
    shutil.rmtree(tmp, ignore_errors=True)

    frames = []
    for line in lines:
        try:
            frames.append(json.loads(line))
        except json.JSONDecodeError:
            fail(f"non-JSON line on stdout: {line!r}\n--- stderr ---\n{''.join(stderr_chunks)}")
    by_id = {f.get("id"): f for f in frames if "id" in f}
    if 1 not in by_id or by_id[1].get("result", {}).get("serverInfo", {}).get("name") != "kintsugi":
        fail(f"no valid initialize response: {by_id.get(1)}")
    tools = by_id.get(2, {}).get("result", {}).get("tools", [])
    if not tools:
        fail(f"tools/list returned no tools: {by_id.get(2)}")
    print(
        f"OK: kintsugi {kintsugi.__version__} from {location.parent}; {len(tools)} tools over stdio"
    )


if __name__ == "__main__":
    main()
