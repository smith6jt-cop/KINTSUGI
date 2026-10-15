"""End-to-end smoke test of the MCP server over its real stdio transport.

Spawns ``kintsugi mcp start`` (through ``python -m kintsugi.cli`` so the CLI
code path is exercised, not just the server module), drives it with raw
JSON-RPC frames, and asserts that *every* byte on stdout is protocol. This is
the regression test for the two historical stdout leaks (the Rich start banner
and ``KintsugiProject.load()`` prints), which no in-process test can catch.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import threading
from pathlib import Path

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE
from kintsugi.mcp.tool_specs import TOOL_SPECS
from kintsugi.project import PROJECT_CONFIG_FILE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")


def _frame(obj: dict) -> str:
    return json.dumps(obj) + "\n"


def _read_until_id(proc: subprocess.Popen, *, wanted_id: int, timeout: float) -> list[str]:
    """Read stdout lines until a frame with ``wanted_id`` arrives (or timeout)."""
    import queue
    import time

    assert proc.stdout is not None
    lines: list[str] = []
    q: queue.Queue[str | None] = queue.Queue()

    def pump() -> None:
        for line in proc.stdout:  # type: ignore[union-attr]
            q.put(line)
        q.put(None)

    threading.Thread(target=pump, daemon=True).start()
    deadline = time.monotonic() + timeout
    while True:
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            pytest.fail(f"timed out waiting for response id={wanted_id}; got {lines!r}")
        try:
            line = q.get(timeout=remaining)
        except queue.Empty:
            continue
        if line is None:
            break
        lines.append(line)
        try:
            if json.loads(line).get("id") == wanted_id:
                break
        except json.JSONDecodeError:
            # Returned to the caller, which fails the test with the stderr context.
            break
    return lines


@pytest.fixture
def stdio_project(tmp_path: Path) -> Path:
    """A minimal project with a real kintsugi_project.json and one raw TIFF."""
    tifffile = pytest.importorskip("tifffile")
    (tmp_path / PROJECT_CONFIG_FILE).write_text(
        json.dumps({"config": {"name": "stdio-smoke"}, "kintsugi_path": ""})
    )
    raw = tmp_path / "data" / "raw" / "cyc001"
    raw.mkdir(parents=True)
    tifffile.imwrite(raw / "CD3.tif", np.zeros((8, 8), dtype=np.uint16))
    return tmp_path


_INIT = [
    {
        "jsonrpc": "2.0",
        "id": 1,
        "method": "initialize",
        "params": {
            "protocolVersion": "2025-06-18",
            "capabilities": {},
            "clientInfo": {"name": "kintsugi-stdio-test", "version": "0"},
        },
    },
    {"jsonrpc": "2.0", "method": "notifications/initialized"},
    {"jsonrpc": "2.0", "id": 2, "method": "tools/list", "params": {}},
]


def _converse(extra_args: list[str], cwd: Path, requests: list[dict]) -> dict:
    """Run ``kintsugi mcp start`` with ``extra_args``, send ``requests``, return frames by id.

    The child gets no project environment variable, so only ``--project-dir`` (or
    the working directory) can tell it its project, as with most MCP clients.
    """
    env = {k: v for k, v in os.environ.items() if k not in ("KINTSUGI_PROJECT_DIR",)}
    env["PYTHONUNBUFFERED"] = "1"
    proc = subprocess.Popen(
        [sys.executable, "-m", "kintsugi.cli", "mcp", "start", *extra_args],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        cwd=cwd,
        env=env,
    )
    assert proc.stdin is not None and proc.stdout is not None and proc.stderr is not None
    stderr_chunks: list[str] = []
    stderr_pump = threading.Thread(
        target=lambda: stderr_chunks.extend(proc.stderr),
        daemon=True,  # type: ignore[arg-type]
    )
    stderr_pump.start()
    stdout_lines: list[str] = []
    try:
        proc.stdin.write("".join(_frame(r) for r in requests))
        proc.stdin.flush()
        # Keep stdin open until the last response arrives: the server shuts
        # down on EOF and would otherwise drop the in-flight tools/call.
        stdout_lines = _read_until_id(proc, wanted_id=requests[-1]["id"], timeout=180)
    finally:
        proc.stdin.close()
        try:
            proc.wait(timeout=60)
        except subprocess.TimeoutExpired:
            proc.kill()
            proc.wait()
        stderr_pump.join(timeout=10)

    stdout = "".join(stdout_lines)
    stderr = "".join(stderr_chunks)

    # Every stdout line must be a JSON-RPC frame; stderr is free-form.
    lines = [line for line in stdout.splitlines() if line.strip()]
    assert lines, f"no stdout frames; stderr:\n{stderr}"
    assert stdout.lstrip().startswith("{"), f"stdout starts with non-JSON:\n{stdout[:200]}"
    frames = []
    for line in lines:
        try:
            frames.append(json.loads(line))
        except json.JSONDecodeError:
            pytest.fail(f"non-JSON line on stdout: {line!r}\nstderr:\n{stderr}")
    return {f.get("id"): f for f in frames if "id" in f}


def _list_channels_call(arguments: dict) -> dict:
    return {
        "jsonrpc": "2.0",
        "id": 3,
        "method": "tools/call",
        "params": {"name": "list_channels", "arguments": arguments},
    }


def test_stdio_stream_is_pure_json_rpc(stdio_project: Path, tmp_path_factory):
    """--project-dir alone (no env var, unrelated cwd) gives tools their default project."""
    elsewhere = tmp_path_factory.mktemp("elsewhere")
    by_id = _converse(
        ["--project-dir", str(stdio_project)], elsewhere, [*_INIT, _list_channels_call({})]
    )

    init = by_id[1]["result"]
    assert init["serverInfo"]["name"] == "kintsugi"
    assert init.get("instructions"), "server instructions missing"

    tools = by_id[2]["result"]["tools"]
    assert {t["name"] for t in tools} == {spec["name"] for spec in TOOL_SPECS}

    call = by_id[3]["result"]
    assert not call.get("isError"), call
    assert call["structuredContent"]["status"] == "success"
    assert call["structuredContent"]["project"] == str(stdio_project.resolve())
    names = {ch["channel"] for ch in call["structuredContent"]["channels"]}
    assert "CD3" in names


def test_stdio_without_a_project_lists_tools_and_explains(tmp_path: Path):
    """No --project-dir, no env var, a non-project cwd: tools work, defaults are refused."""
    by_id = _converse([], tmp_path, [*_INIT, _list_channels_call({})])

    assert {t["name"] for t in by_id[2]["result"]["tools"]} == {s["name"] for s in TOOL_SPECS}
    call = by_id[3]["result"]
    assert call.get("isError"), call
    assert "No project configured" in call["content"][0]["text"]
    assert list(tmp_path.iterdir()) == [], "the server wrote into its working directory"
