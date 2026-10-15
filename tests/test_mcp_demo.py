"""The recorded MCP session (evals/mcp/demo.py) runs end to end without errors."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")

EVALS = Path(__file__).resolve().parents[1] / "evals" / "mcp"


@pytest.fixture(scope="module")
def transcript() -> str:
    pytest.importorskip("dask_image")
    pytest.importorskip("tifffile")
    sys.path.insert(0, str(EVALS))
    try:
        import demo
    finally:
        sys.path.remove(str(EVALS))
    return demo.run()


def test_every_step_succeeds(transcript):
    assert "isError: true" not in transcript, transcript
    for tool in ("load_channel", "subtract_blank", "assess_quality", "approve_and_learn"):
        assert f"`{tool}`" in transcript


def test_thumbnail_is_an_image_block(transcript):
    assert "content: image block (image/png" in transcript


def test_transcript_is_deterministic_and_path_free(transcript):
    sys.path.insert(0, str(EVALS))
    try:
        import demo
    finally:
        sys.path.remove(str(EVALS))
    assert "/tmp/" not in transcript and "kintsugi-mcp-demo" not in transcript
    assert demo.run() == transcript


def test_docs_transcript_is_current(transcript):
    sys.path.insert(0, str(EVALS))
    try:
        import demo
    finally:
        sys.path.remove(str(EVALS))
    docs = Path(__file__).resolve().parents[1] / "docs" / "mcp.md"
    assert demo.is_current(docs, transcript), (
        "docs/mcp.md transcript is stale; run: python evals/mcp/demo.py --write docs/mcp.md"
    )
