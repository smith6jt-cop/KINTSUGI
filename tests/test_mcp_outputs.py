"""Protocol-level tests for visible thumbnails and typed tool outputs.

Through the SDK's in-process client these check exactly what a connected
Claude Code session receives:

* ``get_thumbnail`` delivers a renderable ``image`` content block (not only a
  base64 string buried in JSON) and keeps its metadata as structured content;
  the large ``data_uri`` is opt-in.
* The tools annotated with the ``TypedDict`` shapes in ``kintsugi.mcp.schemas``
  advertise a real ``outputSchema`` (with ``properties``) and their real results
  pass the SDK's output validation, i.e. are not turned into opaque errors.

In-band error results are deliberately not asserted here (their ``isError``
semantics are owned by the server tests).

Run with: pytest tests/test_mcp_outputs.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")

# Tools whose return annotation is a TypedDict from kintsugi.mcp.schemas.
TYPED_OUTPUT_TOOLS = ("get_image_stats", "assess_quality", "suggest_with_learning")

CHANNEL = "cyc001_CD3"
BLANK = "cyc001_blank"
SHAPE = (96, 128)


@pytest.fixture
def seeded_channels():
    """Seed a small signal + blank image pair into the shared MCP session state."""
    from kintsugi.mcp.tools.signal_isolation import (
        _loaded_images,
        _processing_history,
        _store_image,
    )

    _loaded_images.clear()
    _processing_history.clear()

    rng = np.random.default_rng(7)
    signal = rng.integers(100, 500, SHAPE, dtype=np.uint16)
    signal[20:40, 30:60] = 4000  # a bright structure so stats are not flat
    blank = rng.integers(50, 400, SHAPE, dtype=np.uint16)
    _store_image(CHANNEL, signal, {"project": None})
    _store_image(BLANK, blank, {"project": None})

    yield CHANNEL

    _loaded_images.clear()
    _processing_history.clear()


@pytest.fixture
def isolated_learning(tmp_path, monkeypatch):
    """Keep the learning engine's global database out of the real home directory."""
    from kintsugi.mcp.tools import learning

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    learning._learning_engines.clear()
    yield
    learning._learning_engines.clear()


async def _call(name: str, arguments: dict):
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        return await client.call_tool(name, arguments)


# =============================================================================
# get_thumbnail: image content block + structured metadata
# =============================================================================


async def test_thumbnail_arrives_as_image_block(seeded_channels):
    pytest.importorskip("PIL")

    result = await _call("get_thumbnail", {"channel": seeded_channels, "max_size": 64})

    assert not result.is_error
    assert result.content[0].type == "image"
    assert result.content[0].mime_type == "image/png"
    assert result.content[0].data  # base64 payload present
    assert result.content[1].type == "text"

    meta = result.structured_content
    assert meta["status"] == "success"
    assert meta["channel"] == seeded_channels
    assert meta["original_shape"] == list(SHAPE)
    thumb = meta["thumbnail"]
    assert thumb["format"] == "png"
    assert thumb["width"] <= 64 and thumb["height"] <= 64
    assert "data_uri" not in thumb, "data_uri must be opt-in"


async def test_thumbnail_data_uri_is_opt_in(seeded_channels):
    pytest.importorskip("PIL")

    result = await _call(
        "get_thumbnail",
        {"channel": seeded_channels, "max_size": 32, "include_data_uri": True},
    )

    assert not result.is_error
    assert result.content[0].type == "image"
    data_uri = result.structured_content["thumbnail"]["data_uri"]
    assert data_uri.startswith("data:image/png;base64,")
    # The embedded payload is the same bytes as the image block.
    assert data_uri.split(",", 1)[1] == result.content[0].data


async def test_thumbnail_jpeg_sets_matching_mime_type(seeded_channels):
    pytest.importorskip("PIL")

    result = await _call("get_thumbnail", {"channel": seeded_channels, "format": "jpg"})

    assert not result.is_error
    assert result.content[0].type == "image"
    assert result.content[0].mime_type == "image/jpeg"
    assert result.structured_content["thumbnail"]["format"] == "jpeg"


def test_render_thumbnail_helper_is_pure():
    """The encode helper works on a bare array without any session state."""
    pytest.importorskip("PIL")
    from kintsugi.mcp.tools.visualization import _render_thumbnail

    data = np.zeros((300, 200), dtype=np.uint16)
    data[100:200, 50:150] = 10_000

    encoded, meta = _render_thumbnail(data, 100, True, "png")

    assert encoded.startswith(b"\x89PNG")
    assert meta == {"width": 66, "height": 100, "format": "png"}


# =============================================================================
# Typed outputs: advertised schema + validated real results
# =============================================================================


async def test_typed_tools_advertise_output_properties():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}

    for name in TYPED_OUTPUT_TOOLS:
        schema = tools[name].output_schema
        assert schema is not None, name
        props = schema.get("properties")
        assert props, f"{name} should advertise properties, got {schema}"
        assert "status" in props, name
        assert "error" in props, name

    assert "statistics" in tools["get_image_stats"].output_schema["properties"]
    assert "quality_score" in tools["assess_quality"].output_schema["properties"]
    assert "operations" in tools["suggest_with_learning"].output_schema["properties"]


async def test_get_image_stats_round_trip(seeded_channels):
    result = await _call("get_image_stats", {"channel": seeded_channels})

    assert not result.is_error, result.content[0].text
    meta = result.structured_content
    assert meta["status"] == "success"
    assert meta["channel"] == seeded_channels
    stats = meta["statistics"]
    assert isinstance(stats["mean"], float)
    assert stats["shape"] == list(SHAPE)
    assert isinstance(stats["histogram"]["counts"][0], int)
    # Keys beyond the schema's named ones survive (ImageStatistics allows extras).
    assert isinstance(stats["zero_pixels"], int)


async def test_assess_quality_round_trip(seeded_channels):
    result = await _call("assess_quality", {"channel": seeded_channels})

    assert not result.is_error, result.content[0].text
    meta = result.structured_content
    assert meta["status"] == "success"
    assert isinstance(meta["quality_score"], float)
    assert 0.0 <= meta["quality_score"] <= 1.0
    assert isinstance(meta["requires_review"], bool)
    assert isinstance(meta["metrics"]["snr"], float)
    assert isinstance(meta["recommendations"], list)
    assert meta["threshold_used"] == 0.85


async def test_suggest_with_learning_round_trip(seeded_channels, isolated_learning, tmp_path):
    result = await _call(
        "suggest_with_learning",
        {
            "channel": seeded_channels,
            "tissue_type": "tonsil",
            "marker_name": "CD3",
            "blank_channel": BLANK,
            "project_path": str(tmp_path / "proj"),
        },
    )

    assert not result.is_error, result.content[0].text
    meta = result.structured_content
    assert meta["status"] == "success"
    assert meta["channel"] == seeded_channels
    assert meta["marker_name_normalized"]
    assert isinstance(meta["image_characteristics"]["snr"], float)
    assert set(meta["operations"]) >= {"blank_subtraction", "denoise", "clahe"}
    op = meta["operations"]["blank_subtraction"]
    assert isinstance(op["confidence"], float)
    assert isinstance(op["learned_records"], int)
    assert op["source"] in {"heuristic", "blended", "learned"}
    assert isinstance(op["recommended_parameters"], dict)
    assert meta["suggested_order"] and all(isinstance(s, str) for s in meta["suggested_order"])
