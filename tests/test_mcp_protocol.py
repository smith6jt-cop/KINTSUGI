"""Protocol-level tests of the KINTSUGI MCP server through the SDK's in-process client.

These exercise exactly what a connected Claude Code session sees: the advertised
tool list (names, titles, behaviour annotations, input and output schemas), the
server instructions, and call results (structured content, validation errors).
"""

from __future__ import annotations

import inspect
import typing

import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE
from kintsugi.mcp.tool_specs import TOOL_SPECS

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")

HINT_KEYS = ("readOnlyHint", "destructiveHint", "idempotentHint", "openWorldHint")


@pytest.fixture(scope="module")
def surface():
    """Everything the client learns at initialize + tools/list, read once per module."""
    import asyncio

    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async def _collect():
        async with Client(create_server()) as client:
            tools = (await client.list_tools()).tools
            return {
                "tools": {t.name: t for t in tools},
                "instructions": client.instructions,
                "server_info": client.server_info,
            }

    return asyncio.run(_collect())


def _handler(spec):
    import importlib

    module = importlib.import_module(f"kintsugi.mcp.tools.{spec['module']}")
    return getattr(module, spec["attr"])


def test_initialize_exposes_identity_and_instructions(surface):
    assert surface["server_info"].name == "kintsugi"
    assert surface["server_info"].version
    assert "load_channel" in surface["instructions"]
    assert "approve_and_learn" in surface["instructions"]


def test_tool_list_matches_registry(surface):
    assert set(surface["tools"]) == {spec["name"] for spec in TOOL_SPECS}


@pytest.mark.parametrize("spec", TOOL_SPECS, ids=lambda s: s["name"])
def test_tool_metadata_on_the_wire(surface, spec):
    tool = surface["tools"][spec["name"]]
    assert tool.title == spec["title"]
    assert tool.description == spec["description"]
    assert tool.output_schema is not None, "structured output must be declared"
    hints = tool.annotations.model_dump(by_alias=True, exclude_none=True)
    assert {k: hints.get(k) for k in HINT_KEYS} == spec["annotations"]


@pytest.mark.parametrize("spec", TOOL_SPECS, ids=lambda s: s["name"])
def test_input_schema_matches_signature(surface, spec):
    """The advertised schema is derived from the handler and fully described."""
    tool = surface["tools"][spec["name"]]
    schema = tool.input_schema
    props = schema.get("properties", {})

    handler = _handler(spec)
    params = [
        name
        for name, p in inspect.signature(handler).parameters.items()
        if name != "ctx" and p.kind is not inspect.Parameter.VAR_KEYWORD
    ]
    assert set(props) == set(params), f"{spec['name']}: schema/signature drift"

    required = set(schema.get("required", []))
    sig = inspect.signature(handler)
    for name in params:
        prop = props[name]
        assert prop.get("description", "").strip(), f"{spec['name']}.{name} lacks a description"
        has_default = sig.parameters[name].default is not inspect.Parameter.empty
        assert (name in required) == (not has_default), f"{spec['name']}.{name} required-ness"
        if has_default and sig.parameters[name].default is not None:
            assert prop.get("default") == sig.parameters[name].default, f"{spec['name']}.{name}"

    # Literal parameters are advertised as enums.
    hints = typing.get_type_hints(handler, include_extras=True)
    for name in params:
        inner = hints[name]
        if typing.get_origin(inner) is typing.Annotated:
            inner = typing.get_args(inner)[0]
        if typing.get_origin(inner) is typing.Literal:
            assert props[name].get("enum") == list(typing.get_args(inner)), f"{spec['name']}.{name}"


async def test_enum_violation_is_a_tool_error():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        result = await client.call_tool("denoise", {"channel": "x", "method": "gaussian"})
    assert result.is_error is True
    assert "percentile" in result.content[0].text


async def test_range_violation_is_a_tool_error():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        result = await client.call_tool("denoise", {"channel": "x", "filter_size": 0})
    assert result.is_error is True
    assert "filter_size" in result.content[0].text


async def test_missing_required_argument_is_a_tool_error():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        result = await client.call_tool("load_channel", {"cycle": "1"})
    assert result.is_error is True
    assert "channel" in result.content[0].text


async def test_structured_content_round_trip():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        result = await client.call_tool("get_processing_history", {"channel": "nope"})
    assert result.is_error is False
    assert result.structured_content["status"] == "success"
    assert result.structured_content["channel"] == "nope"
    assert result.content[0].type == "text"


# =============================================================================
# Errors, resources, prompts, session tools
# =============================================================================


async def test_resources_and_templates_are_advertised_and_readable(tmp_path, monkeypatch):
    import json

    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    monkeypatch.setenv("KINTSUGI_PROJECT_DIR", str(tmp_path))
    async with Client(create_server()) as client:
        uris = {r.uri for r in (await client.list_resources()).resources}
        templates = {
            t.uri_template for t in (await client.list_resource_templates()).resource_templates
        }
        images = await client.read_resource("kintsugi://session/images")
        history = await client.read_resource("kintsugi://session/history/CD3")
        channels = await client.read_resource("kintsugi://project/channels")

    assert uris == {
        "kintsugi://session/images",
        "kintsugi://learning/stats",
        "kintsugi://project/channels",
    }
    assert templates == {"kintsugi://session/history/{channel}"}
    state = json.loads(images.contents[0].text)
    assert state["project_root"] == str(tmp_path.resolve())
    assert json.loads(history.contents[0].text)["channel"] == "CD3"
    assert json.loads(channels.contents[0].text)["project"] == str(tmp_path.resolve())


async def test_prompt_is_advertised_and_renders():
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        prompts = {p.name: p for p in (await client.list_prompts()).prompts}
        rendered = await client.get_prompt(
            "signal_isolation_workflow",
            {"channel": "CD3", "blank_channel": "Blank", "tissue_type": "tonsil"},
        )

    prompt = prompts["signal_isolation_workflow"]
    assert prompt.title == "Guided signal isolation"
    assert {a.name for a in prompt.arguments} == {
        "channel",
        "blank_channel",
        "tissue_type",
        "project_path",
    }
    text = rendered.messages[0].content.text
    assert "subtract_blank" in text and "approve_and_learn" in text and "tonsil" in text


async def test_session_tools_round_trip(sample_numpy_image, clear_session_state):
    from mcp.client import Client

    from kintsugi.mcp.server import create_server
    from kintsugi.mcp.tools.signal_isolation import _store_image

    _store_image("CD3", sample_numpy_image)
    async with Client(create_server()) as client:
        state = await client.call_tool("get_session_state", {})
        unloaded = await client.call_tool("unload_channel", {"channel": "CD3"})
        after = await client.call_tool("get_session_state", {})

    assert state.is_error is False
    assert "CD3" in state.structured_content["loaded_images"]
    assert state.structured_content["memory_usage"]["total_mb"] > 0
    assert unloaded.structured_content["unloaded"] == ["CD3"]
    assert after.structured_content["total_loaded"] == 0


def test_session_state_has_a_typed_output_schema(surface):
    schema = surface["tools"]["get_session_state"].output_schema
    assert {"status", "project_root", "loaded_images", "memory_usage"} <= set(schema["properties"])


@pytest.fixture
def sample_numpy_image():
    import numpy as np

    return np.random.default_rng(0).integers(0, 4096, (64, 64), dtype=np.uint16)


@pytest.fixture
def clear_session_state():
    from kintsugi.mcp.tools.signal_isolation import clear_session

    clear_session()
    yield
    clear_session()


# =============================================================================
# In-band errors are MCP tool errors on the wire
# =============================================================================


async def test_in_band_error_is_reported_as_tool_error():
    """Handlers return {"error": ...}; on the wire that becomes isError with the JSON text."""
    import json

    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        result = await client.call_tool("unload_channel", {"channel": "never-loaded"})
    assert result.is_error is True
    text = result.content[0].text
    payload = json.loads(text[text.index("{") :])
    assert payload["error"] == "Channel not loaded: never-loaded"


async def test_direct_call_keeps_in_band_error_contract():
    """The same handler, called from Python, still returns the plain error dict."""
    from kintsugi.mcp.tools.workflow import unload_channel

    assert (await unload_channel(channel="never-loaded")) == {
        "error": "Channel not loaded: never-loaded"
    }
