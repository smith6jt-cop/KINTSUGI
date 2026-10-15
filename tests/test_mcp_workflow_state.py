"""The MCP tools keep a persistent per-project workflow state (``<project>/.kintsugi/``)."""

from __future__ import annotations

import contextlib
import io
import json
from pathlib import Path

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")


@pytest.fixture
def project(tmp_path, monkeypatch):
    """A KINTSUGI project with cycle 1 CD3 + Blank, an isolated HOME and a clean session."""
    import tifffile

    from kintsugi.mcp.session import clear_session, set_project_root
    from kintsugi.mcp.tools import learning
    from kintsugi.project import KintsugiProject

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    with contextlib.redirect_stdout(io.StringIO()):
        proj = KintsugiProject.create(tmp_path / "proj", name="wf", setup_notebooks=False)
    raw = proj.paths.raw / "cyc001"
    raw.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    tifffile.imwrite(raw / "CD3.tif", rng.integers(500, 4000, (64, 64)).astype(np.uint16))
    tifffile.imwrite(raw / "Blank.tif", rng.integers(300, 900, (64, 64)).astype(np.uint16))
    set_project_root(proj.paths.root)
    clear_session()
    learning._learning_engines.clear()
    yield proj.paths.root
    clear_session()
    learning._learning_engines.clear()


def _restart():
    """What a client restarting the server does to in-process state."""
    from kintsugi.mcp.session import clear_session
    from kintsugi.mcp.workflow_store import clear_workflow_stores

    clear_session()
    clear_workflow_stores()


async def _process(project_root: Path) -> None:
    from kintsugi.mcp.tools import learning, quality_assessment, signal_isolation, workflow

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    await signal_isolation.load_channel(cycle="1", channel="Blank")
    sub = await signal_isolation.subtract_blank(
        signal_channel="cyc001_CD3", blank_channel="cyc001_Blank"
    )
    assert sub["status"] == "success", sub
    await quality_assessment.assess_quality(channel=sub["output_name"])
    saved = await workflow.save_processed(channel=sub["output_name"], output_name="CD3_iso")
    assert saved["status"] == "success", saved
    approved = await learning.approve_and_learn(
        channel=sub["output_name"],
        tissue_type="tonsil",
        marker_name="CD3",
        operations_params={"subtract_blank": {"blank_scale_factor": 0.8}},
        quality_before=0.3,
        quality_after=0.8,
        notes="clean",
    )
    assert approved["workflow_channel"] == "cyc001_CD3"


async def test_a_full_session_persists_and_survives_a_restart(project):
    from kintsugi.mcp.tools import workflow

    await _process(project)
    _restart()

    status = await workflow.get_workflow_status()
    assert status["status"] == "success"
    assert status["project_root"] == str(project.resolve())
    cd3 = status["channels"]["cyc001_CD3"]
    assert cd3["status"] == "approved" and cd3["quality_score"] == 0.8
    assert cd3["last_operation"] == "saved" and cd3["steps"] == 2  # subtract + saved
    assert status["channels"]["cyc001_Blank"]["status"] == "pending"
    assert status["approved"] == ["cyc001_CD3"]

    # The derived image found its project through lineage, and the steps are
    # available after the restart from the persisted state.
    history = await workflow.get_processing_history(channel="cyc001_CD3")
    assert history["source"] == "persisted"
    assert [s["operation"] for s in history["steps"]] == ["subtract_blank", "saved"]
    assert history["steps"][0]["image"].startswith("cyc001_CD3")


async def test_review_decision_on_an_unloaded_channel_after_a_restart(project):
    from kintsugi.mcp.tools import signal_isolation, workflow

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    _restart()

    result = await workflow.set_review_decision(
        channel="cyc001_CD3", decision="needs_adjustment", notes="halo around follicles"
    )
    assert result["status"] == "success"
    status = await workflow.get_workflow_status()
    assert status["channels"]["cyc001_CD3"]["status"] == "needs_adjustment"
    assert status["processing_queue"] == ["cyc001_CD3"]

    missing = await workflow.set_review_decision(channel="cyc009_X", decision="rejected")
    assert "not in the workflow state" in missing["error"]


async def test_read_only_tools_do_not_write_state(project):
    from kintsugi.mcp.tools import quality_assessment, signal_isolation, visualization, workflow

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    state_file = project / ".kintsugi" / "workflow_state.json"
    before = state_file.read_text()

    await quality_assessment.assess_quality(channel="cyc001_CD3")
    await quality_assessment.compute_snr(channel="cyc001_CD3")
    await visualization.get_image_stats(channel="cyc001_CD3")
    await workflow.get_workflow_status()
    await workflow.get_session_state()
    assert state_file.read_text() == before

    from kintsugi.mcp.session import _loaded_images, _processing_history

    # assess_quality annotates the session metadata instead of the history.
    assert "quality_score" in _loaded_images["cyc001_CD3"]["metadata"]
    assert _processing_history["cyc001_CD3"] == []


async def test_unloading_keeps_persisted_state(project):
    from kintsugi.mcp.tools import signal_isolation, workflow

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    await workflow.unload_channel(clear_all=True)
    status = await workflow.get_workflow_status()
    assert list(status["channels"]) == ["cyc001_CD3"]


async def test_session_state_includes_the_workflow_summary(project):
    from kintsugi.mcp.tools import signal_isolation, workflow

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    state = await workflow.get_session_state()
    assert state["workflow"]["total_channels"] == 1
    assert state["workflow"]["status_breakdown"] == {"pending": 1}


async def test_without_a_project_nothing_is_tracked(tmp_path, monkeypatch):
    from kintsugi.mcp.session import _store_image, clear_session
    from kintsugi.mcp.tools import workflow

    monkeypatch.chdir(tmp_path)
    clear_session()
    _store_image("ch", np.ones((4, 4), dtype=np.uint16))
    assert "No project configured" in (await workflow.get_workflow_status())["error"]
    assert (
        "No project configured"
        in (await workflow.set_review_decision(channel="ch", decision="rejected"))["error"]
    )
    assert "workflow" not in await workflow.get_session_state()
    not_a_project = tmp_path / "plain"
    not_a_project.mkdir()
    result = await workflow.get_workflow_status(project_path=str(not_a_project))
    assert "not a KINTSUGI project" in result["error"]
    assert list(tmp_path.rglob(".kintsugi")) == []
    clear_session()


async def test_lineage_is_inherited_by_derived_images(project):
    from kintsugi.mcp.session import _loaded_images
    from kintsugi.mcp.tools import signal_isolation

    await signal_isolation.load_channel(cycle="1", channel="CD3")
    await signal_isolation.load_channel(cycle="1", channel="Blank")
    sub = await signal_isolation.subtract_blank(
        signal_channel="cyc001_CD3", blank_channel="cyc001_Blank"
    )
    den = await signal_isolation.denoise(channel=sub["output_name"])
    meta = _loaded_images[den["output_name"]]["metadata"]
    assert meta["root"] == "cyc001_CD3" and meta["cycle"] == "cyc001"
    assert meta["project"] == str(project)


async def test_resource_and_typed_output(project):
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        await client.call_tool("load_channel", {"cycle": "1", "channel": "CD3"})
        resource = await client.read_resource("kintsugi://project/workflow")
        call = await client.call_tool("get_workflow_status", {})
        tools = {t.name: t for t in (await client.list_tools()).tools}

    data = json.loads(resource.contents[0].text)
    assert data["channels"]["cyc001_CD3"]["status"] == "pending"
    assert call.structured_content["channels"] == data["channels"]
    props = tools["get_workflow_status"].output_schema["properties"]
    assert {"channels", "review_queue", "processing_queue", "project_root"} <= set(props)
    assert tools["set_review_decision"].annotations.read_only_hint is False


async def test_default_subtract_blank_survives_serialisation(project):
    """numpy values in a tool result (here numpy.bool_ parameters) used to fail the call."""
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        await client.call_tool("load_channel", {"cycle": "1", "channel": "CD3"})
        await client.call_tool("load_channel", {"cycle": "1", "channel": "Blank"})
        result = await client.call_tool(
            "subtract_blank", {"signal_channel": "cyc001_CD3", "blank_channel": "cyc001_Blank"}
        )
    assert not result.is_error, result.content[0].text
    assert result.structured_content["status"] == "success"
