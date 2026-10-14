"""Progress notifications, worker-thread execution and the default project root.

These tests drive the server through the SDK's in-memory client so that the
``ctx`` injection, ``ctx.report_progress`` and ``asyncio.to_thread`` paths in
``kintsugi.mcp.tools.signal_isolation`` are exercised end to end, exactly as
Claude Code would exercise them over the protocol.

Run with: pytest tests/test_mcp_progress.py -v
"""

from __future__ import annotations

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE

pytestmark = pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")

# Handlers that take an injected ``ctx`` and report progress.
PROGRESS_TOOLS = {
    "subtract_blank",
    "denoise_advanced",
    "optimize_parameters",
    "cluster_channels",
    "propagate_parameters",
}


# =============================================================================
# Fixtures and helpers
# =============================================================================


@pytest.fixture
def seeded_images():
    """Seed a small signal/blank pair into the shared session and clear it afterwards."""
    from kintsugi.mcp.tools.signal_isolation import (
        _loaded_images,
        _processing_history,
        _store_image,
    )

    _loaded_images.clear()
    _processing_history.clear()

    rng = np.random.default_rng(7)
    signal = rng.integers(100, 2000, (64, 64), dtype=np.uint16)
    signal[16:32, 16:32] = rng.integers(3000, 8000, (16, 16), dtype=np.uint16)
    blank = rng.integers(50, 400, (64, 64), dtype=np.uint16)
    _store_image("sig", signal)
    _store_image("blk", blank)

    yield

    _loaded_images.clear()
    _processing_history.clear()


async def _call_with_progress(name: str, arguments: dict):
    """Call ``name`` over the in-memory client, collecting progress notifications."""
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    events: list[tuple[float, float | None, str | None]] = []

    async def on_progress(progress: float, total: float | None, message: str | None) -> None:
        events.append((progress, total, message))

    async with Client(create_server()) as client:
        result = await client.call_tool(name, arguments, progress_callback=on_progress)

    return result, events


def _assert_progress_sequence(events):
    """Progress arrived, carries a total, and finishes at total/total."""
    assert events, "no progress notification arrived"
    assert all(total is not None for _, total, _ in events)
    last_progress, last_total, _ = events[-1]
    assert last_progress == last_total


# =============================================================================
# Schema: ctx is injected, never exposed
# =============================================================================


async def test_ctx_is_not_advertised_in_input_schemas():
    """The injected context parameter must not leak into the tool's inputSchema."""
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}

    for name in PROGRESS_TOOLS:
        assert "ctx" not in tools[name].input_schema["properties"], name
        assert "ctx" not in (tools[name].input_schema.get("required") or []), name


async def test_load_channel_project_path_is_optional():
    """project_path defaults to the Claude Code project, so it is no longer required."""
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        tools = {t.name: t for t in (await client.list_tools()).tools}

    schema = tools["load_channel"].input_schema
    assert "project_path" in schema["properties"]
    assert "project_path" not in schema["required"]
    assert set(schema["required"]) == {"cycle", "channel"}


# =============================================================================
# subtract_blank
# =============================================================================


async def test_subtract_blank_global_reports_progress(seeded_images):
    pytest.importorskip("dask_image")

    result, events = await _call_with_progress(
        "subtract_blank",
        {
            "signal_channel": "sig",
            "blank_channel": "blk",
            "blank_clip_factor": 0,
            "blank_scale_factor": 1.0,
            "auto_suggest": False,
        },
    )

    assert result.is_error is False
    assert result.structured_content["status"] == "success"
    assert result.structured_content["method"] == "global"
    _assert_progress_sequence(events)
    # Start (0/total) is always the first notification.
    assert events[0][0] == 0


async def test_subtract_blank_weighted_reports_progress(seeded_images):
    result, events = await _call_with_progress(
        "subtract_blank",
        {
            "signal_channel": "sig",
            "blank_channel": "blk",
            "method": "weighted",
            "blank_clip_factor": 0,
            "blank_scale_factor": 1.0,
        },
    )

    assert result.is_error is False
    assert result.structured_content["status"] == "success"
    assert result.structured_content["method"] == "weighted"
    _assert_progress_sequence(events)
    assert events[0][0] == 0


async def test_subtract_blank_direct_call_needs_no_ctx(seeded_images):
    """Direct (in-process) calls pass no ctx and behave exactly as before."""
    pytest.importorskip("dask_image")
    from kintsugi.mcp.tools.signal_isolation import subtract_blank

    result = await subtract_blank(
        signal_channel="sig",
        blank_channel="blk",
        blank_clip_factor=0,
        blank_scale_factor=1.0,
        auto_suggest=False,
    )
    assert result["status"] == "success"
    assert result["output_name"] == "sig_subtracted"


# =============================================================================
# denoise_advanced
# =============================================================================


async def test_denoise_advanced_reports_progress(seeded_images):
    result, events = await _call_with_progress(
        "denoise_advanced", {"channel": "sig", "method": "adaptive", "strength": "light"}
    )

    assert result.is_error is False
    assert result.structured_content["status"] == "success"
    assert result.structured_content["output_name"] == "sig_adaptive_denoised"
    _assert_progress_sequence(events)
    assert any("denoising" in (message or "") for _, _, message in events)


# =============================================================================
# optimize_parameters
# =============================================================================


async def test_optimize_parameters_reports_progress(seeded_images):
    pytest.importorskip("optuna")

    result, events = await _call_with_progress(
        "optimize_parameters",
        {
            "signal_channel": "sig",
            "blank_channel": "blk",
            "n_trials": 2,
            "timeout": 30,
            "optimize_clean": False,
        },
    )

    assert result.is_error is False, result.content
    assert result.structured_content["status"] == "success"
    assert result.structured_content["n_trials_completed"] >= 1
    _assert_progress_sequence(events)
    # The total is the trial budget and the run starts at trial 0.
    assert events[0] == (0, 2, "trial 0/2")
    # Per-trial notifications are scheduled from Optuna's worker thread back
    # onto the event loop; at least one must have been delivered.
    assert any(
        progress >= 1 and (message or "").startswith("trial ") for progress, _, message in events
    ), events


async def test_optimizer_trial_callback_counts_trials():
    """optimize_signal_isolation forwards a (completed, n_trials) callback to Optuna."""
    pytest.importorskip("optuna")
    from kintsugi.signal.optimizer import optimize_signal_isolation

    rng = np.random.default_rng(3)
    signal = rng.integers(100, 2000, (32, 32), dtype=np.uint16)
    blank = rng.integers(50, 400, (32, 32), dtype=np.uint16)

    seen: list[tuple[int, int]] = []
    result = optimize_signal_isolation(
        signal,
        blank,
        n_trials=3,
        timeout=30,
        optimize_clean=False,
        show_progress=False,
        trial_callback=lambda completed, total: seen.append((completed, total)),
    )

    assert result["status"] == "success"
    assert seen == [(1, 3), (2, 3), (3, 3)]


# =============================================================================
# cluster_channels / propagate_parameters
# =============================================================================


async def test_cluster_and_propagate_report_progress(tmp_path, seeded_images):
    pytest.importorskip("sklearn")
    tifffile = pytest.importorskip("tifffile")

    project = tmp_path / "proj"
    edf = project / "data" / "processed" / "edf"
    edf.mkdir(parents=True)
    rng = np.random.default_rng(1)
    for name, level in (("CD3", 200), ("CD4", 220), ("DAPI", 3000), ("CD8", 210)):
        tifffile.imwrite(edf / f"{name}.tif", rng.integers(0, level, (64, 64), dtype=np.uint16))

    clustered, cluster_events = await _call_with_progress(
        "cluster_channels", {"project_path": str(project), "n_clusters": 2}
    )
    assert clustered.is_error is False, clustered.content
    assert clustered.structured_content["status"] == "success"
    _assert_progress_sequence(cluster_events)

    cid = clustered.structured_content["assignments"]["CD3"]
    propagated, propagate_events = await _call_with_progress(
        "propagate_parameters",
        {
            "project_path": str(project),
            "cluster_id": cid,
            "params": {"clean_params": {"background_threshold": 10}},
        },
    )
    assert propagated.is_error is False, propagated.content
    assert propagated.structured_content["status"] == "success"
    _assert_progress_sequence(propagate_events)


# =============================================================================
# Default project root
# =============================================================================


async def test_load_channel_defaults_to_claude_project_dir(tmp_path, monkeypatch, seeded_images):
    """Without project_path, load_channel reads CLAUDE_PROJECT_DIR."""
    pytest.importorskip("zarr")
    tifffile = pytest.importorskip("tifffile")
    from kintsugi.mcp.tools.signal_isolation import _loaded_images, load_channel

    project = tmp_path / "proj"
    raw = project / "data" / "raw" / "cyc001"
    raw.mkdir(parents=True)
    rng = np.random.default_rng(5)
    tifffile.imwrite(raw / "CD3_cyc001.tif", rng.integers(0, 500, (64, 64), dtype=np.uint16))

    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(project))
    result = await load_channel(cycle="1", channel="CD3")

    assert result.get("status") == "loaded", result
    assert result["name"] == "cyc001_CD3"
    assert "cyc001_CD3" in _loaded_images
    assert _loaded_images["cyc001_CD3"]["metadata"]["project"] == str(project.resolve())
