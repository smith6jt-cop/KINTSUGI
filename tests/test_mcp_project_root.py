"""How the MCP server decides which project ``project_path`` defaults to.

Clients differ in the working directory and environment they give a server, so
the project is configured explicitly (``kintsugi mcp start --project-dir``, then
``KINTSUGI_PROJECT_DIR``) and the working directory counts only when it is a
KINTSUGI project. Without a project, tools say so and write nothing.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from kintsugi.mcp._sdk import MCP_AVAILABLE
from kintsugi.mcp.session import (
    PROJECT_MARKER,
    ProjectNotConfigured,
    configured_project_root,
    learning_project_root,
    resolve_project_path,
    set_project_root,
)


def _project(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True)
    (path / PROJECT_MARKER).write_text(json.dumps({"config": {"name": path.name}}))
    return path


@pytest.fixture
def bare_cwd(tmp_path, monkeypatch):
    """A working directory that is not a KINTSUGI project."""
    cwd = tmp_path / "cwd"
    cwd.mkdir()
    monkeypatch.chdir(cwd)
    return cwd


def test_precedence_explicit_then_env_then_project_cwd(tmp_path, monkeypatch, bare_cwd):
    assert configured_project_root() is None

    cwd_project = _project(tmp_path / "cwd_project")
    monkeypatch.chdir(cwd_project)
    assert configured_project_root() == cwd_project.resolve()

    env_dir = tmp_path / "env"  # taken as given, like any explicit setting
    monkeypatch.setenv("KINTSUGI_PROJECT_DIR", str(env_dir))
    assert configured_project_root() == env_dir.resolve()

    explicit = _project(tmp_path / "explicit")
    set_project_root(explicit)
    assert configured_project_root() == explicit.resolve()

    set_project_root(None)
    monkeypatch.setenv("KINTSUGI_PROJECT_DIR", "  ")  # blank counts as unset
    assert configured_project_root() == cwd_project.resolve()


def test_claude_project_dir_is_not_read(tmp_path, monkeypatch, bare_cwd):
    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(_project(tmp_path / "p")))
    assert configured_project_root() is None
    with pytest.raises(ProjectNotConfigured, match="--project-dir"):
        resolve_project_path(None)
    assert resolve_project_path("/explicit") == "/explicit"


def test_learning_default_needs_a_real_project(tmp_path, monkeypatch, bare_cwd):
    monkeypatch.setenv("KINTSUGI_PROJECT_DIR", str(tmp_path / "not_a_project"))
    assert learning_project_root() is None
    project = _project(tmp_path / "real")
    monkeypatch.setenv("KINTSUGI_PROJECT_DIR", str(project))
    assert learning_project_root() == str(project.resolve())


def test_learning_engine_defaults_to_the_configured_project(tmp_path, monkeypatch, bare_cwd):
    from kintsugi.mcp.tools import learning

    monkeypatch.setenv("HOME", str(tmp_path / "home"))
    monkeypatch.setenv("USERPROFILE", str(tmp_path / "home"))
    learning._learning_engines.clear()
    try:
        assert learning._get_learning_engine().project_path is None  # user-wide only
        project = _project(tmp_path / "proj")
        set_project_root(project)
        engine = learning._get_learning_engine()
        assert engine.project_path is not None
        assert Path(engine.project_path).resolve() == project.resolve()
    finally:
        learning._learning_engines.clear()


@pytest.mark.parametrize("tool", ["load_channel", "list_channels"])
async def test_filesystem_tools_refuse_without_a_project(tool, bare_cwd):
    from kintsugi.mcp.tools import signal_isolation, workflow

    if tool == "load_channel":
        result = await signal_isolation.load_channel(cycle="1", channel="CD3")
    else:
        result = await workflow.list_channels()
    assert "No project configured" in result["error"]
    assert list(bare_cwd.iterdir()) == []


async def test_save_processed_never_writes_under_cwd(tmp_path, bare_cwd):
    from kintsugi.mcp.session import _store_image, clear_session
    from kintsugi.mcp.tools.workflow import save_processed

    clear_session()
    try:
        _store_image("ch", np.ones((8, 8), dtype=np.uint16))
        result = await save_processed(channel="ch", output_name="out")
        assert "No project configured" in result["error"]
        assert list(bare_cwd.iterdir()) == []

        project = _project(tmp_path / "proj")
        set_project_root(project)
        saved = await save_processed(channel="ch", output_name="out")
        assert saved["status"] == "success", saved
        assert Path(saved["output_path"]).is_relative_to(project.resolve())
    finally:
        clear_session()


async def test_session_state_reports_null_root_without_a_project(bare_cwd):
    from kintsugi.mcp.tools.workflow import get_session_state

    assert (await get_session_state())["project_root"] is None


def test_model_roots_ignore_claude_project_dir(tmp_path, monkeypatch, bare_cwd):
    from kintsugi.mcp.path_safety import allowed_model_roots

    monkeypatch.setenv("CLAUDE_PROJECT_DIR", str(tmp_path / "claude"))
    monkeypatch.setenv("KINTSUGI_MODEL_DIR", str(tmp_path / "models"))
    assert allowed_model_roots(None) == [tmp_path / "models"]


def test_mcp_start_rejects_a_directory_that_is_not_a_project(tmp_path):
    from click.testing import CliRunner

    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "start", "--project-dir", str(tmp_path)])
    assert result.exit_code == 2
    assert "is not a KINTSUGI project" in " ".join(result.output.split())  # Rich wraps lines
    assert configured_project_root() != tmp_path.resolve()


@pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")
async def test_resources_without_a_project(bare_cwd):
    from mcp.client import Client

    from kintsugi.mcp.server import create_server

    async with Client(create_server()) as client:
        images = await client.read_resource("kintsugi://session/images")
        channels = await client.read_resource("kintsugi://project/channels")
        call = await client.call_tool("list_channels", {})

    assert json.loads(images.contents[0].text)["project_root"] is None
    assert "No project configured" in json.loads(channels.contents[0].text)["error"]
    assert call.is_error and "No project configured" in call.content[0].text


@pytest.mark.skipif(not MCP_AVAILABLE, reason="mcp 2.x not installed")
async def test_learning_reads_create_no_files_and_approval_creates_state(
    tmp_path, monkeypatch, bare_cwd
):
    from mcp.client import Client

    from kintsugi.mcp.server import create_server
    from kintsugi.mcp.session import _store_image, clear_session
    from kintsugi.mcp.tools import learning

    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    project = _project(tmp_path / "proj")
    set_project_root(project)
    learning._learning_engines.clear()
    clear_session()
    try:
        _store_image("cyc001_CD3", np.random.default_rng(0).integers(0, 4000, (32, 32)))
        lookup = {"tissue_type": "tonsil", "marker_name": "CD3"}
        async with Client(create_server()) as client:
            reads = [
                await client.call_tool(
                    "get_learned_parameters", {**lookup, "operation": "denoise"}
                ),
                await client.call_tool("get_learning_statistics", {}),
                await client.call_tool(
                    "suggest_with_learning", {"channel": "cyc001_CD3", **lookup}
                ),
            ]
            await client.read_resource("kintsugi://learning/stats")
            assert not any(r.is_error for r in reads), [r.content[0].text for r in reads]
            assert not (project / ".kintsugi").exists()
            assert list(home.iterdir()) == []

            approved = await client.call_tool(
                "approve_and_learn",
                {
                    "channel": "cyc001_CD3",
                    **lookup,
                    "operations_params": {"denoise": {"filter_size": 3}},
                    "quality_before": 0.3,
                    "quality_after": 0.8,
                },
            )
        assert not approved.is_error, approved.content[0].text
        assert (project / ".kintsugi" / "parameter_learning.db").exists()
        assert (home / ".kintsugi" / "parameter_learning.db").exists()
    finally:
        learning._learning_engines.clear()
        clear_session()
