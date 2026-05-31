"""
Tests for MCP server P0 hardening: path-safety helpers, tool-level input
validation, and config portability.

Run with: pytest tests/test_mcp_path_safety.py -v
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from kintsugi.mcp.path_safety import (
    PathSafetyError,
    clamp_float,
    clamp_int,
    ensure_within,
    safe_filename,
    validate_choice,
)

# =============================================================================
# safe_filename
# =============================================================================


@pytest.mark.parametrize(
    "name",
    ["CD3_20260101", "DAPI", "foo.bar", "channel_subtracted_v2"],
)
def test_safe_filename_accepts_bare_names(name):
    assert safe_filename(name) == name


@pytest.mark.parametrize(
    "name",
    ["../evil", "a/b", "/etc/passwd", "..", ".", "a\\b", "x\x00y", "", "   "],
)
def test_safe_filename_rejects_traversal_and_separators(name):
    with pytest.raises(PathSafetyError):
        safe_filename(name)


def test_safe_filename_uses_param_in_message():
    with pytest.raises(PathSafetyError, match="output_name"):
        safe_filename("../x", "output_name")


# =============================================================================
# ensure_within
# =============================================================================


def test_ensure_within_accepts_child(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    child = root / "data" / "out.tif"
    assert ensure_within(child, root) == child.resolve()


def test_ensure_within_accepts_root_itself(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    assert ensure_within(root, root) == root.resolve()


def test_ensure_within_rejects_escape(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    with pytest.raises(PathSafetyError):
        ensure_within(root / ".." / "secret.txt", root)


def test_ensure_within_rejects_symlink_escape(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    secret = outside / "secret.txt"
    secret.write_text("classified")

    link = root / "link"
    try:
        link.symlink_to(secret)
    except (OSError, NotImplementedError):
        pytest.skip("symlinks not supported on this platform/filesystem")

    # resolve() follows the symlink to the outside target -> must be rejected.
    with pytest.raises(PathSafetyError):
        ensure_within(link, root)


# =============================================================================
# validate_choice / clamp_int / clamp_float
# =============================================================================


def test_validate_choice_accepts_member():
    assert validate_choice("global", {"global", "weighted"}, "method") == "global"


def test_validate_choice_rejects_nonmember():
    with pytest.raises(ValueError, match="method"):
        validate_choice("bogus", {"global", "weighted"}, "method")


@pytest.mark.parametrize("value", [1, 5, 99])
def test_clamp_int_accepts_in_range(value):
    assert clamp_int(value, "filter_size", minimum=1, maximum=99) == value


@pytest.mark.parametrize("value", [0, 100, -3])
def test_clamp_int_rejects_out_of_range(value):
    with pytest.raises(ValueError):
        clamp_int(value, "filter_size", minimum=1, maximum=99)


def test_clamp_int_rejects_bool_and_nonint():
    with pytest.raises(ValueError):
        clamp_int(True, "n", minimum=0, maximum=10)
    with pytest.raises(ValueError):
        clamp_int(3.5, "n", minimum=0, maximum=10)


def test_clamp_float_bounds():
    assert clamp_float(0.5, "clip_limit", minimum=0.0, maximum=1.0) == 0.5
    with pytest.raises(ValueError):
        clamp_float(1.5, "clip_limit", minimum=0.0, maximum=1.0)


# =============================================================================
# Tool-level validation (handlers return in-band {"error": ...})
# =============================================================================


@pytest.fixture
def loaded_channel():
    """Load a small image into the shared MCP session state and clean up after."""
    from kintsugi.mcp.tools import signal_isolation as si

    name = "cyc001_TEST"
    img = np.random.randint(100, 500, (64, 64), dtype=np.uint16)
    si._store_image(name, img, {"project": None})
    yield name
    si.clear_session()


async def test_save_processed_rejects_traversal_output_name(tmp_path, loaded_channel):
    from kintsugi.mcp.tools import workflow

    project = tmp_path / "proj"
    (project / "data" / "processed").mkdir(parents=True)

    result = await workflow.save_processed(
        channel=loaded_channel,
        output_name="../../evil",
        project_path=str(project),
    )

    assert "error" in result
    # Nothing should have been written outside the project tree.
    assert not (tmp_path / "evil.tiff").exists()


async def test_save_processed_rejects_bad_format(tmp_path, loaded_channel):
    from kintsugi.mcp.tools import workflow

    project = tmp_path / "proj"
    (project / "data" / "processed").mkdir(parents=True)

    result = await workflow.save_processed(
        channel=loaded_channel,
        output_name="ok_name",
        format="exe",
        project_path=str(project),
    )
    assert "error" in result
    assert "format" in result["error"]


async def test_denoise_rejects_out_of_range_filter_size(loaded_channel):
    from kintsugi.mcp.tools import signal_isolation as si

    result = await si.denoise(channel=loaded_channel, method="median", filter_size=0)
    assert "error" in result
    assert "filter_size" in result["error"]


async def test_denoise_rejects_unknown_method(loaded_channel):
    from kintsugi.mcp.tools import signal_isolation as si

    result = await si.denoise(channel=loaded_channel, method="wavelet")
    assert "error" in result
    assert "method" in result["error"]


async def test_optimize_parameters_rejects_huge_n_trials(loaded_channel):
    from kintsugi.mcp.tools import signal_isolation as si

    result = await si.optimize_parameters(
        signal_channel=loaded_channel,
        blank_channel=loaded_channel,
        n_trials=10_000_000,
    )
    assert "error" in result
    assert "n_trials" in result["error"]


async def test_denoise_advanced_rejects_missing_model_path(loaded_channel):
    from kintsugi.mcp.tools import signal_isolation as si

    result = await si.denoise_advanced(
        channel=loaded_channel,
        method="n2v",
        model_path="/nonexistent/path/model.pt",
    )
    assert "error" in result
    assert "model_path" in result["error"]


# =============================================================================
# Config portability
# =============================================================================


def test_resolve_executable_returns_nonempty():
    from kintsugi.project import _resolve_kintsugi_executable

    cmd = _resolve_kintsugi_executable()
    assert isinstance(cmd, str) and cmd


def test_committed_mcp_json_is_portable():
    """The repo's own .mcp.json must not hardcode a machine-specific path."""
    repo_root = Path(__file__).resolve().parents[1]
    mcp_file = repo_root / ".mcp.json"
    if not mcp_file.exists():
        pytest.skip(".mcp.json not present")

    config = json.loads(mcp_file.read_text())
    server = config["mcpServers"]["kintsugi"]

    assert server["command"] == "kintsugi", "command should be a PATH-resolvable name"
    assert not Path(server["command"]).is_absolute()
    # No hardcoded absolute cwd that ties the repo to one machine/account.
    assert "cwd" not in server or not Path(server["cwd"]).is_absolute()
