"""Tests for the generated Claude Code MCP configuration and stdout hygiene."""

from __future__ import annotations

import json
from pathlib import Path

from click.testing import CliRunner

from kintsugi.project import PROJECT_CONFIG_FILE, KintsugiProject, write_claude_config


def _read(path: Path) -> dict:
    return json.loads(path.read_text())


def test_write_claude_config_creates_both_files(tmp_path):
    result = write_claude_config(tmp_path)

    mcp_json = _read(tmp_path / ".mcp.json")
    server = mcp_json["mcpServers"]["kintsugi"]
    assert server["type"] == "stdio"
    # The project travels in args: clients differ in the working directory and
    # environment they give a server, and several ignore a cwd field.
    assert server["args"] == ["mcp", "start", "--project-dir", str(tmp_path.resolve())]
    assert isinstance(server["command"], str) and server["command"]
    assert "cwd" not in server

    settings = _read(tmp_path / ".claude" / "settings.local.json")
    assert settings == {"enabledMcpjsonServers": ["kintsugi"]}
    # Approving every present and future project server is over-broad.
    assert "enableAllProjectMcpServers" not in settings

    assert result["mcp_json"] == tmp_path / ".mcp.json"
    assert result["server"] == server


def test_write_claude_config_merges_and_is_idempotent(tmp_path):
    (tmp_path / ".mcp.json").write_text(
        json.dumps({"mcpServers": {"other": {"command": "other-server"}}})
    )
    (tmp_path / ".claude").mkdir()
    (tmp_path / ".claude" / "settings.local.json").write_text(
        json.dumps({"permissions": {"allow": ["Bash(ls)"]}, "enabledMcpjsonServers": ["other"]})
    )

    write_claude_config(tmp_path)
    first = (tmp_path / ".mcp.json").read_text()
    write_claude_config(tmp_path)
    assert (tmp_path / ".mcp.json").read_text() == first

    mcp_json = _read(tmp_path / ".mcp.json")
    assert set(mcp_json["mcpServers"]) == {"other", "kintsugi"}
    assert mcp_json["mcpServers"]["other"] == {"command": "other-server"}

    settings = _read(tmp_path / ".claude" / "settings.local.json")
    assert settings["permissions"] == {"allow": ["Bash(ls)"]}
    assert settings["enabledMcpjsonServers"] == ["other", "kintsugi"]


def test_write_claude_config_force_rewrites_entry(tmp_path):
    (tmp_path / ".mcp.json").write_text(
        json.dumps({"mcpServers": {"kintsugi": {"command": "stale", "cwd": "/old"}}})
    )
    write_claude_config(tmp_path)
    assert _read(tmp_path / ".mcp.json")["mcpServers"]["kintsugi"]["command"] == "stale"

    write_claude_config(tmp_path, force=True)
    server = _read(tmp_path / ".mcp.json")["mcpServers"]["kintsugi"]
    assert server["command"] != "stale"
    assert "cwd" not in server


def test_mcp_config_cli_print_only_writes_nothing(tmp_path):
    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "config", str(tmp_path), "--print-only"])
    assert result.exit_code == 0, result.output
    assert '"kintsugi"' in result.output
    assert not (tmp_path / ".mcp.json").exists()


def test_mcp_config_cli_writes_files(tmp_path):
    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "config", str(tmp_path)])
    assert result.exit_code == 0, result.output
    assert (tmp_path / ".mcp.json").exists()
    assert (tmp_path / ".claude" / "settings.local.json").exists()
    assert "claude mcp add --scope project" in result.output


def test_project_load_writes_nothing_to_stdout(tmp_path, capsys):
    """KintsugiProject.load() runs inside the MCP server; stdout is the protocol stream."""
    (tmp_path / PROJECT_CONFIG_FILE).write_text(
        json.dumps({"config": {"name": "quiet"}, "kintsugi_path": ""})
    )
    project = KintsugiProject.load(tmp_path)
    assert project.config.name == "quiet"
    assert capsys.readouterr().out == ""


def test_write_claude_config_gitignores_the_local_settings(tmp_path):
    (tmp_path / ".gitignore").write_text("data/raw/\n")

    first = write_claude_config(tmp_path)
    write_claude_config(tmp_path)  # idempotent
    write_claude_config(tmp_path, force=True)

    text = (tmp_path / ".gitignore").read_text()
    assert first["gitignore"] == tmp_path / ".gitignore"
    assert text.startswith("data/raw/\n")  # existing rules untouched
    assert text.count(".claude/settings.local.json") == 1


def test_write_claude_config_respects_an_existing_equivalent_rule(tmp_path):
    (tmp_path / ".gitignore").write_text("/.claude/settings.local.json\n")

    write_claude_config(tmp_path)

    assert (tmp_path / ".gitignore").read_text() == "/.claude/settings.local.json\n"
