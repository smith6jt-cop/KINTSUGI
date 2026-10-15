"""Tests for the generated MCP client configuration (``write_mcp_config``) and stdout hygiene."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from click.testing import CliRunner

from kintsugi.project import PROJECT_CONFIG_FILE, KintsugiProject, write_mcp_config


def _read(path: Path) -> dict:
    return json.loads(path.read_text())


@pytest.fixture
def project(tmp_path) -> Path:
    root = tmp_path / "proj"
    root.mkdir()
    (root / PROJECT_CONFIG_FILE).write_text(json.dumps({"config": {"name": "Tonsil run"}}))
    return root


def test_claude_code_files(project):
    result = write_mcp_config(project, ["claude-code"])

    server = _read(project / ".mcp.json")["mcpServers"]["kintsugi"]
    assert server["type"] == "stdio"
    # The project travels in args: clients differ in the working directory and
    # environment they give a server, and several ignore a cwd field.
    assert server["args"] == ["mcp", "start", "--project-dir", str(project.resolve())]
    assert isinstance(server["command"], str) and server["command"]
    assert "cwd" not in server

    settings = _read(project / ".claude" / "settings.local.json")
    assert settings == {"enabledMcpjsonServers": ["kintsugi"]}
    # Approving every present and future project server is over-broad.
    assert "enableAllProjectMcpServers" not in settings

    assert project / ".mcp.json" in result["written"]
    assert result["server"]["args"] == server["args"]


def test_merges_and_is_idempotent(project):
    (project / ".mcp.json").write_text(
        json.dumps({"mcpServers": {"other": {"command": "other-server"}}})
    )
    (project / ".claude").mkdir()
    (project / ".claude" / "settings.local.json").write_text(
        json.dumps({"permissions": {"allow": ["Bash(ls)"]}, "enabledMcpjsonServers": ["other"]})
    )

    write_mcp_config(project, ["claude-code"])
    first = (project / ".mcp.json").read_text()
    again = write_mcp_config(project, ["claude-code"])
    assert (project / ".mcp.json").read_text() == first
    assert again["written"] == [] and project / ".mcp.json" in again["unchanged"]

    mcp_json = _read(project / ".mcp.json")
    assert set(mcp_json["mcpServers"]) == {"other", "kintsugi"}
    assert mcp_json["mcpServers"]["other"] == {"command": "other-server"}

    settings = _read(project / ".claude" / "settings.local.json")
    assert settings["permissions"] == {"allow": ["Bash(ls)"]}
    assert settings["enabledMcpjsonServers"] == ["other", "kintsugi"]


def test_force_rewrites_only_our_entry(project):
    (project / ".mcp.json").write_text(
        json.dumps(
            {"mcpServers": {"kintsugi": {"command": "stale", "cwd": "/old"}, "other": {"x": 1}}}
        )
    )
    write_mcp_config(project, ["claude-code"])
    assert _read(project / ".mcp.json")["mcpServers"]["kintsugi"]["command"] == "stale"

    write_mcp_config(project, ["claude-code"], force=True)
    servers = _read(project / ".mcp.json")["mcpServers"]
    assert servers["kintsugi"]["command"] != "stale" and "cwd" not in servers["kintsugi"]
    assert servers["other"] == {"x": 1}


def test_refuses_a_directory_that_is_not_a_project(tmp_path):
    with pytest.raises(ValueError, match="not a KINTSUGI project"):
        write_mcp_config(tmp_path)
    assert list(tmp_path.iterdir()) == []
    write_mcp_config(tmp_path, ["cursor"], require_project=False)
    assert (tmp_path / ".cursor" / "mcp.json").exists()


def test_cli_print_only_writes_nothing_and_shows_project_dir(project):
    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "config", str(project), "--print-only"])
    assert result.exit_code == 0, result.output
    flat = " ".join(result.output.split())
    assert '"kintsugi"' in flat and "--project-dir" in flat
    assert "[mcp_servers.kintsugi]" in flat  # Codex (TOML)
    assert sorted(p.name for p in project.iterdir()) == [PROJECT_CONFIG_FILE]


def test_cli_writes_every_project_client_and_prints_snippets(project):
    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "config", str(project)])
    assert result.exit_code == 0, result.output
    for rel in (
        ".mcp.json",
        ".claude/settings.local.json",
        ".vscode/mcp.json",
        ".cursor/mcp.json",
        ".gemini/settings.json",
        ".agents/mcp_config.json",
        ".codex/config.toml",
        ".zed/settings.json",
        "AGENTS.md",
    ):
        assert (project / rel).exists(), rel
    flat = " ".join(result.output.split())
    for display in ("Claude Desktop", "Windsurf", "Cline", "JetBrains"):
        assert display in flat  # snippet panels for user / IDE scoped clients


def test_cli_client_option_restricts_and_rejects_unknown_names(project):
    from kintsugi.cli import main

    runner = CliRunner()
    result = runner.invoke(main, ["mcp", "config", str(project), "--client", "cursor"])
    assert result.exit_code == 0, result.output
    assert (project / ".cursor" / "mcp.json").exists()
    assert not (project / ".mcp.json").exists()

    bad = runner.invoke(main, ["mcp", "config", str(project), "--client", "notaclient"])
    assert bad.exit_code == 2 and "unknown MCP client" in " ".join(bad.output.split())


def test_cli_refuses_a_non_project_without_force(tmp_path):
    from kintsugi.cli import main

    runner = CliRunner()
    result = runner.invoke(main, ["mcp", "config", str(tmp_path)])
    assert result.exit_code == 2 and "not a KINTSUGI project" in " ".join(result.output.split())
    assert list(tmp_path.iterdir()) == []
    assert runner.invoke(main, ["mcp", "config", str(tmp_path), "--force"]).exit_code == 0
    assert (tmp_path / ".mcp.json").exists()


def test_project_load_writes_nothing_to_stdout(tmp_path, capsys):
    """KintsugiProject.load() runs inside the MCP server; stdout is the protocol stream."""
    (tmp_path / PROJECT_CONFIG_FILE).write_text(
        json.dumps({"config": {"name": "quiet"}, "kintsugi_path": ""})
    )
    project = KintsugiProject.load(tmp_path)
    assert project.config.name == "quiet"
    assert capsys.readouterr().out == ""


def test_gitignores_generated_files_once_under_one_header(project):
    (project / ".gitignore").write_text("data/raw/\n")

    write_mcp_config(project)
    write_mcp_config(project)  # idempotent
    write_mcp_config(project, force=True)

    text = (project / ".gitignore").read_text()
    assert text.startswith("data/raw/\n")  # existing rules untouched
    for rel in (".claude/settings.local.json", ".mcp.json", ".codex/config.toml"):
        assert text.splitlines().count(rel) == 1, rel
    assert text.count("MCP client configs") == 1


def test_respects_an_existing_equivalent_rule(project):
    (project / ".gitignore").write_text("/.claude/settings.local.json\n/.mcp.json\n")

    write_mcp_config(project, ["claude-code"])

    assert (project / ".gitignore").read_text() == "/.claude/settings.local.json\n/.mcp.json\n"
