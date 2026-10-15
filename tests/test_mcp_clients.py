"""Tests for the MCP client registry (``kintsugi.mcp.clients``) and the generated AGENTS.md."""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from kintsugi.mcp import clients as mc
from kintsugi.mcp.docs_sync import render_tool_list
from kintsugi.project import PROJECT_CONFIG_FILE, write_mcp_config

try:  # Python 3.11+
    import tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10
    tomllib = pytest.importorskip("tomli")

COMMAND = "/opt/env/bin/kintsugi"


def _parse(client: mc.ClientSpec, text: str) -> dict:
    return tomllib.loads(text) if client.format == "toml" else json.loads(text)


@pytest.fixture
def project(tmp_path) -> Path:
    root = tmp_path / "proj"
    root.mkdir()
    (root / PROJECT_CONFIG_FILE).write_text(json.dumps({"config": {"name": "Spleen"}}))
    return root


@pytest.mark.parametrize("client", mc.project_clients(), ids=lambda c: c.name)
def test_every_project_client_renders_a_parseable_entry(client, project):
    data = _parse(client, mc.render(client, project, COMMAND))
    entry = data[client.key]["kintsugi"]
    assert entry["command"] == COMMAND
    assert entry["args"] == ["mcp", "start", "--project-dir", str(project)]
    if "cwd" in client.extra:
        assert entry["cwd"] == str(project)


@pytest.mark.parametrize("client", mc.project_clients(), ids=lambda c: c.name)
def test_merge_keeps_other_servers_and_is_idempotent(client, project):
    if client.format == "toml":
        existing = '[mcp_servers.other]\ncommand = "other"\n\n[profile]\nmodel = "m"\n'
    else:
        existing = json.dumps({client.key: {"other": {"command": "other"}}, "theme": "dark"})
    merged = mc.render(client, project, COMMAND, existing)
    data = _parse(client, merged)
    assert data[client.key]["other"]["command"] == "other"
    assert data[client.key]["kintsugi"]["command"] == COMMAND
    assert mc.render(client, project, COMMAND, merged) is None  # nothing to change


def test_codex_force_replaces_only_our_section(project):
    codex = mc.by_name()["codex"]
    existing = (
        '[mcp_servers.kintsugi]\ncommand = "stale"\n\n'
        '[mcp_servers.kintsugi.env]\nX = "1"\n\n'
        '[mcp_servers.other]\ncommand = "other"\n'
    )
    assert mc.render(codex, project, COMMAND, existing) is None  # present, no force
    data = tomllib.loads(mc.render(codex, project, COMMAND, existing, force=True))
    assert data["mcp_servers"]["kintsugi"]["command"] == COMMAND
    assert "env" not in data["mcp_servers"]["kintsugi"]
    assert data["mcp_servers"]["other"] == {"command": "other"}


def test_gemini_keeps_existing_context_file_names(project):
    gemini = mc.by_name()["gemini-cli"]
    existing = json.dumps({"context": {"fileName": ["GEMINI.md", "NOTES.md"]}})
    data = json.loads(mc.render(gemini, project, COMMAND, existing))
    assert data["context"]["fileName"] == ["GEMINI.md", "NOTES.md", "AGENTS.md"]


def test_unparseable_json_is_left_untouched(project):
    (project / ".cursor").mkdir()
    (project / ".cursor" / "mcp.json").write_text("{not json")
    with pytest.raises(mc.ClientConfigError, match="left unchanged"):
        write_mcp_config(project, ["cursor"])
    assert (project / ".cursor" / "mcp.json").read_text() == "{not json"


def test_user_and_ide_clients_only_print_snippets(project, tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    names = [c.name for c in mc.CLIENTS if c.scope != "project"]

    result = write_mcp_config(project, names)

    assert result["written"] == [] and set(result["snippets"]) == set(names)
    for text in result["snippets"].values():
        assert json.loads(text)["mcpServers"]["kintsugi"]["args"][-1] == str(project.resolve())
    assert list(home.iterdir()) == []
    assert sorted(p.name for p in project.iterdir()) == ["AGENTS.md", PROJECT_CONFIG_FILE]


def test_resolve_clients():
    assert [c.name for c in mc.resolve_clients(None)] == [c.name for c in mc.project_clients()]
    assert mc.resolve_clients(["none"]) == []
    assert [c.name for c in mc.resolve_clients(["zed", "zed"])] == ["zed"]
    with pytest.raises(mc.ClientConfigError, match="choose from"):
        mc.resolve_clients(["nope"])


def test_agents_md_is_created_with_a_current_tool_block(project):
    write_mcp_config(project, ["cursor"])
    text = (project / "AGENTS.md").read_text()
    assert text.startswith("# Spleen: guidance for AI assistants")
    assert render_tool_list() in text
    assert str(project) not in text  # committable: no machine-specific paths


def test_existing_agents_md_keeps_user_text_and_refreshes_the_block(project):
    stale = (
        "# Mine\n\nKeep this.\n\n<!-- mcp-tools:list:begin -->\nold\n<!-- mcp-tools:list:end -->\n"
    )
    (project / "AGENTS.md").write_text(stale)
    write_mcp_config(project, ["cursor"])
    text = (project / "AGENTS.md").read_text()
    assert text.startswith("# Mine\n\nKeep this.") and "old" not in text
    assert render_tool_list() in text

    plain = "# Mine, no generated block\n"
    (project / "AGENTS.md").write_text(plain)
    write_mcp_config(project, ["cursor"])
    assert (project / "AGENTS.md").read_text() == plain


def test_never_writes_files_that_shadow_agents_md(project):
    write_mcp_config(project)
    for name in mc.NEVER_GENERATED:
        assert not (project / name).exists()


def test_init_writes_agents_md_and_the_selected_clients(tmp_path, monkeypatch):
    from click.testing import CliRunner

    from kintsugi.cli import main
    from kintsugi.project import KintsugiProject

    monkeypatch.setattr(KintsugiProject, "setup_notebooks", lambda self, **kw: None)
    runner = CliRunner()

    every = tmp_path / "every"
    assert runner.invoke(main, ["init", str(every), "--force"]).exit_code == 0
    for client in mc.project_clients():
        assert (every / client.path).exists(), client.name
    assert (every / "AGENTS.md").exists()

    one = tmp_path / "one"
    result = runner.invoke(main, ["init", str(one), "--force", "--mcp-client", "cursor"])
    assert result.exit_code == 0, result.output
    assert (one / ".cursor" / "mcp.json").exists() and not (one / ".mcp.json").exists()

    none = tmp_path / "none"
    assert (
        runner.invoke(main, ["init", str(none), "--force", "--mcp-client", "none"]).exit_code == 0
    )
    assert not (none / ".mcp.json").exists() and not (none / "AGENTS.md").exists()

    bad = runner.invoke(main, ["init", str(tmp_path / "bad"), "--mcp-client", "nope"])
    assert bad.exit_code == 2 and not (tmp_path / "bad").exists()


def test_gemini_single_string_file_name_becomes_a_list(project):
    gemini = mc.by_name()["gemini-cli"]
    existing = json.dumps({"context": {"fileName": "GEMINI.md"}})
    data = json.loads(mc.render(gemini, project, COMMAND, existing))
    assert data["context"]["fileName"] == ["GEMINI.md", "AGENTS.md"]


def test_toml_strings_outside_the_bmp_stay_valid(tmp_path):
    codex = mc.by_name()["codex"]
    root = tmp_path / 'proj "\U0001f52c" \\ dir'
    data = tomllib.loads(mc.render(codex, root, COMMAND))
    assert data["mcp_servers"]["kintsugi"]["args"][-1] == str(root)
    assert data["mcp_servers"]["kintsugi"]["cwd"] == str(root)


@pytest.mark.parametrize(
    "header", ['[mcp_servers."kintsugi"]', "[ mcp_servers . 'kintsugi' ]"], ids=["basic", "literal"]
)
def test_codex_recognises_equivalent_header_keys(project, header):
    codex = mc.by_name()["codex"]
    existing = (
        f'{header}\ncommand = "stale"\n\n'
        '[mcp_servers."kintsugi".env]\nX = "1"\n\n'
        '[mcp_servers.other]\ncommand = "other"\n'
    )
    assert mc.render(codex, project, COMMAND, existing) is None  # present, no force
    data = tomllib.loads(mc.render(codex, project, COMMAND, existing, force=True))
    assert data["mcp_servers"]["kintsugi"]["command"] == COMMAND
    assert "env" not in data["mcp_servers"]["kintsugi"]
    assert data["mcp_servers"]["other"] == {"command": "other"}


def test_invalid_existing_toml_is_left_untouched(project):
    (project / ".codex").mkdir()
    (project / ".codex" / "config.toml").write_text("[mcp_servers\nbroken = \n")
    with pytest.raises(mc.ClientConfigError, match="not valid TOML"):
        write_mcp_config(project, ["codex"])
    assert (project / ".codex" / "config.toml").read_text() == "[mcp_servers\nbroken = \n"


def test_init_prints_snippets_for_user_scoped_clients(tmp_path, monkeypatch):
    from click.testing import CliRunner

    from kintsugi.cli import main
    from kintsugi.project import KintsugiProject

    monkeypatch.setattr(KintsugiProject, "setup_notebooks", lambda self, **kw: None)
    root = tmp_path / "desk"
    result = CliRunner().invoke(
        main, ["init", str(root), "--force", "--mcp-client", "claude-desktop"]
    )
    assert result.exit_code == 0, result.output
    assert "Claude Desktop: paste into" in result.output
    assert '"mcpServers"' in result.output and "--project-dir" in result.output
    assert not (root / ".mcp.json").exists()
