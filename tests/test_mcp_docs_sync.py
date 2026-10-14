"""Tests for the generated MCP tool reference (docs_sync + `kintsugi mcp tools`)."""

from __future__ import annotations

import pytest
from click.testing import CliRunner

from kintsugi.mcp import docs_sync
from kintsugi.mcp.tool_specs import TOOL_SPECS

TEMPLATE = """# Heading

Intro text.

<!-- mcp-tools:table:begin -->
stale
<!-- mcp-tools:table:end -->

Outro text.
"""


def test_table_lists_every_tool():
    table = docs_sync.render_tool_table()
    for spec in TOOL_SPECS:
        assert f"`{spec['name']}`" in table
        assert spec["title"] in table
    assert f"{len(TOOL_SPECS)} tools." in table


def test_list_groups_every_tool():
    listing = docs_sync.render_tool_list()
    for spec in TOOL_SPECS:
        assert f"`{spec['name']}`" in listing
    assert "**Parameter learning:**" in listing


def test_splice_is_idempotent_and_preserves_surroundings():
    once = docs_sync.splice_block(TEMPLATE, "table")
    twice = docs_sync.splice_block(once, "table")
    assert once == twice
    assert once.startswith("# Heading\n\nIntro text.\n\n<!-- mcp-tools:table:begin -->\n")
    assert once.endswith("<!-- mcp-tools:table:end -->\n\nOutro text.\n")
    assert "stale" not in once
    assert "`load_channel`" in once


def test_splice_requires_markers():
    with pytest.raises(docs_sync.DocsSyncError):
        docs_sync.splice_block("no markers here", "table")


def test_sync_file_reports_stale_then_current(tmp_path):
    doc = tmp_path / "doc.md"
    doc.write_text(TEMPLATE)
    assert docs_sync.sync_file(doc, write=False) == {"table": False}
    assert "stale" in doc.read_text()

    assert docs_sync.sync_file(doc, write=True) == {"table": False}
    assert docs_sync.sync_file(doc, write=False) == {"table": True}
    assert docs_sync.check_files([doc]) == []

    # An edited row counts as drift.
    doc.write_text(doc.read_text().replace("`load_channel`", "`load_chanel`", 1))
    assert docs_sync.check_files([doc]) == [f"{doc}: mcp-tools:table block is stale"]


def test_cli_markdown_write_and_check(tmp_path):
    from kintsugi.cli import main

    runner = CliRunner()
    result = runner.invoke(main, ["mcp", "tools", "--markdown"])
    assert result.exit_code == 0, result.output
    assert "| `load_channel` |" in result.output

    doc = tmp_path / "doc.md"
    doc.write_text(TEMPLATE)
    result = runner.invoke(main, ["mcp", "tools", "--check", str(doc)])
    assert result.exit_code == 1
    assert "stale" in result.output

    result = runner.invoke(main, ["mcp", "tools", "--write", str(doc)])
    assert result.exit_code == 0, result.output
    assert "`load_channel`" in doc.read_text()

    result = runner.invoke(main, ["mcp", "tools", "--check", str(doc)])
    assert result.exit_code == 0, result.output


def test_cli_table_lists_titles():
    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["mcp", "tools"])
    assert result.exit_code == 0, result.output
    assert "load_channel" in result.output
    assert f"{len(TOOL_SPECS)} tools" in result.output
