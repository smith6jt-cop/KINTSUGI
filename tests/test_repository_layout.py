"""The repository ships assistant guidance for every platform and no client configuration."""

from __future__ import annotations

import subprocess
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
GUIDANCE_DIRS = ("", "src/kintsugi", "src/kintsugi/signal", "workflow", "notebooks")


def _tracked(*patterns: str) -> list[str]:
    try:
        out = subprocess.run(
            ["git", "ls-files", "--", *patterns],
            cwd=REPO,
            capture_output=True,
            text=True,
            check=True,
        ).stdout
    except (OSError, subprocess.CalledProcessError):
        pytest.skip("not a git checkout")
    return out.split()


@pytest.mark.parametrize("directory", GUIDANCE_DIRS)
def test_agents_md_is_canonical_and_claude_md_imports_it(directory):
    folder = REPO / directory
    assert (folder / "AGENTS.md").read_text().startswith("#")
    assert (folder / "CLAUDE.md").read_text() == "@AGENTS.md\n"


def test_root_agents_md_names_the_claude_md_stub_and_no_duplicate_files():
    """Guards against a blind CLAUDE.md -> AGENTS.md replacement (it broke both once)."""
    import re

    text = (REPO / "AGENTS.md").read_text()
    assert "Claude Code reads it through the one-line `CLAUDE.md` (`@AGENTS.md`)" in text
    for line in text.splitlines():
        files = re.findall(r"(?<![\w/.-])[\w./-]+\.md\b", line)
        if line.lstrip().startswith("kintsugi "):
            assert len(files) == len(set(files)), line


def test_repository_tracks_no_mcp_client_configuration():
    """Client files are personal: contributors keep their own, untracked."""
    tracked = _tracked(
        ".mcp.json",
        ".claude",
        ".cursor",
        ".gemini",
        ".agents",
        ".codex",
        ".zed",
        ".windsurf",
        ".devin",
        ".vscode/mcp.json",
    )
    assert tracked == []
