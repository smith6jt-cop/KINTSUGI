"""Configuration for the MCP clients that can run the KINTSUGI server.

Every client starts the same stdio server, ``kintsugi mcp start --project-dir
<project>``, with the absolute path of the ``kintsugi`` executable: clients
differ in the working directory and environment they give a server (several
ignore a ``cwd`` field, some pass on almost no environment variables, GUI clients
may not inherit ``PATH``), so the project travels in ``args``.

Clients whose configuration lives in the project (``scope="project"``) get a
file written or merged into the project; clients configured per user or in an
IDE dialog get a snippet to paste, and nothing is ever written to the user's
home directory. Existing files are merged: other servers and settings are kept,
and the ``kintsugi`` entry is replaced only with ``force``. Only the standard
library and the tool registry are imported, so this works without the ``mcp``
extra.
"""

from __future__ import annotations

import json
import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Literal

SERVER_NAME = "kintsugi"
AGENTS_FILE = "AGENTS.md"


class ClientConfigError(ValueError):
    """Raised for an unknown client name or a file that cannot be merged safely."""


@dataclass(frozen=True)
class ClientSpec:
    """How one MCP client is configured to start the KINTSUGI server."""

    name: str
    display: str
    scope: Literal["project", "user", "ide"]
    # Project scope: path relative to the project root. Otherwise: where the user's
    # configuration lives, as text for the snippet header.
    path: str
    key: str = "mcpServers"
    format: Literal["json", "toml"] = "json"
    # Fields added to the server entry beyond command/args; "{root}" is substituted.
    extra: dict[str, Any] = field(default_factory=dict)
    # Extra top-level settings merged into the file (lists are merged by union).
    settings: dict[str, Any] = field(default_factory=dict)
    instructions: str = AGENTS_FILE
    docs_url: str = ""
    next_steps: str = ""


CLIENTS: tuple[ClientSpec, ...] = (
    ClientSpec(
        name="claude-code",
        display="Claude Code",
        scope="project",
        path=".mcp.json",
        extra={"type": "stdio"},
        instructions="CLAUDE.md, else AGENTS.md",
        docs_url="https://code.claude.com/docs/en/mcp",
        next_steps="Open the project folder, accept the trust dialog, then run /mcp.",
    ),
    ClientSpec(
        name="vscode",
        display="VS Code (GitHub Copilot)",
        scope="project",
        path=".vscode/mcp.json",
        key="servers",
        extra={"type": "stdio", "cwd": "{root}"},
        docs_url="https://code.visualstudio.com/docs/copilot/customization/mcp-servers",
        next_steps="Open the folder, trust the workspace, and start 'kintsugi' from "
        "MCP: List Servers.",
    ),
    ClientSpec(
        name="cursor",
        display="Cursor",
        scope="project",
        path=".cursor/mcp.json",
        extra={"type": "stdio"},
        docs_url="https://cursor.com/docs/context/mcp",
        next_steps="Open the folder; enable 'kintsugi' under Settings > MCP if prompted.",
    ),
    ClientSpec(
        name="gemini-cli",
        display="Gemini CLI",
        scope="project",
        path=".gemini/settings.json",
        extra={"cwd": "{root}"},
        settings={"context": {"fileName": ["AGENTS.md", "GEMINI.md"]}},
        docs_url="https://geminicli.com/docs/tools/mcp-server/",
        next_steps="Run 'gemini' in the project folder (trust it with 'gemini trust' if asked), "
        "then /mcp.",
    ),
    ClientSpec(
        name="antigravity",
        display="Antigravity",
        scope="project",
        path=".agents/mcp_config.json",
        extra={"cwd": "{root}"},
        docs_url="https://antigravity.google/docs/mcp",
        next_steps="Open the folder as the workspace; the server appears under MCP servers.",
    ),
    ClientSpec(
        name="codex",
        display="OpenAI Codex CLI",
        scope="project",
        path=".codex/config.toml",
        key="mcp_servers",
        format="toml",
        extra={"cwd": "{root}", "startup_timeout_sec": 60, "tool_timeout_sec": 600},
        docs_url="https://developers.openai.com/codex/mcp",
        next_steps="Mark the project trusted (Codex loads .codex/config.toml only in trusted "
        "projects), then run 'codex' there and /mcp.",
    ),
    ClientSpec(
        name="zed",
        display="Zed",
        scope="project",
        path=".zed/settings.json",
        key="context_servers",
        docs_url="https://zed.dev/docs/ai/mcp",
        next_steps="Open the folder (trust the worktree); the server is listed in the "
        "Agent Panel settings.",
    ),
    ClientSpec(
        name="claude-desktop",
        display="Claude Desktop",
        scope="user",
        path="claude_desktop_config.json (macOS: ~/Library/Application Support/Claude/; "
        "Windows: %APPDATA%\\Claude\\)",
        instructions="(none)",
        docs_url="https://modelcontextprotocol.io/docs/develop/connect-local-servers",
        next_steps="Merge the entry into mcpServers and restart Claude Desktop.",
    ),
    ClientSpec(
        name="windsurf",
        display="Windsurf / Devin Desktop",
        scope="user",
        path="mcp_config.json (~/.config/devin/ or ~/.codeium/; Windows: %APPDATA%\\devin\\)",
        docs_url="https://docs.windsurf.com/windsurf/cascade/mcp",
        next_steps="Merge the entry into mcpServers, then refresh the MCP servers list.",
    ),
    ClientSpec(
        name="cline",
        display="Cline",
        scope="user",
        path="cline_mcp_settings.json (open it from Cline > MCP Servers > Configure)",
        docs_url="https://docs.cline.bot/mcp/configuring-mcp-servers",
        next_steps="Merge the entry into mcpServers and save; Cline restarts the server.",
    ),
    ClientSpec(
        name="jetbrains",
        display="JetBrains AI Assistant",
        scope="ide",
        path="Settings > Tools > AI Assistant > Model Context Protocol (MCP) > Add",
        docs_url="https://www.jetbrains.com/help/ai-assistant/mcp.html",
        next_steps="Choose STDIO, paste the JSON, and set the working directory to the "
        "project folder.",
    ),
)

# Instruction files other tools prefer over AGENTS.md; never generated, because
# creating one would hide AGENTS.md from clients that read only the first match.
NEVER_GENERATED = (".rules", ".cursorrules", ".windsurfrules", ".clinerules")


def by_name() -> dict[str, ClientSpec]:
    return {c.name: c for c in CLIENTS}


def project_clients() -> list[ClientSpec]:
    return [c for c in CLIENTS if c.scope == "project"]


def resolve_clients(names: Iterable[str] | None) -> list[ClientSpec]:
    """``None`` means every project-scoped client; ``"none"`` (alone) means no client."""
    if names is None:
        return project_clients()
    names = list(names)
    if names == ["none"]:
        return []
    registry = by_name()
    unknown = [n for n in names if n not in registry]
    if unknown:
        raise ClientConfigError(
            f"unknown MCP client {', '.join(map(repr, unknown))}; "
            f"choose from {', '.join(registry)} or 'none'"
        )
    return [registry[n] for n in dict.fromkeys(names)]


def server_entry(client: ClientSpec, root: Path, command: str) -> dict[str, Any]:
    """The ``kintsugi`` server entry for ``client``."""
    entry: dict[str, Any] = {}
    if "type" in client.extra:
        entry["type"] = client.extra["type"]
    entry["command"] = command
    entry["args"] = ["mcp", "start", "--project-dir", str(root)]
    for key, value in client.extra.items():
        if key != "type":
            entry[key] = str(root) if value == "{root}" else value
    return entry


# --------------------------------------------------------------------------- JSON


def _union(existing: Any, wanted: Any) -> Any:
    """Merge ``wanted`` into ``existing``: dicts recursively, lists by union."""
    if isinstance(existing, dict) and isinstance(wanted, dict):
        merged = dict(existing)
        for key, value in wanted.items():
            merged[key] = _union(existing[key], value) if key in existing else value
        return merged
    if isinstance(wanted, list):
        # A setting that accepts "a string or a list" (Gemini's context.fileName)
        # may hold a single string; keep it as the first element.
        if isinstance(existing, str):
            existing = [existing]
        if isinstance(existing, list):
            return existing + [v for v in wanted if v not in existing]
    return existing


def _merge_json(
    client: ClientSpec, text: str | None, entry: dict[str, Any], force: bool, path: Path
) -> str | None:
    """Return the new file text, or ``None`` when nothing changes."""
    try:
        data = json.loads(text) if text and text.strip() else {}
    except json.JSONDecodeError as exc:
        raise ClientConfigError(f"{path} is not valid JSON ({exc}); left unchanged") from exc
    if not isinstance(data, dict):
        raise ClientConfigError(f"{path} does not hold a JSON object; left unchanged")
    original = json.dumps(data, sort_keys=True)
    servers = data.setdefault(client.key, {})
    if not isinstance(servers, dict):
        raise ClientConfigError(f"{path}: {client.key!r} is not an object; left unchanged")
    if force or SERVER_NAME not in servers:
        servers[SERVER_NAME] = entry
    data = _union(data, client.settings) if client.settings else data
    if text is not None and json.dumps(data, sort_keys=True) == original:
        return None
    return json.dumps(data, indent=2) + "\n"


# --------------------------------------------------------------------------- TOML

_TOML_HEADER = re.compile(r"^\s*\[\[?(.+?)\]\]?\s*(?:#.*)?$")
_TOML_KEY_PART = re.compile(r"""\s*(?:"((?:[^"\\]|\\.)*)"|'([^']*)'|([A-Za-z0-9_-]+))\s*(\.|$)""")

try:  # Python 3.11+
    import tomllib as _tomllib
except ModuleNotFoundError:  # pragma: no cover - Python 3.10 without tomli
    try:
        import tomli as _tomllib  # type: ignore[no-redef]
    except ModuleNotFoundError:
        _tomllib = None  # type: ignore[assignment]


def _toml_string(value: str) -> str:
    """A TOML basic string: escape quote, backslash and control characters only.

    Everything else (including characters outside the BMP) stays literal UTF-8;
    JSON-style surrogate-pair escapes are not valid TOML.
    """
    out = []
    for ch in value:
        if ch in ('"', "\\"):
            out.append("\\" + ch)
        elif ord(ch) < 0x20 or ord(ch) == 0x7F:
            out.append(f"\\u{ord(ch):04X}")
        else:
            out.append(ch)
    return '"' + "".join(out) + '"'


def _toml_value(value: Any) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, int | float):
        return str(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(v) for v in value) + "]"
    return _toml_string(str(value))


def _toml_key_parts(key: str) -> tuple[str, ...] | None:
    """Split a TOML dotted key into its parts (quotes removed); None if unparseable."""
    parts: list[str] = []
    pos = 0
    while pos < len(key):
        match = _TOML_KEY_PART.match(key, pos)
        if not match or match.end() == pos:
            return None
        basic, literal, bare, _sep = match.groups()
        if basic is not None:
            try:
                parts.append(json.loads(f'"{basic}"'))
            except json.JSONDecodeError:
                return None
        else:
            parts.append(literal if literal is not None else bare)
        pos = match.end()
    return tuple(parts) or None


def _toml_check(text: str, what: str) -> None:
    """Raise ClientConfigError if ``text`` is not valid TOML (when a parser is available)."""
    if _tomllib is None:
        return
    try:
        _tomllib.loads(text)
    except _tomllib.TOMLDecodeError as exc:
        raise ClientConfigError(f"{what} is not valid TOML ({exc}); left unchanged") from exc


def _render_toml(client: ClientSpec, entry: dict[str, Any]) -> str:
    lines = [f"[{client.key}.{SERVER_NAME}]"]
    lines += [f"{key} = {_toml_value(value)}" for key, value in entry.items()]
    return "\n".join(lines) + "\n"


def _merge_toml(
    client: ClientSpec, text: str | None, entry: dict[str, Any], force: bool
) -> str | None:
    section = _render_toml(client, entry)
    if not text:
        return section
    _toml_check(text, f"the existing {client.path}")
    ours = (client.key, SERVER_NAME)
    lines = text.splitlines(keepends=True)
    start = end = None
    for i, line in enumerate(lines):
        match = _TOML_HEADER.match(line)
        if not match:
            continue
        parts = _toml_key_parts(match.group(1))
        if start is None and parts == ours:
            start = i
        elif start is not None and (parts is None or parts[: len(ours)] != ours):
            end = i
            break
    if start is None:
        sep = "" if text.endswith("\n\n") else ("\n" if text.endswith("\n") else "\n\n")
        merged = text + sep + section
    elif not force:
        return None
    else:
        end = len(lines) if end is None else end
        tail = "".join(lines[end:])
        merged = "".join(lines[:start]) + section + ("\n" + tail if tail else "")
    _toml_check(merged, f"the merged {client.path}")
    return merged


# --------------------------------------------------------------------------- API


def render(
    client: ClientSpec,
    root: Path,
    command: str,
    existing: str | None = None,
    *,
    force: bool = False,
) -> str | None:
    """Content for ``client``'s file (merged into ``existing``), or ``None`` if unchanged."""
    entry = server_entry(client, root, command)
    if client.format == "toml":
        return _merge_toml(client, existing, entry, force)
    return _merge_json(client, existing, entry, force, root / client.path)


def snippet(client: ClientSpec, root: Path, command: str) -> str:
    """Text a user pastes into a user- or IDE-scoped configuration."""
    entry = server_entry(client, root, command)
    if client.format == "toml":
        return _render_toml(client, entry)
    return json.dumps({client.key: {SERVER_NAME: entry}}, indent=2)


def write_project_configs(
    root: Path, clients: Sequence[ClientSpec], command: str, *, force: bool = False
) -> tuple[list[Path], list[Path]]:
    """Write or merge each project-scoped client file; return ``(written, unchanged)``."""
    written: list[Path] = []
    unchanged: list[Path] = []
    for client in clients:
        if client.scope != "project":
            continue
        path = root / client.path
        existing = path.read_text() if path.exists() else None
        content = render(client, root, command, existing, force=force)
        if content is None:
            unchanged.append(path)
            continue
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        written.append(path)
    return written, unchanged


# --------------------------------------------------------------------------- AGENTS.md

AGENTS_TEMPLATE = """\
# {name}: guidance for AI assistants

This is a KINTSUGI project (multiplex immunofluorescence image processing). The
KINTSUGI MCP server gives an assistant tools to load channels, isolate signal,
assess quality and learn parameters for this project.

## The MCP server

`kintsugi mcp config .` writes the MCP client configuration for this project (run
it again after moving the project or the Python environment; collaborators run it
on their own machine, because the generated files hold absolute paths and are not
committed). Every client starts the same server:

```bash
kintsugi mcp start --project-dir /path/to/this/project
```

`project_path` arguments default to this project. See
https://github.com/smith6jt-cop/KINTSUGI/blob/main/docs/mcp.md for every tool,
resource and the security model.

## Tools

<!-- mcp-tools:list:begin -->
<!-- mcp-tools:list:end -->

## Signal isolation, one marker at a time

1. `list_channels`, then `load_channel` for the marker and its blank.
2. `suggest_with_learning` for starting parameters; preview weighted subtraction with
   `analyze_weighted_subtraction`.
3. Apply `subtract_blank`, then `denoise` / `apply_clahe` / `clean_background` as needed;
   check each step with `assess_quality` and `get_thumbnail`.
4. Show the result to the user; `save_processed` only after they agree.
5. `approve_and_learn` records what worked so the next marker starts closer.

## Layout

- `data/raw/` input cycles (never modified), `data/processed/` outputs.
- `meta/experiment.json` microscope parameters, `meta/CHANNELNAMES.txt` marker names.
- `.kintsugi/` project state (learning database, workflow state); not committed.
"""


def write_agents_md(root: Path, name: str) -> Path | None:
    """Create ``AGENTS.md``, or refresh only its tool block; return the path if written."""
    from kintsugi.mcp.docs_sync import find_blocks, splice_block

    path = root / AGENTS_FILE
    if path.exists():
        text = path.read_text()
        if "list" not in find_blocks(text):
            return None  # the user's own file: never replaced
        updated = splice_block(text, "list")
    else:
        updated = splice_block(AGENTS_TEMPLATE.format(name=name), "list")
        text = None
    if updated == text:
        return None
    path.write_text(updated)
    return path
