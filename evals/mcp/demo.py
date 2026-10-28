"""Record a real KINTSUGI MCP session as a Markdown transcript.

Creates a throw-away project with synthetic signal/blank channels, drives the
live server through the SDK's in-process client exactly as an MCP client would,
and prints what every tool returned. The output is deterministic (seeded data,
fixed names, temp paths normalised) so it can be embedded in the docs between
``<!-- mcp-demo:begin -->`` / ``<!-- mcp-demo:end -->`` markers and checked in CI.

Run::

    python evals/mcp/demo.py                       # print the transcript
    python evals/mcp/demo.py --write docs/mcp.md   # splice it into the docs
    python evals/mcp/demo.py --check docs/mcp.md   # fail if the docs are stale
"""

from __future__ import annotations

import argparse
import asyncio
import contextlib
import io
import json
import os
import re
import shutil
import sys
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

BEGIN = "<!-- mcp-demo:begin -->"
END = "<!-- mcp-demo:end -->"

STEPS: list[tuple[str, dict[str, Any]]] = [
    ("list_channels", {"cycle": "1"}),
    ("load_channel", {"cycle": "1", "channel": "CD3"}),
    ("load_channel", {"cycle": "1", "channel": "Blank"}),
    ("get_image_stats", {"channel": "cyc001_CD3"}),
    ("get_thumbnail", {"channel": "cyc001_CD3", "max_size": 64}),
    (
        "suggest_with_learning",
        {
            "channel": "cyc001_CD3",
            "tissue_type": "tonsil",
            "marker_name": "CD3",
            "blank_channel": "cyc001_Blank",
        },
    ),
    (
        "analyze_weighted_subtraction",
        {"signal_channel": "cyc001_CD3", "blank_channel": "cyc001_Blank", "n_ranges": 3},
    ),
    (
        "subtract_blank",
        {
            "signal_channel": "cyc001_CD3",
            "blank_channel": "cyc001_Blank",
            "method": "weighted",
            "n_ranges": 3,
            "tissue_type": "tonsil",
            "marker_name": "CD3",
        },
    ),
    ("assess_quality", {"channel": "cyc001_CD3_subtracted"}),
    ("save_processed", {"channel": "cyc001_CD3_subtracted", "output_name": "CD3_isolated"}),
    (
        "approve_and_learn",
        {
            "channel": "cyc001_CD3_subtracted",
            "tissue_type": "tonsil",
            "marker_name": "CD3",
            "operations_params": {"blank_subtraction": {"method": "weighted", "n_ranges": 3}},
            "quality_before": 0.4,
            "quality_after": 0.8,
        },
    ),
    (
        "get_learned_parameters",
        {"tissue_type": "tonsil", "marker_name": "CD3", "operation": "blank_subtraction"},
    ),
    ("get_session_state", {}),
    ("get_workflow_status", {}),
]


def make_project(root: Path) -> Path:
    """A minimal KINTSUGI project with one cycle of synthetic CD3 + Blank images."""
    from kintsugi.project import KintsugiProject

    with contextlib.redirect_stdout(io.StringIO()):
        project = KintsugiProject.create(root, name="mcp-demo", setup_notebooks=False)
    raw = project.paths.raw / "cyc001"
    raw.mkdir(parents=True, exist_ok=True)

    import tifffile

    rng = np.random.default_rng(42)
    yy, xx = np.mgrid[0:256, 0:256]
    cells = np.zeros((256, 256))
    for cy, cx in rng.integers(16, 240, size=(40, 2)):
        cells += 3000 * np.exp(-((yy - cy) ** 2 + (xx - cx) ** 2) / (2 * 4.0**2))
    autofluorescence = 400 + 300 * np.exp(-((yy - 128) ** 2 + (xx - 128) ** 2) / (2 * 80.0**2))
    noise = rng.normal(0, 40, (256, 256))
    signal = np.clip(cells + autofluorescence + noise, 0, 65535).astype(np.uint16)
    blank = np.clip(autofluorescence + rng.normal(0, 40, (256, 256)), 0, 65535).astype(np.uint16)
    tifffile.imwrite(raw / "CD3.tif", signal)
    tifffile.imwrite(raw / "Blank.tif", blank)
    return project.paths.root


def _normalise(value: Any, root: Path) -> Any:
    """Make results deterministic and compact: round floats, shorten lists, hide temp paths."""
    if isinstance(value, dict):
        return {k: _normalise(v, root) for k, v in value.items()}
    if isinstance(value, list):
        items = [_normalise(v, root) for v in value[:6]]
        return items + ([f"... ({len(value)} items)"] if len(value) > 6 else [])
    if isinstance(value, float):
        return round(value, 3)
    if isinstance(value, str):
        value = value.replace(str(root), "<project>")
        value = re.sub(
            r"data:image/\w+;base64,[A-Za-z0-9+/=]+", "data:image/...;base64,<elided>", value
        )
        value = re.sub(r"\d{8}_\d{6}_[0-9a-f]{6}", "<session>", value)
        value = re.sub(r"\d{4}-\d{2}-\d{2}T[\d:.]+", "<time>", value)
        value = re.sub(r"\d{8}_\d{6}", "<timestamp>", value)
        return value if len(value) <= 200 else value[:200] + "..."
    return value


def _render_result(result: Any, root: Path) -> list[str]:
    lines: list[str] = []
    blocks = [c.type for c in result.content]
    if "image" in blocks:
        image = next(c for c in result.content if c.type == "image")
        lines.append(
            f"content: image block ({image.mime_type}, {len(image.data)} base64 chars) + text"
        )
    structured = _normalise(result.structured_content, root)
    if result.is_error:
        lines.append("isError: true")
        text = result.content[0].text
        lines.append(text[text.index("{") :] if "{" in text else text)
    else:
        lines.append("```json")
        lines.append(json.dumps(structured, indent=2, default=str))
        lines.append("```")
    return lines


async def record(root: Path) -> str:
    from mcp.client import Client

    from kintsugi.mcp.server import create_server
    from kintsugi.mcp.session import clear_session
    from kintsugi.mcp.tools import learning
    from kintsugi.mcp.workflow_store import clear_workflow_stores

    # Fresh session, learning engines and workflow stores (they are cached per project,
    # and the user-wide engine would otherwise keep pointing at a previous run's HOME).
    clear_session()
    learning._learning_engines.clear()
    clear_workflow_stores()
    out: list[str] = []
    async with Client(create_server()) as client:
        out.append(
            f"Server `{client.server_info.name}` {client.server_info.version}, "
            f"{len((await client.list_tools()).tools)} tools, project `<project>`."
        )
        out.append("")
        for i, (tool, args) in enumerate(STEPS, 1):
            out.append(f"**{i}. `{tool}`** `{json.dumps(args)}`")
            out.append("")
            result = await client.call_tool(tool, args)
            out.extend(_render_result(result, root))
            out.append("")
    clear_session()
    learning._learning_engines.clear()
    clear_workflow_stores()
    return "\n".join(out).rstrip() + "\n"


def run() -> str:
    """Record the transcript in an isolated HOME/project and return it."""
    tmp = Path(tempfile.mkdtemp(prefix="kintsugi-mcp-demo-"))
    from kintsugi.mcp.session import set_project_root

    env_backup = {k: os.environ.get(k) for k in ("HOME", "USERPROFILE", "KINTSUGI_PROJECT_DIR")}
    try:
        home = tmp / "home"
        home.mkdir()
        os.environ["HOME"] = str(home)
        os.environ["USERPROFILE"] = str(home)
        root = make_project(tmp / "demo_project")
        os.environ.pop("KINTSUGI_PROJECT_DIR", None)
        # What `kintsugi mcp start --project-dir <root>` does.
        set_project_root(root)
        return asyncio.run(record(root))
    finally:
        set_project_root(None)
        for key, value in env_backup.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value
        shutil.rmtree(tmp, ignore_errors=True)


def _pattern() -> re.Pattern[str]:
    return re.compile(rf"({re.escape(BEGIN)}\n)(.*?)(\n?{re.escape(END)})", re.DOTALL)


def splice(path: Path, transcript: str) -> bool:
    """Replace the demo block in ``path``; return True if the file changed."""
    text = path.read_text()
    if not _pattern().search(text):
        raise SystemExit(f"{path}: no {BEGIN} / {END} markers")
    new = _pattern().sub(lambda m: f"{m.group(1)}{transcript.rstrip()}\n{END}", text, count=1)
    changed = new != text
    if changed:
        path.write_text(new)
    return changed


def is_current(path: Path, transcript: str) -> bool:
    match = _pattern().search(path.read_text())
    return bool(match) and match.group(2).strip() == transcript.strip()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--write", metavar="FILE", help="Splice the transcript into FILE.")
    parser.add_argument("--check", metavar="FILE", help="Fail if FILE's transcript is stale.")
    args = parser.parse_args(argv)

    transcript = run()
    if args.check:
        path = Path(args.check)
        if not is_current(path, transcript):
            print(
                f"{path}: recorded MCP transcript is stale; run: "
                f"python evals/mcp/demo.py --write {path}",
                file=sys.stderr,
            )
            return 1
        print(f"{path}: recorded MCP transcript is current")
        return 0
    if args.write:
        changed = splice(Path(args.write), transcript)
        print(f"{args.write}: {'updated' if changed else 'already current'}")
        return 0
    print(transcript)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
