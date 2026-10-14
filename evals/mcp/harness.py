"""Evaluation harness for the KINTSUGI MCP tool surface.

This follows the ``mcp-builder`` philosophy: measure tool quality with realistic
tasks instead of guessing. It has two modes.

**Structural** (default, no network / API key):
  * Token cost of the tool definitions injected into every agent context.
  * Schema-quality checks (missing descriptions, terse/overlong descriptions).
  * Flags curated confusable tool groups (consolidation candidates).
  * Confirms every eval task's expected tool exists on the live server.

**Agentic** (``--agentic``, requires ``ANTHROPIC_API_KEY`` + ``anthropic`` SDK):
  * Presents the live tool schemas to Claude for each task with
    ``tool_choice=auto`` and records the model's *first* tool call.
  * Scores tool **selection** (did it pick an acceptable tool?). Because we only
    measure selection, no image data or handler execution is required.

Run::

    python evals/mcp/harness.py                 # structural only
    ANTHROPIC_API_KEY=... python evals/mcp/harness.py --agentic
    python evals/mcp/harness.py --json report.json
"""

from __future__ import annotations

import argparse
import asyncio
import json
import os
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from tasks import TASKS  # noqa: E402  (local sibling module)

# Curated groups whose members are easy to confuse and are candidates for
# consolidation or sharper "Use this when ..." descriptions. Flagged, not failed.
CONFUSABLE_GROUPS: list[tuple[str, list[str]]] = [
    ("denoising", ["denoise", "denoise_advanced"]),
    ("background removal", ["subtract_blank", "gaussian_subtract", "estimate_background"]),
    ("quality", ["assess_quality", "compute_snr"]),
    (
        "parameter suggestion",
        ["suggest_parameters", "suggest_with_learning", "get_learned_parameters"],
    ),
    ("record learning", ["record_successful_parameters", "approve_and_learn"]),
]

DEFAULT_MODEL = "claude-sonnet-4-6"


def estimate_tokens(text: str) -> int:
    """Rough token estimate (~4 chars/token). Approximate; for relative sizing."""
    return max(1, round(len(text) / 4))


def load_tools() -> list[dict[str, Any]]:
    """Build the live MCP server and return its advertised tool definitions.

    Tools are read through the SDK's in-process client, so what the harness
    measures is exactly what a connected Claude Code session receives.
    """
    try:
        from kintsugi.mcp._sdk import MCP_AVAILABLE, unavailable_message
        from kintsugi.mcp.server import create_server
    except Exception as exc:  # pragma: no cover
        raise SystemExit(f"Could not import the MCP server: {exc}")
    if not MCP_AVAILABLE:
        raise SystemExit(unavailable_message())

    async def _list() -> list[Any]:
        from mcp.client import Client

        async with Client(create_server()) as client:
            return list((await client.list_tools()).tools)

    return [
        {
            "name": t.name,
            "title": t.title,
            "description": t.description or "",
            "input_schema": t.input_schema or {},
            "output_schema": t.output_schema,
            "annotations": (
                t.annotations.model_dump(by_alias=True, exclude_none=True)
                if t.annotations
                else None
            ),
        }
        for t in tools_sorted(await_(_list()))
    ]


def await_(coro: Any) -> Any:
    """Run a coroutine to completion from synchronous code."""
    return asyncio.run(coro)


def tools_sorted(tools: list[Any]) -> list[Any]:
    """Stable order for reports."""
    return sorted(tools, key=lambda t: t.name)


def run_structural(tools: list[dict[str, Any]], tasks: list[dict[str, Any]]) -> dict[str, Any]:
    """Compute and print structural metrics for the tool surface."""
    names = {t["name"] for t in tools}

    per_tool_tokens: dict[str, int] = {}
    for t in tools:
        blob = t["name"] + "\n" + t["description"] + "\n" + json.dumps(t["input_schema"])
        per_tool_tokens[t["name"]] = estimate_tokens(blob)
    total_tokens = sum(per_tool_tokens.values())

    # Schema quality.
    missing_tool_desc = [t["name"] for t in tools if len(t["description"].strip()) < 20]
    missing_title = [t["name"] for t in tools if not (t.get("title") or "").strip()]
    missing_annotations = [t["name"] for t in tools if not t.get("annotations")]
    params_missing_desc: list[str] = []
    have_output_schema = 0
    for t in tools:
        if t["output_schema"] is not None:
            have_output_schema += 1
        props = (t["input_schema"] or {}).get("properties", {})
        for pname, pschema in props.items():
            if not (pschema.get("description") or "").strip():
                params_missing_desc.append(f"{t['name']}.{pname}")

    # Task coverage, both directions: tasks that name absent tools, and tools no
    # task ever targets (which the agentic pass can therefore never measure).
    coverage_gaps = [
        {"task": task["id"], "missing": sorted(set(task["expected"]) - names)}
        for task in tasks
        if not set(task["expected"]) & names
    ]
    targeted = set().union(*(set(task["expected"]) for task in tasks)) if tasks else set()
    untargeted_tools = sorted(names - targeted)

    report = {
        "tool_count": len(tools),
        "total_definition_tokens_est": total_tokens,
        "mean_tokens_per_tool_est": round(total_tokens / max(1, len(tools)), 1),
        "tools_with_output_schema": have_output_schema,
        "top_tools_by_tokens": sorted(per_tool_tokens.items(), key=lambda kv: kv[1], reverse=True)[
            :5
        ],
        "tools_missing_description": missing_tool_desc,
        "tools_missing_title": missing_title,
        "tools_missing_annotations": missing_annotations,
        "params_missing_description_count": len(params_missing_desc),
        "untargeted_tools": untargeted_tools,
        "confusable_groups": [
            {"theme": theme, "tools": [m for m in members if m in names]}
            for theme, members in CONFUSABLE_GROUPS
        ],
        "task_coverage_gaps": coverage_gaps,
    }

    print("=" * 70)
    print("STRUCTURAL EVALUATION — KINTSUGI MCP tool surface")
    print("=" * 70)
    print(f"Tools registered:                 {report['tool_count']}")
    print(
        f"Tools with output schema:         {report['tools_with_output_schema']}/{report['tool_count']}"
    )
    print(
        f"Tool-definition tokens (est.):    ~{total_tokens:,} "
        f"(~{report['mean_tokens_per_tool_est']}/tool, injected every call)"
    )
    print(f"Params missing a description:     {report['params_missing_description_count']}")
    print(f"Tools with weak/no description:   {missing_tool_desc or 'none'}")
    print(f"Tools missing a title:            {missing_title or 'none'}")
    print(f"Tools missing annotations:        {missing_annotations or 'none'}")
    print(f"Tools no eval task targets:       {untargeted_tools or 'none'}")
    print("\nHeaviest tools (token est.):")
    for name, tok in report["top_tools_by_tokens"]:
        print(f"   {tok:>5}  {name}")
    print("\nConfusable groups (consolidation / 'Use this when ...' candidates):")
    for grp in report["confusable_groups"]:
        print(f"   [{grp['theme']}]  {', '.join(grp['tools'])}")
    if coverage_gaps:
        print("\nWARNING — eval tasks referencing absent tools:")
        for gap in coverage_gaps:
            print(f"   {gap['task']}: missing {gap['missing']}")
    else:
        print(f"\nTask coverage: all {len(tasks)} tasks reference existing tools. OK")
    return report


def run_agentic(
    tools: list[dict[str, Any]], tasks: list[dict[str, Any]], model: str
) -> dict[str, Any]:
    """Score tool *selection*: does the model pick an acceptable tool per task?"""
    try:
        import anthropic
    except ImportError:
        raise SystemExit("Agentic mode needs the SDK: pip install anthropic")
    if not os.environ.get("ANTHROPIC_API_KEY"):
        raise SystemExit("Agentic mode needs ANTHROPIC_API_KEY in the environment.")

    client = anthropic.Anthropic()
    api_tools = [
        {"name": t["name"], "description": t["description"], "input_schema": t["input_schema"]}
        for t in tools
    ]

    results = []
    passed = 0
    print("\n" + "=" * 70)
    print(f"AGENTIC TOOL-SELECTION EVALUATION  (model={model})")
    print("=" * 70)
    for task in tasks:
        msg = client.messages.create(
            model=model,
            max_tokens=512,
            tools=api_tools,
            tool_choice={"type": "auto"},
            messages=[{"role": "user", "content": task["prompt"]}],
        )
        chosen = next((b.name for b in msg.content if getattr(b, "type", "") == "tool_use"), None)
        ok = chosen in task["expected"]
        passed += int(ok)
        results.append(
            {"task": task["id"], "chosen": chosen, "expected": sorted(task["expected"]), "pass": ok}
        )
        mark = "PASS" if ok else "FAIL"
        print(
            f"   [{mark}] {task['id']:<28} chose={chosen!s:<28} expected={sorted(task['expected'])}"
        )

    accuracy = passed / max(1, len(tasks))
    print(f"\nTool-selection accuracy: {passed}/{len(tasks)} = {accuracy:.0%}")
    return {
        "model": model,
        "accuracy": accuracy,
        "passed": passed,
        "total": len(tasks),
        "results": results,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate the KINTSUGI MCP tool surface.")
    parser.add_argument(
        "--agentic",
        action="store_true",
        help="Run the live tool-selection eval (needs ANTHROPIC_API_KEY).",
    )
    parser.add_argument("--model", default=DEFAULT_MODEL, help=f"Model (default: {DEFAULT_MODEL}).")
    parser.add_argument("--json", metavar="PATH", help="Write the full report as JSON.")
    parser.add_argument(
        "--strict",
        action="store_true",
        help=(
            "Exit non-zero on structural problems (coverage gaps, untargeted tools, "
            "missing descriptions/titles/annotations)."
        ),
    )
    args = parser.parse_args(argv)

    tools = load_tools()
    structural = run_structural(tools, TASKS)
    report: dict[str, Any] = {"structural": structural}
    if args.agentic:
        report["agentic"] = run_agentic(tools, TASKS, args.model)

    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2))
        print(f"\nWrote report to {args.json}")

    if args.strict:
        problems = {
            "task_coverage_gaps": structural["task_coverage_gaps"],
            "tools_missing_description": structural["tools_missing_description"],
            "tools_missing_title": structural["tools_missing_title"],
            "tools_missing_annotations": structural["tools_missing_annotations"],
            "params_missing_description_count": structural["params_missing_description_count"],
            "untargeted_tools": structural["untargeted_tools"],
        }
        failing = {k: v for k, v in problems.items() if v}
        if failing:
            print(f"\nstrict: structural problems detected — failing: {failing}", file=sys.stderr)
            return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
