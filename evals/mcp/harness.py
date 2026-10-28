"""Evaluation harness for the KINTSUGI MCP tool surface.

This follows the ``mcp-builder`` philosophy: measure tool quality with realistic
tasks instead of guessing. It has two modes.

**Structural** (default, no network / API key):
  * Token cost of the tool definitions injected into every agent context.
  * Schema-quality checks (missing descriptions, terse/overlong descriptions).
  * Flags curated confusable tool groups (consolidation candidates).
  * Confirms every eval task's expected tool exists on the live server.

**Agentic** (``--agentic --provider NAME --model MODEL``, see ``providers.py``):
  * Presents the live tool schemas to the model for each task with automatic
    tool choice and records the model's *first* tool call.
  * Scores tool **selection** (did it pick an acceptable tool?). Because we only
    measure selection, no image data or handler execution is required.
  * ``--provider anthropic`` uses the Anthropic Messages API; ``--provider openai``
    uses the OpenAI Chat Completions API, which other platforms and local servers
    (Ollama, vLLM) also serve (``OPENAI_BASE_URL`` / ``--base-url``).
  * There is no default provider or model. In a terminal the harness asks for a
    missing one; otherwise it exits with a usage error.

Run::

    python evals/mcp/harness.py                 # structural only
    ANTHROPIC_API_KEY=... python evals/mcp/harness.py --agentic --provider anthropic --model MODEL
    OPENAI_API_KEY=... python evals/mcp/harness.py --agentic --provider openai --model MODEL
    OPENAI_BASE_URL=http://localhost:11434/v1 python evals/mcp/harness.py \
        --agentic --provider openai --model MODEL      # local model (Ollama)
    python evals/mcp/harness.py --json report.json
"""

from __future__ import annotations

import argparse
import asyncio
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
import providers  # noqa: E402  (local sibling module)
from tasks import TASKS  # noqa: E402  (local sibling module)

# Curated groups whose members are easy to confuse and are candidates for
# consolidation or sharper "Use this when ..." descriptions. Flagged, not failed.
CONFUSABLE_GROUPS: list[tuple[str, list[str]]] = [
    ("denoising", ["denoise", "denoise_advanced"]),
    ("background removal", ["subtract_blank", "gaussian_subtract", "estimate_background"]),
    ("quality", ["assess_quality", "compute_snr", "batch_assess_quality"]),
    (
        "parameter suggestion",
        ["suggest_parameters", "suggest_with_learning", "get_learned_parameters"],
    ),
    ("record learning", ["record_successful_parameters", "approve_and_learn"]),
    (
        "session vs project state",
        ["get_session_state", "get_workflow_status", "get_processing_history"],
    ),
    ("review outcome", ["approve_and_learn", "set_review_decision"]),
]


def estimate_tokens(text: str) -> int:
    """Rough token estimate (~4 chars/token). Approximate; for relative sizing."""
    return max(1, round(len(text) / 4))


def load_tools() -> list[dict[str, Any]]:
    """Build the live MCP server and return its advertised tool definitions.

    Tools are read through the SDK's in-process client, so what the harness
    measures is exactly what a connected MCP client receives.
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
    tools: list[dict[str, Any]],
    tasks: list[dict[str, Any]],
    *,
    provider: str,
    model: str,
    base_url: str | None = None,
) -> dict[str, Any]:
    """Score tool *selection*: does the model pick an acceptable tool per task?"""
    host = providers.endpoint_host(provider, base_url)
    results = []
    passed = 0
    print("\n" + "=" * 70)
    print(f"AGENTIC TOOL-SELECTION EVALUATION  (provider={provider}, model={model}, host={host})")
    print("=" * 70)
    for task in tasks:
        try:
            chosen = providers.select_tool(
                provider, model, task["prompt"], tools, base_url=base_url
            )
        except providers.ProviderError as exc:
            raise SystemExit(f"Agentic eval stopped at task {task['id']}: {exc}") from None
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
        "provider": provider,
        "model": model,
        "host": host,
        "accuracy": accuracy,
        "passed": passed,
        "total": len(tasks),
        "results": results,
    }


DOCS_WITH_TOOL_BLOCKS = ("docs/mcp.md", "docs/cli.md", "README.md", "AGENTS.md")


def _docs_drift() -> list[str]:
    """Documented tool tables that no longer match the registry (see kintsugi.mcp.docs_sync)."""
    from kintsugi.mcp import docs_sync

    repo = Path(__file__).resolve().parents[2]
    present = [repo / rel for rel in DOCS_WITH_TOOL_BLOCKS if (repo / rel).exists()]
    return docs_sync.check_files(present)


def _interactive() -> bool:
    """True when a person can answer prompts (stdin and stdout are terminals)."""
    return sys.stdin.isatty() and sys.stdout.isatty()


def _choose_provider_and_model(args: argparse.Namespace) -> None:
    """Fill ``args.provider`` and ``args.model``: never a default, asked for in a terminal."""
    if args.provider and args.model:
        return
    usage = (
        f"--agentic requires --provider {{{','.join(providers.PROVIDERS)}}} and --model MODEL "
        "(there is no default provider or model)"
    )
    if not _interactive():
        raise ValueError(usage)
    try:
        while not args.provider:
            answer = input(f"Model provider ({', '.join(providers.PROVIDERS)}): ").strip()
            if answer in providers.PROVIDERS:
                args.provider = answer
            elif not answer:
                raise ValueError(usage)
            else:
                print(f"Unknown provider {answer!r}.")
        if not args.model:
            args.model = input(f"Model identifier at {args.provider}: ").strip()
            if not args.model:
                raise ValueError(usage)
    except EOFError:
        raise ValueError(usage) from None


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Evaluate the KINTSUGI MCP tool surface.")
    parser.add_argument(
        "--agentic",
        action="store_true",
        help="Run the live tool-selection eval against a model (needs --provider and --model).",
    )
    parser.add_argument(
        "--provider",
        choices=providers.PROVIDERS,
        help=(
            "API format: 'anthropic' (Messages API) or 'openai' (Chat Completions, also "
            "served by other platforms and local servers). Required with --agentic; "
            "asked for in a terminal."
        ),
    )
    parser.add_argument(
        "--model",
        help="Model identifier at the provider. Required with --agentic; asked for in a terminal.",
    )
    parser.add_argument(
        "--base-url",
        help="Endpoint base URL (default: ANTHROPIC_BASE_URL / OPENAI_BASE_URL, else the vendor's).",
    )
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
    if args.agentic:
        try:
            _choose_provider_and_model(args)
        except ValueError as exc:
            parser.error(str(exc))

    tools = load_tools()
    structural = run_structural(tools, TASKS)
    report: dict[str, Any] = {"structural": structural}
    if args.agentic:
        report["agentic"] = run_agentic(
            tools, TASKS, provider=args.provider, model=args.model, base_url=args.base_url
        )

    if args.json:
        Path(args.json).write_text(json.dumps(report, indent=2))
        print(f"\nWrote report to {args.json}")

    if args.strict:
        docs_drift = _docs_drift()
        problems = {
            "docs_drift": docs_drift,
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
