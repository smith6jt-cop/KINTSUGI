# KINTSUGI MCP — Evaluation Harness

Measures the quality of the KINTSUGI MCP **tool surface** with realistic tasks
instead of intuition, following the `mcp-builder` philosophy. Use it to drive and
verify MCP improvements (tool consolidation, description quality, resources/prompts).

## Why

The server exposes 32 tools to Claude. Two questions decide whether that surface
is *good*:

1. **Is it cheap and clear?** (structural) — how many tokens do the tool
   definitions cost, and are any tools/params under-described, untitled,
   unannotated or confusable?
2. **Does the model use it correctly?** (agentic) — given a realistic request,
   does Claude pick the right tool?

## Run it

```bash
# Structural only — no network, no API key. Needs mcp>=2.0 (pip install -e ".[claude]").
python evals/mcp/harness.py

# Same, but exit non-zero on any structural problem (the CI gate).
python evals/mcp/harness.py --strict

# Agentic tool-selection eval — needs an API key + the anthropic SDK.
pip install anthropic
ANTHROPIC_API_KEY=... python evals/mcp/harness.py --agentic

# Save the full report (structural + agentic) as JSON.
python evals/mcp/harness.py --agentic --json mcp_eval_report.json
```

The harness reads the tools **through the in-process MCP client**, so it sees
exactly what Claude Code sees (titles, annotations, the signature-derived input
schemas, output schemas). The agentic mode measures **tool selection only** — it
records the model's first `tool_use` per task and never executes the handler — so
it needs **no image data or project**, just the live tool schemas.

## Files

| File | Purpose |
|------|---------|
| `tasks.py` | 29 realistic natural-language tasks with `expected` acceptable tool(s); every registered tool is the target of at least one task. Includes deliberate confusable pairs (`denoise` vs `denoise_advanced`; the three background-removal tools; `assess_quality` vs `compute_snr`; the parameter-suggestion and record-learning groups). |
| `harness.py` | Builds the live server, runs the structural metrics, and (opt-in) the agentic tool-selection loop. |

Add a task by appending an `EvalTask` to `TASKS`; the structural pass verifies that
its `expected` tools exist on the server, and `--strict` fails if any tool is left
without a task.

## Strict checks (`--strict`)

`--strict` exits non-zero when any of these is non-empty, which makes it a
deterministic PR gate (no API key needed; `.github/workflows/mcp-eval.yml`):

| Check | Fails when |
|-------|-----------|
| `task_coverage_gaps` | a task's `expected` tool does not exist on the server |
| `tools_missing_description` | a tool has no or a trivially short description |
| `tools_missing_title` | a tool is advertised without a `title` |
| `tools_missing_annotations` | a tool carries no behaviour annotations (readOnly/destructive/idempotent/openWorld) |
| `params_missing_description_count` | any input parameter lacks a description |
| `untargeted_tools` | a registered tool is not the target of any eval task |

### Docs-drift check

The documented tool surface is generated, not hand-written: `kintsugi mcp tools
--write FILE...` renders the registry (`kintsugi.mcp.tool_specs.TOOL_SPECS`) into
the `<!-- mcp-tools:table -->` / `<!-- mcp-tools:list -->` blocks of `docs/mcp.md`,
`docs/cli.md`, `README.md` and `CLAUDE.md`, and `kintsugi mcp tools --check FILE...`
exits non-zero when a block is stale. CI runs the check, so adding, renaming or
re-describing a tool without regenerating the docs fails the build
(`tests/test_mcp_docs_sync.py` covers the generator itself).

## Findings (baseline, structural pass)

`python evals/mcp/harness.py --json baseline.json` on the current server:

- **32 tools**, **32/32 carry an output schema** (structured output on every tool;
  `get_image_stats`, `assess_quality`, `suggest_with_learning` and
  `get_session_state` declare typed schemas).
- **~8.0k tokens** of tool definitions injected on every call (~250/tool). This is
  up from ~3.9k before the schema work because 25 previously hidden parameters are
  now advertised (with descriptions, enums and ranges) and every tool carries a
  title and annotations. Heaviest: `subtract_blank` (~1.0k, the weighted-method
  parameters), `record_successful_parameters`, `approve_and_learn`,
  `denoise_advanced`, `optimize_parameters`. Still far below e.g. GitHub's MCP
  server (≈55k for 43 tools); no urgent token pressure.
- **0 params missing a description**, **0 tools with weak descriptions**,
  **0 tools missing a title**, **0 tools missing annotations**.
- **0 untargeted tools**; **all 29 tasks reference existing tools**.
- Earlier eval-driven fix, kept for the record: the harness caught
  `analyze_weighted_subtraction` — documented and implemented, but never
  registered. Now wired up and targeted.

### Open items the eval points at (future PRs)

- **Confusable groups** worth either consolidating or giving sharper
  "Use this when …" descriptions, to be confirmed by the agentic pass:
  - background removal: `subtract_blank` / `gaussian_subtract` / `estimate_background`
  - denoising: `denoise` / `denoise_advanced`
  - quality: `assess_quality` / `compute_snr`
  - parameter suggestion: `suggest_parameters` / `suggest_with_learning` / `get_learned_parameters`
  - record learning: `record_successful_parameters` / `approve_and_learn`
- **Resources and the prompt** (`kintsugi://session/images`,
  `kintsugi://session/history/{channel}`, `kintsugi://learning/stats`,
  `kintsugi://project/channels`; `signal_isolation_workflow`) exist but are not yet
  covered by tasks — extend the harness with resource-read and prompt tasks.
- **`subtract_blank` token cost** — the weighted-method parameters make it the
  heaviest definition; consider trimming descriptions once the agentic pass shows
  they are not needed for correct selection.

Run the agentic pass with a key to turn the "confusable groups" hypotheses into
measured tool-selection accuracy, then act on the tools that get mis-selected.

## Recorded session

`python evals/mcp/demo.py` creates a synthetic project, drives the live server
through the in-process client (list, load, stats, thumbnail, suggest, preview,
weighted subtraction, assess, save, approve, recall, session state) and prints a
deterministic Markdown transcript. `--write docs/mcp.md` splices it into the
docs between `<!-- mcp-demo:begin/end -->` markers and `tests/test_mcp_demo.py`
fails when that block is stale, so the documented session is always the real one.
