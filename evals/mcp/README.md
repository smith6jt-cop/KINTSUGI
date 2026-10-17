# KINTSUGI MCP — Evaluation Harness

Measures the quality of the KINTSUGI MCP **tool surface** with realistic tasks
instead of intuition, following the `mcp-builder` philosophy. Use it to drive and
verify MCP improvements (tool consolidation, description quality, resources/prompts).

## Why

The server exposes 34 tools to whichever AI assistant connects to it. Two questions decide whether that surface
is *good*:

1. **Is it cheap and clear?** (structural) — how many tokens do the tool
   definitions cost, and are any tools/params under-described, untitled,
   unannotated or confusable?
2. **Does the model use it correctly?** (agentic) — given a realistic request,
   does the model pick the right tool?

## Run it

```bash
# Structural only — no network, no API key. Needs mcp>=2.0 (pip install -e ".[mcp]").
python evals/mcp/harness.py

# Same, but exit non-zero on any structural problem (the CI gate).
python evals/mcp/harness.py --strict

# Agentic tool-selection eval — a model provider and a model (no SDK needed).
ANTHROPIC_API_KEY=... python evals/mcp/harness.py --agentic --provider anthropic --model MODEL
OPENAI_API_KEY=... python evals/mcp/harness.py --agentic --provider openai --model MODEL

# Any OpenAI-compatible endpoint, for example a local model served by Ollama.
OPENAI_BASE_URL=http://localhost:11434/v1 \
  python evals/mcp/harness.py --agentic --provider openai --model MODEL

# Save the full report (structural + agentic) as JSON.
python evals/mcp/harness.py --agentic --provider NAME --model MODEL --json mcp_eval_report.json
```

The harness reads the tools **through the in-process MCP client**, so it sees
exactly what a connected client sees (titles, annotations, the signature-derived
input schemas, output schemas). The agentic mode measures **tool selection only** —
it records the model's first tool call per task and never executes the handler — so
it needs **no image data or project**, just the live tool schemas.

### Providers

`providers.py` speaks two wire formats with the standard library (one JSON POST
per task, no vendor SDK). There is no default provider or model, and none is
inferred from the API keys that happen to be set: pass `--provider` and `--model`,
or answer the prompts the harness shows in a terminal. Without a terminal, a missing
choice is a usage error (exit 2).

| `--provider` | API | Key | Base URL (`--base-url` overrides) |
|---|---|---|---|
| `anthropic` | Messages (`POST /v1/messages`) | `ANTHROPIC_API_KEY` | `ANTHROPIC_BASE_URL`, else `https://api.anthropic.com` |
| `openai` | Chat Completions (`POST /chat/completions`) | `OPENAI_API_KEY` (optional for a custom base URL) | `OPENAI_BASE_URL`, else `https://api.openai.com/v1` |

`openai` covers every service and server that implements Chat Completions: point
`OPENAI_BASE_URL` at it (Ollama `http://localhost:11434/v1`, a vLLM server, or a
hosted platform's OpenAI-compatible endpoint, with that platform's key in
`OPENAI_API_KEY`). For this format the input schemas are sent without annotation-only
keywords (`title`, `default`, `format`) and with `anyOf: [T, null]` reduced to `T`;
constraints (`enum`, bounds, `required`, `additionalProperties`) are kept. Rate-limit and overload responses are retried
twice; any other error stops the run with the HTTP status and the API's message
(never the key). The report's `agentic` block records `provider`, `model` and the
endpoint `host`.

The CI job (`.github/workflows/mcp-eval.yml`) runs the agentic pass only with an
explicit choice: on manual dispatch, the `provider` (`none`, `anthropic`, `openai`),
`model` and optional `base_url` inputs; on the weekly schedule, the repository
variables `EVAL_PROVIDER`, `EVAL_MODEL` and optional `EVAL_BASE_URL`. The key comes
from the secret `ANTHROPIC_API_KEY` or `OPENAI_API_KEY`. With neither a provider nor
a model the pass is skipped with a notice; with only one of them the job fails. The
run writes `mcp_eval_report_<provider>.json`.

## Files

| File | Purpose |
|------|---------|
| `tasks.py` | 31 realistic natural-language tasks with `expected` acceptable tool(s); every registered tool is the target of at least one task. Includes deliberate confusable pairs (`denoise` vs `denoise_advanced`; the three background-removal tools; `assess_quality` vs `compute_snr`; the parameter-suggestion and record-learning groups). |
| `harness.py` | Builds the live server, runs the structural metrics, and (opt-in) the agentic tool-selection loop. |
| `providers.py` | Standard-library adapters for the agentic loop: Anthropic Messages and OpenAI-compatible Chat Completions. |

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
`docs/cli.md`, `README.md` and `AGENTS.md`, and `kintsugi mcp tools --check FILE...`
exits non-zero when a block is stale. CI runs the check, so adding, renaming or
re-describing a tool without regenerating the docs fails the build
(`tests/test_mcp_docs_sync.py` covers the generator itself).

## Findings (baseline, structural pass)

`python evals/mcp/harness.py --json baseline.json` on the current server:

- **34 tools**, **34/34 carry an output schema** (structured output on every tool;
  `get_image_stats`, `assess_quality`, `suggest_with_learning`,
  `get_session_state` and `get_workflow_status` declare typed schemas).
- **~8.6k tokens** of tool definitions injected on every call (~254/tool). This is
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
  - quality: `assess_quality` / `compute_snr` / `batch_assess_quality`
  - parameter suggestion: `suggest_parameters` / `suggest_with_learning` / `get_learned_parameters`
  - record learning: `record_successful_parameters` / `approve_and_learn`
- **Resources and the prompt** (`kintsugi://session/images`,
  `kintsugi://session/history/{channel}`, `kintsugi://learning/stats`,
  `kintsugi://project/channels`; `signal_isolation_workflow`) exist but are not yet
  covered by tasks — extend the harness with resource-read and prompt tasks.
- **`subtract_blank` token cost** — the weighted-method parameters make it the
  heaviest definition; consider trimming descriptions once the agentic pass shows
  they are not needed for correct selection.

Run the agentic pass with any provider to turn the "confusable groups" hypotheses into
measured tool-selection accuracy, then act on the tools that get mis-selected.

## Recorded session

`python evals/mcp/demo.py` creates a synthetic project, drives the live server
through the in-process client (list, load, stats, thumbnail, suggest, preview,
weighted subtraction, assess, save, approve, recall, session state) and prints a
deterministic Markdown transcript. `--write docs/mcp.md` splices it into the
docs between `<!-- mcp-demo:begin/end -->` markers and `tests/test_mcp_demo.py`
fails when that block is stale, so the documented session is always the real one.
