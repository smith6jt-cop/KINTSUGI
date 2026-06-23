# KINTSUGI MCP — Evaluation Harness

Measures the quality of the KINTSUGI MCP **tool surface** with realistic tasks
instead of intuition, following the `mcp-builder` philosophy. Use it to drive and
verify MCP improvements (tool consolidation, description quality, resources/prompts).

## Why

The server exposes 27 tools to Claude. Two questions decide whether that surface
is *good*:

1. **Is it cheap and clear?** (structural) — how many tokens do the tool
   definitions cost, and are any tools/params under-described or confusable?
2. **Does the model use it correctly?** (agentic) — given a realistic request,
   does Claude pick the right tool?

## Run it

```bash
# Structural only — no network, no API key. Runs anywhere mcp is installed.
python evals/mcp/harness.py

# Agentic tool-selection eval — needs an API key + the anthropic SDK.
pip install anthropic
ANTHROPIC_API_KEY=sk-... python evals/mcp/harness.py --agentic --model claude-sonnet-4-6

# Save the full report (structural + agentic) as JSON.
python evals/mcp/harness.py --agentic --json mcp_eval_report.json
```

The agentic mode measures **tool selection only** — it records the model's first
`tool_use` per task and never executes the handler — so it needs **no image data
or project**, just the live tool schemas.

## Files

| File | Purpose |
|------|---------|
| `tasks.py` | 15 realistic natural-language tasks with `expected` acceptable tool(s). Includes deliberate confusable pairs (`denoise` vs `denoise_advanced`; the three background-removal tools; `assess_quality` vs `compute_snr`). |
| `harness.py` | Builds the live FastMCP server, runs the structural metrics, and (opt-in) the agentic tool-selection loop. |

Add a task by appending an `EvalTask` to `TASKS`; the structural pass will verify
its `expected` tools exist on the server.

## Findings (baseline, structural pass)

Run against the post-migration server:

- **27 tools**, **27/27 carry an output schema** (structured output is live).
- **~3.9k tokens** of tool definitions injected on every call (~144/tool) — modest;
  no urgent token pressure (compare: GitHub's MCP server ≈ 55k for 43 tools).
- **0 params missing a description**, **0 tools with weak descriptions** — the
  curated input schemas held up.
- **First eval-driven fix:** the harness caught `analyze_weighted_subtraction` —
  documented in `CLAUDE.md` and implemented, but **never registered**. Now wired up.

### Open items the eval points at (future PRs)

- **Confusable groups** worth either consolidating or giving sharper
  "Use this when …" descriptions, to be confirmed by the agentic pass:
  - background removal: `subtract_blank` / `gaussian_subtract` / `estimate_background`
  - denoising: `denoise` / `denoise_advanced`
  - parameter suggestion: `suggest_parameters` / `suggest_with_learning` / `get_learned_parameters`
- **Resources/prompts** (deferred P1 items) — once added, extend the harness with
  resource-read and prompt tasks.

Run the agentic pass with a key to turn the "confusable groups" hypotheses into
measured tool-selection accuracy, then act on the tools that get mis-selected.
