---
name: mcp-fastmcp-structured-output
description: "Migrate a low-level MCP server to FastMCP with structured output, preserve curated input schemas, and evaluate the tool surface. Trigger: FastMCP migration, adding outputSchema/structuredContent, 'None is not of type' output-validation errors, building MCP tool-selection evals."
author: KINTSUGI Team
date: 2026-06-22
---

# MCP Server: FastMCP Migration + Structured Output + Evaluation

## Experiment Overview
| Item | Details |
|------|---------|
| **Date** | 2026-06-22 |
| **Goal** | Migrate `kintsugi.mcp` from low-level `mcp.server.Server` to FastMCP, deliver structured tool output, and stand up an evaluation/improvement workflow |
| **Environment** | `mcp==1.27.2` (floor `>=1.10`), Python 3.11/3.12, FastMCP (`mcp.server.fastmcp`) |
| **Status** | Success |

## Context
The MCP server hand-wrote a `Tool(...)` list + an `elif` dispatch in `call_tool`,
returned only `TextContent` (a JSON string), and exposed no output schemas. Goal:
move to FastMCP so each tool returns machine-readable `structuredContent`, on a
foundation that makes later resources/prompts/elicitation cheap.

## What Worked (verified approach)

1. **Floor `mcp>=1.10.0,<2`.** Structured output (`Tool.outputSchema`,
   `CallToolResult.structuredContent`, FastMCP `structured_output`) and
   elicitation **do not exist below 1.10** (the SDK release after spec 2025-06-18).

2. **Register existing handlers; keep them FastMCP-agnostic.** A declarative
   registry (`tool_specs.py`: name → module/attr + description + input schema)
   plus, in `server.py`:
   ```python
   server = FastMCP("kintsugi")
   handler = getattr(module, spec["attr"])           # plain async dict-returning fn
   server.add_tool(handler, name=spec["name"],
                   description=spec["description"], structured_output=True)
   ```
   Handlers stay plain `async` functions returning `dict[str, Any]`, so the unit
   tests keep calling them directly.

3. **Returns: `dict[str, Any]` + `structured_output=True`.** FastMCP then emits
   the dict as `structuredContent` **and** a JSON text block, for both the
   success and the in-band `{"error": ...}` branches — no per-tool models needed.

4. **Preserve curated input schemas (descriptions + enums) by overriding** the
   advertised schema after registration; validation still runs from the handler
   signature:
   ```python
   # FastMCP's signature-derived schema drops per-param descriptions + enums.
   server._tool_manager.get_tool(name).parameters = curated_input_schema  # wrap in try/except
   ```

5. **Evaluation harness** (`evals/mcp/`): a structural pass (tool-definition token
   cost, schema-quality checks, confusable-group flags, task coverage) + an
   agentic pass (`tool_choice="auto"`, score the model's first `tool_use`).
   Tool **selection** needs only the schemas — **no image data or handler
   execution** — so it runs with just an API key. Wire agentic into CI behind an
   `ANTHROPIC_API_KEY` secret (manual + scheduled); gate PRs with the structural
   `--strict` pass. The harness immediately found `analyze_weighted_subtraction`
   documented + implemented but **never registered**.

## What Failed (do not repeat)

| Attempt | Result |
|---------|--------|
| Strict `TypedDict` returns (`total=False`) to get rich per-field outputSchema | **FAILED.** FastMCP materializes absent keys as `None`, then validates against the non-nullable field type → `Output validation error: None is not of type 'string'`, `structuredContent=None`, `isError=True` on **both** success and error returns |
| `Optional`-field `TypedDict`s (to dodge the above) | Validates, but **pollutes every result** with `null` keys for absent fields (noisy/misleading) → use plain `dict[str, Any]` instead |
| Rely on docstrings for per-parameter descriptions under FastMCP | **FAILED.** FastMCP 1.27 does **not** parse docstring args; descriptions are dropped. Use `Annotated[T, Field(description=...)]` or override `.parameters` |
| Raise the `mcp` floor only in `pyproject.toml` extras | **Incomplete.** `src/kintsugi/deps.py` `OPTIONAL_GROUPS` is a parallel source of truth used by `kintsugi install`/`check` — update both (`claude` *and* `dev`) |

## Verification
- `pytest tests/test_mcp_tools.py tests/test_mcp_path_safety.py` — direct-call
  tests stay green (handlers unchanged); add an in-memory client test
  (`mcp.shared.memory.create_connected_server_and_client_session`) asserting
  `outputSchema` present + `structuredContent` populated + text block + error
  branch clean.
- `python evals/mcp/harness.py --strict` — structural gate (fails on coverage
  gaps / missing descriptions).
- Preserve the `create_server()` `ImportError` contract by simulating
  `MCP_AVAILABLE = False` (don't swallow the import in a bare `except`).

## Key Files (reference implementation)
- `src/kintsugi/mcp/server.py` — FastMCP wiring, `_advertise_schema` override
- `src/kintsugi/mcp/tool_specs.py` — declarative registry (mcp-free, also feeds `kintsugi mcp tools`)
- `evals/mcp/{tasks.py,harness.py,README.md}` — evaluation harness
- `.github/workflows/mcp-eval.yml` — keyed agentic eval + structural PR gate
