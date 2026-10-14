---
name: mcp-fastmcp-structured-output
description: "Build or migrate a Python MCP server on the 2.x SDK (MCPServer) with signature-derived, enforced input schemas, TypedDict output schemas, titles/annotations, and an in-process client for tests and evals. Trigger: mcp 2.x migration, 'mcp.server.fastmcp' ImportError, adding outputSchema/structuredContent, 'None is not of type' output-validation errors, building MCP tool-selection evals."
author: KINTSUGI Team
date: 2026-10-14
---

# MCP Server on the 2.x SDK: MCPServer + Structured Output + Evaluation

## Experiment Overview
| Item | Details |
|------|---------|
| **Date** | 2026-06-22 (FastMCP migration), updated 2026-10-14 (mcp 2.x) |
| **Goal** | Run `kintsugi.mcp` on the official 2.x SDK, derive and enforce input schemas from the handlers, deliver typed structured output, and keep an evaluation/docs-drift gate |
| **Environment** | `mcp>=2.0,<3` (verified with 2.3.0), Python 3.10–3.12, `mcp.server.mcpserver.MCPServer` |
| **Status** | Success |

## Context
The first migration (2026-06) moved from the low-level `mcp.server.Server` to
`mcp.server.fastmcp.FastMCP` on the 1.x line. The 2.x SDK that `pip install mcp`
resolves today renamed that class to `MCPServer` and removed `mcp.server.fastmcp`
and the in-memory test helper, so the server failed to import while the project
still pinned `mcp<2`. This skill records what holds on 2.x.

## What Worked (verified approach)

1. **Pin `mcp>=2.0,<3` everywhere** (`pyproject.toml` extra, `deps.py`
   `OPTIONAL_GROUPS` for both `claude` *and* `dev`, env files) and gate the import
   in one module (`kintsugi/mcp/_sdk.py`) that distinguishes "mcp is not
   installed" from "an unsupported mcp 1.x is installed" and exposes
   `MCPServer`, `Context`, `ToolAnnotations`, `ImageContent`, `CallToolResult`.
   Install text must read `pip install -e ".[claude]"` — `pip install
   kintsugi[claude]` installs an unrelated PyPI package.

2. **`MCPServer(name, instructions=..., version=..., website_url=...)`** and
   register each handler with a **title** and **`ToolAnnotations`**
   (readOnly/destructive/idempotent/openWorld) from a declarative, SDK-free
   registry (`tool_specs.py`: name, handler, title, description, annotations).
   Handlers stay plain `async` functions, so unit tests call them directly.

3. **Signature-derived schemas replace the private override.** Put parameter
   descriptions, enums and ranges in the signature —
   `Annotated[int, Field(description=..., ge=1, le=99)]`, `Literal["global",
   "weighted"]` — and the SDK both advertises and **enforces** them: a bad enum or
   out-of-range value comes back as `isError` before the handler runs. The 1.x
   trick of overwriting `server._tool_manager.get_tool(name).parameters` with a
   curated schema is gone, and so is the 884-line JSON-Schema copy that drifted.

4. **`TypedDict` outputs work on 2.x.** Return `TypedDict(total=False)` classes
   (`schemas.py`: `ImageStatsResult`, `QualityResult`, `LearningSuggestion`,
   session state) and the SDK emits them as `outputSchema` + `structuredContent`
   without materialising absent keys as `None` (the 1.x failure below). Plain
   `dict[str, Any]` is still fine for tools that do not need a typed schema.

5. **Test and evaluate through the in-process client.**
   `from mcp.client import Client` replaces the removed
   `mcp.shared.memory.create_connected_server_and_client_session`:
   ```python
   from mcp.client import Client
   async with Client(create_server()) as client:
       tools = (await client.list_tools()).tools      # titles, annotations, schemas on the wire
       result = await client.call_tool("denoise", {"channel": "x", "method": "gaussian"})
       assert result.is_error                           # schema enforced by the SDK
   ```
   `tests/test_mcp_protocol.py` checks registry parity, metadata on the wire,
   schema/signature parity, validation errors as `isError` and structured round
   trips; `tests/test_mcp_stdio.py` drives a real subprocess and asserts stdout is
   pure JSON-RPC (the CLI banner goes to stderr).

6. **Generate the docs from the registry.** `kintsugi mcp tools
   --markdown/--write/--check` (`docs_sync.py`) splices the tool table/list into
   `<!-- mcp-tools:* -->` blocks; `--check` fails on drift in CI.
   `evals/mcp/harness.py --strict` fails on missing descriptions, titles,
   annotations, untargeted tools and task-coverage gaps.

## What Failed (do not repeat)

| Attempt | Result |
|---------|--------|
| Importing `mcp.server.fastmcp.FastMCP` under a `mcp<2` pin | **FAILED on 2.x.** The module no longer exists; the server did not start. Target `mcp.server.mcpserver.MCPServer` and pin `>=2.0,<3` |
| `mcp.shared.memory.create_connected_server_and_client_session` for tests | **Removed in 2.x.** Use `from mcp.client import Client` (`async with Client(server)`) |
| Overriding `server._tool_manager...parameters` to keep curated descriptions/enums | **Obsolete and private.** Put them in `Annotated[..., Field(...)]` / `Literal[...]`; the SDK enforces them |
| Strict `TypedDict` returns on **1.x** FastMCP | Absent keys became `None` → `Output validation error: None is not of type 'string'`, `isError=True`. Fixed on 2.x; if you must support 1.x, return `dict[str, Any]` |
| Rely on docstrings for per-parameter descriptions | **Dropped.** Only `Annotated` / `Field(description=...)` reaches the schema |
| Printing the startup banner to stdout in `kintsugi mcp start` | **Corrupts the stdio transport.** Everything but JSON-RPC goes to stderr |
| Raising the `mcp` floor only in `pyproject.toml` extras | **Incomplete.** `deps.py` `OPTIONAL_GROUPS` is a parallel source of truth for `kintsugi install`/`check` |

## Verification
- `pytest tests/test_mcp_tools.py tests/test_mcp_path_safety.py tests/test_mcp_protocol.py tests/test_mcp_config.py tests/test_mcp_stdio.py tests/test_mcp_docs_sync.py` (`KINTSUGI_REQUIRE_MCP=1` makes MCP skips fail).
- `python evals/mcp/harness.py --strict` and `kintsugi mcp tools --check docs/mcp.md docs/cli.md README.md CLAUDE.md`.
- Preserve the `create_server()` `ImportError` contract by simulating
  `MCP_AVAILABLE = False` in `_sdk.py` (don't swallow the import in a bare `except`).

## Key Files (reference implementation)
- `src/kintsugi/mcp/_sdk.py` — 2.x import gate and version message
- `src/kintsugi/mcp/server.py` — `MCPServer` wiring: tools, resources, prompt, instructions
- `src/kintsugi/mcp/tool_specs.py` — declarative registry (mcp-free, also feeds `kintsugi mcp tools`)
- `src/kintsugi/mcp/schemas.py` — `TypedDict` output schemas
- `src/kintsugi/mcp/docs_sync.py` — generated tool blocks (`--write` / `--check`)
- `evals/mcp/{tasks.py,harness.py,README.md}` — evaluation harness; `.github/workflows/mcp-eval.yml` — keyed agentic eval + structural PR gate
- `docs/mcp.md` — canonical user-facing guide
