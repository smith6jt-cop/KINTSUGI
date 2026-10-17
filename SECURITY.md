# Security Policy

## Reporting a vulnerability

Please report suspected vulnerabilities privately, not in a public issue:

- GitHub private vulnerability reporting:
  <https://github.com/smith6jt-cop/KINTSUGI/security/advisories/new>
- or by email to the maintainer listed in `pyproject.toml`.

You will get an acknowledgement within a week. Fixes ship in the next patch
release of the current minor version; older minors are not maintained.

## Scope and trust model of the MCP server

`kintsugi mcp start` is a **local, single-user stdio server**: an MCP client
spawns it as a subprocess, talks JSON-RPC over stdin/stdout, and it has no
network listener, no authentication, and no multi-tenant state. It runs with the
permissions of the user who started the client. Its project comes from
`--project-dir` or `KINTSUGI_PROJECT_DIR` (or a working directory that is a
KINTSUGI project); without one, tools that need a project refuse rather than
write under an arbitrary directory. Everything the
model sends as a tool argument is treated as untrusted input:

- **Path traversal** - filenames must be bare names (`safe_filename`), and
  every resolved path must stay inside the project directory
  (`ensure_within`, which follows symlinks so a symlink escape is rejected).
- **Model checkpoints** - `denoise_advanced` and `predict_parameters` only load
  files under a KINTSUGI `project_path`, the server's project, or
  `~/.kintsugi/models` (`allowed_model_roots`); PyTorch checkpoints are loaded
  with `weights_only=True`.
- **Resource limits** - numeric arguments carry schema bounds that the SDK
  enforces, plus stricter in-body clamps (`clamp_int`/`clamp_float` reject
  booleans, NaN and infinities).
- **Approval scope** - the generated `.claude/settings.local.json` approves only
  the `kintsugi` server (`enabledMcpjsonServers`), never every project server.

The security helpers live in `src/kintsugi/mcp/path_safety.py` and are covered
by `tests/test_mcp_path_safety.py`. See `docs/mcp.md` for the full design.
