# Contributing to KINTSUGI

The full contributor guide lives in the documentation:
[docs/contributing.md](docs/contributing.md)
(rendered at <https://kintsugi.readthedocs.io/en/latest/contributing.html>).

Quick start:

```bash
pip install -e ".[dev,mcp]"
pre-commit install && pre-commit install --hook-type commit-msg
pytest tests/ -q
ruff check src/ tests/ evals/ && ruff format --check src/ tests/ evals/
```

Commits follow [Conventional Commits](https://www.conventionalcommits.org/)
(`feat:`, `fix:`, `docs:`, ...). Changes to the MCP server must keep
`python evals/mcp/harness.py --strict` and
`kintsugi mcp tools --check docs/mcp.md docs/cli.md README.md AGENTS.md` green.

Guidance for AI coding assistants lives in `AGENTS.md` files (the root one and one in
`src/kintsugi/`, `src/kintsugi/signal/`, `workflow/` and `notebooks/`); each `CLAUDE.md`
next to them is the single line `@AGENTS.md`, so every assistant reads the same text.
Edit the `AGENTS.md` files, never the stubs.

The repository ships no MCP client configuration (`.mcp.json`, `.claude/`, `.cursor/`
and the like are gitignored), so files you keep for your own assistant stay local. To
try the MCP server, create a project under the gitignored `test_data/` and configure
your client for it:

```bash
kintsugi init test_data/demo --name demo --mcp-client NAME   # 'kintsugi mcp clients' lists them
```
