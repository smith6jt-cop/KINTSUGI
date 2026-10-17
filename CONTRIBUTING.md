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
`kintsugi mcp tools --check docs/mcp.md docs/cli.md README.md CLAUDE.md AGENTS.md` green.
`AGENTS.md` is the cross-tool entry point for AI assistants; `CLAUDE.md` holds the
same project knowledge for Claude Code.

The repository's own `.mcp.json` starts `kintsugi mcp start` without `--project-dir`:
the repository root is not a KINTSUGI project, so that server runs project-less and
tools that read or write project files need an explicit `project_path` (for example a
project created with `kintsugi init` under the gitignored `test_data/`).
