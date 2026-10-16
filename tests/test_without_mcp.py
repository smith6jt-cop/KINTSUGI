"""KINTSUGI works without the optional MCP SDK and without any AI assistant configured."""

from __future__ import annotations

import json
import subprocess
import sys
import textwrap

# Runs in a fresh interpreter where ``import mcp`` fails, as on an install without
# the ``mcp`` extra; stdin is not a terminal, so nothing is prompted.
_SCRIPT = textwrap.dedent(
    """
    import importlib.metadata as md
    import json, sys

    # The optional SDK is not installed: no module and no distribution metadata.
    sys.modules["mcp"] = None
    _version = md.version

    def version(name):
        if name == "mcp":
            raise md.PackageNotFoundError(name)
        return _version(name)

    md.version = version
    import kintsugi.project
    kintsugi.project.KintsugiProject.setup_notebooks = lambda self, **kw: None
    from kintsugi.assist import ParameterLearningEngine, WorkflowState  # noqa: F401
    from kintsugi.cli import main

    def run(argv):
        try:
            main(argv, standalone_mode=False)
        except SystemExit as exc:
            return exc.code or 0
        except Exception as exc:  # click errors in non-standalone mode
            return type(exc).__name__
        return 0

    root = sys.argv[1]
    print(json.dumps({
        "init": run(["init", root, "--force", "--name", "no-mcp"]),
        "start": run(["mcp", "start", "--project-dir", root]),
    }))
    """
)


def test_init_works_without_the_mcp_sdk(tmp_path):
    root = tmp_path / "proj"
    proc = subprocess.run(
        [sys.executable, "-c", _SCRIPT, str(root)],
        capture_output=True,
        text=True,
        stdin=subprocess.DEVNULL,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stderr
    codes = json.loads(proc.stdout.strip().splitlines()[-1])
    assert codes["init"] == 0, proc.stdout + proc.stderr
    assert (root / "kintsugi_project.json").exists()
    assert (root / "meta" / "experiment.json").exists()
    for name in (".mcp.json", ".claude", ".cursor", ".vscode/mcp.json", "AGENTS.md"):
        assert not (root / name).exists(), name
    # The MCP server, and only it, needs the SDK: it says how to install it.
    assert codes["start"] not in (0, None)
    output = " ".join((proc.stdout + proc.stderr).split())
    assert "The 'mcp' package is not installed." in output
    assert 'Install with: pip install -e ".[mcp]"' in output
