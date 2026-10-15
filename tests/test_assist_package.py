"""The ``kintsugi.assist`` package: lazy exports and the renamed install group."""

from __future__ import annotations

import subprocess
import sys

import pytest


def test_exports_resolve_lazily():
    """Importing the learning engine must not import the Jupyter tuner module."""
    code = (
        "import sys\n"
        "from kintsugi.assist import ParameterLearningEngine\n"
        "assert ParameterLearningEngine.__module__ == 'kintsugi.assist.parameter_learning'\n"
        "assert 'kintsugi.assist.param_tuner' not in sys.modules, 'tuners imported eagerly'\n"
    )
    subprocess.run([sys.executable, "-c", code], check=True)


def test_every_export_and_submodule_imports():
    import kintsugi.assist as assist

    for name in assist.__all__:
        assert getattr(assist, name) is not None
        assert name in dir(assist)
    from kintsugi.assist import param_tuner, workflow_state  # submodule fallback

    assert param_tuner.blank_subtraction_tuner is assist.blank_subtraction_tuner
    assert workflow_state.WorkflowState is assist.WorkflowState


def test_unknown_name_raises_attribute_error():
    import kintsugi.assist as assist

    with pytest.raises(AttributeError):
        assist.does_not_exist  # noqa: B018


def test_old_package_name_is_gone():
    with pytest.raises(ModuleNotFoundError):
        __import__("kintsugi.claude")


def test_install_claude_points_to_mcp():
    from click.testing import CliRunner

    from kintsugi.cli import main

    result = CliRunner().invoke(main, ["install", "claude"])
    assert result.exit_code == 1
    assert "renamed to 'mcp'" in result.output
    assert "kintsugi install mcp" in result.output
