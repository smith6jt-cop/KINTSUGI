"""
Tests for project initialization and template refresh behavior.
"""

from click.testing import CliRunner

from kintsugi.project import KintsugiProject, ProjectConfig


def _create_fake_repo(repo_path):
    """Create a minimal fake repo structure with notebooks for testing."""
    notebooks_dir = repo_path / "notebooks"
    notebooks_dir.mkdir(parents=True)
    (repo_path / "pyproject.toml").write_text("[tool.poetry]\nname = 'kintsugi'\n")

    # Core workflow notebook
    (notebooks_dir / "2_Cycle_Processing.ipynb").write_text("v1")
    # Supporting module and file
    support_dir = notebooks_dir / "Kreg"
    support_dir.mkdir()
    (support_dir / "module.py").write_text("v1")
    (notebooks_dir / "Kutils.py").write_text("v1")

    return notebooks_dir


def test_setup_notebooks_overwrite_updates_files(tmp_path):
    """Ensure overwrite=True refreshes notebooks and support files."""
    repo_path = tmp_path / "repo"
    notebooks_dir = _create_fake_repo(repo_path)

    project_root = tmp_path / "project"
    project = KintsugiProject(
        project_root, config=ProjectConfig(name="test"), kintsugi_path=repo_path
    )
    project.paths.create_all()

    project.setup_notebooks(overwrite=True)
    dest_nb = project.paths.notebooks / "2_Cycle_Processing.ipynb"
    dest_support = project.paths.notebooks / "Kreg" / "module.py"
    dest_file = project.paths.notebooks / "Kutils.py"

    assert dest_nb.read_text() == "v1"
    assert dest_support.read_text() == "v1"
    assert dest_file.read_text() == "v1"

    # Update source templates
    (notebooks_dir / "2_Cycle_Processing.ipynb").write_text("v2")
    (notebooks_dir / "Kreg" / "module.py").write_text("v2")
    (notebooks_dir / "Kutils.py").write_text("v2")

    # Default run should keep existing copies
    project.setup_notebooks()
    assert dest_nb.read_text() == "v1"

    # Overwrite should refresh everything
    project.setup_notebooks(overwrite=True)
    assert dest_nb.read_text() == "v2"
    assert dest_support.read_text() == "v2"
    assert dest_file.read_text() == "v2"


def test_init_force_updates_existing_project(tmp_path, monkeypatch):
    """kintsugi init --force should refresh notebooks for existing projects."""
    repo_path = tmp_path / "repo"
    notebooks_dir = _create_fake_repo(repo_path)

    # Ensure CLI uses our fake repo path
    monkeypatch.setattr(
        KintsugiProject,
        "_detect_kintsugi_path",
        staticmethod(lambda: repo_path),
    )

    project_root = tmp_path / "project"
    runner = CliRunner()

    from kintsugi.cli import main

    # Initial creation
    result = runner.invoke(main, ["init", str(project_root), "--force"])
    assert result.exit_code == 0
    dest_nb = project_root / "notebooks" / "2_Cycle_Processing.ipynb"
    assert dest_nb.read_text() == "v1"

    # Update source and rerun with --force to refresh
    (notebooks_dir / "2_Cycle_Processing.ipynb").write_text("v2")
    result = runner.invoke(main, ["init", str(project_root), "--force"])
    assert result.exit_code == 0
    assert dest_nb.read_text() == "v2"
