"""
Tests for project initialization and template refresh behavior.
"""

import json

from click.testing import CliRunner

from kintsugi.project import KintsugiProject, ProjectConfig, init_project


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


# =========================================================================
# Issue 2: init_project() name update on existing projects
# =========================================================================


def test_init_project_updates_name_on_existing(tmp_path):
    """init_project() with a new name should update and persist it."""
    project_root = tmp_path / "proj"
    # Create initial project
    project = KintsugiProject.create(project_root, name="OldName")
    project.save()
    assert project.config.name == "OldName"

    # Reload with a different name
    reloaded = init_project(project_root, name="NewName", interactive=False, check_processed=False)
    assert reloaded.config.name == "NewName"

    # Verify it was persisted to disk
    reloaded2 = KintsugiProject.load(project_root)
    assert reloaded2.config.name == "NewName"


def test_init_project_preserves_name_when_none(tmp_path):
    """init_project() without a name should keep the existing name."""
    project_root = tmp_path / "proj"
    project = KintsugiProject.create(project_root, name="KeepMe")
    project.save()

    reloaded = init_project(project_root, name=None, interactive=False, check_processed=False)
    assert reloaded.config.name == "KeepMe"


def test_init_project_no_update_when_same_name(tmp_path):
    """init_project() with the same name should not trigger a save."""
    project_root = tmp_path / "proj"
    project = KintsugiProject.create(project_root, name="SameName")
    project.save()

    reloaded = init_project(project_root, name="SameName", interactive=False, check_processed=False)
    assert reloaded.config.name == "SameName"


# =========================================================================
# Issue 3: --force updates default_parameters.json, name, experiment.json
# =========================================================================


def test_force_updates_project_name(tmp_path, monkeypatch):
    """kintsugi init --force --name should update project name."""
    repo_path = tmp_path / "repo"
    _create_fake_repo(repo_path)
    monkeypatch.setattr(KintsugiProject, "_detect_kintsugi_path", staticmethod(lambda: repo_path))

    project_root = tmp_path / "project"
    runner = CliRunner()
    from kintsugi.cli import main

    # Create project with initial name
    result = runner.invoke(main, ["init", str(project_root), "--force", "--name", "Original"])
    assert result.exit_code == 0

    # Force update with new name
    result = runner.invoke(main, ["init", str(project_root), "--force", "--name", "Updated"])
    assert result.exit_code == 0
    assert "Updated project name" in result.output

    # Verify persisted
    project = KintsugiProject.load(project_root)
    assert project.config.name == "Updated"


def test_force_regenerates_default_parameters(tmp_path, monkeypatch):
    """kintsugi init --force should regenerate default_parameters.json."""
    repo_path = tmp_path / "repo"
    _create_fake_repo(repo_path)
    monkeypatch.setattr(KintsugiProject, "_detect_kintsugi_path", staticmethod(lambda: repo_path))

    project_root = tmp_path / "project"
    runner = CliRunner()
    from kintsugi.cli import main

    # Create project
    result = runner.invoke(main, ["init", str(project_root), "--force"])
    assert result.exit_code == 0

    params_file = project_root / "configs" / "default_parameters.json"
    assert params_file.exists()

    # Corrupt the file
    params_file.write_text("{}")

    # Force should regenerate it
    result = runner.invoke(main, ["init", str(project_root), "--force"])
    assert result.exit_code == 0
    assert "Regenerated default_parameters.json" in result.output

    with open(params_file) as f:
        params = json.load(f)
    assert "illumination_correction" in params
    assert "stitching" in params


def test_force_updates_experiment_json_selectively(tmp_path, monkeypatch):
    """kintsugi init --force with microscope params should update only those fields."""
    repo_path = tmp_path / "repo"
    _create_fake_repo(repo_path)
    monkeypatch.setattr(KintsugiProject, "_detect_kintsugi_path", staticmethod(lambda: repo_path))

    project_root = tmp_path / "project"
    runner = CliRunner()
    from kintsugi.cli import main

    # Create project with defaults
    result = runner.invoke(main, ["init", str(project_root), "--force"])
    assert result.exit_code == 0

    exp_file = project_root / "meta" / "experiment.json"
    assert exp_file.exists()

    with open(exp_file) as f:
        original = json.load(f)
    original_na = original["numerical_aperture"]

    # Force update with only xy_pixel_size changed
    result = runner.invoke(main, ["init", str(project_root), "--force", "--xy-pixel-size", "500.0"])
    assert result.exit_code == 0
    assert "Updated experiment.json" in result.output

    with open(exp_file) as f:
        updated = json.load(f)
    assert updated["xy_pixel_size"] == 500.0
    # Other fields should be preserved
    assert updated["numerical_aperture"] == original_na


def test_create_default_params_file_standalone(tmp_path):
    """_create_default_params_file() can be called independently."""
    project_root = tmp_path / "proj"
    project = KintsugiProject.create(project_root, name="test")

    # Verify it was created during create()
    params_file = project.paths.configs / "default_parameters.json"
    assert params_file.exists()

    # Delete and recreate
    params_file.unlink()
    assert not params_file.exists()

    result_path = project._create_default_params_file()
    assert result_path.exists()
    with open(result_path) as f:
        params = json.load(f)
    assert "illumination_correction" in params
    assert "registration" in params
