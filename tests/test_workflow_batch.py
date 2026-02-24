"""Tests for kintsugi workflow batch command and supporting functions."""

from pathlib import Path

import pytest
import yaml


# ---------------------------------------------------------------------------
# Helper: create a minimal project directory structure
# ---------------------------------------------------------------------------
def _make_project(
    tmp_path: Path,
    name: str,
    staged: bool = True,
    config: bool = True,
    registered: bool = False,
    signal_isolated: bool = False,
    stitched_cycles: list[int] | None = None,
    decon_cycles: list[int] | None = None,
    edf_cycles: list[int] | None = None,
) -> Path:
    """Create a minimal project directory with sentinel files."""
    proj = tmp_path / name
    proj.mkdir()

    if config:
        wf = proj / "workflow"
        wf.mkdir(parents=True)
        cfg = {
            "project_dir": str(proj),
            "kintsugi_dir": "/dummy",
            "cycles": [1, 2, 3],
            "channels": [1, 2, 3, 4],
            "n_zplanes": 20,
            "resources": {
                "accounts": [
                    {
                        "name": "testacct",
                        "partition_gpu": "hpg-b200",
                        "partition_cpu": "hpg-default",
                        "gpu_slots": 3,
                        "cpu_slots": 8,
                    }
                ],
                "total_gpu_slots": 3,
                "total_cpu_slots": 8,
                "total_slots": 11,
                "cpu_cpus_per_task": 8,
                "cpu_utilization_cap": 0.85,
            },
        }
        (wf / "config.yaml").write_text(yaml.dump(cfg))

    if staged:
        raw = proj / "data" / "raw"
        raw.mkdir(parents=True)
        (raw / ".staged").touch()

    if registered:
        reg = proj / "data" / "processed" / "registered"
        reg.mkdir(parents=True)
        (reg / ".snakemake_complete").touch()

    if signal_isolated:
        si = proj / "data" / "processed" / "signal_isolated"
        si.mkdir(parents=True)
        (si / ".snakemake_complete").touch()

    for cyc_list, stage_name in [
        (stitched_cycles, "stitched"),
        (decon_cycles, "deconvolved"),
        (edf_cycles, "edf"),
    ]:
        if cyc_list:
            for c in cyc_list:
                d = proj / "data" / "processed" / stage_name / f"cyc{c:02d}"
                d.mkdir(parents=True)
                (d / ".snakemake_complete").touch()

    return proj


# ===========================================================================
# Test _discover_batch_projects
# ===========================================================================
class TestDiscoverBatchProjects:
    def test_basic_discovery(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "proj_A")
        _make_project(tmp_path, "proj_B")
        _make_project(tmp_path, "proj_C", staged=False)  # no .staged

        result = _discover_batch_projects(tmp_path)
        names = [p.name for p in result]
        assert "proj_A" in names
        assert "proj_B" in names
        assert "proj_C" not in names  # no .staged

    def test_excludes_completed(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "proj_done", registered=True, signal_isolated=True)
        _make_project(tmp_path, "proj_todo")

        result = _discover_batch_projects(tmp_path)
        names = [p.name for p in result]
        assert "proj_done" not in names
        assert "proj_todo" in names

    def test_registered_but_not_isolated_is_eligible(self, tmp_path):
        """Projects with registration done but no signal isolation are still eligible."""
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "proj_reg_only", registered=True)

        result = _discover_batch_projects(tmp_path)
        names = [p.name for p in result]
        assert "proj_reg_only" in names

    def test_force_includes_completed(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "proj_done", registered=True, signal_isolated=True)

        result = _discover_batch_projects(tmp_path, force=True)
        names = [p.name for p in result]
        assert "proj_done" in names

    def test_dataset_filter(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "CX_19-001_spleen")
        _make_project(tmp_path, "CX_19-002_lymph")
        _make_project(tmp_path, "CX_19-003_spleen")

        result = _discover_batch_projects(tmp_path, dataset="19-002")
        assert len(result) == 1
        assert result[0].name == "CX_19-002_lymph"

    def test_dataset_filter_case_insensitive(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "CX_19-001_Spleen")

        result = _discover_batch_projects(tmp_path, dataset="spleen")
        assert len(result) == 1

    def test_no_config(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "no_config", config=False)

        result = _discover_batch_projects(tmp_path)
        assert len(result) == 0

    def test_sorted_by_name(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "zzz_proj")
        _make_project(tmp_path, "aaa_proj")
        _make_project(tmp_path, "mmm_proj")

        result = _discover_batch_projects(tmp_path)
        names = [p.name for p in result]
        assert names == sorted(names)

    def test_skips_files(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        _make_project(tmp_path, "real_project")
        (tmp_path / "not_a_dir.txt").write_text("hello")

        result = _discover_batch_projects(tmp_path)
        assert len(result) == 1

    def test_empty_directory(self, tmp_path):
        from kintsugi.cli import _discover_batch_projects

        result = _discover_batch_projects(tmp_path)
        assert result == []


# ===========================================================================
# Test _detect_project_stage
# ===========================================================================
class TestDetectProjectStage:
    def test_signal_isolated(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", registered=True, signal_isolated=True)
        assert "signal_isolated" in _detect_project_stage(proj)

    def test_registered(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", registered=True)
        assert "registered" in _detect_project_stage(proj)

    def test_edf_partial(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", edf_cycles=[1, 2])
        assert "edf" in _detect_project_stage(proj)

    def test_decon_partial(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", decon_cycles=[1])
        assert "deconvolved" in _detect_project_stage(proj)

    def test_stitched_partial(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", stitched_cycles=[1, 2, 3])
        assert "stitched" in _detect_project_stage(proj)

    def test_staged_not_started(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj")
        assert "staged" in _detect_project_stage(proj)

    def test_not_staged(self, tmp_path):
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(tmp_path, "proj", staged=False)
        assert "not staged" in _detect_project_stage(proj)

    def test_stage_priority(self, tmp_path):
        """Signal isolated takes priority over all partial stages."""
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(
            tmp_path, "proj", registered=True, signal_isolated=True,
            edf_cycles=[1], decon_cycles=[1],
        )
        assert "signal_isolated" in _detect_project_stage(proj)

    def test_registered_priority_over_partial(self, tmp_path):
        """Registered takes priority over partial stages."""
        from kintsugi.cli import _detect_project_stage

        proj = _make_project(
            tmp_path, "proj", registered=True, edf_cycles=[1], decon_cycles=[1]
        )
        assert "registered" in _detect_project_stage(proj)


# ===========================================================================
# Test CLI batch command (Click integration)
# ===========================================================================
class TestBatchCLI:
    def test_batch_help(self):
        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "batch", "--help"])
        assert result.exit_code == 0
        assert "Run processing pipeline across multiple datasets" in result.output

    def test_stop_help(self):
        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "stop", "--help"])
        assert result.exit_code == 0
        assert "Stop a running batch process" in result.output

    def test_batch_no_eligible(self, tmp_path):
        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "batch", str(tmp_path), "--dry-run"])
        assert "No eligible projects found" in result.output

    def test_batch_dry_run_shows_table(self, tmp_path, monkeypatch):
        from click.testing import CliRunner

        from kintsugi.cli import main

        # Create eligible projects
        _make_project(tmp_path, "proj_A")
        _make_project(tmp_path, "proj_B")

        # Mock detect_live_multi_account to avoid needing SLURM
        def mock_detect(*args, **kwargs):
            return {
                "total_gpu_slots": 5,
                "total_gpu_avail": 5,
                "total_cpu_slots": 19,
                "total_avail": 24,
                "accounts": [],
            }

        monkeypatch.setattr(
            "kintsugi.hpc.detect_live_multi_account", mock_detect
        )

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "batch", str(tmp_path), "--dry-run"])
        assert result.exit_code == 0
        assert "proj_A" in result.output
        assert "proj_B" in result.output
        assert "Eligible Datasets" in result.output
        assert "Dry run" in result.output

    def test_stop_no_pidfile(self, tmp_path):
        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "stop", str(tmp_path)])
        assert "No batch PID file found" in result.output


# ===========================================================================
# Test GPU guard in workflow_run
# ===========================================================================
class TestGPUGuard:
    def test_no_gpu_slots_error(self, tmp_path, monkeypatch):
        """workflow run should fail if no GPU slots are detected."""
        from click.testing import CliRunner

        from kintsugi.cli import main

        # Create minimal project
        proj = tmp_path / "test_proj"
        proj.mkdir()
        wf = proj / "workflow"
        wf.mkdir(parents=True)
        cfg = {
            "project_dir": str(proj),
            "kintsugi_dir": "/dummy",
            "cycles": [1],
            "channels": [1],
            "n_zplanes": 10,
            "resources": {
                "accounts": [
                    {
                        "name": "testacct",
                        "partition_gpu": "hpg-b200",
                        "partition_cpu": "hpg-default",
                        "gpu_slots": 0,
                        "cpu_slots": 8,
                    }
                ],
                "total_gpu_slots": 0,
                "total_cpu_slots": 8,
            },
        }
        (wf / "config.yaml").write_text(yaml.dump(cfg))
        (wf / "Snakefile").write_text("# dummy")

        # Create SLURM profile so we reach the GPU check
        profiles = wf / "profiles" / "slurm"
        profiles.mkdir(parents=True)
        (profiles / "config.yaml").write_text("# dummy profile")

        # Mock detect_live_multi_account to return 0 GPU slots
        def mock_detect(*args, **kwargs):
            return {
                "total_gpu_slots": 0,
                "total_gpu_avail": 0,
                "total_cpu_slots": 8,
                "total_avail": 8,
                "accounts": [],
            }

        monkeypatch.setattr(
            "kintsugi.hpc.detect_live_multi_account", mock_detect
        )

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "run", str(proj)])
        assert result.exit_code != 0
        assert "No GPU slots detected" in result.output

    def test_no_profile_no_local_errors(self, tmp_path):
        """workflow run without --local should error when no SLURM profile."""
        from click.testing import CliRunner

        from kintsugi.cli import main

        # Create project with config but NO profiles directory
        proj = tmp_path / "test_proj"
        proj.mkdir()
        wf = proj / "workflow"
        wf.mkdir(parents=True)
        cfg = {
            "project_dir": str(proj),
            "kintsugi_dir": "/dummy",
            "cycles": [1],
            "channels": [1],
            "n_zplanes": 10,
            "resources": {},
        }
        (wf / "config.yaml").write_text(yaml.dump(cfg))
        (wf / "Snakefile").write_text("# dummy")

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "run", str(proj)])
        assert result.exit_code != 0
        assert "No SLURM profile found" in result.output
