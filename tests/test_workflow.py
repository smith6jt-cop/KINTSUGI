"""
Tests for the Snakemake workflow implementation.

Covers:
  - CLI commands: ``kintsugi workflow config`` and ``kintsugi workflow run``
  - Snakefile parsing and config validation
  - Workflow script helper functions (quality gate, cycle-dir detection, etc.)
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import yaml
from click.testing import CliRunner

from kintsugi.cli import main

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def project_dir(tmp_path):
    """Create a minimal KINTSUGI project with meta/experiment.json."""
    proj = tmp_path / "project"
    (proj / "meta").mkdir(parents=True)
    (proj / "data" / "raw" / "cyc01").mkdir(parents=True)
    (proj / "data" / "raw" / "cyc02").mkdir(parents=True)
    (proj / "data" / "processed").mkdir(parents=True)
    (proj / "slurm").mkdir(parents=True)

    experiment = {
        "tile_rows": 3,
        "tile_cols": 3,
        "tile_overlap": 0.1,
        "xy_pixel_size": 377.0,
        "z_step_size": 1500.0,
        "numerical_aperture": 0.75,
        "tissue_refractive_index": 1.44,
        "wavelengths": {
            "1": [358.0, 461.0],
            "2": [753.0, 775.0],
            "3": [560.0, 575.0],
            "4": [648.0, 668.0],
        },
        "channels_per_cycle": 4,
        "n_cycles": 2,
        "n_zplanes": 5,
    }
    (proj / "meta" / "experiment.json").write_text(json.dumps(experiment, indent=2))
    return proj


@pytest.fixture
def project_with_slurm_config(project_dir):
    """Project that also has slurm/config.sh with account/partition."""
    config_sh = textwrap.dedent("""\
        #!/usr/bin/env bash
        export ACCOUNT="my_acct"
        export PARTITION="gpu-a100"
    """)
    (project_dir / "slurm" / "config.sh").write_text(config_sh)
    return project_dir


@pytest.fixture
def project_with_workflow(project_dir):
    """Project that already has a workflow directory configured."""
    wf = project_dir / "workflow"
    wf.mkdir(parents=True, exist_ok=True)

    # Minimal Snakefile
    (wf / "Snakefile").write_text('configfile: "config.yaml"\nrule all:\n    input: []\n')

    # Minimal config
    config = {
        "project_dir": str(project_dir),
        "kintsugi_dir": "/fake/kintsugi",
        "cycles": [1, 2],
        "channels": [1, 2, 3, 4],
    }
    (wf / "config.yaml").write_text(yaml.dump(config))

    # Profile
    profile_dir = wf / "profiles" / "slurm"
    profile_dir.mkdir(parents=True)
    (profile_dir / "config.yaml").write_text("executor: slurm\njobs: 8\n")

    return project_dir


# ============================================================================
# _resolve_project_dir helper
# ============================================================================


class TestResolveProjectDir:
    """Test the _resolve_project_dir helper function."""

    def test_returns_same_path_for_project_root(self, project_dir):
        """Normal project dir is returned as-is (resolved)."""
        from kintsugi.cli import _resolve_project_dir

        result = _resolve_project_dir(str(project_dir))
        assert result == project_dir.resolve()

    def test_detects_workflow_subdir(self, project_with_workflow):
        """When given a workflow/ dir, returns the parent project dir."""
        from kintsugi.cli import _resolve_project_dir

        wf_dir = project_with_workflow / "workflow"
        result = _resolve_project_dir(str(wf_dir))
        assert result == project_with_workflow.resolve()

    def test_requires_config_yaml(self, project_dir):
        """A dir named 'workflow' without config.yaml is not auto-detected."""
        from kintsugi.cli import _resolve_project_dir

        wf_dir = project_dir / "workflow"
        wf_dir.mkdir(exist_ok=True)
        # No config.yaml inside
        result = _resolve_project_dir(str(wf_dir))
        assert result == wf_dir.resolve()

    def test_requires_parent_meta_dir(self, tmp_path):
        """A workflow/ dir without meta/ in parent is not auto-detected."""
        from kintsugi.cli import _resolve_project_dir

        # Create workflow dir with config.yaml but no meta/ in parent
        wf_dir = tmp_path / "fake_project" / "workflow"
        wf_dir.mkdir(parents=True)
        (wf_dir / "config.yaml").write_text("cycles: [1]")
        # No meta/ dir in parent
        result = _resolve_project_dir(str(wf_dir))
        assert result == wf_dir.resolve()

    def test_non_workflow_dir_passthrough(self, tmp_path):
        """Arbitrary directory named 'something' is returned unchanged."""
        from kintsugi.cli import _resolve_project_dir

        some_dir = tmp_path / "something"
        some_dir.mkdir()
        result = _resolve_project_dir(str(some_dir))
        assert result == some_dir.resolve()


# ============================================================================
# CLI: install command
# ============================================================================


class TestInstallCommand:
    """Test ``kintsugi install`` command."""

    def test_install_list_shows_all_groups(self, runner):
        """--list should show all groups from deps.py (excluding 'full')."""
        from kintsugi.deps import OPTIONAL_GROUPS

        result = runner.invoke(main, ["install", "--list"])
        assert result.exit_code == 0
        for name in OPTIONAL_GROUPS:
            if name == "full":
                assert name not in result.output
            else:
                assert name in result.output

    def test_install_unknown_group_fails(self, runner):
        """Unknown group should produce error."""
        result = runner.invoke(main, ["install", "nonexistent_group_xyz"])
        assert result.exit_code != 0
        assert "Unknown group" in result.output

    def test_install_no_arg_shows_list(self, runner):
        """No argument should show the group listing."""
        result = runner.invoke(main, ["install"])
        assert result.exit_code == 0
        assert "gpu" in result.output
        assert "analysis" in result.output

    @patch("subprocess.run")
    def test_install_conda_flag(self, mock_run, runner):
        """--conda flag should use conda_cmd when available."""

        mock_run.return_value = MagicMock(returncode=0)
        result = runner.invoke(main, ["install", "gpu", "--conda"])
        assert result.exit_code == 0
        # Should have called the conda command
        call_args = mock_run.call_args
        assert "conda" in call_args[0][0] or "conda" in call_args.kwargs.get("args", [""])[0]

    @patch("subprocess.run")
    @patch("kintsugi.cli._find_project_root")
    def test_install_all_uses_extras(self, mock_root, mock_run, runner, tmp_path):
        """'all' with project root should use single pip extras invocation."""
        mock_root.return_value = tmp_path
        mock_run.return_value = MagicMock(returncode=0)
        result = runner.invoke(main, ["install", "all"])
        assert result.exit_code == 0
        # Should NOT have printed "Installing full:" or sequential per-group lines
        assert "Installing full:" not in result.output
        # Single pip call with extras string
        assert mock_run.call_count == 1
        cmd_arg = (
            mock_run.call_args[0][0]
            if mock_run.call_args[0]
            else mock_run.call_args.kwargs.get("args", "")
        )
        assert "pip install -e" in cmd_arg
        assert "gpu" in cmd_arg
        assert "rapids" not in cmd_arg
        assert "full" not in cmd_arg
        # rapids exclusion note
        assert "rapids" in result.output

    @patch("subprocess.run")
    @patch("kintsugi.cli._find_project_root")
    def test_install_all_fallback_sequential(self, mock_root, mock_run, runner):
        """'all' without project root should fall back to sequential + repair."""
        mock_root.return_value = None
        mock_run.return_value = MagicMock(returncode=0)
        result = runner.invoke(main, ["install", "all"])
        assert result.exit_code == 0
        assert "pyproject.toml not found" in result.output
        # Should have multiple subprocess calls (one per group + repair)
        assert mock_run.call_count > 1
        # Repair step runs 'pip install -e .'
        repair_call = mock_run.call_args_list[-1]
        assert "pip install -e ." in repair_call[0][0]

    @patch("subprocess.run")
    @patch("kintsugi.cli._find_project_root")
    def test_install_all_skips_full_and_rapids(self, mock_root, mock_run, runner, tmp_path):
        """'all' should skip both 'full' composite and 'rapids'."""
        mock_root.return_value = tmp_path
        mock_run.return_value = MagicMock(returncode=0)
        result = runner.invoke(main, ["install", "all"])
        assert result.exit_code == 0
        cmd_arg = mock_run.call_args[0][0]
        # Neither full nor rapids in the extras string
        extras_part = cmd_arg.split("[")[1].split("]")[0]
        assert "full" not in extras_part
        assert "rapids" not in extras_part

    @patch("subprocess.run")
    def test_install_prints_note(self, mock_run, runner):
        """Groups with notes should print the note after install."""
        mock_run.return_value = MagicMock(returncode=0)
        result = runner.invoke(main, ["install", "kronos"])
        assert result.exit_code == 0
        assert "KRONOS" in result.output


# ============================================================================
# CLI: workflow config
# ============================================================================


class TestWorkflowConfig:
    """Test ``kintsugi workflow config`` command."""

    def test_config_generates_yaml(self, runner, project_dir):
        """Config command creates workflow/config.yaml from experiment.json."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0, result.output

        config_file = project_dir / "workflow" / "config.yaml"
        assert config_file.exists()

        cfg = yaml.safe_load(config_file.read_text())
        assert cfg["project_dir"].replace("\\", "/") == project_dir.as_posix()
        assert cfg["cycles"] == [1, 2]
        assert cfg["channels"] == [1, 2, 3, 4]
        assert cfg["tile_rows"] == 3
        assert cfg["tile_cols"] == 3
        assert cfg["tile_overlap"] == 0.1

    def test_config_reads_experiment_json(self, runner, project_dir):
        """Cycles, channels, and tile grid come from experiment.json."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_dir / "workflow" / "config.yaml").read_text())
        assert cfg["cycles"] == [1, 2]  # n_cycles=2
        assert cfg["channels"] == [1, 2, 3, 4]  # channels_per_cycle=4

    def test_config_defaults_without_experiment_json(self, runner, tmp_path):
        """CLI now requires experiment.json and exits with error if absent."""
        bare_proj = tmp_path / "bare"
        bare_proj.mkdir()

        result = runner.invoke(main, ["workflow", "config", str(bare_proj)])
        assert result.exit_code != 0, "Should fail without experiment.json"
        # Check for an error message about missing experiment metadata
        output_lower = result.output.lower()
        assert "experiment" in output_lower or "meta" in output_lower or "not found" in output_lower

    def test_config_reads_slurm_config_sh(self, runner, project_with_slurm_config):
        """Config generates resources section (multi-account architecture)."""
        result = runner.invoke(main, ["workflow", "config", str(project_with_slurm_config)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_with_slurm_config / "workflow" / "config.yaml").read_text())
        # Multi-account architecture uses resources.accounts list, not flat account/partition
        assert "resources" in cfg
        res = cfg["resources"]
        # Should have standard resource keys (memory, time limits, etc.)
        assert any(
            key in res for key in ["accounts", "total_slots", "mem_stitch", "mem_decon", "mem_edf"]
        )

    def test_config_print_only(self, runner, project_dir):
        """--print-only outputs to stdout without writing a file."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir), "--print-only"])
        assert result.exit_code == 0
        assert "project_dir:" in result.output
        assert "cycles:" in result.output

        # Should NOT have created the file
        assert not (project_dir / "workflow" / "config.yaml").exists()

    def test_config_includes_edf_params(self, runner, project_dir):
        """Config includes EDF parameters with correct defaults."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_dir / "workflow" / "config.yaml").read_text())
        assert "edf" in cfg
        assert cfg["edf"]["radius_x"] == 2
        assert cfg["edf"]["radius_y"] == 2
        assert cfg["edf"]["sigma"] == 10.0

    def test_config_includes_quality_gate(self, runner, project_dir):
        """Config includes quality gate thresholds."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_dir / "workflow" / "config.yaml").read_text())
        assert "quality_gate" in cfg
        assert cfg["quality_gate"]["max_zero_pct"] == 0.10
        assert cfg["quality_gate"]["max_sat_pct"] == 0.05
        assert cfg["quality_gate"]["min_valid_slices"] == 3

    def test_config_includes_stitching_params(self, runner, project_dir):
        """Config includes BaSiC and stitching parameters."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_dir / "workflow" / "config.yaml").read_text())
        assert cfg["ncc_threshold"] == 0.078
        assert cfg["pou"] == 0.5
        assert cfg["blend_sigma"] == 15.0
        assert cfg["basic_flatfield_min"] == 0.1
        assert cfg["basic_max_iterations"] == 500

    def test_config_includes_resource_defaults(self, runner, project_dir):
        """Config has SLURM resource entries for each step."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        cfg = yaml.safe_load((project_dir / "workflow" / "config.yaml").read_text())
        res = cfg["resources"]
        for key in ["mem_stitch", "mem_decon", "mem_edf", "time_stitch", "time_decon", "time_edf"]:
            assert key in res, f"Missing resource key: {key}"

    def test_config_copies_snakefile(self, runner, project_dir):
        """Config copies Snakefile and scripts from the KINTSUGI installation."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        wf = project_dir / "workflow"
        assert (wf / "Snakefile").exists()
        assert (wf / "scripts" / "stitch.py").exists()
        assert (wf / "scripts" / "deconvolve.py").exists()
        assert (wf / "scripts" / "edf.py").exists()

    def test_config_copies_profile(self, runner, project_dir):
        """Config copies the SLURM profile directory."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0

        profile = project_dir / "workflow" / "profiles" / "slurm" / "config.yaml"
        assert profile.exists()
        profile_cfg = yaml.safe_load(profile.read_text())
        assert profile_cfg["executor"] == "slurm"

    def test_config_always_overwrites_snakefile(self, runner, project_with_workflow):
        """Config always overwrites Snakefile (by design, ensures latest rules)."""
        snakefile = project_with_workflow / "workflow" / "Snakefile"
        original = snakefile.read_text()

        result = runner.invoke(main, ["workflow", "config", str(project_with_workflow)])
        assert result.exit_code == 0

        # Snakefile SHOULD be overwritten with the full Snakefile from installation
        updated = snakefile.read_text()
        assert updated != original, "Snakefile should have been overwritten"
        # The overwritten Snakefile should contain real rule definitions
        assert "rule stitch:" in updated or "rule all:" in updated

    def test_config_next_steps_no_cd(self, runner, project_dir):
        """Next steps should NOT include 'cd workflow/' (path bug fix)."""
        result = runner.invoke(main, ["workflow", "config", str(project_dir)])
        assert result.exit_code == 0, result.output
        assert "cd " not in result.output or "cd " + str(project_dir) not in result.output
        # Should NOT tell user to cd into workflow dir
        assert "cd " + str(project_dir / "workflow") not in result.output
        # Should have proper kintsugi workflow commands
        assert "kintsugi workflow check" in result.output
        assert "kintsugi workflow run" in result.output
        # Should NOT suggest raw snakemake -n (use kintsugi wrapper instead)
        assert "snakemake -n" not in result.output


# ============================================================================
# CLI: workflow run
# ============================================================================


def _mock_gpu_pool():
    """Return a mock GPU resource pool for SLURM-mode tests.

    The CLI computes -j as: total_gpu_avail + min(4, total_cpu_slots).
    With 12 GPU + 4 CPU → 12 + 4 = 16, matching test_run_jobs_flag's -j 16.
    """
    return {
        "accounts": [
            {
                "name": "test_acct",
                "partition_gpu": "gpu",
                "partition_cpu": "hpg-default",
                "gpu_slots": 12,
                "cpu_slots": 4,
            }
        ],
        "total_gpu_slots": 12,
        "total_gpu_avail": 12,
        "total_cpu_slots": 4,
        "total_slots": 16,
    }


class TestWorkflowRun:
    """Test ``kintsugi workflow run`` command."""

    def test_run_fails_without_config_and_experiment(self, runner, tmp_path):
        """Run exits with error when both config.yaml and experiment.json are missing."""
        bare_proj = tmp_path / "bare"
        bare_proj.mkdir()
        result = runner.invoke(main, ["workflow", "run", str(bare_proj)])
        assert result.exit_code != 0
        assert "No workflow/config.yaml" in result.output
        assert "experiment.json" in result.output

    def test_run_fails_without_snakefile(self, runner, project_dir):
        """Run exits with error when Snakefile is missing."""
        wf = project_dir / "workflow"
        wf.mkdir(parents=True, exist_ok=True)
        (wf / "config.yaml").write_text("project_dir: /tmp\ncycles: [1]\n")

        result = runner.invoke(main, ["workflow", "run", str(project_dir)])
        assert result.exit_code != 0
        assert "No workflow/Snakefile" in result.output

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_invokes_snakemake(self, mock_run, mock_gpu, runner, project_with_workflow):
        """Run calls snakemake with correct arguments."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(main, ["workflow", "run", str(project_with_workflow)])
        assert result.exit_code == 0

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "snakemake"
        assert "--directory" in cmd

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_dry_run_flag(self, mock_run, mock_gpu, runner, project_with_workflow):
        """--dry-run passes -n to snakemake."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(main, ["workflow", "run", str(project_with_workflow), "--dry-run"])
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        assert "-n" in cmd

    @patch("subprocess.run")
    def test_run_local_mode(self, mock_run, runner, project_with_workflow):
        """--local runs with --cores instead of SLURM profile."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(
            main, ["workflow", "run", str(project_with_workflow), "--local", "--cores", "8"]
        )
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        assert "--cores" in cmd
        idx = cmd.index("--cores")
        assert cmd[idx + 1] == "8"
        assert "--profile" not in cmd

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_uses_slurm_profile(self, mock_run, mock_gpu, runner, project_with_workflow):
        """Default run uses the SLURM profile."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(main, ["workflow", "run", str(project_with_workflow)])
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        assert "--profile" in cmd
        profile_idx = cmd.index("--profile")
        assert "slurm" in cmd[profile_idx + 1]

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_forcerun(self, mock_run, mock_gpu, runner, project_with_workflow):
        """--forcerun passes through to snakemake."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(
            main, ["workflow", "run", str(project_with_workflow), "--forcerun", "stitch"]
        )
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        assert "--forcerun" in cmd
        assert "stitch" in cmd

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_cycle_range(self, mock_run, mock_gpu, runner, project_with_workflow):
        """--cycles '1-2' builds target paths for those cycles."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(
            main, ["workflow", "run", str(project_with_workflow), "--cycles", "1-2"]
        )
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        # Should contain target paths for cyc01 and cyc02
        target_args = [a for a in cmd if "cyc01" in a or "cyc02" in a]
        assert len(target_args) == 2

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_cycle_list(self, mock_run, mock_gpu, runner, project_with_workflow):
        """--cycles '1,2' builds target paths for those cycles."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(
            main, ["workflow", "run", str(project_with_workflow), "--cycles", "1,2"]
        )
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        target_args = [a for a in cmd if ".snakemake_complete" in a]
        assert len(target_args) == 2

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run", side_effect=FileNotFoundError)
    def test_run_snakemake_not_installed(self, mock_run, mock_gpu, runner, project_with_workflow):
        """Helpful error message when snakemake is not installed."""
        result = runner.invoke(main, ["workflow", "run", str(project_with_workflow)])
        assert result.exit_code != 0
        assert "snakemake not found" in result.output

    @patch("kintsugi.hpc.detect_live_multi_account", return_value=_mock_gpu_pool())
    @patch("subprocess.run")
    def test_run_jobs_flag(self, mock_run, mock_gpu, runner, project_with_workflow):
        """-j flag sets max concurrent SLURM jobs."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(main, ["workflow", "run", str(project_with_workflow), "-j", "16"])
        assert result.exit_code == 0

        cmd = mock_run.call_args[0][0]
        idx = cmd.index("-j")
        assert cmd[idx + 1] == "16"

    @patch("subprocess.run")
    @patch("kintsugi.cli.generate_workflow_config")
    def test_run_auto_generates_config(self, mock_gen, mock_run, runner, project_dir):
        """Run auto-generates config when experiment.json exists but config.yaml is missing."""
        wf_dir = project_dir / "workflow"

        # generate_workflow_config should create config + Snakefile
        def _fake_gen(proj_dir, print_only=False):
            wf_dir.mkdir(parents=True, exist_ok=True)
            (wf_dir / "config.yaml").write_text(
                yaml.dump({"project_dir": str(proj_dir), "cycles": [1, 2]})
            )
            (wf_dir / "Snakefile").write_text("rule all:\n    input: []\n")
            return wf_dir / "config.yaml"

        mock_gen.side_effect = _fake_gen
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(main, ["workflow", "run", str(project_dir), "--local"])
        assert result.exit_code == 0, result.output
        assert "Auto-generating" in result.output
        mock_gen.assert_called_once()
        mock_run.assert_called_once()

    def test_run_auto_config_fails_without_experiment(self, runner, tmp_path):
        """No config and no experiment.json gives a clear error."""
        bare_proj = tmp_path / "bare"
        bare_proj.mkdir()
        result = runner.invoke(main, ["workflow", "run", str(bare_proj)])
        assert result.exit_code != 0
        assert "experiment.json" in result.output

    @patch("subprocess.Popen")
    @patch("kintsugi.dashboard.scan_project_progress")
    @patch("kintsugi.dashboard.render_dashboard")
    def test_run_dashboard_uses_popen(
        self, mock_render, mock_scan, mock_popen, runner, project_with_workflow
    ):
        """--dashboard flag launches snakemake via Popen, not subprocess.run."""
        # Simulate snakemake finishing immediately
        mock_proc = MagicMock()
        mock_proc.pid = 12345
        mock_proc.poll.return_value = 0  # finished immediately
        mock_proc.stderr = MagicMock()
        mock_proc.stderr.read.return_value = b""
        mock_popen.return_value = mock_proc
        mock_scan.return_value = MagicMock()

        runner.invoke(
            main, ["workflow", "run", str(project_with_workflow), "--dashboard", "--local"]
        )
        # Popen should have been called (not subprocess.run)
        mock_popen.assert_called_once()
        # scan_project_progress should have been called at least once (final render)
        assert mock_scan.called

    @patch("subprocess.run")
    def test_run_dashboard_dry_run_falls_through(self, mock_run, runner, project_with_workflow):
        """--dashboard with --dry-run falls through to blocking subprocess.run."""
        mock_run.return_value = MagicMock(returncode=0)

        result = runner.invoke(
            main,
            [
                "workflow",
                "run",
                str(project_with_workflow),
                "--dashboard",
                "--dry-run",
                "--local",
            ],
        )
        assert result.exit_code == 0
        # Should use subprocess.run (blocking), not Popen
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert "-n" in cmd


# ============================================================================
# CLI: workflow help
# ============================================================================


class TestWorkflowHelp:
    """Basic CLI integration tests."""

    def test_workflow_group_help(self, runner):
        """workflow --help works."""
        result = runner.invoke(main, ["workflow", "--help"])
        assert result.exit_code == 0
        assert "config" in result.output
        assert "run" in result.output

    def test_workflow_config_help(self, runner):
        """workflow config --help works."""
        result = runner.invoke(main, ["workflow", "config", "--help"])
        assert result.exit_code == 0
        assert "PROJECT_DIR" in result.output

    def test_workflow_run_help(self, runner):
        """workflow run --help works."""
        result = runner.invoke(main, ["workflow", "run", "--help"])
        assert result.exit_code == 0
        assert "--dry-run" in result.output
        assert "--local" in result.output
        assert "--dashboard" in result.output


# ============================================================================
# Snakefile Validation
# ============================================================================


class TestSnakefileValidation:
    """Validate Snakefile structure and config.yaml schema."""

    def test_snakefile_exists(self):
        """workflow/Snakefile exists in the repo."""
        snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
        assert snakefile.exists()

    def test_snakefile_has_required_rules(self):
        """Snakefile defines stitch, deconvolve, edf, and all rules."""
        snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
        content = snakefile.read_text()
        for rule_name in ["rule all:", "rule stitch:", "rule deconvolve:", "rule edf:"]:
            assert rule_name in content, f"Missing {rule_name}"

    def test_snakefile_sentinel_dependencies(self):
        """Each rule's input depends on the prior step's sentinel."""
        snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
        content = snakefile.read_text()

        # deconvolve depends on stitch sentinel
        assert "stitched/cyc{cycle}/.snakemake_complete" in content
        # edf depends on deconvolve sentinel
        assert "deconvolved/cyc{cycle}/.snakemake_complete" in content

    def test_snakefile_has_cyc_fmt_helper(self):
        """Snakefile includes the cyc_fmt() helper for zero-padded cycles."""
        snakefile = Path(__file__).parent.parent / "workflow" / "Snakefile"
        content = snakefile.read_text()
        assert "def cyc_fmt" in content

    def test_config_yaml_template_is_valid(self):
        """workflow/config.yaml is valid YAML with expected keys."""
        config_file = Path(__file__).parent.parent / "workflow" / "config.yaml"
        assert config_file.exists()

        cfg = yaml.safe_load(config_file.read_text())
        expected_keys = [
            "project_dir",
            "kintsugi_dir",
            "cycles",
            "channels",
            "tile_rows",
            "tile_cols",
            "edf",
            "quality_gate",
            "resources",
        ]
        for key in expected_keys:
            assert key in cfg, f"Missing key in config.yaml: {key}"

    def test_config_edf_section(self):
        """EDF config section has all required fields."""
        config_file = Path(__file__).parent.parent / "workflow" / "config.yaml"
        cfg = yaml.safe_load(config_file.read_text())

        edf = cfg["edf"]
        for key in ["radius_x", "radius_y", "sigma", "tiles", "blend_depth", "z_smooth_sigma"]:
            assert key in edf, f"Missing edf key: {key}"

    def test_config_quality_gate_section(self):
        """Quality gate section has all threshold fields."""
        config_file = Path(__file__).parent.parent / "workflow" / "config.yaml"
        cfg = yaml.safe_load(config_file.read_text())

        qg = cfg["quality_gate"]
        for key in ["max_zero_pct", "max_sat_pct", "min_valid_slices"]:
            assert key in qg, f"Missing quality_gate key: {key}"

    def test_slurm_profile_is_valid(self):
        """SLURM profile is valid YAML with expected settings."""
        profile = Path(__file__).parent.parent / "workflow" / "profiles" / "slurm" / "config.yaml"
        assert profile.exists()

        cfg = yaml.safe_load(profile.read_text())
        assert cfg["executor"] == "slurm"
        assert "default-resources" in cfg
        assert cfg["latency-wait"] >= 60
        assert cfg["retries"] >= 1

    def test_slurm_profile_no_redundant_j_flag(self):
        """Profile header no longer includes -j 8 in usage examples."""
        profile = Path(__file__).parent.parent / "workflow" / "profiles" / "slurm" / "config.yaml"
        content = profile.read_text()
        # The CoPilot review fix removed "-j 8" from comments
        assert "-j 8" not in content

    def test_conda_env_yaml_exists(self):
        """Conda environment spec exists (GPU variant)."""
        env_file = Path(__file__).parent.parent / "workflow" / "envs" / "kintsugi_gpu.yaml"
        assert env_file.exists()

        cfg = yaml.safe_load(env_file.read_text())
        assert "dependencies" in cfg


# ============================================================================
# Workflow Script Helpers
# ============================================================================


class TestWorkflowScripts:
    """Tests for helper functions extracted from workflow scripts."""

    def test_scripts_exist(self):
        """All three workflow scripts exist."""
        scripts_dir = Path(__file__).parent.parent / "workflow" / "scripts"
        for name in ["stitch.py", "deconvolve.py", "edf.py"]:
            assert (scripts_dir / name).exists(), f"Missing {name}"

    def test_stitch_script_references_snakemake(self):
        """Stitch script uses snakemake.params, wildcards, output."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "stitch.py"
        content = script.read_text()
        assert "snakemake.params.project_dir" in content
        assert "snakemake.wildcards.cycle" in content
        assert "snakemake.output.sentinel" in content

    def test_deconvolve_script_references_snakemake(self):
        """Deconvolve script uses snakemake interface."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "deconvolve.py"
        content = script.read_text()
        assert "snakemake.params.project_dir" in content
        assert "snakemake.wildcards.cycle" in content
        assert "snakemake.output.sentinel" in content

    def test_edf_script_references_snakemake(self):
        """EDF script uses snakemake interface."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "edf.py"
        content = script.read_text()
        assert "snakemake.params.project_dir" in content
        assert "snakemake.wildcards.cycle" in content
        assert "snakemake.output.sentinel" in content

    def test_scripts_no_wait_for_input(self):
        """No script calls wait_for_input() (Snakemake handles deps).

        Docstring/comment mentions are fine — only actual function calls
        or imports are disallowed.
        """
        import re as _re

        scripts_dir = Path(__file__).parent.parent / "workflow" / "scripts"
        for name in ["stitch.py", "deconvolve.py", "edf.py"]:
            content = (scripts_dir / name).read_text()
            # Strip comments and docstrings: search only non-string, non-comment code
            # Simple approach: check for lines that call or define wait_for_input
            # excluding lines inside triple-quoted strings or starting with #
            lines = content.splitlines()
            code_lines = []
            in_docstring = False
            for line in lines:
                stripped = line.strip()
                if stripped.startswith('"""') or stripped.startswith("'''"):
                    in_docstring = not in_docstring
                    continue
                if in_docstring or stripped.startswith("#"):
                    continue
                code_lines.append(line)
            code_text = "\n".join(code_lines)
            calls = _re.findall(r"\bwait_for_input\s*\(", code_text)
            imports = _re.findall(
                r"^\s*(?:from|import)\b.*wait_for_input", code_text, _re.MULTILINE
            )
            assert not calls and not imports, f"{name} still calls wait_for_input()"

    def test_scripts_no_skip_existing(self):
        """No script has skip-existing logic (Snakemake handles re-runs)."""
        scripts_dir = Path(__file__).parent.parent / "workflow" / "scripts"
        for name in ["stitch.py", "deconvolve.py", "edf.py"]:
            content = (scripts_dir / name).read_text()
            assert "skip_existing" not in content.lower(), f"{name} has skip-existing"
            assert "already exists" not in content.lower(), f"{name} has skip-existing check"

    def test_scripts_write_sentinel(self):
        """Each script writes a sentinel file with stage metadata."""
        scripts_dir = Path(__file__).parent.parent / "workflow" / "scripts"
        for name in ["stitch.py", "deconvolve.py", "edf.py"]:
            content = (scripts_dir / name).read_text()
            assert "sentinel" in content, f"{name} missing sentinel write"
            assert "snakemake.output.sentinel" in content

    def test_stitch_imports_correct_modules(self):
        """Stitch script imports from verified KINTSUGI modules."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "stitch.py"
        content = script.read_text()
        assert "from Kstitch.stitching import stitch_images" in content
        assert "from kintsugi.kcorrect_gpu import KCorrectGPU" in content
        assert "from kintsugi.stitch_blend import stitch_with_blending" in content

    def test_deconvolve_imports_correct_modules(self):
        """Deconvolve script imports from verified KINTSUGI modules."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "deconvolve.py"
        content = script.read_text()
        assert "from Kdecon import decon" in content

    def test_edf_imports_correct_modules(self):
        """EDF script imports from verified KINTSUGI modules."""
        script = Path(__file__).parent.parent / "workflow" / "scripts" / "edf.py"
        content = script.read_text()
        assert "from kintsugi.edf import EDFProcessor" in content


# ============================================================================
# Quality Gate Logic (extracted from edf.py for testability)
# ============================================================================


class TestQualityGate:
    """Test quality gate logic from edf.py (tested via pure numpy)."""

    def _quality_gate_zstack(self, stack, max_zero_pct=0.10, max_sat_pct=0.05):
        """Reference implementation matching edf.py's quality_gate_zstack."""
        valid_z_indices = []
        excluded_reasons = []
        for z_idx in range(stack.shape[0]):
            slice_data = stack[z_idx]
            total_pixels = slice_data.size
            pct_zero = np.sum(slice_data == 0) / total_pixels
            pct_sat = np.sum(slice_data >= 64000) / total_pixels
            if pct_zero > max_zero_pct:
                excluded_reasons.append(f"z={z_idx + 1}: {pct_zero * 100:.1f}% zeros")
                continue
            if pct_sat > max_sat_pct:
                excluded_reasons.append(f"z={z_idx + 1}: {pct_sat * 100:.1f}% saturated")
                continue
            valid_z_indices.append(z_idx)
        return valid_z_indices, excluded_reasons

    def test_all_clean_slices(self):
        """All slices pass when there are no artifacts."""
        rng = np.random.default_rng(42)
        stack = rng.integers(100, 30000, size=(10, 64, 64), dtype=np.uint16)
        valid, excluded = self._quality_gate_zstack(stack)
        assert len(valid) == 10
        assert len(excluded) == 0

    def test_zero_slices_excluded(self):
        """Slices with >10% zeros are excluded."""
        rng = np.random.default_rng(42)
        stack = rng.integers(100, 30000, size=(5, 100, 100), dtype=np.uint16)
        # Make slice 2 mostly zeros
        stack[2, :, :] = 0
        valid, excluded = self._quality_gate_zstack(stack)
        assert 2 not in valid
        assert len(excluded) == 1

    def test_saturated_slices_excluded(self):
        """Slices with >5% saturation are excluded."""
        rng = np.random.default_rng(42)
        stack = rng.integers(100, 30000, size=(5, 100, 100), dtype=np.uint16)
        # Make slice 3 mostly saturated
        stack[3, :, :] = 65000
        valid, excluded = self._quality_gate_zstack(stack)
        assert 3 not in valid
        assert len(excluded) == 1

    def test_multiple_bad_slices(self):
        """Multiple bad slices are independently excluded."""
        rng = np.random.default_rng(42)
        stack = rng.integers(100, 30000, size=(6, 100, 100), dtype=np.uint16)
        stack[0, :, :] = 0  # zeros
        stack[4, :, :] = 65000  # saturated
        valid, excluded = self._quality_gate_zstack(stack)
        assert 0 not in valid
        assert 4 not in valid
        assert len(valid) == 4
        assert len(excluded) == 2

    def test_borderline_zero_passes(self):
        """Slice at exactly max_zero_pct threshold passes."""
        stack = np.ones((3, 100, 100), dtype=np.uint16) * 500
        # Set exactly 10% of pixels to zero (threshold is > 0.10, so == passes)
        stack[1, :10, :] = 0  # 10 rows * 100 cols = 1000 pixels = 10%
        valid, excluded = self._quality_gate_zstack(stack, max_zero_pct=0.10)
        assert 1 in valid  # exactly 10% passes (> not >=)


# ============================================================================
# Find Cycle Directory Logic (extracted from stitch.py)
# ============================================================================


class TestFindCycleDir:
    """Test cycle directory detection logic from stitch.py."""

    def _find_cycle_dir(self, raw_dir, cycle_num):
        """Reference implementation matching stitch.py's find_cycle_dir."""

        # List actual directory entries so we match by real (case-preserving)
        # name.  Using Path.exists() directly is unreliable on case-insensitive
        # filesystems (macOS, Windows) where "cyc01" matches "Cyc01".
        entries = {p.name: p for p in raw_dir.iterdir() if p.is_dir()}

        exact_patterns = [
            f"cyc{cycle_num:03d}",
            f"cyc{cycle_num:02d}",
            f"Cyc{cycle_num:02d}",
        ]
        for pattern in exact_patterns:
            if pattern in entries:
                return entries[pattern]

        glob_patterns = [
            f"cyc{cycle_num:03d}_*",
            f"cyc{cycle_num:02d}_*",
            f"Cyc{cycle_num:02d}_*",
        ]
        for pattern in glob_patterns:
            matches = list(raw_dir.glob(pattern))
            if matches:
                return matches[0]

        return None

    def test_short_format(self, tmp_path):
        """Finds cyc01/ style directories."""
        (tmp_path / "cyc01").mkdir()
        assert self._find_cycle_dir(tmp_path, 1).name == "cyc01"

    def test_long_format(self, tmp_path):
        """Finds cyc001/ style directories."""
        (tmp_path / "cyc001").mkdir()
        assert self._find_cycle_dir(tmp_path, 1).name == "cyc001"

    def test_long_form_with_suffix(self, tmp_path):
        """Finds cyc001_DAPI_CD3/ style directories."""
        (tmp_path / "cyc001_DAPI_CD3_CD8").mkdir()
        result = self._find_cycle_dir(tmp_path, 1)
        assert result is not None
        assert result.name.startswith("cyc001")

    def test_uppercase_format(self, tmp_path):
        """Finds Cyc01/ style directories."""
        (tmp_path / "Cyc01").mkdir()
        assert self._find_cycle_dir(tmp_path, 1).name == "Cyc01"

    def test_not_found(self, tmp_path):
        """Returns None when cycle directory doesn't exist."""
        assert self._find_cycle_dir(tmp_path, 99) is None

    def test_prefers_exact_over_glob(self, tmp_path):
        """Exact directory name is preferred over glob with suffix."""
        (tmp_path / "cyc001_markers").mkdir()
        (tmp_path / "cyc001").mkdir()
        result = self._find_cycle_dir(tmp_path, 1)
        # Exact pattern cyc001 matches first (avoids case-insensitive glob issues)
        assert result.name == "cyc001"

    def test_falls_back_to_glob(self, tmp_path):
        """Falls back to glob when no exact match exists."""
        (tmp_path / "cyc001_DAPI_CD3").mkdir()
        result = self._find_cycle_dir(tmp_path, 1)
        assert result is not None
        assert result.name == "cyc001_DAPI_CD3"


# ============================================================================
# Sentinel File Format
# ============================================================================


class TestSentinelFormat:
    """Test sentinel file metadata format."""

    def _parse_sentinel(self, text):
        """Parse key=value sentinel content."""
        data = {}
        for line in text.strip().splitlines():
            if "=" in line:
                key, value = line.split("=", 1)
                data[key.strip()] = value.strip()
        return data

    def test_stitch_sentinel_format(self):
        """Stitch sentinel has expected metadata fields."""
        sentinel_text = (
            "stage=stitch\n"
            "cycle=1\n"
            "completed=2026-01-15T10:30:00\n"
            "channels=1-4\n"
            "zplanes=15\n"
            "duration_minutes=12.3\n"
        )
        data = self._parse_sentinel(sentinel_text)
        assert data["stage"] == "stitch"
        assert data["cycle"] == "1"
        assert "completed" in data
        assert data["channels"] == "1-4"

    def test_decon_sentinel_format(self):
        """Decon sentinel has expected metadata fields."""
        sentinel_text = (
            "stage=decon\n"
            "cycle=3\n"
            "completed=2026-01-15T11:00:00\n"
            "channels=1-4\n"
            "successful=4\n"
            "duration_minutes=45.2\n"
        )
        data = self._parse_sentinel(sentinel_text)
        assert data["stage"] == "decon"
        assert data["successful"] == "4"

    def test_edf_sentinel_format(self):
        """EDF sentinel has expected metadata fields."""
        sentinel_text = (
            "stage=edf\n"
            "cycle=1\n"
            "completed=2026-01-15T12:00:00\n"
            "channels=1-4\n"
            "successful=4\n"
            "duration_minutes=5.1\n"
        )
        data = self._parse_sentinel(sentinel_text)
        assert data["stage"] == "edf"
        assert data["successful"] == "4"

    def test_sentinel_scripts_have_matching_fields(self):
        """Each script writes the stage field matching its expected value."""
        scripts_dir = Path(__file__).parent.parent / "workflow" / "scripts"
        expected_stages = {
            "stitch.py": "stitch",
            "deconvolve.py": "decon",
            "edf.py": "edf",
        }
        for script_name, expected_stage in expected_stages.items():
            content = (scripts_dir / script_name).read_text()
            assert f"stage={expected_stage}" in content, (
                f"{script_name} missing stage={expected_stage}"
            )
