"""
Tests for the Snakemake progress dashboard module.

Covers:
  - Project progress scanning (sentinel detection, stage statuses)
  - Log timing extraction
  - Completion time estimation
  - SLURM job name parsing
  - Memory formatting
  - Elapsed time parsing
  - Multi-project scanning
  - Dashboard rendering (smoke test)
"""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from kintsugi.dashboard import (
    CycleStatus,
    HardwareStatus,
    JobInfo,
    ProjectStatus,
    _check_registration_status,
    _check_signal_isolation_status,
    _check_stage_status,
    _format_mem,
    _normalize_rule,
    _parse_elapsed,
    _parse_job_name,
    _parse_log_duration,
    estimate_completion,
    render_dashboard,
    scan_project_progress,
)

# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def project_dir(tmp_path):
    """Create a minimal project with workflow config and processed data."""
    proj = tmp_path / "project"
    (proj / "meta").mkdir(parents=True)
    (proj / "data" / "raw" / "cyc01").mkdir(parents=True)
    (proj / "data" / "raw" / "cyc02").mkdir(parents=True)
    (proj / "data" / "raw" / "cyc03").mkdir(parents=True)
    (proj / "data" / "processed").mkdir(parents=True)
    (proj / "workflow").mkdir(parents=True)
    (proj / "slurm" / "logs" / "snakemake").mkdir(parents=True)
    (proj / "qc_plots").mkdir(parents=True)

    config = {
        "project_dir": str(proj),
        "project_name": "TestProject",
        "kintsugi_dir": "/opt/KINTSUGI",
        "cycles": [1, 2, 3],
        "channels": [1, 2, 3, 4],
        "n_zplanes": 10,
        "channel_names": {1: ["DAPI", "Blank", "Blank", "Blank"]},
        "resources": {},
    }
    with open(proj / "workflow" / "config.yaml", "w") as f:
        yaml.dump(config, f)

    return proj


@pytest.fixture
def completed_project(project_dir):
    """Project with all stages complete for cycles 1-2, pending for cycle 3."""
    for cyc in ["01", "02"]:
        for stage_dir in ["stitched", "deconvolved", "edf"]:
            d = project_dir / "data" / "processed" / stage_dir / f"cyc{cyc}"
            d.mkdir(parents=True, exist_ok=True)
            (d / ".snakemake_complete").touch()
            (d / "CH1.tif").touch()

    # Cycle 03: stitch done, decon partial, edf pending
    d = project_dir / "data" / "processed" / "stitched" / "cyc03"
    d.mkdir(parents=True, exist_ok=True)
    (d / ".snakemake_complete").touch()

    d = project_dir / "data" / "processed" / "deconvolved" / "cyc03"
    d.mkdir(parents=True, exist_ok=True)
    (d / "CH1.tif").touch()  # Partial — no sentinel

    # Registration done
    reg_dir = project_dir / "data" / "processed" / "registered"
    reg_dir.mkdir(parents=True, exist_ok=True)
    (reg_dir / ".snakemake_complete").touch()

    # Signal isolation done
    si_dir = project_dir / "data" / "processed" / "signal_isolated"
    si_dir.mkdir(parents=True, exist_ok=True)
    (si_dir / ".snakemake_complete").touch()

    # QC done for stitch and decon
    (project_dir / "qc_plots" / ".snakemake_complete_stitch").touch()
    (project_dir / "qc_plots" / ".snakemake_complete_decon").touch()

    return project_dir


# ============================================================================
# Stage status detection
# ============================================================================


class TestStageStatus:
    def test_done_with_sentinel(self, project_dir):
        d = project_dir / "data" / "processed" / "stitched" / "cyc01"
        d.mkdir(parents=True)
        (d / ".snakemake_complete").touch()
        assert _check_stage_status(project_dir, "stitch", "01") == "done"

    def test_partial_without_sentinel(self, project_dir):
        d = project_dir / "data" / "processed" / "deconvolved" / "cyc01"
        d.mkdir(parents=True)
        (d / "CH1.tif").touch()
        assert _check_stage_status(project_dir, "deconvolve", "01") == "partial"

    def test_pending_no_directory(self, project_dir):
        assert _check_stage_status(project_dir, "edf", "01") == "pending"

    def test_pending_empty_directory(self, project_dir):
        d = project_dir / "data" / "processed" / "edf" / "cyc01"
        d.mkdir(parents=True)
        assert _check_stage_status(project_dir, "edf", "01") == "pending"


class TestRegistrationStatus:
    def test_done(self, project_dir):
        d = project_dir / "data" / "processed" / "registered"
        d.mkdir(parents=True)
        (d / ".snakemake_complete").touch()
        assert _check_registration_status(project_dir) == "done"

    def test_partial(self, project_dir):
        d = project_dir / "data" / "processed" / "registered" / "cyc01"
        d.mkdir(parents=True)
        assert _check_registration_status(project_dir) == "partial"

    def test_pending(self, project_dir):
        assert _check_registration_status(project_dir) == "pending"


class TestSignalIsolationStatus:
    def test_done(self, project_dir):
        d = project_dir / "data" / "processed" / "signal_isolated"
        d.mkdir(parents=True)
        (d / ".snakemake_complete").touch()
        assert _check_signal_isolation_status(project_dir) == "done"

    def test_partial_with_manifest(self, project_dir):
        d = project_dir / "data" / "processed" / "signal_isolated"
        d.mkdir(parents=True)
        (d / "signal_isolation_manifest.json").write_text("{}")
        assert _check_signal_isolation_status(project_dir) == "partial"

    def test_partial_with_tifs(self, project_dir):
        d = project_dir / "data" / "processed" / "signal_isolated"
        d.mkdir(parents=True)
        (d / "CD3e.tif").touch()
        assert _check_signal_isolation_status(project_dir) == "partial"

    def test_pending(self, project_dir):
        assert _check_signal_isolation_status(project_dir) == "pending"


# ============================================================================
# Job name parsing
# ============================================================================


class TestJobNameParsing:
    def test_smk_stitch(self):
        rule, cycle = _parse_job_name("smk-stitch-cyc03")
        assert rule == "stitch"
        assert cycle == "03"

    def test_snakejob_decon(self):
        rule, cycle = _parse_job_name("snakejob.deconvolve.12")
        assert rule == "deconvolve"
        assert cycle == "12"

    def test_registration_no_cycle(self):
        rule, cycle = _parse_job_name("smk-registration")
        assert rule == "registration"
        assert cycle == ""

    def test_group_edf(self):
        rule, cycle = _parse_job_name("group_edf_05")
        assert rule == "edf"
        assert cycle == "05"

    def test_signal_isolation(self):
        rule, cycle = _parse_job_name("smk-signal_isolation")
        assert rule == "signal_isolation"
        assert cycle == ""

    def test_qc_signal_isolation(self):
        rule, cycle = _parse_job_name("smk-qc_signal_isolation")
        assert rule == "qc_signal_isolation"
        assert cycle == ""

    def test_qc_rule(self):
        rule, cycle = _parse_job_name("smk-qc_stitch")
        assert rule == "qc_stitch"
        assert cycle == ""

    def test_unrelated_job(self):
        rule, cycle = _parse_job_name("jupyter-notebook")
        assert rule == ""
        assert cycle == ""


class TestNormalizeRule:
    def test_decon_alias(self):
        assert _normalize_rule("decon") == "deconvolve"

    def test_identity(self):
        assert _normalize_rule("stitch") == "stitch"
        assert _normalize_rule("edf") == "edf"
        assert _normalize_rule("registration") == "registration"

    def test_signal_isolation(self):
        assert _normalize_rule("signal_isolation") == "signal_isolation"

    def test_unknown_passthrough(self):
        assert _normalize_rule("vessel3d") == "vessel3d"


# ============================================================================
# Elapsed time parsing
# ============================================================================


class TestParseElapsed:
    def test_hhmmss(self):
        assert _parse_elapsed("01:23:45") == 5025.0

    def test_mmss(self):
        assert _parse_elapsed("05:30") == 330.0

    def test_with_days(self):
        assert _parse_elapsed("1-02:00:00") == 93600.0

    def test_empty(self):
        assert _parse_elapsed("") == 0.0

    def test_invalid(self):
        assert _parse_elapsed("invalid") == 0.0

    def test_non_numeric_days(self):
        """Non-numeric day part should return 0.0 instead of raising."""
        assert _parse_elapsed("abc-01:00:00") == 0.0

    def test_non_numeric_hours(self):
        assert _parse_elapsed("xx:30:00") == 0.0

    def test_non_numeric_in_day_section(self):
        """Day section like '1.5-02:00:00' should handle gracefully."""
        assert _parse_elapsed("1.5-02:00:00") == 0.0


# ============================================================================
# Log timing extraction
# ============================================================================


class TestLogDuration:
    def test_duration_pattern(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text(
            "KINTSUGI stitch - 2026-02-22 10:00:00\n"
            "Processing...\n"
            "Duration: 12m 34s\n"
            "Exit code: 0\n"
        )
        assert _parse_log_duration(log) == 754.0

    def test_timestamp_fallback(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text(
            "KINTSUGI stitch - 2026-02-22 10:00:00\nProcessing...\nCompleted: 2026-02-22 10:08:30\n"
        )
        assert _parse_log_duration(log) == 510.0

    def test_empty_log(self, tmp_path):
        log = tmp_path / "test.log"
        log.write_text("")
        assert _parse_log_duration(log) == 0.0

    def test_missing_file(self, tmp_path):
        log = tmp_path / "nonexistent.log"
        assert _parse_log_duration(log) == 0.0


# ============================================================================
# Memory formatting
# ============================================================================


class TestFormatMem:
    def test_megabytes(self):
        assert _format_mem("4000M") == "3.9G"

    def test_gigabytes(self):
        assert _format_mem("48G") == "48.0G"

    def test_raw_number(self):
        assert _format_mem("4000") == "3.9G"

    def test_small_mb(self):
        assert _format_mem("512M") == "512M"

    def test_empty(self):
        assert _format_mem("") == ""

    def test_kilobytes(self):
        assert _format_mem("2048K") == "2M"


# ============================================================================
# Completion time estimation
# ============================================================================


class TestEstimateCompletion:
    def test_all_done(self, completed_project):
        status = _make_fully_completed_status(completed_project)
        est = estimate_completion(status)
        # Registration is done, all cycles stitch/decon/edf done
        assert est["wall_clock_seconds"] == 0
        assert est["eta"] is None

    def test_pending_stages(self, project_dir):
        status = _make_empty_status(project_dir)
        est = estimate_completion(status)
        assert est["wall_clock_seconds"] > 0
        assert est["eta"] is not None
        assert est["confidence"] == "low"  # No historical data

    def test_with_timing_history(self, project_dir):
        status = _make_empty_status(project_dir)
        # Add timing history
        status.cycle_statuses["01"].stages["stitch"] = "done"
        status.cycle_statuses["01"].timings["stitch"] = 300.0
        est = estimate_completion(status)
        assert est["confidence"] == "medium"  # One stage has history

    def test_parallelism_from_running(self, project_dir):
        status = _make_empty_status(project_dir)
        # Mark two jobs as running
        status.cycle_statuses["01"].stages["stitch"] = "running"
        status.cycle_statuses["02"].stages["stitch"] = "running"
        status.cycle_statuses["01"].jobs["stitch"] = JobInfo(
            job_id="123", state="RUNNING", elapsed="00:05:00"
        )
        status.cycle_statuses["02"].jobs["stitch"] = JobInfo(
            job_id="124", state="RUNNING", elapsed="00:03:00"
        )
        est = estimate_completion(status)
        assert est["parallelism"] == 2


# ============================================================================
# Full project scanning
# ============================================================================


class TestScanProject:
    @patch("kintsugi.dashboard._attach_slurm_jobs")
    @patch("kintsugi.dashboard._attach_log_timings")
    def test_scan_basic(self, mock_timings, mock_jobs, completed_project):
        status = scan_project_progress(completed_project)

        assert status.project_name == "TestProject"
        assert len(status.cycles) == 3
        assert status.cycle_statuses["01"].stages["stitch"] == "done"
        assert status.cycle_statuses["01"].stages["deconvolve"] == "done"
        assert status.cycle_statuses["01"].stages["edf"] == "done"
        assert status.cycle_statuses["03"].stages["stitch"] == "done"
        assert status.cycle_statuses["03"].stages["deconvolve"] == "partial"
        assert status.cycle_statuses["03"].stages["edf"] == "pending"
        assert status.registration_status == "done"
        assert status.signal_isolation_status == "done"
        assert status.qc_statuses["stitch"] == "done"
        assert status.qc_statuses["decon"] == "done"
        assert status.qc_statuses["edf"] == "pending"
        assert status.qc_statuses["registration"] == "pending"
        assert status.qc_statuses["signal_isolation"] == "pending"

    @patch("kintsugi.dashboard._attach_slurm_jobs")
    @patch("kintsugi.dashboard._attach_log_timings")
    def test_scan_empty_project(self, mock_timings, mock_jobs, project_dir):
        status = scan_project_progress(project_dir)

        assert len(status.cycles) == 3
        for cyc in ["01", "02", "03"]:
            for stage in ("stitch", "deconvolve", "edf"):
                assert status.cycle_statuses[cyc].stages[stage] == "pending"
        assert status.registration_status == "pending"

    @patch("kintsugi.dashboard._attach_slurm_jobs")
    @patch("kintsugi.dashboard._attach_log_timings")
    def test_scan_non_numeric_cycles(self, mock_timings, mock_jobs, tmp_path):
        """Non-numeric cycle values in config should not crash scanner."""
        proj = tmp_path / "bad_cycles"
        (proj / "workflow").mkdir(parents=True)
        (proj / "data" / "processed").mkdir(parents=True)
        config = {
            "project_name": "BadCycles",
            "cycles": ["abc", "def"],
        }
        (proj / "workflow" / "config.yaml").write_text(yaml.dump(config))
        status = scan_project_progress(proj)
        assert len(status.cycles) == 0


# ============================================================================
# Dashboard rendering (smoke test)
# ============================================================================


class TestRenderDashboard:
    def test_render_no_crash(self, completed_project):
        """Verify render_dashboard doesn't raise on valid input."""
        status = _make_completed_status(completed_project)
        # Should not raise
        render_dashboard(
            [status],
            show_hardware=False,
            show_estimates=True,
            show_jobs=True,
        )

    def test_render_empty_project(self, project_dir):
        status = _make_empty_status(project_dir)
        render_dashboard(
            [status],
            show_hardware=False,
            show_estimates=True,
            show_jobs=False,
        )

    def test_render_with_hardware(self, project_dir):
        status = _make_empty_status(project_dir)
        status.hardware = HardwareStatus(
            accounts=[
                {
                    "name": "maigan",
                    "gpu_alloc": 3,
                    "gpu_used": 2,
                    "gpu_avail": 1,
                    "cpu_alloc": 104,
                    "cpu_used": 32,
                    "mem_alloc_gb": 812,
                    "mem_used_gb": 256,
                },
            ],
            total_gpu_alloc=3,
            total_gpu_used=2,
            total_gpu_avail=1,
            total_cpu_alloc=104,
            total_cpu_used=32,
            total_mem_alloc_gb=812,
            total_mem_used_gb=256,
        )
        render_dashboard(
            [status],
            show_hardware=True,
            show_estimates=False,
            show_jobs=False,
        )


# ============================================================================
# CLI integration
# ============================================================================


class TestStatusCLI:
    @patch("kintsugi.dashboard._attach_slurm_jobs")
    @patch("kintsugi.dashboard._attach_log_timings")
    def test_json_output(self, mock_timings, mock_jobs, completed_project):
        import json

        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "status", str(completed_project), "--json"])
        assert result.exit_code == 0
        data = json.loads(result.output)
        assert len(data) == 1
        assert data[0]["project_name"] == "TestProject"
        assert "cycle_statuses" in data[0]
        assert "estimate" in data[0]

    def test_missing_config(self, tmp_path):
        from click.testing import CliRunner

        from kintsugi.cli import main

        runner = CliRunner()
        result = runner.invoke(main, ["workflow", "status", str(tmp_path)])
        assert result.exit_code != 0
        assert "config.yaml" in result.output


# ============================================================================
# Helpers
# ============================================================================


def _make_fully_completed_status(project_dir: Path) -> ProjectStatus:
    """Build a ProjectStatus where ALL cycles, registration, and signal isolation are done."""
    status = ProjectStatus(
        project_dir=project_dir,
        project_name="TestProject",
        cycles=["01", "02", "03"],
    )

    for cyc in ["01", "02", "03"]:
        cs = CycleStatus(cycle=cyc)
        for stage in ("stitch", "deconvolve", "edf"):
            cs.stages[stage] = "done"
            cs.timings[stage] = 300.0
        status.cycle_statuses[cyc] = cs

    status.registration_status = "done"
    status.registration_timing = 1200.0
    status.signal_isolation_status = "done"
    status.signal_isolation_timing = 600.0
    status.qc_statuses = {
        "stitch": "done",
        "decon": "done",
        "edf": "done",
        "registration": "done",
        "signal_isolation": "done",
    }

    return status


def _make_completed_status(project_dir: Path) -> ProjectStatus:
    """Build a ProjectStatus where cycles 1-2 are done, 3 is partial."""
    status = ProjectStatus(
        project_dir=project_dir,
        project_name="TestProject",
        cycles=["01", "02", "03"],
    )

    for cyc in ["01", "02"]:
        cs = CycleStatus(cycle=cyc)
        for stage in ("stitch", "deconvolve", "edf"):
            cs.stages[stage] = "done"
            cs.timings[stage] = 300.0
        status.cycle_statuses[cyc] = cs

    cs3 = CycleStatus(cycle="03")
    cs3.stages["stitch"] = "done"
    cs3.stages["deconvolve"] = "partial"
    cs3.stages["edf"] = "pending"
    cs3.timings["stitch"] = 280.0
    status.cycle_statuses["03"] = cs3

    status.registration_status = "done"
    status.registration_timing = 1200.0
    status.signal_isolation_status = "done"
    status.signal_isolation_timing = 600.0
    status.qc_statuses = {
        "stitch": "done",
        "decon": "done",
        "edf": "pending",
        "registration": "pending",
        "signal_isolation": "pending",
    }

    return status


def _make_empty_status(project_dir: Path) -> ProjectStatus:
    """Build a ProjectStatus where everything is pending."""
    status = ProjectStatus(
        project_dir=project_dir,
        project_name="TestProject",
        cycles=["01", "02", "03"],
    )

    for cyc in ["01", "02", "03"]:
        cs = CycleStatus(cycle=cyc)
        for stage in ("stitch", "deconvolve", "edf"):
            cs.stages[stage] = "pending"
        status.cycle_statuses[cyc] = cs

    status.registration_status = "pending"
    status.signal_isolation_status = "pending"
    status.qc_statuses = {
        "stitch": "pending",
        "decon": "pending",
        "edf": "pending",
        "registration": "pending",
        "signal_isolation": "pending",
    }

    return status
