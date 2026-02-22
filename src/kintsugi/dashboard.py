"""
Snakemake progress dashboard for KINTSUGI pipeline on HiPerGator.

Provides real-time tracking of:
- Per-cycle, per-stage completion status across all project locations
- SLURM job states (running, pending, failed) with node assignments
- Memory and GPU hardware utilization from sacct/squeue
- Estimated completion times based on historical log timing data

Usage:
    kintsugi workflow status .              # One-shot status
    kintsugi workflow status . --watch      # Auto-refresh dashboard
    kintsugi workflow status . --watch -i 30  # Refresh every 30s
"""

from __future__ import annotations

import os
import re
import subprocess
import time
from dataclasses import dataclass, field
from datetime import datetime, timedelta
from pathlib import Path
from typing import Any  # noqa: I001

# ---------------------------------------------------------------------------
# Data models
# ---------------------------------------------------------------------------

STAGES = ("stitch", "deconvolve", "edf", "registration")
STAGE_DIRS = {
    "stitch": "stitched",
    "deconvolve": "deconvolved",
    "edf": "edf",
    "registration": "registered",
}
QC_STAGES = ("stitch", "decon", "edf", "registration")


@dataclass
class JobInfo:
    """SLURM job metadata for a single job."""

    job_id: str
    state: str  # RUNNING, PENDING, COMPLETED, FAILED, TIMEOUT, etc.
    node: str = ""
    partition: str = ""
    account: str = ""
    elapsed: str = ""  # HH:MM:SS
    time_limit: str = ""
    mem_requested: str = ""
    mem_used: str = ""  # from sacct MaxRSS
    gpu_id: str = ""
    rule: str = ""
    cycle: str = ""


@dataclass
class CycleStatus:
    """Progress status for a single cycle across all stages."""

    cycle: str  # zero-padded, e.g. "01"
    stages: dict[str, str] = field(default_factory=dict)  # stage -> status
    jobs: dict[str, JobInfo | None] = field(default_factory=dict)  # stage -> job
    timings: dict[str, float] = field(default_factory=dict)  # stage -> seconds


@dataclass
class ProjectStatus:
    """Full progress status for a project."""

    project_dir: Path
    project_name: str
    cycles: list[str]  # zero-padded cycle strings
    cycle_statuses: dict[str, CycleStatus] = field(default_factory=dict)
    registration_status: str = "pending"
    registration_job: JobInfo | None = None
    registration_timing: float = 0.0
    qc_statuses: dict[str, str] = field(default_factory=dict)  # qc_stage -> status
    hardware: HardwareStatus | None = None
    scan_time: datetime = field(default_factory=datetime.now)


@dataclass
class HardwareStatus:
    """Live HiPerGator hardware utilization."""

    accounts: list[dict[str, Any]] = field(default_factory=list)
    total_gpu_alloc: int = 0
    total_gpu_used: int = 0
    total_gpu_avail: int = 0
    total_cpu_alloc: int = 0
    total_cpu_used: int = 0
    total_mem_alloc_gb: int = 0
    total_mem_used_gb: int = 0
    kintsugi_jobs: list[JobInfo] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Progress scanning
# ---------------------------------------------------------------------------


def scan_project_progress(
    project_dir: str | Path,
    config: dict[str, Any] | None = None,
) -> ProjectStatus:
    """Scan a project directory for pipeline progress.

    Checks sentinel files, counts output files per cycle/stage, queries
    SLURM for active jobs, and parses log files for timing data.

    Parameters
    ----------
    project_dir : str or Path
        Path to the KINTSUGI project directory.
    config : dict, optional
        Parsed workflow/config.yaml. If None, will attempt to load it.

    Returns
    -------
    ProjectStatus
        Complete progress snapshot.
    """
    project_dir = Path(project_dir).resolve()

    if config is None:
        config = _load_config(project_dir)

    try:
        cycles = [f"{int(c):02d}" for c in config.get("cycles", [])]
    except (ValueError, TypeError):
        cycles = []
    project_name = config.get("project_name", project_dir.name)

    status = ProjectStatus(
        project_dir=project_dir,
        project_name=project_name,
        cycles=cycles,
    )

    # Scan per-cycle stages
    for cyc in cycles:
        cs = CycleStatus(cycle=cyc)
        for stage in ("stitch", "deconvolve", "edf"):
            cs.stages[stage] = _check_stage_status(project_dir, stage, cyc)
        status.cycle_statuses[cyc] = cs

    # Registration (single job, all cycles)
    status.registration_status = _check_registration_status(project_dir)

    # QC reports
    for qc_stage in QC_STAGES:
        sentinel = project_dir / "qc_plots" / f".snakemake_complete_{qc_stage}"
        status.qc_statuses[qc_stage] = "done" if sentinel.exists() else "pending"

    # Parse SLURM jobs
    _attach_slurm_jobs(status)

    # Parse log timings
    _attach_log_timings(status)

    # Hardware status
    if config.get("resources", {}).get("accounts"):
        account_names = [a["name"] for a in config["resources"]["accounts"]]
        status.hardware = _scan_hardware(account_names, config.get("resources", {}))

    status.scan_time = datetime.now()
    return status


def _load_config(project_dir: Path) -> dict[str, Any]:
    """Load workflow/config.yaml from a project directory."""
    config_path = project_dir / "workflow" / "config.yaml"
    if not config_path.exists():
        return {}
    try:
        import yaml

        with open(config_path) as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _check_stage_status(project_dir: Path, stage: str, cycle: str) -> str:
    """Check completion status of a processing stage for a cycle.

    Returns one of: 'done', 'partial', 'pending'.
    """
    stage_dir_name = STAGE_DIRS.get(stage, stage)
    stage_dir = project_dir / "data" / "processed" / stage_dir_name / f"cyc{cycle}"
    sentinel = stage_dir / ".snakemake_complete"

    if sentinel.exists():
        return "done"

    # Check for partial output (some files exist but sentinel missing)
    if stage_dir.exists():
        tifs = list(stage_dir.glob("*.tif"))
        if tifs:
            return "partial"

    return "pending"


def _check_registration_status(project_dir: Path) -> str:
    """Check registration completion status."""
    reg_dir = project_dir / "data" / "processed" / "registered"
    sentinel = reg_dir / ".snakemake_complete"

    if sentinel.exists():
        return "done"

    if reg_dir.exists():
        # Check for any registered output
        cyc_dirs = [d for d in reg_dir.iterdir() if d.is_dir() and d.name.startswith("cyc")]
        if cyc_dirs:
            return "partial"

    return "pending"


# ---------------------------------------------------------------------------
# SLURM job queries
# ---------------------------------------------------------------------------

# Regex to parse Snakemake job names: smk-<rule>-<cycle> or smk-<rule>
_SMK_JOB_RE = re.compile(
    r"(?:smk-|snakejob\.|group_)"
    r"(stitch|deconvolve|decon|edf|registration|qc_\w+|vessel3d|spillover)"
    r"(?:[._-](?:cyc)?(\d+))?"
)

# Also match SLURM log file paths from --output/--error flags
_LOG_CYCLE_RE = re.compile(r"(stitch|decon|edf|registration|qc_\w+)_cyc(\d+)")


def _attach_slurm_jobs(status: ProjectStatus) -> None:
    """Query squeue for KINTSUGI jobs and attach to cycle statuses."""
    user = os.environ.get("USER", "")
    if not user:
        return

    try:
        proc = subprocess.run(
            [
                "squeue",
                "-u",
                user,
                "-h",
                "-o",
                "%i|%j|%T|%N|%P|%a|%M|%l|%m|%b",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if proc.returncode != 0:
            return
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return

    for line in proc.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.strip().split("|")
        if len(parts) < 10:
            continue

        job_id, name, state, node, partition, account, elapsed, limit, mem, gres = parts[:10]

        # Parse rule and cycle from job name
        rule, cycle = _parse_job_name(name)
        if not rule:
            continue

        job = JobInfo(
            job_id=job_id,
            state=state,
            node=node,
            partition=partition,
            account=account,
            elapsed=elapsed,
            time_limit=limit,
            mem_requested=mem,
            gpu_id=gres if gres and gres != "N/A" else "",
            rule=rule,
            cycle=cycle,
        )

        # Map rule name to canonical stage
        stage = _normalize_rule(rule)

        if stage == "registration":
            status.registration_job = job
            if state == "RUNNING":
                status.registration_status = "running"
            elif state == "PENDING":
                status.registration_status = "queued"
        elif cycle and cycle in status.cycle_statuses:
            cs = status.cycle_statuses[cycle]
            cs.jobs[stage] = job
            if state == "RUNNING" and cs.stages.get(stage) != "done":
                cs.stages[stage] = "running"
            elif state == "PENDING" and cs.stages.get(stage) == "pending":
                cs.stages[stage] = "queued"

    # Also check sacct for recently completed/failed jobs
    _attach_recent_sacct(status)


def _parse_job_name(name: str) -> tuple[str, str]:
    """Extract rule and cycle from a Snakemake SLURM job name.

    Returns (rule, cycle) where cycle is zero-padded or empty string.
    """
    m = _SMK_JOB_RE.search(name)
    if m:
        rule = m.group(1)
        cycle = f"{int(m.group(2)):02d}" if m.group(2) else ""
        return rule, cycle
    return "", ""


def _normalize_rule(rule: str) -> str:
    """Map various rule names to canonical stage names."""
    mapping = {
        "stitch": "stitch",
        "deconvolve": "deconvolve",
        "decon": "deconvolve",
        "edf": "edf",
        "registration": "registration",
    }
    return mapping.get(rule, rule)


def _attach_recent_sacct(status: ProjectStatus) -> None:
    """Query sacct for recently failed/completed jobs (last 2 hours)."""
    user = os.environ.get("USER", "")
    if not user:
        return

    start_time = (datetime.now() - timedelta(hours=2)).strftime("%Y-%m-%dT%H:%M")

    try:
        proc = subprocess.run(
            [
                "sacct",
                "-u",
                user,
                "-S",
                start_time,
                "--format=JobID,JobName%50,State,Elapsed,MaxRSS,ExitCode",
                "-n",
                "-P",
                "--noconvert",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if proc.returncode != 0:
            return
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return

    for line in proc.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.strip().split("|")
        if len(parts) < 6:
            continue

        job_id, name, state, elapsed, max_rss, exit_code = parts[:6]

        # Skip sub-steps (batch, extern)
        if "." in job_id:
            # Use MaxRSS from batch sub-step to update parent job
            parent_id = job_id.split(".")[0]
            if max_rss:
                _update_job_mem(status, parent_id, max_rss)
            continue

        rule, cycle = _parse_job_name(name)
        if not rule:
            continue

        stage = _normalize_rule(rule)

        # Only attach FAILED/TIMEOUT states (COMPLETED already tracked by sentinels)
        if state in ("FAILED", "TIMEOUT", "OUT_OF_MEMORY", "CANCELLED"):
            if stage == "registration":
                if status.registration_status not in ("done", "running"):
                    status.registration_status = "failed"
            elif cycle and cycle in status.cycle_statuses:
                cs = status.cycle_statuses[cycle]
                if cs.stages.get(stage) not in ("done", "running"):
                    cs.stages[stage] = "failed"


def _update_job_mem(status: ProjectStatus, parent_id: str, max_rss: str) -> None:
    """Update a job's actual memory usage from sacct MaxRSS."""
    # Search through all attached jobs
    if status.registration_job and status.registration_job.job_id == parent_id:
        status.registration_job.mem_used = max_rss
        return

    for cs in status.cycle_statuses.values():
        for job in cs.jobs.values():
            if job and job.job_id == parent_id:
                job.mem_used = max_rss
                return


# ---------------------------------------------------------------------------
# Log timing extraction
# ---------------------------------------------------------------------------


def _attach_log_timings(status: ProjectStatus) -> None:
    """Parse Snakemake log files for actual processing durations."""
    log_dir = status.project_dir / "slurm" / "logs" / "snakemake"
    if not log_dir.exists():
        return

    for log_file in log_dir.glob("*.log"):
        m = _LOG_CYCLE_RE.search(log_file.stem)
        if not m:
            # Check for registration (no cycle)
            if "registration" in log_file.stem:
                elapsed = _parse_log_duration(log_file)
                if elapsed > 0:
                    status.registration_timing = elapsed
            continue

        rule_name = m.group(1)
        cycle = f"{int(m.group(2)):02d}"
        stage = _normalize_rule(rule_name)

        if cycle in status.cycle_statuses:
            elapsed = _parse_log_duration(log_file)
            if elapsed > 0:
                status.cycle_statuses[cycle].timings[stage] = elapsed


def _parse_log_duration(log_path: Path) -> float:
    """Extract processing duration from a KINTSUGI log file.

    Looks for the log_utils summary_after "Duration: XmYs" pattern,
    or falls back to computing from log_header/log_footer timestamps.

    Returns elapsed seconds, or 0 if unparseable.
    """
    try:
        content = log_path.read_text(errors="replace")
    except OSError:
        return 0.0

    # Pattern 1: "Duration: 12m 34s"
    m = re.search(r"Duration:\s*(\d+)m\s*(\d+)s", content)
    if m:
        return int(m.group(1)) * 60 + int(m.group(2))

    # Pattern 2: Parse timestamps from header and footer
    timestamps = re.findall(r"(\d{4}-\d{2}-\d{2}\s+\d{2}:\d{2}:\d{2})", content)
    if len(timestamps) >= 2:
        try:
            fmt = "%Y-%m-%d %H:%M:%S"
            t_start = datetime.strptime(timestamps[0], fmt)
            t_end = datetime.strptime(timestamps[-1], fmt)
            delta = (t_end - t_start).total_seconds()
            if delta > 0:
                return delta
        except ValueError:
            pass

    return 0.0


# ---------------------------------------------------------------------------
# Hardware / resource monitoring
# ---------------------------------------------------------------------------


def _scan_hardware(
    account_names: list[str],
    resources_config: dict[str, Any],
) -> HardwareStatus:
    """Query HiPerGator for live hardware utilization."""
    from kintsugi.hpc import detect_account_resources, detect_current_usage

    hw = HardwareStatus()

    for acct_name in account_names:
        alloc = detect_account_resources(acct_name)
        usage = detect_current_usage(acct_name)

        gpu_avail = max(0, alloc["gpus"] - usage["gpus"])

        hw.accounts.append(
            {
                "name": acct_name,
                "gpu_alloc": alloc["gpus"],
                "gpu_used": usage["gpus"],
                "gpu_avail": gpu_avail,
                "cpu_alloc": alloc["cpus"],
                "cpu_used": usage["cpus"],
                "mem_alloc_gb": alloc["mem_gb"],
                "mem_used_gb": usage["mem_gb"],
            }
        )
        hw.total_gpu_alloc += alloc["gpus"]
        hw.total_gpu_used += usage["gpus"]
        hw.total_gpu_avail += gpu_avail
        hw.total_cpu_alloc += alloc["cpus"]
        hw.total_cpu_used += usage["cpus"]
        hw.total_mem_alloc_gb += alloc["mem_gb"]
        hw.total_mem_used_gb += usage["mem_gb"]

    # Gather all KINTSUGI jobs for hardware detail
    hw.kintsugi_jobs = _get_kintsugi_jobs()

    return hw


def _get_kintsugi_jobs() -> list[JobInfo]:
    """Get all KINTSUGI-related SLURM jobs for current user."""
    user = os.environ.get("USER", "")
    if not user:
        return []

    try:
        proc = subprocess.run(
            [
                "squeue",
                "-u",
                user,
                "-h",
                "-o",
                "%i|%j|%T|%N|%P|%a|%M|%l|%m|%b",
            ],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if proc.returncode != 0:
            return []
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return []

    jobs = []
    for line in proc.stdout.strip().split("\n"):
        if not line.strip():
            continue
        parts = line.strip().split("|")
        if len(parts) < 10:
            continue

        job_id, name, state, node, partition, account, elapsed, limit, mem, gres = parts[:10]

        rule, cycle = _parse_job_name(name)
        if not rule:
            continue

        jobs.append(
            JobInfo(
                job_id=job_id,
                state=state,
                node=node,
                partition=partition,
                account=account,
                elapsed=elapsed,
                time_limit=limit,
                mem_requested=mem,
                gpu_id=gres if gres and gres != "N/A" else "",
                rule=rule,
                cycle=cycle,
            )
        )

    return jobs


# ---------------------------------------------------------------------------
# Completion time estimation
# ---------------------------------------------------------------------------


def estimate_completion(status: ProjectStatus) -> dict[str, Any]:
    """Estimate remaining time to pipeline completion.

    Uses historical timing from completed stages to predict remaining
    work. Falls back to default estimates if no history is available.

    Returns
    -------
    dict
        {
            "total_remaining_seconds": float,
            "eta": datetime or None,
            "per_stage_remaining": dict[str, float],
            "confidence": str,  # "high", "medium", "low"
        }
    """
    # Collect historical timings per stage
    stage_avg: dict[str, list[float]] = {}
    for cs in status.cycle_statuses.values():
        for stage, elapsed in cs.timings.items():
            if elapsed > 0:
                stage_avg.setdefault(stage, []).append(elapsed)

    # Compute averages
    avg_times: dict[str, float] = {}
    for stage, times in stage_avg.items():
        avg_times[stage] = sum(times) / len(times)

    # Default estimates (seconds) for HiPerGator GPU
    # Based on measured timings from CX_19-003 (9x7 grid, 4ch, 20z)
    defaults = {
        "stitch": 480,  # ~8 min
        "deconvolve": 720,  # ~12 min
        "edf": 90,  # ~1.5 min
        "registration": 1800,  # ~30 min
    }

    confidence = "high" if len(avg_times) >= 2 else ("medium" if avg_times else "low")

    # Calculate remaining work
    remaining_per_stage: dict[str, float] = {}
    total_remaining = 0.0

    for _cyc, cs in status.cycle_statuses.items():
        for stage in ("stitch", "deconvolve", "edf"):
            if cs.stages.get(stage) in ("done",):
                continue
            if cs.stages.get(stage) == "running":
                # Estimate remaining portion of running job
                job = cs.jobs.get(stage)
                elapsed_sec = _parse_elapsed(job.elapsed) if job else 0
                expected = avg_times.get(stage, defaults.get(stage, 600))
                remaining = max(0, expected - elapsed_sec)
            else:
                remaining = avg_times.get(stage, defaults.get(stage, 600))

            remaining_per_stage[stage] = remaining_per_stage.get(stage, 0) + remaining
            total_remaining += remaining

    # Registration (sequential, happens once after all EDF)
    if status.registration_status not in ("done",):
        if status.registration_status == "running" and status.registration_job:
            elapsed_sec = _parse_elapsed(status.registration_job.elapsed)
            expected = avg_times.get("registration", defaults["registration"])
            remaining = max(0, expected - elapsed_sec)
        else:
            remaining = avg_times.get("registration", defaults["registration"])
        remaining_per_stage["registration"] = remaining
        total_remaining += remaining

    # Account for parallelism: GPU-only scheduling means jobs run N at a time
    # where N is total GPU slots. This is a rough estimate.
    n_running = sum(
        1
        for cs in status.cycle_statuses.values()
        for stage in ("stitch", "deconvolve", "edf")
        if cs.stages.get(stage) == "running"
    )
    parallelism = max(1, n_running) if n_running > 0 else 1

    # Rough wall-clock estimate: divide non-registration work by parallelism
    non_reg_remaining = total_remaining - remaining_per_stage.get("registration", 0)
    wall_clock = (non_reg_remaining / parallelism) + remaining_per_stage.get("registration", 0)

    eta = datetime.now() + timedelta(seconds=wall_clock) if wall_clock > 0 else None

    return {
        "total_remaining_seconds": total_remaining,
        "wall_clock_seconds": wall_clock,
        "eta": eta,
        "per_stage_remaining": remaining_per_stage,
        "confidence": confidence,
        "parallelism": parallelism,
    }


def _parse_elapsed(elapsed_str: str) -> float:
    """Parse SLURM elapsed time string to seconds.

    Handles formats: HH:MM:SS, D-HH:MM:SS, MM:SS
    """
    if not elapsed_str:
        return 0.0

    try:
        days = 0
        if "-" in elapsed_str:
            day_part, elapsed_str = elapsed_str.split("-", 1)
            days = int(day_part)

        parts = elapsed_str.split(":")
        if len(parts) == 3:
            h, m, s = int(parts[0]), int(parts[1]), int(parts[2])
            return days * 86400 + h * 3600 + m * 60 + s
        elif len(parts) == 2:
            m, s = int(parts[0]), int(parts[1])
            return days * 86400 + m * 60 + s
    except ValueError:
        pass
    return 0.0


# ---------------------------------------------------------------------------
# Multi-project scanning
# ---------------------------------------------------------------------------


def scan_all_projects(
    base_dir: str | Path,
    project_dirs: list[str | Path] | None = None,
) -> list[ProjectStatus]:
    """Scan multiple project directories for pipeline progress.

    Parameters
    ----------
    base_dir : str or Path
        Base directory to search for projects (if project_dirs is None).
    project_dirs : list, optional
        Explicit list of project directories. If None, auto-discovers
        projects under base_dir by looking for workflow/config.yaml.

    Returns
    -------
    list[ProjectStatus]
        Status for each discovered project.
    """
    base_dir = Path(base_dir).resolve()

    if project_dirs is None:
        project_dirs = []
        # Auto-discover: look for workflow/config.yaml
        for config_file in base_dir.rglob("workflow/config.yaml"):
            project_dirs.append(config_file.parent.parent)
        # Also check if base_dir itself is a project
        if (base_dir / "workflow" / "config.yaml").exists():
            if base_dir not in project_dirs:
                project_dirs = [base_dir] + list(project_dirs)

    results = []
    for pdir in project_dirs:
        pdir = Path(pdir).resolve()
        try:
            results.append(scan_project_progress(pdir))
        except Exception:
            # Skip projects that fail to scan
            pass

    return results


# ---------------------------------------------------------------------------
# Rich terminal rendering
# ---------------------------------------------------------------------------

# Status symbols and colors
_STATUS_STYLE = {
    "done": ("[bold green]", "OK"),
    "running": ("[bold cyan]", ">>"),
    "queued": ("[yellow]", ".."),
    "partial": ("[yellow]", "~~"),
    "pending": ("[dim]", "--"),
    "failed": ("[bold red]", "XX"),
}


def render_dashboard(
    statuses: list[ProjectStatus],
    show_hardware: bool = True,
    show_estimates: bool = True,
    show_jobs: bool = True,
) -> None:
    """Render the progress dashboard to the terminal using Rich.

    Parameters
    ----------
    statuses : list[ProjectStatus]
        One or more project status snapshots to display.
    show_hardware : bool
        Include hardware utilization section.
    show_estimates : bool
        Include completion time estimates.
    show_jobs : bool
        Include active SLURM job details.
    """
    from rich.console import Console
    from rich.panel import Panel

    console = Console()

    # Header
    now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    console.print()
    console.print(
        Panel(
            f"[bold]KINTSUGI Pipeline Dashboard[/bold]  |  {now}",
            style="blue",
            expand=True,
        )
    )

    for status in statuses:
        _render_project(console, status, show_estimates, show_jobs)

    # Hardware section (shared across projects — uses first project's hardware data)
    if show_hardware:
        hw = next((s.hardware for s in statuses if s.hardware), None)
        if hw:
            _render_hardware(console, hw)

    console.print()


def _render_project(
    console: Any,
    status: ProjectStatus,
    show_estimates: bool,
    show_jobs: bool,
) -> None:
    """Render a single project's progress table."""
    from rich.table import Table
    from rich.text import Text

    # Summary line
    n_cycles = len(status.cycles)
    done_count = sum(
        1
        for cs in status.cycle_statuses.values()
        if all(cs.stages.get(s) == "done" for s in ("stitch", "deconvolve", "edf"))
    )
    console.print()
    console.print(
        f"[bold]{status.project_name}[/bold]  "
        f"({status.project_dir})  "
        f"[dim]{done_count}/{n_cycles} cycles complete[/dim]"
    )

    # Progress table: rows = cycles, columns = stages
    table = Table(
        show_header=True,
        header_style="bold",
        show_lines=False,
        pad_edge=False,
        expand=False,
    )
    table.add_column("Cycle", style="bold", width=6)
    table.add_column("Stitch", width=10, justify="center")
    table.add_column("Decon", width=10, justify="center")
    table.add_column("EDF", width=10, justify="center")
    table.add_column("Timing", width=16, justify="right")
    if show_jobs:
        table.add_column("Job / Node", width=28)

    for cyc in status.cycles:
        cs = status.cycle_statuses.get(cyc)
        if not cs:
            continue

        row = [f"cyc{cyc}"]

        for stage in ("stitch", "deconvolve", "edf"):
            st = cs.stages.get(stage, "pending")
            style, symbol = _STATUS_STYLE.get(st, ("[dim]", "??"))
            cell = Text(f" {symbol} ", style=style.strip("[]"))
            row.append(cell)

        # Timing column: show cumulative time for completed stages
        timing_parts = []
        for stage in ("stitch", "deconvolve", "edf"):
            t = cs.timings.get(stage, 0)
            if t > 0:
                timing_parts.append(f"{t / 60:.0f}m")
        timing_str = " + ".join(timing_parts) if timing_parts else ""
        if timing_parts:
            total_t = sum(cs.timings.get(s, 0) for s in ("stitch", "deconvolve", "edf"))
            timing_str += f" = {total_t / 60:.0f}m"
        row.append(timing_str)

        if show_jobs:
            # Show the most relevant active job
            job_str = ""
            for stage in ("stitch", "deconvolve", "edf"):
                job = cs.jobs.get(stage)
                if job and job.state in ("RUNNING", "PENDING"):
                    mem_info = ""
                    if job.mem_used:
                        mem_info = f" [{_format_mem(job.mem_used)}]"
                    job_str = (
                        f"{job.job_id} {job.state[:3]} "
                        f"{job.node or '(queued)'}"
                        f" {job.elapsed}{mem_info}"
                    )
                    break
            row.append(job_str)

        table.add_row(*row)

    # Registration row
    reg_style, reg_symbol = _STATUS_STYLE.get(status.registration_status, ("[dim]", "??"))
    reg_timing = f"{status.registration_timing / 60:.0f}m" if status.registration_timing > 0 else ""
    reg_job_str = ""
    if show_jobs and status.registration_job:
        rj = status.registration_job
        if rj.state in ("RUNNING", "PENDING"):
            mem_info = ""
            if rj.mem_used:
                mem_info = f" [{_format_mem(rj.mem_used)}]"
            reg_job_str = (
                f"{rj.job_id} {rj.state[:3]} {rj.node or '(queued)'} {rj.elapsed}{mem_info}"
            )

    reg_cell = Text(f" {reg_symbol} ", style=reg_style.strip("[]"))
    reg_row = ["REG", reg_cell, "", "", reg_timing]
    if show_jobs:
        reg_row.append(reg_job_str)
    table.add_row(*reg_row)

    console.print(table)

    # QC status summary
    qc_parts = []
    for qc_stage in QC_STAGES:
        st = status.qc_statuses.get(qc_stage, "pending")
        style, symbol = _STATUS_STYLE.get(st, ("[dim]", "??"))
        qc_parts.append(f"{style}{symbol}[/] {qc_stage}")
    console.print(f"  QC: {' | '.join(qc_parts)}")

    # Completion estimate
    if show_estimates:
        est = estimate_completion(status)
        if est["wall_clock_seconds"] > 0:
            wall = est["wall_clock_seconds"]
            hours = int(wall // 3600)
            mins = int((wall % 3600) // 60)
            eta = est["eta"]
            eta_str = eta.strftime("%H:%M") if eta else "N/A"
            confidence = est["confidence"]
            par = est["parallelism"]

            console.print(
                f"  ETA: ~{hours}h{mins:02d}m remaining "
                f"(ETA {eta_str}, {par}x parallel, "
                f"confidence: {confidence})"
            )
        else:
            console.print("  [green]Pipeline complete[/green]")


def _render_hardware(console: Any, hw: HardwareStatus) -> None:
    """Render hardware utilization section."""
    from rich.table import Table

    console.print()
    console.print("[bold]HiPerGator Resources[/bold]")

    table = Table(
        show_header=True,
        header_style="bold",
        show_lines=False,
        pad_edge=False,
        expand=False,
    )
    table.add_column("Account", style="bold", width=12)
    table.add_column("GPUs", width=14, justify="center")
    table.add_column("CPUs", width=14, justify="center")
    table.add_column("Memory", width=18, justify="center")

    for acct in hw.accounts:
        gpu_str = f"{acct['gpu_used']}/{acct['gpu_alloc']}"
        cpu_str = f"{acct['cpu_used']}/{acct['cpu_alloc']}"
        mem_str = f"{acct['mem_used_gb']}/{acct['mem_alloc_gb']} GB"

        # Color code utilization
        gpu_pct = acct["gpu_used"] / acct["gpu_alloc"] if acct["gpu_alloc"] > 0 else 0
        gpu_style = "green" if gpu_pct < 0.5 else ("yellow" if gpu_pct < 0.9 else "red")

        table.add_row(
            acct["name"],
            f"[{gpu_style}]{gpu_str}[/{gpu_style}]",
            cpu_str,
            mem_str,
        )

    # Totals row
    total_gpu_str = f"{hw.total_gpu_used}/{hw.total_gpu_alloc}"
    total_cpu_str = f"{hw.total_cpu_used}/{hw.total_cpu_alloc}"
    total_mem_str = f"{hw.total_mem_used_gb}/{hw.total_mem_alloc_gb} GB"
    table.add_row(
        "[bold]Total[/bold]",
        f"[bold]{total_gpu_str}[/bold]",
        f"[bold]{total_cpu_str}[/bold]",
        f"[bold]{total_mem_str}[/bold]",
    )

    console.print(table)

    # Active jobs detail
    if hw.kintsugi_jobs:
        console.print()
        console.print(f"  Active KINTSUGI jobs: {len(hw.kintsugi_jobs)}")

        job_table = Table(
            show_header=True,
            header_style="bold dim",
            show_lines=False,
            pad_edge=False,
            expand=False,
        )
        job_table.add_column("JobID", width=10)
        job_table.add_column("Rule", width=12)
        job_table.add_column("Cycle", width=6)
        job_table.add_column("State", width=8)
        job_table.add_column("Node", width=14)
        job_table.add_column("Elapsed", width=10)
        job_table.add_column("Mem", width=10)
        job_table.add_column("GPU", width=10)

        for job in hw.kintsugi_jobs:
            state_color = {
                "RUNNING": "cyan",
                "PENDING": "yellow",
                "COMPLETING": "green",
            }.get(job.state, "dim")

            job_table.add_row(
                job.job_id,
                job.rule,
                job.cycle or "-",
                f"[{state_color}]{job.state[:7]}[/{state_color}]",
                job.node or "(queued)",
                job.elapsed,
                _format_mem(job.mem_requested),
                job.gpu_id or "-",
            )

        console.print(job_table)


def _format_mem(mem_str: str) -> str:
    """Format a memory string for display.

    Handles SLURM formats like '4000M', '48G', '64000', raw bytes, etc.
    """
    if not mem_str:
        return ""

    mem_str = mem_str.strip()

    # Already human-readable with unit
    m = re.match(r"(\d+(?:\.\d+)?)\s*([KMGT]?)i?[Bb]?", mem_str, re.IGNORECASE)
    if m:
        val = float(m.group(1))
        unit = m.group(2).upper()
        if unit == "K":
            return f"{val / 1024:.0f}M" if val >= 1024 else f"{val:.0f}K"
        elif unit == "M":
            return f"{val / 1024:.1f}G" if val >= 1024 else f"{val:.0f}M"
        elif unit == "G":
            return f"{val:.1f}G"
        elif unit == "T":
            return f"{val:.1f}T"
        else:
            # Raw number — assume MB (SLURM default)
            if val >= 1024:
                return f"{val / 1024:.1f}G"
            return f"{val:.0f}M"

    return mem_str


# ---------------------------------------------------------------------------
# Watch mode (auto-refresh)
# ---------------------------------------------------------------------------


def watch_dashboard(
    project_dirs: list[str | Path],
    interval: int = 30,
    show_hardware: bool = True,
    show_estimates: bool = True,
    show_jobs: bool = True,
) -> None:
    """Run the dashboard in watch mode with auto-refresh.

    Clears the terminal and re-renders at each interval.
    Press Ctrl+C to exit.

    Parameters
    ----------
    project_dirs : list
        Project directories to monitor.
    interval : int
        Refresh interval in seconds (default: 30).
    show_hardware : bool
        Include hardware utilization section.
    show_estimates : bool
        Include completion time estimates.
    show_jobs : bool
        Include active SLURM job details.
    """
    import sys

    from rich.console import Console

    console = Console()

    try:
        while True:
            console.clear()

            statuses = []
            for pdir in project_dirs:
                try:
                    statuses.append(scan_project_progress(pdir))
                except Exception as e:
                    console.print(f"[red]Error scanning {pdir}: {e}[/red]")

            if statuses:
                render_dashboard(
                    statuses,
                    show_hardware=show_hardware,
                    show_estimates=show_estimates,
                    show_jobs=show_jobs,
                )

            console.print(f"\n[dim]Refreshing every {interval}s. Press Ctrl+C to exit.[/dim]")

            time.sleep(interval)
    except KeyboardInterrupt:
        print("\nDashboard stopped.")
        sys.exit(0)
