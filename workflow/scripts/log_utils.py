"""
Python logging utilities for KINTSUGI Snakemake pipeline.

Replicates the structured logging from slurm/utils.sh:
  - log_header / log_footer  → job boundary markers
  - summary_before / summary_after → data state snapshots
  - log_info / log_warn / log_error → timestamped messages
  - Directory summarization (file counts, sizes)

These are designed for headless SLURM execution where structured log
output is the primary debugging interface.
"""

import os
import shutil
import subprocess
import time
from datetime import datetime
from pathlib import Path


def log_header(step_name: str, cycle: int, project_dir: str | Path) -> None:
    """Print a structured header at the start of a processing step."""
    project_dir = Path(project_dir)
    print("=" * 60)
    print(f"KINTSUGI {step_name} - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"Job ID:      {os.environ.get('SLURM_JOB_ID', 'local')}")
    print(f"Array Task:  {cycle}")
    print(f"Project:     {project_dir}")
    print(f"Node:        {os.environ.get('SLURMD_NODENAME', _hostname())}")
    print(f"GPUs:        {os.environ.get('CUDA_VISIBLE_DEVICES', 'none')}")
    print("=" * 60)


def log_footer(exit_code: int = 0) -> None:
    """Print a structured footer at the end of a processing step."""
    print()
    print("=" * 60)
    print(f"Completed: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Exit code: {exit_code}")
    print("=" * 60)


def log_info(msg: str) -> None:
    """Print a timestamped info message."""
    print(f"[INFO] {datetime.now().strftime('%H:%M:%S')} {msg}")


def log_warn(msg: str) -> None:
    """Print a timestamped warning message."""
    print(f"[WARN] {datetime.now().strftime('%H:%M:%S')} {msg}")


def log_error(msg: str) -> None:
    """Print a timestamped error message."""
    print(f"[ERROR] {datetime.now().strftime('%H:%M:%S')} {msg}")


def summarize_directory(directory: str | Path, label: str) -> str:
    """
    Generate a text summary of a data directory.

    Returns a multi-line string with file counts and total size.
    Equivalent to the summarize_directory() function in utils.sh.
    """
    directory = Path(directory)
    lines = [f"{label}:"]

    if not directory.exists():
        lines.append("  (not found)")
        return "\n".join(lines)

    lines.append(f"  Path: {directory}")

    tiff_count = len(list(directory.rglob("*.tif"))) + len(
        list(directory.rglob("*.tiff"))
    )
    zarr_count = len(
        [d for d in directory.rglob("*.zarr") if d.is_dir()]
    )
    lines.append(f"  TIFFs: {tiff_count}")
    lines.append(f"  Zarrs: {zarr_count}")

    total_size = _dir_size_human(directory)
    lines.append(f"  Total size: {total_size}")

    return "\n".join(lines)


def summary_before(step_name: str, project_dir: str | Path) -> None:
    """
    Print a before-processing data summary.

    Mirrors summary_before() in utils.sh — shows input data state
    and GPU memory before processing starts.
    """
    project_dir = Path(project_dir)
    print()
    print("=" * 60)
    print(f"BEFORE: {step_name}")
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)

    input_dirs = {
        "stitch": project_dir / "data" / "raw",
        "decon": project_dir / "data" / "processed" / "stitched",
        "edf": project_dir / "data" / "processed" / "deconvolved",
    }

    input_dir = input_dirs.get(step_name)
    if input_dir:
        label_map = {
            "stitch": "Raw Data",
            "decon": "Stitched Data",
            "edf": "Deconvolved Data",
        }
        print(summarize_directory(input_dir, label_map.get(step_name, "Input")))

    print()
    print("GPU Memory (before):")
    print(_gpu_memory_summary())
    print()


def summary_after(
    step_name: str,
    project_dir: str | Path,
    elapsed_seconds: float,
    exit_code: int = 0,
) -> None:
    """
    Print an after-processing data summary.

    Mirrors summary_after() in utils.sh — shows output data state,
    timing, and GPU memory after processing completes.
    """
    project_dir = Path(project_dir)
    elapsed_min = int(elapsed_seconds // 60)
    elapsed_sec = int(elapsed_seconds % 60)

    print()
    print("=" * 60)
    print(f"AFTER: {step_name}")
    print(f"Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print("=" * 60)
    print(f"Exit code: {exit_code}")
    print(f"Duration: {elapsed_min}m {elapsed_sec}s")
    print()

    output_dirs = {
        "stitch": project_dir / "data" / "processed" / "stitched",
        "decon": project_dir / "data" / "processed" / "deconvolved",
        "edf": project_dir / "data" / "processed" / "edf",
    }

    output_dir = output_dirs.get(step_name)
    if output_dir:
        label_map = {
            "stitch": "Stitched Data",
            "decon": "Deconvolved Data",
            "edf": "EDF Data",
        }
        print(summarize_directory(output_dir, label_map.get(step_name, "Output")))

    print()
    print("GPU Memory (after):")
    print(_gpu_memory_summary())
    print()


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------
def _hostname() -> str:
    """Get hostname, with fallback."""
    try:
        import socket

        return socket.gethostname()
    except Exception:
        return "unknown"


def _dir_size_human(path: Path) -> str:
    """Get human-readable directory size."""
    try:
        total = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
        for unit in ["B", "KB", "MB", "GB", "TB"]:
            if total < 1024:
                return f"{total:.1f}{unit}"
            total /= 1024
        return f"{total:.1f}PB"
    except Exception:
        return "unknown"


def _gpu_memory_summary() -> str:
    """Query GPU memory via nvidia-smi."""
    try:
        result = subprocess.run(
            [
                "nvidia-smi",
                "--query-gpu=index,memory.used,memory.total",
                "--format=csv,noheader",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0 and result.stdout.strip():
            return "  " + result.stdout.strip().replace("\n", "\n  ")
        return "  N/A"
    except Exception:
        return "  nvidia-smi not available"
