"""
HPC (High-Performance Computing) detection and configuration utilities.

Provides auto-detection of SLURM cluster settings from the environment,
including account names, partitions, GPU types, and conda environments.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path
from typing import Any

# ============================================================================
# Account blocklist — these accounts are NEVER used for any KINTSUGI jobs.
# Enforced in get_user_slurm_accounts(), detect_multi_account_resources(),
# and detect_live_multi_account().  To add/remove accounts, edit this set.
#
# - brusko: permanently blocked (legacy)
# - clive:  blocked Apr 2026 — QOS pool throttled to 1 GPU / 312.5 GB and
#           regularly saturated by other group members; KINTSUGI jobs end up
#           stuck on QOSGrpMemLimit indefinitely. Use maigan exclusively.
# ============================================================================
BLOCKED_ACCOUNTS: frozenset[str] = frozenset({"brusko", "clive"})


def _filter_blocked(accounts: list[str]) -> list[str]:
    """Remove blocked accounts from a list. Single enforcement point."""
    return [a for a in accounts if a not in BLOCKED_ACCOUNTS]


def detect_hpc_settings() -> dict[str, Any]:
    """
    Auto-detect HPC cluster settings from environment.

    Detects:
    - SLURM account from SLURM_ACCOUNT or SBATCH_ACCOUNT environment variables
    - Partition from SLURM_PARTITION environment variable
    - GPU type from nvidia-smi query
    - Conda environment from CONDA_DEFAULT_ENV
    - QOS from environment or sane defaults

    Returns
    -------
    dict
        Dictionary with detected settings:
        - account: str or None
        - partition: str or None
        - gpu_type: str or None
        - conda_env: str or None
        - qos: str or None
    """
    settings: dict[str, Any] = {}

    # SLURM account - check multiple environment variables
    for env_var in ["SLURM_ACCOUNT", "SBATCH_ACCOUNT", "SLURM_JOB_ACCOUNT"]:
        if env_var in os.environ:
            settings["account"] = os.environ[env_var]
            break

    # Partition
    if "SLURM_PARTITION" in os.environ:
        settings["partition"] = os.environ["SLURM_PARTITION"]

    # QOS (Quality of Service)
    if "SLURM_QOS" in os.environ:
        settings["qos"] = os.environ["SLURM_QOS"]

    # GPU type from nvidia-smi
    gpu_type = _detect_gpu_type()
    if gpu_type:
        settings["gpu_type"] = gpu_type

    # Conda environment
    if "CONDA_DEFAULT_ENV" in os.environ:
        settings["conda_env"] = os.environ["CONDA_DEFAULT_ENV"]

    return settings


def _detect_gpu_type() -> str | None:
    """
    Detect GPU type using nvidia-smi.

    Returns
    -------
    str or None
        Normalized GPU type name (e.g., 'a100', 'v100', 'b200') or None if detection fails
    """
    try:
        result = subprocess.run(
            ["nvidia-smi", "--query-gpu=name", "--format=csv,noheader"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            gpu_name = result.stdout.strip().split("\n")[0]
            return _normalize_gpu_type(gpu_name)
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    return None


def _normalize_gpu_type(gpu_name: str) -> str:
    """
    Normalize GPU name to SLURM constraint format.

    Parameters
    ----------
    gpu_name : str
        Full GPU name from nvidia-smi (e.g., "NVIDIA A100-SXM4-80GB")

    Returns
    -------
    str
        Normalized GPU type for SLURM (e.g., "a100")
    """
    gpu_name_lower = gpu_name.lower()

    # Map common GPU names to SLURM constraint names
    gpu_mappings = {
        "a100": "a100",
        "a40": "a40",
        "v100": "v100",
        "h100": "h100",
        "b200": "b200",
        "turin": "geforce",  # hpg-turin partition uses GeForce GPUs
        "rtx 3090": "rtx3090",
        "rtx 4090": "rtx4090",
        "rtx 6000": "rtx6000",
        "quadro": "quadro",
        "titan": "titan",
        "geforce": "geforce",
    }

    for pattern, normalized in gpu_mappings.items():
        if pattern in gpu_name_lower:
            return normalized

    # Fallback: extract alphanumeric model name
    import re

    match = re.search(r"([a-z]\d+)", gpu_name_lower)
    if match:
        return match.group(1)

    return "gpu"


def get_slurm_defaults() -> dict[str, Any]:
    """
    Get sensible default SLURM settings for KINTSUGI workloads.

    Returns
    -------
    dict
        Default SLURM configuration values
    """
    return {
        # HPC resources (HiPerGator RHEL9+)
        # GPU Partition Selection:
        #   - hpg-b200:  NVIDIA B200 GPUs (180GB VRAM, 8/node)
        #   - hpg-turin: NVIDIA L4 GPUs (24GB VRAM, 3/node)
        "partition": "hpg-b200",
        "qos": "",
        "account": "",
        "gpus_per_node": 2,
        # Fallback partition settings (for high-demand partitions like hpg-b200)
        "partition_fallback": "hpg-turin",
        "qos_fallback": "",
        # Account chain settings (multi-account with GPU/CPU fallback)
        "account_chain": "",  # Comma-separated: GPU accounts first, then burst
        "cpu_only_accounts": "",  # Accounts that cannot request GPUs
        "partition_cpu": "hpg-default",  # Partition for CPU-only accounts
        # CPU resource adjustments (for burst accounts)
        "cpu_time_multiplier": 5,  # CPU processing is slower
        "cpu_cpus_per_task": 8,
        "cpu_mem_correction": 16,
        "cpu_mem_stitch": 48,
        "cpu_mem_decon": 128,
        "cpu_mem_edf": 96,
        # Memory per job (GB) - GPU mode (CuPy processes in GPU memory)
        "mem_correction": 16,
        "mem_stitch": 48,
        "mem_decon": 48,
        "mem_edf": 48,
        # Time limits
        "time_correction": "02:00:00",
        "time_stitch": "04:00:00",
        "time_decon": "04:00:00",
        "time_edf": "01:00:00",
        # Processing parameters
        "start_cycle": 1,
        "end_cycle": 1,
        "start_channel": 1,
        "end_channel": 4,
        "output_format": "zarr",
        # Tile grid
        "tile_rows": 5,
        "tile_cols": 5,
        "tile_overlap": 0.3,
        # Microscope parameters
        "xy_vox": 377,
        "z_vox": 1500,
        "mic_na": 0.75,
        "tissue_ri": 1.44,
        "wavelengths": "1:358:461,2:753:775,3:560:575,4:648:668",
        # Conda
        "conda_env": "kintsugi",
        # Email
        "email": "",
        "mail_type": "END,FAIL",
    }


def is_slurm_available() -> bool:
    """
    Check if SLURM is available on the system.

    Returns
    -------
    bool
        True if sbatch command is available
    """
    try:
        result = subprocess.run(
            ["which", "sbatch"],
            capture_output=True,
            timeout=5,
        )
        return result.returncode == 0
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return False


def get_user_slurm_accounts() -> list[str]:
    """
    Get list of SLURM accounts the current user has access to.

    Returns
    -------
    list
        List of account names, empty if detection fails
    """
    try:
        result = subprocess.run(
            [
                "sacctmgr",
                "show",
                "associations",
                f"user={os.environ.get('USER', '')}",
                "format=account",
                "-n",
                "-P",
            ],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            accounts = list(
                dict.fromkeys(
                    line.strip() for line in result.stdout.strip().split("\n") if line.strip()
                )
            )
            return _filter_blocked(accounts)
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    return []


def get_account_chain_from_user() -> str:
    """
    Get comma-separated account chain, with GPU accounts first.

    Orders accounts with GPU-enabled accounts first, followed by
    CPU-only burst accounts (accounts ending in '-b').

    Returns
    -------
    str
        Comma-separated account chain, or empty string if detection fails
    """
    accounts = get_user_slurm_accounts()
    if accounts:
        # GPU accounts first (no -b suffix), then burst accounts (-b suffix)
        gpu_accounts = [a for a in accounts if not a.endswith("-b")]
        burst_accounts = [a for a in accounts if a.endswith("-b")]
        return ",".join(gpu_accounts + burst_accounts)
    return ""


def get_cpu_only_accounts_from_user() -> str:
    """
    Get comma-separated list of CPU-only accounts.

    Detects burst accounts (ending in '-b') which typically don't have GPU access.

    Returns
    -------
    str
        Comma-separated list of CPU-only accounts
    """
    accounts = get_user_slurm_accounts()
    if accounts:
        burst_accounts = [a for a in accounts if a.endswith("-b")]
        return ",".join(burst_accounts)
    return ""


def detect_account_resources(account: str) -> dict[str, int]:
    """
    Detect SLURM account resource limits via sacctmgr.

    Queries the account's QOS GrpTRES limits to determine available
    CPUs, GPUs, and memory.

    Parameters
    ----------
    account : str
        SLURM account name (e.g., 'maigan', 'clive')

    Returns
    -------
    dict
        {"cpus": int, "gpus": int, "mem_gb": int}
        Falls back to {0, 0, 0} if sacctmgr is unavailable.
    """
    result = {"cpus": 0, "gpus": 0, "mem_gb": 0}
    if not account:
        return result

    try:
        proc = subprocess.run(
            ["sacctmgr", "show", "qos", account, "format=GrpTRES%80", "-n", "-P"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if proc.returncode != 0 or not proc.stdout.strip():
            return result

        # Parse GrpTRES string like "cpu=104,gres/gpu=3,mem=812G"
        tres_str = proc.stdout.strip().split("|")[0].strip()
        if not tres_str:
            return result

        import re

        for part in tres_str.split(","):
            part = part.strip()
            if part.startswith("cpu="):
                result["cpus"] = int(part.split("=", 1)[1])
            elif part.startswith("gres/gpu="):
                result["gpus"] = int(part.split("=", 1)[1])
            elif part.startswith("mem="):
                mem_str = part.split("=", 1)[1]
                m = re.match(r"(\d+)\s*([GMgm]?)", mem_str)
                if m:
                    val = int(m.group(1))
                    unit = m.group(2).upper()
                    if unit == "M":
                        result["mem_gb"] = val // 1024
                    else:
                        result["mem_gb"] = val

    except (subprocess.SubprocessError, FileNotFoundError, OSError, ValueError):
        pass

    return result


def detect_current_usage(account: str) -> dict[str, int]:
    """
    Detect current resource usage on a SLURM account.

    Uses slurmInfo (HiPerGator) to get investment QOS usage,
    falling back to squeue if unavailable.

    Parameters
    ----------
    account : str
        SLURM account name

    Returns
    -------
    dict
        {"cpus": int, "gpus": int, "mem_gb": int}
        Includes running + pending jobs from all users on the account.
    """
    usage = {"cpus": 0, "gpus": 0, "mem_gb": 0}
    if not account:
        return usage

    import re

    # Try slurmInfo first (HiPerGator tool, provides clean investment usage)
    try:
        proc = subprocess.run(
            ["slurmInfo", "-g", account],
            capture_output=True,
            text=True,
            timeout=15,
        )
        if proc.returncode == 0:
            # Parse "Investment (<account>):  CPU  MEM(GB)  CPU  MEM(GB)  CPU  MEM(GB)"
            # The "Total" columns (last pair) include running + pending
            for line in proc.stdout.split("\n"):
                # Match investment line: "Investment (maigan):     4       32     0        0     4       32"
                m = re.match(
                    rf"\s*Investment\s+\({re.escape(account)}\):"
                    r"\s+(\d+)\s+(\d+)"  # running: cpu, mem
                    r"\s+(\d+)\s+(\d+)"  # pending: cpu, mem
                    r"\s+(\d+)\s+(\d+)",  # total:   cpu, mem
                    line,
                )
                if m:
                    usage["cpus"] = int(m.group(5))  # total CPUs
                    usage["mem_gb"] = int(m.group(6))  # total MEM(GB)
                    break
            # slurmInfo doesn't show GPU usage per account, estimate from squeue
            usage["gpus"] = _count_gpu_usage_squeue(account)
            return usage
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass

    # Fallback: parse squeue directly
    try:
        proc = subprocess.run(
            ["squeue", "-A", account, "--state=RUNNING,PENDING", "-h", "-o", "%C|%m|%b"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if proc.returncode != 0:
            return usage

        for line in proc.stdout.strip().split("\n"):
            if not line.strip():
                continue
            parts = line.strip().split("|")
            if len(parts) < 3:
                continue

            try:
                usage["cpus"] += int(parts[0])
            except ValueError:
                pass

            mem_str = parts[1].strip()
            m = re.match(r"(\d+)\s*([GMgm]?)", mem_str)
            if m:
                val = int(m.group(1))
                unit = m.group(2).upper()
                if unit == "M":
                    usage["mem_gb"] += val // 1024
                else:
                    usage["mem_gb"] += val

            gres_str = parts[2].strip()
            if gres_str and gres_str != "N/A":
                gm = re.search(r"gpu(?::[^:]+)?:(\d+)", gres_str)
                if gm:
                    usage["gpus"] += int(gm.group(1))

    except (subprocess.SubprocessError, FileNotFoundError, OSError, ValueError):
        pass

    return usage


def _count_gpu_usage_squeue(account: str) -> int:
    """Count GPUs in use on an account via squeue GRES field."""
    try:
        proc = subprocess.run(
            ["squeue", "-A", account, "--state=RUNNING,PENDING", "-h", "-o", "%b"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if proc.returncode != 0:
            return 0

        import re

        total = 0
        for line in proc.stdout.strip().split("\n"):
            line = line.strip()
            if line and line != "N/A":
                m = re.search(r"gpu(?::[^:]+)?:(\d+)", line)
                if m:
                    total += int(m.group(1))
        return total
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        return 0


def _compute_slots(
    alloc: dict, usage: dict, cpus_per_job: int, mem_per_job: int, has_gpus: bool
) -> int:
    """Compute available job slots from allocation minus current usage."""
    avail_cpus = max(0, alloc["cpus"] - usage["cpus"])
    avail_mem = max(0, alloc["mem_gb"] - usage["mem_gb"])

    if has_gpus:
        avail_gpus = max(0, alloc["gpus"] - usage["gpus"])
        if avail_gpus == 0:
            return 0
        slots = [avail_gpus]
        if cpus_per_job > 0:
            slots.append(avail_cpus // cpus_per_job)
        if mem_per_job > 0:
            slots.append(avail_mem // mem_per_job)
        return min(slots)
    else:
        if cpus_per_job <= 0:
            return 0
        by_cpus = avail_cpus // cpus_per_job
        by_mem = avail_mem // mem_per_job if mem_per_job > 0 else by_cpus
        return min(by_cpus, by_mem)


def detect_dual_pool_resources(
    gpu_account: str,
    cpu_account: str,
    resources_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Detect combined GPU+CPU pool resources from two SLURM accounts.

    Reports both **max** slots (from QOS allocation) and **available** slots
    (allocation minus current usage by all users on the account).

    Parameters
    ----------
    gpu_account : str
        Account with GPU access (e.g., 'maigan')
    cpu_account : str
        Account for CPU-only jobs (e.g., 'clive')
    resources_config : dict, optional
        Resource config with per-job requirements:
        - cpus_per_gpu_job (default 4)
        - mem_per_gpu_job (default 48 GB)
        - cpus_per_cpu_job (default 8)
        - mem_per_cpu_job (default 48 GB)

    Returns
    -------
    dict
        {
            "gpu_slots": int,       # max GPU slots from allocation
            "cpu_slots": int,       # max CPU slots from allocation
            "total_slots": int,     # gpu_slots + cpu_slots (max)
            "gpu_avail": int,       # currently available GPU slots
            "cpu_avail": int,       # currently available CPU slots
            "total_avail": int,     # currently available total
            "gpu_alloc": {...},     # raw QOS allocation
            "cpu_alloc": {...},     # raw QOS allocation
            "gpu_usage": {...},     # current usage on GPU account
            "cpu_usage": {...},     # current usage on CPU account
        }
    """
    rc = resources_config or {}
    cpus_per_gpu = rc.get("cpus_per_gpu_job", 4)
    mem_per_gpu = rc.get("mem_per_gpu_job", 48)
    cpus_per_cpu = rc.get("cpus_per_cpu_job", 8)
    mem_per_cpu = rc.get("mem_per_cpu_job", 48)

    gpu_alloc = detect_account_resources(gpu_account)
    has_cpu_pool = cpu_account and cpu_account != gpu_account
    cpu_alloc = (
        detect_account_resources(cpu_account)
        if has_cpu_pool
        else {"cpus": 0, "gpus": 0, "mem_gb": 0}
    )

    gpu_usage = detect_current_usage(gpu_account)
    cpu_usage = (
        detect_current_usage(cpu_account) if has_cpu_pool else {"cpus": 0, "gpus": 0, "mem_gb": 0}
    )

    zero_usage = {"cpus": 0, "gpus": 0, "mem_gb": 0}

    # Max slots (from full allocation, ignoring current usage)
    gpu_slots_max = _compute_slots(gpu_alloc, zero_usage, cpus_per_gpu, mem_per_gpu, has_gpus=True)
    cpu_slots_max = _compute_slots(cpu_alloc, zero_usage, cpus_per_cpu, mem_per_cpu, has_gpus=False)

    # Available slots (allocation minus current usage)
    gpu_slots_avail = _compute_slots(gpu_alloc, gpu_usage, cpus_per_gpu, mem_per_gpu, has_gpus=True)
    cpu_slots_avail = _compute_slots(
        cpu_alloc, cpu_usage, cpus_per_cpu, mem_per_cpu, has_gpus=False
    )

    return {
        "gpu_slots": gpu_slots_max,
        "cpu_slots": cpu_slots_max,
        "total_slots": gpu_slots_max + cpu_slots_max,
        "gpu_avail": gpu_slots_avail,
        "cpu_avail": cpu_slots_avail,
        "total_avail": gpu_slots_avail + cpu_slots_avail,
        "gpu_alloc": gpu_alloc,
        "cpu_alloc": cpu_alloc,
        "gpu_usage": gpu_usage,
        "cpu_usage": cpu_usage,
    }


def detect_multi_account_resources(
    accounts: list[str],
    resources_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Detect resources across multiple SLURM accounts.

    For each account, computes:
    - GPU slots: direct from QOS gres/gpu allocation
    - CPU slots: floor(utilization_cap * account_cpus / cpus_per_cpu_job)

    GPU and CPU partitions are separate resource pools on HiPerGator,
    so GPU jobs don't consume from the CPU allocation.

    Parameters
    ----------
    accounts : list[str]
        SLURM account names (e.g., ['maigan', 'clive'])
    resources_config : dict, optional
        Resource config with:
        - cpus_per_cpu_job (default 8)
        - cpu_utilization_cap (default 0.85)
        - partition_gpu (default "hpg-b200,hpg-turin")
        - partition_cpu (default "hpg-default")

    Returns
    -------
    dict
        {
            "accounts": [
                {"name": "maigan", "gpu_slots": 2, "cpu_slots": 8, "alloc": {...}},
                {"name": "clive",  "gpu_slots": 3, "cpu_slots": 11, "alloc": {...}},
            ],
            "total_gpu_slots": 5,
            "total_cpu_slots": 19,
            "total_slots": 24,
        }
    """
    import math

    rc = resources_config or {}
    cpus_per_cpu = rc.get("cpus_per_cpu_job", 8)
    util_cap = rc.get("cpu_utilization_cap", 0.85)
    partition_gpu = rc.get("partition_gpu", "hpg-b200,hpg-turin")
    partition_cpu = rc.get("partition_cpu", "hpg-default")

    # Safety net: always filter blocked accounts even if caller forgot
    accounts = _filter_blocked(accounts)

    account_details = []
    total_gpu = 0
    total_cpu = 0

    for acct_name in accounts:
        if not acct_name or acct_name.endswith("-b"):
            continue

        alloc = detect_account_resources(acct_name)
        gpu_slots = alloc["gpus"]
        cpu_slots = math.floor(util_cap * alloc["cpus"] / cpus_per_cpu) if cpus_per_cpu > 0 else 0

        account_details.append(
            {
                "name": acct_name,
                "gpu_slots": gpu_slots,
                "cpu_slots": cpu_slots,
                "partition_gpu": partition_gpu,
                "partition_cpu": partition_cpu,
                "alloc": alloc,
            }
        )
        total_gpu += gpu_slots
        total_cpu += cpu_slots

    return {
        "accounts": account_details,
        "total_gpu_slots": total_gpu,
        "total_cpu_slots": total_cpu,
        "total_slots": total_gpu + total_cpu,
    }


def detect_live_multi_account(
    accounts: list[str],
    resources_config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Detect resources with live usage across multiple SLURM accounts.

    Like detect_multi_account_resources() but also queries current usage
    to compute available slots.

    Parameters
    ----------
    accounts : list[str]
        SLURM account names
    resources_config : dict, optional
        Same as detect_multi_account_resources()

    Returns
    -------
    dict
        Same as detect_multi_account_resources() plus:
        - Each account entry has "usage", "gpu_avail", "cpu_avail"
        - Top-level "total_gpu_avail", "total_cpu_avail", "total_avail"
    """
    import math

    rc = resources_config or {}
    cpus_per_cpu = rc.get("cpus_per_cpu_job", 8)
    mem_per_cpu = rc.get("mem_per_cpu_job", 48)
    util_cap = rc.get("cpu_utilization_cap", 0.85)
    partition_gpu = rc.get("partition_gpu", "hpg-b200,hpg-turin")
    partition_cpu = rc.get("partition_cpu", "hpg-default")

    # Safety net: always filter blocked accounts even if caller forgot
    accounts = _filter_blocked(accounts)

    account_details = []
    total_gpu = 0
    total_cpu = 0
    total_gpu_avail = 0
    total_cpu_avail = 0

    for acct_name in accounts:
        if not acct_name or acct_name.endswith("-b"):
            continue

        alloc = detect_account_resources(acct_name)
        usage = detect_current_usage(acct_name)

        gpu_slots = alloc["gpus"]
        cpu_slots = math.floor(util_cap * alloc["cpus"] / cpus_per_cpu) if cpus_per_cpu > 0 else 0

        # Available GPU slots
        gpu_avail = max(0, alloc["gpus"] - usage["gpus"])

        # Available CPU slots (from remaining capacity after usage, with util cap)
        avail_cpus = max(0, math.floor(util_cap * alloc["cpus"]) - usage["cpus"])
        avail_mem = max(0, alloc["mem_gb"] - usage["mem_gb"])
        cpu_avail_by_cpus = avail_cpus // cpus_per_cpu if cpus_per_cpu > 0 else 0
        cpu_avail_by_mem = avail_mem // mem_per_cpu if mem_per_cpu > 0 else cpu_avail_by_cpus
        cpu_avail = min(cpu_avail_by_cpus, cpu_avail_by_mem)

        account_details.append(
            {
                "name": acct_name,
                "gpu_slots": gpu_slots,
                "cpu_slots": cpu_slots,
                "gpu_avail": gpu_avail,
                "cpu_avail": cpu_avail,
                "partition_gpu": partition_gpu,
                "partition_cpu": partition_cpu,
                "alloc": alloc,
                "usage": usage,
            }
        )
        total_gpu += gpu_slots
        total_cpu += cpu_slots
        total_gpu_avail += gpu_avail
        total_cpu_avail += cpu_avail

    return {
        "accounts": account_details,
        "total_gpu_slots": total_gpu,
        "total_cpu_slots": total_cpu,
        "total_slots": total_gpu + total_cpu,
        "total_gpu_avail": total_gpu_avail,
        "total_cpu_avail": total_cpu_avail,
        "total_avail": total_gpu_avail + total_cpu_avail,
    }


def generate_slurm_config(
    project_name: str,
    project_dir: str | Path,
    kintsugi_dir: str | Path,
    detected_settings: dict[str, Any] | None = None,
    user_settings: dict[str, Any] | None = None,
    experiment_config: dict[str, Any] | None = None,
) -> str:
    """
    Generate SLURM configuration file content.

    Priority order: hardcoded defaults < experiment_config < detected HPC < user CLI overrides

    Parameters
    ----------
    project_name : str
        Name of the project
    project_dir : str or Path
        Path to the project directory
    kintsugi_dir : str or Path
        Path to the KINTSUGI installation
    detected_settings : dict, optional
        Auto-detected settings from detect_hpc_settings()
    user_settings : dict, optional
        User-provided settings to override defaults
    experiment_config : dict, optional
        Experiment config from /meta/experiment.json (ExperimentConfig.to_dict()).
        Maps experiment.json keys to config.sh keys.

    Returns
    -------
    str
        Content for config.sh file
    """
    # Start with defaults
    config = get_slurm_defaults()

    # Apply experiment config (from /meta/experiment.json)
    if experiment_config:
        # Map experiment.json field names to config.sh field names
        field_map = {
            "n_cycles": "end_cycle",
            "channels_per_cycle": "end_channel",
            "tile_rows": "tile_rows",
            "tile_cols": "tile_cols",
            "tile_overlap": "tile_overlap",
            "xy_pixel_size": "xy_vox",
            "z_step_size": "z_vox",
            "numerical_aperture": "mic_na",
            "tissue_refractive_index": "tissue_ri",
        }
        for exp_key, config_key in field_map.items():
            value = experiment_config.get(exp_key)
            if value is not None:
                config[config_key] = value

        # Convert wavelengths dict to string format
        wl = experiment_config.get("wavelengths")
        if wl and isinstance(wl, dict):
            parts = []
            for ch in sorted(wl.keys(), key=lambda k: int(k)):
                ex, em = wl[str(ch)] if str(ch) in wl else wl[ch]
                parts.append(f"{ch}:{int(ex)}:{int(em)}")
            if parts:
                config["wavelengths"] = ",".join(parts)

    # Apply detected HPC settings (account, partition, GPU type, etc.)
    if detected_settings:
        for key, value in detected_settings.items():
            if value is not None:
                config[key] = value

    # Apply user overrides
    if user_settings:
        for key, value in user_settings.items():
            if value is not None:
                config[key] = value

    # Generate config file content
    content = f"""#!/bin/bash
# =============================================================================
# KINTSUGI SLURM Configuration
# Generated by: kintsugi init --slurm
# =============================================================================

# =============================================================================
# PROJECT PATHS
# =============================================================================
export PROJECT_NAME="{project_name}"
export PROJECT_DIR="{project_dir}"
export KINTSUGI_DIR="{kintsugi_dir}"

# Derived paths (usually don't need to change)
export NOTEBOOK_DIR="${{PROJECT_DIR}}/notebooks"
export DATA_DIR="${{PROJECT_DIR}}/data"
export LOG_DIR="${{PROJECT_DIR}}/slurm/logs"
export RAW_DIR="${{DATA_DIR}}/raw"
export PROCESSED_DIR="${{DATA_DIR}}/processed"

# =============================================================================
# PROCESSING PARAMETERS
# =============================================================================
# Cycles to process
export START_CYCLE={config["start_cycle"]}
export END_CYCLE={config["end_cycle"]}

# Channels
export START_CHANNEL={config["start_channel"]}
export END_CHANNEL={config["end_channel"]}

# Data format: 'tiff' or 'zarr'
export OUTPUT_FORMAT="{config["output_format"]}"

# Tile grid (for stitching)
export TILE_ROWS={config["tile_rows"]}
export TILE_COLS={config["tile_cols"]}
export TILE_OVERLAP={config["tile_overlap"]}

# =============================================================================
# MICROSCOPE PARAMETERS (for deconvolution)
# =============================================================================
export XY_VOX={config["xy_vox"]}          # nm per pixel XY
export Z_VOX={config["z_vox"]}          # nm per z-slice
export MIC_NA={config["mic_na"]}         # Numerical aperture
export TISSUE_RI={config["tissue_ri"]}      # Refractive index

# Wavelengths: channel:excitation:emission (nm)
export WAVELENGTHS="{config["wavelengths"]}"

# =============================================================================
# HPC RESOURCES
# =============================================================================
# GPU Partition Selection (HiPerGator RHEL9+):
#   - hpg-b200:  NVIDIA B200 GPUs (180GB VRAM, 8/node) - for large models/extreme performance
#   - hpg-turin: NVIDIA L4 GPUs (24GB VRAM, 3/node) - for standard GPU workloads
# See: https://docs.rc.ufl.edu/scheduler/gpu_access/
export PARTITION="{config["partition"]}"
export QOS="{config["qos"]}"
export ACCOUNT="{config["account"]}"
export GPUS_PER_NODE={config["gpus_per_node"]}

# Fallback partition settings (used when primary partition is unavailable)
# Leave empty to disable fallback and fail if primary not available
export PARTITION_FALLBACK="{config.get("partition_fallback", "")}"
export QOS_FALLBACK="{config.get("qos_fallback", "")}"

# =============================================================================
# ACCOUNT CHAIN (for multi-account GPU/CPU fallback)
# =============================================================================
# Comma-separated list of accounts to try in order
# GPU accounts first, then CPU-only burst accounts for fallback
export ACCOUNT_CHAIN="{config.get("account_chain", "")}"

# CPU-only accounts (cannot request GPUs)
# Auto-detected: accounts ending in -b are assumed CPU-only if this is empty
export CPU_ONLY_ACCOUNTS="{config.get("cpu_only_accounts", "")}"

# CPU-only partition (used with burst accounts)
export PARTITION_CPU="{config.get("partition_cpu", "hpg-default")}"

# CPU resource adjustments (CPU processing is slower than GPU)
export CPU_TIME_MULTIPLIER={config.get("cpu_time_multiplier", 5)}
export CPU_CPUS_PER_TASK={config.get("cpu_cpus_per_task", 8)}
export CPU_MEM_CORRECTION={config.get("cpu_mem_correction", 16)}
export CPU_MEM_STITCH={config.get("cpu_mem_stitch", 48)}
export CPU_MEM_DECON={config.get("cpu_mem_decon", 48)}
export CPU_MEM_EDF={config.get("cpu_mem_edf", 16)}

# Allocation limits (auto-detected from SLURM if not set)
# Uncomment and set to override auto-detection:
# export ALLOC_CPUS=104
# export ALLOC_MEM=812
# export ALLOC_GPUS=3

# Memory per job (GB)
export MEM_CORRECTION={config["mem_correction"]}
export MEM_STITCH={config["mem_stitch"]}
export MEM_DECON={config["mem_decon"]}
export MEM_EDF={config["mem_edf"]}

# Time limits (HH:MM:SS)
export TIME_CORRECTION="{config["time_correction"]}"
export TIME_STITCH="{config["time_stitch"]}"
export TIME_DECON="{config["time_decon"]}"
export TIME_EDF="{config["time_edf"]}"

# =============================================================================
# CONDA ENVIRONMENT
# =============================================================================
export CONDA_ENV="{config["conda_env"]}"

# =============================================================================
# EMAIL NOTIFICATIONS (optional)
# =============================================================================
export EMAIL="{config["email"]}"
export MAIL_TYPE="{config["mail_type"]}"
"""
    return content


def generate_slurm_readme(project_name: str, kintsugi_dir: str | Path) -> str:
    """
    Generate README content for the project's slurm directory.

    Parameters
    ----------
    project_name : str
        Name of the project
    kintsugi_dir : str or Path
        Path to the KINTSUGI installation

    Returns
    -------
    str
        Content for README.md file
    """
    return f"""# SLURM Processing Workflow for {project_name}

This guide walks you through the complete workflow from raw data to processed images.

---

## Step 1: Prerequisites

Before starting, ensure you have:

- [ ] Project created with `kintsugi init --slurm`
- [ ] SLURM access configured (account, partition)
- [ ] Raw microscopy data ready to copy

---

## Step 2: Copy Raw Data

Copy your cycle folders to the `data/raw/` directory:

```bash
# From your data source:
cp -r /path/to/cyc001 data/raw/
cp -r /path/to/cyc002 data/raw/
# ... continue for all cycles
```

**Naming conventions** (both supported):
- Short form: `cyc001/`, `cyc002/`, etc.
- Long form: `cyc001_DAPI_Blank_Blank_Blank/`, etc.

**Expected structure:**
```
data/raw/
├── cyc001/
│   └── *.tif (tile images)
├── cyc002/
│   └── *.tif
└── ...
```

---

## Step 3: Create Channel Names File

Create `meta/CHANNELNAMES.txt` with your marker names:

```bash
nano meta/CHANNELNAMES.txt
```

**Simple list format** (one marker per line, cycles in order):
```
DAPI-01
Blank
Blank
Blank
DAPI-02
CD31
CD8
CD45
```

**Cycle-prefixed format** (alternative):
```
1: DAPI, Blank, Blank, Blank
2: DAPI, CD31, CD8, CD45
```

---

## Step 4: Configure Experiment Metadata

Edit `meta/experiment.json` with your microscope parameters:

```bash
nano meta/experiment.json
```

**Key parameters to verify:**

| Parameter | Description | Example |
|-----------|-------------|---------|
| `tile_rows`, `tile_cols` | Tile grid dimensions | 5, 5 |
| `xy_pixel_size` | XY voxel size (nm) | 377.0 |
| `z_step_size` | Z step size (nm) | 1500.0 |
| `numerical_aperture` | Objective NA | 0.75 |
| `tissue_refractive_index` | Tissue RI | 1.44 |
| `wavelengths` | Channel wavelengths | See file |
| `n_cycles` | Number of cycles | Auto-detected |
| `n_zplanes` | Z-planes per stack | Auto-detected |

---

## Step 5: Review SLURM Configuration (Optional)

The `slurm/config.sh` file contains HPC-specific settings. Most parameters are
auto-detected, but you may want to review:

```bash
nano slurm/config.sh
```

**Common settings to check:**
- `ACCOUNT` - Your SLURM account
- `PARTITION` - GPU partition (e.g., hpg-b200)
- `TIME_*` - Time limits for each step
- `MEM_*` - Memory allocation

---

## Step 6: Preview Jobs (Recommended)

Before submitting, preview what will run:

```bash
kintsugi slurm submit . --dry-run
```

This shows:
- Which cycles will be processed
- Job parameters (partition, memory, time)
- Commands that will execute

**Verify:**
- [ ] Correct number of cycles detected
- [ ] Proper tile grid (e.g., 5x5)
- [ ] Expected wavelengths

---

## Step 7: Submit Processing Jobs

Submit all processing steps:

```bash
kintsugi slurm submit .
```

**Options:**
```bash
# Submit specific steps only
kintsugi slurm submit . --steps correction,stitching

# Submit specific cycles
kintsugi slurm submit . --cycles 1-3

# Submit single step for testing
kintsugi slurm submit . --steps correction --cycles 1
```

---

## Step 8: Monitor Progress

**Check job status:**
```bash
squeue -u $USER
# Or filter by project:
squeue -u $USER -n "kintsugi_*_{project_name}"
```

**View logs:**
```bash
# Find latest run directory
ls -lt slurm/runs/

# View stdout/stderr
tail -f slurm/runs/<timestamp>/logs/*.out
tail -f slurm/runs/<timestamp>/logs/*.err
```

**Check completion:**
```bash
kintsugi slurm status .
```

---

## Step 9: Review QC Images

After each step completes, review QC outputs:

**Location:** `slurm/runs/<timestamp>/qc/`

| Step | QC Images | What to Check |
|------|-----------|---------------|
| Correction | Flatfield images | Smooth illumination profiles |
| Stitching | Stitched overviews | Tile alignment, no gaps |
| Deconvolution | Before/after comparisons | Detail enhancement, no artifacts |
| EDF | Focus projections | Sharp features, proper z-selection |

---

## Step 10: Evaluate Results & Next Steps

**Output locations:**

| Step | Output Directory |
|------|------------------|
| Correction | `data/processed/corrected/` |
| Stitching | `data/processed/stitched/` |
| Deconvolution | `data/processed/deconvolved/` |
| EDF | `data/processed/edf/` |

**Next steps:**
1. Open `3_Signal_Isolation_QC.ipynb` for signal isolation and QC
2. Use Claude MCP tools for interactive processing assistance
3. Proceed to segmentation with `4_Segmentation_Analysis.ipynb`

---

## Quick Reference

### Processing Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_correction.sh` | Illumination correction (BaSiC) |
| 2 | `02_stitching.sh` | Tile stitching |
| 3 | `03_deconvolution.sh` | Richardson-Lucy deconvolution |
| 4 | `04_edf.sh` | Extended depth of focus |

### Common Commands

```bash
# Submit all steps
kintsugi slurm submit .

# Submit specific steps
kintsugi slurm submit . --steps decon,edf

# Submit specific cycles
kintsugi slurm submit . --cycles 1-3

# Preview without submitting
kintsugi slurm submit . --dry-run

# Check status
kintsugi slurm status .
```

---

## Advanced Configuration

### Partition Fallback

If the primary partition is busy, configure a fallback:

```bash
# In slurm/config.sh:
export PARTITION_FALLBACK="hpg-turin"
export QOS_FALLBACK=""
```

### Multi-Account GPU/CPU Fallback

For users with GPU and CPU-only accounts:

```bash
# In slurm/config.sh:
export ACCOUNT_CHAIN="maigan,clive,maigan-b,clive-b"
export CPU_ONLY_ACCOUNTS="maigan-b,clive-b"
export PARTITION_CPU="hpg-default"
```

CPU mode automatically adjusts:
- No GPUs requested
- More CPUs allocated
- Extended time limits

---

## Troubleshooting

| Issue | Solution |
|-------|----------|
| Job fails immediately | Check `slurm/runs/*/logs/*.err` |
| Out of memory | Increase `MEM_*` in config.sh |
| Time limit exceeded | Increase `TIME_*` in config.sh |
| Partition not found | Verify `PARTITION` setting |
| Resource unavailable | Configure fallback partition |
| Wrong parameters | Verify `meta/experiment.json` |
| Missing channels | Check `meta/CHANNELNAMES.txt` |

**Debug tip:** Always run `--dry-run` first to verify configuration.
"""
