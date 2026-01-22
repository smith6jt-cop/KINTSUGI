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
        # HPC resources
        "partition": "gpu",
        "qos": "",
        "account": "",
        "gpu_type": "a100",
        "gpus_per_node": 2,
        # Memory per job (GB)
        "mem_correction": 64,
        "mem_stitch": 128,
        "mem_decon": 192,
        "mem_edf": 64,
        # Time limits
        "time_correction": "04:00:00",
        "time_stitch": "06:00:00",
        "time_decon": "08:00:00",
        "time_edf": "02:00:00",
        # Processing parameters
        "start_cycle": 1,
        "end_cycle": 1,
        "start_channel": 1,
        "end_channel": 4,
        "output_format": "zarr",
        # Tile grid
        "tile_rows": 5,
        "tile_cols": 5,
        "tile_overlap": 0.1,
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
            ["sacctmgr", "show", "user", os.environ.get("USER", ""), "-n", "-p", "format=account"],
            capture_output=True,
            text=True,
            timeout=10,
        )
        if result.returncode == 0:
            accounts = [
                line.strip().rstrip("|")
                for line in result.stdout.strip().split("\n")
                if line.strip()
            ]
            return accounts
    except (subprocess.SubprocessError, FileNotFoundError, OSError):
        pass
    return []


def generate_slurm_config(
    project_name: str,
    project_dir: str | Path,
    kintsugi_dir: str | Path,
    detected_settings: dict[str, Any] | None = None,
    user_settings: dict[str, Any] | None = None,
) -> str:
    """
    Generate SLURM configuration file content.

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

    Returns
    -------
    str
        Content for config.sh file
    """
    # Start with defaults
    config = get_slurm_defaults()

    # Apply detected settings
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
    content = f'''#!/bin/bash
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
export PARTITION="{config["partition"]}"
export QOS="{config["qos"]}"
export ACCOUNT="{config["account"]}"
export GPU_TYPE="{config["gpu_type"]}"
export GPUS_PER_NODE={config["gpus_per_node"]}

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
'''
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
    return f'''# SLURM Job Submission

## Quick Start

1. Edit configuration:
   ```bash
   nano slurm/config.sh
   ```

2. Submit all steps:
   ```bash
   kintsugi slurm submit .
   # Or directly:
   {kintsugi_dir}/slurm/submit.sh --project .
   ```

3. Monitor jobs:
   ```bash
   squeue -u $USER -n "kintsugi_*_{project_name}"
   ```

## Configuration

Edit `slurm/config.sh` to set:
- `START_CYCLE`, `END_CYCLE` - Cycles to process
- `TILE_ROWS`, `TILE_COLS` - Tile grid dimensions
- `XY_VOX`, `Z_VOX` - Voxel sizes (nm)
- `WAVELENGTHS` - Channel wavelengths

## Processing Steps

| Step | Script | Description |
|------|--------|-------------|
| 1 | `01_correction.sh` | Illumination correction (BaSiC) |
| 2 | `02_stitching.sh` | Tile stitching |
| 3 | `03_deconvolution.sh` | Richardson-Lucy deconvolution |
| 4 | `04_edf.sh` | Extended depth of focus |

## Outputs

Logs and QC images are saved to `slurm/runs/<timestamp>/`.

## Commands

```bash
# Submit all steps
kintsugi slurm submit .

# Submit specific steps
kintsugi slurm submit . --steps decon,edf

# Submit for specific cycles
kintsugi slurm submit . --cycles 1-3

# Preview without submitting
kintsugi slurm submit . --dry-run
```

## Troubleshooting

- **Job fails immediately**: Check `slurm/runs/*/logs/*.err` for error details
- **Out of memory**: Increase `MEM_*` values in config.sh
- **Time limit exceeded**: Increase `TIME_*` values in config.sh
- **GPU not found**: Verify `GPU_TYPE` matches available hardware
'''
