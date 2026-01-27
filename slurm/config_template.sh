#!/bin/bash
# =============================================================================
# KINTSUGI SLURM Configuration Template
# Copy this to your project and customize
# =============================================================================

# =============================================================================
# PROJECT PATHS (REQUIRED - customize these)
# =============================================================================
export PROJECT_NAME="my_project"
export PROJECT_DIR="/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/${PROJECT_NAME}"
export KINTSUGI_DIR="/blue/maigan/smith6jt/KINTSUGI"

# Derived paths (usually don't need to change)
export NOTEBOOK_DIR="${PROJECT_DIR}/notebooks"
export DATA_DIR="${PROJECT_DIR}/data"
export LOG_DIR="${PROJECT_DIR}/slurm/logs"
export RAW_DIR="${DATA_DIR}/raw"
export PROCESSED_DIR="${DATA_DIR}/processed"

# =============================================================================
# PROCESSING PARAMETERS
# =============================================================================
# Cycles to process
export START_CYCLE=1
export END_CYCLE=7

# Channels
export START_CHANNEL=1
export END_CHANNEL=4

# Data format: 'tiff' or 'zarr'
export OUTPUT_FORMAT="zarr"

# Tile grid (for stitching)
export TILE_ROWS=5
export TILE_COLS=5
export TILE_OVERLAP=0.1

# =============================================================================
# MICROSCOPE PARAMETERS (for deconvolution)
# =============================================================================
export XY_VOX=377          # nm per pixel XY
export Z_VOX=1500          # nm per z-slice
export MIC_NA=0.75         # Numerical aperture
export TISSUE_RI=1.44      # Refractive index

# Wavelengths: channel:excitation:emission (nm)
export WAVELENGTHS="1:358:461,2:753:775,3:560:575,4:648:668"

# =============================================================================
# HPC RESOURCES
# =============================================================================
# GPU Partition Selection (HiPerGator RHEL9+):
#   - hpg-b200:  NVIDIA B200 GPUs (180GB VRAM, 8/node) - for large models/extreme performance
#   - hpg-turin: NVIDIA L4 GPUs (24GB VRAM, 3/node) - for standard GPU workloads
# See: https://docs.rc.ufl.edu/scheduler/gpu_access/
export PARTITION="hpg-b200"
export QOS="maigan-b"
export ACCOUNT="maigan"

# Fallback GPU settings (used when primary GPU partition has no available nodes)
# When set: script checks availability and switches to fallback if primary unavailable
export PARTITION_FALLBACK="hpg-turin"
export QOS_FALLBACK=""
# Leave fallback values empty to disable automatic fallback retries; jobs will
# continue using the primary settings and may remain queued per Slurm behavior.

# -----------------------------------------------------------------------------
# ALLOCATION LIMITS (shared allocation - total resources available to you)
# These are used to automatically calculate safe concurrency levels.
# Set these to match your actual SLURM allocation to prevent over-subscription.
# -----------------------------------------------------------------------------
export ALLOC_CPUS=64       # Total CPUs in your allocation
export ALLOC_MEM=600       # Total memory (GB) in your allocation
export ALLOC_GPUS=2        # Total GPUs in your allocation

# -----------------------------------------------------------------------------
# Per-Job Resource Requests
# The script automatically calculates max concurrent jobs based on:
#   max_concurrent = min(ALLOC_CPUS/CPUS_PER_TASK, ALLOC_MEM/max_mem, ALLOC_GPUS/GPUS_PER_NODE)
# -----------------------------------------------------------------------------

# CPUs per job (must request at least 1 CPU per GPU on HiPerGator)
export CPUS_PER_TASK=4

# GPUs per job (GPU type is determined by PARTITION setting)
export GPUS_PER_NODE=1

# Memory per job (GB) - per processing step
export MEM_CORRECTION=32
export MEM_STITCH=64
export MEM_DECON=96
export MEM_EDF=32

# Manual override for max concurrent cycles (optional)
# Leave empty or set to "auto" for automatic calculation based on allocation limits
# Set to a number to force a specific concurrency level
export MAX_CONCURRENT_CYCLES="auto"

# Time limits (HH:MM:SS)
export TIME_CORRECTION="04:00:00"
export TIME_STITCH="06:00:00"
export TIME_DECON="08:00:00"
export TIME_EDF="02:00:00"

# =============================================================================
# CONDA ENVIRONMENT
# =============================================================================
export CONDA_ENV="kintsugi"

# =============================================================================
# EMAIL NOTIFICATIONS (optional)
# =============================================================================
export EMAIL=""
export MAIL_TYPE="END,FAIL"
