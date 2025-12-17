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
export PARTITION="gpu"
export QOS="maigan-b"
export ACCOUNT="maigan"
export GPU_TYPE="b200"
export GPUS_PER_NODE=2

# Memory per job (GB)
export MEM_CORRECTION=64
export MEM_STITCH=128
export MEM_DECON=192
export MEM_EDF=64

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
