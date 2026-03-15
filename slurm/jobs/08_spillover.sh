#!/bin/bash
# =============================================================================
# KINTSUGI REDSEA Spillover Correction Job
# CPU-only post-segmentation spillover compensation
# =============================================================================
#
# Prerequisites: Signal-isolated multichannel TIFF + cell segmentation mask
# Output: Corrected and uncorrected per-cell quantification CSVs
#
# Usage:
#   sbatch --export=ALL,PROJECT_DIR=/path/to/project 08_spillover.sh
#
# Reference: Bai et al., Front. Immunol. 12:652631 (2021)
# =============================================================================

#SBATCH --job-name=kintsugi_spillover
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64gb
#SBATCH --time=02:00:00
#SBATCH --qos=${QOS:-${SLURM_ACCOUNT}-b}

# Source utilities (use KINTSUGI_DIR from environment, not relative path)
KINTSUGI_SLURM="${KINTSUGI_DIR}/slurm"
if [ -f "${KINTSUGI_SLURM}/utils.sh" ]; then
    source "${KINTSUGI_SLURM}/utils.sh"
else
    # Fallback: define minimal stubs if utils.sh not found
    log_info() { echo "[INFO] $*"; }
    log_error() { echo "[ERROR] $*" >&2; }
    init_logging() { :; }
    summary_before() { :; }
    summary_after() { :; }
    log_footer() { :; }
    generate_run_id() { date '+%Y%m%d_%H%M%S'; }
fi

# Initialize logging
init_logging "spillover" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "spillover"

echo ""
log_info "Starting REDSEA spillover correction..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${KINTSUGI_DIR}/notebooks:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
import tifffile

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'spillover'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.spillover import (
    run_redsea_correction,
    estimate_element_size,
    identify_surface_markers,
    compute_spillover_metrics,
    plot_spillover_comparison,
    plot_spillover_summary_heatmap,
)

# =============================================================================
# METADATA LOADING
# =============================================================================

def load_experiment_config(project_dir):
    """Load experiment configuration from /meta/experiment.json."""
    config_path = project_dir / 'meta' / 'experiment.json'
    if config_path.exists():
        try:
            with open(config_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load experiment.json: {e}")
    return {}

experiment_config = load_experiment_config(PROJECT_DIR)

print("=" * 60)
print("REDSEA Spillover Correction")
print("=" * 60)

# =============================================================================
# CONFIGURATION
# =============================================================================

# REDSEA parameters (from environment or defaults)
ELEMENT_SHAPE = os.environ.get('REDSEA_ELEMENT_SHAPE', 'star')
ELEMENT_SIZE = int(os.environ.get('REDSEA_ELEMENT_SIZE', 0))  # 0 = auto-estimate

# Input paths
DATA_DIR = PROJECT_DIR / 'data'
SIG_DIR = DATA_DIR / 'processed' / 'signal_isolated'
SEG_DIR = DATA_DIR / 'processed' / 'segmented'
OUTPUT_DIR = DATA_DIR / 'processed' / 'spillover_corrected'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Find multichannel image
image_candidates = list(SIG_DIR.glob("*.ome.tif")) + list(SIG_DIR.glob("*.tif"))
if not image_candidates:
    print(f"ERROR: No signal-isolated images found in {SIG_DIR}")
    sys.exit(1)
image_path = image_candidates[0]
print(f"Image: {image_path}")

# Find segmentation mask
mask_candidates = list(SEG_DIR.glob("*cell*mask*.tif")) + list(SEG_DIR.glob("*cell*.tif"))
if not mask_candidates:
    # Try alternative naming
    mask_candidates = list(SEG_DIR.glob("*.tif"))
if not mask_candidates:
    print(f"ERROR: No segmentation masks found in {SEG_DIR}")
    sys.exit(1)
mask_path = mask_candidates[0]
print(f"Mask: {mask_path}")

# Load channel names
from Kio import load_channel_names
channels_per_cycle = experiment_config.get('channels_per_cycle', 4)
channel_name_dict = load_channel_names(PROJECT_DIR / 'meta', channels_per_cycle=channels_per_cycle)

# Flatten to a single list of unique marker names
marker_names = []
if channel_name_dict:
    for cycle_num in sorted(channel_name_dict.keys()):
        for name in channel_name_dict[cycle_num]:
            if name not in marker_names and name.lower() not in {'blank', 'empty', 'dapi'}:
                marker_names.append(name)
    # Add DAPI at the beginning (if present)
    if any('dapi' in str(n).lower() for names in channel_name_dict.values() for n in names):
        marker_names.insert(0, "DAPI")

if not marker_names:
    print("ERROR: Could not load marker names from CHANNELNAMES.txt")
    sys.exit(1)

print(f"Markers: {len(marker_names)} ({marker_names[:5]}...)")

# =============================================================================
# LOAD DATA
# =============================================================================
start_time = time.time()

image = tifffile.imread(str(image_path))
if image.ndim == 3 and image.shape[0] < image.shape[-1]:
    image = image.transpose(1, 2, 0)  # (C, Y, X) -> (Y, X, C)
print(f"Image shape: {image.shape}")

mask = tifffile.imread(str(mask_path))
print(f"Mask shape: {mask.shape}, {np.unique(mask).max()} cells")

# Validate marker count matches channels
n_channels = image.shape[-1] if image.ndim == 3 else 1
if len(marker_names) != n_channels:
    print(f"WARNING: {len(marker_names)} marker names but {n_channels} channels")
    # Truncate or pad as needed
    if len(marker_names) > n_channels:
        marker_names = marker_names[:n_channels]
    else:
        while len(marker_names) < n_channels:
            marker_names.append(f"CH{len(marker_names)+1}")

# =============================================================================
# PARAMETER ESTIMATION
# =============================================================================
if ELEMENT_SIZE == 0:
    ELEMENT_SIZE = estimate_element_size(mask)
    print(f"Auto-estimated element_size: {ELEMENT_SIZE}")

surface_markers = identify_surface_markers(marker_names)
print(f"Surface markers: {len(surface_markers)}")

# =============================================================================
# RUN REDSEA
# =============================================================================
print(f"\nRunning REDSEA correction...")
print(f"  element_shape: {ELEMENT_SHAPE}")
print(f"  element_size: {ELEMENT_SIZE}")
print(f"  markers_of_interest: {surface_markers}")

result = run_redsea_correction(
    multichannel_image=image,
    segmentation_mask=mask,
    marker_names=marker_names,
    markers_of_interest=surface_markers if surface_markers else None,
    element_shape=ELEMENT_SHAPE,
    element_size=ELEMENT_SIZE,
    output_dir=OUTPUT_DIR,
)

print(f"\nCorrection complete: {result.n_cells} cells x {result.n_markers} markers")

# =============================================================================
# QC
# =============================================================================
try:
    exclusive_pairs = [
        ("CD3", "CD20"), ("CD3e", "CD20"),
        ("CD4", "CD8"), ("CD4", "CD8a"),
    ]
    available_pairs = [
        (a, b) for a, b in exclusive_pairs
        if a in result.corrected.columns and b in result.corrected.columns
    ]

    if available_pairs:
        metrics = compute_spillover_metrics(
            result.corrected, result.uncorrected, available_pairs
        )
        print("\nDouble-positive reduction:")
        for pair, m in metrics["pair_metrics"].items():
            print(f"  {pair[0]}/{pair[1]}: "
                  f"{m['uncorrected_dp_pct']:.1f}% -> {m['corrected_dp_pct']:.1f}%")

        for pair in available_pairs:
            plot_spillover_comparison(
                result.corrected, result.uncorrected, pair,
                save_path=QC_DIR / f"qc_{pair[0]}_{pair[1]}.pdf",
            )

    plot_spillover_summary_heatmap(
        result.corrected, result.uncorrected, marker_names,
        save_path=QC_DIR / "qc_summary_heatmap.pdf",
    )
    print(f"QC plots saved to: {QC_DIR}")
except Exception as e:
    print(f"WARNING: QC visualization failed: {e}")

# =============================================================================
# COMPLETION
# =============================================================================
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"Spillover Correction Complete")
print(f"{'='*60}")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {OUTPUT_DIR}")

# Write completion marker
from datetime import datetime
marker_file = OUTPUT_DIR / ".complete"
with open(marker_file, 'w') as f:
    f.write(f"stage=spillover_correction\n")
    f.write(f"completed={datetime.now().isoformat()}\n")
    f.write(f"n_cells={result.n_cells}\n")
    f.write(f"n_markers={result.n_markers}\n")
    f.write(f"element_shape={ELEMENT_SHAPE}\n")
    f.write(f"element_size={ELEMENT_SIZE}\n")
    f.write(f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
    f.write(f"duration_minutes={elapsed:.1f}\n")

os.sync()
print(f"Completion marker written: {marker_file}")
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Spillover correction finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "spillover" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
