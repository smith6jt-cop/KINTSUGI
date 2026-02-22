#!/bin/bash
# =============================================================================
# KINTSUGI 3D Vessel Segmentation Job
# Frangi vesselness + skeletonization + graph morphometry (automated/batch)
#
# Pipeline position: After deconvolution (step 3), parallel with EDF (step 4)
# Input: data/processed/deconvolved/cyc##/CH#/*.tif
# Output: data/processed/vessel_3d/
# =============================================================================

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
init_logging "vessel3d" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "vessel3d"

echo ""
log_info "Starting 3D Vessel Segmentation..."

module load conda
conda activate ${CONDA_ENV:-KINTSUGI}

export PYTHONPATH="${KINTSUGI_DIR}:${KINTSUGI_DIR}/notebooks:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import json
import time
import gc
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor

import numpy as np

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'vessel3d'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.vessel3d import (
    VesselSpacing,
    segment_vessels_3d,
    segment_vessels_multichannel,
    export_vessel_results,
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
print("Metadata Loading")
print("=" * 60)
if experiment_config:
    print(f"  Loaded experiment.json")
else:
    print(f"  experiment.json not found, using defaults")

# Physical spacing
spacing = VesselSpacing.from_experiment(experiment_config)
print(f"  Spacing: XY={spacing.xy} um, Z={spacing.z} um, ratio={spacing.ratio:.2f}")

# Device mode — force CPU if GPU VRAM < 40 GB (L4 = 23 GB, too small for isotropic volumes)
DEVICE_MODE = os.environ.get('KINTSUGI_DEVICE_MODE', 'gpu')
try:
    import subprocess
    nvidia_out = subprocess.check_output(
        ['nvidia-smi', '--query-gpu=memory.total', '--format=csv,noheader,nounits'],
        text=True, timeout=10,
    ).strip()
    gpu_vram_mb = int(nvidia_out.split('\n')[0])
    if gpu_vram_mb < 40000 and DEVICE_MODE != 'cpu':
        print(f"  WARNING: GPU VRAM = {gpu_vram_mb} MB (< 40 GB). Forcing device='cpu'.")
        DEVICE_MODE = 'cpu'
except Exception as e:
    print(f"  Could not query GPU VRAM ({e}), using {DEVICE_MODE}")
print(f"  Device mode: {DEVICE_MODE}")

# Vessel marker parameters from environment
VESSEL_CYCLE = int(os.environ.get('VESSEL_CYCLE', 2))
VESSEL_CHANNEL = int(os.environ.get('VESSEL_CHANNEL', 2))
VESSEL_MARKER = os.environ.get('VESSEL_MARKER', 'CD31')

# Sensitivity tuning
VESSEL_PRESET = os.environ.get('VESSEL_PRESET', '') or None

# Multi-channel combination (e.g., "CD31+SMA")
VESSEL_COMBINE = os.environ.get('VESSEL_COMBINE', '') or None
VESSEL_COMBINE_CYCLE = os.environ.get('VESSEL_COMBINE_CYCLE', '')
VESSEL_COMBINE_CHANNEL = os.environ.get('VESSEL_COMBINE_CHANNEL', '')
VESSEL_COMBINE_MARKER = os.environ.get('VESSEL_COMBINE_MARKER', '')

# Build kwargs — only include params that were explicitly set via env vars.
# Params left out stay None in segment_vessels_3d() so presets can fill them.
seg_kwargs = dict(
    spacing=spacing,
    make_isotropic=True,
    device=DEVICE_MODE,
    marker_name=VESSEL_MARKER,
    preset=VESSEL_PRESET,
)
if 'VESSEL_SIGMAS' in os.environ:
    seg_kwargs['sigmas'] = [float(s) for s in os.environ['VESSEL_SIGMAS'].split(',')]
if 'VESSEL_ALPHA' in os.environ:
    seg_kwargs['alpha'] = float(os.environ['VESSEL_ALPHA'])
if 'VESSEL_BETA' in os.environ:
    seg_kwargs['beta'] = float(os.environ['VESSEL_BETA'])
if 'VESSEL_THRESHOLD' in os.environ and os.environ['VESSEL_THRESHOLD']:
    seg_kwargs['threshold'] = float(os.environ['VESSEL_THRESHOLD'])
if 'VESSEL_MIN_SIZE' in os.environ:
    seg_kwargs['min_size'] = int(os.environ['VESSEL_MIN_SIZE'])
if 'VESSEL_CLOSING_RADIUS' in os.environ:
    seg_kwargs['closing_radius'] = int(os.environ['VESSEL_CLOSING_RADIUS'])
if 'VESSEL_DENOISE_SIGMA' in os.environ:
    seg_kwargs['denoise_sigma'] = float(os.environ['VESSEL_DENOISE_SIGMA'])
if 'VESSEL_PRUNE_MIN_UM' in os.environ:
    seg_kwargs['prune_min_length_um'] = float(os.environ['VESSEL_PRUNE_MIN_UM'])

print(f"  Vessel marker: Cycle {VESSEL_CYCLE}, CH{VESSEL_CHANNEL} ({VESSEL_MARKER})")
if VESSEL_PRESET:
    print(f"  Preset: {VESSEL_PRESET}")
print(f"  Explicit overrides: {[k for k in seg_kwargs if k not in ('spacing', 'make_isotropic', 'device', 'marker_name', 'preset')]}")
if VESSEL_COMBINE:
    print(f"  Multi-channel: {VESSEL_COMBINE} (cycle {VESSEL_COMBINE_CYCLE}, "
          f"CH{VESSEL_COMBINE_CHANNEL}, {VESSEL_COMBINE_MARKER})")

# =============================================================================
# INPUT VALIDATION
# =============================================================================

DATA_DIR = PROJECT_DIR / 'data'
DECON_DIR = DATA_DIR / 'processed' / 'deconvolved'
OUTPUT_DIR = DATA_DIR / 'processed' / 'vessel_3d'
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

input_dir = DECON_DIR / f"cyc{VESSEL_CYCLE:02d}" / f"CH{VESSEL_CHANNEL}"
if not input_dir.exists():
    print(f"ERROR: Deconvolution output not found: {input_dir}")
    sys.exit(1)

# Check if already processed
output_csv = OUTPUT_DIR / f"vessel_features_{VESSEL_MARKER}.csv"
if output_csv.exists():
    print(f"\nSKIP: Vessel segmentation already complete for {VESSEL_MARKER}")
    print(f"Output exists: {output_csv}")
    print("To reprocess, delete the output files first.")
    sys.exit(0)

# =============================================================================
# LOAD Z-STACK
# =============================================================================

print(f"\n{'='*60}")
print(f"3D Vessel Segmentation - {VESSEL_MARKER}")
print(f"{'='*60}")

from skimage.io import imread
from natsort import natsorted

tiff_files = natsorted(list(input_dir.glob('*.tif')))
print(f"Loading {len(tiff_files)} z-planes from {input_dir}...")

# Parallel loading
def load_single(f):
    return imread(str(f))

with ThreadPoolExecutor(max_workers=4) as executor:
    slices = list(executor.map(load_single, tiff_files))

volume = np.stack(slices, axis=0)
del slices
gc.collect()

print(f"Volume: {volume.shape}, {volume.dtype}, {volume.nbytes / 1e9:.2f} GB")

# =============================================================================
# RUN PIPELINE
# =============================================================================

start_time = time.time()

if VESSEL_COMBINE:
    # Multi-channel mode: load second volume and run combined pipeline
    combine_cycle = int(VESSEL_COMBINE_CYCLE)
    combine_channel = int(VESSEL_COMBINE_CHANNEL)
    combine_marker = VESSEL_COMBINE_MARKER

    input_dir2 = DECON_DIR / f"cyc{combine_cycle:02d}" / f"CH{combine_channel}"
    if not input_dir2.exists():
        print(f"ERROR: Second channel not found: {input_dir2}")
        sys.exit(1)

    tiff_files2 = natsorted(list(input_dir2.glob('*.tif')))
    print(f"Loading {len(tiff_files2)} z-planes for {combine_marker} from {input_dir2}...")
    with ThreadPoolExecutor(max_workers=4) as executor:
        slices2 = list(executor.map(load_single, tiff_files2))
    volume2 = np.stack(slices2, axis=0)
    del slices2
    gc.collect()
    print(f"Volume 2: {volume2.shape}, {volume2.dtype}, {volume2.nbytes / 1e9:.2f} GB")

    volumes = {VESSEL_MARKER: volume, combine_marker: volume2}
    result = segment_vessels_multichannel(
        volumes,
        combine_method='union',
        **seg_kwargs,
    )
    del volume, volume2
    gc.collect()

    # Update marker name for output files
    VESSEL_MARKER = result.marker_name
else:
    # Single-channel mode
    result = segment_vessels_3d(volume, **seg_kwargs)
    del volume
    gc.collect()

# =============================================================================
# EXPORT
# =============================================================================

outputs = export_vessel_results(
    result,
    output_dir=OUTPUT_DIR,
    save_vesselness=False,
)

elapsed = (time.time() - start_time) / 60

# =============================================================================
# QC IMAGE
# =============================================================================

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    # MIP of mask
    mip = result.binary_mask.max(axis=0).astype(np.float32)
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.imshow(mip, cmap='gray')
    ax.set_title(f"{VESSEL_MARKER} Vessel Mask MIP")
    ax.axis('off')
    qc_path = QC_DIR / f"{VESSEL_MARKER}_vessel_mask_mip.png"
    plt.savefig(str(qc_path), dpi=100, bbox_inches='tight')
    plt.close()
    print(f"QC saved: {qc_path}")
except Exception as e:
    print(f"QC save failed: {e}")

# =============================================================================
# SUMMARY
# =============================================================================

print(f"\n{'='*60}")
print(f"Vessel Segmentation Complete")
print(f"{'='*60}")
print(f"Marker: {VESSEL_MARKER}")
print(f"Time: {elapsed:.1f} minutes")
s = result.summary()
if s:
    print(f"Total segments: {s['n_segments']}")
    print(f"  Junction-junction: {s.get('n_junctions', 'N/A')}")
    print(f"  Endpoint-endpoint: {s.get('n_endpoints', 'N/A')}")
    print(f"Total network length: {s['total_length_um']:.0f} um")
    print(f"Total network volume: {s.get('total_volume_um3', 0):.0f} um^3")
    print(f"Mean diameter: {s['mean_diameter_um']:.1f} um")
    print(f"Mean tortuosity: {s['mean_tortuosity']:.2f}")
    ba = s.get('mean_branching_angle_deg')
    print(f"Mean branching angle: {ba:.1f} deg" if ba is not None else "Mean branching angle: N/A")
elif result.features is not None:
    print(f"Total segments: {len(result.features)}")
print(f"Output: {OUTPUT_DIR}")

# Write completion marker
from datetime import datetime
marker_file = OUTPUT_DIR / ".complete"
try:
    with open(marker_file, 'w') as f:
        f.write(f"stage=vessel3d\n")
        f.write(f"marker={VESSEL_MARKER}\n")
        f.write(f"completed={datetime.now().isoformat()}\n")
        f.write(f"duration_minutes={elapsed:.1f}\n")
        f.write(f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
    os.sync()
    print(f"Completion marker written: {marker_file}")
except Exception as e:
    print(f"WARNING: Failed to write completion marker: {e}")

PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Vessel segmentation finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "vessel3d" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
