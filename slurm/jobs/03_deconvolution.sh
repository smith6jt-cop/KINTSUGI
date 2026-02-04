#!/bin/bash
# =============================================================================
# KINTSUGI Deconvolution Job
# Multi-GPU deconvolution for one cycle with logging and QC
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
init_logging "decon" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "decon"

echo ""
log_info "Starting deconvolution processing..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

# Run deconvolution
python << 'PYTHON_SCRIPT'
import sys
import os
import json
import time
from pathlib import Path

# Setup paths
PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'decon'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from Kdecon import decon
from kintsugi.gpu import cleanup_gpu_memory
import numpy as np

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

# Load metadata
experiment_config = load_experiment_config(PROJECT_DIR)

print("=" * 60)
print("Metadata Loading")
print("=" * 60)
if experiment_config:
    print(f"  Loaded experiment.json")
else:
    print(f"  experiment.json not found, using environment defaults")

# Check device mode from environment (set by submit.sh based on account type)
DEVICE_MODE = os.environ.get('KINTSUGI_DEVICE_MODE', 'gpu')
print(f"  Device mode: {DEVICE_MODE}")

# Get configuration: experiment.json -> environment variable -> default
# CYCLE can be passed via --export or derived from SLURM_ARRAY_TASK_ID
CYCLE = int(os.environ.get('CYCLE', os.environ.get('SLURM_ARRAY_TASK_ID', 1)))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL',
    experiment_config.get('channels_per_cycle', 4)))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'zarr')

# Microscope parameters (from experiment.json or environment)
XY_VOX = float(os.environ.get('XY_VOX',
    experiment_config.get('xy_pixel_size', 377)))
Z_VOX = float(os.environ.get('Z_VOX',
    experiment_config.get('z_step_size', 1500)))
MIC_NA = float(os.environ.get('MIC_NA',
    experiment_config.get('numerical_aperture', 0.75)))
TISSUE_RI = float(os.environ.get('TISSUE_RI',
    experiment_config.get('tissue_refractive_index', 1.44)))

# Parse wavelengths: experiment.json -> environment -> defaults
WAVELENGTHS = {}
if 'wavelengths' in experiment_config:
    # Load from experiment.json (keys are strings in JSON)
    for ch_str, (ex, em) in experiment_config['wavelengths'].items():
        WAVELENGTHS[int(ch_str)] = (float(ex), float(em))
    print(f"  Wavelengths loaded from experiment.json")
else:
    # Parse from environment variable
    wl_str = os.environ.get('WAVELENGTHS', '1:358:461,2:753:775,3:560:575,4:648:668')
    for item in wl_str.split(','):
        parts = item.split(':')
        if len(parts) == 3:
            ch, ex, em = int(parts[0]), float(parts[1]), float(parts[2])
            WAVELENGTHS[ch] = (ex, em)
    print(f"  Wavelengths from environment: {wl_str}")

print(f"  XY pixel size: {XY_VOX} nm")
print(f"  Z step size: {Z_VOX} nm")
print(f"  NA: {MIC_NA}, RI: {TISSUE_RI}")
print("=" * 60)

# Paths
DATA_DIR = PROJECT_DIR / 'data'
STITCH_DIR = DATA_DIR / 'processed' / 'stitched'
DECON_DIR = DATA_DIR / 'processed' / 'deconvolved'
DECON_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Initialize CUDA if in GPU mode
if DEVICE_MODE != 'cpu':
    try:
        import cupy as cp
        # Initialize CUDA context on device 0 (SLURM sets CUDA_VISIBLE_DEVICES)
        cp.cuda.Device(0).use()
        # Verify GPU is accessible with a simple operation
        _ = cp.zeros(1)
        n_gpus = 1  # SLURM allocates specific GPUs, we see them as device 0
        GPU_IDS = [0]
        print(f"CUDA initialized successfully on device 0")
        print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
    except Exception as e:
        print(f"WARNING: CUDA initialization failed: {e}")
        print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
        print("Falling back to CPU processing")
        DEVICE_MODE = 'cpu'
        n_gpus = 0
        GPU_IDS = []
else:
    print("Running in CPU mode (CPU-only burst account)")
    n_gpus = 0
    GPU_IDS = []

print(f"\n{'='*60}")
print(f"Deconvolution - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Device mode: {DEVICE_MODE}")
if DEVICE_MODE == 'gpu':
    print(f"GPUs: {n_gpus} devices")
else:
    print("Running on CPU (burst account)")
print(f"Output: {OUTPUT_FORMAT}")
print(f"Voxel size: {XY_VOX}x{XY_VOX}x{Z_VOX} nm")
print(f"QC output: {QC_DIR}")

# =============================================================================
# SKIP-EXISTING CHECK
# =============================================================================
def check_cycle_complete(decon_dir, cycle, start_ch, end_ch):
    """Check if deconvolution output already exists and is complete."""
    cycle_out = decon_dir / f"cyc{cycle:02d}"
    if not cycle_out.exists():
        return False, "output directory missing"

    for ch in range(start_ch, end_ch + 1):
        ch_dir = cycle_out / f"CH{ch}"
        if not ch_dir.exists():
            return False, f"CH{ch} directory missing"

        # Check for at least one output file
        tif_files = list(ch_dir.glob("*.tif"))
        if len(tif_files) == 0:
            return False, f"CH{ch} has no output files"

    return True, "complete"

# Check if this cycle is already processed
is_complete, status = check_cycle_complete(DECON_DIR, CYCLE, START_CHANNEL, END_CHANNEL)
if is_complete:
    print(f"\n{'='*60}")
    print(f"SKIP: Cycle {CYCLE} deconvolution already complete ({status})")
    print(f"{'='*60}")
    print(f"Output exists: {DECON_DIR / f'cyc{CYCLE:02d}'}")
    print("To reprocess, delete the output directory first.")
    sys.exit(0)
else:
    print(f"\nCycle {CYCLE} needs deconvolution: {status}")

channels = list(range(START_CHANNEL, END_CHANNEL + 1))
results = []

def save_qc_slice(data, output_path, title=""):
    """Save middle slice as QC image."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        if data.ndim == 3:
            slice_data = data[data.shape[0]//2]
        else:
            slice_data = data

        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(slice_data, cmap='gray')
        ax.set_title(title)
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"  QC saved: {output_path}")
    except Exception as e:
        print(f"  QC save failed: {e}")

def process_channel(ch):
    """Deconvolve one channel (device auto-detected from CUDA context)."""
    print(f"\n  Channel {ch} starting...")

    try:
        lambda_ex, lambda_em = WAVELENGTHS.get(ch, (560, 575))

        # Use device mode from environment (respects account type)
        decon(
            base_dir=str(PROJECT_DIR),
            stitch_dir=str(STITCH_DIR),
            dec_cycle=CYCLE,
            dec_channel=ch,
            xy_vox=XY_VOX,
            z_vox=Z_VOX,
            iterations=25,
            mic_NA=MIC_NA,
            tissue_RI=TISSUE_RI,
            damping=0,
            stop_criterion=5.0,
            device=DEVICE_MODE,
            device_id=None,  # Let decon auto-detect
            wavelengths=WAVELENGTHS,
            decon_dir=str(DECON_DIR)
        )

        # Generate QC image
        output_files = list((DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{ch}").glob("*.tif"))
        if output_files:
            from skimage.io import imread
            # Load middle z-plane for QC
            mid_idx = len(output_files) // 2
            data = imread(str(output_files[mid_idx]))
            qc_path = QC_DIR / f"cyc{CYCLE:02d}_CH{ch}_decon.png"
            save_qc_slice(data, str(qc_path), f"Cycle {CYCLE} Channel {ch} - Deconvolved")

        print(f"  Channel {ch} complete")
        return ch, True

    except Exception as e:
        print(f"  Channel {ch} FAILED: {e}")
        import traceback
        traceback.print_exc()
        return ch, False
    finally:
        cleanup_gpu_memory()

# Process channels sequentially (one GPU per SLURM job)
start_time = time.time()

print(f"\nProcessing {len(channels)} channels sequentially")
results = [process_channel(ch) for ch in channels]

# Summary
successful = sum(1 for _, ok in results if ok)
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"Deconvolution Complete")
print(f"{'='*60}")
print(f"Success: {successful}/{len(channels)} channels")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {DECON_DIR}")
print(f"QC images: {QC_DIR}")

# Exit with error if any channels failed
if successful < len(channels):
    sys.exit(1)
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Deconvolution finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "decon" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
