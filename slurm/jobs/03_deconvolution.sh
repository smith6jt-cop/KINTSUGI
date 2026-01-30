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

# Get configuration from environment
CYCLE = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL', 4))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'zarr')

# Microscope parameters
XY_VOX = float(os.environ.get('XY_VOX', 377))
Z_VOX = float(os.environ.get('Z_VOX', 1500))
MIC_NA = float(os.environ.get('MIC_NA', 0.75))
TISSUE_RI = float(os.environ.get('TISSUE_RI', 1.44))

# Parse wavelengths from environment
WAVELENGTHS = {}
wl_str = os.environ.get('WAVELENGTHS', '1:358:461,2:753:775,3:560:575,4:648:668')
for item in wl_str.split(','):
    parts = item.split(':')
    if len(parts) == 3:
        ch, ex, em = int(parts[0]), float(parts[1]), float(parts[2])
        WAVELENGTHS[ch] = (ex, em)

# Paths
DATA_DIR = PROJECT_DIR / 'data'
STITCH_DIR = DATA_DIR / 'processed' / 'stitched'
DECON_DIR = DATA_DIR / 'processed' / 'deconvolved'
DECON_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Initialize CUDA properly - let CuPy auto-detect and initialize
# Don't explicitly set device_id to avoid context conflicts
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
    n_gpus = 0
    GPU_IDS = []

print(f"\n{'='*60}")
print(f"Deconvolution - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"GPUs: {n_gpus} devices")
print(f"Output: {OUTPUT_FORMAT}")
print(f"Voxel size: {XY_VOX}x{XY_VOX}x{Z_VOX} nm")
print(f"QC output: {QC_DIR}")

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

        # Don't pass device_id - let decon auto-detect from initialized context
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
            device='gpu' if n_gpus > 0 else 'cpu',
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
