#!/bin/bash
# =============================================================================
# KINTSUGI Extended Depth of Focus (EDF) Job
# Multi-GPU EDF processing with logging and QC
# =============================================================================

# Source utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils.sh"

# Initialize logging
init_logging "edf" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "edf"

echo ""
log_info "Starting EDF processing..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import time
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
from itertools import cycle as iter_cycle

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'edf'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.edf import EDFProcessor
from kintsugi.gpu import cleanup_gpu_memory
import cupy as cp
import numpy as np

CYCLE = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL', 4))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'zarr')

DATA_DIR = PROJECT_DIR / 'data'
DECON_DIR = DATA_DIR / 'processed' / 'deconvolved'
EDF_DIR = DATA_DIR / 'processed' / 'edf'
EDF_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

n_gpus = cp.cuda.runtime.getDeviceCount()
GPU_IDS = list(range(n_gpus))

print(f"\n{'='*60}")
print(f"EDF Processing - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"GPUs: {n_gpus}")
print(f"Input: {DECON_DIR}")
print(f"Output: {EDF_DIR}")
print(f"QC output: {QC_DIR}")

channels = list(range(START_CHANNEL, END_CHANNEL + 1))

def save_qc_image(data, output_path, title=""):
    """Save EDF result as QC image."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 8))
        # Auto-contrast
        vmin, vmax = np.percentile(data, [1, 99])
        im = ax.imshow(data, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"  QC saved: {output_path}")
    except Exception as e:
        print(f"  QC save failed: {e}")

def process_edf(args):
    """Process EDF for one channel."""
    ch, gpu_id = args
    print(f"\n  [GPU{gpu_id}] Channel {ch} starting...")

    try:
        cp.cuda.Device(gpu_id).use()

        processor = EDFProcessor(device='gpu')

        input_path = DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{ch}"
        output_path = EDF_DIR / f"cyc{CYCLE:02d}"
        output_path.mkdir(parents=True, exist_ok=True)

        # Find input files
        input_files = sorted(input_path.glob("*.tif"))
        if not input_files:
            print(f"  [GPU{gpu_id}] No input files found for channel {ch}")
            return ch, False

        # Load z-stack
        from skimage.io import imread, imsave
        stack = np.array([imread(str(f)) for f in input_files])

        # Run EDF
        edf_result = processor.process(stack)

        # Save result
        output_file = output_path / f"CH{ch}_edf.tif"
        imsave(str(output_file), edf_result.astype(np.uint16))

        # Generate QC image
        qc_path = QC_DIR / f"cyc{CYCLE:02d}_CH{ch}_edf.png"
        save_qc_image(edf_result, str(qc_path), f"Cycle {CYCLE} Channel {ch} - EDF")

        print(f"  [GPU{gpu_id}] Channel {ch} complete -> {output_file}")
        return ch, True

    except Exception as e:
        print(f"  [GPU{gpu_id}] Channel {ch} FAILED: {e}")
        import traceback
        traceback.print_exc()
        return ch, False
    finally:
        cleanup_gpu_memory()

# Process channels
start_time = time.time()

if n_gpus >= 2 and len(channels) >= 2:
    print(f"\n[MULTI-GPU] {len(channels)} channels across {n_gpus} GPUs")
    pairs = list(zip(channels, iter_cycle(GPU_IDS)))
    with ThreadPoolExecutor(max_workers=n_gpus) as executor:
        results = list(executor.map(process_edf, pairs))
else:
    results = [process_edf((ch, 0)) for ch in channels]

successful = sum(1 for _, ok in results if ok)
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"EDF Complete")
print(f"{'='*60}")
print(f"Success: {successful}/{len(channels)} channels")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {EDF_DIR}")
print(f"QC images: {QC_DIR}")

if successful < len(channels):
    sys.exit(1)
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "EDF processing finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "edf" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
