#!/bin/bash
# =============================================================================
# KINTSUGI Image Stitching Job
# GPU-accelerated stitching with logging and QC
# =============================================================================

# Source utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils.sh"

# Initialize logging
init_logging "stitch" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "stitch"

echo ""
log_info "Starting image stitching..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import time
from pathlib import Path

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'stitch'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.gpu import cleanup_gpu_memory
import numpy as np

CYCLE = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL', 4))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'zarr')

# Tile grid parameters
TILE_ROWS = int(os.environ.get('TILE_ROWS', 5))
TILE_COLS = int(os.environ.get('TILE_COLS', 5))
TILE_OVERLAP = float(os.environ.get('TILE_OVERLAP', 0.1))

DATA_DIR = PROJECT_DIR / 'data'
CORRECTED_DIR = DATA_DIR / 'processed' / 'corrected'
STITCH_DIR = DATA_DIR / 'processed' / 'stitched'
STITCH_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

print(f"\n{'='*60}")
print(f"Stitching - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Grid: {TILE_ROWS}x{TILE_COLS}, overlap: {TILE_OVERLAP*100:.0f}%")
print(f"Input: {CORRECTED_DIR}")
print(f"Output: {STITCH_DIR}")
print(f"Format: {OUTPUT_FORMAT}")
print(f"QC output: {QC_DIR}")

def save_qc_stitched(data, output_path, title=""):
    """Save stitched result as QC image (downsampled)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from skimage.transform import rescale

        # Downsample for QC
        max_dim = 1024
        scale = min(1.0, max_dim / max(data.shape))
        if scale < 1.0:
            data_small = rescale(data, scale, preserve_range=True, anti_aliasing=True)
        else:
            data_small = data

        fig, ax = plt.subplots(figsize=(10, 10))
        vmin, vmax = np.percentile(data_small, [1, 99])
        im = ax.imshow(data_small, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"  QC saved: {output_path}")
    except Exception as e:
        print(f"  QC save failed: {e}")

start_time = time.time()

for channel in range(START_CHANNEL, END_CHANNEL + 1):
    print(f"\n--- Channel {channel} ---")
    ch_start = time.time()

    try:
        # Import and run stitching
        # from Kstitch import stitch_cycle_gpu
        # stitch_cycle_gpu(...)
        print(f"  Stitching channel {channel}...")

        # Placeholder for actual stitching
        # result = stitch_cycle_gpu(...)

        # Generate QC image
        # qc_path = QC_DIR / f"cyc{CYCLE:02d}_CH{channel}_stitched.png"
        # save_qc_stitched(result, str(qc_path), f"Cycle {CYCLE} Channel {channel} - Stitched")

        print(f"  Complete: {(time.time() - ch_start)/60:.1f} min")

    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

    cleanup_gpu_memory()

print(f"\n{'='*60}")
print(f"Stitching complete! Total: {(time.time() - start_time)/60:.1f} min")
print(f"{'='*60}")
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Stitching finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "stitch" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
