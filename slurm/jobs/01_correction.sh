#!/bin/bash
# =============================================================================
# KINTSUGI Illumination Correction Job (BaSiC)
# GPU-accelerated correction with logging and QC
# =============================================================================

# Source utilities
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/../utils.sh"

# Initialize logging
init_logging "correction" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "correction"

echo ""
log_info "Starting illumination correction..."

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
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'correction'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.gpu import cleanup_gpu_memory
import numpy as np

CYCLE = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL', 4))

DATA_DIR = PROJECT_DIR / 'data'
RAW_DIR = DATA_DIR / 'raw'
CORRECTED_DIR = DATA_DIR / 'processed' / 'corrected'
CORRECTED_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

print(f"\n{'='*60}")
print(f"Illumination Correction - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Raw: {RAW_DIR}")
print(f"Output: {CORRECTED_DIR}")
print(f"QC output: {QC_DIR}")

# Find cycle directory
cycle_patterns = [f"cyc{CYCLE:03d}", f"cyc{CYCLE:02d}", f"Cyc{CYCLE:02d}"]
cycle_dir = None
for pattern in cycle_patterns:
    test_dir = RAW_DIR / pattern
    if test_dir.exists():
        cycle_dir = test_dir
        break

if cycle_dir is None:
    print(f"ERROR: Cycle directory not found in {RAW_DIR}")
    sys.exit(1)

print(f"Found: {cycle_dir}")

def save_qc_comparison(before, after, output_path, title=""):
    """Save before/after comparison as QC image."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(1, 2, figsize=(12, 6))

        vmin = min(np.percentile(before, 1), np.percentile(after, 1))
        vmax = max(np.percentile(before, 99), np.percentile(after, 99))

        axes[0].imshow(before, cmap='gray', vmin=vmin, vmax=vmax)
        axes[0].set_title('Before Correction')
        axes[0].axis('off')

        axes[1].imshow(after, cmap='gray', vmin=vmin, vmax=vmax)
        axes[1].set_title('After Correction')
        axes[1].axis('off')

        plt.suptitle(title)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
        print(f"  QC saved: {output_path}")
    except Exception as e:
        print(f"  QC save failed: {e}")

corrector = KCorrectGPU(use_gpu=True, verbose=True)

for channel in range(START_CHANNEL, END_CHANNEL + 1):
    print(f"\n--- Channel {channel} ---")
    start = time.time()

    try:
        # Find channel images - adjust patterns for your naming convention
        print(f"  Processing channel {channel}...")
        # corrector.process(input_dir, output_dir, ...)
        print(f"  Complete: {time.time() - start:.1f}s")

        # Generate QC image (placeholder - actual implementation depends on data)
        # qc_path = QC_DIR / f"cyc{CYCLE:02d}_CH{channel}_correction.png"
        # save_qc_comparison(before, after, str(qc_path), f"Cycle {CYCLE} Channel {channel}")

    except Exception as e:
        print(f"  ERROR: {e}")
        import traceback
        traceback.print_exc()

    cleanup_gpu_memory()

print(f"\n{'='*60}")
print("Correction complete!")
print(f"{'='*60}")
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Illumination correction finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "correction" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
