#!/bin/bash
# =============================================================================
# KINTSUGI Signal Isolation Job
# CPU-only autofluorescence subtraction on registered images
# =============================================================================
#
# Prerequisites: Registration complete (data/processed/registered/)
# Output: Signal-isolated images in data/processed/signal_isolated/
#
# Usage:
#   sbatch --export=ALL,PROJECT_DIR=/path/to/project 06_signal_isolation.sh
#
# =============================================================================

#SBATCH --job-name=kintsugi_signal_isolation
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
init_logging "signal_isolation" "${RUN_ID:-$(generate_run_id)}"

# Record start time
START_TIME=$(date +%s)

# Generate before summary
summary_before "signal_isolation"

echo ""
log_info "Starting signal isolation (autofluorescence subtraction)..."

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${KINTSUGI_DIR}/notebooks:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
from pathlib import Path

PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

# Validate registration is complete
REG_DIR = PROJECT_DIR / 'data' / 'processed' / 'registered'
if not REG_DIR.exists() or not any(REG_DIR.glob("cyc*/*.tif")):
    print("ERROR: Registration output not found. Run registration first.")
    sys.exit(1)

print(f"Registration output validated: {REG_DIR}")

# Use the kintsugi CLI batch signal isolation
from kintsugi.signal.batch import run_batch_signal_isolation

try:
    run_batch_signal_isolation(PROJECT_DIR)
except Exception as e:
    print(f"ERROR: Signal isolation failed: {e}")
    import traceback
    traceback.print_exc()
    sys.exit(1)

print("Signal isolation complete.")
PYTHON_SCRIPT

EXIT_CODE=$?

echo ""
log_info "Signal isolation finished with exit code: ${EXIT_CODE}"

# Generate after summary
summary_after "signal_isolation" "${EXIT_CODE}" "${START_TIME}"

# Log footer
log_footer "${EXIT_CODE}"

exit ${EXIT_CODE}
