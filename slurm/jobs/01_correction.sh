#!/bin/bash
# =============================================================================
# KINTSUGI Illumination Correction Job
# NOTE: Correction is now integrated into 02_stitching.sh for efficiency.
# This script is kept for backward compatibility but simply passes through.
# =============================================================================

# Get task ID - CYCLE can be passed via --export or derived from SLURM_ARRAY_TASK_ID
TASK_ID="${CYCLE:-${SLURM_ARRAY_TASK_ID:-0}}"

echo "============================================================"
echo "Illumination Correction - Cycle ${TASK_ID}"
echo "Job ID: ${SLURM_JOB_ID:-local} (Array Job: ${SLURM_ARRAY_JOB_ID:-N/A})"
echo "============================================================"
echo "NOTE: Correction is integrated into the stitching step."
echo "This job is a passthrough for pipeline compatibility."
echo "Actual correction happens in 02_stitching.sh"
echo "============================================================"

# Small delay to ensure SLURM accounting properly registers the array task
# This prevents "Badly formatted array jobid with task_id=-2" errors from seff
# when the job completes too quickly for SLURM's accounting daemon
sleep 2

exit 0
