#!/bin/bash
# =============================================================================
# KINTSUGI Illumination Correction Job
# NOTE: Correction is now integrated into 02_stitching.sh for efficiency.
# This script is kept for backward compatibility but simply passes through.
# =============================================================================

echo "============================================================"
echo "Illumination Correction - Cycle ${SLURM_ARRAY_TASK_ID}"
echo "============================================================"
echo "NOTE: Correction is integrated into the stitching step."
echo "This job is a passthrough for pipeline compatibility."
echo "Actual correction happens in 02_stitching.sh"
echo "============================================================"

exit 0
