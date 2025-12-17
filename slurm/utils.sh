#!/bin/bash
# =============================================================================
# KINTSUGI SLURM Utilities
# Shared logging, summary, and QC functions
# =============================================================================

# =============================================================================
# LOGGING
# =============================================================================

# Initialize logging for a job step
init_logging() {
    local step_name=$1
    local run_id=$2

    export STEP_NAME="${step_name}"
    export RUN_ID="${run_id}"
    export RUN_DIR="${PROJECT_DIR}/slurm/runs/${RUN_ID}"
    export LOG_FILE="${RUN_DIR}/logs/${step_name}_cyc${SLURM_ARRAY_TASK_ID:-0}.log"
    export QC_DIR="${RUN_DIR}/qc/${step_name}"
    export SUMMARY_DIR="${RUN_DIR}/summaries"

    mkdir -p "$(dirname "${LOG_FILE}")"
    mkdir -p "${QC_DIR}"
    mkdir -p "${SUMMARY_DIR}"

    # Start logging
    exec > >(tee -a "${LOG_FILE}") 2>&1

    log_header
}

log_header() {
    echo "============================================================"
    echo "KINTSUGI ${STEP_NAME} - $(date)"
    echo "============================================================"
    echo "Run ID:      ${RUN_ID}"
    echo "Job ID:      ${SLURM_JOB_ID:-local}"
    echo "Array Task:  ${SLURM_ARRAY_TASK_ID:-N/A}"
    echo "Project:     ${PROJECT_DIR}"
    echo "Node:        ${SLURMD_NODENAME:-$(hostname)}"
    echo "GPUs:        ${CUDA_VISIBLE_DEVICES:-none}"
    echo "============================================================"
}

log_footer() {
    local exit_code=$1
    echo ""
    echo "============================================================"
    echo "Completed: $(date)"
    echo "Exit code: ${exit_code}"
    echo "Log file:  ${LOG_FILE}"
    echo "============================================================"
}

log_info() {
    echo "[INFO] $(date '+%H:%M:%S') $*"
}

log_warn() {
    echo "[WARN] $(date '+%H:%M:%S') $*" >&2
}

log_error() {
    echo "[ERROR] $(date '+%H:%M:%S') $*" >&2
}

# =============================================================================
# DATA SUMMARY - Before/After State
# =============================================================================

# Count files and get total size
count_files() {
    local dir=$1
    local pattern=${2:-"*"}

    if [ -d "${dir}" ]; then
        local count=$(find "${dir}" -name "${pattern}" -type f 2>/dev/null | wc -l)
        local size=$(du -sh "${dir}" 2>/dev/null | cut -f1)
        echo "${count} files (${size})"
    else
        echo "0 files (directory not found)"
    fi
}

# Generate data summary for a directory
summarize_directory() {
    local dir=$1
    local label=$2

    echo "${label}:"
    if [ -d "${dir}" ]; then
        echo "  Path: ${dir}"
        echo "  TIFFs: $(find "${dir}" -name "*.tif" -o -name "*.tiff" 2>/dev/null | wc -l)"
        echo "  Zarrs: $(find "${dir}" -name "*.zarr" -type d 2>/dev/null | wc -l)"
        echo "  Total size: $(du -sh "${dir}" 2>/dev/null | cut -f1)"
    else
        echo "  (not found)"
    fi
}

# Generate before-processing summary
summary_before() {
    local step_name=$1
    local summary_file="${SUMMARY_DIR}/${step_name}_cyc${SLURM_ARRAY_TASK_ID:-0}_before.txt"

    {
        echo "============================================================"
        echo "BEFORE: ${step_name} - Cycle ${SLURM_ARRAY_TASK_ID:-all}"
        echo "Time: $(date)"
        echo "============================================================"
        echo ""

        case ${step_name} in
            correction|stitch)
                summarize_directory "${PROJECT_DIR}/data/raw" "Raw Data"
                ;;
            decon)
                summarize_directory "${PROJECT_DIR}/data/processed/stitched" "Stitched Data"
                ;;
            edf)
                summarize_directory "${PROJECT_DIR}/data/processed/deconvolved" "Deconvolved Data"
                ;;
            registration)
                summarize_directory "${PROJECT_DIR}/data/processed/edf" "EDF Data"
                ;;
        esac

        echo ""
        echo "GPU Memory (before):"
        if command -v nvidia-smi &> /dev/null; then
            nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "  N/A"
        else
            echo "  nvidia-smi not available"
        fi

    } > "${summary_file}"

    log_info "Before summary: ${summary_file}"
    cat "${summary_file}"
}

# Generate after-processing summary
summary_after() {
    local step_name=$1
    local exit_code=$2
    local start_time=$3
    local summary_file="${SUMMARY_DIR}/${step_name}_cyc${SLURM_ARRAY_TASK_ID:-0}_after.txt"

    local end_time=$(date +%s)
    local duration=$((end_time - start_time))
    local duration_min=$((duration / 60))
    local duration_sec=$((duration % 60))

    {
        echo "============================================================"
        echo "AFTER: ${step_name} - Cycle ${SLURM_ARRAY_TASK_ID:-all}"
        echo "Time: $(date)"
        echo "============================================================"
        echo ""
        echo "Exit code: ${exit_code}"
        echo "Duration: ${duration_min}m ${duration_sec}s"
        echo ""

        case ${step_name} in
            correction)
                summarize_directory "${PROJECT_DIR}/data/processed/corrected" "Corrected Data"
                ;;
            stitch)
                summarize_directory "${PROJECT_DIR}/data/processed/stitched" "Stitched Data"
                ;;
            decon)
                summarize_directory "${PROJECT_DIR}/data/processed/deconvolved" "Deconvolved Data"
                ;;
            edf)
                summarize_directory "${PROJECT_DIR}/data/processed/edf" "EDF Data"
                ;;
            registration)
                summarize_directory "${PROJECT_DIR}/data/processed/registered" "Registered Data"
                ;;
        esac

        echo ""
        echo "GPU Memory (after):"
        if command -v nvidia-smi &> /dev/null; then
            nvidia-smi --query-gpu=index,memory.used,memory.total --format=csv,noheader 2>/dev/null || echo "  N/A"
        else
            echo "  nvidia-smi not available"
        fi

        echo ""
        echo "Output files created:"
        case ${step_name} in
            correction)
                find "${PROJECT_DIR}/data/processed/corrected" -name "*.tif" -newer "${summary_file%_after.txt}_before.txt" 2>/dev/null | head -20
                ;;
            stitch)
                find "${PROJECT_DIR}/data/processed/stitched" -name "*.tif" -newer "${summary_file%_after.txt}_before.txt" 2>/dev/null | head -20
                ;;
            decon)
                find "${PROJECT_DIR}/data/processed/deconvolved" \( -name "*.tif" -o -name "*.zarr" \) -newer "${summary_file%_after.txt}_before.txt" 2>/dev/null | head -20
                ;;
            edf)
                find "${PROJECT_DIR}/data/processed/edf" -name "*.tif" -newer "${summary_file%_after.txt}_before.txt" 2>/dev/null | head -20
                ;;
        esac

    } > "${summary_file}"

    log_info "After summary: ${summary_file}"
    cat "${summary_file}"
}

# =============================================================================
# QC IMAGE GENERATION
# =============================================================================

# Save a QC thumbnail (requires Python with PIL/skimage)
save_qc_thumbnail() {
    local input_file=$1
    local output_name=$2
    local max_dim=${3:-512}

    local output_file="${QC_DIR}/${output_name}"

    python3 << PYEOF
import sys
from pathlib import Path

try:
    import numpy as np
    from skimage.io import imread, imsave
    from skimage.transform import rescale

    img = imread("${input_file}")

    # Handle 3D stacks - take middle slice or max projection
    if img.ndim == 3:
        img = img[img.shape[0]//2]  # Middle slice

    # Rescale if needed
    scale = min(1.0, ${max_dim} / max(img.shape))
    if scale < 1.0:
        img = rescale(img, scale, preserve_range=True, anti_aliasing=True)

    # Normalize to 8-bit for PNG
    img = img.astype(np.float32)
    img = (img - img.min()) / (img.max() - img.min() + 1e-10) * 255
    img = img.astype(np.uint8)

    imsave("${output_file}", img)
    print(f"Saved: ${output_file}")

except Exception as e:
    print(f"QC thumbnail failed: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

# Generate comparison QC image (before/after side by side)
save_qc_comparison() {
    local before_file=$1
    local after_file=$2
    local output_name=$3

    local output_file="${QC_DIR}/${output_name}"

    python3 << PYEOF
import sys
from pathlib import Path

try:
    import numpy as np
    from skimage.io import imread, imsave
    from skimage.transform import rescale
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    def load_and_normalize(path, max_dim=400):
        img = imread(path)
        if img.ndim == 3:
            img = img[img.shape[0]//2]
        scale = min(1.0, max_dim / max(img.shape))
        if scale < 1.0:
            img = rescale(img, scale, preserve_range=True, anti_aliasing=True)
        img = img.astype(np.float32)
        img = (img - img.min()) / (img.max() - img.min() + 1e-10)
        return img

    before = load_and_normalize("${before_file}")
    after = load_and_normalize("${after_file}")

    fig, axes = plt.subplots(1, 2, figsize=(10, 5))
    axes[0].imshow(before, cmap='gray')
    axes[0].set_title('Before')
    axes[0].axis('off')
    axes[1].imshow(after, cmap='gray')
    axes[1].set_title('After')
    axes[1].axis('off')

    plt.tight_layout()
    plt.savefig("${output_file}", dpi=100, bbox_inches='tight')
    plt.close()

    print(f"Saved: ${output_file}")

except Exception as e:
    print(f"QC comparison failed: {e}", file=sys.stderr)
    sys.exit(1)
PYEOF
}

# =============================================================================
# RUN MANAGEMENT
# =============================================================================

# Generate unique run ID
generate_run_id() {
    date +"%Y%m%d_%H%M%S"
}

# Create run directory structure
create_run_structure() {
    local run_id=$1
    local run_dir="${PROJECT_DIR}/slurm/runs/${run_id}"

    mkdir -p "${run_dir}/logs"
    mkdir -p "${run_dir}/qc/correction"
    mkdir -p "${run_dir}/qc/stitch"
    mkdir -p "${run_dir}/qc/decon"
    mkdir -p "${run_dir}/qc/edf"
    mkdir -p "${run_dir}/qc/registration"
    mkdir -p "${run_dir}/summaries"

    # Create run info file
    cat > "${run_dir}/run_info.txt" << EOF
KINTSUGI Processing Run
=======================
Run ID:      ${run_id}
Project:     ${PROJECT_DIR}
Started:     $(date)
User:        ${USER}
Host:        $(hostname)
SLURM Job:   ${SLURM_JOB_ID:-N/A}

Configuration:
  Cycles:    ${START_CYCLE:-?} - ${END_CYCLE:-?}
  Channels:  ${START_CHANNEL:-?} - ${END_CHANNEL:-?}
  Format:    ${OUTPUT_FORMAT:-tiff}
  GPU Type:  ${GPU_TYPE:-auto}
EOF

    echo "${run_dir}"
}

# Finalize run with summary
finalize_run() {
    local run_id=$1
    local run_dir="${PROJECT_DIR}/slurm/runs/${run_id}"

    # Generate final summary
    cat > "${run_dir}/summary.txt" << EOF
KINTSUGI Processing Summary
===========================
Run ID:      ${run_id}
Completed:   $(date)

Output Locations:
  Corrected:    ${PROJECT_DIR}/data/processed/corrected/
  Stitched:     ${PROJECT_DIR}/data/processed/stitched/
  Deconvolved:  ${PROJECT_DIR}/data/processed/deconvolved/
  EDF:          ${PROJECT_DIR}/data/processed/edf/
  Registered:   ${PROJECT_DIR}/data/processed/registered/

QC Images:      ${run_dir}/qc/
Logs:           ${run_dir}/logs/
Step Summaries: ${run_dir}/summaries/

Data Summary:
$(summarize_directory "${PROJECT_DIR}/data/processed" "Processed Data")

EOF

    log_info "Run summary: ${run_dir}/summary.txt"
}
