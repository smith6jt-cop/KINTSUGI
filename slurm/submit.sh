#!/bin/bash
# =============================================================================
# KINTSUGI Universal Pipeline Submission
# Submit processing jobs for any KINTSUGI project
# =============================================================================

set -e

KINTSUGI_SLURM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KINTSUGI_DIR="$(dirname "${KINTSUGI_SLURM}")"

# Source utilities
source "${KINTSUGI_SLURM}/utils.sh"

# Defaults
STEPS="all"
DRY_RUN=false
CYCLES_OVERRIDE=""
PROJECT_DIR=""
CONFIG_FILE=""
RUN_ID=""

usage() {
    cat << EOF
KINTSUGI SLURM Pipeline Submission

Usage: $0 --project <PROJECT_DIR> [OPTIONS]

Required:
  --project DIR       Path to KINTSUGI project directory
                      (must contain data/ and notebooks/ subdirectories)

Options:
  --config FILE       Custom config file (default: PROJECT_DIR/slurm/config.sh)
  --steps STEPS       Comma-separated: correction,stitch,decon,edf,all (default: all)
  --cycles RANGE      Override cycles: '1-7' or '1,2,5' (default: from config)
  --dry-run           Show commands without submitting
  --help              Show this help

Examples:
  # Full pipeline for a project
  $0 --project /blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L

  # Just deconvolution for cycles 1-3
  $0 --project ~/my_project --steps decon --cycles 1-3

  # Preview without submitting
  $0 --project ~/my_project --dry-run

  # Multiple projects (run separately)
  $0 --project /path/to/project1 &
  $0 --project /path/to/project2 &
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project|-p)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --config|-c)
            CONFIG_FILE="$2"
            shift 2
            ;;
        --steps|-s)
            STEPS="$2"
            shift 2
            ;;
        --cycles)
            CYCLES_OVERRIDE="$2"
            shift 2
            ;;
        --dry-run|-n)
            DRY_RUN=true
            shift
            ;;
        --help|-h)
            usage
            exit 0
            ;;
        *)
            echo "ERROR: Unknown option: $1"
            usage
            exit 1
            ;;
    esac
done

# Validate project directory
if [ -z "${PROJECT_DIR}" ]; then
    echo "ERROR: --project is required"
    usage
    exit 1
fi

PROJECT_DIR="$(realpath "${PROJECT_DIR}")"

if [ ! -d "${PROJECT_DIR}" ]; then
    echo "ERROR: Project directory does not exist: ${PROJECT_DIR}"
    exit 1
fi

if [ ! -d "${PROJECT_DIR}/data" ]; then
    echo "ERROR: Not a valid KINTSUGI project (missing data/): ${PROJECT_DIR}"
    exit 1
fi

# Find or create config
if [ -z "${CONFIG_FILE}" ]; then
    CONFIG_FILE="${PROJECT_DIR}/slurm/config.sh"
fi

if [ ! -f "${CONFIG_FILE}" ]; then
    echo "Config not found. Creating default config at: ${CONFIG_FILE}"
    mkdir -p "$(dirname "${CONFIG_FILE}")"

    # Generate config from template
    PROJECT_NAME="$(basename "${PROJECT_DIR}")"
    sed "s|PROJECT_NAME=\"my_project\"|PROJECT_NAME=\"${PROJECT_NAME}\"|g; \
         s|/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/\${PROJECT_NAME}|${PROJECT_DIR}|g" \
         "${KINTSUGI_SLURM}/config_template.sh" > "${CONFIG_FILE}"

    echo "Please review and edit: ${CONFIG_FILE}"
    echo "Then re-run this command."
    exit 0
fi

# Source configuration
source "${CONFIG_FILE}"

# Override with command line
export PROJECT_DIR
export KINTSUGI_DIR
export NOTEBOOK_DIR="${PROJECT_DIR}/notebooks"
export DATA_DIR="${PROJECT_DIR}/data"

# Generate unique run ID and create run structure
RUN_ID=$(generate_run_id)
export RUN_ID
RUN_DIR=$(create_run_structure "${RUN_ID}")
export RUN_DIR
export LOG_DIR="${RUN_DIR}/logs"
export QC_DIR="${RUN_DIR}/qc"
export SUMMARY_DIR="${RUN_DIR}/summaries"

# Determine cycle array
if [ -n "${CYCLES_OVERRIDE}" ]; then
    ARRAY_SPEC="${CYCLES_OVERRIDE}"
else
    ARRAY_SPEC="${START_CYCLE}-${END_CYCLE}"
fi

PROJECT_NAME="$(basename "${PROJECT_DIR}")"

echo "============================================================"
echo "KINTSUGI Pipeline Submission"
echo "============================================================"
echo "Project:  ${PROJECT_NAME}"
echo "Path:     ${PROJECT_DIR}"
echo "Run ID:   ${RUN_ID}"
echo "Config:   ${CONFIG_FILE}"
echo "Cycles:   ${ARRAY_SPEC}"
echo "Steps:    ${STEPS}"
echo "Format:   ${OUTPUT_FORMAT}"
echo "Dry run:  ${DRY_RUN}"
echo "============================================================"
echo ""
echo "Output structure:"
echo "  Logs:       ${LOG_DIR}/"
echo "  QC Images:  ${QC_DIR}/"
echo "  Summaries:  ${SUMMARY_DIR}/"
echo "============================================================"

# Submit job function
submit_job() {
    local step_name=$1
    local script=$2
    local mem=$3
    local time=$4
    local dep_type=$5
    local dep_jobid=$6

    local job_name="kintsugi_${step_name}_${PROJECT_NAME}"

    # Ensure log and QC directories exist (sbatch validates paths at submission time)
    mkdir -p "${LOG_DIR}"
    mkdir -p "${QC_DIR}/${step_name}"

    local cmd="sbatch"
    cmd="${cmd} --job-name=${job_name}"
    cmd="${cmd} --output=${LOG_DIR}/${step_name}_%A_%a.out"
    cmd="${cmd} --error=${LOG_DIR}/${step_name}_%A_%a.err"
    cmd="${cmd} --partition=${PARTITION}"
    if [ -n "${QOS}" ]; then
        cmd="${cmd} --qos=${QOS}"
    fi
    cmd="${cmd} --account=${ACCOUNT}"
    cmd="${cmd} --nodes=1"
    cmd="${cmd} --ntasks=1"
    cmd="${cmd} --cpus-per-task=8"
    cmd="${cmd} --mem=${mem}gb"
    cmd="${cmd} --gpus=${GPU_TYPE}:${GPUS_PER_NODE}"
    cmd="${cmd} --time=${time}"
    cmd="${cmd} --array=${ARRAY_SPEC}"

    # Export all config variables including RUN_ID for unified logging
    cmd="${cmd} --export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},RUN_ID=${RUN_ID},QC_DIR=${QC_DIR}/${step_name}"

    if [ -n "${dep_jobid}" ]; then
        cmd="${cmd} --dependency=${dep_type}:${dep_jobid}"
    fi

    if [ -n "${EMAIL}" ]; then
        cmd="${cmd} --mail-user=${EMAIL} --mail-type=${MAIL_TYPE}"
    fi

    cmd="${cmd} ${script}"

    if [ "${DRY_RUN}" = true ]; then
        # Print to stderr so it's visible (stdout is captured for job ID)
        echo "[DRY RUN] ${cmd}" >&2
        echo "DRY_${step_name}_JOB"
    else
        echo "Submitting ${step_name}..." >&2
        result=$(${cmd} 2>&1)
        if [ $? -ne 0 ]; then
            echo "ERROR: sbatch failed: ${result}" >&2
            echo ""
            return 1
        fi
        jobid=$(echo ${result} | grep -oP '\d+$')
        echo "  Job ID: ${jobid}" >&2
        echo "${jobid}"
    fi
}

# Job IDs for dependencies
JOB_CORRECTION=""
JOB_STITCH=""
JOB_DECON=""
JOB_EDF=""

# Run steps
run_step() {
    local step=$1
    echo ""
    case ${step} in
        correction)
            echo "--- Step 1: Illumination Correction ---"
            JOB_CORRECTION=$(submit_job "correct" "${KINTSUGI_SLURM}/jobs/01_correction.sh" \
                "${MEM_CORRECTION}" "${TIME_CORRECTION}" "" "")
            ;;
        stitch)
            echo "--- Step 2: Stitching ---"
            JOB_STITCH=$(submit_job "stitch" "${KINTSUGI_SLURM}/jobs/02_stitching.sh" \
                "${MEM_STITCH}" "${TIME_STITCH}" "afterok" "${JOB_CORRECTION}")
            ;;
        decon)
            echo "--- Step 3: Deconvolution ---"
            JOB_DECON=$(submit_job "decon" "${KINTSUGI_SLURM}/jobs/03_deconvolution.sh" \
                "${MEM_DECON}" "${TIME_DECON}" "afterok" "${JOB_STITCH}")
            ;;
        edf)
            echo "--- Step 4: Extended Depth of Focus ---"
            JOB_EDF=$(submit_job "edf" "${KINTSUGI_SLURM}/jobs/04_edf.sh" \
                "${MEM_EDF}" "${TIME_EDF}" "afterok" "${JOB_DECON}")
            ;;
    esac
}

if [ "${STEPS}" = "all" ]; then
    run_step "correction"
    run_step "stitch"
    run_step "decon"
    run_step "edf"
else
    IFS=',' read -ra STEP_ARRAY <<< "${STEPS}"
    for step in "${STEP_ARRAY[@]}"; do
        run_step "${step}"
    done
fi

echo ""
echo "============================================================"
echo "Jobs submitted for: ${PROJECT_NAME}"
echo "============================================================"
echo ""
echo "Run ID:     ${RUN_ID}"
echo "Run Dir:    ${RUN_DIR}"
echo ""
echo "Monitor jobs:  squeue -u \$USER -n kintsugi_*_${PROJECT_NAME}"
echo "View logs:     ls ${LOG_DIR}/"
echo "View QC:       ls ${QC_DIR}/"
echo "View summary:  cat ${RUN_DIR}/run_info.txt"
echo ""
echo "After completion:"
echo "  All QC images: ${QC_DIR}/"
echo "  Summaries:     ${SUMMARY_DIR}/"
echo ""
