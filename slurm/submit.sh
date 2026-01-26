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

# =============================================================================
# TERMINAL LOGGING
# =============================================================================

# Colors for terminal output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Terminal logging functions
term_log() {
    local level=$1
    shift
    local timestamp=$(date '+%Y-%m-%d %H:%M:%S')
    case ${level} in
        INFO)
            echo -e "${BLUE}[${timestamp}]${NC} ${GREEN}[INFO]${NC} $*" >&2
            ;;
        WARN)
            echo -e "${BLUE}[${timestamp}]${NC} ${YELLOW}[WARN]${NC} $*" >&2
            ;;
        ERROR)
            echo -e "${BLUE}[${timestamp}]${NC} ${RED}[ERROR]${NC} $*" >&2
            ;;
        *)
            echo -e "${BLUE}[${timestamp}]${NC} $*" >&2
            ;;
    esac
}

term_info() { term_log INFO "$@"; }
term_warn() { term_log WARN "$@"; }
term_error() { term_log ERROR "$@"; }

# =============================================================================
# GPU RESOURCE CHECKING
# =============================================================================

# Check if a GPU type is available in a partition
# Returns: 0 if available, 1 if not available
check_gpu_availability() {
    local partition=$1
    local gpu_type=$2

    if [ "${DRY_RUN}" = true ]; then
        term_info "DRY RUN: Skipping GPU availability check for ${gpu_type} on ${partition}"
        return 0
    fi

    # Check if sinfo is available
    if ! command -v sinfo &> /dev/null; then
        term_warn "sinfo command not found - cannot verify GPU availability"
        return 0  # Proceed anyway, let sbatch handle it
    fi

    term_info "Checking availability of ${gpu_type} GPUs on partition ${partition}..."

    # Query SLURM for available GPUs of the specified type
    # Using sinfo to check if the partition has idle/mix nodes with the GPU
    local gres_info
    gres_info=$(sinfo -p "${partition}" -N -o "%N %G %t" 2>/dev/null | grep -i "${gpu_type}" | grep -E "(idle|mix)" || true)

    if [ -n "${gres_info}" ]; then
        local available_nodes=$(echo "${gres_info}" | wc -l)
        term_info "Found ${available_nodes} node(s) with ${gpu_type} GPUs available/mixed on ${partition}"
        return 0
    fi

    # Also check via squeue for pending jobs which might indicate resource contention
    local queue_count=0
    if command -v squeue &> /dev/null; then
        queue_count=$(squeue -p "${partition}" --gres="gpu:${gpu_type}" -h 2>/dev/null | wc -l || echo "0")

        if [ "${queue_count}" -gt 10 ]; then
            term_warn "High queue depth (${queue_count} jobs) for ${gpu_type} on ${partition} - jobs may wait"
        fi
    else
        term_warn "squeue command not found - skipping queue depth check for ${gpu_type} on ${partition}"
    fi

    # Check if any nodes with this GPU type exist at all
    local total_nodes
    total_nodes=$(sinfo -p "${partition}" -N -o "%N %G" 2>/dev/null | grep -i "${gpu_type}" | wc -l || echo "0")

    if [ "${total_nodes}" -eq 0 ]; then
        term_error "No nodes with ${gpu_type} GPUs found on partition ${partition}"
        return 1
    fi

    term_warn "All ${gpu_type} nodes on ${partition} appear busy - job will be queued"
    return 2
}

# Get current GPU resource status
show_gpu_status() {
    local partition=$1
    local gpu_type=$2

    if ! command -v sinfo &> /dev/null; then
        return
    fi

    term_info "Current GPU status for ${gpu_type} on ${partition}:"
    sinfo -p "${partition}" -N -o "  %N %G %t %C" 2>/dev/null | \
        grep -i "${gpu_type}" | head -10 >&2 || \
        echo "  No specific GPU info available" >&2
}

# =============================================================================
# RESOURCE CALCULATION
# =============================================================================

# Calculate optimal max concurrent cycles based on allocation limits
# Sets COMPUTED_MAX_CONCURRENT as output (global)
calculate_max_concurrent() {
    local cpus_per_task=${CPUS_PER_TASK:-4}
    local gpus_per_node=${GPUS_PER_NODE:-1}
    local alloc_cpus=${ALLOC_CPUS:-64}
    local alloc_mem=${ALLOC_MEM:-600}
    local alloc_gpus=${ALLOC_GPUS:-2}

    # Find max memory requirement across all steps
    local max_mem=${MEM_DECON:-96}
    [ "${MEM_STITCH:-64}" -gt "${max_mem}" ] && max_mem=${MEM_STITCH:-64}
    [ "${MEM_CORRECTION:-32}" -gt "${max_mem}" ] && max_mem=${MEM_CORRECTION:-32}
    [ "${MEM_EDF:-32}" -gt "${max_mem}" ] && max_mem=${MEM_EDF:-32}

    # Calculate limits from each resource
    local max_by_cpu=$((alloc_cpus / cpus_per_task))
    local max_by_mem=$((alloc_mem / max_mem))
    local max_by_gpu=$((alloc_gpus / gpus_per_node))

    # Find minimum (limiting factor)
    COMPUTED_MAX_CONCURRENT=${max_by_cpu}
    LIMITING_RESOURCE="CPUs"

    if [ ${max_by_mem} -lt ${COMPUTED_MAX_CONCURRENT} ]; then
        COMPUTED_MAX_CONCURRENT=${max_by_mem}
        LIMITING_RESOURCE="Memory"
    fi

    if [ ${max_by_gpu} -lt ${COMPUTED_MAX_CONCURRENT} ]; then
        COMPUTED_MAX_CONCURRENT=${max_by_gpu}
        LIMITING_RESOURCE="GPUs"
    fi

    # Store calculation details for display
    RESOURCE_CALC_DETAILS="CPU: ${alloc_cpus}/${cpus_per_task}=${max_by_cpu}, MEM: ${alloc_mem}/${max_mem}=${max_by_mem}, GPU: ${alloc_gpus}/${gpus_per_node}=${max_by_gpu}"

    # Validate we can run at least 1 job
    if [ ${COMPUTED_MAX_CONCURRENT} -lt 1 ]; then
        term_error "============================================================"
        term_error "RESOURCE VALIDATION FAILED"
        term_error "============================================================"
        term_error "Cannot fit even 1 job within allocation limits!"
        term_error ""
        term_error "Allocation limits:"
        term_error "  CPUs:   ${alloc_cpus}"
        term_error "  Memory: ${alloc_mem} GB"
        term_error "  GPUs:   ${alloc_gpus}"
        term_error ""
        term_error "Per-job requirements:"
        term_error "  CPUs:   ${cpus_per_task}"
        term_error "  Memory: ${max_mem} GB (max across steps)"
        term_error "  GPUs:   ${gpus_per_node}"
        term_error ""
        term_error "Options:"
        term_error "  1. Reduce per-job resources in config.sh"
        term_error "  2. Increase allocation limits (if possible)"
        term_error "  3. Request a larger allocation from your HPC admin"
        term_error "============================================================"
        return 1
    fi

    return 0
}

# Resolve MAX_CONCURRENT_CYCLES (auto or manual)
# Sets EFFECTIVE_MAX_CONCURRENT as output (global)
resolve_max_concurrent() {
    local configured=${MAX_CONCURRENT_CYCLES:-auto}

    # Calculate optimal value
    if ! calculate_max_concurrent; then
        return 1
    fi

    if [ "${configured}" = "auto" ] || [ -z "${configured}" ]; then
        EFFECTIVE_MAX_CONCURRENT=${COMPUTED_MAX_CONCURRENT}
        term_info "Auto-calculated max concurrent cycles: ${EFFECTIVE_MAX_CONCURRENT}"
        term_info "  Calculation: ${RESOURCE_CALC_DETAILS}"
        term_info "  Limiting factor: ${LIMITING_RESOURCE}"
    else
        # Manual override - validate it
        EFFECTIVE_MAX_CONCURRENT=${configured}
        if [ ${EFFECTIVE_MAX_CONCURRENT} -gt ${COMPUTED_MAX_CONCURRENT} ]; then
            term_warn "Manual MAX_CONCURRENT_CYCLES=${EFFECTIVE_MAX_CONCURRENT} exceeds safe limit of ${COMPUTED_MAX_CONCURRENT}"
            term_warn "  Calculation: ${RESOURCE_CALC_DETAILS}"
            term_warn "  This may cause job failures due to resource over-subscription!"
            term_warn ""
            if [ "${DRY_RUN}" != true ] && [ -t 0 ]; then
                read -p "Continue with potentially unsafe concurrency? [y/N] " -n 1 -r
                echo ""
                if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                    term_info "Submission cancelled. Set MAX_CONCURRENT_CYCLES=auto for safe defaults."
                    return 1
                fi
            fi
        else
            term_info "Using configured MAX_CONCURRENT_CYCLES=${EFFECTIVE_MAX_CONCURRENT} (safe limit: ${COMPUTED_MAX_CONCURRENT})"
        fi
    fi

    # Calculate actual resource usage with this concurrency
    local cpus_per_task=${CPUS_PER_TASK:-4}
    local gpus_per_node=${GPUS_PER_NODE:-1}
    local max_mem=${MEM_DECON:-96}

    TOTAL_CPUS_USED=$((EFFECTIVE_MAX_CONCURRENT * cpus_per_task))
    TOTAL_MEM_USED=$((EFFECTIVE_MAX_CONCURRENT * max_mem))
    TOTAL_GPUS_USED=$((EFFECTIVE_MAX_CONCURRENT * gpus_per_node))

    return 0
}

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

echo ""
term_info "============================================================"
term_info "KINTSUGI Pipeline Submission"
term_info "============================================================"
term_info "Project:  ${PROJECT_NAME}"
term_info "Path:     ${PROJECT_DIR}"
term_info "Run ID:   ${RUN_ID}"
term_info "Config:   ${CONFIG_FILE}"
term_info "Cycles:   ${ARRAY_SPEC}"
term_info "Steps:    ${STEPS}"
term_info "Format:   ${OUTPUT_FORMAT}"
term_info "Dry run:  ${DRY_RUN}"
term_info "------------------------------------------------------------"
term_info "GPU Configuration:"
term_info "  Primary:  ${GPU_TYPE} on ${PARTITION}"
if [ -n "${GPU_TYPE_FALLBACK}" ] && [ -n "${PARTITION_FALLBACK}" ]; then
    term_info "  Fallback: ${GPU_TYPE_FALLBACK} on ${PARTITION_FALLBACK}"
else
    term_info "  Fallback: (not configured)"
fi
term_info "------------------------------------------------------------"

# Calculate optimal concurrency based on allocation limits
term_info "Calculating resource allocation..."
if ! resolve_max_concurrent; then
    exit 1
fi

term_info "------------------------------------------------------------"
term_info "Resource Allocation:"
term_info "  Allocation limits: ${ALLOC_CPUS:-64} CPUs, ${ALLOC_MEM:-600}GB mem, ${ALLOC_GPUS:-2} GPUs"
term_info "  Per-job request:   ${CPUS_PER_TASK:-4} CPUs, ${MEM_DECON:-96}GB mem, ${GPUS_PER_NODE:-1} GPU"
term_info "  Max concurrent:    ${EFFECTIVE_MAX_CONCURRENT} cycles (limited by ${LIMITING_RESOURCE})"
term_info "  Total usage:       ${TOTAL_CPUS_USED} CPUs, ${TOTAL_MEM_USED}GB mem, ${TOTAL_GPUS_USED} GPUs"
term_info "============================================================"
echo ""

# Initial resource availability check
term_info "Checking cluster resource availability..."
if ! check_gpu_availability "${PARTITION}" "${GPU_TYPE}"; then
    term_warn "Primary GPU resources may not be immediately available"
    if [ -n "${GPU_TYPE_FALLBACK}" ] && [ -n "${PARTITION_FALLBACK}" ]; then
        if check_gpu_availability "${PARTITION_FALLBACK}" "${GPU_TYPE_FALLBACK}"; then
            term_info "Fallback GPU resources are available - will use if primary fails"
        else
            term_error "Neither primary nor fallback GPU resources appear available"
            term_error "Jobs will likely fail or wait indefinitely"
            if [ "${DRY_RUN}" != true ]; then
                # Only prompt if running interactively (stdin is a TTY)
                if [ -t 0 ]; then
                    echo ""
                    read -p "Continue anyway? [y/N] " -n 1 -r
                    echo ""
                    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                        term_info "Submission cancelled"
                        exit 1
                    fi
                else
                    # Non-interactive mode: abort by default for safety
                    term_error "Non-interactive mode: aborting due to unavailable resources"
                    term_error "Use --dry-run to preview or ensure resources are available"
                    exit 1
                fi
            fi
        fi
    else
        term_warn "No fallback GPU configured - jobs will wait for primary resources"
    fi
else
    term_info "Primary GPU resources appear available"
fi

term_info "Output structure:"
term_info "  Logs:       ${LOG_DIR}/"
term_info "  QC Images:  ${QC_DIR}/"
term_info "  Summaries:  ${SUMMARY_DIR}/"
term_info "============================================================"

# Build sbatch command array with specified GPU settings
# Sets SBATCH_CMD array as output (global)
build_sbatch_cmd_array() {
    local step_name=$1
    local script=$2
    local mem=$3
    local time=$4
    local dep_type=$5
    local dep_jobid=$6
    local use_partition=$7
    local use_gpu_type=$8
    local use_qos=$9

    local job_name="kintsugi_${step_name}_${PROJECT_NAME}"

    SBATCH_CMD=(
        sbatch
        "--job-name=${job_name}"
        "--output=${LOG_DIR}/${step_name}_%A_%a.out"
        "--error=${LOG_DIR}/${step_name}_%A_%a.err"
        "--partition=${use_partition}"
    )

    if [ -n "${use_qos}" ]; then
        SBATCH_CMD+=("--qos=${use_qos}")
    fi

    SBATCH_CMD+=(
        "--account=${ACCOUNT}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${CPUS_PER_TASK:-4}"
        "--mem=${mem}gb"
        "--gpus=${use_gpu_type}:${GPUS_PER_NODE:-1}"
        "--time=${time}"
        "--array=${ARRAY_SPEC}%${EFFECTIVE_MAX_CONCURRENT}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},RUN_ID=${RUN_ID},QC_DIR=${QC_DIR}/${step_name}"
    )

    if [ -n "${dep_jobid}" ]; then
        # Validate dependency job ID is numeric to prevent "Badly formatted array jobid" errors
        if ! [[ "${dep_jobid}" =~ ^[0-9]+$ ]]; then
            echo "ERROR: Invalid dependency job ID '${dep_jobid}' - must be numeric" >&2
            return 1
        fi
        SBATCH_CMD+=("--dependency=${dep_type}:${dep_jobid}")
    fi

    if [ -n "${EMAIL}" ]; then
        SBATCH_CMD+=("--mail-user=${EMAIL}" "--mail-type=${MAIL_TYPE}")
    fi

    SBATCH_CMD+=("${script}")
}

# Format command array as string for display
format_cmd_for_display() {
    printf '%q ' "${SBATCH_CMD[@]}"
}

# Submit job function with fallback support
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

    term_info "Preparing job submission for ${step_name}..."
    term_info "  Partition: ${PARTITION}, GPU: ${GPU_TYPE}, Memory: ${mem}GB, Time: ${time}"

    # Determine which GPU configuration to use
    local use_partition="${PARTITION}"
    local use_gpu_type="${GPU_TYPE}"
    local use_qos="${QOS}"
    local using_fallback=false

    # Check primary GPU availability
    if ! check_gpu_availability "${PARTITION}" "${GPU_TYPE}"; then
        term_warn "Primary GPU (${GPU_TYPE}) not available on partition ${PARTITION}"

        # Check if fallback is configured
        if [ -n "${GPU_TYPE_FALLBACK}" ] && [ -n "${PARTITION_FALLBACK}" ]; then
            term_info "Checking fallback: ${GPU_TYPE_FALLBACK} on partition ${PARTITION_FALLBACK}..."
            if check_gpu_availability "${PARTITION_FALLBACK}" "${GPU_TYPE_FALLBACK}"; then
                term_info "Switching to fallback GPU configuration"
                use_partition="${PARTITION_FALLBACK}"
                use_gpu_type="${GPU_TYPE_FALLBACK}"
                use_qos="${QOS_FALLBACK}"
                using_fallback=true
            else
                term_error "Fallback GPU (${GPU_TYPE_FALLBACK}) also not available on ${PARTITION_FALLBACK}"
                term_error "No GPU resources available. Aborting submission."
                show_gpu_status "${PARTITION}" "${GPU_TYPE}"
                return 1
            fi
        else
            term_error "No fallback GPU configured and primary GPU unavailable. Aborting."
            term_error "Configure GPU_TYPE_FALLBACK and PARTITION_FALLBACK in config.sh to enable fallback."
            return 1
        fi
    fi

    # Build command with selected GPU configuration
    build_sbatch_cmd_array "${step_name}" "${script}" "${mem}" "${time}" "${dep_type}" "${dep_jobid}" \
        "${use_partition}" "${use_gpu_type}" "${use_qos}"

    if [ "${DRY_RUN}" = true ]; then
        # Print to stderr so it's visible (stdout is captured for job ID)
        if [ "${using_fallback}" = true ]; then
            term_info "DRY RUN - would execute (using fallback configuration):"
        else
            term_info "DRY RUN - would execute:"
        fi
        echo "  $(format_cmd_for_display)" >&2
        echo "DRY_${step_name}_JOB"
        return 0
    fi

    # Submit job
    if [ "${using_fallback}" = true ]; then
        term_info "Submitting ${step_name} to ${use_partition} with ${use_gpu_type} GPU (fallback)..."
    else
        term_info "Submitting ${step_name} to ${use_partition} with ${use_gpu_type} GPU..."
    fi
    local result
    local exit_code
    result=$("${SBATCH_CMD[@]}" 2>&1)
    exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        local jobid
        # Extract job ID from "Submitted batch job XXXXXXXX" message
        # More robust pattern that handles array job output and edge cases
        jobid=$(echo "${result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)

        # Fallback to simpler pattern if grep -P not available
        if [ -z "${jobid}" ]; then
            jobid=$(echo "${result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
        fi

        # Validate job ID is numeric and non-empty
        if [ -z "${jobid}" ] || ! [[ "${jobid}" =~ ^[0-9]+$ ]]; then
            term_error "Failed to extract job ID from sbatch output:"
            term_error "  Output: ${result}"
            return 1
        fi

        term_info "Successfully submitted ${step_name} - Job ID: ${jobid}"
        echo "${jobid}"
        return 0
    fi

    # Primary submission failed - preserve error for reporting
    local primary_error="${result}"
    term_warn "Primary submission failed: ${primary_error}"

    # Common resource-related error patterns
    if echo "${primary_error}" | grep -qiE "(invalid|unavailable|not available|no nodes|cannot satisfy|partition.*invalid|gres.*invalid|constraint.*invalid)"; then
        term_warn "Submission failed due to resource unavailability"

        # Try fallback if configured and not already using it
        if [ "${using_fallback}" = false ] && [ -n "${GPU_TYPE_FALLBACK}" ] && [ -n "${PARTITION_FALLBACK}" ]; then
            term_info "Attempting fallback submission to ${PARTITION_FALLBACK} with ${GPU_TYPE_FALLBACK} GPU..."

            build_sbatch_cmd_array "${step_name}" "${script}" "${mem}" "${time}" "${dep_type}" "${dep_jobid}" \
                "${PARTITION_FALLBACK}" "${GPU_TYPE_FALLBACK}" "${QOS_FALLBACK}"

            local fallback_result
            fallback_result=$("${SBATCH_CMD[@]}" 2>&1)
            exit_code=$?

            if [ ${exit_code} -eq 0 ]; then
                local jobid
                # Extract job ID using same robust pattern as primary submission
                jobid=$(echo "${fallback_result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
                if [ -z "${jobid}" ]; then
                    jobid=$(echo "${fallback_result}" | grep "Submitted batch job" | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
                fi

                if [ -z "${jobid}" ] || ! [[ "${jobid}" =~ ^[0-9]+$ ]]; then
                    term_error "Failed to extract job ID from fallback sbatch output:"
                    term_error "  Output: ${fallback_result}"
                    return 1
                fi

                term_info "Fallback submission successful - Job ID: ${jobid}"
                term_info "NOTE: Using ${GPU_TYPE_FALLBACK} GPU instead of ${GPU_TYPE}"
                echo "${jobid}"
                return 0
            else
                term_error "Fallback submission also failed: ${fallback_result}"
            fi
        elif [ "${using_fallback}" = true ]; then
            term_error "Already using fallback configuration, no further options available"
        else
            term_error "No fallback GPU configured. Set GPU_TYPE_FALLBACK and PARTITION_FALLBACK in config.sh"
        fi
    fi

    # Submission failed
    term_error "Job submission failed for ${step_name}"
    term_error "Error: $(echo "${primary_error}" | head -1)"
    term_error ""
    term_error "Troubleshooting:"
    term_error "  1. Check available partitions: sinfo"
    term_error "  2. Check your account access: sacctmgr show user \$USER"
    term_error "  3. Check GPU availability: sinfo -p ${PARTITION} -o '%N %G %t'"
    term_error "  4. Check queue status: squeue -p ${PARTITION}"
    term_error ""

    return 1
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
            if [ -z "${JOB_CORRECTION}" ]; then
                term_error "Failed to submit correction job - aborting pipeline"
                exit 1
            fi
            ;;
        stitch)
            echo "--- Step 2: Stitching ---"
            JOB_STITCH=$(submit_job "stitch" "${KINTSUGI_SLURM}/jobs/02_stitching.sh" \
                "${MEM_STITCH}" "${TIME_STITCH}" "afterok" "${JOB_CORRECTION}")
            if [ -z "${JOB_STITCH}" ]; then
                term_error "Failed to submit stitching job - aborting pipeline"
                exit 1
            fi
            ;;
        decon)
            echo "--- Step 3: Deconvolution ---"
            JOB_DECON=$(submit_job "decon" "${KINTSUGI_SLURM}/jobs/03_deconvolution.sh" \
                "${MEM_DECON}" "${TIME_DECON}" "afterok" "${JOB_STITCH}")
            if [ -z "${JOB_DECON}" ]; then
                term_error "Failed to submit deconvolution job - aborting pipeline"
                exit 1
            fi
            ;;
        edf)
            echo "--- Step 4: Extended Depth of Focus ---"
            JOB_EDF=$(submit_job "edf" "${KINTSUGI_SLURM}/jobs/04_edf.sh" \
                "${MEM_EDF}" "${TIME_EDF}" "afterok" "${JOB_DECON}")
            if [ -z "${JOB_EDF}" ]; then
                term_error "Failed to submit EDF job - aborting pipeline"
                exit 1
            fi
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
term_info "============================================================"
term_info "Jobs submitted for: ${PROJECT_NAME}"
term_info "============================================================"
echo ""
term_info "Run ID:     ${RUN_ID}"
term_info "Run Dir:    ${RUN_DIR}"
echo ""
term_info "Monitor jobs:  squeue -u \$USER -n kintsugi_*_${PROJECT_NAME}"
term_info "View logs:     ls ${LOG_DIR}/"
term_info "View QC:       ls ${QC_DIR}/"
term_info "View summary:  cat ${RUN_DIR}/run_info.txt"
echo ""
term_info "After completion:"
term_info "  All QC images: ${QC_DIR}/"
term_info "  Summaries:     ${SUMMARY_DIR}/"
echo ""

# Final status summary
if [ "${DRY_RUN}" = true ]; then
    term_info "DRY RUN completed - no jobs were actually submitted"
else
    term_info "All jobs submitted successfully"
    term_info "Use 'squeue -u \$USER' to monitor job status"
fi
