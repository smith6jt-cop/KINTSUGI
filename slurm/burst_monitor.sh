#!/bin/bash
# =============================================================================
# KINTSUGI Job Auto-Promotion Monitor
# Monitors pending jobs and promotes them to better resources when available
# =============================================================================
#
# This script implements auto-promotion for two job types:
#   1. Burst -> Allocated: Preemptible GPU jobs to guaranteed GPU resources
#   2. CPU -> GPU: CPU-only jobs to GPU resources when GPUs become available
#
# When allocated GPU resources become available while jobs are pending:
#   1. Detects pending burst jobs and CPU jobs
#   2. Checks allocated partition availability
#   3. Cancels pending job and resubmits to GPU partition
#
# Promotion priority: Burst jobs first (already GPU-ready), then CPU jobs
#
# Usage:
#   ./burst_monitor.sh --project /path/to/project [--run-id RUN_ID] [--dry-run]
#
# The monitor runs as a background process and exits when:
#   - All burst/CPU jobs complete (success)
#   - No pending jobs remain
#   - User cancels (Ctrl+C)
#
# =============================================================================

set -e

# Default values
PROJECT_DIR=""
RUN_ID=""
DRY_RUN=false
CHECK_INTERVAL=60  # seconds between checks
MAX_RUNTIME=86400  # 24 hours max runtime

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m'

log_info() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ${GREEN}[MONITOR]${NC} $*"
}

log_warn() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ${YELLOW}[MONITOR]${NC} $*"
}

log_error() {
    echo -e "${BLUE}[$(date '+%Y-%m-%d %H:%M:%S')]${NC} ${RED}[MONITOR]${NC} $*"
}

usage() {
    cat << EOF
KINTSUGI Burst-to-Allocated Auto-Promotion Monitor

Usage: $0 --project <PROJECT_DIR> [OPTIONS]

Required:
  --project DIR       Path to KINTSUGI project directory

Options:
  --run-id ID         Specific run ID to monitor (default: all kintsugi_* jobs)
  --interval SECS     Check interval in seconds (default: 60)
  --dry-run           Show what would be done without taking action
  --help              Show this help

Examples:
  # Monitor all kintsugi burst jobs for a project
  $0 --project /path/to/project

  # Monitor specific run
  $0 --project /path/to/project --run-id 20260204_130539

  # Preview mode
  $0 --project /path/to/project --dry-run
EOF
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --project|-p)
            PROJECT_DIR="$2"
            shift 2
            ;;
        --run-id|-r)
            RUN_ID="$2"
            shift 2
            ;;
        --interval|-i)
            CHECK_INTERVAL="$2"
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
    log_error "--project is required"
    usage
    exit 1
fi

PROJECT_DIR="$(realpath "${PROJECT_DIR}")"

if [ ! -d "${PROJECT_DIR}" ]; then
    log_error "Project directory does not exist: ${PROJECT_DIR}"
    exit 1
fi

# Source project config
CONFIG_FILE="${PROJECT_DIR}/slurm/config.sh"
if [ ! -f "${CONFIG_FILE}" ]; then
    log_error "Config not found: ${CONFIG_FILE}"
    exit 1
fi

source "${CONFIG_FILE}"

PROJECT_NAME="$(basename "${PROJECT_DIR}")"

# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

# Get pending burst jobs for this project
# Returns: job_id step_name cycle_id (tab-separated)
get_pending_burst_jobs() {
    local job_pattern="kintsugi_*_burst_${PROJECT_NAME}"

    if [ -n "${RUN_ID}" ]; then
        # If run ID specified, filter by it
        job_pattern="kintsugi_*_burst_${PROJECT_NAME}"
    fi

    # Get pending burst jobs: job_id, job_name, state
    squeue -u "${USER}" --name="${job_pattern}" --state=PD -h -o "%i %j %t" 2>/dev/null || true
}

# Get pending CPU jobs for this project (these are also candidates for promotion)
get_pending_cpu_jobs() {
    local job_pattern="kintsugi_*_cpu_${PROJECT_NAME}"
    squeue -u "${USER}" --name="${job_pattern}" --state=PD -h -o "%i %j %t" 2>/dev/null || true
}

# Check if allocated GPU partition has idle resources
# Returns: number of idle nodes (0 if none)
check_allocated_availability() {
    local partition="${PARTITION:-hpg-b200}"

    # Get node states: A/I/O/T format (Allocated/Idle/Other/Total)
    local node_info
    node_info=$(sinfo -p "${partition}" -h -o "%F" 2>/dev/null | head -1)

    if [ -z "${node_info}" ]; then
        echo "0"
        return
    fi

    # Extract idle count (second field in A/I/O/T)
    local idle_nodes
    idle_nodes=$(echo "${node_info}" | cut -d'/' -f2)

    echo "${idle_nodes:-0}"
}

# Extract step name from job name (e.g., "kintsugi_stitch_burst_project" -> "stitch")
extract_step_name() {
    local job_name=$1
    echo "${job_name}" | sed -E 's/kintsugi_([^_]+)_.*/\1/'
}

# Get the original sbatch parameters from a running/pending job
# This uses scontrol to retrieve job details
get_job_details() {
    local job_id=$1
    scontrol show job "${job_id}" 2>/dev/null
}

# Cancel a job
cancel_job() {
    local job_id=$1

    if [ "${DRY_RUN}" = true ]; then
        log_info "DRY RUN: Would cancel job ${job_id}"
        return 0
    fi

    scancel "${job_id}" 2>/dev/null
    return $?
}

# Resubmit a job to allocated partition
# Parameters: job_id, step_name
promote_job() {
    local job_id=$1
    local job_name=$2
    local step_name
    step_name=$(extract_step_name "${job_name}")

    # Determine the correct job script
    local script
    case ${step_name} in
        correct)
            script="${KINTSUGI_DIR}/slurm/jobs/01_correction.sh"
            local mem=${MEM_CORRECTION:-64}
            local time=${TIME_CORRECTION:-04:00:00}
            ;;
        stitch)
            script="${KINTSUGI_DIR}/slurm/jobs/02_stitching.sh"
            local mem=${MEM_STITCH:-128}
            local time=${TIME_STITCH:-06:00:00}
            ;;
        decon)
            script="${KINTSUGI_DIR}/slurm/jobs/03_deconvolution.sh"
            local mem=${MEM_DECON:-180}
            local time=${TIME_DECON:-08:00:00}
            ;;
        edf)
            script="${KINTSUGI_DIR}/slurm/jobs/04_edf.sh"
            local mem=${MEM_EDF:-64}
            local time=${TIME_EDF:-02:00:00}
            ;;
        *)
            log_warn "Unknown step: ${step_name} - skipping promotion"
            return 1
            ;;
    esac

    # Get job array info
    local job_details
    job_details=$(get_job_details "${job_id}")

    # Extract array task info (if this is an array job)
    local array_spec=""
    if echo "${job_details}" | grep -q "ArrayTaskId"; then
        local array_task
        array_task=$(echo "${job_details}" | grep "ArrayTaskId" | sed 's/.*ArrayTaskId=//' | cut -d' ' -f1)
        array_spec="--array=${array_task}"
    fi

    log_info "Promoting ${job_name} (${job_id}) to allocated partition ${PARTITION}"

    # Build new sbatch command for allocated partition
    local new_job_name="${job_name/_burst_/_promoted_}"

    local sbatch_cmd=(
        sbatch
        "--job-name=${new_job_name}"
        "--partition=${PARTITION}"
        "--account=${ACCOUNT}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${CPUS_PER_TASK:-8}"
        "--mem=${mem}gb"
        "--gpus=${GPUS_PER_NODE:-1}"
        "--time=${time}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},KINTSUGI_DEVICE_MODE=gpu"
    )

    if [ -n "${array_spec}" ]; then
        sbatch_cmd+=("${array_spec}")
    fi

    if [ -n "${QOS}" ]; then
        sbatch_cmd+=("--qos=${QOS}")
    fi

    sbatch_cmd+=("${script}")

    if [ "${DRY_RUN}" = true ]; then
        log_info "DRY RUN: Would execute:"
        echo "  ${sbatch_cmd[*]}"
        return 0
    fi

    # Cancel the burst job first
    log_info "Cancelling burst job ${job_id}..."
    if ! cancel_job "${job_id}"; then
        log_warn "Failed to cancel job ${job_id} - may have already started"
        return 1
    fi

    # Wait a moment for cancel to propagate
    sleep 2

    # Submit to allocated partition
    log_info "Submitting to allocated partition..."
    local result
    result=$("${sbatch_cmd[@]}" 2>&1)
    local exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        local new_job_id
        new_job_id=$(echo "${result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
        if [ -z "${new_job_id}" ]; then
            new_job_id=$(echo "${result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
        fi
        log_info "SUCCESS: Promoted ${job_id} -> ${new_job_id} on ${PARTITION}"
        return 0
    else
        log_error "Failed to submit promoted job: ${result}"
        return 1
    fi
}

# Promote a CPU job to GPU when GPU resources become available
# Parameters: job_id, job_name
promote_cpu_to_gpu() {
    local job_id=$1
    local job_name=$2
    local step_name
    step_name=$(extract_step_name "${job_name}")

    # Determine the correct job script and resources
    local script mem time
    case ${step_name} in
        correct)
            script="${KINTSUGI_DIR}/slurm/jobs/01_correction.sh"
            mem=${MEM_CORRECTION:-64}
            time=${TIME_CORRECTION:-04:00:00}
            ;;
        stitch)
            script="${KINTSUGI_DIR}/slurm/jobs/02_stitching.sh"
            mem=${MEM_STITCH:-128}
            time=${TIME_STITCH:-06:00:00}
            ;;
        decon)
            script="${KINTSUGI_DIR}/slurm/jobs/03_deconvolution.sh"
            mem=${MEM_DECON:-180}
            time=${TIME_DECON:-08:00:00}
            ;;
        edf)
            script="${KINTSUGI_DIR}/slurm/jobs/04_edf.sh"
            mem=${MEM_EDF:-64}
            time=${TIME_EDF:-02:00:00}
            ;;
        *)
            log_warn "Unknown step: ${step_name} - skipping CPU->GPU promotion"
            return 1
            ;;
    esac

    # Get job array info
    local job_details
    job_details=$(get_job_details "${job_id}")

    # Extract array task info (if this is an array job)
    local array_spec=""
    if echo "${job_details}" | grep -q "ArrayTaskId"; then
        local array_task
        array_task=$(echo "${job_details}" | grep "ArrayTaskId" | sed 's/.*ArrayTaskId=//' | cut -d' ' -f1)
        array_spec="--array=${array_task}"
    fi

    log_info "Promoting CPU job ${job_name} (${job_id}) to GPU partition ${PARTITION}"

    # Build new sbatch command for GPU execution
    local new_job_name="${job_name/_cpu_/_gpu_promoted_}"

    local sbatch_cmd=(
        sbatch
        "--job-name=${new_job_name}"
        "--partition=${PARTITION}"
        "--account=${ACCOUNT}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${CPUS_PER_TASK:-8}"
        "--mem=${mem}gb"
        "--gpus=${GPUS_PER_NODE:-1}"
        "--time=${time}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},KINTSUGI_DEVICE_MODE=gpu"
    )

    if [ -n "${array_spec}" ]; then
        sbatch_cmd+=("${array_spec}")
    fi

    if [ -n "${QOS}" ]; then
        sbatch_cmd+=("--qos=${QOS}")
    fi

    sbatch_cmd+=("${script}")

    if [ "${DRY_RUN}" = true ]; then
        log_info "DRY RUN: Would execute (CPU->GPU promotion):"
        echo "  ${sbatch_cmd[*]}"
        return 0
    fi

    # Cancel the CPU job first
    log_info "Cancelling CPU job ${job_id}..."
    if ! cancel_job "${job_id}"; then
        log_warn "Failed to cancel job ${job_id} - may have already started"
        return 1
    fi

    # Wait a moment for cancel to propagate
    sleep 2

    # Submit to GPU partition
    log_info "Submitting to GPU partition..."
    local result
    result=$("${sbatch_cmd[@]}" 2>&1)
    local exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        local new_job_id
        new_job_id=$(echo "${result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
        if [ -z "${new_job_id}" ]; then
            new_job_id=$(echo "${result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
        fi
        log_info "SUCCESS: CPU ${job_id} -> GPU ${new_job_id} on ${PARTITION}"
        return 0
    else
        log_error "Failed to submit GPU job: ${result}"
        return 1
    fi
}

# =============================================================================
# MAIN MONITORING LOOP
# =============================================================================

log_info "============================================================"
log_info "KINTSUGI Burst-to-Allocated Auto-Promotion Monitor"
log_info "============================================================"
log_info "Project:           ${PROJECT_NAME}"
log_info "Allocated partition: ${PARTITION}"
log_info "Fallback partition:  ${PARTITION_FALLBACK:-none}"
log_info "Check interval:    ${CHECK_INTERVAL}s"
log_info "Dry run:           ${DRY_RUN}"
log_info "============================================================"
log_info ""
log_info "Monitoring for pending burst jobs..."
log_info "Press Ctrl+C to stop"
log_info ""

start_time=$(date +%s)
promoted_count=0
burst_promoted_count=0
cpu_promoted_count=0

# Trap Ctrl+C
trap 'log_info "Monitor stopped by user"; exit 0' INT TERM

while true; do
    current_time=$(date +%s)
    elapsed=$((current_time - start_time))

    # Check max runtime
    if [ ${elapsed} -gt ${MAX_RUNTIME} ]; then
        log_warn "Max runtime (${MAX_RUNTIME}s) exceeded - stopping monitor"
        break
    fi

    # Get pending burst jobs
    pending_burst=$(get_pending_burst_jobs)
    pending_cpu=$(get_pending_cpu_jobs)

    # Count pending jobs (filter empty lines, handle empty output)
    if [ -n "${pending_burst}" ]; then
        burst_count=$(echo "${pending_burst}" | grep -c "^[0-9]" 2>/dev/null) || burst_count=0
    else
        burst_count=0
    fi
    if [ -n "${pending_cpu}" ]; then
        cpu_count=$(echo "${pending_cpu}" | grep -c "^[0-9]" 2>/dev/null) || cpu_count=0
    else
        cpu_count=0
    fi

    if [ "${burst_count}" -eq 0 ] && [ "${cpu_count}" -eq 0 ]; then
        log_info "No pending burst/CPU jobs found - all jobs running or complete"
        break
    fi

    log_info "Found ${burst_count} pending burst jobs, ${cpu_count} pending CPU jobs"

    # Check if allocated resources are available
    idle_nodes=$(check_allocated_availability)

    if [ "${idle_nodes}" -gt 0 ]; then
        log_info "Allocated partition ${PARTITION} has ${idle_nodes} idle node(s) - attempting promotion"

        # Process burst jobs first (they're preemptible GPU jobs -> guaranteed GPU)
        if [ -n "${pending_burst}" ]; then
            log_info "Processing pending burst jobs (burst -> allocated)..."
            while IFS= read -r line; do
                if [ -z "${line}" ]; then
                    continue
                fi

                job_id=$(echo "${line}" | awk '{print $1}')
                job_name=$(echo "${line}" | awk '{print $2}')

                if [ -n "${job_id}" ] && [ -n "${job_name}" ]; then
                    if promote_job "${job_id}" "${job_name}"; then
                        ((burst_promoted_count++))
                        ((promoted_count++))
                    fi

                    # Re-check availability after each promotion
                    idle_nodes=$(check_allocated_availability)
                    if [ "${idle_nodes}" -eq 0 ]; then
                        log_info "No more idle nodes - waiting for next cycle"
                        break
                    fi
                fi
            done <<< "${pending_burst}"
        fi

        # Process CPU jobs (CPU -> GPU when GPUs become available)
        # Re-check availability in case burst promotions consumed some nodes
        idle_nodes=$(check_allocated_availability)
        if [ "${idle_nodes}" -gt 0 ] && [ -n "${pending_cpu}" ]; then
            log_info "Processing pending CPU jobs (CPU -> GPU)..."
            while IFS= read -r line; do
                if [ -z "${line}" ]; then
                    continue
                fi

                job_id=$(echo "${line}" | awk '{print $1}')
                job_name=$(echo "${line}" | awk '{print $2}')

                if [ -n "${job_id}" ] && [ -n "${job_name}" ]; then
                    if promote_cpu_to_gpu "${job_id}" "${job_name}"; then
                        ((cpu_promoted_count++))
                        ((promoted_count++))
                    fi

                    # Re-check availability after each promotion
                    idle_nodes=$(check_allocated_availability)
                    if [ "${idle_nodes}" -eq 0 ]; then
                        log_info "No more idle nodes - waiting for next cycle"
                        break
                    fi
                fi
            done <<< "${pending_cpu}"
        fi
    else
        log_info "No idle nodes on ${PARTITION} - waiting..."
    fi

    # Sleep before next check
    log_info "Sleeping ${CHECK_INTERVAL}s before next check..."
    sleep "${CHECK_INTERVAL}"
done

log_info "============================================================"
log_info "Monitor Summary"
log_info "============================================================"
log_info "Burst -> Allocated: ${burst_promoted_count}"
log_info "CPU -> GPU:         ${cpu_promoted_count}"
log_info "Total promoted:     ${promoted_count}"
log_info "Runtime:            ${elapsed}s"
log_info "============================================================"
