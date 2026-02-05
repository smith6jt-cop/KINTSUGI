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
USE_BURST=false

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

# Check if user has hit QOS limits on a partition
# Returns: 0 if no QOS issues, 1 if QOS limits reached
check_qos_limits() {
    local partition=$1
    local account=${2:-${ACCOUNT}}

    if ! command -v squeue &> /dev/null; then
        return 0  # Can't check, assume OK
    fi

    # Check for pending jobs with QOS-related reasons
    # QOSGrpGRES = GPU quota reached, QOSGrpCpuLimit = CPU quota reached
    local qos_blocked
    qos_blocked=$(squeue -u "${USER}" -p "${partition}" -A "${account}" -h -t PD \
        -o "%r" 2>/dev/null | grep -c "QOSGrp" || echo "0")

    if [ "${qos_blocked}" -gt 0 ]; then
        term_warn "Account ${account} has ${qos_blocked} jobs blocked by QOS limits on ${partition}"
        return 1
    fi

    # Also check running+pending GPU count against typical QOS limits
    # This catches the case where we'd exceed limits with new submission
    local running_gpus
    running_gpus=$(squeue -u "${USER}" -p "${partition}" -A "${account}" -h -t R,PD \
        -o "%b" 2>/dev/null | grep -o "gpu:[0-9]*" | cut -d: -f2 | \
        awk '{sum+=$1} END {print sum+0}' || echo "0")

    # Get QOS GPU limit if available (requires sacctmgr)
    local qos_gpu_limit=""
    if command -v sacctmgr &> /dev/null; then
        # Try to get the GrpTRES limit for GPUs from the account's QOS
        local qos_name="${account}"
        qos_gpu_limit=$(sacctmgr show qos "${qos_name}" format=GrpTRES -n -P 2>/dev/null | \
            grep -o "gres/gpu=[0-9]*" | cut -d= -f2 || echo "")
    fi

    if [ -n "${qos_gpu_limit}" ] && [ "${qos_gpu_limit}" -gt 0 ] 2>/dev/null; then
        if [ "${running_gpus}" -ge "${qos_gpu_limit}" ]; then
            term_warn "Account ${account} using ${running_gpus}/${qos_gpu_limit} GPUs (at limit)"
            return 1
        fi
        term_info "Account ${account} using ${running_gpus}/${qos_gpu_limit} GPUs"
    fi

    return 0
}

# Check if GPU resources are available in a partition
# Returns: 0 if idle GPUs available, 1 if partition invalid, 2 if resources busy, 3 if QOS limits reached
check_gpu_availability() {
    local partition=$1
    local account=${2:-${ACCOUNT}}

    if [ "${DRY_RUN}" = true ]; then
        term_info "DRY RUN: Skipping GPU availability check for partition ${partition}"
        return 0
    fi

    # Check if sinfo is available
    if ! command -v sinfo &> /dev/null; then
        term_warn "sinfo command not found - cannot verify GPU availability"
        return 0  # Proceed anyway, let sbatch handle it
    fi

    term_info "Checking GPU availability on partition ${partition}..."

    # FIRST: Check QOS limits - this takes priority over node availability
    # Even if nodes are idle, QOS limits will block jobs
    if ! check_qos_limits "${partition}" "${account}"; then
        term_warn "QOS limits reached for account ${account} on ${partition}"
        return 3
    fi

    # Check if partition exists and get idle node count
    # Format: NODES(A/I/O/T) where I = idle nodes
    local partition_info
    partition_info=$(sinfo -p "${partition}" -h -o "%F" 2>/dev/null | head -1)

    if [ -z "${partition_info}" ]; then
        term_error "Partition ${partition} not found or has no nodes"
        return 1
    fi

    # Extract idle count from A/I/O/T format (second field)
    local idle_nodes
    idle_nodes=$(echo "${partition_info}" | cut -d'/' -f2)

    if [ "${idle_nodes}" -gt 0 ] 2>/dev/null; then
        term_info "Found ${idle_nodes} idle node(s) on ${partition}"
        return 0
    fi

    # No idle nodes - check queue depth
    local queue_count=0
    if command -v squeue &> /dev/null; then
        queue_count=$(squeue -p "${partition}" -h 2>/dev/null | wc -l || echo "0")

        if [ "${queue_count}" -gt 10 ]; then
            term_warn "High queue depth (${queue_count} jobs) on ${partition}"
        fi
    fi

    term_warn "No idle nodes on ${partition} (0 idle, ${queue_count} jobs queued) - partition busy"
    return 2
}

# Get current GPU resource status
show_gpu_status() {
    local partition=$1

    if ! command -v sinfo &> /dev/null; then
        return
    fi

    term_info "Current node status on partition ${partition}:"
    sinfo -p "${partition}" -N -o "  %N %G %t %C" 2>/dev/null | head -10 >&2 || \
        echo "  No node info available" >&2
}

# =============================================================================
# TIME UTILITIES
# =============================================================================

# Apply time multiplier to a time string (HH:MM:SS)
# Usage: apply_time_multiplier "02:00:00" 5 -> "10:00:00"
apply_time_multiplier() {
    local time_str=$1
    local multiplier=${2:-1}

    # Parse HH:MM:SS
    local hours minutes seconds
    IFS=':' read -r hours minutes seconds <<< "${time_str}"

    # Convert to total seconds
    local total_seconds=$(( (10#${hours} * 3600) + (10#${minutes} * 60) + 10#${seconds} ))

    # Apply multiplier
    total_seconds=$(( total_seconds * multiplier ))

    # Convert back to HH:MM:SS
    hours=$(( total_seconds / 3600 ))
    minutes=$(( (total_seconds % 3600) / 60 ))
    seconds=$(( total_seconds % 60 ))

    printf "%02d:%02d:%02d" ${hours} ${minutes} ${seconds}
}

# =============================================================================
# ACCOUNT CHAIN SELECTION
# =============================================================================

# Check if account is CPU-only (burst account)
# Returns: 0 if CPU-only, 1 if GPU-enabled
is_cpu_only_account() {
    local account=$1
    local cpu_only="${CPU_ONLY_ACCOUNTS:-}"

    # If CPU_ONLY_ACCOUNTS is set, check against it
    if [ -n "${cpu_only}" ]; then
        IFS=',' read -ra CPU_ACCOUNTS <<< "${cpu_only}"
        for cpu_acct in "${CPU_ACCOUNTS[@]}"; do
            cpu_acct=$(echo "${cpu_acct}" | xargs)
            [ "${account}" = "${cpu_acct}" ] && return 0
        done
        return 1
    fi

    # Default heuristic: accounts ending in -b are CPU-only burst accounts
    [[ "${account}" == *-b ]] && return 0 || return 1
}

# Get account chain (backward compatible with single ACCOUNT)
get_account_chain() {
    if [ -n "${ACCOUNT_CHAIN:-}" ]; then
        echo "${ACCOUNT_CHAIN}"
    elif [ -n "${ACCOUNT:-}" ]; then
        echo "${ACCOUNT}"
    else
        echo ""
    fi
}

# Select first available account from chain
# Sets global: SELECTED_ACCOUNT, SELECTED_MODE
select_account() {
    local chain
    chain=$(get_account_chain)

    if [ -z "${chain}" ]; then
        term_error "No accounts configured. Set ACCOUNT or ACCOUNT_CHAIN in config.sh"
        return 1
    fi

    IFS=',' read -ra ACCOUNTS <<< "${chain}"

    # For now, select the first account in the chain
    # Future enhancement: check account availability with sacctmgr
    for account in "${ACCOUNTS[@]}"; do
        account=$(echo "${account}" | xargs)  # Trim whitespace

        if is_cpu_only_account "${account}"; then
            SELECTED_ACCOUNT="${account}"
            SELECTED_MODE="cpu"
            term_info "Selected account: ${account} (CPU-only burst)"
            return 0
        else
            SELECTED_ACCOUNT="${account}"
            SELECTED_MODE="gpu"
            term_info "Selected account: ${account} (GPU-enabled)"
            return 0
        fi
    done

    term_error "No available accounts in chain: ${chain}"
    return 1
}

# Get burst QOS name for an account
# Adds "-b" suffix if not already present
get_burst_qos() {
    local account=$1
    if [[ "${account}" == *-b ]]; then
        echo "${account}"
    else
        echo "${account}-b"
    fi
}

# Check if burst QOS exists for current account
# Returns 0 if burst QOS is available, 1 otherwise
burst_qos_available() {
    local account_to_check="${1:-${SELECTED_ACCOUNT:-${ACCOUNT:-}}}"

    if [ -z "${account_to_check}" ]; then
        return 1
    fi

    local burst_qos
    burst_qos=$(get_burst_qos "${account_to_check}")

    if [ "${DRY_RUN}" = true ]; then
        # In dry run, assume burst is available
        return 0
    fi

    # Check if QOS exists using sacctmgr
    if command -v sacctmgr &> /dev/null; then
        local qos_list
        qos_list=$(sacctmgr show qos format=name -n 2>/dev/null || true)
        if [ -n "${qos_list}" ] && echo "${qos_list}" | grep -q "^${burst_qos}$"; then
            return 0
        fi
    fi

    # Fallback: check if user has access to burst QOS via association
    if command -v sacctmgr &> /dev/null; then
        local assoc_qos
        assoc_qos=$(sacctmgr show assoc where user="${USER}" format=qos -n 2>/dev/null || true)
        if [ -n "${assoc_qos}" ] && echo "${assoc_qos}" | grep -q "${burst_qos}"; then
            return 0
        fi
    fi

    return 1
}

# Calculate effective resources based on GPU or CPU mode
# Sets globals: EFFECTIVE_GPUS, EFFECTIVE_CPUS_PER_TASK, EFFECTIVE_MEM_*, EFFECTIVE_TIME_*
calculate_resources_for_mode() {
    local mode=$1

    if [ "${mode}" = "cpu" ]; then
        # CPU mode: no GPUs, more CPUs, adjusted memory
        EFFECTIVE_GPUS=0
        EFFECTIVE_CPUS_PER_TASK=${CPU_CPUS_PER_TASK:-8}
        EFFECTIVE_MEM_CORRECTION=${CPU_MEM_CORRECTION:-32}
        EFFECTIVE_MEM_STITCH=${CPU_MEM_STITCH:-48}
        EFFECTIVE_MEM_DECON=${CPU_MEM_DECON:-64}
        EFFECTIVE_MEM_EDF=${CPU_MEM_EDF:-24}
        EFFECTIVE_PARTITION=${PARTITION_CPU:-hpg-default}

        # Apply time multiplier for CPU mode
        local multiplier=${CPU_TIME_MULTIPLIER:-5}
        EFFECTIVE_TIME_CORRECTION=$(multiply_time "${TIME_CORRECTION}" "${multiplier}")
        EFFECTIVE_TIME_STITCH=$(multiply_time "${TIME_STITCH}" "${multiplier}")
        EFFECTIVE_TIME_DECON=$(multiply_time "${TIME_DECON}" "${multiplier}")
        EFFECTIVE_TIME_EDF=$(multiply_time "${TIME_EDF}" "${multiplier}")

        term_info "CPU mode resources:"
        term_info "  CPUs per task: ${EFFECTIVE_CPUS_PER_TASK}"
        term_info "  Partition: ${EFFECTIVE_PARTITION}"
        term_info "  Time multiplier: ${multiplier}x"
    else
        # GPU mode: use standard settings
        EFFECTIVE_GPUS=${GPUS_PER_NODE:-1}
        EFFECTIVE_CPUS_PER_TASK=${CPUS_PER_TASK:-4}
        EFFECTIVE_MEM_CORRECTION=${MEM_CORRECTION:-32}
        EFFECTIVE_MEM_STITCH=${MEM_STITCH:-64}
        EFFECTIVE_MEM_DECON=${MEM_DECON:-96}
        EFFECTIVE_MEM_EDF=${MEM_EDF:-32}
        EFFECTIVE_PARTITION=${PARTITION}
        EFFECTIVE_TIME_CORRECTION=${TIME_CORRECTION}
        EFFECTIVE_TIME_STITCH=${TIME_STITCH}
        EFFECTIVE_TIME_DECON=${TIME_DECON}
        EFFECTIVE_TIME_EDF=${TIME_EDF}
    fi
}

# Multiply time string by a factor
# Usage: multiply_time "04:00:00" 5 -> "20:00:00"
multiply_time() {
    local time_str=$1
    local multiplier=$2

    # Parse HH:MM:SS
    local hours=$(echo "${time_str}" | cut -d: -f1 | sed 's/^0*//')
    local mins=$(echo "${time_str}" | cut -d: -f2 | sed 's/^0*//')
    local secs=$(echo "${time_str}" | cut -d: -f3 | sed 's/^0*//')

    # Handle empty values
    [ -z "${hours}" ] && hours=0
    [ -z "${mins}" ] && mins=0
    [ -z "${secs}" ] && secs=0

    # Calculate total seconds and multiply
    local total_secs=$(( (hours * 3600 + mins * 60 + secs) * multiplier ))

    # Cap at 7 days (SLURM limit for most queues)
    local max_secs=$((7 * 24 * 3600))
    [ ${total_secs} -gt ${max_secs} ] && total_secs=${max_secs}

    # Convert back to HH:MM:SS
    local new_hours=$((total_secs / 3600))
    local new_mins=$(((total_secs % 3600) / 60))
    local new_secs=$((total_secs % 60))

    printf "%02d:%02d:%02d" ${new_hours} ${new_mins} ${new_secs}
}

# =============================================================================
# RESOURCE CALCULATION
# =============================================================================

# Calculate optimal max concurrent jobs based on allocation limits
# Uses a dual-pool architecture: GPU slots + CPU slots = total concurrent
# Sets COMPUTED_MAX_CONCURRENT, GPU_SLOTS, CPU_SLOTS as output (globals)
calculate_max_concurrent() {
    local cpus_per_task=${CPUS_PER_TASK:-4}
    local gpus_per_node=${GPUS_PER_NODE:-1}
    local alloc_cpus=${ALLOC_CPUS:-64}
    local alloc_mem=${ALLOC_MEM:-600}
    local alloc_gpus=${ALLOC_GPUS:-2}

    # CPU job resource requirements
    local cpu_cpus=${CPU_CPUS_PER_TASK:-8}
    local cpu_mem=${CPU_MEM_DECON:-48}

    # Find max memory requirement for GPU jobs across all steps
    local gpu_max_mem=${MEM_DECON:-96}
    [ "${MEM_STITCH:-64}" -gt "${gpu_max_mem}" ] && gpu_max_mem=${MEM_STITCH:-64}
    [ "${MEM_CORRECTION:-32}" -gt "${gpu_max_mem}" ] && gpu_max_mem=${MEM_CORRECTION:-32}
    [ "${MEM_EDF:-32}" -gt "${gpu_max_mem}" ] && gpu_max_mem=${MEM_EDF:-32}

    # ==========================================================================
    # GPU Pool: Limited by allocated GPUs
    # ==========================================================================
    local gpu_slots=$((alloc_gpus / gpus_per_node))

    # Also check GPU pool against CPU and memory limits
    local gpu_slots_by_cpu=$((alloc_cpus / cpus_per_task))
    local gpu_slots_by_mem=$((alloc_mem / gpu_max_mem))

    # GPU pool is minimum of GPU, CPU, and memory limits
    [ "${gpu_slots_by_cpu}" -lt "${gpu_slots}" ] && gpu_slots=${gpu_slots_by_cpu}
    [ "${gpu_slots_by_mem}" -lt "${gpu_slots}" ] && gpu_slots=${gpu_slots_by_mem}
    [ "${gpu_slots}" -lt 0 ] && gpu_slots=0

    # ==========================================================================
    # CPU Pool: Limited by remaining resources after GPU allocation
    # ==========================================================================
    # CPUs remaining after GPU jobs
    local cpus_used_by_gpu=$((gpu_slots * cpus_per_task))
    local cpus_for_cpu_jobs=$((alloc_cpus - cpus_used_by_gpu))
    local cpu_slots_by_cpu=$((cpus_for_cpu_jobs / cpu_cpus))

    # Memory remaining after GPU jobs
    local mem_used_by_gpu=$((gpu_slots * gpu_max_mem))
    local mem_for_cpu_jobs=$((alloc_mem - mem_used_by_gpu))
    local cpu_slots_by_mem=$((mem_for_cpu_jobs / cpu_mem))

    # CPU pool is minimum of CPU and memory limits
    local cpu_slots=${cpu_slots_by_cpu}
    [ "${cpu_slots_by_mem}" -lt "${cpu_slots}" ] && cpu_slots=${cpu_slots_by_mem}
    [ "${cpu_slots}" -lt 0 ] && cpu_slots=0

    # ==========================================================================
    # Total concurrent = GPU slots + CPU slots
    # ==========================================================================
    COMPUTED_MAX_CONCURRENT=$((gpu_slots + cpu_slots))
    GPU_SLOTS=${gpu_slots}
    CPU_SLOTS=${cpu_slots}

    # Export for use by burst_monitor.sh
    export GPU_SLOTS
    export CPU_SLOTS

    # Store calculation details for display
    RESOURCE_CALC_DETAILS="GPU pool: ${gpu_slots} (${alloc_gpus} GPUs), CPU pool: ${cpu_slots} (${cpus_for_cpu_jobs} CPUs, ${mem_for_cpu_jobs}GB remaining)"

    # Determine limiting factor for display
    if [ "${gpu_slots}" -eq 0 ]; then
        LIMITING_RESOURCE="GPUs (no GPU slots available)"
    elif [ "${cpu_slots}" -eq 0 ]; then
        LIMITING_RESOURCE="remaining resources (no CPU slots after GPU allocation)"
    elif [ "${cpu_slots_by_mem}" -lt "${cpu_slots_by_cpu}" ]; then
        LIMITING_RESOURCE="memory for CPU pool"
    else
        LIMITING_RESOURCE="CPUs for CPU pool"
    fi

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
        term_error "Per-job requirements (GPU):"
        term_error "  CPUs:   ${cpus_per_task}"
        term_error "  Memory: ${gpu_max_mem} GB (max across steps)"
        term_error "  GPUs:   ${gpus_per_node}"
        term_error ""
        term_error "Per-job requirements (CPU):"
        term_error "  CPUs:   ${cpu_cpus}"
        term_error "  Memory: ${cpu_mem} GB"
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
        term_info "Resource pool calculation:"
        term_info "  GPU job slots: ${GPU_SLOTS} (from ${ALLOC_GPUS:-2} GPUs)"
        term_info "  CPU job slots: ${CPU_SLOTS} (from remaining resources)"
        term_info "  Total concurrent jobs: ${EFFECTIVE_MAX_CONCURRENT}"
        term_info "  ${RESOURCE_CALC_DETAILS}"
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
  --use-burst         Also submit overflow cycles to burst QOS for faster processing
                      (burst jobs are preemptible but utilize idle cluster resources)
  --dry-run           Show commands without submitting
  --help              Show this help

Examples:
  # Full pipeline for a project
  $0 --project /blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L

  # Just deconvolution for cycles 1-3
  $0 --project ~/my_project --steps decon --cycles 1-3

  # Preview without submitting
  $0 --project ~/my_project --dry-run

  # Use burst resources for faster processing (preemptible overflow)
  $0 --project ~/my_project --use-burst

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
        --use-burst|-b)
            USE_BURST=true
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
term_info "Burst:    ${USE_BURST}"
term_info "------------------------------------------------------------"
term_info "GPU Configuration:"
term_info "  Primary partition:  ${PARTITION}"
if [ -n "${PARTITION_FALLBACK}" ]; then
    term_info "  Fallback partition: ${PARTITION_FALLBACK}"
else
    term_info "  Fallback partition: (not configured)"
fi
if [ "${USE_BURST}" = true ]; then
    term_info "  Burst partition:    ${PARTITION_BURST:-hpg-default} (preemptible)"
fi
term_info "------------------------------------------------------------"

# Calculate optimal concurrency based on allocation limits
term_info "Calculating resource allocation..."
if ! resolve_max_concurrent; then
    exit 1
fi

term_info "------------------------------------------------------------"
term_info "Resource Allocation (Dual-Pool Architecture):"
term_info "  Allocation limits: ${ALLOC_CPUS:-64} CPUs, ${ALLOC_MEM:-600}GB mem, ${ALLOC_GPUS:-2} GPUs"
term_info "  GPU jobs: ${CPUS_PER_TASK:-4} CPUs, ${MEM_DECON:-96}GB mem, ${GPUS_PER_NODE:-1} GPU each"
term_info "  CPU jobs: ${CPU_CPUS_PER_TASK:-8} CPUs, ${CPU_MEM_DECON:-48}GB mem each"
term_info "  GPU slots: ${GPU_SLOTS}, CPU slots: ${CPU_SLOTS}"
term_info "  Total concurrent: ${EFFECTIVE_MAX_CONCURRENT} jobs"
term_info "  Total usage: ${TOTAL_CPUS_USED} CPUs, ${TOTAL_MEM_USED}GB mem, ${TOTAL_GPUS_USED} GPUs"
term_info "============================================================"
echo ""

# Initial resource availability check
term_info "Checking cluster resource availability..."
check_gpu_availability "${PARTITION}" "${ACCOUNT}"
initial_check=$?
if [ ${initial_check} -ne 0 ]; then
    if [ ${initial_check} -eq 3 ]; then
        term_warn "QOS limits reached on primary partition ${PARTITION}"
    else
        term_warn "Primary partition resources may not be immediately available"
    fi
    if [ -n "${PARTITION_FALLBACK}" ]; then
        check_gpu_availability "${PARTITION_FALLBACK}" "${ACCOUNT}"
        fallback_check=$?
        if [ ${fallback_check} -eq 0 ]; then
            term_info "Fallback partition ${PARTITION_FALLBACK} has resources - will use automatically"
        elif [ ${fallback_check} -eq 2 ]; then
            term_info "Fallback partition ${PARTITION_FALLBACK} busy but available - will queue there"
        else
            # fallback_check is 1 (invalid) or 3 (QOS limits)
            if [ ${fallback_check} -eq 3 ]; then
                term_error "QOS limits also reached on fallback partition ${PARTITION_FALLBACK}"
            else
                term_error "Fallback partition ${PARTITION_FALLBACK} not available"
            fi
            term_error "Neither primary nor fallback partition resources appear available"
            if [ "${DRY_RUN}" != true ]; then
                if [ -t 0 ]; then
                    echo ""
                    read -p "Continue anyway? [y/N] " -n 1 -r
                    echo ""
                    if [[ ! $REPLY =~ ^[Yy]$ ]]; then
                        term_info "Submission cancelled"
                        exit 1
                    fi
                else
                    term_error "Non-interactive mode: aborting due to unavailable resources"
                    term_error "Use --dry-run to preview or ensure resources are available"
                    exit 1
                fi
            fi
        fi
    else
        term_warn "No fallback partition configured - jobs will wait for primary resources"
    fi
else
    term_info "Primary partition resources appear available"
fi

term_info "Output structure:"
term_info "  Logs:       ${LOG_DIR}/"
term_info "  QC Images:  ${QC_DIR}/"
term_info "  Summaries:  ${SUMMARY_DIR}/"
term_info "============================================================"

# Build sbatch command array with specified partition settings
# Sets SBATCH_CMD array as output (global)
# Note: GPU type is determined by partition (hpg-b200 = B200, hpg-turin = L4)
build_sbatch_cmd_array() {
    local step_name=$1
    local script=$2
    local mem=$3
    local time=$4
    local dep_type=$5
    local dep_jobid=$6
    local use_partition=$7
    local use_qos=$8
    local use_account=$9
    local use_mode=${10:-gpu}
    local use_cpus=${11:-${CPUS_PER_TASK:-4}}
    local use_gpus=${12:-${GPUS_PER_NODE:-1}}

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
        "--account=${use_account}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${use_cpus}"
        "--mem=${mem}gb"
    )

    # Only request GPUs in GPU mode
    if [ "${use_mode}" = "gpu" ] && [ "${use_gpus}" -gt 0 ]; then
        SBATCH_CMD+=("--gpus=${use_gpus}")
    fi

    # Use appropriate array limit based on job mode
    local array_limit
    if [ "${use_mode}" = "gpu" ]; then
        array_limit=${GPU_SLOTS:-${EFFECTIVE_MAX_CONCURRENT}}
    else
        array_limit=${CPU_SLOTS:-${EFFECTIVE_MAX_CONCURRENT}}
    fi

    SBATCH_CMD+=(
        "--time=${time}"
        "--array=${ARRAY_SPEC}%${array_limit}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},RUN_ID=${RUN_ID},QC_DIR=${QC_DIR}/${step_name},KINTSUGI_DEVICE_MODE=${use_mode}"
    )

    if [ -n "${dep_jobid}" ]; then
        # In dry-run mode, accept placeholder IDs (DRY_*_JOB); in real mode, validate numeric
        if [[ "${dep_jobid}" =~ ^DRY_ ]]; then
            # Dry-run placeholder - add dependency for display purposes
            SBATCH_CMD+=("--dependency=${dep_type}:${dep_jobid}")
        elif [[ "${dep_jobid}" =~ ^[0-9]+$ ]]; then
            # Real numeric job ID
            SBATCH_CMD+=("--dependency=${dep_type}:${dep_jobid}")
        else
            echo "ERROR: Invalid dependency job ID '${dep_jobid}' - must be numeric" >&2
            return 1
        fi
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

# Submit job function with account chain and fallback support
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

    # Select account from chain (sets SELECTED_ACCOUNT and SELECTED_MODE)
    if ! select_account; then
        term_error "Failed to select account - aborting"
        return 1
    fi

    # Calculate effective resources based on mode
    calculate_resources_for_mode "${SELECTED_MODE}"

    # Get step-specific memory and time based on mode
    local effective_mem effective_time
    case ${step_name} in
        correct)
            effective_mem=${EFFECTIVE_MEM_CORRECTION}
            effective_time=${EFFECTIVE_TIME_CORRECTION}
            ;;
        stitch)
            effective_mem=${EFFECTIVE_MEM_STITCH}
            effective_time=${EFFECTIVE_TIME_STITCH}
            ;;
        decon)
            effective_mem=${EFFECTIVE_MEM_DECON}
            effective_time=${EFFECTIVE_TIME_DECON}
            ;;
        edf)
            effective_mem=${EFFECTIVE_MEM_EDF}
            effective_time=${EFFECTIVE_TIME_EDF}
            ;;
        *)
            effective_mem=${mem}
            effective_time=${time}
            ;;
    esac

    term_info "Preparing job submission for ${step_name}..."
    term_info "  Account: ${SELECTED_ACCOUNT} (${SELECTED_MODE} mode)"
    term_info "  Partition: ${EFFECTIVE_PARTITION}, Memory: ${effective_mem}GB, Time: ${effective_time}"
    if [ "${SELECTED_MODE}" = "gpu" ]; then
        term_info "  GPUs: ${EFFECTIVE_GPUS}, CPUs: ${EFFECTIVE_CPUS_PER_TASK}"
    else
        term_info "  CPUs: ${EFFECTIVE_CPUS_PER_TASK} (no GPUs - CPU-only account)"
    fi

    # Determine which partition to use
    local use_partition="${EFFECTIVE_PARTITION}"
    local use_qos="${QOS}"
    local using_fallback=false

    # For GPU mode, check availability and potentially fallback
    # Strategy: Submit to whichever partition has idle GPUs, or both if both busy
    local submit_to_fallback_too=false

    if [ "${SELECTED_MODE}" = "gpu" ]; then
        check_gpu_availability "${use_partition}" "${SELECTED_ACCOUNT}"
        local gpu_check_status=$?

        if [ ${gpu_check_status} -eq 1 ]; then
            # Partition doesn't exist or has no nodes - try fallback
            term_warn "Primary partition ${use_partition} not available"

            # Check if fallback is configured
            if [ -n "${PARTITION_FALLBACK}" ]; then
                term_info "Checking fallback partition: ${PARTITION_FALLBACK}..."
                check_gpu_availability "${PARTITION_FALLBACK}" "${SELECTED_ACCOUNT}"
                local fallback_status=$?
                # Accept fallback if GPUs are available (0) or busy (2)
                if [ ${fallback_status} -eq 0 ] || [ ${fallback_status} -eq 2 ]; then
                    term_info "Switching to fallback partition"
                    use_partition="${PARTITION_FALLBACK}"
                    use_qos="${QOS_FALLBACK}"
                    using_fallback=true
                else
                    term_error "Fallback partition ${PARTITION_FALLBACK} also not available"
                    term_error "No GPU resources available. Aborting submission."
                    show_gpu_status "${EFFECTIVE_PARTITION}"
                    return 1
                fi
            else
                term_error "No fallback partition configured and primary partition unavailable. Aborting."
                term_error "Configure PARTITION_FALLBACK in config.sh to enable fallback."
                return 1
            fi
        elif [ ${gpu_check_status} -eq 3 ]; then
            # QOS limits reached on primary partition - switch to fallback
            term_warn "QOS limits reached on ${use_partition} - switching to fallback"

            if [ -n "${PARTITION_FALLBACK}" ]; then
                term_info "Checking fallback partition: ${PARTITION_FALLBACK}..."
                check_gpu_availability "${PARTITION_FALLBACK}" "${SELECTED_ACCOUNT}"
                local fallback_status=$?

                if [ ${fallback_status} -eq 0 ] || [ ${fallback_status} -eq 2 ]; then
                    # Fallback available - switch to it (don't submit to primary at all)
                    term_info "Switching to fallback partition ${PARTITION_FALLBACK}"
                    use_partition="${PARTITION_FALLBACK}"
                    use_qos="${QOS_FALLBACK}"
                    using_fallback=true
                elif [ ${fallback_status} -eq 3 ]; then
                    # QOS limits also reached on fallback - submit to fallback anyway (will queue)
                    term_warn "QOS limits also reached on fallback - jobs will queue on ${PARTITION_FALLBACK}"
                    use_partition="${PARTITION_FALLBACK}"
                    use_qos="${QOS_FALLBACK}"
                    using_fallback=true
                else
                    term_error "Fallback partition ${PARTITION_FALLBACK} not available"
                    term_error "QOS limits reached on primary and fallback unavailable. Aborting."
                    return 1
                fi
            else
                term_error "QOS limits reached on ${use_partition} and no fallback configured."
                term_error "Configure PARTITION_FALLBACK in config.sh to enable fallback."
                return 1
            fi
        elif [ ${gpu_check_status} -eq 2 ]; then
            # Primary GPUs are busy - always queue on primary, optionally also submit to fallback
            term_info "Primary partition ${use_partition} is busy - jobs will queue"

            if [ -n "${PARTITION_FALLBACK}" ]; then
                term_info "Checking fallback partition: ${PARTITION_FALLBACK}..."
                check_gpu_availability "${PARTITION_FALLBACK}" "${SELECTED_ACCOUNT}"
                local fallback_status=$?

                if [ ${fallback_status} -eq 0 ]; then
                    # Fallback has idle GPUs - submit to BOTH (primary queues, fallback may start faster)
                    term_info "Fallback partition has idle GPUs - will ALSO submit to ${PARTITION_FALLBACK}"
                    submit_to_fallback_too=true
                elif [ ${fallback_status} -eq 2 ]; then
                    # Both partitions busy - submit to BOTH, first to get resources wins
                    term_info "Both partitions busy - will submit to BOTH (first to start wins)"
                    submit_to_fallback_too=true
                elif [ ${fallback_status} -eq 3 ]; then
                    # QOS limits on fallback but primary is just busy - stick with primary
                    term_info "Fallback has QOS limits - queuing on primary ${use_partition} only"
                else
                    # Fallback not available, stick with primary only
                    term_info "Fallback partition not available - queuing on primary ${use_partition} only"
                fi
            else
                term_info "No fallback configured - queuing on primary ${use_partition}"
            fi
        fi
    fi

    # Build command with selected partition configuration and mode
    build_sbatch_cmd_array "${step_name}" "${script}" "${effective_mem}" "${effective_time}" \
        "${dep_type}" "${dep_jobid}" "${use_partition}" "${use_qos}" "${SELECTED_ACCOUNT}" \
        "${SELECTED_MODE}" "${EFFECTIVE_CPUS_PER_TASK}" "${EFFECTIVE_GPUS}"

    if [ "${DRY_RUN}" = true ]; then
        # Print to stderr so it's visible (stdout is captured for job ID)
        if [ "${using_fallback}" = true ]; then
            term_info "DRY RUN - would execute (using fallback partition):"
        else
            term_info "DRY RUN - would execute:"
        fi
        echo "  $(format_cmd_for_display)" >&2

        # Show fallback submission if both partitions busy
        if [ "${submit_to_fallback_too}" = true ]; then
            term_info "DRY RUN - would ALSO submit to fallback partition (both busy):"
            build_sbatch_cmd_array "${step_name}" "${script}" "${effective_mem}" "${effective_time}" \
                "${dep_type}" "${dep_jobid}" "${PARTITION_FALLBACK}" "${QOS_FALLBACK}" \
                "${SELECTED_ACCOUNT}" "${SELECTED_MODE}" "${EFFECTIVE_CPUS_PER_TASK}" "${EFFECTIVE_GPUS}"
            echo "  $(format_cmd_for_display)" >&2
            term_info "NOTE: Both jobs would race - first to get resources wins"
        fi

        echo "DRY_${step_name}_JOB"
        return 0
    fi

    # Submit job
    if [ "${using_fallback}" = true ]; then
        term_info "Submitting ${step_name} to ${use_partition} (fallback)..."
    else
        term_info "Submitting ${step_name} to ${use_partition}..."
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

        # If both partitions are busy, also submit to fallback (first to start wins)
        if [ "${submit_to_fallback_too}" = true ]; then
            term_info "Also submitting to fallback partition ${PARTITION_FALLBACK}..."
            build_sbatch_cmd_array "${step_name}" "${script}" "${effective_mem}" "${effective_time}" \
                "${dep_type}" "${dep_jobid}" "${PARTITION_FALLBACK}" "${QOS_FALLBACK}" \
                "${SELECTED_ACCOUNT}" "${SELECTED_MODE}" "${EFFECTIVE_CPUS_PER_TASK}" "${EFFECTIVE_GPUS}"

            local fallback_result
            fallback_result=$("${SBATCH_CMD[@]}" 2>&1)
            if [ $? -eq 0 ]; then
                local fallback_jobid
                fallback_jobid=$(echo "${fallback_result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
                if [ -z "${fallback_jobid}" ]; then
                    fallback_jobid=$(echo "${fallback_result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
                fi
                term_info "Fallback job submitted: ${fallback_jobid} (on ${PARTITION_FALLBACK})"
                term_info "NOTE: Jobs ${jobid} and ${fallback_jobid} will race - first to get resources wins"
                term_info "      Cancel the other with: scancel <jobid> once one starts"
            else
                term_warn "Fallback submission failed (non-critical): ${fallback_result}"
            fi
        fi

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
        if [ "${using_fallback}" = false ] && [ -n "${PARTITION_FALLBACK}" ]; then
            term_info "Attempting fallback submission to ${PARTITION_FALLBACK}..."

            build_sbatch_cmd_array "${step_name}" "${script}" "${effective_mem}" "${effective_time}" \
                "${dep_type}" "${dep_jobid}" "${PARTITION_FALLBACK}" "${QOS_FALLBACK}" \
                "${SELECTED_ACCOUNT}" "${SELECTED_MODE}" "${EFFECTIVE_CPUS_PER_TASK}" "${EFFECTIVE_GPUS}"

            local fallback_result
            fallback_result=$("${SBATCH_CMD[@]}" 2>&1)
            exit_code=$?

            if [ ${exit_code} -eq 0 ]; then
                local jobid
                # Extract job ID using same robust pattern as primary submission
                jobid=$(echo "${fallback_result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
                if [ -z "${jobid}" ]; then
                    jobid=$(echo "${fallback_result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
                fi

                if [ -z "${jobid}" ] || ! [[ "${jobid}" =~ ^[0-9]+$ ]]; then
                    term_error "Failed to extract job ID from fallback sbatch output:"
                    term_error "  Output: ${fallback_result}"
                    return 1
                fi

                term_info "Fallback submission successful - Job ID: ${jobid}"
                term_info "NOTE: Using fallback partition ${PARTITION_FALLBACK} instead of ${PARTITION}"
                echo "${jobid}"
                return 0
            else
                term_error "Fallback submission also failed: ${fallback_result}"
            fi
        elif [ "${using_fallback}" = true ]; then
            term_error "Already using fallback partition, no further options available"
        else
            term_error "No fallback partition configured. Set PARTITION_FALLBACK in config.sh"
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

# Submit burst job (preemptible, opportunistic resources)
# Called after the main allocated job when --use-burst is enabled
submit_burst_job() {
    local step_name=$1
    local script=$2
    local mem=$3
    local time=$4
    local dep_type=$5
    local dep_jobid=$6

    # Ensure account is selected (may have been cleared after submit_job)
    if [ -z "${SELECTED_ACCOUNT}" ]; then
        if ! select_account; then
            term_warn "Failed to select account for burst job"
            return 0
        fi
    fi

    # Get burst QOS from the primary account (not a -b account)
    local base_account="${SELECTED_ACCOUNT}"
    # If current account ends in -b, strip it to get base
    if [[ "${base_account}" == *-b ]]; then
        base_account="${base_account%-b}"
    fi

    local burst_qos
    burst_qos=$(get_burst_qos "${base_account}")

    # Check if burst QOS is available
    if ! burst_qos_available; then
        term_warn "Burst QOS ${burst_qos} not available - skipping burst submission"
        return 0
    fi

    local job_name="kintsugi_${step_name}_burst_${PROJECT_NAME}"

    # Burst jobs use GPU mode (can request GPUs on burst, they're just preemptible)
    local burst_mode="gpu"
    local burst_gpus="${GPUS_PER_NODE:-1}"
    local burst_cpus="${CPUS_PER_TASK:-4}"

    # For GPU burst jobs, use the primary GPU partition with burst QOS
    # This allows preemptible access to GPU resources on the same partition
    local burst_partition="${PARTITION}"

    term_info "Submitting BURST job for ${step_name}..."
    term_info "  Account: ${base_account}, QOS: ${burst_qos}"
    term_info "  Partition: ${burst_partition} (preemptible GPU)"
    term_info "  Note: Burst jobs may be preempted and automatically requeued"

    # Build burst sbatch command (similar to regular but with burst QOS and --requeue)
    SBATCH_CMD=(
        sbatch
        "--job-name=${job_name}"
        "--output=${LOG_DIR}/${step_name}_burst_%A_%a.out"
        "--error=${LOG_DIR}/${step_name}_burst_%A_%a.err"
        "--partition=${burst_partition}"
        "--qos=${burst_qos}"
        "--account=${base_account}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${burst_cpus}"
        "--mem=${mem}gb"
        "--requeue"  # Automatically requeue if preempted
    )

    # Request GPUs for burst jobs (preemptible but available when cluster is idle)
    if [ "${burst_gpus}" -gt 0 ]; then
        SBATCH_CMD+=("--gpus=${burst_gpus}")
    fi

    # Burst jobs are GPU jobs, use GPU_SLOTS limit
    local burst_array_limit=${GPU_SLOTS:-${EFFECTIVE_MAX_CONCURRENT}}

    SBATCH_CMD+=(
        "--time=${time}"
        "--array=${ARRAY_SPEC}%${burst_array_limit}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},RUN_ID=${RUN_ID},QC_DIR=${QC_DIR}/${step_name},KINTSUGI_DEVICE_MODE=${burst_mode}"
    )

    if [ -n "${dep_jobid}" ]; then
        # Accept DRY_ placeholders for dry-run mode, numeric IDs for real runs
        if [[ "${dep_jobid}" =~ ^DRY_ ]] || [[ "${dep_jobid}" =~ ^[0-9]+$ ]]; then
            SBATCH_CMD+=("--dependency=${dep_type}:${dep_jobid}")
        fi
    fi

    SBATCH_CMD+=("${script}")

    if [ "${DRY_RUN}" = true ]; then
        term_info "DRY RUN - would execute (burst):"
        echo "  $(printf '%q ' "${SBATCH_CMD[@]}")" >&2
        echo "DRY_${step_name}_BURST_JOB"
        return 0
    fi

    # Submit burst job
    local result
    result=$("${SBATCH_CMD[@]}" 2>&1)
    local exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        local jobid
        jobid=$(echo "${result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
        if [ -z "${jobid}" ]; then
            jobid=$(echo "${result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
        fi

        if [ -n "${jobid}" ] && [[ "${jobid}" =~ ^[0-9]+$ ]]; then
            term_info "Burst job submitted: ${jobid} (preemptible)"
            echo "${jobid}"
            return 0
        fi
    fi

    term_warn "Burst job submission failed (non-critical): ${result}"
    return 0  # Don't fail the pipeline if burst submission fails
}

# Submit concurrent CPU job (to supplement GPU processing)
# Uses burst QOS on CPU partition for parallel processing with GPU jobs
submit_cpu_job() {
    local step_name=$1
    local script=$2
    local mem=$3
    local time=$4
    local dep_type=$5
    local dep_jobid=$6

    # Get base account (strip -b suffix if present)
    local base_account="${ACCOUNT:-maigan}"
    if [[ "${base_account}" == *-b ]]; then
        base_account="${base_account%-b}"
    fi

    # Get burst QOS for CPU jobs
    local cpu_qos
    cpu_qos=$(get_burst_qos "${base_account}")

    # Check if burst QOS is available
    if ! burst_qos_available "${base_account}"; then
        # No burst QOS available - skip CPU job
        return 0
    fi

    # Use base account with burst QOS (not a separate account)
    local cpu_account="${base_account}"

    local job_name="kintsugi_${step_name}_cpu_${PROJECT_NAME}"

    # Get step-specific CPU memory
    local cpu_mem
    case ${step_name} in
        correct) cpu_mem="${CPU_MEM_CORRECTION:-${mem}}" ;;
        stitch)  cpu_mem="${CPU_MEM_STITCH:-${mem}}" ;;
        decon)   cpu_mem="${CPU_MEM_DECON:-${mem}}" ;;
        edf)     cpu_mem="${CPU_MEM_EDF:-${mem}}" ;;
        *)       cpu_mem="${mem}" ;;
    esac

    # Apply time multiplier for CPU processing
    local cpu_time
    cpu_time=$(apply_time_multiplier "${time}" "${CPU_TIME_MULTIPLIER:-5}")

    local cpu_cpus="${CPU_CPUS_PER_TASK:-8}"
    local cpu_partition="${PARTITION_CPU:-hpg-default}"

    term_info "Submitting CPU job for ${step_name} (concurrent processing)..."
    term_info "  Account: ${cpu_account}, QOS: ${cpu_qos}, Partition: ${cpu_partition}"
    term_info "  CPUs: ${cpu_cpus}, Memory: ${cpu_mem}GB, Time: ${cpu_time}"
    term_info "  Note: CPU jobs run concurrently with GPU jobs (preemptible)"

    # Build CPU sbatch command with burst QOS
    SBATCH_CMD=(
        sbatch
        "--job-name=${job_name}"
        "--output=${LOG_DIR}/${step_name}_cpu_%A_%a.out"
        "--error=${LOG_DIR}/${step_name}_cpu_%A_%a.err"
        "--partition=${cpu_partition}"
        "--qos=${cpu_qos}"
        "--account=${cpu_account}"
        "--nodes=1"
        "--ntasks=1"
        "--cpus-per-task=${cpu_cpus}"
        "--mem=${cpu_mem}gb"
        "--time=${cpu_time}"
        "--requeue"
        "--array=${ARRAY_SPEC}%${CPU_SLOTS:-${EFFECTIVE_MAX_CONCURRENT}}"
        "--export=ALL,PROJECT_DIR=${PROJECT_DIR},KINTSUGI_DIR=${KINTSUGI_DIR},RUN_ID=${RUN_ID},QC_DIR=${QC_DIR}/${step_name},KINTSUGI_DEVICE_MODE=cpu"
    )

    # No GPUs for CPU jobs

    if [ -n "${dep_jobid}" ]; then
        # Accept DRY_ placeholders for dry-run mode, numeric IDs for real runs
        if [[ "${dep_jobid}" =~ ^DRY_ ]] || [[ "${dep_jobid}" =~ ^[0-9]+$ ]]; then
            SBATCH_CMD+=("--dependency=${dep_type}:${dep_jobid}")
        fi
    fi

    SBATCH_CMD+=("${script}")

    if [ "${DRY_RUN}" = true ]; then
        term_info "DRY RUN - would execute (CPU):"
        echo "  $(printf '%q ' "${SBATCH_CMD[@]}")" >&2
        echo "DRY_${step_name}_CPU_JOB"
        return 0
    fi

    # Submit CPU job
    local result
    result=$("${SBATCH_CMD[@]}" 2>&1)
    local exit_code=$?

    if [ ${exit_code} -eq 0 ]; then
        local jobid
        jobid=$(echo "${result}" | grep -oP 'Submitted batch job \K[0-9]+' | head -1)
        if [ -z "${jobid}" ]; then
            jobid=$(echo "${result}" | grep "Submitted batch job" | head -1 | sed 's/.*Submitted batch job //' | tr -d '[:space:]')
        fi

        if [ -n "${jobid}" ] && [[ "${jobid}" =~ ^[0-9]+$ ]]; then
            term_info "CPU job submitted: ${jobid}"
            echo "${jobid}"
            return 0
        fi
    fi

    term_warn "CPU job submission failed (non-critical): ${result}"
    return 0  # Don't fail the pipeline if CPU submission fails
}

# Job IDs for dependencies
JOB_CORRECTION=""
JOB_STITCH=""
JOB_DECON=""
JOB_EDF=""
JOB_CORRECTION_BURST=""
JOB_STITCH_BURST=""
JOB_DECON_BURST=""
JOB_EDF_BURST=""
JOB_CORRECTION_CPU=""
JOB_STITCH_CPU=""
JOB_DECON_CPU=""
JOB_EDF_CPU=""

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
            # Submit concurrent CPU job (always, if CPU accounts configured)
            JOB_CORRECTION_CPU=$(submit_cpu_job "correct" "${KINTSUGI_SLURM}/jobs/01_correction.sh" \
                "${MEM_CORRECTION}" "${TIME_CORRECTION}" "" "")
            # Submit burst job if enabled
            if [ "${USE_BURST}" = true ]; then
                JOB_CORRECTION_BURST=$(submit_burst_job "correct" "${KINTSUGI_SLURM}/jobs/01_correction.sh" \
                    "${MEM_CORRECTION}" "${TIME_CORRECTION}" "" "")
            fi
            ;;
        stitch)
            echo "--- Step 2: Stitching ---"
            JOB_STITCH=$(submit_job "stitch" "${KINTSUGI_SLURM}/jobs/02_stitching.sh" \
                "${MEM_STITCH}" "${TIME_STITCH}" "aftercorr" "${JOB_CORRECTION}")
            if [ -z "${JOB_STITCH}" ]; then
                term_error "Failed to submit stitching job - aborting pipeline"
                exit 1
            fi
            # Submit concurrent CPU job (depends on CPU correction if it exists)
            local cpu_dep="${JOB_CORRECTION_CPU:-${JOB_CORRECTION}}"
            JOB_STITCH_CPU=$(submit_cpu_job "stitch" "${KINTSUGI_SLURM}/jobs/02_stitching.sh" \
                "${MEM_STITCH}" "${TIME_STITCH}" "aftercorr" "${cpu_dep}")
            # Submit burst job if enabled (depends on burst correction if it exists)
            if [ "${USE_BURST}" = true ]; then
                local burst_dep="${JOB_CORRECTION_BURST:-${JOB_CORRECTION}}"
                JOB_STITCH_BURST=$(submit_burst_job "stitch" "${KINTSUGI_SLURM}/jobs/02_stitching.sh" \
                    "${MEM_STITCH}" "${TIME_STITCH}" "aftercorr" "${burst_dep}")
            fi
            ;;
        decon)
            echo "--- Step 3: Deconvolution ---"
            JOB_DECON=$(submit_job "decon" "${KINTSUGI_SLURM}/jobs/03_deconvolution.sh" \
                "${MEM_DECON}" "${TIME_DECON}" "aftercorr" "${JOB_STITCH}")
            if [ -z "${JOB_DECON}" ]; then
                term_error "Failed to submit deconvolution job - aborting pipeline"
                exit 1
            fi
            # Submit concurrent CPU job
            local cpu_dep="${JOB_STITCH_CPU:-${JOB_STITCH}}"
            JOB_DECON_CPU=$(submit_cpu_job "decon" "${KINTSUGI_SLURM}/jobs/03_deconvolution.sh" \
                "${MEM_DECON}" "${TIME_DECON}" "aftercorr" "${cpu_dep}")
            # Submit burst job if enabled
            if [ "${USE_BURST}" = true ]; then
                local burst_dep="${JOB_STITCH_BURST:-${JOB_STITCH}}"
                JOB_DECON_BURST=$(submit_burst_job "decon" "${KINTSUGI_SLURM}/jobs/03_deconvolution.sh" \
                    "${MEM_DECON}" "${TIME_DECON}" "aftercorr" "${burst_dep}")
            fi
            ;;
        edf)
            echo "--- Step 4: Extended Depth of Focus ---"
            JOB_EDF=$(submit_job "edf" "${KINTSUGI_SLURM}/jobs/04_edf.sh" \
                "${MEM_EDF}" "${TIME_EDF}" "aftercorr" "${JOB_DECON}")
            if [ -z "${JOB_EDF}" ]; then
                term_error "Failed to submit EDF job - aborting pipeline"
                exit 1
            fi
            # Submit concurrent CPU job
            local cpu_dep="${JOB_DECON_CPU:-${JOB_DECON}}"
            JOB_EDF_CPU=$(submit_cpu_job "edf" "${KINTSUGI_SLURM}/jobs/04_edf.sh" \
                "${MEM_EDF}" "${TIME_EDF}" "aftercorr" "${cpu_dep}")
            # Submit burst job if enabled
            if [ "${USE_BURST}" = true ]; then
                local burst_dep="${JOB_DECON_BURST:-${JOB_DECON}}"
                JOB_EDF_BURST=$(submit_burst_job "edf" "${KINTSUGI_SLURM}/jobs/04_edf.sh" \
                    "${MEM_EDF}" "${TIME_EDF}" "aftercorr" "${burst_dep}")
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
