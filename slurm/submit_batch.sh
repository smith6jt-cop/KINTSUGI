#!/bin/bash
# =============================================================================
# KINTSUGI Batch Submission - Process Multiple Datasets
# =============================================================================

KINTSUGI_SLURM="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

usage() {
    cat << EOF
KINTSUGI Batch Pipeline Submission

Process multiple datasets in parallel or sequence.

Usage:
  $0 <project_list.txt>              # Submit all projects in file
  $0 --list <dir1> <dir2> ...        # Submit listed directories
  $0 --find <parent_dir>             # Find and submit all projects under parent

Options:
  --steps STEPS       Processing steps (default: all)
  --cycles RANGE      Override cycle range for all projects
  --sequential        Run projects one after another (default: parallel)
  --dry-run           Preview without submitting
  --help              Show this help

Project List File Format (one project per line):
  /path/to/project1
  /path/to/project2
  # Comments start with #

Examples:
  # Submit from list file
  $0 ~/my_projects.txt

  # Submit specific projects
  $0 --list /blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L \\
            /blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC2-1L

  # Find all projects under a directory
  $0 --find /blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN

  # Run sequentially (project2 waits for project1)
  $0 --sequential ~/my_projects.txt
EOF
}

# Defaults
STEPS="all"
CYCLES=""
DRY_RUN=false
SEQUENTIAL=false
MODE=""
PROJECTS=()

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --list)
            MODE="list"
            shift
            while [[ $# -gt 0 && ! "$1" =~ ^-- ]]; do
                PROJECTS+=("$1")
                shift
            done
            ;;
        --find)
            MODE="find"
            PARENT_DIR="$2"
            shift 2
            ;;
        --steps)
            STEPS="$2"
            shift 2
            ;;
        --cycles)
            CYCLES="$2"
            shift 2
            ;;
        --sequential)
            SEQUENTIAL=true
            shift
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
            if [ -z "${MODE}" ] && [ -f "$1" ]; then
                MODE="file"
                PROJECT_FILE="$1"
            else
                echo "Unknown option: $1"
                usage
                exit 1
            fi
            shift
            ;;
    esac
done

# Collect projects based on mode
case ${MODE} in
    file)
        if [ ! -f "${PROJECT_FILE}" ]; then
            echo "ERROR: File not found: ${PROJECT_FILE}"
            exit 1
        fi
        while IFS= read -r line; do
            # Skip comments and empty lines
            [[ "$line" =~ ^#.*$ ]] && continue
            [[ -z "$line" ]] && continue
            PROJECTS+=("$line")
        done < "${PROJECT_FILE}"
        ;;

    find)
        if [ ! -d "${PARENT_DIR}" ]; then
            echo "ERROR: Directory not found: ${PARENT_DIR}"
            exit 1
        fi
        # Find directories containing data/ subfolder (KINTSUGI projects)
        while IFS= read -r dir; do
            PROJECTS+=("$dir")
        done < <(find "${PARENT_DIR}" -maxdepth 2 -type d -name "data" -exec dirname {} \;)
        ;;

    list)
        # Projects already populated from command line
        ;;

    *)
        echo "ERROR: No projects specified"
        usage
        exit 1
        ;;
esac

if [ ${#PROJECTS[@]} -eq 0 ]; then
    echo "ERROR: No projects found"
    exit 1
fi

echo "============================================================"
echo "KINTSUGI Batch Submission"
echo "============================================================"
echo "Projects: ${#PROJECTS[@]}"
echo "Steps: ${STEPS}"
echo "Sequential: ${SEQUENTIAL}"
echo "Dry run: ${DRY_RUN}"
echo "============================================================"

# Build common arguments
SUBMIT_ARGS=""
[ -n "${STEPS}" ] && SUBMIT_ARGS="${SUBMIT_ARGS} --steps ${STEPS}"
[ -n "${CYCLES}" ] && SUBMIT_ARGS="${SUBMIT_ARGS} --cycles ${CYCLES}"
[ "${DRY_RUN}" = true ] && SUBMIT_ARGS="${SUBMIT_ARGS} --dry-run"

LAST_JOB=""

for project in "${PROJECTS[@]}"; do
    echo ""
    echo ">>> Project: $(basename "${project}")"
    echo "    Path: ${project}"

    if [ ! -d "${project}/data" ]; then
        echo "    SKIPPED: Not a valid KINTSUGI project"
        continue
    fi

    cmd="${KINTSUGI_SLURM}/submit.sh --project ${project} ${SUBMIT_ARGS}"

    if [ "${SEQUENTIAL}" = true ] && [ -n "${LAST_JOB}" ]; then
        # TODO: Add dependency on last job
        echo "    (waiting for previous project)"
    fi

    if [ "${DRY_RUN}" = true ]; then
        echo "    [DRY RUN] ${cmd}"
    else
        echo "    Submitting..."
        ${cmd}
    fi
done

echo ""
echo "============================================================"
echo "Batch submission complete!"
echo "============================================================"
echo ""
echo "Monitor all jobs: squeue -u \$USER"
