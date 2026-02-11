#!/bin/bash
# setup_all_projects.sh - Create KINTSUGI projects from manifest
set -euo pipefail

MANIFEST="/blue/maigan/smith6jt/dataset_manifest.csv"
BLUE_BASE="/blue/maigan/smith6jt/KINTSUGI_Projects"
KINTSUGI_DIR="/blue/maigan/smith6jt/KINTSUGI"
LOG="${BLUE_BASE}/setup.log"

mkdir -p "${BLUE_BASE}"
echo "Setup started: $(date)" | tee "$LOG"

success=0
failed=0

# Read manifest (skip header)
while IFS=',' read -r dataset orange_path \
    n_cycles n_cycle_dirs channels tile_rows tile_cols n_zplanes \
    xy_pixel z_step na ri ch_file raw_gb first_cycle rest; do

    PROJECT_DIR="${BLUE_BASE}/${dataset}"

    echo "========================================" | tee -a "$LOG"
    echo "Setting up: ${dataset}" | tee -a "$LOG"
    echo "  Source: ${orange_path}" | tee -a "$LOG"
    echo "  Cycles: ${n_cycles}, Grid: ${tile_rows}x${tile_cols}" | tee -a "$LOG"
    echo "========================================" | tee -a "$LOG"

    # 1. Create project directory structure
    mkdir -p "${PROJECT_DIR}/data/raw"
    mkdir -p "${PROJECT_DIR}/data/processed"
    mkdir -p "${PROJECT_DIR}/meta"
    mkdir -p "${PROJECT_DIR}/notebooks"
    mkdir -p "${PROJECT_DIR}/workflow"
    mkdir -p "${PROJECT_DIR}/configs"
    mkdir -p "${PROJECT_DIR}/logs"

    # 2. Copy metadata files from orange
    if [ -f "${orange_path}/experiment.json" ]; then
        echo "  Copying experiment.json..." | tee -a "$LOG"
        cp "${orange_path}/experiment.json" "${PROJECT_DIR}/meta/"
    else
        echo "  ERROR: experiment.json not found at ${orange_path}" | tee -a "$LOG"
        failed=$((failed + 1))
        continue
    fi

    if [ "${ch_file}" != "MISSING" ]; then
        # Check both orange_path and parent (wrapper) dir for channel names file
        if [ -f "${orange_path}/${ch_file}" ]; then
            echo "  Copying ${ch_file} from data root..." | tee -a "$LOG"
            cp "${orange_path}/${ch_file}" "${PROJECT_DIR}/meta/CHANNELNAMES.txt"
        elif [ -f "$(dirname "${orange_path}")/${ch_file}" ]; then
            echo "  Copying ${ch_file} from wrapper dir..." | tee -a "$LOG"
            cp "$(dirname "${orange_path}")/${ch_file}" "${PROJECT_DIR}/meta/CHANNELNAMES.txt"
        else
            echo "  WARNING: ${ch_file} not found at expected locations" | tee -a "$LOG"
        fi
    else
        echo "  WARNING: No channel names file found for ${dataset}" | tee -a "$LOG"
    fi

    # 3. Initialize KINTSUGI project
    # Default ? values to known constants for this instrument
    [[ "${ri}" == "?" ]] && ri="1.44"
    [[ "${na}" == "?" ]] && na="0.75"

    if [ ! -f "${PROJECT_DIR}/kintsugi_project.json" ]; then
        echo "  Running kintsugi init..." | tee -a "$LOG"
        if kintsugi init "${PROJECT_DIR}" \
            --name "${dataset}" \
            --tile-rows "${tile_rows}" \
            --tile-cols "${tile_cols}" \
            --xy-pixel-size "${xy_pixel}" \
            --z-step-size "${z_step}" \
            --numerical-aperture "${na}" \
            --tissue-ri "${ri}" \
            --slurm \
            --force 2>&1 | tee -a "$LOG"; then
            echo "  kintsugi init: OK" | tee -a "$LOG"
            success=$((success + 1))
        else
            echo "  kintsugi init: FAILED (exit code $?)" | tee -a "$LOG"
            failed=$((failed + 1))
        fi
    else
        echo "  Project already initialized, skipping init" | tee -a "$LOG"
        success=$((success + 1))
    fi

    echo "  Done: ${dataset}" | tee -a "$LOG"
    echo "" | tee -a "$LOG"

done < <(tail -n +2 "$MANIFEST")

echo "========================================" | tee -a "$LOG"
echo "Setup complete: $(date)" | tee -a "$LOG"
echo "  Success: ${success}" | tee -a "$LOG"
echo "  Failed:  ${failed}" | tee -a "$LOG"
echo "  Log: ${LOG}" | tee -a "$LOG"
echo "All projects created. Next: stage raw data." | tee -a "$LOG"
