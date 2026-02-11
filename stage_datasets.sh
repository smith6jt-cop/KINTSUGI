#!/bin/bash
# stage_datasets.sh - Submit rsync staging jobs to SLURM
set -euo pipefail

MANIFEST="/blue/maigan/smith6jt/dataset_manifest.csv"
BLUE_BASE="/blue/maigan/smith6jt/KINTSUGI_Projects"
MAX_CONCURRENT=5   # limit active staging jobs

count=0
tail -n +2 "$MANIFEST" | while IFS=',' read -r dataset orange_path \
    n_cycles n_cycle_dirs channels tile_rows tile_cols n_zplanes \
    xy_pixel z_step na ri ch_file raw_gb first_cycle rest; do

    PROJECT_DIR="${BLUE_BASE}/${dataset}"
    RAW_DIR="${PROJECT_DIR}/data/raw"

    # Skip if already staged (check for .staged marker)
    if [ -f "${RAW_DIR}/.staged" ]; then
        echo "SKIP (already staged): ${dataset}"
        continue
    fi

    # Throttle: wait if too many staging jobs are running
    while [ $(squeue -u $USER -n 'kintsugi_stage*' -h | wc -l) \
            -ge ${MAX_CONCURRENT} ]; do
        echo "Waiting for staging slot (${MAX_CONCURRENT} active)..."
        sleep 60
    done

    echo "Staging: ${dataset} (${raw_gb} GB)"

    sbatch --job-name="kintsugi_stage_${dataset}" \
        --partition=hpg-default --account=maigan \
        --mem=4gb --cpus-per-task=2 --time=06:00:00 \
        --output="${PROJECT_DIR}/logs/stage_%j.out" \
        --wrap="
            rsync -av '${orange_path}'/*/  '${RAW_DIR}/' \
              --include='*/' --include='*.tif' --exclude='*' && \
            echo 'staged from ${orange_path}' > '${RAW_DIR}/.staged' && \
            echo 'Staging complete: ${dataset}'
        "

    count=$((count + 1))
done

echo "Submitted ${count} staging jobs."
echo "Monitor with: squeue -u \$USER -n 'kintsugi_stage*'"
