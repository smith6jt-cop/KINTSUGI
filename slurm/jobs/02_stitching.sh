#!/bin/bash
# =============================================================================
# KINTSUGI Correction + Stitching Job
# GPU-accelerated BaSiC correction and image stitching
# Processes one cycle, all channels, all z-planes
# =============================================================================

set -e

module load conda
conda activate ${CONDA_ENV:-kintsugi}

export PYTHONPATH="${KINTSUGI_DIR}:${PROJECT_DIR}/notebooks:${PYTHONPATH}"

python << 'PYTHON_SCRIPT'
import sys
import os
import gc
import time
from glob import glob
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

import numpy as np
import pandas as pd
from skimage.io import imread, imsave

# Setup paths
PROJECT_DIR = Path(os.environ['PROJECT_DIR'])
KINTSUGI_DIR = Path(os.environ['KINTSUGI_DIR'])
QC_DIR = Path(os.environ.get('QC_DIR', PROJECT_DIR / 'slurm' / 'qc' / 'stitch'))
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR))

from Kstitch.stitching import stitch_images
from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending
from kintsugi.gpu import cleanup_gpu_memory

# Get configuration from environment
CYCLE = int(os.environ.get('SLURM_ARRAY_TASK_ID', 1))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL', 4))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'tiff')

# Tile grid parameters
TILE_ROWS = int(os.environ.get('TILE_ROWS', 5))
TILE_COLS = int(os.environ.get('TILE_COLS', 5))
TILE_OVERLAP = float(os.environ.get('TILE_OVERLAP', 0.1))

# BaSiC correction parameters
BASIC_FLATFIELD_MIN = 0.1
BASIC_MAX_ITERATIONS = 500
BASIC_OPTIMIZATION_TOLERANCE = 1e-6

# Stitching parameters (read from environment or use notebook-tuned defaults)
OVERLAP_PERCENTAGE = TILE_OVERLAP * 100
INITIAL_NCC_THRESHOLD = float(os.environ.get('NCC_THRESHOLD', 0.078))
POU = float(os.environ.get('POU', 0.5))
BLEND_SIGMA = 15.0

# Paths
DATA_DIR = PROJECT_DIR / 'data'
RAW_DIR = DATA_DIR / 'raw'
STITCH_DIR = DATA_DIR / 'processed' / 'stitched'
STITCH_DIR.mkdir(parents=True, exist_ok=True)
QC_DIR.mkdir(parents=True, exist_ok=True)

# Find cycle directory (supports long-form names like cyc001_reg001_200210_170925)
def find_cycle_dir(raw_dir, cycle_num):
    """Find cycle directory supporting both short and long-form names."""
    # Try glob patterns for long-form names first (e.g., cyc001_reg001_200210_170925)
    glob_patterns = [
        f"cyc{cycle_num:03d}_*",  # cyc001_... (3-digit with suffix)
        f"cyc{cycle_num:02d}_*",  # cyc01_... (2-digit with suffix)
        f"Cyc{cycle_num:02d}_*",  # Cyc01_... (capitalized)
    ]
    for pattern in glob_patterns:
        matches = list(raw_dir.glob(pattern))
        if matches:
            # Return first match (should be unique per cycle)
            return matches[0]

    # Fallback to exact short-form names (e.g., cyc001, cyc01)
    exact_patterns = [f"cyc{cycle_num:03d}", f"cyc{cycle_num:02d}", f"Cyc{cycle_num:02d}"]
    for pattern in exact_patterns:
        test_dir = raw_dir / pattern
        if test_dir.exists():
            return test_dir

    return None

cycle_dir = find_cycle_dir(RAW_DIR, CYCLE)

if cycle_dir is None:
    print(f"ERROR: Cycle directory not found in {RAW_DIR}")
    print(f"  Looking for cycle {CYCLE} (tried cyc{CYCLE:03d}*, cyc{CYCLE:02d}*, etc.)")
    print(f"  Available directories: {[d.name for d in RAW_DIR.iterdir() if d.is_dir()][:10]}")
    sys.exit(1)

print(f"Found cycle directory: {cycle_dir.name}")

# Detect z-planes from first channel
sample_files = sorted(cycle_dir.glob("*_Z*_CH1.tif"))
if not sample_files:
    print(f"ERROR: No files found matching *_Z*_CH1.tif in {cycle_dir}")
    sys.exit(1)

# Extract z-plane numbers
import re
z_planes = set()
for f in sample_files:
    match = re.search(r'_Z(\d+)_', f.name)
    if match:
        z_planes.add(int(match.group(1)))
n_zplanes = len(z_planes)
n_tiles = TILE_ROWS * TILE_COLS

print(f"\n{'='*60}")
print(f"Correction + Stitching - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Input: {cycle_dir}")
print(f"Output: {STITCH_DIR}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Z-planes: {n_zplanes}")
print(f"Tiles: {TILE_ROWS}x{TILE_COLS} = {n_tiles}")
print(f"Overlap: {TILE_OVERLAP*100:.0f}%")
print(f"QC output: {QC_DIR}")

# Generate row/col indices for snake pattern
rows = []
cols = []
for r in range(TILE_ROWS):
    for c in range(TILE_COLS):
        rows.append(r)
        if r % 2 == 0:
            cols.append(c)
        else:
            cols.append(TILE_COLS - 1 - c)

def alphanumeric_key(s):
    """Natural sorting key."""
    return [int(c) if c.isdigit() else c for c in re.split('([0-9]+)', str(s))]

def load_tiles_parallel(file_list, max_workers=4):
    """Load tiles in parallel."""
    def load_one(f):
        return imread(str(f))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        tiles = list(executor.map(load_one, file_list))
    return np.array(tiles)

def save_qc_image(data, output_path, title=""):
    """Save stitched result as QC image (downsampled)."""
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        from skimage.transform import rescale

        max_dim = 1024
        scale = min(1.0, max_dim / max(data.shape))
        if scale < 1.0:
            data_small = rescale(data, scale, preserve_range=True, anti_aliasing=True)
        else:
            data_small = data

        fig, ax = plt.subplots(figsize=(10, 10))
        vmin, vmax = np.percentile(data_small, [1, 99])
        ax.imshow(data_small, cmap='gray', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis('off')
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches='tight')
        plt.close()
    except Exception as e:
        print(f"  QC save failed: {e}")

def process_zplane(cycle, channel, zplane, device_id=0):
    """Process a single z-plane: load tiles, correct, stitch, save."""

    # Find tile files for this z-plane
    pattern = f'*_Z{zplane:03d}_CH{channel}.tif'
    tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)

    if not tile_files:
        # Try 2-digit z format
        pattern = f'*_Z{zplane:02d}_CH{channel}.tif'
        tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)

    if not tile_files:
        return None, f"No files for Z{zplane} CH{channel}"

    if len(tile_files) != n_tiles:
        print(f"  Warning: Expected {n_tiles} tiles, found {len(tile_files)}")

    # Load tiles
    tiles = load_tiles_parallel(tile_files)

    # Normalize to float
    dtype_max = np.iinfo(tiles.dtype).max if np.issubdtype(tiles.dtype, np.integer) else 1.0
    tiles_float = tiles.astype(np.float64) / dtype_max

    # BaSiC illumination correction
    corrector = KCorrectGPU(use_gpu=True, verbose=False, device_id=device_id)
    flatfield, darkfield = corrector.fit(
        tiles_float,
        if_darkfield=True,
        max_iterations=BASIC_MAX_ITERATIONS,
        optimization_tolerance=BASIC_OPTIMIZATION_TOLERANCE
    )

    # Apply correction with flatfield clamping
    flatfield_safe = np.clip(flatfield, BASIC_FLATFIELD_MIN, None)
    corrected = (tiles_float - darkfield) / flatfield_safe
    corrected = np.clip(corrected, 0, 1)
    corrected = (corrected * dtype_max).astype(np.uint16)

    # Get or compute stitch model
    output_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    output_dir.mkdir(parents=True, exist_ok=True)
    pkl_path = output_dir / "result_df.pkl"

    ref_zplane = n_zplanes // 2
    ch1_pkl = STITCH_DIR / f"cyc{CYCLE:02d}" / "CH1" / "result_df.pkl"

    if channel == 1 and zplane == ref_zplane:
        # Compute stitch model from CH1 reference z-plane
        result_df, _ = stitch_images(
            corrected, rows, cols,
            initial_ncc_threshold=INITIAL_NCC_THRESHOLD,
            overlap_percentage=OVERLAP_PERCENTAGE,
            pou=POU,
            max_cores=4,
            use_gpu=True
        )
        result_df.to_pickle(str(pkl_path))
    elif ch1_pkl.exists():
        result_df = pd.read_pickle(str(ch1_pkl))
    else:
        return None, f"No stitch model at {ch1_pkl}"

    # Normalize positions
    result_df = result_df.copy()
    result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
    result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

    # Stitch with blending
    overlap_fraction = (TILE_OVERLAP, TILE_OVERLAP)
    stitched = stitch_with_blending(
        corrected,
        result_df,
        blend=True,
        sigma=BLEND_SIGMA,
        overlap_fraction=overlap_fraction,
        output_dtype=np.uint16
    )

    # Save stitched result
    output_file = output_dir / f"{zplane:02d}.tif"
    imsave(str(output_file), stitched, check_contrast=False)

    return str(output_file), None

# Process all channels and z-planes
start_time = time.time()
channels = list(range(START_CHANNEL, END_CHANNEL + 1))
ref_zplane = n_zplanes // 2

for channel in channels:
    print(f"\n--- Channel {channel}/{END_CHANNEL} ---")
    ch_start = time.time()

    output_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Process reference z-plane first for CH1 (stitch model)
    if channel == 1:
        print(f"  Computing stitch model from Z{ref_zplane}...")
        result, error = process_zplane(CYCLE, channel, ref_zplane)
        if error:
            print(f"  ERROR: {error}")
            sys.exit(1)

    # Process all z-planes
    processed = 0
    errors = []

    for zplane in range(1, n_zplanes + 1):
        if channel == 1 and zplane == ref_zplane:
            processed += 1
            continue  # Already done

        result, error = process_zplane(CYCLE, channel, zplane)
        if error:
            errors.append(f"Z{zplane}: {error}")
        else:
            processed += 1

        # Progress update every 5 z-planes
        if processed % 5 == 0:
            print(f"  Progress: {processed}/{n_zplanes} z-planes")

    ch_time = time.time() - ch_start
    print(f"  Channel {channel} complete: {processed}/{n_zplanes} z-planes in {ch_time:.1f}s")

    if errors:
        print(f"  Errors: {len(errors)}")
        for e in errors[:3]:
            print(f"    - {e}")

    # Generate QC image from middle z-plane
    qc_file = output_dir / f"{ref_zplane:02d}.tif"
    if qc_file.exists():
        qc_data = imread(str(qc_file))
        qc_path = QC_DIR / f"cyc{CYCLE:02d}_CH{channel}_stitched.png"
        save_qc_image(qc_data, str(qc_path), f"Cycle {CYCLE} CH{channel} Z{ref_zplane}")
        print(f"  QC saved: {qc_path.name}")

    # GPU cleanup between channels
    cleanup_gpu_memory()
    gc.collect()

total_time = time.time() - start_time

print(f"\n{'='*60}")
print(f"Correction + Stitching Complete")
print(f"{'='*60}")
print(f"Cycle: {CYCLE}")
print(f"Channels: {len(channels)}")
print(f"Z-planes per channel: {n_zplanes}")
print(f"Total time: {total_time/60:.1f} minutes")
print(f"Output: {STITCH_DIR}")
print(f"{'='*60}")
PYTHON_SCRIPT

exit $?
