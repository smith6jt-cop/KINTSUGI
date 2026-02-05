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
import json
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

# =============================================================================
# METADATA LOADING
# =============================================================================

def load_experiment_config(project_dir):
    """Load experiment configuration from /meta/experiment.json."""
    config_path = project_dir / 'meta' / 'experiment.json'
    if config_path.exists():
        try:
            with open(config_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load experiment.json: {e}")
    return {}

def load_project_config(project_dir):
    """Load project configuration from kintsugi_project.json."""
    config_path = project_dir / 'kintsugi_project.json'
    if config_path.exists():
        try:
            with open(config_path) as f:
                return json.load(f)
        except (json.JSONDecodeError, IOError) as e:
            print(f"Warning: Could not load kintsugi_project.json: {e}")
    return {}

# Load metadata files
experiment_config = load_experiment_config(PROJECT_DIR)
project_config = load_project_config(PROJECT_DIR)

# Check device mode from environment (set by submit.sh based on account type)
DEVICE_MODE = os.environ.get('KINTSUGI_DEVICE_MODE', 'gpu')

print("=" * 60)
print("Metadata Loading")
print("=" * 60)
if experiment_config:
    print(f"  Loaded experiment.json")
else:
    print(f"  experiment.json not found, using environment defaults")
print(f"  Device mode: {DEVICE_MODE}")

# Get configuration: experiment.json -> environment variable -> default
# CYCLE can be passed via --export or derived from SLURM_ARRAY_TASK_ID
CYCLE = int(os.environ.get('CYCLE', os.environ.get('SLURM_ARRAY_TASK_ID', 1)))
START_CHANNEL = int(os.environ.get('START_CHANNEL', 1))
END_CHANNEL = int(os.environ.get('END_CHANNEL',
    experiment_config.get('channels_per_cycle', 4)))
OUTPUT_FORMAT = os.environ.get('OUTPUT_FORMAT', 'tiff')

# Tile grid parameters (from experiment.json or environment)
TILE_ROWS = int(os.environ.get('TILE_ROWS',
    experiment_config.get('tile_rows', 5)))
TILE_COLS = int(os.environ.get('TILE_COLS',
    experiment_config.get('tile_cols', 5)))
TILE_OVERLAP = float(os.environ.get('TILE_OVERLAP',
    experiment_config.get('tile_overlap', 0.1)))

# BaSiC correction parameters
BASIC_FLATFIELD_MIN = 0.1
BASIC_MAX_ITERATIONS = 500
BASIC_OPTIMIZATION_TOLERANCE = 1e-6

# Stitching parameters (from experiment.json or environment)
OVERLAP_PERCENTAGE = TILE_OVERLAP * 100
INITIAL_NCC_THRESHOLD = float(os.environ.get('NCC_THRESHOLD',
    experiment_config.get('ncc_threshold', 0.078)))
POU = float(os.environ.get('POU',
    experiment_config.get('pou', 0.5)))
BLEND_SIGMA = 15.0

print(f"  Tile grid: {TILE_ROWS}x{TILE_COLS}")
print(f"  Tile overlap: {TILE_OVERLAP*100:.0f}%")
print(f"  NCC threshold: {INITIAL_NCC_THRESHOLD}")
print("=" * 60)

# Initialize GPU with fallback to CPU
if DEVICE_MODE != 'cpu':
    try:
        import cupy as cp
        # Initialize CUDA context on device 0 (SLURM sets CUDA_VISIBLE_DEVICES)
        cp.cuda.Device(0).use()
        # Verify GPU is accessible with a simple operation
        _ = cp.zeros(1)
        print(f"CUDA initialized successfully on device 0")
        print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
    except Exception as e:
        print(f"WARNING: CUDA initialization failed: {e}")
        print(f"CUDA_VISIBLE_DEVICES: {os.environ.get('CUDA_VISIBLE_DEVICES', 'not set')}")
        print("Falling back to CPU processing")
        DEVICE_MODE = 'cpu'
else:
    print("Running in CPU mode (CPU-only burst account)")

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

# =============================================================================
# SKIP-EXISTING CHECK
# =============================================================================
def check_cycle_complete(stitch_dir, cycle, start_ch, end_ch, n_zplanes):
    """Check if cycle output already exists and is complete."""
    cycle_out = stitch_dir / f"cyc{cycle:02d}"
    if not cycle_out.exists():
        return False, "output directory missing"

    for ch in range(start_ch, end_ch + 1):
        ch_dir = cycle_out / f"CH{ch}"
        if not ch_dir.exists():
            return False, f"CH{ch} directory missing"

        # Check for expected z-plane files
        tif_files = list(ch_dir.glob("*.tif"))
        if len(tif_files) < n_zplanes:
            return False, f"CH{ch} has {len(tif_files)}/{n_zplanes} z-planes"

    return True, "complete"

# Check if this cycle is already processed
is_complete, status = check_cycle_complete(STITCH_DIR, CYCLE, START_CHANNEL, END_CHANNEL, n_zplanes)
if is_complete:
    print(f"\n{'='*60}")
    print(f"SKIP: Cycle {CYCLE} already processed ({status})")
    print(f"{'='*60}")
    print(f"Output exists: {STITCH_DIR / f'cyc{CYCLE:02d}'}")
    print("To reprocess, delete the output directory first.")
    sys.exit(0)
else:
    print(f"\nCycle {CYCLE} needs processing: {status}")

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

    # BaSiC illumination correction (use GPU if available)
    use_gpu = (DEVICE_MODE == 'gpu')
    corrector = KCorrectGPU(use_gpu=use_gpu, verbose=False, device_id=device_id if use_gpu else None)
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
            use_gpu=use_gpu
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

# =============================================================================
# WRITE COMPLETION MARKER
# =============================================================================
# This marker signals to downstream jobs (deconvolution) that stitching is complete
# and all output files have been written and flushed to disk/NFS.

output_cycle_dir = STITCH_DIR / f"cyc{CYCLE:02d}"
marker_file = output_cycle_dir / ".complete"

try:
    with open(marker_file, 'w') as f:
        f.write(f"stage=stitch\n")
        f.write(f"cycle={CYCLE}\n")
        f.write(f"completed={datetime.now().isoformat()}\n")
        f.write(f"channels={START_CHANNEL}-{END_CHANNEL}\n")
        f.write(f"zplanes={n_zplanes}\n")
        f.write(f"job_id={os.environ.get('SLURM_JOB_ID', 'local')}\n")
        f.write(f"array_task_id={os.environ.get('SLURM_ARRAY_TASK_ID', '0')}\n")
        f.write(f"duration_minutes={total_time/60:.1f}\n")

    # Force NFS sync to ensure marker is visible to other nodes
    os.sync()
    print(f"Completion marker written: {marker_file}")
except Exception as e:
    print(f"WARNING: Failed to write completion marker: {e}")
    # Don't fail the job if marker write fails - the data is still valid
PYTHON_SCRIPT

exit $?
