"""
Snakemake wrapper for KINTSUGI stitching (replaces 02_stitching.sh).

Performs BaSiC illumination correction and GPU-accelerated tile stitching
for one cycle (all channels, all z-planes).

Snakemake guarantees:
  - Input raw directory exists before this script runs
  - Output sentinel is only written if this script exits successfully
  - No need for skip-existing checks (Snakemake handles re-runs)
  - No need for wait_for_input() polling (Snakemake tracks dependencies)
"""

import gc
import json
import re
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd
from skimage.io import imread, imsave

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = list(snakemake.params.channels)

# Setup Python path (same as old job scripts)
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities (replicates slurm/utils.sh structured logging)
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_header, log_footer, summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration from experiment.json and workflow config
# ---------------------------------------------------------------------------
experiment_config = {}
config_path = PROJECT_DIR / "meta" / "experiment.json"
if config_path.exists():
    with open(config_path) as f:
        experiment_config = json.load(f)

wf_config = snakemake.config

TILE_ROWS = int(wf_config.get("tile_rows", experiment_config.get("tile_rows", 5)))
TILE_COLS = int(wf_config.get("tile_cols", experiment_config.get("tile_cols", 5)))
TILE_OVERLAP = float(
    wf_config.get("tile_overlap", experiment_config.get("tile_overlap", 0.1))
)
START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)

# BaSiC parameters
BASIC_FLATFIELD_MIN = float(wf_config.get("basic_flatfield_min", 0.1))
BASIC_MAX_ITERATIONS = int(wf_config.get("basic_max_iterations", 500))
BASIC_OPTIMIZATION_TOLERANCE = float(
    wf_config.get("basic_optimization_tolerance", 1e-6)
)

# Stitching parameters
OVERLAP_PERCENTAGE = TILE_OVERLAP * 100
INITIAL_NCC_THRESHOLD = float(
    wf_config.get("ncc_threshold", experiment_config.get("ncc_threshold", 0.078))
)
POU = float(wf_config.get("pou", experiment_config.get("pou", 0.5)))
BLEND_SIGMA = float(wf_config.get("blend_sigma", 15.0))

# ---------------------------------------------------------------------------
# GPU initialization
# ---------------------------------------------------------------------------
use_gpu = True
try:
    import cupy as cp

    cp.cuda.Device(0).use()
    _ = cp.zeros(1)
    print(f"CUDA initialized successfully on device 0")
except Exception as e:
    print(f"WARNING: CUDA initialization failed: {e}")
    print("Falling back to CPU processing")
    use_gpu = False

# ---------------------------------------------------------------------------
# Imports (after PYTHONPATH setup)
# ---------------------------------------------------------------------------
from Kstitch.stitching import stitch_images
from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending
from kintsugi.gpu import cleanup_gpu_memory

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data"
RAW_DIR = DATA_DIR / "raw"
STITCH_DIR = DATA_DIR / "processed" / "stitched"
STITCH_DIR.mkdir(parents=True, exist_ok=True)


# ---------------------------------------------------------------------------
# Find cycle directory (supports long-form names)
# ---------------------------------------------------------------------------
def find_cycle_dir(raw_dir, cycle_num):
    """Find cycle directory supporting both short and long-form names."""
    glob_patterns = [
        f"cyc{cycle_num:03d}_*",
        f"cyc{cycle_num:02d}_*",
        f"Cyc{cycle_num:02d}_*",
    ]
    for pattern in glob_patterns:
        matches = list(raw_dir.glob(pattern))
        if matches:
            return matches[0]

    exact_patterns = [
        f"cyc{cycle_num:03d}",
        f"cyc{cycle_num:02d}",
        f"Cyc{cycle_num:02d}",
    ]
    for pattern in exact_patterns:
        test_dir = raw_dir / pattern
        if test_dir.exists():
            return test_dir
    return None


cycle_dir = find_cycle_dir(RAW_DIR, CYCLE)
if cycle_dir is None:
    print(f"ERROR: Cycle directory not found in {RAW_DIR}")
    print(
        f"  Looking for cycle {CYCLE} (tried cyc{CYCLE:03d}*, cyc{CYCLE:02d}*, etc.)"
    )
    print(
        f"  Available: {[d.name for d in RAW_DIR.iterdir() if d.is_dir()][:10]}"
    )
    sys.exit(1)

print(f"Found cycle directory: {cycle_dir.name}")

# ---------------------------------------------------------------------------
# Detect z-planes
# ---------------------------------------------------------------------------
sample_files = sorted(cycle_dir.glob("*_Z*_CH1.tif"))
if not sample_files:
    print(f"ERROR: No files found matching *_Z*_CH1.tif in {cycle_dir}")
    sys.exit(1)

z_planes = set()
for f in sample_files:
    match = re.search(r"_Z(\d+)_", f.name)
    if match:
        z_planes.add(int(match.group(1)))
n_zplanes = len(z_planes)
n_tiles = TILE_ROWS * TILE_COLS

# Structured logging header (mirrors utils.sh log_header)
log_header("stitch", CYCLE, PROJECT_DIR)
summary_before("stitch", PROJECT_DIR)

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
print(f"GPU: {use_gpu}")

# ---------------------------------------------------------------------------
# Tile loading helpers
# ---------------------------------------------------------------------------
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
    return [int(c) if c.isdigit() else c for c in re.split("([0-9]+)", str(s))]


def load_tiles_parallel(file_list, max_workers=4):
    """Load tiles in parallel."""

    def load_one(f):
        return imread(str(f))

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        tiles = list(executor.map(load_one, file_list))
    return np.array(tiles)


# ---------------------------------------------------------------------------
# QC helper
# ---------------------------------------------------------------------------
def save_qc_image(data, output_path, title=""):
    """Save stitched result as QC image (downsampled)."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from skimage.transform import rescale

        max_dim = 1024
        scale = min(1.0, max_dim / max(data.shape))
        if scale < 1.0:
            data_small = rescale(
                data, scale, preserve_range=True, anti_aliasing=True
            )
        else:
            data_small = data

        fig, ax = plt.subplots(figsize=(10, 10))
        vmin, vmax = np.percentile(data_small, [1, 99])
        ax.imshow(data_small, cmap="gray", vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis("off")
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches="tight")
        plt.close()
        print(f"  QC saved: {Path(output_path).name}")
    except Exception as e:
        print(f"  QC save failed: {e}")


# ---------------------------------------------------------------------------
# Per-z-plane processing
# ---------------------------------------------------------------------------
def process_zplane(cycle, channel, zplane, device_id=0):
    """Process a single z-plane: load tiles, correct, stitch, save."""
    pattern = f"*_Z{zplane:03d}_CH{channel}.tif"
    tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)

    if not tile_files:
        pattern = f"*_Z{zplane:02d}_CH{channel}.tif"
        tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)

    if not tile_files:
        return None, f"No files for Z{zplane} CH{channel}"

    if len(tile_files) != n_tiles:
        print(f"  Warning: Expected {n_tiles} tiles, found {len(tile_files)}")

    # Load tiles
    tiles = load_tiles_parallel(tile_files)

    # Normalize to float
    dtype_max = (
        np.iinfo(tiles.dtype).max
        if np.issubdtype(tiles.dtype, np.integer)
        else 1.0
    )
    tiles_float = tiles.astype(np.float64) / dtype_max

    # BaSiC illumination correction
    corrector = KCorrectGPU(
        use_gpu=use_gpu, verbose=False, device_id=device_id if use_gpu else None
    )
    flatfield, darkfield = corrector.fit(
        tiles_float,
        if_darkfield=True,
        max_iterations=BASIC_MAX_ITERATIONS,
        optimization_tolerance=BASIC_OPTIMIZATION_TOLERANCE,
    )

    flatfield_safe = np.clip(flatfield, BASIC_FLATFIELD_MIN, None)
    corrected = (tiles_float - darkfield) / flatfield_safe
    corrected = np.clip(corrected, 0, 1)
    corrected = (corrected * dtype_max).astype(np.uint16)

    # Get or compute stitch model
    output_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_zplane = n_zplanes // 2
    ch1_pkl = STITCH_DIR / f"cyc{CYCLE:02d}" / "CH1" / "result_df.pkl"

    if channel == 1 and zplane == ref_zplane:
        result_df, _ = stitch_images(
            corrected,
            rows,
            cols,
            initial_ncc_threshold=INITIAL_NCC_THRESHOLD,
            overlap_percentage=OVERLAP_PERCENTAGE,
            pou=POU,
            max_cores=4,
            use_gpu=use_gpu,
        )
        pkl_path = output_dir / "result_df.pkl"
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
        output_dtype=np.uint16,
    )

    # Save
    output_file = output_dir / f"{zplane:02d}.tif"
    imsave(str(output_file), stitched, check_contrast=False)

    return str(output_file), None


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
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

        if processed % 5 == 0:
            print(f"  Progress: {processed}/{n_zplanes} z-planes")

    ch_time = time.time() - ch_start
    print(
        f"  Channel {channel} complete: {processed}/{n_zplanes} z-planes "
        f"in {ch_time:.1f}s"
    )

    if errors:
        print(f"  Errors: {len(errors)}")
        for e in errors[:3]:
            print(f"    - {e}")

    # Generate QC image from middle z-plane
    qc_file = output_dir / f"{ref_zplane:02d}.tif"
    if qc_file.exists():
        log_dir = Path(snakemake.log[0]).parent
        log_dir.mkdir(parents=True, exist_ok=True)
        qc_data = imread(str(qc_file))
        qc_path = log_dir / f"cyc{CYCLE:02d}_CH{channel}_stitched.png"
        save_qc_image(qc_data, str(qc_path), f"Cycle {CYCLE} CH{channel} Z{ref_zplane}")

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

# Structured logging after-summary (mirrors utils.sh summary_after)
summary_after("stitch", PROJECT_DIR, total_time, exit_code=0)

# ---------------------------------------------------------------------------
# Write sentinel (replaces .complete marker)
# ---------------------------------------------------------------------------
sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(
    f"stage=stitch\n"
    f"cycle={CYCLE}\n"
    f"completed={datetime.now().isoformat()}\n"
    f"channels={START_CHANNEL}-{END_CHANNEL}\n"
    f"zplanes={n_zplanes}\n"
    f"duration_minutes={total_time/60:.1f}\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
