"""
Snakemake wrapper for KINTSUGI stitching (replaces 02_stitching.sh).

Performs BaSiC illumination correction and GPU-accelerated tile stitching
for one cycle (all channels, all z-planes).

Snakemake guarantees:
  - Input raw directory exists before this script runs
  - Output sentinel is only written if this script exits successfully
  - No need for wait_for_input() polling (Snakemake tracks dependencies)

Per-channel skip-existing: If a sentinel is missing but some channels are
already complete (e.g. interrupted job), completed channels are skipped.
A channel is considered complete only when ALL expected z-plane TIFs exist.
"""

import gc
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
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities (replicates slurm/utils.sh structured logging)
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_header, log_footer, summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration from config.yaml (single source of truth)
# ---------------------------------------------------------------------------
wf_config = snakemake.config

# Tile grid — no fallback defaults; config.yaml must provide these
TILE_ROWS = int(wf_config["tile_rows"])
TILE_COLS = int(wf_config["tile_cols"])
TILE_OVERLAP = float(wf_config["tile_overlap"])
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
INITIAL_NCC_THRESHOLD = float(wf_config["ncc_threshold"])
POU = float(wf_config["pou"])
BLEND_SIGMA = float(wf_config.get("blend_sigma", 15.0))

# Boundary QC + multi-z fallback (detects and corrects rare misalignments)
BOUNDARY_QC_ENABLED = bool(wf_config.get("boundary_qc_enabled", True))
BOUNDARY_QC_THRESHOLD = float(wf_config.get("boundary_qc_threshold", 0.65))
BOUNDARY_QC_MAX_ALTERNATES = int(wf_config.get("boundary_qc_max_alternates", 4))

# Number of CPU cores available (from SLURM allocation)
CPUS = int(getattr(snakemake.resources, "cpus_per_task", 4))

# ---------------------------------------------------------------------------
# GPU initialization (respects device_mode from Snakemake rule)
# ---------------------------------------------------------------------------
REQUESTED_DEVICE_MODE = getattr(snakemake.params, "device_mode", "gpu")

use_gpu = False
if REQUESTED_DEVICE_MODE == "gpu":
    try:
        import cupy as cp

        cp.cuda.Device(0).use()
        _ = cp.zeros(1)
        print(f"CUDA initialized successfully on device 0")
        use_gpu = True
    except Exception as e:
        print(f"WARNING: CUDA initialization failed: {e}")
        print("Falling back to CPU processing")
else:
    print(f"CPU mode (requested by Snakemake rule)")

# ---------------------------------------------------------------------------
# Imports (after PYTHONPATH setup)
# ---------------------------------------------------------------------------
from Kstitch.stitching import stitch_images
from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending
from kintsugi.gpu import cleanup_gpu_memory
from kintsugi.qc.boundary_check import (
    check_boundary_quality,
    select_alternate_zplanes,
)

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


def load_tiles_parallel(file_list, max_workers=None):
    """Load tiles in parallel."""
    if max_workers is None:
        max_workers = CPUS

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
def _load_and_correct_tiles(zplane, channel, device_id=0):
    """Load + BaSiC-correct all tiles for one (zplane, channel).

    Returns (corrected_uint16, error_or_None).
    """
    pattern = f"*_Z{zplane:03d}_CH{channel}.tif"
    tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)
    if not tile_files:
        pattern = f"*_Z{zplane:02d}_CH{channel}.tif"
        tile_files = sorted(cycle_dir.glob(pattern), key=alphanumeric_key)
    if not tile_files:
        return None, f"No files for Z{zplane} CH{channel}"
    if len(tile_files) != n_tiles:
        print(f"  Warning: Expected {n_tiles} tiles, found {len(tile_files)}")

    tiles = load_tiles_parallel(tile_files)
    dtype_max = (
        np.iinfo(tiles.dtype).max
        if np.issubdtype(tiles.dtype, np.integer)
        else 1.0
    )
    tiles_float = tiles.astype(np.float64) / dtype_max

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
    return corrected, None


def _stitch_with_model(corrected, result_df):
    """Apply an existing stitch model to corrected tiles, return stitched image."""
    result_df = result_df.copy()
    result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
    result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()
    overlap_fraction = (TILE_OVERLAP, TILE_OVERLAP)
    return stitch_with_blending(
        corrected,
        result_df,
        blend=True,
        sigma=BLEND_SIGMA,
        overlap_fraction=overlap_fraction,
        output_dtype=np.uint16,
    )


def _compute_stitch_model(corrected):
    """Run phase correlation + spanning tree to produce a stitch model."""
    result_df, _ = stitch_images(
        corrected,
        rows,
        cols,
        initial_ncc_threshold=INITIAL_NCC_THRESHOLD,
        overlap_percentage=OVERLAP_PERCENTAGE,
        pou=POU,
        max_cores=CPUS,
        use_gpu=use_gpu,
    )
    return result_df


def process_zplane(cycle, channel, zplane, device_id=0, model_df=None):
    """Process a single z-plane: load tiles, correct, stitch, save.

    If ``model_df`` is provided, that stitch model is used (no recomputation).
    Otherwise, for (channel=1, zplane=ref_zplane) a model is computed and
    saved; for other z-planes the saved CH1 model is loaded.
    """
    corrected, error = _load_and_correct_tiles(zplane, channel, device_id=device_id)
    if error:
        return None, error

    output_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    output_dir.mkdir(parents=True, exist_ok=True)

    ref_zplane = n_zplanes // 2
    ch1_pkl = STITCH_DIR / f"cyc{CYCLE:02d}" / "CH1" / "result_df.pkl"

    if model_df is not None:
        result_df = model_df
    elif channel == 1 and zplane == ref_zplane:
        result_df = _compute_stitch_model(corrected)
        pkl_path = output_dir / "result_df.pkl"
        result_df.to_pickle(str(pkl_path))
    elif ch1_pkl.exists():
        result_df = pd.read_pickle(str(ch1_pkl))
    else:
        return None, f"No stitch model at {ch1_pkl}"

    stitched = _stitch_with_model(corrected, result_df)

    output_file = output_dir / f"{zplane:02d}.tif"
    imsave(str(output_file), stitched, check_contrast=False)

    return str(output_file), None


def _run_boundary_qc_with_fallback(cycle):
    """Check reference-z-plane boundaries; try alternate z-planes if needed.

    Must be called AFTER the initial ref_zplane has been stitched with the
    initial model. Modifies the saved result_df.pkl if a better model is
    found, and re-stitches the reference z-plane.

    Writes a JSON log to data/processed/stitched/cyc{NN}/CH1/boundary_qc.json.
    """
    import json

    ref_zplane = n_zplanes // 2
    ch1_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / "CH1"
    ch1_pkl = ch1_dir / "result_df.pkl"
    ref_tif = ch1_dir / f"{ref_zplane:02d}.tif"

    if not ref_tif.exists() or not ch1_pkl.exists():
        print("  BOUNDARY QC: skipped (reference stitch output missing)")
        return

    # Load the reference stitched image and the initial model
    initial_df = pd.read_pickle(str(ch1_pkl))
    ref_stitched = imread(str(ref_tif))

    # Need tile dimensions — get them from a raw tile
    sample_files = sorted(
        cycle_dir.glob(f"*_Z{ref_zplane:03d}_CH1.tif"), key=alphanumeric_key
    )
    if not sample_files:
        sample_files = sorted(
            cycle_dir.glob(f"*_Z{ref_zplane:02d}_CH1.tif"), key=alphanumeric_key
        )
    if not sample_files:
        print("  BOUNDARY QC: skipped (cannot determine tile size)")
        return
    sample_tile = imread(str(sample_files[0]))
    tile_h_px, tile_w_px = sample_tile.shape

    initial_passed, initial_flagged, initial_scores = check_boundary_quality(
        ref_stitched, initial_df, tile_h_px, tile_w_px,
        threshold=BOUNDARY_QC_THRESHOLD,
    )

    qc_log = {
        "ref_zplane": ref_zplane,
        "threshold": BOUNDARY_QC_THRESHOLD,
        "initial": {
            "n_flagged": len(initial_flagged),
            "worst_score": min(
                [f["score"] for f in initial_flagged], default=None,
            ),
            "flagged": [
                {
                    "orientation": f["orientation"],
                    "tile_a": list(f["tile_a"]),
                    "tile_b": list(f["tile_b"]),
                    "score": f["score"],
                } for f in initial_flagged[:10]
            ],
        },
        "alternates": [],
        "final_model_zplane": ref_zplane,
        "final_n_flagged": len(initial_flagged),
    }

    if initial_passed:
        print(f"  BOUNDARY QC: all {len(initial_scores)} boundaries passed")
        with open(ch1_dir / "boundary_qc.json", "w") as fh:
            json.dump(qc_log, fh, indent=2)
        return

    print(
        f"  BOUNDARY QC: {len(initial_flagged)} boundaries flagged "
        f"(worst ratio={initial_flagged[0]['score']:.3f}); trying alternate z-planes"
    )
    for f in initial_flagged[:5]:
        print(
            f"    {f['orientation']:10s} {f['tile_a']} | {f['tile_b']}  "
            f"ratio={f['score']:.3f}"
        )

    best_df = initial_df
    best_n_flagged = len(initial_flagged)
    best_worst_score = initial_flagged[0]["score"]
    best_zplane = ref_zplane

    alternates = select_alternate_zplanes(
        ref_zplane, n_zplanes, BOUNDARY_QC_MAX_ALTERNATES,
    )
    for alt_z in alternates:
        print(f"  Trying alternate Z{alt_z}...")
        alt_corrected, err = _load_and_correct_tiles(alt_z, channel=1)
        if err:
            print(f"    skipped: {err}")
            continue
        alt_df = _compute_stitch_model(alt_corrected)
        alt_stitched = _stitch_with_model(alt_corrected, alt_df)
        alt_passed, alt_flagged, alt_scores = check_boundary_quality(
            alt_stitched, alt_df, tile_h_px, tile_w_px,
            threshold=BOUNDARY_QC_THRESHOLD,
        )
        alt_worst = min(
            [f["score"] for f in alt_flagged], default=1.0,
        )
        qc_log["alternates"].append({
            "zplane": alt_z,
            "passed": alt_passed,
            "n_flagged": len(alt_flagged),
            "worst_score": alt_worst,
        })
        print(
            f"    Z{alt_z}: {len(alt_flagged)} flagged, "
            f"worst={alt_worst:.3f}"
        )
        # Adopt alt if it has strictly fewer flagged, or equal flagged but
        # higher worst score (= less severe)
        if (
            len(alt_flagged) < best_n_flagged
            or (len(alt_flagged) == best_n_flagged and alt_worst > best_worst_score)
        ):
            best_df = alt_df
            best_n_flagged = len(alt_flagged)
            best_worst_score = alt_worst
            best_zplane = alt_z
            if alt_passed:
                print(f"    → Z{alt_z} passes cleanly, stopping search")
                break

    qc_log["final_model_zplane"] = best_zplane
    qc_log["final_n_flagged"] = best_n_flagged

    if best_zplane != ref_zplane:
        print(
            f"  BOUNDARY QC: adopting model from Z{best_zplane} "
            f"({best_n_flagged} flagged vs initial {len(initial_flagged)})"
        )
        best_df.to_pickle(str(ch1_pkl))
        # Re-stitch the reference z-plane with the new model
        result, error = process_zplane(CYCLE, 1, ref_zplane, model_df=best_df)
        if error:
            print(f"  WARNING: re-stitch of Z{ref_zplane} failed: {error}")
    else:
        print(
            f"  BOUNDARY QC: no alternate improved on initial model; "
            f"keeping Z{ref_zplane}"
        )

    with open(ch1_dir / "boundary_qc.json", "w") as fh:
        json.dump(qc_log, fh, indent=2)


# ---------------------------------------------------------------------------
# Skip-existing helper
# ---------------------------------------------------------------------------
def channel_complete(channel):
    """Check if all z-planes are already stitched for this channel."""
    ch_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    if not ch_dir.exists():
        return False
    for z in range(1, n_zplanes + 1):
        if not (ch_dir / f"{z:02d}.tif").exists():
            return False
    # CH1 also needs the stitch model pickle
    if channel == 1 and not (ch_dir / "result_df.pkl").exists():
        return False
    return True


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
start_time = time.time()
channels = list(range(START_CHANNEL, END_CHANNEL + 1))
ref_zplane = n_zplanes // 2
skipped_channels = 0

for channel in channels:
    if channel_complete(channel):
        print(f"\n--- Channel {channel}/{END_CHANNEL} --- SKIPPED (all {n_zplanes} z-planes exist)")
        skipped_channels += 1
        continue

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

        # Post-stitch boundary QC: detect misaligned boundaries and try
        # alternate z-planes if the initial model produces visible ghosting.
        if BOUNDARY_QC_ENABLED and n_zplanes >= 2:
            _run_boundary_qc_with_fallback(CYCLE)

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

processed_channels = len(channels) - skipped_channels

print(f"\n{'='*60}")
print(f"Correction + Stitching Complete")
print(f"{'='*60}")
print(f"Cycle: {CYCLE}")
print(f"Channels: {len(channels)} ({processed_channels} processed, {skipped_channels} skipped)")
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
    f"skipped={skipped_channels}\n"
    f"duration_minutes={total_time/60:.1f}\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
