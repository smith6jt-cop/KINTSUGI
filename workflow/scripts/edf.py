"""
Snakemake wrapper for KINTSUGI EDF (replaces 04_edf.sh).

Performs extended depth of focus (variance-based z-projection) with
quality gating for one cycle (all channels).

Snakemake guarantees:
  - Deconvolution output exists before this script runs
  - No need for wait_for_input() polling

Per-channel skip-existing: If a sentinel is missing but some channels are
already complete (e.g. interrupted job), completed channels are skipped.
A channel is considered complete when its marker-named EDF output file exists.
"""

import gc
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime
from itertools import cycle as iter_cycle
from pathlib import Path

import numpy as np
from skimage.io import imread, imsave

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = list(snakemake.params.channels)

# Setup Python path
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

START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)

# EDF parameters from workflow config
edf_cfg = wf_config.get("edf", {})
EDF_PARAMS = {
    "radius_x": int(edf_cfg.get("radius_x", 2)),
    "radius_y": int(edf_cfg.get("radius_y", 2)),
    "sigma": float(edf_cfg.get("sigma", 10.0)),
    "tiles": tuple(edf_cfg.get("tiles", [3, 3])),
    "blend_depth": int(edf_cfg.get("blend_depth", 0)),
    "z_smooth_sigma": float(edf_cfg.get("z_smooth_sigma", 1.0)),
}

# Quality gate thresholds
qg_cfg = wf_config.get("quality_gate", {})
QUALITY_GATE = {
    "max_zero_pct": float(qg_cfg.get("max_zero_pct", 0.10)),
    "max_sat_pct": float(qg_cfg.get("max_sat_pct", 0.05)),
    "min_valid_slices": int(qg_cfg.get("min_valid_slices", 3)),
}

# Number of CPU cores available (from SLURM allocation)
CPUS = int(getattr(snakemake.resources, "cpus_per_task", 4))

# ---------------------------------------------------------------------------
# GPU initialization (respects device_mode from Snakemake rule)
# ---------------------------------------------------------------------------
REQUESTED_DEVICE_MODE = getattr(snakemake.params, "device_mode", "gpu")

DEVICE_MODE = "cpu"
n_gpus = 0
GPU_IDS = []
if REQUESTED_DEVICE_MODE == "gpu":
    try:
        import cupy as cp

        cp.cuda.Device(0).use()
        _ = cp.zeros(1)
        n_gpus = cp.cuda.runtime.getDeviceCount()
        GPU_IDS = list(range(n_gpus))
        print(f"CUDA initialized: {n_gpus} GPU(s) available")
        DEVICE_MODE = "gpu"
    except Exception as e:
        print(f"WARNING: CUDA initialization failed: {e}")
        print("Falling back to CPU mode")
else:
    print(f"CPU mode (requested by Snakemake rule)")

EDF_PARAMS["backend"] = "numpy" if DEVICE_MODE == "cpu" else "auto"
EDF_PARAMS["device"] = DEVICE_MODE

# ---------------------------------------------------------------------------
# Channel names for marker-named output
# ---------------------------------------------------------------------------
CHANNEL_NAMES_CONFIG = getattr(snakemake.params, "channel_names", {})


def get_channel_output_name(cycle, channel):
    """Get marker-named output filename, falling back to CH#_edf.tif."""
    cyc_key = int(cycle)
    ch_idx = int(channel) - 1
    names = CHANNEL_NAMES_CONFIG.get(cyc_key, CHANNEL_NAMES_CONFIG.get(str(cyc_key), []))
    if names and ch_idx < len(names):
        return names[ch_idx].replace("/", "_").replace(" ", "_") + ".tif"
    return f"CH{channel}_edf.tif"


# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from kintsugi.edf import EDFProcessor
from kintsugi.gpu import cleanup_gpu_memory

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data"
DECON_DIR = DATA_DIR / "processed" / "deconvolved"
EDF_DIR = DATA_DIR / "processed" / "edf"
EDF_DIR.mkdir(parents=True, exist_ok=True)

# Structured logging header (mirrors utils.sh)
log_header("edf", CYCLE, PROJECT_DIR)
summary_before("edf", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"EDF Processing - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Device mode: {DEVICE_MODE}")
print(f"Input: {DECON_DIR}")
print(f"Output: {EDF_DIR}")
print(f"EDF params: radius={EDF_PARAMS['radius_x']}, sigma={EDF_PARAMS['sigma']}")


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def load_stack_parallel(tiff_files):
    """Load z-stack from TIFF files with parallel I/O."""

    def load_single(f):
        return imread(str(f))

    with ThreadPoolExecutor(max_workers=CPUS) as executor:
        slices = list(executor.map(load_single, tiff_files))
    return np.stack(slices, axis=0)


def quality_gate_zstack(stack, channel_name):
    """Quality gate: detect and exclude corrupted z-slices."""
    valid_z_indices = []
    excluded_reasons = []

    for z_idx in range(stack.shape[0]):
        slice_data = stack[z_idx]
        total_pixels = slice_data.size

        pct_zero = np.sum(slice_data == 0) / total_pixels
        pct_sat = np.sum(slice_data >= 64000) / total_pixels

        if pct_zero > QUALITY_GATE["max_zero_pct"]:
            excluded_reasons.append(f"z={z_idx+1}: {pct_zero*100:.1f}% zeros")
            continue

        if pct_sat > QUALITY_GATE["max_sat_pct"]:
            excluded_reasons.append(f"z={z_idx+1}: {pct_sat*100:.1f}% saturated")
            continue

        valid_z_indices.append(z_idx)

    return valid_z_indices, excluded_reasons


def save_qc_image(data, output_path, title=""):
    """Save EDF result as QC image."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(8, 8))
        vmin, vmax = np.percentile(data, [1, 99])
        im = ax.imshow(data, cmap="gray", vmin=vmin, vmax=vmax)
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches="tight")
        plt.close()
        print(f"    QC saved: {Path(output_path).name}")
    except Exception as e:
        print(f"    QC save failed: {e}")


# ---------------------------------------------------------------------------
# Skip-existing helper
# ---------------------------------------------------------------------------
def channel_edf_complete(ch):
    """Check if EDF output file exists for this channel."""
    output_path = EDF_DIR / f"cyc{CYCLE:02d}"
    output_file = output_path / get_channel_output_name(CYCLE, ch)
    return output_file.exists()


# ---------------------------------------------------------------------------
# Per-channel EDF processing
# ---------------------------------------------------------------------------
def process_edf(args):
    """Process EDF for one channel with quality gate."""
    ch, gpu_id = args
    label = f"GPU{gpu_id}" if DEVICE_MODE == "gpu" else "CPU"
    print(f"\n  [{label}] Channel {ch} starting...")

    try:
        if DEVICE_MODE == "gpu" and n_gpus > 0:
            import cupy as cp

            cp.cuda.Device(gpu_id).use()

        processor = EDFProcessor(backend=EDF_PARAMS["backend"], method="variance")

        input_path = DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{ch}"
        output_path = EDF_DIR / f"cyc{CYCLE:02d}"
        output_path.mkdir(parents=True, exist_ok=True)

        input_files = sorted(input_path.glob("*.tif"))
        if not input_files:
            print(f"    No input files found for channel {ch}")
            return ch, False, "No input files"

        # Load z-stack
        print(f"    Loading {len(input_files)} z-slices...")
        stack = load_stack_parallel(input_files)
        n_zplanes_local = stack.shape[0]
        print(f"    Stack shape: {stack.shape}")

        # Quality gate
        valid_z_indices, excluded_reasons = quality_gate_zstack(stack, f"CH{ch}")

        if excluded_reasons:
            print(f"    Quality gate excluded z-slices:")
            for reason in excluded_reasons:
                print(f"      - {reason}")

        if len(valid_z_indices) < QUALITY_GATE["min_valid_slices"]:
            msg = (
                f"Only {len(valid_z_indices)} valid z-slices, "
                f"need at least {QUALITY_GATE['min_valid_slices']}"
            )
            print(f"    ERROR: {msg}")
            return ch, False, msg

        stack_filtered = stack[valid_z_indices]
        print(f"    Valid z-slices: {len(valid_z_indices)}/{n_zplanes_local}")

        # Z bounds
        if len(valid_z_indices) > 4:
            z_start = 2
            z_end = len(valid_z_indices) - 1
        else:
            z_start = 1
            z_end = len(valid_z_indices)

        print(f"    Using z_start={z_start}, z_end={z_end}")

        # Run EDF
        edf_result = processor.process(
            stack_filtered,
            radius_x=EDF_PARAMS["radius_x"],
            radius_y=EDF_PARAMS["radius_y"],
            sigma=EDF_PARAMS["sigma"],
            z_start=z_start,
            z_end=z_end,
            tiles=EDF_PARAMS["tiles"],
            device=EDF_PARAMS["device"],
            blend_depth=EDF_PARAMS.get("blend_depth", 0),
            z_smooth_sigma=EDF_PARAMS.get("z_smooth_sigma", 0.0),
        )

        # Save (marker-named if channel names available)
        output_file = output_path / get_channel_output_name(CYCLE, ch)
        imsave(str(output_file), edf_result.astype(np.uint16), check_contrast=False)

        pct_zero = 100 * np.sum(edf_result == 0) / edf_result.size
        print(f"    Output: {pct_zero:.2f}% zeros (should be <5%)")
        if pct_zero > 5:
            print(f"    WARNING: High zero percentage in output")

        # QC image
        log_dir = Path(snakemake.log[0]).parent
        qc_path = log_dir / f"cyc{CYCLE:02d}_CH{ch}_edf.png"
        save_qc_image(
            edf_result, str(qc_path), f"Cycle {CYCLE} CH{ch} EDF ({pct_zero:.1f}% zeros)"
        )

        print(f"  [{label}] Channel {ch} complete -> {output_file.name}")

        del stack, stack_filtered, edf_result
        gc.collect()

        return ch, True, None

    except Exception as e:
        print(f"  [{label}] Channel {ch} FAILED: {e}")
        import traceback

        traceback.print_exc()
        return ch, False, str(e)
    finally:
        cleanup_gpu_memory()


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------
start_time = time.time()

# Filter out channels that are already complete
channels_to_process = []
skipped_channels = []
for ch in CHANNELS:
    if channel_edf_complete(ch):
        print(f"  Channel {ch} SKIPPED ({get_channel_output_name(CYCLE, ch)} exists)")
        skipped_channels.append(ch)
    else:
        channels_to_process.append(ch)

if not channels_to_process:
    print(f"\nAll {len(CHANNELS)} channels already complete — nothing to do")
    results = []
elif DEVICE_MODE == "cpu":
    print(f"\n[CPU] Processing {len(channels_to_process)} channels sequentially")
    results = [process_edf((ch, 0)) for ch in channels_to_process]
elif n_gpus >= 2 and len(channels_to_process) >= 2:
    print(f"\n[MULTI-GPU] {len(channels_to_process)} channels across {n_gpus} GPUs")
    pairs = list(zip(channels_to_process, iter_cycle(GPU_IDS)))
    with ThreadPoolExecutor(max_workers=n_gpus) as executor:
        results = list(executor.map(process_edf, pairs))
else:
    print(f"\n[GPU] Processing {len(channels_to_process)} channels on GPU 0")
    results = [process_edf((ch, 0)) for ch in channels_to_process]

successful = sum(1 for _, ok, _ in results if ok) + len(skipped_channels)
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"EDF Complete")
print(f"{'='*60}")
print(f"Success: {successful}/{len(CHANNELS)} channels ({len(skipped_channels)} skipped)")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {EDF_DIR}")

failed = [(ch, err) for ch, ok, err in results if not ok]
if failed:
    print(f"\nFailed channels:")
    for ch, err in failed:
        print(f"  CH{ch}: {err}")
    summary_after("edf", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    sys.exit(1)

# Structured logging after-summary (mirrors utils.sh summary_after)
summary_after("edf", PROJECT_DIR, elapsed * 60, exit_code=0)

# ---------------------------------------------------------------------------
# Write sentinel
# ---------------------------------------------------------------------------
sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(
    f"stage=edf\n"
    f"cycle={CYCLE}\n"
    f"completed={datetime.now().isoformat()}\n"
    f"channels={START_CHANNEL}-{END_CHANNEL}\n"
    f"successful={successful}\n"
    f"skipped={len(skipped_channels)}\n"
    f"duration_minutes={elapsed:.1f}\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
