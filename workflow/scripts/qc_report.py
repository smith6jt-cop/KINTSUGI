"""
Snakemake wrapper for KINTSUGI QC report generation.

Dispatches to Kprocess QC functions based on snakemake.params.stage:
  - stitch -> run_stitched_qc()
  - decon  -> run_decon_qc()
  - edf    -> run_edf_qc()

Produces summary heatmaps, z-plane profile PDFs, and cross-stage
comparison plots. Each stage caches its statistics DataFrame as a
pickle for downstream comparison.
"""

# Headless rendering — must precede any matplotlib/Kprocess import
import matplotlib
matplotlib.use("Agg")

import os
import pickle
import sys
import time
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
STAGE = snakemake.params.stage
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLES = list(snakemake.params.cycles)
CHANNELS = list(snakemake.params.channels)
N_ZPLANES = int(snakemake.params.n_zplanes)
CHANNEL_NAMES = dict(snakemake.params.channel_names)

START_CYCLE = min(CYCLES)
END_CYCLE = max(CYCLES)
START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)

# Setup Python path (same as other wrapper scripts)
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_footer, log_info, log_warn

# ---------------------------------------------------------------------------
# Output directories
# ---------------------------------------------------------------------------
QC_DIR = PROJECT_DIR / "qc_plots"
QC_DIR.mkdir(parents=True, exist_ok=True)

CACHE_DIR = PROJECT_DIR / "cache"
CACHE_DIR.mkdir(parents=True, exist_ok=True)

DATA_DIR = PROJECT_DIR / "data" / "processed"

# ---------------------------------------------------------------------------
# GPU initialization (with CPU fallback)
# ---------------------------------------------------------------------------
DEVICE_MODE = "cpu"
GPU_IDS = []
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

# ---------------------------------------------------------------------------
# Import Kprocess QC functions
# ---------------------------------------------------------------------------
from Kprocess import run_stitched_qc, run_decon_qc, run_edf_qc, run_registration_qc

# ---------------------------------------------------------------------------
# Cache helpers
# ---------------------------------------------------------------------------
CACHE_NAMES = {
    "stitch": "stitch_stats.pkl",
    "decon": "decon_stats.pkl",
    "edf": "edf_stats.pkl",
    "registration": "registration_stats.pkl",
}


def load_cache(stage_name):
    """Load cached stats DataFrame for a previous stage, or None."""
    cache_path = CACHE_DIR / CACHE_NAMES.get(stage_name, "")
    if cache_path.exists():
        log_info(f"Loading cached stats from {cache_path}")
        with open(cache_path, "rb") as f:
            return pickle.load(f)
    return None


# ---------------------------------------------------------------------------
# Build channel_name_dict (int keys -> list of names)
# ---------------------------------------------------------------------------
channel_name_dict = {}
for k, v in CHANNEL_NAMES.items():
    channel_name_dict[int(k)] = list(v) if v else []

# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
print(f"\n{'='*60}")
print(f"QC Report - Stage: {STAGE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Cycles: {START_CYCLE}-{END_CYCLE}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Z-planes: {N_ZPLANES}")
print(f"Device mode: {DEVICE_MODE}")
print(f"Output: {QC_DIR}")

start_time = time.time()

gpu_ids = GPU_IDS if GPU_IDS else None

if STAGE == "stitch":
    stitch_dir = DATA_DIR / "stitched"
    cache_file = str(CACHE_DIR / CACHE_NAMES["stitch"])

    log_info(f"Running stitched QC on {stitch_dir}")
    stats_df = run_stitched_qc(
        stitch_dir=str(stitch_dir),
        cache_file=cache_file,
        start_cycle=START_CYCLE,
        end_cycle=END_CYCLE,
        start_channel=START_CHANNEL,
        end_channel=END_CHANNEL,
        n_zplanes=N_ZPLANES,
        qc_output_dir=str(QC_DIR),
        raw_stats_df=None,
        gpu_device_ids=gpu_ids,
        channel_name_dict=channel_name_dict,
    )

elif STAGE == "decon":
    decon_dir = DATA_DIR / "deconvolved"
    cache_file = str(CACHE_DIR / CACHE_NAMES["decon"])

    # Load stitched stats for cross-stage comparison
    stitched_stats_df = load_cache("stitch")

    log_info(f"Running deconvolved QC on {decon_dir}")
    stats_df = run_decon_qc(
        decon_dir=str(decon_dir),
        cache_file=cache_file,
        start_cycle=START_CYCLE,
        end_cycle=END_CYCLE,
        start_channel=START_CHANNEL,
        end_channel=END_CHANNEL,
        n_zplanes=N_ZPLANES,
        qc_output_dir=str(QC_DIR),
        stitched_stats_df=stitched_stats_df,
        gpu_device_ids=gpu_ids,
        channel_name_dict=channel_name_dict,
    )

elif STAGE == "edf":
    cache_file = str(CACHE_DIR / CACHE_NAMES["edf"])

    # Load decon stats for cross-stage comparison
    decon_stats_df = load_cache("decon")

    log_info(f"Running EDF QC")
    stats_df = run_edf_qc(
        project_dir=str(PROJECT_DIR),
        cache_file=cache_file,
        start_cycle=START_CYCLE,
        end_cycle=END_CYCLE,
        n_zplanes=N_ZPLANES,
        qc_output_dir=str(QC_DIR),
        edf_dir=str(DATA_DIR / "edf"),
        decon_stats_df=decon_stats_df,
        gpu_device_ids=gpu_ids,
        channel_name_dict=channel_name_dict,
    )

elif STAGE == "registration":
    cache_file = str(CACHE_DIR / CACHE_NAMES["registration"])

    log_info("Running registration QC")
    stats_df = run_registration_qc(
        project_dir=str(PROJECT_DIR),
        cache_file=cache_file,
        start_cycle=START_CYCLE,
        end_cycle=END_CYCLE,
        qc_output_dir=str(QC_DIR),
        channel_name_dict=channel_name_dict,
    )

else:
    print(f"ERROR: Unknown QC stage '{STAGE}'")
    log_footer(1)
    sys.exit(1)

elapsed = (time.time() - start_time) / 60

# ---------------------------------------------------------------------------
# Write sentinel
# ---------------------------------------------------------------------------
sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(
    f"stage=qc_{STAGE}\n"
    f"completed={datetime.now().isoformat()}\n"
    f"cycles={START_CYCLE}-{END_CYCLE}\n"
    f"channels={START_CHANNEL}-{END_CHANNEL}\n"
    f"duration_minutes={elapsed:.1f}\n"
)

print(f"\n{'='*60}")
print(f"QC Report Complete - {STAGE}")
print(f"{'='*60}")
print(f"Time: {elapsed:.1f} minutes")
print(f"Sentinel: {sentinel}")
log_footer(0)
