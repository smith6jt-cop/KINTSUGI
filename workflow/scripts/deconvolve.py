"""
Snakemake wrapper for KINTSUGI deconvolution (replaces 03_deconvolution.sh).

Performs Richardson-Lucy deconvolution for one cycle (all channels).

Snakemake guarantees:
  - Stitching output exists before this script runs (via sentinel dependency)
  - No need for wait_for_input() polling
  - No need for skip-existing checks
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = list(snakemake.params.channels)

# Setup Python path
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities (replicates slurm/utils.sh structured logging)
sys.path.insert(0, str(Path(snakemake.workflow.basedir) / "scripts"))
from log_utils import log_header, log_footer, summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
experiment_config = {}
config_path = PROJECT_DIR / "meta" / "experiment.json"
if config_path.exists():
    with open(config_path) as f:
        experiment_config = json.load(f)

START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)
OUTPUT_FORMAT = snakemake.config.get("output_format", "zarr")

# Microscope parameters (experiment.json -> defaults)
XY_VOX = float(experiment_config.get("xy_pixel_size", 377))
Z_VOX = float(experiment_config.get("z_step_size", 1500))
MIC_NA = float(experiment_config.get("numerical_aperture", 0.75))
TISSUE_RI = float(experiment_config.get("tissue_refractive_index", 1.44))

# Parse wavelengths
WAVELENGTHS = {}
if "wavelengths" in experiment_config:
    for ch_str, (ex, em) in experiment_config["wavelengths"].items():
        WAVELENGTHS[int(ch_str)] = (float(ex), float(em))
else:
    # Default wavelengths
    WAVELENGTHS = {
        1: (358.0, 461.0),
        2: (753.0, 775.0),
        3: (560.0, 575.0),
        4: (648.0, 668.0),
    }

# ---------------------------------------------------------------------------
# GPU initialization
# ---------------------------------------------------------------------------
DEVICE_MODE = "gpu"
try:
    import cupy as cp

    cp.cuda.Device(0).use()
    _ = cp.zeros(1)
    print(f"CUDA initialized successfully on device 0")
except Exception as e:
    print(f"WARNING: CUDA initialization failed: {e}")
    print("Falling back to CPU processing")
    DEVICE_MODE = "cpu"

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
from Kdecon import decon
from kintsugi.gpu import cleanup_gpu_memory

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data"
STITCH_DIR = DATA_DIR / "processed" / "stitched"
DECON_DIR = DATA_DIR / "processed" / "deconvolved"
DECON_DIR.mkdir(parents=True, exist_ok=True)

# Structured logging header (mirrors utils.sh)
log_header("decon", CYCLE, PROJECT_DIR)
summary_before("decon", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"Deconvolution - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Device mode: {DEVICE_MODE}")
print(f"Output: {OUTPUT_FORMAT}")
print(f"Voxel size: {XY_VOX}x{XY_VOX}x{Z_VOX} nm")
print(f"NA: {MIC_NA}, RI: {TISSUE_RI}")


# ---------------------------------------------------------------------------
# QC helper
# ---------------------------------------------------------------------------
def save_qc_slice(data, output_path, title=""):
    """Save middle slice as QC image."""
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        if data.ndim == 3:
            slice_data = data[data.shape[0] // 2]
        else:
            slice_data = data

        fig, ax = plt.subplots(figsize=(8, 8))
        im = ax.imshow(slice_data, cmap="gray")
        ax.set_title(title)
        ax.axis("off")
        plt.colorbar(im, ax=ax, fraction=0.046)
        plt.tight_layout()
        plt.savefig(output_path, dpi=100, bbox_inches="tight")
        plt.close()
        print(f"  QC saved: {output_path}")
    except Exception as e:
        print(f"  QC save failed: {e}")


# ---------------------------------------------------------------------------
# Per-channel processing
# ---------------------------------------------------------------------------
def process_channel(ch):
    """Deconvolve one channel."""
    print(f"\n  Channel {ch} starting...")

    try:
        lambda_ex, lambda_em = WAVELENGTHS.get(ch, (560, 575))

        decon(
            base_dir=str(PROJECT_DIR),
            stitch_dir=str(STITCH_DIR),
            dec_cycle=CYCLE,
            dec_channel=ch,
            xy_vox=XY_VOX,
            z_vox=Z_VOX,
            iterations=25,
            mic_NA=MIC_NA,
            tissue_RI=TISSUE_RI,
            damping=0,
            stop_criterion=5.0,
            device=DEVICE_MODE,
            device_id=None,
            wavelengths=WAVELENGTHS,
            decon_dir=str(DECON_DIR),
        )

        # Generate QC image
        output_files = list(
            (DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{ch}").glob("*.tif")
        )
        if output_files:
            from skimage.io import imread

            mid_idx = len(output_files) // 2
            data = imread(str(output_files[mid_idx]))

            # Write QC to the log directory (alongside snakemake log)
            log_dir = Path(snakemake.log[0]).parent
            qc_path = log_dir / f"cyc{CYCLE:02d}_CH{ch}_decon.png"
            save_qc_slice(data, str(qc_path), f"Cycle {CYCLE} Channel {ch} - Deconvolved")

        print(f"  Channel {ch} complete")
        return ch, True

    except Exception as e:
        print(f"  Channel {ch} FAILED: {e}")
        import traceback

        traceback.print_exc()
        return ch, False
    finally:
        cleanup_gpu_memory()


# ---------------------------------------------------------------------------
# Main processing loop
# ---------------------------------------------------------------------------
start_time = time.time()

print(f"\nProcessing {len(CHANNELS)} channels sequentially")
results = [process_channel(ch) for ch in CHANNELS]

successful = sum(1 for _, ok in results if ok)
elapsed = (time.time() - start_time) / 60

print(f"\n{'='*60}")
print(f"Deconvolution Complete")
print(f"{'='*60}")
print(f"Success: {successful}/{len(CHANNELS)} channels")
print(f"Time: {elapsed:.1f} minutes")
print(f"Output: {DECON_DIR}")

if successful < len(CHANNELS):
    failed = [ch for ch, ok in results if not ok]
    print(f"\nFailed channels: {failed}")
    summary_after("decon", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    sys.exit(1)

# Structured logging after-summary (mirrors utils.sh summary_after)
summary_after("decon", PROJECT_DIR, elapsed * 60, exit_code=0)

# ---------------------------------------------------------------------------
# Write sentinel
# ---------------------------------------------------------------------------
sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(
    f"stage=decon\n"
    f"cycle={CYCLE}\n"
    f"completed={datetime.now().isoformat()}\n"
    f"channels={START_CHANNEL}-{END_CHANNEL}\n"
    f"successful={successful}\n"
    f"duration_minutes={elapsed:.1f}\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
