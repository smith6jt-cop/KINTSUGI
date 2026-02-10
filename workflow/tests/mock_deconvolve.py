"""
Mock deconvolve wrapper for testing Snakemake workflow mechanics.

Produces dummy deconvolved TIFFs and QC images in the expected directory
structure without requiring GPU/KINTSUGI dependencies.
"""

import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

try:
    import tifffile

    def _write_tiff(path, data):
        tifffile.imwrite(str(path), data)
except ImportError:
    def _write_tiff(path, data):
        path.write_bytes(data.tobytes())

# ---------------------------------------------------------------------------
# Snakemake interface (identical to real deconvolve.py)
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = list(snakemake.params.channels)

# Logging utilities
sys.path.insert(0, str(Path(snakemake.scriptdir).parent / "scripts"))
from log_utils import log_header, log_footer, summary_before, summary_after

START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)

# ---------------------------------------------------------------------------
# Structured logging
# ---------------------------------------------------------------------------
log_header("decon", CYCLE, PROJECT_DIR)
summary_before("decon", PROJECT_DIR)

STITCH_DIR = PROJECT_DIR / "data" / "processed" / "stitched"
DECON_DIR = PROJECT_DIR / "data" / "processed" / "deconvolved"

print(f"\n{'='*60}")
print(f"[MOCK] Deconvolution - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")

# ---------------------------------------------------------------------------
# Failure injection for testing
# ---------------------------------------------------------------------------
fail_flag = PROJECT_DIR / "slurm" / "logs" / "snakemake" / f"_mock_fail_decon_cyc{CYCLE:02d}"
if fail_flag.exists():
    fail_flag.unlink()
    print("MOCK FAILURE INJECTED")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Produce dummy outputs
# ---------------------------------------------------------------------------
start_time = time.time()

for channel in CHANNELS:
    # Read stitched input to determine z-plane count
    stitch_ch_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    z_files = sorted(stitch_ch_dir.glob("*.tif")) if stitch_ch_dir.exists() else []
    n_zplanes = max(len(z_files), 1)

    decon_ch_dir = DECON_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    decon_ch_dir.mkdir(parents=True, exist_ok=True)

    for z in range(1, n_zplanes + 1):
        dummy = np.random.randint(0, 1000, (64, 64), dtype=np.uint16)
        output_file = decon_ch_dir / f"{z:02d}.tif"
        _write_tiff(output_file, dummy)

    print(f"  Channel {channel}: {n_zplanes} z-planes written")

    # QC image
    log_dir = Path(snakemake.log[0]).parent
    log_dir.mkdir(parents=True, exist_ok=True)
    qc_path = log_dir / f"cyc{CYCLE:02d}_CH{channel}_decon.png"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(np.random.rand(64, 64), cmap="gray")
        ax.set_title(f"[MOCK] Decon Cyc{CYCLE} CH{channel}")
        ax.axis("off")
        plt.savefig(str(qc_path), dpi=50)
        plt.close()
        print(f"  QC saved: {qc_path.name}")
    except Exception as e:
        print(f"  QC save failed: {e}")

elapsed = time.time() - start_time

print(f"\nTime: {elapsed:.1f}s")
summary_after("decon", PROJECT_DIR, elapsed, exit_code=0)

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
    f"successful={END_CHANNEL - START_CHANNEL + 1}\n"
    f"duration_minutes={elapsed/60:.1f}\n"
    f"mock=true\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
