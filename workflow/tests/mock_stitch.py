"""
Mock stitch wrapper for testing Snakemake workflow mechanics.

Produces dummy stitched TIFFs and QC images in the expected directory
structure without requiring GPU/KINTSUGI dependencies.
"""

import re
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
# Snakemake interface (identical to real stitch.py)
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = list(snakemake.params.channels)

# Logging utilities
sys.path.insert(0, str(Path(snakemake.scriptdir).parent / "scripts"))
from log_utils import log_header, log_footer, summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
START_CHANNEL = min(CHANNELS)
END_CHANNEL = max(CHANNELS)

# Detect z-planes from raw data
RAW_DIR = PROJECT_DIR / "data" / "raw"
cycle_dir = None
for pattern in [f"cyc{CYCLE:03d}_*", f"cyc{CYCLE:03d}", f"cyc{CYCLE:02d}"]:
    matches = list(RAW_DIR.glob(pattern))
    if matches:
        cycle_dir = matches[0]
        break
if cycle_dir is None:
    cycle_dir = RAW_DIR / f"cyc{CYCLE:02d}"

sample_files = sorted(cycle_dir.glob("*_Z*_CH1.tif"))
z_planes = set()
for f in sample_files:
    match = re.search(r"_Z(\d+)_", f.name)
    if match:
        z_planes.add(int(match.group(1)))
n_zplanes = max(len(z_planes), 1)

# ---------------------------------------------------------------------------
# Structured logging
# ---------------------------------------------------------------------------
log_header("stitch", CYCLE, PROJECT_DIR)
summary_before("stitch", PROJECT_DIR)

STITCH_DIR = PROJECT_DIR / "data" / "processed" / "stitched"

print(f"\n{'='*60}")
print(f"[MOCK] Correction + Stitching - Cycle {CYCLE}")
print(f"{'='*60}")
print(f"Channels: {START_CHANNEL}-{END_CHANNEL}")
print(f"Z-planes: {n_zplanes}")

# ---------------------------------------------------------------------------
# Failure injection for testing
# ---------------------------------------------------------------------------
fail_flag = PROJECT_DIR / "slurm" / "logs" / "snakemake" / f"_mock_fail_stitch_cyc{CYCLE:02d}"
if fail_flag.exists():
    fail_flag.unlink()
    print("MOCK FAILURE INJECTED")
    sys.exit(1)

# ---------------------------------------------------------------------------
# Produce dummy outputs
# ---------------------------------------------------------------------------
start_time = time.time()

for channel in CHANNELS:
    ch_dir = STITCH_DIR / f"cyc{CYCLE:02d}" / f"CH{channel}"
    ch_dir.mkdir(parents=True, exist_ok=True)

    for z in range(1, n_zplanes + 1):
        dummy = np.random.randint(0, 1000, (64, 64), dtype=np.uint16)
        output_file = ch_dir / f"{z:02d}.tif"
        _write_tiff(output_file, dummy)

    print(f"  Channel {channel}: {n_zplanes} z-planes written")

    # QC image
    log_dir = Path(snakemake.log[0]).parent
    log_dir.mkdir(parents=True, exist_ok=True)
    qc_path = log_dir / f"cyc{CYCLE:02d}_CH{channel}_stitched.png"
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt

        fig, ax = plt.subplots(figsize=(4, 4))
        ax.imshow(np.random.rand(64, 64), cmap="gray")
        ax.set_title(f"[MOCK] Cyc{CYCLE} CH{channel}")
        ax.axis("off")
        plt.savefig(str(qc_path), dpi=50)
        plt.close()
        print(f"  QC saved: {qc_path.name}")
    except Exception as e:
        print(f"  QC save failed: {e}")

total_time = time.time() - start_time

print(f"\nTotal time: {total_time:.1f}s")
summary_after("stitch", PROJECT_DIR, total_time, exit_code=0)

# ---------------------------------------------------------------------------
# Write sentinel
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
    f"mock=true\n"
)
print(f"Sentinel written: {sentinel}")
log_footer(0)
