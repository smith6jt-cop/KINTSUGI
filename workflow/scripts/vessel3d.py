"""
Snakemake wrapper for KINTSUGI 3D vessel segmentation.

Performs Frangi vesselness filtering, binarization, skeletonization, and
graph-based morphometry on deconvolved 3D z-stacks. Targets vascular
markers (CD31, CD34, etc.) to produce binary masks, skeletons, and
per-segment feature tables.

Snakemake guarantees:
  - Deconvolution output exists for this cycle before this script runs

This is a GPU-accelerated step (CuPy Frangi), with CPU fallback.
"""

import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = snakemake.wildcards.cycle
CHANNEL = int(snakemake.params.channel)
MARKER = str(snakemake.params.marker)

# Setup Python path
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_footer, log_header, summary_after, summary_before

# ---------------------------------------------------------------------------
# Configuration from config.yaml
# ---------------------------------------------------------------------------
wf_config = snakemake.config
vessel_cfg = wf_config.get("vessel3d", {})

SIGMAS = vessel_cfg.get("sigmas", [1, 2, 4, 8])
MIN_SIZE = int(vessel_cfg.get("min_size", 500))
PRUNE_MIN_UM = float(vessel_cfg.get("prune_min_um", 5.0))

# Pixel spacing from config (nm -> um)
XY_PIXEL_SIZE = float(wf_config.get("xy_pixel_size", 377.0)) / 1000.0  # nm -> um
Z_STEP_SIZE = float(wf_config.get("z_step_size", 1500.0)) / 1000.0  # nm -> um

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import tifffile

from kintsugi.vessel3d import (
    VesselSpacing,
    export_vessel_results,
    segment_vessels_3d,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DECON_DIR = PROJECT_DIR / "data" / "processed" / "deconvolved" / f"cyc{CYCLE}"
OUTPUT_DIR = PROJECT_DIR / "data" / "processed" / "vessel_3d" / f"cyc{CYCLE}"
SENTINEL = Path(snakemake.output.sentinel)

OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Structured logging
log_header("vessel3d", int(CYCLE), PROJECT_DIR)
summary_before("vessel3d", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"3D Vessel Segmentation - Cycle {CYCLE}, Channel {CHANNEL} ({MARKER})")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Sigmas: {SIGMAS}")
print(f"Min size: {MIN_SIZE} voxels")
print(f"Prune min: {PRUNE_MIN_UM} um")
print(f"Spacing: XY={XY_PIXEL_SIZE:.3f} um, Z={Z_STEP_SIZE:.3f} um")

# ---------------------------------------------------------------------------
# Load deconvolved z-stack for target channel
# ---------------------------------------------------------------------------
start_time = time.time()

# Find z-plane files for the target channel
zplane_files = sorted(DECON_DIR.glob(f"*_CH{CHANNEL}_*.tif"))
if not zplane_files:
    # Try alternative naming
    zplane_files = sorted(DECON_DIR.glob(f"*CH{CHANNEL}*.tif"))

if not zplane_files:
    print(f"ERROR: No z-plane files found for channel {CHANNEL} in {DECON_DIR}")
    print(f"Available files: {list(DECON_DIR.glob('*.tif'))[:10]}")
    log_footer(1)
    sys.exit(1)

print(f"Loading {len(zplane_files)} z-planes...")
planes = []
for f in zplane_files:
    planes.append(tifffile.imread(str(f)))
volume = np.stack(planes, axis=0)
print(f"Volume shape: {volume.shape} ({volume.dtype})")

# ---------------------------------------------------------------------------
# Run segmentation
# ---------------------------------------------------------------------------
spacing = VesselSpacing(xy=XY_PIXEL_SIZE, z=Z_STEP_SIZE)

try:
    result = segment_vessels_3d(
        volume,
        spacing=spacing,
        sigmas=SIGMAS,
        min_size=MIN_SIZE,
        prune_min_length_um=PRUNE_MIN_UM,
        device="auto",
        marker_name=MARKER,
    )

    # Export results
    export_paths = export_vessel_results(
        result,
        output_dir=OUTPUT_DIR,
        marker_name=MARKER,
    )

    elapsed = (time.time() - start_time) / 60

    print(f"\n{'='*60}")
    print(f"Vessel Segmentation Complete")
    print(f"{'='*60}")
    print(f"Vessel volume fraction: {result.binary_mask.sum() / result.binary_mask.size:.4f}")
    if result.features is not None:
        print(f"Segments: {len(result.features)}")
    print(f"Output: {OUTPUT_DIR}")
    print(f"Time: {elapsed:.1f} minutes")

    for name, path in export_paths.items():
        print(f"  {name}: {path}")

    summary_after("vessel3d", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Write sentinel
    SENTINEL.parent.mkdir(parents=True, exist_ok=True)
    SENTINEL.write_text(
        f"stage=vessel3d\n"
        f"cycle={CYCLE}\n"
        f"channel={CHANNEL}\n"
        f"marker={MARKER}\n"
        f"completed={datetime.now().isoformat()}\n"
        f"z_planes={len(zplane_files)}\n"
        f"volume_shape={volume.shape}\n"
        f"segments={len(result.features) if result.features is not None else 0}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {SENTINEL}")
    log_footer(0)

except Exception as e:
    elapsed = (time.time() - start_time) / 60
    print(f"\nERROR: Vessel segmentation failed: {e}")
    import traceback
    traceback.print_exc()
    summary_after("vessel3d", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    sys.exit(1)
