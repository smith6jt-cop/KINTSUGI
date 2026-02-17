"""
Snakemake wrapper for KINTSUGI REDSEA spillover correction.

Performs post-segmentation lateral signal spillover correction using REDSEA
(Bai et al., 2021). Operates on a multichannel image + cell segmentation mask
to produce corrected per-cell quantification tables.

Snakemake guarantees:
  - Segmentation output exists before this script runs
  - Signal-isolated image is available

This is a CPU-only step (no GPU required).
"""

import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
SAMPLE = snakemake.wildcards.sample

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
redsea_cfg = wf_config.get("redsea", {})

ELEMENT_SHAPE = redsea_cfg.get("element_shape", "star")
ELEMENT_SIZE = int(redsea_cfg.get("element_size", 2))
MARKERS_OF_INTEREST = redsea_cfg.get("markers_of_interest", [])

# ---------------------------------------------------------------------------
# Imports
# ---------------------------------------------------------------------------
import tifffile

from kintsugi.spillover import (
    compute_spillover_metrics,
    estimate_element_size,
    identify_surface_markers,
    plot_spillover_comparison,
    plot_spillover_summary_heatmap,
    run_redsea_correction,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
IMAGE_PATH = Path(snakemake.input.image)
MASK_PATH = Path(snakemake.input.mask)
MARKERS_PATH = Path(snakemake.input.markers)

OUTPUT_CORRECTED = Path(snakemake.output.corrected)
OUTPUT_UNCORRECTED = Path(snakemake.output.uncorrected)
SENTINEL = Path(snakemake.output.sentinel)

OUTPUT_DIR = OUTPUT_CORRECTED.parent
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# Structured logging
log_header("spillover", 0, PROJECT_DIR)
summary_before("spillover", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"REDSEA Spillover Correction - {SAMPLE}")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Image: {IMAGE_PATH}")
print(f"Mask: {MASK_PATH}")
print(f"Element shape: {ELEMENT_SHAPE}, size: {ELEMENT_SIZE}")

# ---------------------------------------------------------------------------
# Load inputs
# ---------------------------------------------------------------------------
start_time = time.time()

# Load marker names
marker_names = pd.read_csv(MARKERS_PATH)["marker_name"].tolist()
print(f"Markers: {len(marker_names)} ({marker_names[:5]}...)")

# Load multichannel image
image = tifffile.imread(str(IMAGE_PATH))
if image.ndim == 3 and image.shape[0] < image.shape[-1]:
    # (C, Y, X) -> (Y, X, C) for redseapy
    image = image.transpose(1, 2, 0)
print(f"Image shape: {image.shape}")

# Load segmentation mask
mask = tifffile.imread(str(MASK_PATH))
print(f"Mask shape: {mask.shape}, {np.unique(mask).max()} cells")

# ---------------------------------------------------------------------------
# Parameter estimation
# ---------------------------------------------------------------------------
if ELEMENT_SIZE == 0:
    # Auto-estimate from cell size
    ELEMENT_SIZE = estimate_element_size(mask)
    print(f"Auto-estimated element_size: {ELEMENT_SIZE}")

if not MARKERS_OF_INTEREST:
    # Auto-identify surface markers
    MARKERS_OF_INTEREST = identify_surface_markers(marker_names)
    print(f"Auto-identified {len(MARKERS_OF_INTEREST)} surface markers")

print(f"Correcting markers: {MARKERS_OF_INTEREST}")

# ---------------------------------------------------------------------------
# Run REDSEA
# ---------------------------------------------------------------------------
try:
    result = run_redsea_correction(
        multichannel_image=image,
        segmentation_mask=mask,
        marker_names=marker_names,
        markers_of_interest=MARKERS_OF_INTEREST if MARKERS_OF_INTEREST else None,
        element_shape=ELEMENT_SHAPE,
        element_size=ELEMENT_SIZE,
        output_dir=OUTPUT_DIR,
    )

    # Save corrected and uncorrected CSVs to the expected output paths
    result.corrected.to_csv(OUTPUT_CORRECTED, index=False)
    result.uncorrected.to_csv(OUTPUT_UNCORRECTED, index=False)

    print(f"\nCorrection complete: {result.n_cells} cells x {result.n_markers} markers")
    print(f"Corrected: {OUTPUT_CORRECTED}")
    print(f"Uncorrected: {OUTPUT_UNCORRECTED}")

    # ---------------------------------------------------------------------------
    # QC visualization
    # ---------------------------------------------------------------------------
    try:
        qc_dir = OUTPUT_DIR / "qc"
        qc_dir.mkdir(exist_ok=True)

        # Default mutually exclusive pairs
        exclusive_pairs = [
            ("CD3", "CD20"), ("CD3e", "CD20"),
            ("CD4", "CD8"), ("CD4", "CD8a"),
        ]
        available_pairs = [
            (a, b) for a, b in exclusive_pairs
            if a in result.corrected.columns and b in result.corrected.columns
        ]

        if available_pairs:
            metrics = compute_spillover_metrics(
                result.corrected, result.uncorrected, available_pairs
            )
            print("\nDouble-positive reduction:")
            for pair, m in metrics["pair_metrics"].items():
                print(
                    f"  {pair[0]}/{pair[1]}: "
                    f"{m['uncorrected_dp_pct']:.1f}% -> {m['corrected_dp_pct']:.1f}%"
                )

            for pair in available_pairs:
                plot_spillover_comparison(
                    result.corrected, result.uncorrected, pair,
                    save_path=qc_dir / f"qc_{pair[0]}_{pair[1]}.pdf",
                )

        plot_spillover_summary_heatmap(
            result.corrected, result.uncorrected, marker_names,
            save_path=qc_dir / "qc_summary_heatmap.pdf",
        )
        print(f"QC plots saved to: {qc_dir}")

    except Exception as e:
        print(f"WARNING: QC visualization failed: {e}")
        # Don't fail the job if QC plots fail

    elapsed = (time.time() - start_time) / 60
    print(f"\n{'='*60}")
    print(f"Spillover Correction Complete")
    print(f"{'='*60}")
    print(f"Time: {elapsed:.1f} minutes")

    summary_after("spillover", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Write sentinel
    SENTINEL.parent.mkdir(parents=True, exist_ok=True)
    SENTINEL.write_text(
        f"stage=spillover_correction\n"
        f"sample={SAMPLE}\n"
        f"completed={datetime.now().isoformat()}\n"
        f"n_cells={result.n_cells}\n"
        f"n_markers={result.n_markers}\n"
        f"element_shape={ELEMENT_SHAPE}\n"
        f"element_size={ELEMENT_SIZE}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {SENTINEL}")
    log_footer(0)

except Exception as e:
    elapsed = (time.time() - start_time) / 60
    print(f"\nERROR: Spillover correction failed: {e}")
    import traceback
    traceback.print_exc()
    summary_after("spillover", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    sys.exit(1)
