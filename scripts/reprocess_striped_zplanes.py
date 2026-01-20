#!/usr/bin/env python3
"""
Reprocess z-planes with stripe artifacts using current (fixed) code.

This script reprocesses specific z-planes that were identified as having
stripe artifacts when they were originally created with an older version
of the pipeline.

Usage:
    python scripts/reprocess_striped_zplanes.py

The script will:
1. Load raw tiles
2. Apply BaSiC illumination correction
3. Stitch with sigmoid taper blending
4. Save the corrected images (with backup of originals)
"""

import sys
from pathlib import Path
from glob import glob
import shutil
from datetime import datetime

import numpy as np
import pandas as pd
from skimage.io import imread, imsave
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending
from Kio import load_tiles_parallel, alphanumeric_key


# Configuration
PROJECT_DIR = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
OVERLAP_PERCENTAGE = 30.0
BLEND_SIGMA = 10.0
BASIC_FLATFIELD_MIN = 0.1

# Z-planes to reprocess (cycle, channel, zplane, label)
TARGETS = [
    (1, 3, 13, "cyc01_CH3_Z13"),
    (1, 3, 2, "cyc01_CH3_Z02"),
    (1, 3, 7, "cyc01_CH3_Z07"),  # Comparison
    (2, 3, 2, "cyc02_CH3_Z02"),
    (2, 3, 7, "cyc02_CH3_Z07"),
]


def compute_stripe_metric(img):
    """Compute high-pass column profile std (stripe metric)."""
    img_f = img.astype(np.float32)
    col_profile = img_f.mean(axis=0)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=50)
    return float(col_hp.std())


def reprocess_zplane(raw_dir, stitch_dir, cycle, channel, zplane, use_gpu=True, backup=True):
    """
    Reprocess a single z-plane.

    Parameters
    ----------
    raw_dir : Path
        Directory containing raw tiles
    stitch_dir : Path
        Directory for stitched output
    cycle : int
        Cycle number
    channel : int
        Channel number
    zplane : int
        Z-plane number
    use_gpu : bool
        Use GPU for BaSiC correction
    backup : bool
        Create backup of existing file

    Returns
    -------
    dict : Processing results
    """
    label = f"cyc{cycle:02d}_CH{channel}_Z{zplane:02d}"
    print(f"\n{'='*70}")
    print(f"REPROCESSING: {label}")
    print(f"{'='*70}")

    # Paths
    pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
    files = sorted(
        glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
        key=alphanumeric_key
    )

    if not files:
        print(f"  ERROR: No tiles found for {label}")
        return None

    output_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"
    pkl_path = stitch_dir / f"cyc{cycle:02d}" / "CH1" / "result_df.pkl"

    # Check for stitch model
    if not pkl_path.exists():
        print(f"  ERROR: Stitch model not found: {pkl_path}")
        return None

    # Backup existing file
    if backup and output_path.exists():
        backup_path = output_path.with_suffix(".tif.bak")
        if not backup_path.exists():
            print(f"  Backing up: {output_path.name} -> {backup_path.name}")
            shutil.copy2(output_path, backup_path)

    # Get existing stripe metric for comparison
    existing_metric = None
    if output_path.exists():
        existing = imread(output_path)
        existing_metric = compute_stripe_metric(existing[::4, ::4])  # 4x downsampled

    # Load tiles
    print(f"  Loading {len(files)} tiles...")
    start_time = datetime.now()
    im_array_init = load_tiles_parallel(files, max_workers=4)
    dtype_max = np.iinfo(im_array_init.dtype).max
    im_array = im_array_init.astype(np.float64) / dtype_max
    load_time = (datetime.now() - start_time).total_seconds()
    print(f"    Loaded in {load_time:.1f}s")

    # BaSiC correction
    print(f"  Applying BaSiC correction (GPU={use_gpu})...")
    start_time = datetime.now()
    try:
        corrector = KCorrectGPU(use_gpu=use_gpu, verbose=False)
        flatfield, darkfield = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )
    except Exception as e:
        print(f"    GPU failed ({e}), falling back to CPU...")
        corrector = KCorrectGPU(use_gpu=False, verbose=False)
        flatfield, darkfield = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )

    basic_time = (datetime.now() - start_time).total_seconds()
    print(f"    BaSiC completed in {basic_time:.1f}s")
    print(f"    Flatfield range: [{flatfield.min():.3f}, {flatfield.max():.3f}]")

    # Apply correction with flatfield minimum threshold
    flatfield_safe = np.clip(flatfield, BASIC_FLATFIELD_MIN, None)
    clamped_pct = 100.0 * np.sum(flatfield < BASIC_FLATFIELD_MIN) / flatfield.size
    if clamped_pct > 1.0:
        print(f"    WARNING: {clamped_pct:.1f}% of flatfield clamped")

    corrected = (im_array - darkfield) / flatfield_safe
    corrected = np.clip(corrected, 0, 1)
    corrected_uint16 = (corrected * dtype_max).astype(np.uint16)

    # Load stitch model
    result_df = pd.read_pickle(pkl_path)
    result_df = result_df.copy()
    result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
    result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

    # Stitch with blending
    print(f"  Stitching with sigmoid blending (sigma={BLEND_SIGMA})...")
    start_time = datetime.now()
    overlap_fraction = (OVERLAP_PERCENTAGE / 100.0, OVERLAP_PERCENTAGE / 100.0)
    stitched = stitch_with_blending(
        corrected_uint16,
        result_df,
        blend=True,
        sigma=BLEND_SIGMA,
        overlap_fraction=overlap_fraction,
        output_dtype=np.uint16
    )
    stitch_time = (datetime.now() - start_time).total_seconds()
    print(f"    Stitched in {stitch_time:.1f}s, shape: {stitched.shape}")

    # Save
    print(f"  Saving: {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    imsave(output_path, stitched, check_contrast=False)

    # Calculate new stripe metric
    new_metric = compute_stripe_metric(stitched[::4, ::4])

    # Results
    result = {
        'label': label,
        'cycle': cycle,
        'channel': channel,
        'zplane': zplane,
        'existing_metric': existing_metric,
        'new_metric': new_metric,
        'improvement': existing_metric / new_metric if existing_metric and new_metric > 0 else None,
        'load_time': load_time,
        'basic_time': basic_time,
        'stitch_time': stitch_time,
    }

    print(f"\n  Results:")
    if existing_metric:
        print(f"    Existing stripe metric: {existing_metric:.1f}")
    print(f"    New stripe metric:      {new_metric:.1f}")
    if result['improvement']:
        print(f"    Improvement:            {result['improvement']:.2f}x fewer stripes")

    return result


def main():
    raw_dir = PROJECT_DIR / "data" / "raw"
    stitch_dir = PROJECT_DIR / "data" / "processed" / "stitched"

    print("="*70)
    print("KINTSUGI: Reprocess Striped Z-planes")
    print("="*70)
    print(f"Project: {PROJECT_DIR}")
    print(f"Targets: {len(TARGETS)} z-planes")
    print(f"Parameters:")
    print(f"  Overlap: {OVERLAP_PERCENTAGE}%")
    print(f"  Blend sigma: {BLEND_SIGMA}")
    print(f"  Flatfield min: {BASIC_FLATFIELD_MIN}")

    results = []

    for cycle, channel, zplane, label in TARGETS:
        result = reprocess_zplane(
            raw_dir, stitch_dir,
            cycle, channel, zplane,
            use_gpu=True,
            backup=True
        )
        if result:
            results.append(result)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    print(f"{'Label':<20} {'Existing':<12} {'New':<12} {'Improvement':<12}")
    print("-"*60)

    for r in results:
        existing = f"{r['existing_metric']:.0f}" if r['existing_metric'] else "N/A"
        improvement = f"{r['improvement']:.2f}x" if r['improvement'] else "N/A"
        print(f"{r['label']:<20} {existing:<12} {r['new_metric']:<12.0f} {improvement:<12}")

    # Overall
    improved = [r for r in results if r.get('improvement', 0) and r['improvement'] > 1]
    print(f"\nImproved: {len(improved)}/{len(results)} z-planes")

    if improved:
        avg_improvement = np.mean([r['improvement'] for r in improved])
        print(f"Average improvement: {avg_improvement:.2f}x fewer stripes")


if __name__ == "__main__":
    main()
