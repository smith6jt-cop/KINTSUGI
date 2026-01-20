#!/usr/bin/env python3
"""
Diagnose how sigmoid taper blending creates stripe patterns.

This script investigates why certain z-planes (Z13, Z02) develop 30-pixel
stripe patterns during stitching while others (Z07) don't.

Hypothesis: At low signal levels, the sigmoid taper normalization creates
visible variations at tile boundaries.
"""

import sys
from pathlib import Path
from glob import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import generate_taper, stitch_with_blending
from Kstitch.stitching import stitch_images
from Kio import load_tiles_parallel, alphanumeric_key
import pandas as pd


def compute_stripe_metric(img):
    """Compute high-pass column profile std (stripe metric)."""
    img_f = img.astype(np.float32)
    col_profile = img_f.mean(axis=0)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=50)
    return float(col_hp.std())


def analyze_taper(tile_shape, overlap_fraction, sigma):
    """Analyze properties of the sigmoid taper."""
    taper = generate_taper(tile_shape, overlap_fraction, sigma)

    # Column profile of taper
    taper_col_mean = taper.mean(axis=0)

    return taper, taper_col_mean


def stitch_without_blending(tiles, result_df):
    """Simple tile placement without blending (for comparison)."""
    tile_h, tile_w = tiles.shape[1], tiles.shape[2]

    if "y_pos2" not in result_df.columns:
        result_df = result_df.copy()
        result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
        result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

    out_h = int(result_df["y_pos2"].max() + tile_h)
    out_w = int(result_df["x_pos2"].max() + tile_w)

    stitched = np.zeros((out_h, out_w), dtype=tiles.dtype)

    for i, row in result_df.iterrows():
        y0, x0 = int(row["y_pos2"]), int(row["x_pos2"])
        stitched[y0:y0 + tile_h, x0:x0 + tile_w] = tiles[i]

    return stitched


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    stitch_dir = project_dir / "data" / "processed" / "stitched"
    output_dir = stitch_dir / "blending_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test parameters (from notebook)
    overlap_percentage = 30.0
    BLEND_SIGMA = 10.0

    test_cases = [
        (1, 3, 13, "Z13_problematic"),
        (1, 3, 7, "Z07_comparison"),
    ]

    results = []

    for cycle, channel, zplane, label in test_cases:
        print(f"\n{'='*70}")
        print(f"Processing {label}")
        print(f"{'='*70}")

        # Load tiles
        pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
        files = sorted(
            glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
            key=alphanumeric_key
        )

        if not files:
            print(f"  No files found")
            continue

        print(f"Loading {len(files)} tiles...")
        im_array_init = load_tiles_parallel(files, max_workers=4)
        dtype_max = np.iinfo(im_array_init.dtype).max
        im_array = im_array_init.astype(np.float64) / dtype_max

        # Check raw tile characteristics
        raw_means = [im_array[i].mean() for i in range(len(im_array))]
        raw_stds = [im_array[i].std() for i in range(len(im_array))]
        print(f"  Raw tiles: mean={np.mean(raw_means):.4f}, std={np.mean(raw_stds):.4f}")

        # Apply BaSiC correction
        print("Applying BaSiC correction...")
        corrector = KCorrectGPU(use_gpu=False, verbose=True)
        flatfield, darkfield = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )

        flatfield_min = 0.1
        flatfield_safe = np.clip(flatfield, flatfield_min, None)
        corrected = (im_array - darkfield) / flatfield_safe
        corrected = np.clip(corrected, 0, None)

        # Check corrected tile characteristics
        corr_means = [corrected[i].mean() for i in range(len(corrected))]
        corr_stds = [corrected[i].std() for i in range(len(corrected))]
        print(f"  Corrected tiles: mean={np.mean(corr_means):.4f}, std={np.mean(corr_stds):.4f}")

        # Convert to uint16 for stitching (mimics notebook)
        corrected_uint16 = (corrected * dtype_max).astype(np.uint16)

        # Load existing stitch model
        pkl_path = stitch_dir / f"cyc{cycle:02d}" / "CH1" / "result_df.pkl"
        if pkl_path.exists():
            result_df = pd.read_pickle(pkl_path)
            result_df = result_df.copy()
            result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
            result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()
        else:
            print(f"  No stitch model found at {pkl_path}")
            continue

        # Method 1: Stitch WITH blending (current method)
        print("Stitching WITH sigmoid blending...")
        overlap_fraction = (overlap_percentage / 100.0, overlap_percentage / 100.0)
        stitched_blend = stitch_with_blending(
            corrected_uint16,
            result_df,
            blend=True,
            sigma=BLEND_SIGMA,
            overlap_fraction=overlap_fraction,
            output_dtype=np.uint16
        )

        # Method 2: Stitch WITHOUT blending
        print("Stitching WITHOUT blending...")
        stitched_no_blend = stitch_without_blending(corrected_uint16, result_df)

        # Method 3: Stitch with different sigma values
        stitched_sigma5 = stitch_with_blending(
            corrected_uint16, result_df, blend=True, sigma=5.0,
            overlap_fraction=overlap_fraction, output_dtype=np.uint16
        )
        stitched_sigma20 = stitch_with_blending(
            corrected_uint16, result_df, blend=True, sigma=20.0,
            overlap_fraction=overlap_fraction, output_dtype=np.uint16
        )

        # Analyze stripe metrics for each method
        tile_h, tile_w = corrected_uint16.shape[1], corrected_uint16.shape[2]

        # Sample from center (consistent with previous analysis)
        h, w = stitched_blend.shape
        y_start = (h - tile_h) // 2
        x_start = (w - tile_w) // 2

        crop_blend = stitched_blend[y_start:y_start+tile_h, x_start:x_start+tile_w]
        crop_no_blend = stitched_no_blend[y_start:y_start+tile_h, x_start:x_start+tile_w]
        crop_sigma5 = stitched_sigma5[y_start:y_start+tile_h, x_start:x_start+tile_w]
        crop_sigma20 = stitched_sigma20[y_start:y_start+tile_h, x_start:x_start+tile_w]

        metric_input = compute_stripe_metric(corrected_uint16[0])
        metric_blend = compute_stripe_metric(crop_blend)
        metric_no_blend = compute_stripe_metric(crop_no_blend)
        metric_sigma5 = compute_stripe_metric(crop_sigma5)
        metric_sigma20 = compute_stripe_metric(crop_sigma20)

        print(f"\n  STRIPE METRICS:")
        print(f"    Input tile:        {metric_input:.2f}")
        print(f"    No blending:       {metric_no_blend:.2f} (ratio: {metric_no_blend/metric_input:.2f}x)")
        print(f"    Sigma=5 blending:  {metric_sigma5:.2f} (ratio: {metric_sigma5/metric_input:.2f}x)")
        print(f"    Sigma=10 blending: {metric_blend:.2f} (ratio: {metric_blend/metric_input:.2f}x)")
        print(f"    Sigma=20 blending: {metric_sigma20:.2f} (ratio: {metric_sigma20/metric_input:.2f}x)")

        results.append({
            'label': label,
            'input_mean': np.mean(corr_means),
            'metric_input': metric_input,
            'metric_no_blend': metric_no_blend,
            'metric_sigma5': metric_sigma5,
            'metric_sigma10': metric_blend,
            'metric_sigma20': metric_sigma20,
        })

        # Analyze the weight sum pattern
        print("\n  Analyzing weight distribution...")

        # Recompute stitching manually to get weight_sum
        taper = generate_taper((tile_h, tile_w), overlap_fraction, BLEND_SIGMA)
        out_h = int(result_df["y_pos2"].max() + tile_h)
        out_w = int(result_df["x_pos2"].max() + tile_w)
        weight_sum = np.zeros((out_h, out_w), dtype=np.float64)

        for i, row in result_df.iterrows():
            y0, x0 = int(row["y_pos2"]), int(row["x_pos2"])
            weight_sum[y0:y0 + tile_h, x0:x0 + tile_w] += taper

        # Check weight_sum variations
        weight_crop = weight_sum[y_start:y_start+tile_h, x_start:x_start+tile_w]
        weight_col_profile = weight_crop.mean(axis=0)
        weight_col_hp = weight_col_profile - ndimage.gaussian_filter1d(weight_col_profile, sigma=50)

        print(f"    Weight sum range: [{weight_sum.min():.3f}, {weight_sum.max():.3f}]")
        print(f"    Weight col HP std: {weight_col_hp.std():.6f}")

        # Create visualization
        fig, axes = plt.subplots(3, 4, figsize=(24, 16))
        fig.suptitle(f"Blending Diagnostics: {label}", fontsize=14)

        # Row 0: Images
        im = axes[0, 0].imshow(crop_blend, cmap='gray')
        axes[0, 0].set_title(f"Blended (σ=10)\nStripe: {metric_blend:.1f}")
        plt.colorbar(im, ax=axes[0, 0])

        im = axes[0, 1].imshow(crop_no_blend, cmap='gray')
        axes[0, 1].set_title(f"No blending\nStripe: {metric_no_blend:.1f}")
        plt.colorbar(im, ax=axes[0, 1])

        im = axes[0, 2].imshow(crop_sigma5, cmap='gray')
        axes[0, 2].set_title(f"Blended (σ=5)\nStripe: {metric_sigma5:.1f}")
        plt.colorbar(im, ax=axes[0, 2])

        im = axes[0, 3].imshow(crop_sigma20, cmap='gray')
        axes[0, 3].set_title(f"Blended (σ=20)\nStripe: {metric_sigma20:.1f}")
        plt.colorbar(im, ax=axes[0, 3])

        # Row 1: High-pass filtered
        def show_hp(ax, img, title):
            hp = img.astype(np.float32) - ndimage.gaussian_filter(img.astype(np.float32), sigma=50)
            vmin, vmax = np.percentile(hp, [1, 99])
            im = ax.imshow(hp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
            ax.set_title(title)
            plt.colorbar(im, ax=ax)

        show_hp(axes[1, 0], crop_blend, "HP (σ=10 blend)")
        show_hp(axes[1, 1], crop_no_blend, "HP (no blend)")
        show_hp(axes[1, 2], crop_sigma5, "HP (σ=5 blend)")
        show_hp(axes[1, 3], crop_sigma20, "HP (σ=20 blend)")

        # Row 2: Weight sum and taper analysis
        im = axes[2, 0].imshow(weight_crop, cmap='turbo')
        axes[2, 0].set_title("Weight sum (center crop)")
        plt.colorbar(im, ax=axes[2, 0])

        im = axes[2, 1].imshow(taper, cmap='turbo')
        axes[2, 1].set_title("Single tile taper (σ=10)")
        plt.colorbar(im, ax=axes[2, 1])

        axes[2, 2].plot(weight_col_profile, label='Weight')
        axes[2, 2].set_title("Weight col profile")
        axes[2, 2].set_xlabel("Column")

        axes[2, 3].plot(taper.mean(axis=0), label='Taper')
        axes[2, 3].set_title("Taper col profile")
        axes[2, 3].set_xlabel("Column")

        plt.tight_layout()
        output_path = output_dir / f"blending_diagnosis_{label}.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved: {output_path}")
        plt.close(fig)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for r in results:
        print(f"\n{r['label']}:")
        print(f"  Input mean: {r['input_mean']:.4f}")
        print(f"  No blend vs blend σ=10: {r['metric_no_blend']:.1f} vs {r['metric_sigma10']:.1f}")
        print(f"  Stripe amplification: {r['metric_sigma10']/r['metric_input']:.2f}x")

    # Check if no_blend is better
    print("\n" + "-"*70)
    print("RECOMMENDATION:")
    for r in results:
        if r['metric_no_blend'] < r['metric_sigma10']:
            print(f"  {r['label']}: NO BLENDING is better ({r['metric_no_blend']:.1f} vs {r['metric_sigma10']:.1f})")
        else:
            print(f"  {r['label']}: BLENDING is better ({r['metric_sigma10']:.1f} vs {r['metric_no_blend']:.1f})")


if __name__ == "__main__":
    main()
