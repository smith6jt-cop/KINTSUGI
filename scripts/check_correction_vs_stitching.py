#!/usr/bin/env python3
"""
Determine if stripes are introduced during BaSiC correction or stitching.

This script:
1. Loads raw tiles for Z13 (problematic) and Z07 (comparison)
2. Applies BaSiC correction
3. Compares stripe metrics BEFORE vs AFTER correction
4. Compares individual corrected tiles vs the stitched result
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
from kintsugi.kcorrect_gpu import KCorrectGPU
from Kio import load_tiles_parallel, alphanumeric_key


def compute_stripe_metric(img):
    """Compute high-pass column profile std (stripe metric)."""
    img_f = img.astype(np.float32)
    col_profile = img_f.mean(axis=0)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=50)
    return float(col_hp.std())


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    stitch_dir = project_dir / "data" / "processed" / "stitched"
    output_dir = stitch_dir / "stripe_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    test_cases = [
        (1, 3, 13, "Z13 (problematic)"),
        (1, 3, 7, "Z07 (comparison)"),
    ]

    for cycle, channel, zplane, label in test_cases:
        print(f"\n{'='*70}")
        print(f"Analyzing {label}")
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

        # Compute stripe metrics for raw tiles
        print("Computing stripe metrics for raw tiles...")
        raw_metrics = [compute_stripe_metric(im_array[i]) for i in range(min(10, len(im_array)))]
        raw_avg = np.mean(raw_metrics)
        print(f"  Raw tiles HP std: {raw_avg:.4f} (avg of first 10)")

        # Apply BaSiC correction
        print("Applying BaSiC correction...")
        corrector = KCorrectGPU(use_gpu=False, verbose=True)
        flatfield, darkfield = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )

        # Apply correction to tiles
        flatfield_min = 0.1
        flatfield_safe = np.clip(flatfield, flatfield_min, None)
        corrected = (im_array - darkfield) / flatfield_safe
        corrected = np.clip(corrected, 0, None)

        # Compute stripe metrics for corrected tiles
        print("Computing stripe metrics for corrected tiles...")
        corrected_metrics = [compute_stripe_metric(corrected[i]) for i in range(min(10, len(corrected)))]
        corrected_avg = np.mean(corrected_metrics)
        print(f"  Corrected tiles HP std: {corrected_avg:.4f} (avg of first 10)")

        # Load stitched result
        stitched_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"
        if stitched_path.exists():
            stitched = imread(stitched_path).astype(np.float32) / dtype_max

            # Sample center region (same size as tile)
            h, w = stitched.shape
            tile_h, tile_w = im_array.shape[1], im_array.shape[2]
            y_start = (h - tile_h) // 2
            x_start = (w - tile_w) // 2
            stitched_crop = stitched[y_start:y_start+tile_h, x_start:x_start+tile_w]

            stitched_metric = compute_stripe_metric(stitched_crop)
            print(f"  Stitched (center crop) HP std: {stitched_metric:.4f}")
        else:
            stitched_metric = None
            print(f"  Stitched file not found")

        # Summary
        print(f"\n  SUMMARY for {label}:")
        print(f"    Raw tiles:       {raw_avg:.4f}")
        print(f"    After BaSiC:     {corrected_avg:.4f}  (ratio: {corrected_avg/raw_avg:.2f}x)")
        if stitched_metric:
            print(f"    After stitching: {stitched_metric:.4f}  (ratio: {stitched_metric/corrected_avg:.2f}x)")
            print(f"    Total increase:  {stitched_metric/raw_avg:.2f}x from raw to stitched")

        # Determine where stripes are introduced
        if corrected_avg / raw_avg > 1.5:
            print(f"    => STRIPES INTRODUCED DURING BaSiC CORRECTION")
        elif stitched_metric and stitched_metric / corrected_avg > 1.5:
            print(f"    => STRIPES INTRODUCED DURING STITCHING/BLENDING")
        else:
            print(f"    => No significant stripe amplification")

        # Create diagnostic visualization
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        fig.suptitle(f"Stripe Introduction Analysis: {label}", fontsize=14)

        # Raw tile
        raw_example = im_array[0]
        im = axes[0, 0].imshow(raw_example, cmap='gray')
        axes[0, 0].set_title(f"Raw tile (HP std={raw_metrics[0]:.4f})")
        plt.colorbar(im, ax=axes[0, 0])

        # Flatfield
        im = axes[0, 1].imshow(flatfield, cmap='turbo')
        axes[0, 1].set_title(f"Flatfield (range: [{flatfield.min():.3f}, {flatfield.max():.3f}])")
        plt.colorbar(im, ax=axes[0, 1])

        # Corrected tile
        corrected_example = corrected[0]
        im = axes[0, 2].imshow(corrected_example, cmap='gray')
        axes[0, 2].set_title(f"Corrected tile (HP std={corrected_metrics[0]:.4f})")
        plt.colorbar(im, ax=axes[0, 2])

        # Stitched crop
        if stitched_metric:
            im = axes[0, 3].imshow(stitched_crop, cmap='gray')
            axes[0, 3].set_title(f"Stitched (HP std={stitched_metric:.4f})")
            plt.colorbar(im, ax=axes[0, 3])

        # High-pass filtered versions (bottom row)
        def show_hp(ax, img, title):
            hp = img - ndimage.gaussian_filter(img, sigma=50)
            vmin, vmax = np.percentile(hp, [1, 99])
            im = ax.imshow(hp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
            ax.set_title(title)
            plt.colorbar(im, ax=ax)

        show_hp(axes[1, 0], raw_example, "Raw HP filtered")
        show_hp(axes[1, 1], flatfield, "Flatfield HP filtered")
        show_hp(axes[1, 2], corrected_example, "Corrected HP filtered")
        if stitched_metric:
            show_hp(axes[1, 3], stitched_crop, "Stitched HP filtered")

        plt.tight_layout()

        output_path = output_dir / f"correction_vs_stitching_{label.replace(' ', '_').replace('(', '').replace(')', '')}.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved: {output_path}")
        plt.close(fig)


if __name__ == "__main__":
    main()
