#!/usr/bin/env python3
"""
Compare freshly computed stitched images vs existing ones on disk.

This helps identify if the stripe artifacts are from:
1. The blending algorithm itself
2. Something specific to how the existing files were created
3. A different pipeline version
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy import ndimage


def compute_stripe_metric(img):
    """Compute high-pass column profile std (stripe metric)."""
    img_f = img.astype(np.float32)
    col_profile = img_f.mean(axis=0)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=50)
    return float(col_hp.std())


def analyze_col_fft(img):
    """Analyze column profile FFT to find periodic patterns."""
    img_f = img.astype(np.float32)
    col_profile = img_f.mean(axis=0)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=50)

    fft = np.abs(np.fft.fft(col_hp))
    fft_half = fft[1:len(fft)//2]

    peak_idx = np.argmax(fft_half) + 1
    spacing = len(col_profile) / peak_idx if peak_idx > 0 else 0
    peak_strength = fft_half.max() / fft_half.mean() if fft_half.mean() > 0 else 0

    return {
        'col_hp': col_hp,
        'fft': fft_half,
        'peak_freq': peak_idx,
        'spacing': spacing,
        'peak_strength': peak_strength,
    }


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    stitch_dir = project_dir / "data" / "processed" / "stitched"
    output_dir = stitch_dir / "blending_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    test_cases = [
        (1, 3, 13, "Z13_problematic"),
        (1, 3, 7, "Z07_comparison"),
        (1, 3, 2, "Z02_edge"),
        (1, 3, 5, "Z05_mid"),
    ]

    print("="*70)
    print("EXISTING STITCHED FILE ANALYSIS")
    print("="*70)

    all_results = []

    for cycle, channel, zplane, label in test_cases:
        stitched_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"

        if not stitched_path.exists():
            print(f"\n{label}: File not found")
            continue

        print(f"\n{label} ({stitched_path.name}):")

        img = imread(stitched_path)
        print(f"  Shape: {img.shape}, dtype: {img.dtype}")
        print(f"  Range: [{img.min()}, {img.max()}]")
        print(f"  Mean: {img.mean():.1f}, Std: {img.std():.1f}")

        # Compute stripe metric on full image (4x downsampled)
        ds = 4
        img_ds = img[::ds, ::ds]
        metric_full = compute_stripe_metric(img_ds)

        # Compute on center crop (tile-sized)
        tile_h, tile_w = 1440, 1920
        h, w = img.shape
        y_start = (h - tile_h) // 2
        x_start = (w - tile_w) // 2
        crop = img[y_start:y_start+tile_h, x_start:x_start+tile_w]
        metric_crop = compute_stripe_metric(crop)

        print(f"  Stripe metric (4x DS full): {metric_full:.1f}")
        print(f"  Stripe metric (center crop): {metric_crop:.1f}")

        # FFT analysis
        fft_result = analyze_col_fft(img_ds)
        print(f"  Peak freq: {fft_result['peak_freq']}, spacing: {fft_result['spacing']:.0f}px (at 4x DS)")
        print(f"  Full-res spacing: ~{fft_result['spacing'] * ds:.0f}px")
        print(f"  Peak strength: {fft_result['peak_strength']:.1f}")

        all_results.append({
            'label': label,
            'zplane': zplane,
            'mean': img.mean(),
            'metric_full': metric_full,
            'metric_crop': metric_crop,
            'spacing': fft_result['spacing'] * ds,
            'peak_strength': fft_result['peak_strength'],
        })

    # Summary comparison
    print("\n" + "="*70)
    print("SUMMARY TABLE")
    print("="*70)
    print(f"{'Label':<20} {'Z':<4} {'Mean':<10} {'Stripe(DS)':<12} {'Stripe(crop)':<12} {'Spacing':<10} {'Pk Str':<10}")
    print("-"*80)

    for r in all_results:
        print(f"{r['label']:<20} {r['zplane']:<4} {r['mean']:<10.0f} {r['metric_full']:<12.1f} {r['metric_crop']:<12.1f} {r['spacing']:<10.0f} {r['peak_strength']:<10.1f}")

    # Create comparison visualization
    if len(all_results) >= 2:
        fig, axes = plt.subplots(len(all_results), 4, figsize=(20, 4*len(all_results)))

        for idx, (result, (cycle, channel, zplane, label)) in enumerate(zip(all_results, test_cases)):
            stitched_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"
            img = imread(stitched_path)
            ds = 4
            img_ds = img[::ds, ::ds]

            # Image
            vmin, vmax = np.percentile(img_ds, [1, 99])
            axes[idx, 0].imshow(img_ds, cmap='gray', vmin=vmin, vmax=vmax)
            axes[idx, 0].set_title(f"{label}\nmean={result['mean']:.0f}")

            # HP filtered
            hp = img_ds.astype(np.float32) - ndimage.gaussian_filter(img_ds.astype(np.float32), sigma=50)
            vmin_hp, vmax_hp = np.percentile(hp, [1, 99])
            axes[idx, 1].imshow(hp, cmap='RdBu_r', vmin=vmin_hp, vmax=vmax_hp)
            axes[idx, 1].set_title(f"HP filtered\nstripe={result['metric_full']:.0f}")

            # Column profile
            fft_result = analyze_col_fft(img_ds)
            axes[idx, 2].plot(fft_result['col_hp'])
            axes[idx, 2].set_title(f"Col HP profile\nstd={np.std(fft_result['col_hp']):.1f}")
            axes[idx, 2].set_xlabel("Column")

            # FFT
            axes[idx, 3].semilogy(fft_result['fft'][:200])
            axes[idx, 3].axvline(fft_result['peak_freq'], color='r', linestyle='--')
            axes[idx, 3].set_title(f"Col FFT\npeak={fft_result['peak_freq']}, spacing~{result['spacing']:.0f}px")
            axes[idx, 3].set_xlabel("Frequency")

        plt.tight_layout()
        output_path = output_dir / "existing_stitched_analysis.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved: {output_path}")
        plt.close(fig)

    # Determine if there's a consistent spacing pattern
    print("\n" + "-"*70)
    print("ANALYSIS:")
    spacings = [r['spacing'] for r in all_results]
    if spacings:
        avg_spacing = np.mean(spacings)
        print(f"  Average stripe spacing: {avg_spacing:.0f}px")

        # Check if this matches tile spacing
        tile_width = 1920
        n_cols = 13  # from notebook
        overlap_pct = 30.0
        expected_step = tile_width * (1 - overlap_pct/100)
        print(f"  Expected tile step (30% overlap): {expected_step:.0f}px")
        print(f"  Tile width: {tile_width}px")

        if abs(avg_spacing - expected_step) < 100:
            print(f"  => Stripes appear to correlate with TILE BOUNDARIES")
        elif abs(avg_spacing - tile_width/64) < 50:
            print(f"  => Stripes are ~{tile_width/avg_spacing:.0f} per tile width")


if __name__ == "__main__":
    main()
