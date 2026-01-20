#!/usr/bin/env python3
"""
Check if the 30-pixel stripe pattern exists in raw tiles vs stitched output.
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
from Kio import alphanumeric_key


def analyze_stripe_pattern(img, label):
    """Analyze column stripe pattern in an image."""
    img_f = img.astype(np.float32)

    # Compute column profile
    col_profile = img_f.mean(axis=0)

    # High-pass filter
    hp_sigma = 50
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=hp_sigma)

    # FFT
    col_fft = np.abs(np.fft.fft(col_hp))
    col_fft_half = col_fft[1:len(col_fft)//2]

    # Find peak
    col_peak_idx = np.argmax(col_fft_half) + 1
    col_peak_strength = col_fft_half.max() / col_fft_half.mean() if col_fft_half.mean() > 0 else 0
    col_spacing = len(col_profile) / col_peak_idx if col_peak_idx > 0 else 0

    print(f"  {label}:")
    print(f"    Col HP std: {col_hp.std():.1f}")
    print(f"    Col peak: freq={col_peak_idx}, strength={col_peak_strength:.1f}, spacing={col_spacing:.0f}px")

    return {
        'label': label,
        'col_hp_std': float(col_hp.std()),
        'col_peak_freq': col_peak_idx,
        'col_peak_strength': col_peak_strength,
        'col_spacing': col_spacing,
        'col_profile': col_profile,
        'col_hp': col_hp,
        'col_fft': col_fft_half,
    }


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    stitch_dir = project_dir / "data" / "processed" / "stitched"
    output_dir = stitch_dir / "stripe_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Analyze Z13 (problematic) and Z07 (comparison)
    test_cases = [
        (1, 3, 13, "cyc01_CH3_Z13 (problematic)"),
        (1, 3, 7, "cyc01_CH3_Z07 (comparison)"),
    ]

    for cycle, channel, zplane, label in test_cases:
        print(f"\n{'='*70}")
        print(f"Analyzing {label}")
        print(f"{'='*70}")

        # Load raw tiles
        pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
        raw_files = sorted(
            glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
            key=alphanumeric_key
        )

        if not raw_files:
            print(f"  No raw files found")
            continue

        # Analyze a few raw tiles
        raw_results = []
        for i, f in enumerate(raw_files[:5]):  # First 5 tiles
            img = imread(f)
            r = analyze_stripe_pattern(img, f"Raw tile {i+1}")
            raw_results.append(r)

        # Compute average raw tile metrics
        avg_raw_hp_std = np.mean([r['col_hp_std'] for r in raw_results])
        avg_raw_peak_str = np.mean([r['col_peak_strength'] for r in raw_results])

        # Load stitched image
        stitched_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"
        if stitched_path.exists():
            # Sample a region same size as a tile (for fair comparison)
            stitched_full = imread(stitched_path)

            # Sample from center
            h, w = stitched_full.shape
            tile_h, tile_w = 1440, 1920
            y_start = (h - tile_h) // 2
            x_start = (w - tile_w) // 2
            stitched_crop = stitched_full[y_start:y_start+tile_h, x_start:x_start+tile_w]

            stitched_result = analyze_stripe_pattern(stitched_crop, "Stitched (center crop)")

            # Also analyze full stitched (downsampled)
            ds = 4
            stitched_ds = stitched_full[::ds, ::ds]
            stitched_ds_result = analyze_stripe_pattern(stitched_ds, "Stitched (4x downsampled)")
        else:
            print(f"  Stitched file not found: {stitched_path}")
            continue

        # Create comparison visualization
        fig, axes = plt.subplots(2, 3, figsize=(18, 10))
        fig.suptitle(f"{label}: Raw vs Stitched Stripe Pattern Comparison", fontsize=14)

        # Raw tile example
        raw_tile = imread(raw_files[0])
        axes[0, 0].imshow(raw_tile, cmap='gray')
        axes[0, 0].set_title(f"Raw tile 1 (full)")

        # Raw tile column profile (HP filtered)
        axes[0, 1].plot(raw_results[0]['col_hp'])
        axes[0, 1].set_title(f"Raw tile HP col profile (std={raw_results[0]['col_hp_std']:.1f})")
        axes[0, 1].set_xlabel("Column")

        # Raw tile FFT
        axes[0, 2].semilogy(raw_results[0]['col_fft'][:200])
        axes[0, 2].axvline(raw_results[0]['col_peak_freq'], color='r', linestyle='--')
        axes[0, 2].set_title(f"Raw tile FFT (peak freq={raw_results[0]['col_peak_freq']})")
        axes[0, 2].set_xlabel("Frequency")

        # Stitched crop
        axes[1, 0].imshow(stitched_crop, cmap='gray')
        axes[1, 0].set_title("Stitched (center crop, same size as tile)")

        # Stitched column profile (HP filtered)
        axes[1, 1].plot(stitched_result['col_hp'])
        axes[1, 1].set_title(f"Stitched HP col profile (std={stitched_result['col_hp_std']:.1f})")
        axes[1, 1].set_xlabel("Column")

        # Stitched FFT
        axes[1, 2].semilogy(stitched_result['col_fft'][:200])
        axes[1, 2].axvline(stitched_result['col_peak_freq'], color='r', linestyle='--')
        axes[1, 2].set_title(f"Stitched FFT (peak freq={stitched_result['col_peak_freq']})")
        axes[1, 2].set_xlabel("Frequency")

        plt.tight_layout()

        output_path = output_dir / f"raw_vs_stitched_{label.replace(' ', '_').replace('(', '').replace(')', '')}.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved: {output_path}")
        plt.close(fig)

        # Summary
        print(f"\n  SUMMARY:")
        print(f"    Raw tiles avg HP std: {avg_raw_hp_std:.1f}")
        print(f"    Raw tiles avg peak strength: {avg_raw_peak_str:.1f}")
        print(f"    Stitched crop HP std: {stitched_result['col_hp_std']:.1f}")
        print(f"    Stitched crop peak strength: {stitched_result['col_peak_strength']:.1f}")

        ratio = stitched_result['col_hp_std'] / avg_raw_hp_std if avg_raw_hp_std > 0 else 0
        print(f"    Stitched/Raw HP std ratio: {ratio:.2f}x")

        if ratio > 2:
            print(f"    => STRIPES INTRODUCED DURING STITCHING (ratio > 2)")
        elif ratio > 1.5:
            print(f"    => STRIPES AMPLIFIED DURING STITCHING (ratio 1.5-2)")
        else:
            print(f"    => STRIPES ALREADY IN RAW TILES (ratio < 1.5)")


if __name__ == "__main__":
    main()
