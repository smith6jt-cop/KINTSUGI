#!/usr/bin/env python3
"""
Analyze stripe patterns in stitched images.

This script loads the stitched images for problematic z-planes and analyzes
the stripe patterns to determine their source.
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.io import imread
from scipy import ndimage

# Add notebooks to path
sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))


def analyze_stitched_image(stitched_path: Path, output_dir: Path):
    """Analyze stripe patterns in a stitched image."""

    if not stitched_path.exists():
        print(f"  File not found: {stitched_path}")
        return None

    print(f"\nAnalyzing: {stitched_path}")

    # Load image
    img = imread(stitched_path)
    print(f"  Shape: {img.shape}, dtype: {img.dtype}")
    print(f"  Range: [{img.min()}, {img.max()}]")

    # Downsample for faster analysis
    ds_factor = 4
    img_ds = img[::ds_factor, ::ds_factor].astype(np.float32)

    # Compute row and column profiles
    row_profile = img_ds.mean(axis=1)
    col_profile = img_ds.mean(axis=0)

    # High-pass filter the profiles to remove low-frequency structure
    hp_sigma = 50  # pixels (at downsampled resolution)
    row_hp = row_profile - ndimage.gaussian_filter1d(row_profile, sigma=hp_sigma)
    col_hp = col_profile - ndimage.gaussian_filter1d(col_profile, sigma=hp_sigma)

    # Compute FFT of high-pass filtered profiles
    row_fft = np.abs(np.fft.fft(row_hp))
    col_fft = np.abs(np.fft.fft(col_hp))

    # Skip DC component
    row_fft_half = row_fft[1:len(row_fft)//2]
    col_fft_half = col_fft[1:len(col_fft)//2]

    # Find peaks
    row_peak_idx = np.argmax(row_fft_half) + 1
    col_peak_idx = np.argmax(col_fft_half) + 1

    row_peak_strength = row_fft_half.max() / row_fft_half.mean() if row_fft_half.mean() > 0 else 0
    col_peak_strength = col_fft_half.max() / col_fft_half.mean() if col_fft_half.mean() > 0 else 0

    # Convert peak frequency to physical spacing (pixels)
    row_spacing = len(row_profile) / row_peak_idx if row_peak_idx > 0 else 0
    col_spacing = len(col_profile) / col_peak_idx if col_peak_idx > 0 else 0

    # Scale back to full resolution
    row_spacing_full = row_spacing * ds_factor
    col_spacing_full = col_spacing * ds_factor

    print(f"  Row profile:")
    print(f"    Mean: {row_profile.mean():.1f}, Std: {row_profile.std():.1f}")
    print(f"    HP std: {row_hp.std():.1f}")
    print(f"    FFT peak: freq={row_peak_idx}, strength={row_peak_strength:.1f}")
    print(f"    Peak spacing: {row_spacing_full:.0f} pixels")

    print(f"  Col profile:")
    print(f"    Mean: {col_profile.mean():.1f}, Std: {col_profile.std():.1f}")
    print(f"    HP std: {col_hp.std():.1f}")
    print(f"    FFT peak: freq={col_peak_idx}, strength={col_peak_strength:.1f}")
    print(f"    Peak spacing: {col_spacing_full:.0f} pixels")

    # Create visualization
    fig = plt.figure(figsize=(20, 16))

    # Image with enhanced contrast
    ax1 = fig.add_subplot(3, 3, 1)
    vmin, vmax = np.percentile(img_ds, [1, 99])
    ax1.imshow(img_ds, cmap='gray', vmin=vmin, vmax=vmax)
    ax1.set_title(f"Stitched image (4x downsampled)\nRange: [{img.min()}, {img.max()}]")

    # High-pass filtered image (shows stripes better)
    img_hp = img_ds - ndimage.gaussian_filter(img_ds, sigma=50)
    ax2 = fig.add_subplot(3, 3, 2)
    vmin_hp, vmax_hp = np.percentile(img_hp, [1, 99])
    ax2.imshow(img_hp, cmap='RdBu_r', vmin=vmin_hp, vmax=vmax_hp)
    ax2.set_title("High-pass filtered (sigma=50)\n(shows stripe artifacts)")

    # Horizontal derivative (shows vertical stripes)
    dx = np.diff(img_ds, axis=1)
    ax3 = fig.add_subplot(3, 3, 3)
    vmin_dx, vmax_dx = np.percentile(dx, [1, 99])
    ax3.imshow(dx, cmap='RdBu_r', vmin=vmin_dx, vmax=vmax_dx)
    ax3.set_title("Horizontal derivative\n(shows vertical edges/stripes)")

    # Row profile (raw and HP)
    ax4 = fig.add_subplot(3, 3, 4)
    ax4.plot(row_profile, 'b-', alpha=0.5, label='Raw')
    ax4.plot(row_hp + row_profile.mean(), 'r-', label='HP filtered')
    ax4.set_title(f"Row profile (HP std: {row_hp.std():.1f})")
    ax4.set_xlabel("Row")
    ax4.legend()

    # Column profile (raw and HP)
    ax5 = fig.add_subplot(3, 3, 5)
    ax5.plot(col_profile, 'b-', alpha=0.5, label='Raw')
    ax5.plot(col_hp + col_profile.mean(), 'r-', label='HP filtered')
    ax5.set_title(f"Col profile (HP std: {col_hp.std():.1f})")
    ax5.set_xlabel("Column")
    ax5.legend()

    # Row FFT
    ax6 = fig.add_subplot(3, 3, 6)
    freqs = np.arange(1, len(row_fft_half) + 1)
    ax6.semilogy(freqs[:100], row_fft_half[:100])  # First 100 frequencies
    ax6.axvline(row_peak_idx, color='r', linestyle='--', alpha=0.5)
    ax6.set_title(f"Row FFT (peak={row_peak_idx}, spacing={row_spacing_full:.0f}px)")
    ax6.set_xlabel("Frequency")

    # Column FFT
    ax7 = fig.add_subplot(3, 3, 7)
    ax7.semilogy(freqs[:100], col_fft_half[:100])
    ax7.axvline(col_peak_idx, color='r', linestyle='--', alpha=0.5)
    ax7.set_title(f"Col FFT (peak={col_peak_idx}, spacing={col_spacing_full:.0f}px)")
    ax7.set_xlabel("Frequency")

    # Vertical derivative (shows horizontal stripes)
    dy = np.diff(img_ds, axis=0)
    ax8 = fig.add_subplot(3, 3, 8)
    vmin_dy, vmax_dy = np.percentile(dy, [1, 99])
    ax8.imshow(dy, cmap='RdBu_r', vmin=vmin_dy, vmax=vmax_dy)
    ax8.set_title("Vertical derivative\n(shows horizontal edges/stripes)")

    # 2D FFT (power spectrum)
    fft_2d = np.fft.fftshift(np.abs(np.fft.fft2(img_hp)))
    fft_2d_log = np.log1p(fft_2d)
    ax9 = fig.add_subplot(3, 3, 9)
    ax9.imshow(fft_2d_log, cmap='viridis')
    ax9.set_title("2D FFT power spectrum (log)\n(shows periodic patterns)")

    plt.suptitle(stitched_path.name, fontsize=14)
    plt.tight_layout()

    # Save
    output_path = output_dir / f"stripe_analysis_{stitched_path.stem}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")

    plt.close(fig)

    return {
        'file': str(stitched_path),
        'row_hp_std': float(row_hp.std()),
        'col_hp_std': float(col_hp.std()),
        'row_peak_strength': row_peak_strength,
        'col_peak_strength': col_peak_strength,
        'row_spacing': row_spacing_full,
        'col_spacing': col_spacing_full,
    }


def main():
    # Paths
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    stitch_dir = project_dir / "data" / "processed" / "stitched"
    output_dir = stitch_dir / "stripe_analysis"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Target z-planes (user reported stripe artifacts)
    targets = [
        ("cyc01", "CH3", 13),
        ("cyc02", "CH3", 2),
        ("cyc02", "CH3", 7),
    ]

    # Also analyze nearby "clean" z-planes for comparison
    comparison = [
        ("cyc01", "CH3", 7),
        ("cyc02", "CH3", 10),
    ]

    all_targets = targets + comparison

    print("="*70)
    print("STITCHED IMAGE STRIPE ANALYSIS")
    print("="*70)
    print(f"Stitch dir: {stitch_dir}")
    print(f"Output: {output_dir}")

    results = []
    for cycle, channel, zplane in all_targets:
        stitched_path = stitch_dir / cycle / channel / f"{zplane:02d}.tif"
        result = analyze_stitched_image(stitched_path, output_dir)
        if result:
            result['cycle'] = cycle
            result['channel'] = channel
            result['zplane'] = zplane
            result['is_problematic'] = (cycle, channel, zplane) in targets
            results.append(result)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY COMPARISON")
    print("="*70)
    print(f"{'Cycle':<8} {'CH':<4} {'Z':<4} {'Row HP Std':<12} {'Col HP Std':<12} {'Row Pk Str':<12} {'Col Pk Str':<12} {'Row Sp':<10} {'Col Sp':<10}")
    print("-"*90)

    for r in results:
        marker = "* " if r['is_problematic'] else "  "
        print(f"{marker}{r['cycle']:<6} {r['channel']:<4} Z{r['zplane']:02d}  "
              f"{r['row_hp_std']:<12.1f} {r['col_hp_std']:<12.1f} "
              f"{r['row_peak_strength']:<12.1f} {r['col_peak_strength']:<12.1f} "
              f"{r['row_spacing']:<10.0f} {r['col_spacing']:<10.0f}")

    print("\n* = problematic z-plane (user-reported stripes)")
    print("Row/Col Sp = periodic spacing in pixels (if stripe pattern detected)")
    print("\nCheck PNG files in:", output_dir)


if __name__ == "__main__":
    main()
