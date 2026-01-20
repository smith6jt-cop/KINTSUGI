#!/usr/bin/env python3
"""
Diagnose flatfield/darkfield models for specific z-planes.

This script loads raw tiles for problematic z-planes, computes the BaSiC
illumination correction models, and visualizes them to identify patterns
that may cause stripe artifacts.

Usage:
    python scripts/diagnose_flatfield.py

Target z-planes:
- cyc01, CH3, Z13
- cyc02, CH3, Z2, Z7
"""

import sys
from pathlib import Path
from glob import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Headless backend
import matplotlib.pyplot as plt
from skimage.io import imread

# Add notebooks to path
sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from Kio import load_tiles_parallel, alphanumeric_key


def analyze_flatfield_for_zplane(
    raw_dir: Path,
    cycle: int,
    channel: int,
    zplane: int,
    output_dir: Path,
):
    """Analyze flatfield/darkfield for a single z-plane."""

    # Find tiles
    pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
    files = sorted(
        glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
        key=alphanumeric_key
    )

    if not files:
        print(f"  No files found for cyc{cycle:03d} CH{channel} Z{zplane:03d}")
        return None

    print(f"\nAnalyzing cyc{cycle:02d} CH{channel} Z{zplane:02d}...")
    print(f"  Loading {len(files)} tiles...")

    # Load tiles
    im_array_init = load_tiles_parallel(files, max_workers=4)
    dtype_max = np.iinfo(im_array_init.dtype).max
    im_array = im_array_init.astype(np.float64) / dtype_max

    print(f"  Tile array shape: {im_array.shape}")
    print(f"  Intensity range: [{im_array.min():.4f}, {im_array.max():.4f}]")

    # Compute BaSiC correction with verbose output (CPU for diagnostics)
    corrector = KCorrectGPU(use_gpu=False, verbose=True)
    flatfield, darkfield = corrector.fit(
        im_array,
        if_darkfield=True,
        max_iterations=500,
    )

    # Analyze flatfield
    print(f"\n  Flatfield stats:")
    print(f"    Shape: {flatfield.shape}")
    print(f"    Range: [{flatfield.min():.4f}, {flatfield.max():.4f}]")
    print(f"    Mean: {flatfield.mean():.4f}")
    print(f"    Std: {flatfield.std():.4f}")

    # Check for negative values
    neg_pct = 100.0 * np.sum(flatfield < 0) / flatfield.size
    if neg_pct > 0:
        print(f"    WARNING: {neg_pct:.1f}% negative values!")

    # Check for stripe patterns in flatfield
    # Analyze row and column profiles
    row_profile = flatfield.mean(axis=1)
    col_profile = flatfield.mean(axis=0)

    row_cv = np.std(row_profile) / np.mean(row_profile) * 100
    col_cv = np.std(col_profile) / np.mean(col_profile) * 100

    print(f"    Row CV: {row_cv:.2f}%")
    print(f"    Col CV: {col_cv:.2f}%")

    # Check for periodicity using FFT
    row_fft = np.abs(np.fft.fft(row_profile - np.mean(row_profile)))
    col_fft = np.abs(np.fft.fft(col_profile - np.mean(col_profile)))

    # Find dominant frequencies (skip DC)
    row_fft_half = row_fft[1:len(row_fft)//2]
    col_fft_half = col_fft[1:len(col_fft)//2]

    row_peak_idx = np.argmax(row_fft_half) + 1
    col_peak_idx = np.argmax(col_fft_half) + 1

    row_peak_strength = row_fft_half.max() / row_fft_half.mean() if row_fft_half.mean() > 0 else 0
    col_peak_strength = col_fft_half.max() / col_fft_half.mean() if col_fft_half.mean() > 0 else 0

    print(f"    Row FFT peak: freq={row_peak_idx}, strength={row_peak_strength:.1f}")
    print(f"    Col FFT peak: freq={col_peak_idx}, strength={col_peak_strength:.1f}")

    # Analyze darkfield
    print(f"\n  Darkfield stats:")
    print(f"    Range: [{darkfield.min():.4f}, {darkfield.max():.4f}]")
    print(f"    Mean: {darkfield.mean():.4f}")

    # Create visualization
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))
    fig.suptitle(f"Flatfield/Darkfield Analysis: cyc{cycle:02d} CH{channel} Z{zplane:02d}", fontsize=14)

    # Flatfield image
    im = axes[0, 0].imshow(flatfield, cmap='turbo')
    axes[0, 0].set_title(f"Flatfield (range: [{flatfield.min():.3f}, {flatfield.max():.3f}])")
    plt.colorbar(im, ax=axes[0, 0])

    # Darkfield image
    im = axes[0, 1].imshow(darkfield, cmap='turbo')
    axes[0, 1].set_title(f"Darkfield (range: [{darkfield.min():.3f}, {darkfield.max():.3f}])")
    plt.colorbar(im, ax=axes[0, 1])

    # Example tile (first one)
    example_tile = im_array[0]
    im = axes[0, 2].imshow(example_tile, cmap='gray')
    axes[0, 2].set_title("Example raw tile")
    plt.colorbar(im, ax=axes[0, 2])

    # Row profile
    axes[1, 0].plot(row_profile)
    axes[1, 0].set_title(f"Flatfield row profile (CV: {row_cv:.1f}%)")
    axes[1, 0].set_xlabel("Row")
    axes[1, 0].set_ylabel("Mean intensity")

    # Column profile
    axes[1, 1].plot(col_profile)
    axes[1, 1].set_title(f"Flatfield col profile (CV: {col_cv:.1f}%)")
    axes[1, 1].set_xlabel("Column")
    axes[1, 1].set_ylabel("Mean intensity")

    # Flatfield histogram
    axes[1, 2].hist(flatfield.flatten(), bins=50, edgecolor='black')
    axes[1, 2].axvline(0, color='r', linestyle='--', label='Zero')
    axes[1, 2].axvline(1, color='g', linestyle='--', label='Unity')
    axes[1, 2].set_title("Flatfield histogram")
    axes[1, 2].legend()

    # Row FFT
    freqs_row = np.arange(1, len(row_fft_half) + 1)
    axes[2, 0].semilogy(freqs_row, row_fft_half)
    axes[2, 0].axvline(row_peak_idx, color='r', linestyle='--', alpha=0.5)
    axes[2, 0].set_title(f"Row FFT (peak strength: {row_peak_strength:.1f})")
    axes[2, 0].set_xlabel("Frequency")
    axes[2, 0].set_ylabel("Magnitude")

    # Column FFT
    freqs_col = np.arange(1, len(col_fft_half) + 1)
    axes[2, 1].semilogy(freqs_col, col_fft_half)
    axes[2, 1].axvline(col_peak_idx, color='r', linestyle='--', alpha=0.5)
    axes[2, 1].set_title(f"Col FFT (peak strength: {col_peak_strength:.1f})")
    axes[2, 1].set_xlabel("Frequency")
    axes[2, 1].set_ylabel("Magnitude")

    # Corrected tile
    flatfield_safe = np.clip(flatfield, 0.1, None)
    corrected = (example_tile - darkfield) / flatfield_safe
    corrected = np.clip(corrected, 0, None)
    im = axes[2, 2].imshow(corrected, cmap='gray')
    axes[2, 2].set_title("Corrected tile (clamped flatfield)")
    plt.colorbar(im, ax=axes[2, 2])

    plt.tight_layout()

    # Save
    output_path = output_dir / f"flatfield_analysis_cyc{cycle:02d}_CH{channel}_Z{zplane:02d}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\n  Saved: {output_path}")

    plt.close(fig)

    return {
        'cycle': cycle,
        'channel': channel,
        'zplane': zplane,
        'flatfield': flatfield,
        'darkfield': darkfield,
        'row_cv': row_cv,
        'col_cv': col_cv,
        'row_peak_strength': row_peak_strength,
        'col_peak_strength': col_peak_strength,
        'neg_pct': neg_pct,
    }


def main():
    # Paths
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    output_dir = project_dir / "data" / "processed" / "stitched" / "flatfield_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Target z-planes (user reported stripe artifacts)
    targets = [
        (1, 3, 13),   # cyc01, CH3, Z13
        (2, 3, 2),    # cyc02, CH3, Z2
        (2, 3, 7),    # cyc02, CH3, Z7
    ]

    # Also analyze nearby "clean" z-planes for comparison
    comparison_planes = [
        (1, 3, 7),    # cyc01, CH3, Z7 (mid-stack, hopefully clean)
        (2, 3, 10),   # cyc02, CH3, Z10 (hopefully clean)
    ]

    all_targets = targets + comparison_planes

    print("="*70)
    print("FLATFIELD/DARKFIELD DIAGNOSTIC ANALYSIS")
    print("="*70)
    print(f"Raw data: {raw_dir}")
    print(f"Output: {output_dir}")
    print(f"\nTargets (problematic): {targets}")
    print(f"Comparison (clean): {comparison_planes}")

    results = []
    for cycle, channel, zplane in all_targets:
        result = analyze_flatfield_for_zplane(
            raw_dir, cycle, channel, zplane, output_dir
        )
        if result:
            results.append(result)

    # Summary comparison
    print("\n" + "="*70)
    print("SUMMARY COMPARISON")
    print("="*70)
    print(f"{'Cycle':<6} {'CH':<4} {'Z':<4} {'Row CV':<10} {'Col CV':<10} {'Row Peak':<12} {'Col Peak':<12} {'Neg%':<8}")
    print("-"*70)

    for r in results:
        is_problematic = (r['cycle'], r['channel'], r['zplane']) in targets
        marker = "* " if is_problematic else "  "
        print(f"{marker}cyc{r['cycle']:02d}  CH{r['channel']}  Z{r['zplane']:02d}  "
              f"{r['row_cv']:<10.2f} {r['col_cv']:<10.2f} "
              f"{r['row_peak_strength']:<12.1f} {r['col_peak_strength']:<12.1f} "
              f"{r['neg_pct']:<8.1f}")

    print("\n* = problematic z-plane (user-reported stripes)")
    print("\nDone! Check output images in:", output_dir)


if __name__ == "__main__":
    main()
