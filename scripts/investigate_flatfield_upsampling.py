#!/usr/bin/env python3
"""
Investigate if flatfield upsampling introduces the 30-pixel stripe pattern.

Theory: BaSiC works at 128x128, then resizes to full tile size (1920x1440).
The resize factor of 15 means any 2-pixel pattern at 128 becomes 30px at full resolution.
"""

import sys
from pathlib import Path
from glob import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.io import imread
from skimage.transform import resize
from scipy import ndimage

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from Kio import load_tiles_parallel, alphanumeric_key


def compute_stripe_metric(img, axis='col'):
    """Compute high-pass profile std (stripe metric)."""
    img_f = img.astype(np.float32)
    if axis == 'col':
        profile = img_f.mean(axis=0)
    else:
        profile = img_f.mean(axis=1)
    hp = profile - ndimage.gaussian_filter1d(profile, sigma=50)
    return float(hp.std())


def analyze_fft(profile):
    """Analyze FFT of a 1D profile."""
    hp = profile - ndimage.gaussian_filter1d(profile, sigma=50)
    fft = np.abs(np.fft.fft(hp))
    fft_half = fft[1:len(fft)//2]
    peak_idx = np.argmax(fft_half) + 1
    spacing = len(profile) / peak_idx if peak_idx > 0 else 0
    return peak_idx, spacing, fft_half


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    output_dir = project_dir / "data" / "processed" / "stitched" / "flatfield_upsampling"
    output_dir.mkdir(parents=True, exist_ok=True)

    test_cases = [
        (1, 3, 13, "Z13_problematic"),
        (1, 3, 7, "Z07_comparison"),
    ]

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
            continue

        print(f"Loading {len(files)} tiles...")
        im_array_init = load_tiles_parallel(files, max_workers=4)
        dtype_max = np.iinfo(im_array_init.dtype).max
        im_array = im_array_init.astype(np.float64) / dtype_max

        # Run BaSiC to get flatfield
        print("Computing BaSiC flatfield...")
        corrector = KCorrectGPU(use_gpu=False, verbose=False, working_size=128)

        # First, let's look at the internal working-size images
        images_small = corrector._resize_images(im_array, 128)
        print(f"  Working size images: {images_small.shape}")

        # Fit to get flatfield
        flatfield_full, darkfield_full = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )

        tile_h, tile_w = im_array.shape[1], im_array.shape[2]
        print(f"  Full-size flatfield: {flatfield_full.shape}")

        # Downsample flatfield back to 128x128 to see working-size pattern
        flatfield_128 = resize(flatfield_full, (128, 128), order=1, preserve_range=True)

        # Analyze flatfield at 128x128
        ff_col_profile_128 = flatfield_128.mean(axis=0)
        peak_128, spacing_128, fft_128 = analyze_fft(ff_col_profile_128)
        print(f"\n  Flatfield at 128x128:")
        print(f"    Col profile range: [{ff_col_profile_128.min():.4f}, {ff_col_profile_128.max():.4f}]")
        print(f"    Peak freq: {peak_128}, spacing: {spacing_128:.1f}px")
        print(f"    Col HP std: {compute_stripe_metric(flatfield_128):.6f}")

        # Analyze flatfield at full resolution
        ff_col_profile_full = flatfield_full.mean(axis=0)
        peak_full, spacing_full, fft_full = analyze_fft(ff_col_profile_full)
        print(f"\n  Flatfield at full resolution ({tile_h}x{tile_w}):")
        print(f"    Col profile range: [{ff_col_profile_full.min():.4f}, {ff_col_profile_full.max():.4f}]")
        print(f"    Peak freq: {peak_full}, spacing: {spacing_full:.1f}px")
        print(f"    Col HP std: {compute_stripe_metric(flatfield_full):.6f}")

        # Apply correction to a single tile and check the division
        print("\n  Testing flatfield division on tile 0:")
        tile = im_array[0]
        flatfield_min = 0.1
        flatfield_safe = np.clip(flatfield_full, flatfield_min, None)

        # Before correction
        tile_metric = compute_stripe_metric(tile)
        print(f"    Raw tile stripe metric: {tile_metric:.4f}")

        # After correction
        corrected = (tile - darkfield_full) / flatfield_safe
        corrected = np.clip(corrected, 0, None)
        corrected_metric = compute_stripe_metric(corrected)
        print(f"    Corrected tile stripe metric: {corrected_metric:.4f}")
        print(f"    Ratio: {corrected_metric/tile_metric:.2f}x")

        # Check the flatfield's effect on stripe patterns
        # The division by flatfield amplifies patterns if flatfield has variations
        ff_reciprocal = 1.0 / flatfield_safe
        ff_recip_metric = compute_stripe_metric(ff_reciprocal)
        print(f"    1/flatfield stripe metric: {ff_recip_metric:.6f}")

        # Test different resize methods for flatfield
        print("\n  Testing resize methods on flatfield:")
        from scipy.ndimage import zoom

        # Create a simple test: sinusoidal pattern at 128x128
        test_pattern_128 = np.sin(np.linspace(0, 8*np.pi, 128))  # 4 cycles
        test_profile_128 = test_pattern_128 * 0.01 + 1.0  # Small variation

        # Resize using different methods
        resize_factor = tile_w / 128

        # Method 1: skimage resize with bilinear
        resized_bilinear = resize(test_profile_128.reshape(1, -1), (1, tile_w), order=1)[0]

        # Method 2: scipy zoom with order=1
        resized_zoom1 = zoom(test_profile_128, resize_factor, order=1)

        # Method 3: scipy zoom with order=3 (cubic)
        resized_zoom3 = zoom(test_profile_128, resize_factor, order=3)

        print(f"    Original 128: 4 cycles")
        print(f"    Bilinear resize to {tile_w}: {len(resized_bilinear)} points")

        # Create visualization
        fig, axes = plt.subplots(3, 4, figsize=(20, 14))
        fig.suptitle(f"Flatfield Upsampling Analysis: {label}", fontsize=14)

        # Row 0: Flatfield images
        im = axes[0, 0].imshow(flatfield_128, cmap='turbo')
        axes[0, 0].set_title("Flatfield 128x128")
        plt.colorbar(im, ax=axes[0, 0])

        im = axes[0, 1].imshow(flatfield_full, cmap='turbo')
        axes[0, 1].set_title(f"Flatfield {tile_h}x{tile_w}")
        plt.colorbar(im, ax=axes[0, 1])

        im = axes[0, 2].imshow(ff_reciprocal, cmap='turbo')
        axes[0, 2].set_title("1/Flatfield")
        plt.colorbar(im, ax=axes[0, 2])

        # High-pass of flatfield
        ff_hp = flatfield_full - ndimage.gaussian_filter(flatfield_full, sigma=50)
        vmin, vmax = np.percentile(ff_hp, [1, 99])
        im = axes[0, 3].imshow(ff_hp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        axes[0, 3].set_title("Flatfield HP filtered")
        plt.colorbar(im, ax=axes[0, 3])

        # Row 1: Column profiles
        axes[1, 0].plot(ff_col_profile_128)
        axes[1, 0].set_title(f"Flatfield col profile (128)\nrange=[{ff_col_profile_128.min():.3f}, {ff_col_profile_128.max():.3f}]")

        axes[1, 1].plot(ff_col_profile_full)
        axes[1, 1].set_title(f"Flatfield col profile (full)")

        # HP filtered profiles
        hp_128 = ff_col_profile_128 - ndimage.gaussian_filter1d(ff_col_profile_128, sigma=10)
        hp_full = ff_col_profile_full - ndimage.gaussian_filter1d(ff_col_profile_full, sigma=50)

        axes[1, 2].plot(hp_128)
        axes[1, 2].set_title(f"Flatfield col HP (128)")

        axes[1, 3].plot(hp_full[:500])  # First 500 pixels
        axes[1, 3].set_title(f"Flatfield col HP (full, first 500px)")

        # Row 2: FFT analysis
        axes[2, 0].semilogy(fft_128[:64])
        axes[2, 0].axvline(peak_128, color='r', linestyle='--')
        axes[2, 0].set_title(f"Flatfield 128 FFT (peak={peak_128})")

        axes[2, 1].semilogy(fft_full[:200])
        axes[2, 1].axvline(peak_full, color='r', linestyle='--')
        axes[2, 1].set_title(f"Flatfield full FFT (peak={peak_full}, ~{spacing_full:.0f}px)")

        # Show corrected tile comparison
        im = axes[2, 2].imshow(tile, cmap='gray')
        axes[2, 2].set_title(f"Raw tile (stripe={tile_metric:.4f})")
        plt.colorbar(im, ax=axes[2, 2])

        im = axes[2, 3].imshow(corrected, cmap='gray')
        axes[2, 3].set_title(f"Corrected tile (stripe={corrected_metric:.4f})")
        plt.colorbar(im, ax=axes[2, 3])

        plt.tight_layout()
        output_path = output_dir / f"flatfield_upsampling_{label}.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\n  Saved: {output_path}")
        plt.close(fig)

        # Check if the problem is in how correction is applied at different signal levels
        print("\n  Signal level analysis:")
        tile_mean = tile.mean()
        ff_mean = flatfield_full.mean()
        print(f"    Tile mean: {tile_mean:.4f}")
        print(f"    Flatfield mean: {ff_mean:.4f}")
        print(f"    Flatfield range: [{flatfield_full.min():.4f}, {flatfield_full.max():.4f}]")

        # The key: at low signal, flatfield variations become more visible relative to signal
        snr_raw = tile.mean() / tile.std()
        snr_corrected = corrected.mean() / corrected.std()
        print(f"    Raw SNR: {snr_raw:.2f}")
        print(f"    Corrected SNR: {snr_corrected:.2f}")


if __name__ == "__main__":
    main()
