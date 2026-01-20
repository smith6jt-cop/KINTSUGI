#!/usr/bin/env python3
"""
Compare flatfield frequency content between problematic and clean z-planes.

Hypothesis: The 30-pixel stripe pattern comes from high-frequency content
in the flatfield being amplified during upsampling from 128x128 to full resolution.
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
from skimage.transform import resize

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
from kintsugi.kcorrect_gpu import KCorrectGPU
from Kio import load_tiles_parallel, alphanumeric_key


def analyze_flatfield_frequency(flatfield_128, flatfield_full, label, output_dir):
    """Analyze frequency content of flatfield at working size and full resolution."""

    print(f"\nAnalyzing flatfield: {label}")
    print(f"  Working size: {flatfield_128.shape}")
    print(f"  Full size: {flatfield_full.shape}")

    # 2D FFT of working-size flatfield
    fft_128 = np.fft.fftshift(np.abs(np.fft.fft2(flatfield_128 - flatfield_128.mean())))
    fft_128_log = np.log1p(fft_128)

    # 2D FFT of full-size flatfield
    fft_full = np.fft.fftshift(np.abs(np.fft.fft2(flatfield_full - flatfield_full.mean())))
    fft_full_log = np.log1p(fft_full)

    # Column profiles and FFT
    col_profile_128 = flatfield_128.mean(axis=0)
    col_hp_128 = col_profile_128 - ndimage.gaussian_filter1d(col_profile_128, sigma=10)
    col_fft_128 = np.abs(np.fft.fft(col_hp_128))[1:len(col_hp_128)//2]

    col_profile_full = flatfield_full.mean(axis=0)
    col_hp_full = col_profile_full - ndimage.gaussian_filter1d(col_profile_full, sigma=50)
    col_fft_full = np.abs(np.fft.fft(col_hp_full))[1:len(col_hp_full)//2]

    # Peak detection
    peak_128 = np.argmax(col_fft_128) + 1 if len(col_fft_128) > 0 else 0
    peak_full = np.argmax(col_fft_full) + 1 if len(col_fft_full) > 0 else 0

    spacing_128 = 128 / peak_128 if peak_128 > 0 else 0
    spacing_full = flatfield_full.shape[1] / peak_full if peak_full > 0 else 0

    print(f"  128x128 flatfield:")
    print(f"    Col HP std: {col_hp_128.std():.4f}")
    print(f"    Peak freq: {peak_128}, spacing: {spacing_128:.1f}px")

    print(f"  Full-size flatfield:")
    print(f"    Col HP std: {col_hp_full.std():.4f}")
    print(f"    Peak freq: {peak_full}, spacing: {spacing_full:.1f}px")

    # Check high-frequency content in 128x128 flatfield
    # The "danger zone" is frequencies > 32 (Nyquist/2)
    hf_content_128 = col_fft_128[32:].sum() / col_fft_128.sum() if col_fft_128.sum() > 0 else 0
    print(f"  High-freq content (>32 cycles): {100*hf_content_128:.1f}%")

    # Create visualization
    fig, axes = plt.subplots(3, 3, figsize=(18, 15))
    fig.suptitle(f"Flatfield Frequency Analysis: {label}", fontsize=14)

    # 128x128 flatfield
    im = axes[0, 0].imshow(flatfield_128, cmap='turbo')
    axes[0, 0].set_title("Flatfield (128x128 working size)")
    plt.colorbar(im, ax=axes[0, 0])

    # 128x128 FFT
    im = axes[0, 1].imshow(fft_128_log, cmap='viridis')
    axes[0, 1].set_title("2D FFT (128x128, log scale)")
    plt.colorbar(im, ax=axes[0, 1])

    # 128x128 col profile
    axes[0, 2].plot(col_profile_128, 'b-', alpha=0.5, label='Raw')
    axes[0, 2].plot(col_hp_128 + col_profile_128.mean(), 'r-', label='HP')
    axes[0, 2].set_title(f"128 col profile (HP std={col_hp_128.std():.4f})")
    axes[0, 2].legend()

    # Full-size flatfield
    im = axes[1, 0].imshow(flatfield_full, cmap='turbo')
    axes[1, 0].set_title(f"Flatfield (full size {flatfield_full.shape})")
    plt.colorbar(im, ax=axes[1, 0])

    # Full-size FFT (center crop for visibility)
    center = fft_full_log.shape[0] // 2
    crop = 200
    fft_crop = fft_full_log[center-crop:center+crop, center-crop:center+crop]
    im = axes[1, 1].imshow(fft_crop, cmap='viridis')
    axes[1, 1].set_title("2D FFT (full, center crop, log)")
    plt.colorbar(im, ax=axes[1, 1])

    # Full-size col profile
    # Show a small portion for detail
    mid = len(col_profile_full) // 2
    window = 500
    axes[1, 2].plot(col_profile_full[mid-window:mid+window], 'b-', alpha=0.5)
    axes[1, 2].set_title(f"Full col profile (center 1000px)")

    # 128x128 col FFT
    axes[2, 0].semilogy(col_fft_128)
    axes[2, 0].axvline(peak_128-1, color='r', linestyle='--')
    axes[2, 0].axvline(32, color='orange', linestyle=':', label='Nyquist/2')
    axes[2, 0].set_title(f"128 col FFT (peak={peak_128})")
    axes[2, 0].legend()

    # Full-size col FFT (first 200 frequencies)
    axes[2, 1].semilogy(col_fft_full[:200])
    axes[2, 1].axvline(peak_full-1, color='r', linestyle='--')
    axes[2, 1].axvline(64, color='orange', linestyle=':', label='64 (30px)')
    axes[2, 1].set_title(f"Full col FFT (peak={peak_full})")
    axes[2, 1].legend()

    # Difference from mean (shows local variations)
    flatfield_diff = flatfield_full - flatfield_full.mean()
    im = axes[2, 2].imshow(flatfield_diff, cmap='RdBu_r',
                          vmin=-flatfield_diff.std()*3, vmax=flatfield_diff.std()*3)
    axes[2, 2].set_title("Flatfield deviation from mean")
    plt.colorbar(im, ax=axes[2, 2])

    plt.tight_layout()

    output_path = output_dir / f"flatfield_frequency_{label.replace(' ', '_')}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved: {output_path}")
    plt.close(fig)

    return {
        'label': label,
        'col_hp_std_128': col_hp_128.std(),
        'col_hp_std_full': col_hp_full.std(),
        'peak_freq_128': peak_128,
        'peak_freq_full': peak_full,
        'hf_content_128': hf_content_128,
    }


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    output_dir = project_dir / "data" / "processed" / "stitched" / "flatfield_diagnostics"
    output_dir.mkdir(parents=True, exist_ok=True)

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
            continue

        print(f"Loading {len(files)} tiles...")
        im_array_init = load_tiles_parallel(files, max_workers=4)
        dtype_max = np.iinfo(im_array_init.dtype).max
        im_array = im_array_init.astype(np.float64) / dtype_max

        # Compute BaSiC correction (CPU)
        corrector = KCorrectGPU(use_gpu=False, verbose=True, working_size=128)

        # Get the working-size flatfield by intercepting the resize
        # We need to access the internal _resize_images method
        images_small = corrector._resize_images(im_array, 128)

        # Run the full fit to get flatfield
        flatfield_full, darkfield_full = corrector.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )

        # The working-size flatfield is inside the fit() method
        # Let's recompute it to get the 128x128 version
        images_small_t = np.transpose(images_small, (1, 2, 0))
        mean_image = np.mean(images_small_t, axis=2)
        mean_div = mean_image / (np.mean(mean_image) + 1e-10)

        # This is what the flatfield looks like at 128x128 (before resize)
        # Note: This is an approximation since we can't easily extract
        # the intermediate flatfield from the algorithm
        flatfield_128_approx = resize(
            flatfield_full,
            (128, 128),
            order=1,
            mode='symmetric',
            preserve_range=True,
        )

        result = analyze_flatfield_frequency(
            flatfield_128_approx,
            flatfield_full,
            label,
            output_dir
        )
        results.append(result)

    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    for r in results:
        print(f"\n{r['label']}:")
        print(f"  128px HP std: {r['col_hp_std_128']:.4f}")
        print(f"  Full HP std: {r['col_hp_std_full']:.4f}")
        print(f"  128 peak freq: {r['peak_freq_128']}")
        print(f"  Full peak freq: {r['peak_freq_full']}")
        print(f"  HF content: {100*r['hf_content_128']:.1f}%")


if __name__ == "__main__":
    main()
