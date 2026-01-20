#!/usr/bin/env python3
"""
Reprocess Z13 from scratch and compare with existing stitched file.

This script will:
1. Load raw tiles exactly as the notebook does
2. Apply BaSiC correction
3. Stitch with blending
4. Save the result
5. Compare with existing file pixel-by-pixel
"""

import sys
from pathlib import Path
from glob import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from skimage.io import imread, imsave
from scipy import ndimage
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending
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
    output_dir = stitch_dir / "reprocess_comparison"
    output_dir.mkdir(parents=True, exist_ok=True)

    # Test parameters matching notebook
    cycle = 1
    channel = 3
    zplane = 13
    overlap_percentage = 30.0
    BLEND_SIGMA = 10.0
    BASIC_FLATFIELD_MIN = 0.1

    print(f"="*70)
    print(f"REPROCESSING cyc{cycle} CH{channel} Z{zplane}")
    print(f"="*70)

    # Load tiles exactly as notebook does
    pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
    files = sorted(
        glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
        key=alphanumeric_key
    )

    print(f"\nLoading {len(files)} tiles...")
    print(f"  First file: {files[0]}")
    print(f"  Last file: {files[-1]}")

    im_array_init = load_tiles_parallel(files, max_workers=4)
    dtype_max = np.iinfo(im_array_init.dtype).max
    im_array = im_array_init.astype(np.float64) / dtype_max

    print(f"  Shape: {im_array.shape}")
    print(f"  dtype_max: {dtype_max}")

    # Apply BaSiC correction (using CPU for reproducibility)
    print("\nApplying BaSiC correction (CPU)...")
    corrector = KCorrectGPU(use_gpu=False, verbose=True)
    flatfield, darkfield = corrector.fit(
        im_array,
        if_darkfield=True,
        max_iterations=500,
    )

    # Apply correction
    flatfield_safe = np.clip(flatfield, BASIC_FLATFIELD_MIN, None)
    clamped_pct = 100.0 * np.sum(flatfield < BASIC_FLATFIELD_MIN) / flatfield.size
    if clamped_pct > 1.0:
        print(f"  WARNING: {clamped_pct:.1f}% of flatfield clamped")

    corrected = (im_array - darkfield) / flatfield_safe
    corrected = np.clip(corrected, 0, 1)
    corrected_uint16 = (corrected * dtype_max).astype(np.uint16)

    print(f"  Corrected shape: {corrected_uint16.shape}")
    print(f"  Corrected range: [{corrected_uint16.min()}, {corrected_uint16.max()}]")

    # Load stitch model
    pkl_path = stitch_dir / f"cyc{cycle:02d}" / "CH1" / "result_df.pkl"
    result_df = pd.read_pickle(pkl_path)
    result_df = result_df.copy()
    result_df["y_pos2"] = result_df["y_pos"] - result_df["y_pos"].min()
    result_df["x_pos2"] = result_df["x_pos"] - result_df["x_pos"].min()

    print(f"\nStitch model: {len(result_df)} tiles")
    print(f"  Position range: x=[{result_df['x_pos2'].min():.0f}, {result_df['x_pos2'].max():.0f}], y=[{result_df['y_pos2'].min():.0f}, {result_df['y_pos2'].max():.0f}]")

    # Stitch with blending
    print("\nStitching with sigmoid blending (CPU)...")
    overlap_fraction = (overlap_percentage / 100.0, overlap_percentage / 100.0)
    fresh_stitched = stitch_with_blending(
        corrected_uint16,
        result_df,
        blend=True,
        sigma=BLEND_SIGMA,
        overlap_fraction=overlap_fraction,
        output_dtype=np.uint16
    )

    print(f"  Fresh stitched shape: {fresh_stitched.shape}")
    print(f"  Fresh stitched range: [{fresh_stitched.min()}, {fresh_stitched.max()}]")

    # Save fresh result
    fresh_path = output_dir / f"fresh_cyc{cycle:02d}_CH{channel}_Z{zplane:02d}.tif"
    imsave(fresh_path, fresh_stitched, check_contrast=False)
    print(f"\n  Saved fresh: {fresh_path}")

    # Load existing stitched file
    existing_path = stitch_dir / f"cyc{cycle:02d}" / f"CH{channel}" / f"{zplane:02d}.tif"
    if not existing_path.exists():
        print(f"\nExisting file not found: {existing_path}")
        return

    existing_stitched = imread(existing_path)
    print(f"\nExisting stitched: {existing_path}")
    print(f"  Shape: {existing_stitched.shape}")
    print(f"  Range: [{existing_stitched.min()}, {existing_stitched.max()}]")

    # Compare
    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)

    if fresh_stitched.shape != existing_stitched.shape:
        print(f"  SHAPE MISMATCH: fresh {fresh_stitched.shape} vs existing {existing_stitched.shape}")
        # Use the minimum shape for comparison
        min_h = min(fresh_stitched.shape[0], existing_stitched.shape[0])
        min_w = min(fresh_stitched.shape[1], existing_stitched.shape[1])
        fresh_crop = fresh_stitched[:min_h, :min_w]
        existing_crop = existing_stitched[:min_h, :min_w]
    else:
        fresh_crop = fresh_stitched
        existing_crop = existing_stitched

    # Pixel-by-pixel difference
    diff = fresh_crop.astype(np.float32) - existing_crop.astype(np.float32)
    abs_diff = np.abs(diff)

    print(f"\nPixel differences:")
    print(f"  Mean abs diff: {abs_diff.mean():.2f}")
    print(f"  Max abs diff: {abs_diff.max():.2f}")
    print(f"  Std of diff: {diff.std():.2f}")
    print(f"  Pixels > 100 diff: {(abs_diff > 100).sum()} ({100*(abs_diff > 100).sum()/abs_diff.size:.4f}%)")
    print(f"  Pixels > 1000 diff: {(abs_diff > 1000).sum()} ({100*(abs_diff > 1000).sum()/abs_diff.size:.4f}%)")

    # Stripe metrics
    fresh_metric = compute_stripe_metric(fresh_crop)
    existing_metric = compute_stripe_metric(existing_crop)
    diff_metric = compute_stripe_metric(abs_diff)

    print(f"\nStripe metrics:")
    print(f"  Fresh:    {fresh_metric:.1f}")
    print(f"  Existing: {existing_metric:.1f}")
    print(f"  Diff:     {diff_metric:.1f}")
    print(f"  Ratio (existing/fresh): {existing_metric/fresh_metric:.2f}x")

    # Visualize
    fig, axes = plt.subplots(2, 4, figsize=(24, 12))
    fig.suptitle(f"Fresh vs Existing Comparison: cyc{cycle} CH{channel} Z{zplane}", fontsize=14)

    # Downsampled for visualization
    ds = 4
    fresh_ds = fresh_crop[::ds, ::ds]
    existing_ds = existing_crop[::ds, ::ds]
    diff_ds = diff[::ds, ::ds]

    vmin, vmax = np.percentile(fresh_ds, [1, 99])

    # Row 0: Images
    im = axes[0, 0].imshow(fresh_ds, cmap='gray', vmin=vmin, vmax=vmax)
    axes[0, 0].set_title(f"Fresh (stripe={fresh_metric:.0f})")
    plt.colorbar(im, ax=axes[0, 0])

    im = axes[0, 1].imshow(existing_ds, cmap='gray', vmin=vmin, vmax=vmax)
    axes[0, 1].set_title(f"Existing (stripe={existing_metric:.0f})")
    plt.colorbar(im, ax=axes[0, 1])

    vmin_d, vmax_d = np.percentile(diff_ds, [1, 99])
    im = axes[0, 2].imshow(diff_ds, cmap='RdBu_r', vmin=-max(abs(vmin_d), abs(vmax_d)), vmax=max(abs(vmin_d), abs(vmax_d)))
    axes[0, 2].set_title(f"Difference (mean={abs_diff.mean():.1f})")
    plt.colorbar(im, ax=axes[0, 2])

    # Histogram of differences
    axes[0, 3].hist(diff.flatten(), bins=100, log=True)
    axes[0, 3].set_title("Difference histogram")
    axes[0, 3].set_xlabel("Difference")
    axes[0, 3].axvline(0, color='r', linestyle='--')

    # Row 1: HP filtered
    def show_hp(ax, img, title):
        hp = img.astype(np.float32) - ndimage.gaussian_filter(img.astype(np.float32), sigma=50)
        vmin, vmax = np.percentile(hp, [1, 99])
        im = ax.imshow(hp, cmap='RdBu_r', vmin=vmin, vmax=vmax)
        ax.set_title(title)
        plt.colorbar(im, ax=ax)

    show_hp(axes[1, 0], fresh_ds, "Fresh HP filtered")
    show_hp(axes[1, 1], existing_ds, "Existing HP filtered")

    # Column profiles
    fresh_col = fresh_crop.mean(axis=0)
    existing_col = existing_crop.mean(axis=0)

    axes[1, 2].plot(fresh_col, label='Fresh', alpha=0.7)
    axes[1, 2].plot(existing_col, label='Existing', alpha=0.7)
    axes[1, 2].set_title("Column profiles")
    axes[1, 2].legend()
    axes[1, 2].set_xlabel("Column")

    # HP filtered column profiles
    fresh_col_hp = fresh_col - ndimage.gaussian_filter1d(fresh_col, sigma=50)
    existing_col_hp = existing_col - ndimage.gaussian_filter1d(existing_col, sigma=50)

    axes[1, 3].plot(fresh_col_hp, label='Fresh', alpha=0.7)
    axes[1, 3].plot(existing_col_hp, label='Existing', alpha=0.7)
    axes[1, 3].set_title(f"HP col profiles (fresh std={np.std(fresh_col_hp):.0f}, exist std={np.std(existing_col_hp):.0f})")
    axes[1, 3].legend()
    axes[1, 3].set_xlabel("Column")

    plt.tight_layout()
    output_path = output_dir / f"comparison_cyc{cycle:02d}_CH{channel}_Z{zplane:02d}.png"
    fig.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"\nSaved comparison: {output_path}")
    plt.close(fig)

    # Check if the issue is GPU-specific
    print("\n" + "-"*70)
    print("DIAGNOSIS:")
    if existing_metric > fresh_metric * 2:
        print("  => EXISTING file has 2x+ more stripes than fresh computation")
        print("  => Likely bug in GPU version or earlier pipeline code")
        print("  => RECOMMENDATION: Reprocess with current code")
    elif abs(existing_metric - fresh_metric) / fresh_metric < 0.1:
        print("  => Fresh and existing have similar stripe metrics")
        print("  => Current code reproduces the issue")
    else:
        print(f"  => Some difference: existing/fresh ratio = {existing_metric/fresh_metric:.2f}x")


if __name__ == "__main__":
    main()
