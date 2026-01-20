#!/usr/bin/env python3
"""
Test if GPU BaSiC produces different flatfield than CPU.
"""

import sys
from pathlib import Path
from glob import glob

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))
sys.path.insert(0, str(Path(__file__).parent.parent / "src"))

from kintsugi.kcorrect_gpu import KCorrectGPU
from Kio import load_tiles_parallel, alphanumeric_key


def main():
    project_dir = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
    raw_dir = project_dir / "data" / "raw"
    output_dir = project_dir / "data" / "processed" / "stitched" / "gpu_vs_cpu"
    output_dir.mkdir(parents=True, exist_ok=True)

    cycle = 1
    channel = 3
    zplane = 13

    print(f"Loading tiles for cyc{cycle} CH{channel} Z{zplane}...")
    pattern = f"1_*_Z{zplane:03d}_CH{channel}.tif"
    files = sorted(
        glob(str(raw_dir / f"cyc{cycle:03d}" / pattern)),
        key=alphanumeric_key
    )

    im_array_init = load_tiles_parallel(files, max_workers=4)
    dtype_max = np.iinfo(im_array_init.dtype).max
    im_array = im_array_init.astype(np.float64) / dtype_max

    print(f"  Loaded {len(files)} tiles, shape {im_array.shape}")

    # Run CPU BaSiC
    print("\nRunning CPU BaSiC...")
    corrector_cpu = KCorrectGPU(use_gpu=False, verbose=True)
    flatfield_cpu, darkfield_cpu = corrector_cpu.fit(
        im_array,
        if_darkfield=True,
        max_iterations=500,
    )

    # Try GPU BaSiC
    print("\nRunning GPU BaSiC...")
    try:
        corrector_gpu = KCorrectGPU(use_gpu=True, verbose=True, device_id=0)
        flatfield_gpu, darkfield_gpu = corrector_gpu.fit(
            im_array,
            if_darkfield=True,
            max_iterations=500,
        )
        gpu_available = True
    except Exception as e:
        print(f"  GPU not available: {e}")
        gpu_available = False

    # Compare
    print("\n" + "="*70)
    print("COMPARISON")
    print("="*70)

    print(f"\nCPU Flatfield:")
    print(f"  Range: [{flatfield_cpu.min():.4f}, {flatfield_cpu.max():.4f}]")
    print(f"  Mean: {flatfield_cpu.mean():.4f}")
    print(f"  Std: {flatfield_cpu.std():.4f}")

    print(f"\nCPU Darkfield:")
    print(f"  Range: [{darkfield_cpu.min():.6f}, {darkfield_cpu.max():.6f}]")
    print(f"  Mean: {darkfield_cpu.mean():.6f}")

    if gpu_available:
        print(f"\nGPU Flatfield:")
        print(f"  Range: [{flatfield_gpu.min():.4f}, {flatfield_gpu.max():.4f}]")
        print(f"  Mean: {flatfield_gpu.mean():.4f}")
        print(f"  Std: {flatfield_gpu.std():.4f}")

        print(f"\nGPU Darkfield:")
        print(f"  Range: [{darkfield_gpu.min():.6f}, {darkfield_gpu.max():.6f}]")
        print(f"  Mean: {darkfield_gpu.mean():.6f}")

        # Difference
        ff_diff = np.abs(flatfield_gpu - flatfield_cpu)
        df_diff = np.abs(darkfield_gpu - darkfield_cpu)

        print(f"\nFlatfield difference (GPU-CPU):")
        print(f"  Max: {ff_diff.max():.6f}")
        print(f"  Mean: {ff_diff.mean():.6f}")

        print(f"\nDarkfield difference (GPU-CPU):")
        print(f"  Max: {df_diff.max():.6f}")
        print(f"  Mean: {df_diff.mean():.6f}")

        # Visualize
        fig, axes = plt.subplots(2, 4, figsize=(20, 10))
        fig.suptitle(f"GPU vs CPU BaSiC Comparison: cyc{cycle} CH{channel} Z{zplane}")

        im = axes[0, 0].imshow(flatfield_cpu, cmap='turbo')
        axes[0, 0].set_title("CPU Flatfield")
        plt.colorbar(im, ax=axes[0, 0])

        im = axes[0, 1].imshow(flatfield_gpu, cmap='turbo')
        axes[0, 1].set_title("GPU Flatfield")
        plt.colorbar(im, ax=axes[0, 1])

        im = axes[0, 2].imshow(ff_diff, cmap='hot')
        axes[0, 2].set_title(f"Difference (max={ff_diff.max():.4f})")
        plt.colorbar(im, ax=axes[0, 2])

        axes[0, 3].plot(flatfield_cpu.mean(axis=0), label='CPU')
        axes[0, 3].plot(flatfield_gpu.mean(axis=0), label='GPU')
        axes[0, 3].set_title("Flatfield col profile")
        axes[0, 3].legend()

        im = axes[1, 0].imshow(darkfield_cpu, cmap='turbo')
        axes[1, 0].set_title("CPU Darkfield")
        plt.colorbar(im, ax=axes[1, 0])

        im = axes[1, 1].imshow(darkfield_gpu, cmap='turbo')
        axes[1, 1].set_title("GPU Darkfield")
        plt.colorbar(im, ax=axes[1, 1])

        im = axes[1, 2].imshow(df_diff, cmap='hot')
        axes[1, 2].set_title(f"Difference (max={df_diff.max():.6f})")
        plt.colorbar(im, ax=axes[1, 2])

        axes[1, 3].plot(darkfield_cpu.mean(axis=0), label='CPU')
        axes[1, 3].plot(darkfield_gpu.mean(axis=0), label='GPU')
        axes[1, 3].set_title("Darkfield col profile")
        axes[1, 3].legend()

        plt.tight_layout()
        output_path = output_dir / f"gpu_vs_cpu_basic_Z{zplane}.png"
        fig.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nSaved: {output_path}")
        plt.close(fig)


if __name__ == "__main__":
    main()
