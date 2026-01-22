#!/usr/bin/env python3
"""
Reprocess EDF for affected channels (CD45RO, Podoplan) with quality gate.

This script applies the new EDF parameters that:
1. Exclude edge z-planes (z_start=2, z_end=n_zplanes-1)
2. Auto-detect and exclude corrupted z-slices (>10% zeros or >5% saturation)
"""

import os
import sys
import gc
import numpy as np
from glob import glob
from pathlib import Path
from skimage.io import imread, imsave
from tqdm import tqdm

# Add notebooks to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent / "notebooks"))

from kintsugi.edf import EDFProcessor

# Project paths
PROJECT = Path("/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L")
DECON_DIR = PROJECT / "data/processed/deconvolved"
EDF_DIR = PROJECT / "data/processed/edf"

# EDF parameters (matching notebook settings)
EDF_PARAMS = {
    'radius_x': 2,
    'radius_y': 2,
    'sigma': 10.0,
    'tiles': (3, 3),
    'backend': 'auto',
    'device': 'auto',
    'blend_depth': 0,
    'z_smooth_sigma': 1.0
}

# Affected channels to reprocess
AFFECTED_CHANNELS = [
    {'cycle': 4, 'channel_name': 'Podoplan'},
    {'cycle': 6, 'channel_name': 'CD45RO'},
]


def find_channel_number(cycle: int, channel_name: str) -> int:
    """Find channel number by matching EDF filename to deconvolved directories."""
    edf_path = EDF_DIR / f"cyc{cycle:02d}" / f"{channel_name}.tif"
    if not edf_path.exists():
        raise FileNotFoundError(f"EDF file not found: {edf_path}")

    # Get all EDF files in order and find index
    edf_files = sorted((EDF_DIR / f"cyc{cycle:02d}").glob("*.tif"))
    edf_names = [f.stem for f in edf_files]

    # Get decon channel dirs
    decon_dirs = sorted((DECON_DIR / f"cyc{cycle:02d}").glob("CH*"))

    if len(edf_files) != len(decon_dirs):
        print(f"Warning: EDF count ({len(edf_files)}) != decon count ({len(decon_dirs)})")

    # Match by position (assuming same order)
    try:
        idx = edf_names.index(channel_name)
        return idx + 1  # CH1, CH2, etc.
    except ValueError:
        raise ValueError(f"Channel {channel_name} not found in {edf_names}")


def load_stack_parallel(tiff_files: list) -> np.ndarray:
    """Load z-stack from TIFF files."""
    from concurrent.futures import ThreadPoolExecutor

    def load_single(f):
        return imread(f)

    with ThreadPoolExecutor(max_workers=4) as executor:
        slices = list(executor.map(load_single, tiff_files))

    return np.stack(slices, axis=0)


def process_edf_with_quality_gate(cycle: int, channel: int, channel_name: str):
    """Process EDF with quality gate to exclude corrupted z-slices."""

    print(f"\n{'='*60}")
    print(f"Processing {channel_name} (cyc{cycle:02d}/CH{channel})")
    print(f"{'='*60}")

    # Load deconvolved z-stack
    source_dir = DECON_DIR / f"cyc{cycle:02d}" / f"CH{channel}"
    tiff_files = sorted(source_dir.glob("*.tif"))

    if not tiff_files:
        print(f"ERROR: No files found in {source_dir}")
        return False

    print(f"Loading {len(tiff_files)} z-slices from {source_dir}")
    stack = load_stack_parallel(tiff_files)
    print(f"Stack shape: {stack.shape}, dtype: {stack.dtype}")

    # ======== QUALITY GATE ========
    valid_z_indices = []
    excluded_reasons = []

    print("\nQuality gate analysis:")
    for z_idx in range(stack.shape[0]):
        slice_data = stack[z_idx]
        total_pixels = slice_data.size

        pct_zero = np.sum(slice_data == 0) / total_pixels
        pct_sat = np.sum(slice_data >= 64000) / total_pixels

        status = "✓"
        reason = None

        if pct_zero > 0.10:
            status = "✗"
            reason = f"{pct_zero*100:.1f}% zeros"
            excluded_reasons.append(f"z={z_idx+1}: {reason}")
        elif pct_sat > 0.05:
            status = "✗"
            reason = f"{pct_sat*100:.1f}% saturated"
            excluded_reasons.append(f"z={z_idx+1}: {reason}")
        else:
            valid_z_indices.append(z_idx)

        print(f"  z={z_idx+1:2d}: zeros={pct_zero*100:5.1f}%, sat={pct_sat*100:5.1f}% {status} {reason or ''}")

    print(f"\nValid z-slices: {len(valid_z_indices)}/{stack.shape[0]}")
    if excluded_reasons:
        print("Excluded slices:")
        for reason in excluded_reasons:
            print(f"  - {reason}")

    if len(valid_z_indices) < 3:
        print(f"ERROR: Only {len(valid_z_indices)} valid z-slices, need at least 3")
        return False

    # Filter stack
    stack_filtered = stack[valid_z_indices]
    print(f"Filtered stack shape: {stack_filtered.shape}")

    # Also apply z_start/z_end bounds (exclude first and last of remaining)
    if len(valid_z_indices) > 4:
        # Only exclude edges if we have enough slices
        z_start = 1  # Skip first valid slice
        z_end = len(valid_z_indices) - 1  # Skip last valid slice
        print(f"Applying z_start={z_start+1}, z_end={z_end} (excluding edges of filtered stack)")
    else:
        z_start = 0
        z_end = len(valid_z_indices)
        print(f"Using all {len(valid_z_indices)} valid slices (not enough to exclude edges)")

    # ======== EDF PROCESSING ========
    print("\nRunning EDF...")
    edf_processor = EDFProcessor(backend=EDF_PARAMS['backend'], method='variance')
    print(f"EDF backend: {edf_processor.backend}")

    result = edf_processor.process(
        stack_filtered,
        radius_x=EDF_PARAMS['radius_x'],
        radius_y=EDF_PARAMS['radius_y'],
        sigma=EDF_PARAMS['sigma'],
        z_start=z_start + 1,  # 1-indexed for EDF
        z_end=z_end,
        tiles=EDF_PARAMS['tiles'],
        device=EDF_PARAMS['device'],
        blend_depth=EDF_PARAMS.get('blend_depth', 0),
        z_smooth_sigma=EDF_PARAMS.get('z_smooth_sigma', 0.0)
    )

    # ======== SAVE RESULT ========
    dest_dir = EDF_DIR / f"cyc{cycle:02d}"
    dest_dir.mkdir(parents=True, exist_ok=True)

    # Backup original
    out_path = dest_dir / f"{channel_name}.tif"
    if out_path.exists():
        backup_path = dest_dir / f"{channel_name}_backup.tif"
        print(f"Backing up original to {backup_path.name}")
        out_path.rename(backup_path)

    imsave(out_path, result.astype(np.uint16), check_contrast=False)
    print(f"Saved: {out_path}")

    # ======== VERIFY ========
    pct_zero_result = 100 * np.sum(result == 0) / result.size
    print(f"\nResult statistics:")
    print(f"  Shape: {result.shape}")
    print(f"  Dtype: {result.dtype}")
    print(f"  Min: {result.min()}, Max: {result.max()}")
    print(f"  Percent zeros: {pct_zero_result:.2f}% (should be <5%)")

    if pct_zero_result < 5:
        print(f"  ✓ SUCCESS: Zero percentage is acceptable")
    else:
        print(f"  ⚠ WARNING: Zero percentage still high")

    # Cleanup
    del stack, stack_filtered, result
    gc.collect()

    return True


def main():
    print("=" * 60)
    print("EDF Reprocessing for Affected Channels")
    print("=" * 60)
    print(f"Project: {PROJECT}")
    print(f"Decon dir: {DECON_DIR}")
    print(f"EDF dir: {EDF_DIR}")

    # Initialize EDF processor to check backend
    edf_test = EDFProcessor(backend='auto', method='variance')
    print(f"EDF backend available: {edf_test.backend}")

    results = []
    for channel_info in AFFECTED_CHANNELS:
        cycle = channel_info['cycle']
        channel_name = channel_info['channel_name']

        try:
            # Find channel number
            channel_num = find_channel_number(cycle, channel_name)
            print(f"\nFound {channel_name} at CH{channel_num}")

            # Process
            success = process_edf_with_quality_gate(cycle, channel_num, channel_name)
            results.append((channel_name, success))

        except Exception as e:
            print(f"ERROR processing {channel_name}: {e}")
            import traceback
            traceback.print_exc()
            results.append((channel_name, False))

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    for name, success in results:
        status = "✓ SUCCESS" if success else "✗ FAILED"
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
