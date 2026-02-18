"""
KINTSUGI Processing QC and Statistics Module

This module contains GPU-accelerated functions for computing image statistics
and generating QC plots throughout the KINTSUGI processing pipeline.

Location: notebooks/Kprocess.py

Extracted from 2_Cycle_Processing.ipynb to improve code reusability and
maintainability while preserving processing transparency.

Functions:
---------
Statistics Computation:
    - compute_zplane_stats_gpu: Compute stats for raw tile z-plane
    - compute_stitched_stats_gpu: Compute stats for stitched image
    - compute_decon_stats_gpu: Compute stats for deconvolved image
    - compute_edf_stats_gpu: Compute stats for EDF projection

Parallel Collection:
    - collect_raw_stats_parallel: Collect raw data statistics
    - collect_stitched_stats_parallel: Collect stitched statistics
    - collect_decon_stats_parallel: Collect deconvolved statistics
    - collect_edf_stats_parallel: Collect EDF statistics

Visualization:
    - plot_summary_heatmaps: Create SNR/CV heatmaps by cycle/channel
    - plot_zplane_profiles: Plot statistics across z-planes

Utilities:
    - find_edf_directory: Locate EDF output directory
    - discover_edf_files: Build mapping of EDF files
"""

import os
import pickle
from concurrent.futures import ThreadPoolExecutor, as_completed
from glob import glob
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from skimage.io import imread

from Kio import alphanumeric_key, cleanup_gpu_memory, load_tiles_parallel


# =============================================================================
# CONFIGURATION CLASS
# =============================================================================

class QCConfig:
    """Configuration for QC statistics collection.

    Provides default values that can be overridden when calling functions.
    """
    def __init__(
        self,
        gpu_device_ids: List[int] = None,
        io_workers: int = 4,
        zplanes_per_gpu: int = 2,
        channel_name_dict: Dict[int, List[str]] = None,
        qc_output_dir: Path = None,
    ):
        self.gpu_device_ids = gpu_device_ids or [0]
        self.io_workers = io_workers
        self.zplanes_per_gpu = zplanes_per_gpu
        self.channel_name_dict = channel_name_dict or {}
        self.qc_output_dir = qc_output_dir


# =============================================================================
# RAW DATA STATISTICS
# =============================================================================

def compute_zplane_stats_gpu(
    image_dir: Union[str, Path],
    cycle: int,
    channel: int,
    zplane: int,
    device_id: int = 0,
    io_workers: int = 4,
) -> Optional[Dict[str, float]]:
    """Compute statistics for a single raw z-plane using GPU acceleration.

    Parameters
    ----------
    image_dir : Path
        Directory containing raw images (with cyc### subdirectories)
    cycle : int
        Cycle number
    channel : int
        Channel number (1-indexed)
    zplane : int
        Z-plane number (1-indexed)
    device_id : int
        GPU device ID to use
    io_workers : int
        Number of parallel I/O workers for tile loading

    Returns
    -------
    dict or None
        Dictionary with 'overall_mean', 'overall_std', 'snr', 'row_cv', 'col_cv'
        Returns None if no files found, or dict with 'error' key on failure
    """
    import cupy as cp

    filename_pattern = f'1_*_Z{str(zplane).zfill(3)}_CH{str(channel)}.tif'
    files = sorted(
        glob(os.path.join(image_dir, f'cyc{str(cycle).zfill(3)}', filename_pattern)),
        key=alphanumeric_key
    )

    if not files:
        return None

    try:
        with cp.cuda.Device(device_id):
            im_array = load_tiles_parallel(files, max_workers=io_workers)
            im_gpu = cp.asarray(im_array)

            tile_mean = cp.mean(im_gpu, axis=0)
            row_profile = cp.mean(tile_mean, axis=1)
            col_profile = cp.mean(tile_mean, axis=0)

            row_mean = cp.mean(row_profile)
            col_mean = cp.mean(col_profile)
            row_cv = float(cp.std(row_profile) / row_mean * 100) if float(row_mean) > 0 else 0.0
            col_cv = float(cp.std(col_profile) / col_mean * 100) if float(col_mean) > 0 else 0.0

            overall_mean = float(cp.mean(im_gpu))
            overall_std = float(cp.std(im_gpu))
            snr = overall_mean / overall_std if overall_std > 0 else 0.0

            del im_gpu, tile_mean, row_profile, col_profile
            cp.get_default_memory_pool().free_all_blocks()

        return {
            'overall_mean': overall_mean,
            'overall_std': overall_std,
            'snr': snr,
            'row_cv': row_cv,
            'col_cv': col_cv
        }
    except Exception as e:
        return {'error': str(e)}


def _process_single_raw_stat(args: Tuple) -> Optional[Dict]:
    """Worker function for parallel raw stat collection."""
    image_dir, cycle, channel, zplane, device_id, channel_name_dict, io_workers = args
    stats = compute_zplane_stats_gpu(image_dir, cycle, channel, zplane, device_id, io_workers)
    if stats and 'error' not in stats:
        ch_name = channel_name_dict.get(cycle, ['?']*4)[channel-1] if channel_name_dict else f'CH{channel}'
        return {
            'cycle': cycle,
            'channel': channel,
            'zplane': zplane,
            'channel_name': ch_name,
            **stats
        }
    return None


def collect_raw_stats_parallel(
    image_dir: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Collect raw data statistics using GPU-accelerated parallel processing.

    Parameters
    ----------
    image_dir : Path
        Directory containing raw images
    start_cycle, end_cycle : int
        Cycle range (inclusive)
    start_channel, end_channel : int
        Channel range (inclusive)
    n_zplanes : int
        Number of z-planes
    gpu_device_ids : list
        List of GPU device IDs to use (default: [0])
    io_workers : int
        Number of parallel I/O workers
    zplanes_per_gpu : int
        Number of concurrent z-planes per GPU
    channel_name_dict : dict
        Mapping of cycle -> list of channel names

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: cycle, channel, zplane, channel_name,
        overall_mean, overall_std, snr, row_cv, col_cv
    """
    gpu_ids = gpu_device_ids or [0]
    num_gpus = len(gpu_ids)
    channel_names = channel_name_dict or {}

    # Build work items
    work_items = []
    idx = 0
    for cycle in range(start_cycle, end_cycle + 1):
        for channel in range(start_channel, end_channel + 1):
            for zplane in range(1, n_zplanes + 1):
                device_id = gpu_ids[idx % num_gpus]
                work_items.append((image_dir, cycle, channel, zplane, device_id, channel_names, io_workers))
                idx += 1

    total = len(work_items)
    all_stats = []

    # Use ThreadPoolExecutor with progress tracking
    max_workers = min(io_workers, num_gpus * zplanes_per_gpu)
    print(f"  Using {max_workers} parallel workers across {num_gpus} GPU(s)")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_single_raw_stat, item): item for item in work_items}

        completed = 0
        for future in as_completed(futures):
            result = future.result()
            if result:
                all_stats.append(result)
            completed += 1
            if completed % 20 == 0 or completed == total:
                print(f"\r  Progress: {completed}/{total} ({100*completed/total:.0f}%)", end='', flush=True)

    print()  # Newline after progress
    cleanup_gpu_memory()

    df = pd.DataFrame(all_stats)
    if len(df) > 0:
        df = df.sort_values(['cycle', 'channel', 'zplane']).reset_index(drop=True)
    return df


# =============================================================================
# STITCHED IMAGE STATISTICS
# =============================================================================

def compute_stitched_stats_gpu(
    stitch_dir: Union[str, Path],
    cycle: int,
    channel: int,
    zplane: int,
    device_id: int = 0,
) -> Optional[Dict[str, float]]:
    """Compute statistics for a single stitched image using GPU.

    Parameters
    ----------
    stitch_dir : Path
        Directory containing stitched images
    cycle : int
        Cycle number
    channel : int
        Channel number (1-indexed)
    zplane : int
        Z-plane number (1-indexed)
    device_id : int
        GPU device ID to use

    Returns
    -------
    dict or None
        Statistics dictionary or None if file not found
    """
    import cupy as cp

    stitched_path = os.path.join(stitch_dir, f"cyc{cycle:02d}", f"CH{channel}", f"{zplane:02d}.tif")
    if not os.path.exists(stitched_path):
        return None

    try:
        with cp.cuda.Device(device_id):
            img = imread(stitched_path)
            img_gpu = cp.asarray(img)

            overall_mean = float(cp.mean(img_gpu))
            overall_std = float(cp.std(img_gpu))
            snr = overall_mean / overall_std if overall_std > 0 else 0.0

            # Downsample for profile computation
            step = 4 if img.shape[0] > 2000 else 1
            ds_gpu = img_gpu[::step, ::step]
            row_profile = cp.mean(ds_gpu, axis=1)
            col_profile = cp.mean(ds_gpu, axis=0)

            row_mean = cp.mean(row_profile)
            col_mean = cp.mean(col_profile)
            row_cv = float(cp.std(row_profile) / row_mean * 100) if float(row_mean) > 0 else 0.0
            col_cv = float(cp.std(col_profile) / col_mean * 100) if float(col_mean) > 0 else 0.0

            del img_gpu, ds_gpu, row_profile, col_profile
            cp.get_default_memory_pool().free_all_blocks()

        return {
            'overall_mean': overall_mean,
            'overall_std': overall_std,
            'snr': snr,
            'row_cv': row_cv,
            'col_cv': col_cv
        }
    except Exception as e:
        return {'error': str(e)}


def _process_single_stitched_stat(args: Tuple) -> Optional[Dict]:
    """Worker function for parallel stitched stat collection."""
    stitch_dir, cycle, channel, zplane, device_id, channel_name_dict = args
    stats = compute_stitched_stats_gpu(stitch_dir, cycle, channel, zplane, device_id)
    if stats and 'error' not in stats:
        ch_name = channel_name_dict.get(cycle, ['?']*4)[channel-1] if channel_name_dict else f'CH{channel}'
        return {
            'cycle': cycle,
            'channel': channel,
            'zplane': zplane,
            'channel_name': ch_name,
            **stats
        }
    return None


def collect_stitched_stats_parallel(
    stitch_dir: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Collect stitched statistics using GPU-accelerated parallel processing."""
    gpu_ids = gpu_device_ids or [0]
    num_gpus = len(gpu_ids)
    channel_names = channel_name_dict or {}

    work_items = []
    idx = 0
    for cycle in range(start_cycle, end_cycle + 1):
        for channel in range(start_channel, end_channel + 1):
            for zplane in range(1, n_zplanes + 1):
                device_id = gpu_ids[idx % num_gpus]
                work_items.append((stitch_dir, cycle, channel, zplane, device_id, channel_names))
                idx += 1

    total = len(work_items)
    all_stats = []

    max_workers = min(io_workers, num_gpus * zplanes_per_gpu)
    print(f"  Using {max_workers} parallel workers across {num_gpus} GPU(s)")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_single_stitched_stat, item): item for item in work_items}

        completed = 0
        for future in as_completed(futures):
            result = future.result()
            if result:
                all_stats.append(result)
            completed += 1
            if completed % 20 == 0 or completed == total:
                print(f"\r  Progress: {completed}/{total} ({100*completed/total:.0f}%)", end='', flush=True)

    print()
    cleanup_gpu_memory()

    df = pd.DataFrame(all_stats)
    if len(df) > 0:
        df = df.sort_values(['cycle', 'channel', 'zplane']).reset_index(drop=True)
    return df


# =============================================================================
# DECONVOLVED IMAGE STATISTICS
# =============================================================================

def compute_decon_stats_gpu(
    decon_dir: Union[str, Path],
    cycle: int,
    channel: int,
    zplane: int,
    device_id: int = 0,
) -> Optional[Dict[str, float]]:
    """Compute statistics for a single deconvolved image using GPU.

    Parameters
    ----------
    decon_dir : Path
        Directory containing deconvolved images
    cycle : int
        Cycle number
    channel : int
        Channel number (1-indexed)
    zplane : int
        Z-plane number (1-indexed)
    device_id : int
        GPU device ID to use

    Returns
    -------
    dict or None
        Statistics dictionary or None if file not found
    """
    import cupy as cp

    decon_path = os.path.join(decon_dir, f"cyc{cycle:02d}", f"CH{channel}", f"deconv_{zplane:03d}.tif")
    if not os.path.exists(decon_path):
        return None

    try:
        with cp.cuda.Device(device_id):
            img = imread(decon_path)
            img_gpu = cp.asarray(img)

            overall_mean = float(cp.mean(img_gpu))
            overall_std = float(cp.std(img_gpu))
            snr = overall_mean / overall_std if overall_std > 0 else 0.0

            # Downsample for profile computation
            step = 4 if img.shape[0] > 2000 else 1
            ds_gpu = img_gpu[::step, ::step]
            row_profile = cp.mean(ds_gpu, axis=1)
            col_profile = cp.mean(ds_gpu, axis=0)

            row_mean = cp.mean(row_profile)
            col_mean = cp.mean(col_profile)
            row_cv = float(cp.std(row_profile) / row_mean * 100) if float(row_mean) > 0 else 0.0
            col_cv = float(cp.std(col_profile) / col_mean * 100) if float(col_mean) > 0 else 0.0

            del img_gpu, ds_gpu, row_profile, col_profile
            cp.get_default_memory_pool().free_all_blocks()

        return {
            'overall_mean': overall_mean,
            'overall_std': overall_std,
            'snr': snr,
            'row_cv': row_cv,
            'col_cv': col_cv
        }
    except Exception as e:
        return {'error': str(e)}


def _process_single_decon_stat(args: Tuple) -> Optional[Dict]:
    """Worker function for parallel decon stat collection."""
    decon_dir, cycle, channel, zplane, device_id, channel_name_dict = args
    stats = compute_decon_stats_gpu(decon_dir, cycle, channel, zplane, device_id)
    if stats and 'error' not in stats:
        ch_name = channel_name_dict.get(cycle, ['?']*4)[channel-1] if channel_name_dict else f'CH{channel}'
        return {
            'cycle': cycle,
            'channel': channel,
            'zplane': zplane,
            'channel_name': ch_name,
            **stats
        }
    return None


def collect_decon_stats_parallel(
    decon_dir: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Collect deconvolved statistics using GPU-accelerated parallel processing."""
    gpu_ids = gpu_device_ids or [0]
    num_gpus = len(gpu_ids)
    channel_names = channel_name_dict or {}

    work_items = []
    idx = 0
    for cycle in range(start_cycle, end_cycle + 1):
        for channel in range(start_channel, end_channel + 1):
            for zplane in range(1, n_zplanes + 1):
                device_id = gpu_ids[idx % num_gpus]
                work_items.append((decon_dir, cycle, channel, zplane, device_id, channel_names))
                idx += 1

    total = len(work_items)
    all_stats = []

    max_workers = min(io_workers, num_gpus * zplanes_per_gpu)
    print(f"  Using {max_workers} parallel workers across {num_gpus} GPU(s)")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_single_decon_stat, item): item for item in work_items}

        completed = 0
        for future in as_completed(futures):
            result = future.result()
            if result:
                all_stats.append(result)
            completed += 1
            if completed % 20 == 0 or completed == total:
                print(f"\r  Progress: {completed}/{total} ({100*completed/total:.0f}%)", end='', flush=True)

    print()
    cleanup_gpu_memory()

    df = pd.DataFrame(all_stats)
    if len(df) > 0:
        df = df.sort_values(['cycle', 'channel', 'zplane']).reset_index(drop=True)
    return df


# =============================================================================
# EDF PROJECTION STATISTICS
# =============================================================================

def find_edf_directory(
    project_dir: Union[str, Path],
    default_edf_dir: Union[str, Path] = None,
) -> Path:
    """Find the actual EDF directory (handles different project structures).

    Parameters
    ----------
    project_dir : Path
        Project root directory
    default_edf_dir : Path, optional
        Default EDF directory to try first

    Returns
    -------
    Path
        Path to EDF directory (may not exist if not found)
    """
    project_dir = Path(project_dir)

    # Build list of candidates
    candidates = []
    if default_edf_dir:
        candidates.append(Path(default_edf_dir))
    candidates.extend([
        project_dir / 'data' / 'processed' / 'edf',
        project_dir / 'processed' / 'edf',
    ])

    for candidate in candidates:
        if candidate.exists():
            # Check if it has cycle subdirectories with .tif files
            cycle_dirs = list(candidate.glob("cyc*"))
            if cycle_dirs:
                tif_files = list(candidate.glob("cyc*/*.tif"))
                if tif_files:
                    print(f"  Found EDF directory: {candidate}")
                    print(f"  Contains {len(tif_files)} files across {len(cycle_dirs)} cycles")
                    return candidate

    return default_edf_dir or (project_dir / 'data' / 'processed' / 'edf')


def discover_edf_files(
    edf_directory: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    channel_name_dict: Dict[int, List[str]] = None,
) -> Dict[Tuple[int, int], Tuple[str, str]]:
    """Discover actual EDF files and build a mapping of cycle/channel to file paths.

    Parameters
    ----------
    edf_directory : Path
        Directory containing EDF output
    start_cycle, end_cycle : int
        Cycle range to search
    channel_name_dict : dict, optional
        Mapping of cycle -> list of channel names. Used to determine correct
        channel index from marker name. If not provided, channel index is
        assigned based on alphabetical file order (may be incorrect).

    Returns
    -------
    dict
        Mapping of (cycle, channel) -> (filepath, marker_name)
    """
    file_map = {}

    for cycle in range(start_cycle, end_cycle + 1):
        cycle_dir = Path(edf_directory) / f"cyc{cycle:02d}"
        if not cycle_dir.exists():
            continue

        # Get all .tif files in this cycle directory
        tif_files = sorted(cycle_dir.glob("*.tif"))

        # Get channel names for this cycle if available
        cycle_channels = channel_name_dict.get(cycle, []) if channel_name_dict else []

        for tif_path in tif_files:
            # Extract marker name from filename (remove .tif extension)
            marker_name = tif_path.stem

            # Try to find the actual channel index from channel_name_dict
            channel_idx = None
            if cycle_channels:
                # Look for this marker in the channel list
                for idx, ch_name in enumerate(cycle_channels, start=1):
                    if ch_name == marker_name:
                        channel_idx = idx
                        break

            # If not found in dict, fall back to position in sorted list
            if channel_idx is None:
                channel_idx = list(tif_files).index(tif_path) + 1

            file_map[(cycle, channel_idx)] = (str(tif_path), marker_name)

    return file_map


def compute_edf_stats_gpu(
    edf_path: Union[str, Path],
    device_id: int = 0,
) -> Optional[Dict[str, float]]:
    """Compute statistics for a single EDF projection using GPU.

    Parameters
    ----------
    edf_path : Path
        Path to EDF image file
    device_id : int
        GPU device ID to use

    Returns
    -------
    dict or None
        Statistics dictionary or None if file not found
    """
    import cupy as cp

    if not os.path.exists(edf_path):
        return None

    try:
        with cp.cuda.Device(device_id):
            img = imread(edf_path)
            img_gpu = cp.asarray(img)

            overall_mean = float(cp.mean(img_gpu))
            overall_std = float(cp.std(img_gpu))
            snr = overall_mean / overall_std if overall_std > 0 else 0.0

            # Downsample for profile computation
            step = 4 if img.shape[0] > 2000 else 1
            ds_gpu = img_gpu[::step, ::step]
            row_profile = cp.mean(ds_gpu, axis=1)
            col_profile = cp.mean(ds_gpu, axis=0)

            row_mean = cp.mean(row_profile)
            col_mean = cp.mean(col_profile)
            row_cv = float(cp.std(row_profile) / row_mean * 100) if float(row_mean) > 0 else 0.0
            col_cv = float(cp.std(col_profile) / col_mean * 100) if float(col_mean) > 0 else 0.0

            del img_gpu, ds_gpu, row_profile, col_profile
            cp.get_default_memory_pool().free_all_blocks()

        return {
            'overall_mean': overall_mean,
            'overall_std': overall_std,
            'snr': snr,
            'row_cv': row_cv,
            'col_cv': col_cv
        }
    except Exception as e:
        return {'error': str(e)}


def _process_single_edf_stat(args: Tuple) -> Optional[Dict]:
    """Worker function for parallel EDF stat collection."""
    cycle, channel, edf_path, marker_name, device_id = args
    stats = compute_edf_stats_gpu(edf_path, device_id)
    if stats and 'error' not in stats:
        return {
            'cycle': cycle,
            'channel': channel,
            'channel_name': marker_name,
            **stats
        }
    return None


def collect_edf_stats_parallel(
    file_map: Dict[Tuple[int, int], Tuple[str, str]],
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
) -> pd.DataFrame:
    """Collect EDF statistics using GPU-accelerated parallel processing.

    Parameters
    ----------
    file_map : dict
        Mapping of (cycle, channel) -> (filepath, marker_name)
        from discover_edf_files()
    gpu_device_ids : list
        List of GPU device IDs
    io_workers : int
        Number of parallel I/O workers
    zplanes_per_gpu : int
        Concurrency factor per GPU

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: cycle, channel, channel_name,
        overall_mean, overall_std, snr, row_cv, col_cv
    """
    gpu_ids = gpu_device_ids or [0]
    num_gpus = len(gpu_ids)

    work_items = []
    idx = 0
    for (cycle, channel), (edf_path, marker_name) in file_map.items():
        device_id = gpu_ids[idx % num_gpus]
        work_items.append((cycle, channel, edf_path, marker_name, device_id))
        idx += 1

    total = len(work_items)
    if total == 0:
        return pd.DataFrame()

    all_stats = []

    max_workers = min(io_workers, num_gpus * zplanes_per_gpu)
    print(f"  Using {max_workers} parallel workers across {num_gpus} GPU(s)")

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {executor.submit(_process_single_edf_stat, item): item for item in work_items}

        completed = 0
        for future in as_completed(futures):
            result = future.result()
            if result:
                all_stats.append(result)
            completed += 1
            if completed % 10 == 0 or completed == total:
                print(f"\r  Progress: {completed}/{total} ({100*completed/total:.0f}%)", end='', flush=True)

    print()
    cleanup_gpu_memory()

    df = pd.DataFrame(all_stats)
    if len(df) > 0:
        df = df.sort_values(['cycle', 'channel']).reset_index(drop=True)
    return df


# =============================================================================
# VISUALIZATION FUNCTIONS
# =============================================================================

def plot_summary_heatmaps(
    stats_df: pd.DataFrame,
    stage_name: str = "raw",
    qc_output_dir: Union[str, Path] = None,
    save_pdf: bool = True,
    show: bool = True,
) -> plt.Figure:
    """Create summary heatmaps: SNR and CV by cycle/channel.

    Parameters
    ----------
    stats_df : pd.DataFrame
        DataFrame with statistics (must have cycle, channel, snr, row_cv, overall_mean)
    stage_name : str
        Name for the processing stage (used in title and PDF filename)
    qc_output_dir : Path, optional
        Directory for saving PDF output
    save_pdf : bool
        Whether to save the plot as PDF
    show : bool
        Whether to display the plot

    Returns
    -------
    plt.Figure
        The matplotlib figure
    """
    summary = stats_df.groupby(['cycle', 'channel']).agg({
        'snr': 'mean',
        'row_cv': 'mean',
        'col_cv': 'mean',
        'overall_mean': 'mean'
    }).reset_index()

    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # SNR heatmap
    pivot_snr = summary.pivot(index='cycle', columns='channel', values='snr')
    sns.heatmap(pivot_snr, annot=True, fmt='.2f', cmap='viridis', ax=axes[0])
    axes[0].set_title('Mean SNR by Cycle/Channel')

    # Row CV heatmap
    pivot_cv = summary.pivot(index='cycle', columns='channel', values='row_cv')
    sns.heatmap(pivot_cv, annot=True, fmt='.1f', cmap='RdYlGn_r', ax=axes[1])
    axes[1].set_title('Row CV (%) by Cycle/Channel')

    # Mean intensity heatmap
    pivot_mean = summary.pivot(index='cycle', columns='channel', values='overall_mean')
    sns.heatmap(pivot_mean, annot=True, fmt='.0f', cmap='plasma', ax=axes[2])
    axes[2].set_title('Mean Intensity by Cycle/Channel')

    plt.suptitle(f'{stage_name.upper()} Summary Heatmaps', y=1.02, fontsize=14)
    plt.tight_layout()

    if save_pdf and qc_output_dir:
        qc_output_dir = Path(qc_output_dir)
        qc_output_dir.mkdir(parents=True, exist_ok=True)
        pdf_path = qc_output_dir / f'{stage_name}_summary_heatmaps.pdf'
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=150)
        print(f"  Saved: {pdf_path}")

    if show:
        plt.show()

    return fig


def plot_zplane_profiles(
    stats_df: pd.DataFrame,
    cycle: int,
    channel: int,
    stage_name: str = "raw",
    qc_output_dir: Union[str, Path] = None,
    save_pdf: bool = True,
    show: bool = True,
) -> Optional[plt.Figure]:
    """Plot statistics across all z-planes for one cycle/channel.

    Parameters
    ----------
    stats_df : pd.DataFrame
        DataFrame with statistics (must have cycle, channel, zplane, channel_name,
        overall_mean, snr, row_cv, col_cv, overall_std)
    cycle : int
        Cycle number
    channel : int
        Channel number
    stage_name : str
        Name for the processing stage
    qc_output_dir : Path, optional
        Directory for saving PDF output
    save_pdf : bool
        Whether to save the plot as PDF
    show : bool
        Whether to display the plot

    Returns
    -------
    plt.Figure or None
        The matplotlib figure, or None if no data found
    """
    subset = stats_df[(stats_df['cycle'] == cycle) & (stats_df['channel'] == channel)].sort_values('zplane')
    if len(subset) == 0:
        return None

    ch_name = subset['channel_name'].iloc[0]

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))

    axes[0].plot(subset['zplane'], subset['overall_mean'], 'b-o', markersize=4)
    axes[0].set_xlabel('Z-plane')
    axes[0].set_ylabel('Mean Intensity')
    axes[0].set_title(f'Cyc{cycle} {ch_name}: Intensity')
    axes[0].grid(True, alpha=0.3)

    axes[1].plot(subset['zplane'], subset['snr'], 'g-o', markersize=4)
    axes[1].set_xlabel('Z-plane')
    axes[1].set_ylabel('SNR')
    axes[1].set_title('SNR across Z')
    axes[1].grid(True, alpha=0.3)

    axes[2].plot(subset['zplane'], subset['row_cv'], 'r-o', markersize=4, label='Row CV')
    axes[2].plot(subset['zplane'], subset['col_cv'], 'm-o', markersize=4, label='Col CV')
    axes[2].set_xlabel('Z-plane')
    axes[2].set_ylabel('CV (%)')
    axes[2].set_title('Structure (CV)')
    axes[2].legend()
    axes[2].grid(True, alpha=0.3)

    axes[3].plot(subset['zplane'], subset['overall_std'], color='orange', marker='o', markersize=4)
    axes[3].set_xlabel('Z-plane')
    axes[3].set_ylabel('Std Dev')
    axes[3].set_title('Noise (Std Dev)')
    axes[3].grid(True, alpha=0.3)

    plt.tight_layout()

    if save_pdf and qc_output_dir:
        qc_output_dir = Path(qc_output_dir)
        qc_output_dir.mkdir(parents=True, exist_ok=True)
        pdf_path = qc_output_dir / f'{stage_name}_zprofile_cyc{cycle:02d}_ch{channel}.pdf'
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=150)

    if show:
        plt.show()

    return fig


def plot_edf_summary_heatmaps(
    edf_stats_df: pd.DataFrame,
    qc_output_dir: Union[str, Path] = None,
    save_pdf: bool = True,
    show: bool = True,
) -> plt.Figure:
    """Create EDF-specific summary heatmaps.

    Parameters
    ----------
    edf_stats_df : pd.DataFrame
        DataFrame with EDF statistics
    qc_output_dir : Path, optional
        Directory for saving PDF output
    save_pdf : bool
        Whether to save the plot as PDF
    show : bool
        Whether to display the plot

    Returns
    -------
    plt.Figure
        The matplotlib figure
    """
    fig, axes = plt.subplots(1, 3, figsize=(18, 6))

    # SNR by cycle/channel
    pivot_snr = edf_stats_df.pivot(index='cycle', columns='channel', values='snr')
    sns.heatmap(pivot_snr, annot=True, fmt='.2f', cmap='viridis', ax=axes[0])
    axes[0].set_title('EDF SNR by Cycle/Channel')

    # Row CV
    pivot_cv = edf_stats_df.pivot(index='cycle', columns='channel', values='row_cv')
    sns.heatmap(pivot_cv, annot=True, fmt='.1f', cmap='RdYlGn_r', ax=axes[1])
    axes[1].set_title('EDF Row CV (%) by Cycle/Channel')

    # Mean intensity
    pivot_mean = edf_stats_df.pivot(index='cycle', columns='channel', values='overall_mean')
    sns.heatmap(pivot_mean, annot=True, fmt='.0f', cmap='plasma', ax=axes[2])
    axes[2].set_title('EDF Mean Intensity by Cycle/Channel')

    plt.suptitle("EDF Projections", y=1.02, fontsize=14)
    plt.tight_layout()

    if save_pdf and qc_output_dir:
        qc_output_dir = Path(qc_output_dir)
        qc_output_dir.mkdir(parents=True, exist_ok=True)
        pdf_path = qc_output_dir / 'edf_summary_heatmaps.pdf'
        fig.savefig(pdf_path, format='pdf', bbox_inches='tight', dpi=150)
        print(f"  Saved: {pdf_path}")

    if show:
        plt.show()

    return fig


# =============================================================================
# CONVENIENCE FUNCTIONS
# =============================================================================

def run_raw_qc(
    image_dir: Union[str, Path],
    cache_file: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    qc_output_dir: Union[str, Path] = None,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Run complete raw data QC analysis with caching.

    Parameters
    ----------
    image_dir : Path
        Directory containing raw images
    cache_file : Path
        Path to cache file for statistics
    start_cycle, end_cycle : int
        Cycle range
    start_channel, end_channel : int
        Channel range
    n_zplanes : int
        Number of z-planes
    qc_output_dir : Path, optional
        Directory for QC plot output
    gpu_device_ids : list
        GPU device IDs
    io_workers : int
        I/O worker count
    zplanes_per_gpu : int
        Concurrency per GPU
    channel_name_dict : dict
        Channel name mapping

    Returns
    -------
    pd.DataFrame
        DataFrame with raw statistics
    """
    print("="*70)
    print("COMPREHENSIVE RAW DATA ANALYSIS")
    print("="*70)
    print(f"Analyzing: Cycles {start_cycle}-{end_cycle}, Channels {start_channel}-{end_channel}, Z-planes 1-{n_zplanes}")
    if qc_output_dir:
        print(f"QC plots will be saved to: {qc_output_dir}")
    print()

    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists():
        print(f"Loading cached statistics from {cache_path}")
        with open(cache_path, 'rb') as f:
            stats_df = pickle.load(f)
        print(f"  Loaded {len(stats_df)} cached data points")
    else:
        print("Computing statistics (GPU-accelerated parallel)...")
        stats_df = collect_raw_stats_parallel(
            image_dir, start_cycle, end_cycle, start_channel, end_channel, n_zplanes,
            gpu_device_ids=gpu_device_ids,
            io_workers=io_workers,
            zplanes_per_gpu=zplanes_per_gpu,
            channel_name_dict=channel_name_dict,
        )

        with open(cache_path, 'wb') as f:
            pickle.dump(stats_df, f)
        print(f"  Cached to {cache_path}")

    if len(stats_df) > 0:
        print(f"\nCollected {len(stats_df)} data points")

        # Summary statistics
        print("\n" + "="*70)
        print("OVERALL STATISTICS SUMMARY")
        print("="*70)
        summary = stats_df.groupby(['cycle', 'channel']).agg({
            'overall_mean': ['mean', 'std'],
            'snr': ['mean', 'min', 'max'],
            'row_cv': 'mean',
            'col_cv': 'mean'
        }).round(2)
        print(summary.to_string())

        # Summary heatmaps
        print("\n" + "="*70)
        print("SUMMARY HEATMAPS")
        print("="*70)
        plot_summary_heatmaps(stats_df, stage_name="raw", qc_output_dir=qc_output_dir)

        # Individual channel profiles
        print("\n" + "="*70)
        print("INDIVIDUAL CHANNEL PROFILES (All Z-planes)")
        print("="*70)
        for cycle in range(start_cycle, end_cycle + 1):
            for channel in range(start_channel, end_channel + 1):
                save_pdf = (cycle == start_cycle)
                plot_zplane_profiles(
                    stats_df, cycle, channel,
                    stage_name="raw",
                    qc_output_dir=qc_output_dir,
                    save_pdf=save_pdf
                )
    else:
        print("No data collected - check image_dir and file patterns")

    print("\n" + "="*70)
    print("Raw data analysis complete")
    print("="*70)

    return stats_df


def run_stitched_qc(
    stitch_dir: Union[str, Path],
    cache_file: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    qc_output_dir: Union[str, Path] = None,
    raw_stats_df: pd.DataFrame = None,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Run complete stitched data QC analysis with caching and comparison."""
    print("="*70)
    print("COMPREHENSIVE STITCHED DATA ANALYSIS")
    print("="*70)
    print(f"Analyzing: Cycles {start_cycle}-{end_cycle}, Channels {start_channel}-{end_channel}, Z-planes 1-{n_zplanes}")
    if qc_output_dir:
        print(f"QC plots will be saved to: {qc_output_dir}")
    print()

    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists():
        print(f"Loading cached statistics from {cache_path}")
        with open(cache_path, 'rb') as f:
            stats_df = pickle.load(f)
        print(f"  Loaded {len(stats_df)} cached data points")
    else:
        print("Computing statistics (GPU-accelerated parallel)...")
        stats_df = collect_stitched_stats_parallel(
            stitch_dir, start_cycle, end_cycle, start_channel, end_channel, n_zplanes,
            gpu_device_ids=gpu_device_ids,
            io_workers=io_workers,
            zplanes_per_gpu=zplanes_per_gpu,
            channel_name_dict=channel_name_dict,
        )

        with open(cache_path, 'wb') as f:
            pickle.dump(stats_df, f)
        print(f"  Cached to {cache_path}")

    if len(stats_df) > 0:
        print(f"\nCollected {len(stats_df)} data points")

        # Summary heatmaps
        print("\n" + "="*70)
        print("STITCHED SUMMARY HEATMAPS")
        print("="*70)
        plot_summary_heatmaps(stats_df, stage_name="stitched", qc_output_dir=qc_output_dir)

        # Compare with raw if available
        if raw_stats_df is not None and len(raw_stats_df) > 0:
            print("\n" + "="*70)
            print("RAW vs STITCHED COMPARISON")
            print("="*70)

            raw_agg = raw_stats_df.groupby(['cycle', 'channel']).agg({
                'overall_mean': 'mean', 'snr': 'mean', 'row_cv': 'mean', 'col_cv': 'mean'
            }).add_prefix('raw_')

            stitched_agg = stats_df.groupby(['cycle', 'channel']).agg({
                'overall_mean': 'mean', 'snr': 'mean', 'row_cv': 'mean', 'col_cv': 'mean'
            }).add_prefix('stitched_')

            comparison = raw_agg.join(stitched_agg)
            comparison['snr_change'] = comparison['stitched_snr'] - comparison['raw_snr']
            comparison['cv_change'] = comparison['stitched_row_cv'] - comparison['raw_row_cv']

            print("\nSNR Improvement (Stitched - Raw):")
            print(comparison[['snr_change']].round(2).to_string())

        # Individual channel profiles
        print("\n" + "="*70)
        print("STITCHED CHANNEL PROFILES (All Z-planes)")
        print("="*70)
        for cycle in range(start_cycle, end_cycle + 1):
            for channel in range(start_channel, end_channel + 1):
                plot_zplane_profiles(
                    stats_df, cycle, channel,
                    stage_name="stitched",
                    qc_output_dir=qc_output_dir,
                    save_pdf=True
                )
    else:
        print("No stitched images found - run stitching first")

    print("\n" + "="*70)
    print("Stitched data analysis complete")
    print("="*70)

    return stats_df


def run_decon_qc(
    decon_dir: Union[str, Path],
    cache_file: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    start_channel: int,
    end_channel: int,
    n_zplanes: int,
    qc_output_dir: Union[str, Path] = None,
    stitched_stats_df: pd.DataFrame = None,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Run complete deconvolved data QC analysis with caching and comparison."""
    print("="*70)
    print("COMPREHENSIVE DECONVOLVED DATA ANALYSIS")
    print("="*70)
    print(f"Analyzing: Cycles {start_cycle}-{end_cycle}, Channels {start_channel}-{end_channel}, Z-planes 1-{n_zplanes}")
    if qc_output_dir:
        print(f"QC plots will be saved to: {qc_output_dir}")
    print()

    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    if cache_path.exists():
        print(f"Loading cached statistics from {cache_path}")
        with open(cache_path, 'rb') as f:
            stats_df = pickle.load(f)
        print(f"  Loaded {len(stats_df)} cached data points")
    else:
        print("Computing statistics (GPU-accelerated parallel)...")
        stats_df = collect_decon_stats_parallel(
            decon_dir, start_cycle, end_cycle, start_channel, end_channel, n_zplanes,
            gpu_device_ids=gpu_device_ids,
            io_workers=io_workers,
            zplanes_per_gpu=zplanes_per_gpu,
            channel_name_dict=channel_name_dict,
        )

        with open(cache_path, 'wb') as f:
            pickle.dump(stats_df, f)
        print(f"  Cached to {cache_path}")

    if len(stats_df) > 0:
        print(f"\nCollected {len(stats_df)} data points")

        # Summary heatmaps
        print("\n" + "="*70)
        print("DECONVOLVED SUMMARY HEATMAPS")
        print("="*70)
        plot_summary_heatmaps(stats_df, stage_name="deconvolved", qc_output_dir=qc_output_dir)

        # Compare with stitched if available
        if stitched_stats_df is not None and len(stitched_stats_df) > 0:
            print("\n" + "="*70)
            print("STITCHED vs DECONVOLVED COMPARISON")
            print("="*70)

            stitched_agg = stitched_stats_df.groupby(['cycle', 'channel']).agg({
                'overall_mean': 'mean', 'snr': 'mean', 'row_cv': 'mean'
            }).add_prefix('stitched_')

            decon_agg = stats_df.groupby(['cycle', 'channel']).agg({
                'overall_mean': 'mean', 'snr': 'mean', 'row_cv': 'mean'
            }).add_prefix('decon_')

            comparison = stitched_agg.join(decon_agg)
            comparison['snr_change'] = comparison['decon_snr'] - comparison['stitched_snr']

            print("\nSNR Improvement (Deconvolved - Stitched):")
            print(comparison[['snr_change']].round(2).to_string())

        # Individual channel profiles
        print("\n" + "="*70)
        print("DECONVOLVED CHANNEL PROFILES (All Z-planes)")
        print("="*70)
        for cycle in range(start_cycle, end_cycle + 1):
            for channel in range(start_channel, end_channel + 1):
                plot_zplane_profiles(
                    stats_df, cycle, channel,
                    stage_name="deconvolved",
                    qc_output_dir=qc_output_dir,
                    save_pdf=True
                )
    else:
        print("No deconvolved images found - run deconvolution first")

    print("\n" + "="*70)
    print("Deconvolved data analysis complete")
    print("="*70)

    return stats_df


def run_edf_qc(
    project_dir: Union[str, Path],
    cache_file: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    n_zplanes: int,
    qc_output_dir: Union[str, Path] = None,
    edf_dir: Union[str, Path] = None,
    decon_stats_df: pd.DataFrame = None,
    gpu_device_ids: List[int] = None,
    io_workers: int = 4,
    zplanes_per_gpu: int = 2,
    channel_name_dict: Dict[int, List[str]] = None,
) -> pd.DataFrame:
    """Run complete EDF data QC analysis with caching and comparison."""
    print("="*70)
    print("COMPREHENSIVE EDF PROJECTION ANALYSIS")
    print("="*70)
    print(f"Analyzing: Cycles {start_cycle}-{end_cycle}")
    if qc_output_dir:
        print(f"QC plots will be saved to: {qc_output_dir}")
    print()

    # Find actual EDF directory
    actual_edf_dir = find_edf_directory(project_dir, edf_dir)

    # Discover actual EDF files
    print("Discovering EDF files...")
    file_map = discover_edf_files(actual_edf_dir, start_cycle, end_cycle, channel_name_dict)
    print(f"  Found {len(file_map)} EDF files")

    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)

    use_cache = False
    if cache_path.exists():
        with open(cache_path, 'rb') as f:
            cached_df = pickle.load(f)
        if len(cached_df) > 0:
            print(f"Loading cached statistics from {cache_path}")
            print(f"  Loaded {len(cached_df)} cached data points")
            stats_df = cached_df
            use_cache = True
        else:
            print(f"  Cache is empty but {len(file_map)} files exist - recomputing...")

    if not use_cache:
        if len(file_map) > 0:
            print("Computing statistics (GPU-accelerated parallel)...")
            stats_df = collect_edf_stats_parallel(
                file_map,
                gpu_device_ids=gpu_device_ids,
                io_workers=io_workers,
                zplanes_per_gpu=zplanes_per_gpu,
            )

            if len(stats_df) > 0:
                with open(cache_path, 'wb') as f:
                    pickle.dump(stats_df, f)
                print(f"  Cached to {cache_path}")
        else:
            stats_df = pd.DataFrame()

    if len(stats_df) > 0:
        print(f"\nCollected {len(stats_df)} EDF projections")

        # Summary heatmaps
        print("\n" + "="*70)
        print("EDF SUMMARY HEATMAPS")
        print("="*70)
        plot_edf_summary_heatmaps(stats_df, qc_output_dir=qc_output_dir)

        # Compare with deconvolved middle z-plane if available
        if decon_stats_df is not None and len(decon_stats_df) > 0:
            print("\n" + "="*70)
            print("DECONVOLVED (middle z) vs EDF COMPARISON")
            print("="*70)

            mid_z = n_zplanes // 2
            decon_mid = decon_stats_df[decon_stats_df['zplane'] == mid_z].groupby(['cycle', 'channel']).agg({
                'snr': 'mean'
            }).add_prefix('decon_mid_')

            edf_agg = stats_df.groupby(['cycle', 'channel']).agg({
                'snr': 'mean'
            }).add_prefix('edf_')

            comparison = decon_mid.join(edf_agg)
            comparison['snr_change'] = comparison['edf_snr'] - comparison['decon_mid_snr']

            print("\nSNR Improvement (EDF - Decon middle z):")
            print(comparison[['snr_change']].round(2).to_string())

        # Summary statistics table
        print("\n" + "="*70)
        print("EDF STATISTICS SUMMARY")
        print("="*70)
        print(stats_df[['cycle', 'channel', 'channel_name', 'overall_mean', 'snr', 'row_cv']].to_string(index=False))
    else:
        print("\nNo EDF projections found.")
        print("Run EDF processing first to generate EDF projections.")

    print("\n" + "="*70)
    print("EDF analysis complete")
    print("="*70)

    return stats_df


# ===========================================================================
# Registration QC
# ===========================================================================
# Green/magenta DAPI overlays at full resolution to evaluate nuclei overlap.
# No GPU needed — pyvips streaming I/O + numpy compositing + matplotlib.


def _discover_registered_dapi(registered_dir, start_cycle, end_cycle):
    """Find DAPI TIFFs in registered/cyc##/ directories.

    Handles naming conventions:
      - Cycle 1 (reference): DAPI-01.tif
      - Other cycles: DAPI_cyc02.tif, DAPI_cyc03.tif, etc.

    Returns:
        dict: {cycle_number: Path} mapping, e.g. {1: Path(".../cyc01/DAPI-01.tif")}
    """
    registered_dir = Path(registered_dir)
    dapi_files = {}

    for cyc in range(start_cycle, end_cycle + 1):
        cyc_dir = registered_dir / f"cyc{cyc:02d}"
        if not cyc_dir.exists():
            continue

        # Try common DAPI naming patterns
        candidates = [
            cyc_dir / f"DAPI-{cyc:02d}.tif",       # Reference cycle: DAPI-01.tif
            cyc_dir / f"DAPI_cyc{cyc:02d}.tif",     # Moving cycles: DAPI_cyc02.tif
            cyc_dir / "DAPI.tif",                     # Plain DAPI.tif
        ]

        for candidate in candidates:
            if candidate.exists():
                dapi_files[cyc] = candidate
                break
        else:
            # Fallback: glob for anything with DAPI in the name
            dapi_glob = sorted(cyc_dir.glob("DAPI*.tif"))
            if dapi_glob:
                dapi_files[cyc] = dapi_glob[0]

    return dapi_files


def _select_crop_positions(img_height, img_width, crop_size, n_crops=4):
    """Select deterministic crop positions with 10% edge margin.

    Returns list of (row, col) top-left corners for:
      center, upper-left, upper-right, lower-center
    """
    margin_y = int(img_height * 0.10)
    margin_x = int(img_width * 0.10)

    positions = []
    labels = []

    # Center
    cy = (img_height - crop_size) // 2
    cx = (img_width - crop_size) // 2
    positions.append((cy, cx))
    labels.append("Center")

    # Upper-left (with margin)
    positions.append((margin_y, margin_x))
    labels.append("Upper-left")

    # Upper-right (with margin)
    positions.append((margin_y, img_width - crop_size - margin_x))
    labels.append("Upper-right")

    # Lower-center
    positions.append((img_height - crop_size - margin_y, cx))
    labels.append("Lower-center")

    # Clamp to valid range and limit to n_crops
    valid_positions = []
    valid_labels = []
    for (r, c), label in zip(positions[:n_crops], labels[:n_crops]):
        r = max(0, min(r, img_height - crop_size))
        c = max(0, min(c, img_width - crop_size))
        valid_positions.append((r, c))
        valid_labels.append(label)

    return valid_positions, valid_labels


def _compute_spatial_ncc_grid(ref_path, mov_path, tile_size=512):
    """Compute local NCC between reference and moving images on a tile grid.

    Parameters
    ----------
    ref_path : str or Path
        Path to reference DAPI TIFF.
    mov_path : str or Path
        Path to moving DAPI TIFF.
    tile_size : int
        Side length of each square tile in pixels.

    Returns
    -------
    ncc_grid : np.ndarray
        2D array of NCC values, shape (n_rows, n_cols). NaN for zero-variance tiles.
    img_height : int
        Image height in pixels.
    img_width : int
        Image width in pixels.
    """
    import math

    import pyvips

    ref_img = pyvips.Image.new_from_file(str(ref_path), access="random")
    mov_img = pyvips.Image.new_from_file(str(mov_path), access="random")
    # Use minimum dimensions — registered images may differ slightly
    img_height = min(ref_img.height, mov_img.height)
    img_width = min(ref_img.width, mov_img.width)

    n_rows = math.ceil(img_height / tile_size)
    n_cols = math.ceil(img_width / tile_size)
    ncc_grid = np.full((n_rows, n_cols), np.nan, dtype=np.float64)

    for i in range(n_rows):
        y = i * tile_size
        h = min(tile_size, img_height - y)
        for j in range(n_cols):
            x = j * tile_size
            w = min(tile_size, img_width - x)

            ref_tile = np.ndarray(
                buffer=ref_img.crop(x, y, w, h).write_to_memory(),
                dtype=np.uint16, shape=[h, w],
            ).astype(np.float32)
            mov_tile = np.ndarray(
                buffer=mov_img.crop(x, y, w, h).write_to_memory(),
                dtype=np.uint16, shape=[h, w],
            ).astype(np.float32)

            r = ref_tile - ref_tile.mean()
            m = mov_tile - mov_tile.mean()
            r_norm = np.linalg.norm(r)
            m_norm = np.linalg.norm(m)

            # Skip zero-variance tiles (background)
            if r_norm < 1e-10 or m_norm < 1e-10:
                continue

            ncc_grid[i, j] = np.dot(r.ravel(), m.ravel()) / (r_norm * m_norm)

    return ncc_grid, img_height, img_width


def _plot_spatial_ncc_heatmaps(
    ncc_grids, moving_cycles, ref_cycle, img_height, img_width, tile_size,
    qc_output_dir,
):
    """Render per-cycle NCC heatmaps and save as a single PDF.

    Parameters
    ----------
    ncc_grids : dict
        {cycle_number: ncc_grid} from _compute_spatial_ncc_grid.
    moving_cycles : list of int
        Sorted list of moving cycle numbers.
    ref_cycle : int
        Reference cycle number.
    img_height, img_width : int
        Image dimensions in pixels.
    tile_size : int
        Tile size used for NCC computation.
    qc_output_dir : Path
        Directory for output PDF.

    Returns
    -------
    Path
        Path to the saved PDF.
    """
    import math

    from matplotlib.colors import Normalize
    from mpl_toolkits.axes_grid1 import make_axes_locatable

    n_moving = len(moving_cycles)
    n_cols = min(4, n_moving)
    n_rows = math.ceil(n_moving / n_cols)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(5 * n_cols, 4 * n_rows),
        squeeze=False,
    )

    cmap = plt.cm.RdYlGn.copy()
    cmap.set_bad("0.85")
    norm = Normalize(vmin=0.5, vmax=1.0)

    for idx, mov_cyc in enumerate(moving_cycles):
        row_idx = idx // n_cols
        col_idx = idx % n_cols
        ax = axes[row_idx, col_idx]
        grid = ncc_grids[mov_cyc]

        im = ax.imshow(
            grid, cmap=cmap, norm=norm,
            interpolation="bilinear", origin="upper",
            extent=[0, img_width, img_height, 0],
        )
        min_val = np.nanmin(grid)
        p5_val = np.nanpercentile(grid[~np.isnan(grid)], 5) if np.any(~np.isnan(grid)) else np.nan
        n_poor = int(np.nansum(grid < 0.7))
        ax.set_title(
            f"Cyc {mov_cyc}\u2192{ref_cycle}  "
            f"(min={min_val:.3f}, p5={p5_val:.3f}, poor={n_poor})",
            fontsize=9,
        )
        # Axis ticks at tile boundaries
        x_ticks = np.arange(0, img_width + 1, tile_size * 4)
        y_ticks = np.arange(0, img_height + 1, tile_size * 4)
        ax.set_xticks(x_ticks)
        ax.set_yticks(y_ticks)
        ax.tick_params(labelsize=7)

    # Hide unused axes
    for idx in range(n_moving, n_rows * n_cols):
        row_idx = idx // n_cols
        col_idx = idx % n_cols
        axes[row_idx, col_idx].set_visible(False)

    # Shared colorbar
    fig.subplots_adjust(right=0.92)
    cbar_ax = fig.add_axes([0.94, 0.15, 0.02, 0.7])
    fig.colorbar(
        plt.cm.ScalarMappable(norm=norm, cmap=cmap),
        cax=cbar_ax, label="NCC",
    )

    fig.suptitle(
        f"Spatial NCC Heatmap \u2014 tile={tile_size}px  "
        f"(red=poor, green=good, gray=background)",
        fontsize=12, y=1.02,
    )

    out_path = Path(qc_output_dir) / "registration_qc_spatial_ncc.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _plot_targeted_overlays(
    dapi_files, ref_cycle, ncc_grids, moving_cycles, tile_size, crop_size,
    n_targeted, qc_output_dir,
):
    """Render green/magenta overlays at the worst-NCC tile positions per cycle.

    Parameters
    ----------
    dapi_files : dict
        {cycle_number: Path} mapping to registered DAPI TIFFs.
    ref_cycle : int
        Reference cycle number.
    ncc_grids : dict
        {cycle_number: ncc_grid} from _compute_spatial_ncc_grid.
    moving_cycles : list of int
        Sorted list of moving cycle numbers.
    tile_size : int
        Tile size used for NCC computation.
    crop_size : int
        Size of the overlay crop in pixels.
    n_targeted : int
        Number of worst-NCC tiles to show per cycle.
    qc_output_dir : Path
        Directory for output PDF.

    Returns
    -------
    Path
        Path to the saved PDF.
    """
    import pyvips

    n_rows = len(moving_cycles)
    n_cols = n_targeted

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4 * n_cols, 4 * n_rows),
        squeeze=False,
    )

    ref_img = pyvips.Image.new_from_file(str(dapi_files[ref_cycle]), access="random")
    img_height = ref_img.height
    img_width = ref_img.width

    for row_idx, mov_cyc in enumerate(moving_cycles):
        grid = ncc_grids[mov_cyc]

        # Find n_targeted worst (lowest) non-NaN tiles
        valid_mask = ~np.isnan(grid)
        if not valid_mask.any():
            for col_idx in range(n_cols):
                axes[row_idx, col_idx].text(
                    0.5, 0.5, "No valid tiles", transform=axes[row_idx, col_idx].transAxes,
                    ha="center", va="center", fontsize=9,
                )
                axes[row_idx, col_idx].set_xticks([])
                axes[row_idx, col_idx].set_yticks([])
            continue

        flat_indices = np.argsort(grid[valid_mask])[:n_targeted]
        valid_coords = np.argwhere(valid_mask)
        worst_tiles = valid_coords[flat_indices]

        mov_img = pyvips.Image.new_from_file(str(dapi_files[mov_cyc]), access="random")

        for col_idx in range(n_cols):
            ax = axes[row_idx, col_idx]

            if col_idx >= len(worst_tiles):
                ax.set_visible(False)
                continue

            tile_i, tile_j = worst_tiles[col_idx]
            ncc_val = grid[tile_i, tile_j]

            # Center the crop on the worst tile
            center_y = tile_i * tile_size + tile_size // 2
            center_x = tile_j * tile_size + tile_size // 2
            y = center_y - crop_size // 2
            x = center_x - crop_size // 2

            # Clamp to image bounds
            y = max(0, min(y, img_height - crop_size))
            x = max(0, min(x, img_width - crop_size))

            try:
                ref_crop = _load_crop_pyvips(dapi_files[ref_cycle], y, x, crop_size)
                mov_crop = _load_crop_pyvips(dapi_files[mov_cyc], y, x, crop_size)
                overlay = _make_green_magenta(ref_crop, mov_crop)
                ax.imshow(overlay)
            except Exception as e:
                ax.text(
                    0.5, 0.5, f"Error:\n{e}", transform=ax.transAxes,
                    ha="center", va="center", fontsize=8,
                )

            ax.set_title(f"NCC={ncc_val:.3f} @({x},{y})", fontsize=9)
            ax.set_xticks([])
            ax.set_yticks([])

            if col_idx == 0:
                ax.set_ylabel(f"Cyc {mov_cyc}\u2192{ref_cycle}", fontsize=10)

    fig.suptitle(
        f"Targeted Overlays \u2014 {n_targeted} worst-NCC positions per cycle\n"
        f"Green=Cyc{ref_cycle} (ref), Magenta=moving",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    out_path = Path(qc_output_dir) / "registration_qc_targeted_overlays.pdf"
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    return out_path


def _load_crop_pyvips(tiff_path, row, col, crop_size):
    """Load a crop from a large TIFF via pyvips random access.

    Returns numpy array (crop_size, crop_size) as float32.
    """
    import pyvips

    img = pyvips.Image.new_from_file(str(tiff_path), access="random")
    # Extract crop region
    crop = img.crop(col, row, crop_size, crop_size)
    # Convert to numpy
    data = np.ndarray(
        buffer=crop.write_to_memory(),
        dtype=np.uint16,
        shape=[crop.height, crop.width],
    )
    return data.astype(np.float32)


def _normalize_uint8(img, plow=1, phigh=99.5):
    """Percentile-clip and rescale to uint8 for display."""
    lo = np.percentile(img, plow)
    hi = np.percentile(img, phigh)
    if hi <= lo:
        hi = lo + 1
    clipped = np.clip((img - lo) / (hi - lo), 0, 1)
    return (clipped * 255).astype(np.uint8)


def _make_green_magenta(ref_crop, mov_crop):
    """Create green/magenta overlay: green=reference, magenta=moving.

    White/gray = well-aligned (both signals overlap).
    Green fringe = reference-only signal.
    Magenta fringe = moving-only signal.
    """
    ref_u8 = _normalize_uint8(ref_crop)
    mov_u8 = _normalize_uint8(mov_crop)

    rgb = np.zeros((*ref_u8.shape, 3), dtype=np.uint8)
    rgb[:, :, 0] = mov_u8     # Red channel = moving (magenta = R+B)
    rgb[:, :, 1] = ref_u8     # Green channel = reference
    rgb[:, :, 2] = mov_u8     # Blue channel = moving (magenta = R+B)
    return rgb


def run_registration_qc(
    project_dir,
    cache_file,
    start_cycle,
    end_cycle,
    qc_output_dir=None,
    channel_name_dict=None,
    crop_size=1000,
    n_crops=4,
    spatial_tile_size=512,
    n_targeted_crops=3,
):
    """Generate registration QC: green/magenta DAPI overlays and spatial metrics.

    Produces:
      1. registration_qc_overlay_grid.pdf — Contact sheet of DAPI overlays
         Rows = moving cycles, Columns = crop positions
         Green = reference DAPI, Magenta = moving DAPI
         White/gray = good alignment, color fringing = misalignment

      2. registration_qc_metrics.pdf — Per-cycle spatial NCC quality chart
         showing min NCC and p5 NCC (worst-tile metrics, NOT averages).

      3. registration_qc_spatial_ncc.pdf — Per-cycle heatmap of local NCC
         on a tile grid, showing spatially localized alignment quality.

      4. registration_qc_targeted_overlays.pdf — Green/magenta overlays at
         the worst-NCC tile positions per cycle.

    Parameters
    ----------
    project_dir : str or Path
        Project root directory.
    cache_file : str or Path
        Path to write cache pickle (for consistency with other QC functions).
    start_cycle : int
        First cycle number.
    end_cycle : int
        Last cycle number.
    qc_output_dir : str or Path, optional
        Directory for QC output PDFs. Defaults to project_dir/qc_plots.
    channel_name_dict : dict, optional
        Not used for registration QC (DAPI only), kept for API consistency.
    crop_size : int
        Size of each square crop in pixels. Default 1000.
    n_crops : int
        Number of crop positions. Default 4.
    spatial_tile_size : int
        Tile size for spatial NCC grid computation. Default 512.
    n_targeted_crops : int
        Number of worst-NCC positions to show per cycle. Default 3.

    Returns
    -------
    pd.DataFrame or None
        Summary table with registration metrics per cycle.
    """
    project_dir = Path(project_dir)
    registered_dir = project_dir / "data" / "processed" / "registered"

    if qc_output_dir is None:
        qc_output_dir = project_dir / "qc_plots"
    qc_output_dir = Path(qc_output_dir)
    qc_output_dir.mkdir(parents=True, exist_ok=True)

    n_cycles = end_cycle - start_cycle + 1
    if n_cycles < 2:
        print("Single-cycle project — no registration QC needed.")
        return None

    # Discover DAPI files
    dapi_files = _discover_registered_dapi(registered_dir, start_cycle, end_cycle)
    if len(dapi_files) < 2:
        print(f"WARNING: Found only {len(dapi_files)} registered DAPI files, need >= 2")
        return None

    ref_cycle = start_cycle
    if ref_cycle not in dapi_files:
        print(f"WARNING: Reference cycle {ref_cycle} DAPI not found")
        return None

    print(f"Registration QC: {len(dapi_files)} cycles, reference=cycle {ref_cycle}")
    print(f"  Crop size: {crop_size}x{crop_size}, positions: {n_crops}")

    # -----------------------------------------------------------------------
    # Load VALIS summary CSV for metrics (optional)
    # -----------------------------------------------------------------------
    valis_csv = (
        registered_dir / "registration_data" / "edf" / "data" / "edf_summary.csv"
    )
    valis_df = None
    if valis_csv.exists():
        try:
            valis_df = pd.read_csv(valis_csv)
            print(f"  VALIS summary: {len(valis_df)} rows from {valis_csv.name}")
        except Exception as e:
            print(f"  WARNING: Could not read VALIS summary: {e}")

    # -----------------------------------------------------------------------
    # Determine image dimensions from reference DAPI
    # -----------------------------------------------------------------------
    import pyvips

    ref_img = pyvips.Image.new_from_file(str(dapi_files[ref_cycle]), access="random")
    img_height = ref_img.height
    img_width = ref_img.width
    print(f"  Image size: {img_width}x{img_height}")

    positions, pos_labels = _select_crop_positions(
        img_height, img_width, crop_size, n_crops
    )

    # -----------------------------------------------------------------------
    # Load reference crops (reused for all moving cycles)
    # -----------------------------------------------------------------------
    ref_crops = []
    for (r, c) in positions:
        ref_crops.append(_load_crop_pyvips(dapi_files[ref_cycle], r, c, crop_size))

    # -----------------------------------------------------------------------
    # Build overlay grid: rows=moving cycles, cols=crop positions
    # -----------------------------------------------------------------------
    moving_cycles = sorted([c for c in dapi_files if c != ref_cycle])

    if not moving_cycles:
        print("No moving cycles to compare.")
        return None

    n_rows = len(moving_cycles)
    n_cols = len(positions)

    fig, axes = plt.subplots(
        n_rows, n_cols,
        figsize=(4 * n_cols, 4 * n_rows),
        squeeze=False,
    )

    for row_idx, mov_cyc in enumerate(moving_cycles):
        for col_idx, ((r, c), label) in enumerate(zip(positions, pos_labels)):
            ax = axes[row_idx, col_idx]
            try:
                mov_crop = _load_crop_pyvips(dapi_files[mov_cyc], r, c, crop_size)
                overlay = _make_green_magenta(ref_crops[col_idx], mov_crop)
                ax.imshow(overlay)
            except Exception as e:
                ax.text(
                    0.5, 0.5, f"Error:\n{e}", transform=ax.transAxes,
                    ha='center', va='center', fontsize=8,
                )

            ax.set_xticks([])
            ax.set_yticks([])

            if row_idx == 0:
                ax.set_title(label, fontsize=10)
            if col_idx == 0:
                ax.set_ylabel(f"Cyc {mov_cyc} \u2192 {ref_cycle}", fontsize=10)

    fig.suptitle(
        f"Registration QC \u2014 Green=Cyc{ref_cycle} (ref), Magenta=moving\n"
        f"White/gray=aligned, color fringing=misaligned",
        fontsize=12, y=1.02,
    )
    fig.tight_layout()
    overlay_path = qc_output_dir / "registration_qc_overlay_grid.pdf"
    fig.savefig(overlay_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved overlay grid: {overlay_path}")

    # -----------------------------------------------------------------------
    # Spatial NCC heatmaps — tile-level alignment quality
    # -----------------------------------------------------------------------
    print("\n  Computing spatial alignment heatmaps...")
    ncc_grids = {}
    for mov_cyc in moving_cycles:
        ncc_grid, _, _ = _compute_spatial_ncc_grid(
            dapi_files[ref_cycle], dapi_files[mov_cyc],
            tile_size=spatial_tile_size,
        )
        ncc_grids[mov_cyc] = ncc_grid
        valid = ncc_grid[~np.isnan(ncc_grid)]
        min_ncc = np.nanmin(ncc_grid)
        p5_ncc = np.percentile(valid, 5) if len(valid) > 0 else np.nan
        p25_ncc = np.percentile(valid, 25) if len(valid) > 0 else np.nan
        n_poor = int(np.sum(valid < 0.7))
        print(f"    Cyc {mov_cyc}: min={min_ncc:.4f}, p5={p5_ncc:.4f}, p25={p25_ncc:.4f}, n_poor(<0.7)={n_poor}")

    ncc_path = _plot_spatial_ncc_heatmaps(
        ncc_grids, moving_cycles, ref_cycle,
        img_height, img_width, spatial_tile_size,
        qc_output_dir,
    )
    print(f"  Saved spatial NCC heatmap: {ncc_path}")

    # -----------------------------------------------------------------------
    # Targeted overlays at worst-NCC positions
    # -----------------------------------------------------------------------
    print(f"\n  Generating targeted overlays ({n_targeted_crops} worst per cycle)...")
    targeted_path = _plot_targeted_overlays(
        dapi_files, ref_cycle, ncc_grids, moving_cycles,
        spatial_tile_size, crop_size, n_targeted_crops,
        qc_output_dir,
    )
    print(f"  Saved targeted overlays: {targeted_path}")

    # -----------------------------------------------------------------------
    # Per-cycle spatial NCC summary (replaces average-based VALIS metrics)
    # -----------------------------------------------------------------------
    ncc_threshold = 0.7
    summary_rows = []
    for mov_cyc in moving_cycles:
        grid = ncc_grids[mov_cyc]
        valid = grid[~np.isnan(grid)]
        n_valid = len(valid)
        summary_rows.append({
            "cycle": mov_cyc,
            "min_ncc": np.min(valid) if n_valid > 0 else np.nan,
            "p5_ncc": np.percentile(valid, 5) if n_valid > 0 else np.nan,
            "p25_ncc": np.percentile(valid, 25) if n_valid > 0 else np.nan,
            "median_ncc": np.median(valid) if n_valid > 0 else np.nan,
            "n_poor": int(np.sum(valid < ncc_threshold)) if n_valid > 0 else 0,
            "n_tiles": n_valid,
            "pct_poor": 100.0 * np.sum(valid < ncc_threshold) / n_valid if n_valid > 0 else 0.0,
        })

    summary_df = pd.DataFrame(summary_rows)
    summary_df = summary_df.sort_values("cycle")

    # Print spatial quality table
    print(f"\n{'=' * 85}")
    print(f"SPATIAL REGISTRATION QUALITY (NCC threshold={ncc_threshold})")
    print(f"{'=' * 85}")
    print(
        summary_df[["cycle", "min_ncc", "p5_ncc", "p25_ncc", "median_ncc", "n_poor", "pct_poor"]].to_string(
            index=False,
            float_format=lambda x: f"{x:.4f}" if not np.isnan(x) else "NaN",
        )
    )

    # Flag cycles with poor regions
    worst_cycles = summary_df[summary_df["n_poor"] > 0]
    if len(worst_cycles) > 0:
        print(f"\n  WARNING: {len(worst_cycles)} cycle(s) have tiles with NCC < {ncc_threshold}")
        for _, row in worst_cycles.iterrows():
            print(f"    Cyc {int(row['cycle'])}: {int(row['n_poor'])} poor tiles ({row['pct_poor']:.1f}%), min={row['min_ncc']:.4f}")
    else:
        print(f"\n  All tiles above NCC threshold ({ncc_threshold}) — registration looks uniform")

    # Bar chart: per-cycle min NCC and p5 NCC (spatial quality, not averages)
    if len(summary_df) > 0:
        fig2, ax2 = plt.subplots(figsize=(max(6, len(summary_df) * 0.8), 4))
        x = np.arange(len(summary_df))
        width = 0.35
        bars_min = ax2.bar(
            x - width / 2, summary_df["min_ncc"], width,
            label="min NCC (worst tile)",
            color=["#e74c3c" if v < ncc_threshold else "#f39c12" for v in summary_df["min_ncc"]],
            edgecolor="black", linewidth=0.5,
        )
        bars_p5 = ax2.bar(
            x + width / 2, summary_df["p5_ncc"], width,
            label="p5 NCC (5th percentile)",
            color=["#e74c3c" if v < ncc_threshold else "#2ecc71" for v in summary_df["p5_ncc"]],
            edgecolor="black", linewidth=0.5,
        )
        ax2.axhline(y=ncc_threshold, color="red", linestyle="--", linewidth=1, label=f"NCC={ncc_threshold} threshold")
        ax2.set_xticks(x)
        ax2.set_xticklabels([f"Cyc {int(c)}" for c in summary_df["cycle"]], fontsize=8)
        ax2.set_ylabel("NCC")
        ax2.set_title("Spatial Registration Quality per Cycle (worst-tile metrics)")
        ax2.set_ylim(0, 1.05)
        ax2.legend(fontsize=8)
        fig2.tight_layout()
        metrics_path = qc_output_dir / "registration_qc_metrics.pdf"
        fig2.savefig(metrics_path, dpi=150, bbox_inches="tight")
        plt.close(fig2)
        print(f"  Saved spatial quality chart: {metrics_path}")

    # Cache the summary
    cache_path = Path(cache_file)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "wb") as f:
        pickle.dump(summary_df, f)

    return summary_df
