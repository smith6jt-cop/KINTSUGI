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
) -> Dict[Tuple[int, int], Tuple[str, str]]:
    """Discover actual EDF files and build a mapping of cycle/channel to file paths.

    Parameters
    ----------
    edf_directory : Path
        Directory containing EDF output
    start_cycle, end_cycle : int
        Cycle range to search

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

        for channel_idx, tif_path in enumerate(tif_files, start=1):
            # Extract marker name from filename (remove .tif extension)
            marker_name = tif_path.stem
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
                save_pdf = (cycle == start_cycle)
                plot_zplane_profiles(
                    stats_df, cycle, channel,
                    stage_name="stitched",
                    qc_output_dir=qc_output_dir,
                    save_pdf=save_pdf
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
                save_pdf = (cycle == start_cycle)
                plot_zplane_profiles(
                    stats_df, cycle, channel,
                    stage_name="deconvolved",
                    qc_output_dir=qc_output_dir,
                    save_pdf=save_pdf
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
    file_map = discover_edf_files(actual_edf_dir, start_cycle, end_cycle)
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
