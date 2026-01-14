"""
KINTSUGI I/O, Channel Name, and Parallel Processing Utilities

This module contains functions for:
- Loading and saving channel name configurations
- Extracting channels from OME-TIFF files
- Image resize utilities for processing
- Parallel I/O utilities (tile loading, z-stack loading)
- Progress tracking for long-running operations

Location: notebooks/Kio.py
"""

import gc
import os
import re
import sys
import threading
import time
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from datetime import datetime
from glob import glob
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple, Union

import numpy as np
import pandas as pd
import tifffile as tiff
from skimage.io import imread, imsave
from skimage.transform import resize as skresize


# =============================================================================
# CHANNEL NAME MANAGEMENT
# =============================================================================

def load_channel_names(
    meta_dir: Union[str, Path],
    filename: str = "CHANNELNAMES.txt",
    channels_per_cycle: int = 4,
    make_unique: bool = True
) -> Optional[Dict[int, List[str]]]:
    """
    Load channel names from a text file in the metadata directory.

    Supports multiple formats:
    1. Simple list (one channel per line, with DAPI-01, DAPI-02 marking cycles):
       DAPI-01
       Blank
       Blank
       Blank
       DAPI-02
       CD31
       ...

    2. Cycle-prefixed (colon-separated): "1: DAPI, CD31, CD8, Empty"
    3. Tab-separated: "1\\tDAPI\\tCD31\\tCD8\\tEmpty"
    4. CSV style: "1,DAPI,CD31,CD8,Empty"

    Lines starting with # are treated as comments.

    Parameters
    ----------
    meta_dir : Path
        Path to metadata directory
    filename : str
        Name of channel names file (default: CHANNELNAMES.txt)
    channels_per_cycle : int
        Expected channels per cycle (for validation)
    make_unique : bool
        If True (default), automatically make channel names unique both
        within and across cycles. Set to False to get raw names from file.

    Returns
    -------
    dict or None
        Dictionary mapping cycle number (int) to list of channel names,
        or None if file not found. Names are made unique by default.

    Example
    -------
    >>> channel_dict = load_channel_names(meta_dir)
    >>> print(channel_dict[1])  # Duplicates are made unique
    ['DAPI-01', 'Blank', 'Blank_b', 'Blank_c']
    """
    channel_file = Path(meta_dir) / filename

    if not channel_file.exists():
        # Try alternate names
        alt_names = ["CHANNELNAMES.txt", "channelnames.txt", "channel_names.txt",
                     "channel_names.csv", "channels.txt", "markers.txt"]
        for alt_name in alt_names:
            alt_file = Path(meta_dir) / alt_name
            if alt_file.exists():
                channel_file = alt_file
                break
        else:
            return None

    print(f"Loading channel names from: {channel_file}")

    # Read all non-empty, non-comment lines
    lines = []
    with open(channel_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith('#'):
                lines.append(line)

    if not lines:
        print("  Warning: File is empty")
        return None

    channel_dict = {}

    # Detect format by checking first line
    first_line = lines[0]

    # Check if it's cycle-prefixed format
    is_cycle_prefixed = (
        ':' in first_line or
        '\t' in first_line or
        (first_line.split(',')[0].strip().isdigit() and len(first_line.split(',')) > 2)
    )

    if is_cycle_prefixed:
        # Cycle-prefixed format parsing
        for line_num, line in enumerate(lines, 1):
            try:
                if ':' in line:
                    cycle_str, names_str = line.split(':', 1)
                    cycle = int(cycle_str.strip())
                    names = [n.strip() for n in names_str.split(',')]
                elif '\t' in line:
                    parts = line.split('\t')
                    cycle = int(parts[0].strip())
                    names = [n.strip() for n in parts[1:]]
                else:
                    parts = line.split(',')
                    cycle = int(parts[0].strip())
                    names = [n.strip() for n in parts[1:]]

                channel_dict[cycle] = names

            except (ValueError, IndexError):
                print(f"  Warning: Could not parse line {line_num}: {line}")
                continue
    else:
        # Simple list format (one channel per line)
        # Detect cycles from DAPI-XX pattern
        current_cycle = 0
        cycle_channels = []

        for line in lines:
            # Check if this line indicates a new cycle (e.g., DAPI-01, DAPI-02)
            dapi_match = re.match(r'DAPI[-_]?(\d+)', line, re.IGNORECASE)

            if dapi_match:
                # Save previous cycle if we have channels
                if cycle_channels and current_cycle > 0:
                    channel_dict[current_cycle] = cycle_channels

                # Start new cycle
                current_cycle = int(dapi_match.group(1))
                cycle_channels = [line]  # DAPI is first channel
            elif current_cycle > 0:
                # Add channel to current cycle
                cycle_channels.append(line)

                # If we have enough channels, save and reset
                if len(cycle_channels) == channels_per_cycle:
                    channel_dict[current_cycle] = cycle_channels
                    cycle_channels = []

        # Save final cycle if incomplete
        if cycle_channels and current_cycle > 0:
            channel_dict[current_cycle] = cycle_channels

    if channel_dict:
        # Apply uniquification if requested (default)
        if make_unique:
            channel_dict = make_channel_names_unique(channel_dict)
            unique_suffix = " (with unique names)"
        else:
            unique_suffix = ""

        print(f"  Loaded {len(channel_dict)} cycle(s){unique_suffix}")
        for cyc_num in sorted(channel_dict.keys()):
            print(f"    Cycle {cyc_num}: {', '.join(channel_dict[cyc_num])}")

    return channel_dict


def save_channel_names_template(
    meta_dir: Union[str, Path],
    n_cycles: int = 13,
    n_channels: int = 4,
    filename: str = "channel_names.txt"
) -> Path:
    """
    Save a template channel names file to the metadata directory.

    Parameters
    ----------
    meta_dir : Path
        Path to metadata directory
    n_cycles : int
        Number of cycles to include in template
    n_channels : int
        Number of channels per cycle
    filename : str
        Output filename

    Returns
    -------
    Path
        Path to the created template file

    Example
    -------
    >>> save_channel_names_template(meta_dir, n_cycles=13, n_channels=4)
    """
    template_file = Path(meta_dir) / filename

    with open(template_file, 'w') as f:
        f.write("# Channel Names Configuration\n")
        f.write("# Format: cycle_number: CH1_name, CH2_name, CH3_name, CH4_name\n")
        f.write("# Lines starting with # are comments\n")
        f.write("#\n")
        f.write("# Example:\n")
        f.write("# 1: DAPI, CD31, CD8, Empty\n")
        f.write("# 2: DAPI, CD20, Ki67, CD3e\n")
        f.write("#\n")

        for cyc in range(1, n_cycles + 1):
            channels = ["DAPI"] + [f"CH{i}" for i in range(2, n_channels + 1)]
            f.write(f"{cyc}: {', '.join(channels)}\n")

    print(f"Template saved to: {template_file}")
    return template_file


def make_channel_names_unique(channel_dict: Dict[int, List[str]]) -> Dict[int, List[str]]:
    """
    Ensure all channel names are unique both within and across cycles.

    Naming convention for duplicates (Blank/Empty channels):
    - Format: {name}_{cycle}{channel_letter}
    - cycle: no zero padding (1, 2, ... 9, 10, 11)
    - channel_letter: a=CH2, b=CH3, c=CH4
    - Example: Blank on cycle 9, CH3 -> "Blank_9b"
    - Example: Empty on cycle 6, CH2 -> "Empty_6a"

    DAPI channels after cycle 1 get cycle suffix: DAPI_cyc02, DAPI_cyc03, etc.

    Parameters
    ----------
    channel_dict : dict
        Dictionary mapping cycle numbers (int) to channel name lists

    Returns
    -------
    dict
        Dictionary with unique channel names

    Example
    -------
    >>> channel_dict = {1: ['DAPI-01', 'Blank', 'Blank', 'Blank']}
    >>> unique = make_channel_names_unique(channel_dict)
    >>> print(unique[1])
    ['DAPI-01', 'Blank_1a', 'Blank_1b', 'Blank_1c']
    """
    # Channel index to letter mapping (0-indexed: CH2=0->a, CH3=1->b, CH4=2->c)
    def get_channel_letter(channel_idx):
        """Convert channel index (0-based, excluding CH1) to letter."""
        return chr(ord('a') + channel_idx)

    # First pass: identify all names that appear multiple times (within or across cycles)
    # These need the _Xa suffix format
    global_name_counts = {}
    for cycle_num, cycle_names in channel_dict.items():
        for name in cycle_names:
            if name.upper().startswith('DAPI'):
                continue
            key = name.upper()
            if key not in global_name_counts:
                global_name_counts[key] = 0
            global_name_counts[key] += 1

    # Names that need suffix: appear more than once globally
    needs_suffix = {name for name, count in global_name_counts.items() if count > 1}

    unique_dict = {}

    for cycle_num in sorted(channel_dict.keys()):
        cycle_names = channel_dict[cycle_num]
        new_names = []

        for ch_idx, name in enumerate(cycle_names):
            # DAPI handling: keep first cycle DAPI as-is, suffix others
            if name.upper().startswith('DAPI'):
                if cycle_num > 1:
                    new_names.append(f"DAPI_cyc{cycle_num:02d}")
                else:
                    new_names.append(name)
                continue

            # Check if this name needs a suffix (duplicated anywhere)
            if name.upper() in needs_suffix:
                # Use format: name_Xa where X=cycle, a=channel letter
                # ch_idx is 0-based: CH1=0, CH2=1, CH3=2, CH4=3
                # Letter: CH2->a, CH3->b, CH4->c (ch_idx - 1)
                letter = get_channel_letter(ch_idx - 1) if ch_idx > 0 else 'a'
                unique_name = f"{name}_{cycle_num}{letter}"
                new_names.append(unique_name)
            else:
                # Unique name - keep as-is
                new_names.append(name)

        unique_dict[cycle_num] = new_names

    return unique_dict


# =============================================================================
# OME-TIFF CHANNEL EXTRACTION
# =============================================================================

def extract_channels_from_ome(
    reg_dir: Union[str, Path],
    sig_dir: Union[str, Path],
    start_cycle: int,
    end_cycle: int,
    drop_duplicates: bool = True,
    drop_empty: bool = True,
    dup_name: str = 'DAPI',
    empty_name: str = 'Empty'
) -> List[str]:
    """
    Extract single-channel TIFFs from registered OME-TIFFs.

    Reads OME-XML metadata to get channel names and saves each channel
    as a separate TIFF file for downstream processing.

    Parameters
    ----------
    reg_dir : Path
        Directory containing registered OME-TIFF files
    sig_dir : Path
        Output directory for extracted channel TIFFs
    start_cycle : int
        First cycle to process
    end_cycle : int
        Last cycle to process
    drop_duplicates : bool
        Skip duplicate channels (e.g., DAPI after cycle 1) (default: True)
    drop_empty : bool
        Skip channels marked as 'Empty' (default: True)
    dup_name : str
        Name of duplicate channel to skip (default: 'DAPI')
    empty_name : str
        Name pattern for empty channels (default: 'Empty')

    Returns
    -------
    list
        List of extracted channel names

    Example
    -------
    >>> channels = extract_channels_from_ome(reg_dir, sig_dir, 1, 13)
    >>> print(f"Extracted {len(channels)} channels")
    """
    reg_dir = Path(reg_dir)
    sig_dir = Path(sig_dir)
    sig_dir.mkdir(parents=True, exist_ok=True)

    extracted = []

    for i in range(start_cycle, end_cycle + 1):
        ome_path = reg_dir / f"cyc{str(i).zfill(2)}.ome.tif"

        if not ome_path.exists():
            print(f"Warning: {ome_path} not found, skipping")
            continue

        ome_data = tiff.imread(str(ome_path))

        # Parse OME metadata for channel names
        with tiff.TiffFile(str(ome_path)) as tif:
            xml_desc = tif.ome_metadata

        root = ET.fromstring(xml_desc)
        namespace = {'ome': 'http://www.openmicroscopy.org/Schemas/OME/2016-06'}

        channel_names = []
        for channel in root.findall('.//ome:Channel', namespace):
            channel_names.append(channel.attrib.get('Name', f'Channel_{len(channel_names)}'))

        # Extract each channel
        for idx, name in enumerate(channel_names):
            # Skip duplicate DAPI after cycle 1
            if i > 1 and drop_duplicates and name == dup_name:
                continue

            # Skip empty channels
            if drop_empty and empty_name.lower() in name.lower():
                continue

            # Extract channel data
            if ome_data.ndim > 2:
                channel_data = ome_data[idx, ...]
            else:
                channel_data = ome_data

            # Save
            output_path = sig_dir / f"{name}.tif"
            tiff.imwrite(str(output_path), channel_data)
            extracted.append(name)
            print(f"Cycle {i}: Saved {name}")

    return extracted


# =============================================================================
# IMAGE RESIZE UTILITIES
# =============================================================================

# Resize parameters
RESIZE_ORDER = 1
RESIZE_MODE = "symmetric"
PRESERVE_RANGE = True


def resize_images_list(
    images_list: List[np.ndarray],
    side_size: Optional[float] = None,
    x_side_size: Optional[float] = None,
    y_side_size: Optional[float] = None
) -> List[np.ndarray]:
    """
    Resize a list of images to specified dimensions.

    Used for creating uniform-sized images for illumination correction
    estimation or other processing steps.

    Parameters
    ----------
    images_list : list
        List of 2D numpy arrays (images)
    side_size : float, optional
        Target size for both dimensions (square output)
    x_side_size : float, optional
        Target X dimension size
    y_side_size : float, optional
        Target Y dimension size

    Returns
    -------
    list
        List of resized images

    Example
    -------
    >>> resized = resize_images_list(images, side_size=128)
    >>> print(resized[0].shape)
    (128, 128)
    """
    if side_size is not None:
        y_side_size = x_side_size = side_size

    resized_images_list = []
    for i, im in enumerate(images_list):
        if im.shape[0] != x_side_size or im.shape[1] != y_side_size:
            resized_images_list.append(skresize(
                im,
                (int(x_side_size), int(y_side_size)),
                order=RESIZE_ORDER,
                mode=RESIZE_MODE,
                preserve_range=PRESERVE_RANGE
            ))
        else:
            resized_images_list.append(im)

    return resized_images_list


def resize_stack_for_estimation(
    stack: np.ndarray,
    working_size: int = 128
) -> np.ndarray:
    """
    Resize a stack of images for fast estimation algorithms.

    Converts from (n_images, height, width) to (height, width, n_images)
    format commonly used by illumination correction algorithms.

    Parameters
    ----------
    stack : np.ndarray
        Input stack with shape (n_images, height, width)
    working_size : int
        Target size for height and width (default: 128)

    Returns
    -------
    np.ndarray
        Resized stack with shape (working_size, working_size, n_images)

    Example
    -------
    >>> downsized = resize_stack_for_estimation(stack, working_size=128)
    >>> print(downsized.shape)
    (128, 128, N)
    """
    # Convert to list format
    images_list = [stack[i] for i in range(stack.shape[0])]

    # Resize
    resized_list = resize_images_list(images_list, side_size=working_size)

    # Stack back with transposed dimensions
    downsized = np.stack(resized_list, axis=2)

    return downsized


# =============================================================================
# PARALLEL PROCESSING UTILITIES
# =============================================================================

class ProgressCounter:
    """
    Thread-safe counter for progress tracking with reduced verbosity.

    Used for tracking progress of parallel processing operations.
    Only prints updates on first item, every Nth item, and last item
    to reduce output clutter in HPC environments.

    Parameters
    ----------
    total : int
        Total number of items to process
    desc : str
        Description prefix for progress messages
    print_every : int
        Print progress every N items (default: 4)

    Example
    -------
    >>> progress = ProgressCounter(100, desc="Processing tiles")
    >>> for i in range(100):
    ...     do_work()
    ...     progress.increment(f"Completed tile {i}")
    """

    def __init__(self, total: int, desc: str = "Processing", print_every: int = 4):
        self.count = 0
        self.total = total
        self.desc = desc
        self.lock = threading.Lock()
        self.start_time = datetime.now()
        self.print_every = print_every

    def increment(self, msg: str = "") -> None:
        """Increment counter and optionally print progress."""
        with self.lock:
            self.count += 1
            # Only print on first, every Nth, and last item
            should_print = (
                self.count == 1 or
                self.count == self.total or
                self.count % self.print_every == 0
            )
            if should_print:
                elapsed = (datetime.now() - self.start_time).total_seconds()
                rate = self.count / elapsed if elapsed > 0 else 0
                eta = (self.total - self.count) / rate if rate > 0 else 0
                print(f"[{self.count}/{self.total}] {self.desc} "
                      f"({elapsed:.0f}s elapsed, ~{eta:.0f}s remaining)")
                sys.stdout.flush()


def load_tiles_parallel(
    file_paths: List[str],
    max_workers: Optional[int] = None,
    default_workers: int = 32
) -> np.ndarray:
    """
    Load multiple tile images in parallel using ThreadPoolExecutor.

    Performance: 10-20x faster than sequential loading for 50+ images.
    Uses I/O-bound parallelism which scales well with many workers.

    Parameters
    ----------
    file_paths : list
        List of file paths to load
    max_workers : int, optional
        Maximum number of parallel workers. If None, uses default_workers.
    default_workers : int
        Default number of workers if max_workers not specified (default: 32)

    Returns
    -------
    np.ndarray
        Stacked array of images with shape (n_images, height, width)

    Example
    -------
    >>> files = glob("tiles/*.tif")
    >>> stack = load_tiles_parallel(files, max_workers=64)
    >>> print(stack.shape)
    (117, 2048, 2048)
    """
    if max_workers is None:
        max_workers = default_workers

    n_files = len(file_paths)
    if n_files == 0:
        return np.array([])

    workers = min(max_workers, n_files)

    with ThreadPoolExecutor(max_workers=workers) as executor:
        images = list(executor.map(imread, file_paths))

    return np.stack(images, axis=0)


def load_stack_parallel(
    file_paths: List[str],
    max_workers: Optional[int] = None,
    default_workers: int = 32
) -> np.ndarray:
    """
    Load z-stack images in parallel for EDF processing.

    Identical to load_tiles_parallel but named for clarity when
    loading z-stacks for extended depth of focus processing.

    Parameters
    ----------
    file_paths : list
        List of file paths to load (one per z-plane)
    max_workers : int, optional
        Maximum number of parallel workers
    default_workers : int
        Default number of workers (default: 32)

    Returns
    -------
    np.ndarray
        Stacked array with shape (n_zplanes, height, width)
    """
    return load_tiles_parallel(file_paths, max_workers, default_workers)


def alphanumeric_key(s: str) -> List:
    """
    Key function for natural alphanumeric sorting.

    Splits string into numeric and non-numeric parts for sorting
    so that "tile_2" comes before "tile_10".

    Parameters
    ----------
    s : str
        String to generate sort key for

    Returns
    -------
    list
        List of alternating string and integer parts

    Example
    -------
    >>> sorted(["tile_2", "tile_10", "tile_1"], key=alphanumeric_key)
    ['tile_1', 'tile_2', 'tile_10']
    """
    return [int(c) if c.isdigit() else c.lower() for c in re.split('([0-9]+)', s)]


def cleanup_gpu_memory(
    device_id: int = None,
    device_ids: List[int] = None,
    verbose: bool = False
) -> None:
    """
    Explicitly free GPU memory pools to prevent OOM errors.

    Should be called after processing each channel or batch
    to release memory for subsequent operations.

    Parameters
    ----------
    device_id : int, optional
        Single GPU device ID to clean up (for backwards compatibility)
    device_ids : list of int, optional
        List of GPU device IDs to clean up (preferred for multi-GPU)
    verbose : bool
        Print cleanup status (default: False)

    Examples
    --------
    >>> cleanup_gpu_memory(device_id=0)  # Single GPU
    >>> cleanup_gpu_memory(device_ids=[0, 1, 2, 3])  # Multi-GPU
    """
    # Determine which devices to clean
    if device_ids is not None:
        devices = device_ids
    elif device_id is not None:
        devices = [device_id]
    else:
        devices = [0]  # Default to device 0

    try:
        import cupy as cp
        for dev_id in devices:
            with cp.cuda.Device(dev_id):
                cp.get_default_memory_pool().free_all_blocks()
                cp.get_default_pinned_memory_pool().free_all_blocks()
                if verbose:
                    print(f"  GPU {dev_id}: memory pools cleared")
    except Exception as e:
        if verbose:
            print(f"  Warning: GPU cleanup failed: {e}")
    gc.collect()


# Thread-safe lock for zarr writes (zarr is not fully thread-safe for concurrent writes)
zarr_write_lock = threading.Lock()


@dataclass
class ProcessingConfig:
    """
    Configuration for parallel processing pipeline.

    Centralizes all processing parameters for BaSiC correction,
    stitching, and output formatting. Create an instance in your
    notebook and pass to processing functions.

    Example
    -------
    >>> config = ProcessingConfig(
    ...     n_rows=13, n_cols=9,
    ...     blend_mode='sigmoid', blend_sigma=10.0
    ... )
    >>> process_zplane(..., config=config)
    """
    # Grid configuration
    n_rows: int = 13
    n_cols: int = 9

    # BaSiC illumination correction parameters
    basic_if_darkfield: bool = True
    basic_max_iterations: int = 500
    basic_optimization_tolerance: float = 1e-6
    basic_max_reweight_iterations: int = 25
    basic_reweight_tolerance: float = 1e-3
    basic_flatfield_min: float = 0.1  # Minimum flatfield value to prevent division artifacts

    # Stitching parameters
    initial_ncc_threshold: float = 0.078
    overlap_percentage: float = 30.0
    pou: float = 0.5
    max_cores: int = 16

    # Blend mode
    blend_mode: str = 'sigmoid'  # 'sigmoid' or 'legacy'
    blend_sigma: float = 10.0

    # Output configuration
    output_format: str = 'tiff'  # 'tiff' or 'zarr'
    save_tiff_for_qc: bool = True

    # Parallel processing
    io_workers: int = 32
    zplanes_per_gpu: int = 6

    # GPU settings (populated at runtime)
    gpu_device_ids: List[int] = field(default_factory=lambda: [0])
