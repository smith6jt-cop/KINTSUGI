"""
KINTSUGI I/O and Channel Name Functions

This module contains functions for:
- Loading and saving channel name configurations
- Extracting channels from OME-TIFF files
- Image resize utilities for processing

Location: notebooks/Kio.py
"""

import os
import re
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

import numpy as np
import tifffile as tiff
from skimage.transform import resize as skresize


# =============================================================================
# CHANNEL NAME MANAGEMENT
# =============================================================================

def load_channel_names(
    meta_dir: Union[str, Path],
    filename: str = "CHANNELNAMES.txt",
    channels_per_cycle: int = 4
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

    Returns
    -------
    dict or None
        Dictionary mapping cycle number (int) to list of channel names,
        or None if file not found

    Example
    -------
    >>> channel_dict = load_channel_names(meta_dir)
    >>> print(channel_dict[1])
    ['DAPI-01', 'Blank', 'Blank', 'Blank']
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
        print(f"  Loaded {len(channel_dict)} cycle(s)")
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
    Ensure all channel names across cycles are unique.

    Appends cycle number suffix to duplicate names (except DAPI in cycle 1).

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
    >>> unique_names = make_channel_names_unique(channel_dict)
    """
    seen = {}
    unique_dict = {}

    for cycle_num in sorted(channel_dict.keys()):
        new_names = []

        for name in channel_dict[cycle_num]:
            # Skip DAPI after cycle 1 (usually kept as reference)
            if name.upper().startswith('DAPI') and cycle_num > 1:
                new_names.append(f"DAPI_cyc{cycle_num:02d}")
                continue

            if name in seen:
                # Make unique by appending cycle number
                unique_name = f"{name}_cyc{cycle_num:02d}"
                new_names.append(unique_name)
            else:
                seen[name] = cycle_num
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
