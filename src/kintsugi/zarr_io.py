"""
OME-Zarr I/O utilities for KINTSUGI pipeline.

This module provides a unified interface for reading and writing OME-Zarr format,
enabling chunked, parallel I/O across the entire pipeline.

OME-Zarr Benefits:
- Chunked storage enables parallel read/write
- No file size limits (unlike TIFF 4GB limit)
- Native multi-resolution pyramids
- Cloud-ready format
- Direct GPU memory access (with zarr 3.0+)

Data Structure:
    experiment.zarr/
    ├── .zattrs                    # OME-NGFF metadata
    ├── cyc01/                     # Per-cycle group
    │   ├── raw/                   # Raw tiles (t, z, y, x) per channel
    │   │   ├── CH1/
    │   │   ├── CH2/
    │   │   └── ...
    │   ├── corrected/             # After BaSiC correction
    │   │   └── ...
    │   ├── stitched/              # After stitching (z, c, y, x)
    │   ├── deconvolved/           # After deconvolution
    │   └── edf/                   # After EDF projection (c, y, x)
    ├── cyc02/
    │   └── ...
    └── registered/                # Final registered output
        ├── cyc01/
        └── ...
"""

import os
import json
import numpy as np
from pathlib import Path
from typing import Union, Optional, Tuple, List, Dict, Any, Callable
from datetime import datetime
import warnings

# Zarr imports
try:
    import zarr
    from zarr.storage import DirectoryStore
    ZARR_AVAILABLE = True
except ImportError:
    ZARR_AVAILABLE = False
    warnings.warn("zarr not installed. Install with: pip install zarr")

# OME-Zarr imports
try:
    from ome_zarr.io import parse_url
    from ome_zarr.writer import write_image, write_multiscale
    from ome_zarr.reader import Reader
    OME_ZARR_AVAILABLE = True
except ImportError:
    OME_ZARR_AVAILABLE = False
    warnings.warn("ome-zarr not installed. Install with: pip install ome-zarr")

# Dask for lazy/parallel operations
try:
    import dask.array as da
    from dask import delayed
    DASK_AVAILABLE = True
except ImportError:
    DASK_AVAILABLE = False
    warnings.warn("dask not installed. Install with: pip install dask[array]")


# Default chunking strategies for different data types
DEFAULT_CHUNKS = {
    'raw_tiles': (1, 1, 512, 512),      # (tile, z, y, x)
    'stitched': (1, 1, 1024, 1024),     # (z, c, y, x)
    'edf': (1, 2048, 2048),             # (c, y, x)
    'registered': (1, 2048, 2048),      # (c, y, x)
}

# Compression settings
DEFAULT_COMPRESSOR = zarr.Blosc(cname='zstd', clevel=3, shuffle=zarr.Blosc.BITSHUFFLE)


class KintsugiZarr:
    """
    Main class for managing OME-Zarr storage in KINTSUGI pipeline.

    This class provides a high-level interface for creating, reading, and
    writing OME-Zarr datasets throughout the processing pipeline.

    Parameters
    ----------
    path : str or Path
        Path to the zarr store (directory)
    mode : str
        'r' for read-only, 'w' for write (creates new), 'a' for append/modify

    Attributes
    ----------
    root : zarr.Group
        Root group of the zarr store
    path : Path
        Path to the zarr store

    Examples
    --------
    >>> # Create new zarr store for experiment
    >>> store = KintsugiZarr("experiment.zarr", mode="w")
    >>> store.create_cycle(1, n_channels=4, n_zplanes=17)
    >>>
    >>> # Write stitched data
    >>> store.write_stitched(cycle=1, channel=1, zplane=5, data=stitched_array)
    >>>
    >>> # Read data lazily with dask
    >>> lazy_stack = store.read_stitched_dask(cycle=1, channel=1)
    """

    def __init__(self, path: Union[str, Path], mode: str = 'a'):
        if not ZARR_AVAILABLE:
            raise ImportError("zarr is required. Install with: pip install zarr")

        self.path = Path(path)
        self.mode = mode

        if mode == 'w':
            # Create new store
            self.store = DirectoryStore(str(self.path))
            self.root = zarr.group(store=self.store, overwrite=True)
            self._init_metadata()
        elif mode == 'r':
            # Read-only
            if not self.path.exists():
                raise FileNotFoundError(f"Zarr store not found: {self.path}")
            self.store = DirectoryStore(str(self.path))
            self.root = zarr.open_group(store=self.store, mode='r')
        else:  # 'a' - append/modify
            self.store = DirectoryStore(str(self.path))
            if self.path.exists():
                self.root = zarr.open_group(store=self.store, mode='a')
            else:
                self.root = zarr.group(store=self.store)
                self._init_metadata()

    def _init_metadata(self):
        """Initialize OME-NGFF metadata."""
        self.root.attrs['kintsugi_version'] = '1.0'
        self.root.attrs['created'] = datetime.now().isoformat()
        self.root.attrs['ome_ngff_version'] = '0.4'

    def create_cycle(self,
                     cycle: int,
                     n_channels: int = 4,
                     n_zplanes: int = 17,
                     channel_names: Optional[List[str]] = None,
                     tile_shape: Optional[Tuple[int, int]] = None,
                     n_tiles: Optional[int] = None,
                     dtype: np.dtype = np.uint16) -> zarr.Group:
        """
        Create a cycle group with subgroups for each processing stage.

        Parameters
        ----------
        cycle : int
            Cycle number (1-indexed)
        n_channels : int
            Number of channels
        n_zplanes : int
            Number of z-planes
        channel_names : list, optional
            Names for each channel (e.g., ['DAPI', 'CD31', 'CD8', 'Empty'])
        tile_shape : tuple, optional
            (height, width) of each tile
        n_tiles : int, optional
            Number of tiles
        dtype : dtype
            Data type for arrays

        Returns
        -------
        zarr.Group
            The created cycle group
        """
        cycle_name = f"cyc{cycle:02d}"

        if cycle_name in self.root:
            return self.root[cycle_name]

        cycle_group = self.root.create_group(cycle_name)

        # Store cycle metadata
        cycle_group.attrs['n_channels'] = n_channels
        cycle_group.attrs['n_zplanes'] = n_zplanes
        cycle_group.attrs['dtype'] = str(dtype)

        if channel_names:
            cycle_group.attrs['channel_names'] = channel_names

        # Create stage subgroups
        for stage in ['raw', 'corrected', 'stitched', 'deconvolved', 'edf']:
            cycle_group.create_group(stage)

        return cycle_group

    def write_raw_tiles(self,
                        cycle: int,
                        channel: int,
                        zplane: int,
                        tiles: np.ndarray,
                        tile_positions: Optional[np.ndarray] = None) -> zarr.Array:
        """
        Write raw tiles for a specific cycle/channel/z-plane.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number (1-indexed)
        zplane : int
            Z-plane number (1-indexed)
        tiles : np.ndarray
            Array of tiles with shape (n_tiles, height, width)
        tile_positions : np.ndarray, optional
            Tile positions array with shape (n_tiles, 2) for (row, col)

        Returns
        -------
        zarr.Array
            The created array
        """
        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"
        zplane_name = f"Z{zplane:02d}"

        # Ensure cycle exists
        if cycle_name not in self.root:
            self.create_cycle(cycle)

        raw_group = self.root[cycle_name]['raw']

        # Create channel group if needed
        if channel_name not in raw_group:
            raw_group.create_group(channel_name)

        ch_group = raw_group[channel_name]

        # Write tiles
        arr = ch_group.create_dataset(
            zplane_name,
            data=tiles,
            chunks=(1, tiles.shape[1], tiles.shape[2]),
            compressor=DEFAULT_COMPRESSOR,
            overwrite=True
        )

        if tile_positions is not None:
            arr.attrs['tile_positions'] = tile_positions.tolist()

        return arr

    def write_stitched(self,
                       cycle: int,
                       channel: int,
                       zplane: int,
                       data: np.ndarray,
                       chunks: Optional[Tuple[int, ...]] = None) -> zarr.Array:
        """
        Write stitched image for a specific cycle/channel/z-plane.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number (1-indexed)
        zplane : int
            Z-plane number (1-indexed)
        data : np.ndarray
            2D stitched image
        chunks : tuple, optional
            Chunk size for storage

        Returns
        -------
        zarr.Array
            The created array
        """
        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"

        if cycle_name not in self.root:
            self.create_cycle(cycle)

        stitched_group = self.root[cycle_name]['stitched']

        if channel_name not in stitched_group:
            stitched_group.create_group(channel_name)

        ch_group = stitched_group[channel_name]

        if chunks is None:
            chunks = (min(1024, data.shape[0]), min(1024, data.shape[1]))

        arr = ch_group.create_dataset(
            f"Z{zplane:02d}",
            data=data,
            chunks=chunks,
            compressor=DEFAULT_COMPRESSOR,
            overwrite=True
        )

        return arr

    def read_stitched(self,
                      cycle: int,
                      channel: int,
                      zplane: Optional[int] = None) -> Union[np.ndarray, zarr.Array]:
        """
        Read stitched data.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number
        zplane : int, optional
            Specific z-plane, or None for all

        Returns
        -------
        np.ndarray or zarr.Array
            The requested data
        """
        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"

        ch_group = self.root[cycle_name]['stitched'][channel_name]

        if zplane is not None:
            return ch_group[f"Z{zplane:02d}"][:]
        else:
            # Return all z-planes as 3D array
            zplanes = sorted([k for k in ch_group.keys() if k.startswith('Z')])
            stack = np.stack([ch_group[z][:] for z in zplanes], axis=0)
            return stack

    def read_stitched_dask(self,
                           cycle: int,
                           channel: int) -> 'da.Array':
        """
        Read stitched data as a lazy dask array.

        This enables out-of-core processing and parallel computation.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number

        Returns
        -------
        dask.array.Array
            Lazy dask array for the z-stack
        """
        if not DASK_AVAILABLE:
            raise ImportError("dask is required. Install with: pip install dask[array]")

        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"

        ch_group = self.root[cycle_name]['stitched'][channel_name]
        zplanes = sorted([k for k in ch_group.keys() if k.startswith('Z')])

        # Create lazy dask arrays for each z-plane
        arrays = []
        for z in zplanes:
            zarr_arr = ch_group[z]
            dask_arr = da.from_zarr(zarr_arr)
            arrays.append(dask_arr[np.newaxis, ...])  # Add z dimension

        # Stack along z axis
        return da.concatenate(arrays, axis=0)

    def write_deconvolved(self,
                          cycle: int,
                          channel: int,
                          data: np.ndarray,
                          chunks: Optional[Tuple[int, ...]] = None) -> zarr.Array:
        """
        Write deconvolved z-stack.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number
        data : np.ndarray
            3D deconvolved stack (z, y, x)
        chunks : tuple, optional
            Chunk size

        Returns
        -------
        zarr.Array
            The created array
        """
        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"

        if cycle_name not in self.root:
            self.create_cycle(cycle)

        decon_group = self.root[cycle_name]['deconvolved']

        if chunks is None:
            chunks = (1, min(1024, data.shape[1]), min(1024, data.shape[2]))

        arr = decon_group.create_dataset(
            channel_name,
            data=data,
            chunks=chunks,
            compressor=DEFAULT_COMPRESSOR,
            overwrite=True
        )

        return arr

    def read_deconvolved_dask(self,
                              cycle: int,
                              channel: int) -> 'da.Array':
        """
        Read deconvolved stack as lazy dask array.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number

        Returns
        -------
        dask.array.Array
            Lazy dask array
        """
        if not DASK_AVAILABLE:
            raise ImportError("dask is required")

        cycle_name = f"cyc{cycle:02d}"
        channel_name = f"CH{channel}"

        zarr_arr = self.root[cycle_name]['deconvolved'][channel_name]
        return da.from_zarr(zarr_arr)

    def write_edf(self,
                  cycle: int,
                  channel: int,
                  data: np.ndarray,
                  channel_name: Optional[str] = None) -> zarr.Array:
        """
        Write EDF (2D projection) result.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number
        data : np.ndarray
            2D EDF image
        channel_name : str, optional
            Human-readable channel name (e.g., 'DAPI', 'CD31')

        Returns
        -------
        zarr.Array
            The created array
        """
        cycle_name = f"cyc{cycle:02d}"
        ch_key = f"CH{channel}"

        if cycle_name not in self.root:
            self.create_cycle(cycle)

        edf_group = self.root[cycle_name]['edf']

        chunks = (min(2048, data.shape[0]), min(2048, data.shape[1]))

        arr = edf_group.create_dataset(
            ch_key,
            data=data,
            chunks=chunks,
            compressor=DEFAULT_COMPRESSOR,
            overwrite=True
        )

        if channel_name:
            arr.attrs['channel_name'] = channel_name

        return arr

    def read_edf(self, cycle: int, channel: int) -> np.ndarray:
        """Read EDF result."""
        cycle_name = f"cyc{cycle:02d}"
        return self.root[cycle_name]['edf'][f"CH{channel}"][:]

    def get_all_edf_dask(self, cycle: int) -> 'da.Array':
        """
        Get all EDF channels for a cycle as a lazy dask array.

        Parameters
        ----------
        cycle : int
            Cycle number

        Returns
        -------
        dask.array.Array
            Array with shape (n_channels, height, width)
        """
        if not DASK_AVAILABLE:
            raise ImportError("dask is required")

        cycle_name = f"cyc{cycle:02d}"
        edf_group = self.root[cycle_name]['edf']

        channels = sorted([k for k in edf_group.keys() if k.startswith('CH')])

        arrays = []
        for ch in channels:
            dask_arr = da.from_zarr(edf_group[ch])
            arrays.append(dask_arr[np.newaxis, ...])

        return da.concatenate(arrays, axis=0)

    def create_registered_group(self,
                                pixel_size_um: float = 0.377,
                                channel_names_per_cycle: Optional[Dict[int, List[str]]] = None):
        """
        Create the registered output group with OME-NGFF metadata.

        Parameters
        ----------
        pixel_size_um : float
            Pixel size in micrometers
        channel_names_per_cycle : dict, optional
            Dict mapping cycle number to list of channel names
        """
        if 'registered' not in self.root:
            reg_group = self.root.create_group('registered')
        else:
            reg_group = self.root['registered']

        reg_group.attrs['pixel_size_um'] = pixel_size_um

        if channel_names_per_cycle:
            reg_group.attrs['channel_names'] = channel_names_per_cycle

        return reg_group

    def write_registered(self,
                         cycle: int,
                         data: np.ndarray,
                         channel_names: Optional[List[str]] = None,
                         pixel_size_um: float = 0.377) -> zarr.Array:
        """
        Write registered multi-channel image.

        Parameters
        ----------
        cycle : int
            Cycle number
        data : np.ndarray
            Multi-channel image (c, y, x)
        channel_names : list, optional
            Names for each channel
        pixel_size_um : float
            Pixel size in micrometers

        Returns
        -------
        zarr.Array
            The created array
        """
        if 'registered' not in self.root:
            self.create_registered_group(pixel_size_um)

        reg_group = self.root['registered']
        cycle_name = f"cyc{cycle:02d}"

        chunks = (1, min(2048, data.shape[1]), min(2048, data.shape[2]))

        arr = reg_group.create_dataset(
            cycle_name,
            data=data,
            chunks=chunks,
            compressor=DEFAULT_COMPRESSOR,
            overwrite=True
        )

        if channel_names:
            arr.attrs['channel_names'] = channel_names

        return arr

    def write_ome_zarr(self,
                       cycle: int,
                       data: np.ndarray,
                       axes: str = "cyx",
                       channel_names: Optional[List[str]] = None,
                       pixel_size_um: float = 0.377,
                       generate_pyramid: bool = True,
                       pyramid_levels: int = 4) -> None:
        """
        Write data in full OME-NGFF format with proper metadata.

        Uses ome-zarr-py for compliant OME-NGFF output.

        Parameters
        ----------
        cycle : int
            Cycle number
        data : np.ndarray
            Image data
        axes : str
            Axis order (e.g., "cyx", "zyx", "czyx")
        channel_names : list, optional
            Channel names
        pixel_size_um : float
            Pixel size in micrometers
        generate_pyramid : bool
            Whether to generate multi-resolution pyramid
        pyramid_levels : int
            Number of pyramid levels
        """
        if not OME_ZARR_AVAILABLE:
            raise ImportError("ome-zarr is required. Install with: pip install ome-zarr")

        cycle_name = f"cyc{cycle:02d}"
        output_path = self.path / 'ome_ngff' / f"{cycle_name}.zarr"
        output_path.parent.mkdir(parents=True, exist_ok=True)

        store = parse_url(str(output_path), mode="w").store

        # Build coordinate transformations
        transforms = [
            [{"type": "scale", "scale": [pixel_size_um] * len(axes)}]
        ]

        # Generate pyramid if requested
        if generate_pyramid:
            pyramid = [data]
            for level in range(1, pyramid_levels):
                factor = 2 ** level
                if 'z' in axes:
                    # Don't downsample z
                    z_idx = axes.index('z')
                    slices = [slice(None, None, factor if i != z_idx else 1)
                             for i in range(len(axes))]
                else:
                    slices = [slice(None, None, factor) for _ in range(len(axes))]
                pyramid.append(data[tuple(slices)])
                transforms.append([{"type": "scale",
                                   "scale": [pixel_size_um * factor] * len(axes)}])

            write_multiscale(
                pyramid=pyramid,
                group=zarr.group(store),
                axes=list(axes),
                coordinate_transformations=transforms,
                storage_options=dict(chunks=DEFAULT_CHUNKS.get('registered', (1, 2048, 2048))),
                name=cycle_name
            )
        else:
            write_image(
                image=data,
                group=zarr.group(store),
                axes=list(axes),
                coordinate_transformations=transforms[0],
                storage_options=dict(chunks=DEFAULT_CHUNKS.get('registered', (1, 2048, 2048)))
            )

    def close(self):
        """Close the zarr store."""
        if hasattr(self, 'store'):
            self.store.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False

    def list_cycles(self) -> List[int]:
        """List all cycle numbers in the store."""
        cycles = []
        for key in self.root.keys():
            if key.startswith('cyc'):
                try:
                    cycles.append(int(key[3:]))
                except ValueError:
                    pass
        return sorted(cycles)

    def get_cycle_info(self, cycle: int) -> Dict[str, Any]:
        """Get metadata about a cycle."""
        cycle_name = f"cyc{cycle:02d}"
        if cycle_name not in self.root:
            raise KeyError(f"Cycle {cycle} not found")

        cycle_group = self.root[cycle_name]

        info = dict(cycle_group.attrs)
        info['stages'] = list(cycle_group.keys())

        # Count data in each stage
        for stage in info['stages']:
            stage_group = cycle_group[stage]
            if isinstance(stage_group, zarr.Group):
                info[f'{stage}_contents'] = list(stage_group.keys())

        return info


# Convenience functions for common operations

def convert_tiff_stack_to_zarr(tiff_dir: Union[str, Path],
                                zarr_path: Union[str, Path],
                                cycle: int,
                                channel: int,
                                stage: str = 'stitched',
                                chunks: Optional[Tuple[int, ...]] = None) -> None:
    """
    Convert a directory of TIFF files to zarr format.

    Parameters
    ----------
    tiff_dir : str or Path
        Directory containing TIFF files
    zarr_path : str or Path
        Output zarr store path
    cycle : int
        Cycle number
    channel : int
        Channel number
    stage : str
        Processing stage ('stitched', 'deconvolved', etc.)
    chunks : tuple, optional
        Chunk size for zarr arrays
    """
    from glob import glob
    from skimage.io import imread

    tiff_dir = Path(tiff_dir)
    tiff_files = sorted(glob(str(tiff_dir / "*.tif")) + glob(str(tiff_dir / "*.tiff")))

    if not tiff_files:
        raise ValueError(f"No TIFF files found in {tiff_dir}")

    store = KintsugiZarr(zarr_path, mode='a')

    for i, tiff_file in enumerate(tiff_files):
        img = imread(tiff_file)
        zplane = i + 1

        if stage == 'stitched':
            store.write_stitched(cycle, channel, zplane, img, chunks)
        # Add other stages as needed

    store.close()


def zarr_to_dask_stack(zarr_path: Union[str, Path],
                       cycle: int,
                       channel: int,
                       stage: str = 'stitched') -> 'da.Array':
    """
    Load zarr data as a lazy dask array.

    Parameters
    ----------
    zarr_path : str or Path
        Path to zarr store
    cycle : int
        Cycle number
    channel : int
        Channel number
    stage : str
        Processing stage

    Returns
    -------
    dask.array.Array
        Lazy dask array
    """
    store = KintsugiZarr(zarr_path, mode='r')

    if stage == 'stitched':
        return store.read_stitched_dask(cycle, channel)
    elif stage == 'deconvolved':
        return store.read_deconvolved_dask(cycle, channel)
    else:
        raise ValueError(f"Unknown stage: {stage}")


def process_zarr_parallel(zarr_path: Union[str, Path],
                          cycle: int,
                          channel: int,
                          process_func: Callable,
                          output_stage: str,
                          input_stage: str = 'stitched',
                          **process_kwargs) -> None:
    """
    Process zarr data in parallel using dask.

    Parameters
    ----------
    zarr_path : str or Path
        Path to zarr store
    cycle : int
        Cycle number
    channel : int
        Channel number
    process_func : callable
        Function to apply to each chunk/block
    output_stage : str
        Stage to write output to
    input_stage : str
        Stage to read input from
    **process_kwargs
        Additional arguments passed to process_func
    """
    if not DASK_AVAILABLE:
        raise ImportError("dask is required for parallel processing")

    store = KintsugiZarr(zarr_path, mode='a')

    # Get input as dask array
    if input_stage == 'stitched':
        data = store.read_stitched_dask(cycle, channel)
    elif input_stage == 'deconvolved':
        data = store.read_deconvolved_dask(cycle, channel)
    else:
        raise ValueError(f"Unknown input stage: {input_stage}")

    # Apply function to blocks
    result = data.map_blocks(
        lambda block: process_func(block, **process_kwargs),
        dtype=data.dtype
    )

    # Compute and store result
    result_computed = result.compute()

    if output_stage == 'deconvolved':
        store.write_deconvolved(cycle, channel, result_computed)
    elif output_stage == 'edf':
        store.write_edf(cycle, channel, result_computed)

    store.close()


# =============================================================================
# Unified Output Interface
# =============================================================================

class OutputWriter:
    """
    Unified output writer supporting both TIFF and OME-Zarr formats.

    This class provides a common interface for writing pipeline outputs,
    allowing easy switching between formats without changing processing code.

    TIFF is recommended for:
    - Speed (fastest I/O)
    - Compatibility with existing tools
    - Small to medium datasets

    OME-Zarr is recommended for:
    - Large datasets (>4GB per image)
    - Cloud storage
    - Parallel/chunked access
    - Multi-resolution pyramids

    Parameters
    ----------
    output_dir : str or Path
        Base output directory
    format : str
        Output format: 'tiff' (default) or 'zarr'
    zarr_path : str or Path, optional
        Path for zarr store (required if format='zarr')
    compression : str, optional
        TIFF compression: None, 'zlib', 'lzw' (default: None for speed)

    Examples
    --------
    >>> # TIFF output (fast, compatible)
    >>> writer = OutputWriter("/output/dir", format='tiff')
    >>> writer.write_stitched(cycle=1, channel=2, zplane=5, data=img)
    >>>
    >>> # OME-Zarr output (scalable, cloud-ready)
    >>> writer = OutputWriter("/output/dir", format='zarr', zarr_path="experiment.zarr")
    >>> writer.write_stitched(cycle=1, channel=2, zplane=5, data=img)
    """

    def __init__(self,
                 output_dir: Union[str, Path],
                 format: str = 'tiff',
                 zarr_path: Optional[Union[str, Path]] = None,
                 compression: Optional[str] = None):
        self.output_dir = Path(output_dir)
        self.format = format.lower()
        self.compression = compression

        if self.format not in ('tiff', 'zarr'):
            raise ValueError(f"Unknown format: {format}. Use 'tiff' or 'zarr'")

        if self.format == 'zarr':
            if zarr_path is None:
                zarr_path = self.output_dir / "output.zarr"
            self.zarr_store = KintsugiZarr(zarr_path, mode='a')
        else:
            self.zarr_store = None

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def write_stitched(self,
                       cycle: int,
                       channel: int,
                       zplane: int,
                       data: np.ndarray,
                       channel_name: Optional[str] = None) -> str:
        """
        Write stitched image.

        Parameters
        ----------
        cycle : int
            Cycle number
        channel : int
            Channel number
        zplane : int
            Z-plane number
        data : np.ndarray
            2D stitched image
        channel_name : str, optional
            Human-readable channel name

        Returns
        -------
        str
            Path to the written file/array
        """
        if self.format == 'tiff':
            return self._write_tiff_stitched(cycle, channel, zplane, data)
        else:
            self.zarr_store.write_stitched(cycle, channel, zplane, data)
            return str(self.zarr_store.path)

    def _write_tiff_stitched(self,
                              cycle: int,
                              channel: int,
                              zplane: int,
                              data: np.ndarray) -> str:
        """Write stitched image as TIFF."""
        try:
            import tifffile
            use_tifffile = True
        except ImportError:
            from skimage.io import imsave
            use_tifffile = False

        # Create directory structure
        dest_dir = self.output_dir / f"cyc{cycle:02d}" / f"CH{channel}"
        dest_dir.mkdir(parents=True, exist_ok=True)

        output_path = dest_dir / f"{zplane:02d}.tif"

        if use_tifffile:
            # tifffile is faster and supports more compression options
            tifffile.imwrite(
                str(output_path),
                data,
                compression=self.compression,
                photometric='minisblack'
            )
        else:
            imsave(str(output_path), data, check_contrast=False)

        return str(output_path)

    def write_deconvolved(self,
                          cycle: int,
                          channel: int,
                          zplane: int,
                          data: np.ndarray) -> str:
        """Write deconvolved image."""
        if self.format == 'tiff':
            return self._write_tiff_deconvolved(cycle, channel, zplane, data)
        else:
            # For zarr, we typically write the whole stack at once
            # This method writes individual z-planes for compatibility
            cycle_name = f"cyc{cycle:02d}"
            ch_name = f"CH{channel}"

            if cycle_name not in self.zarr_store.root:
                self.zarr_store.create_cycle(cycle)

            decon_group = self.zarr_store.root[cycle_name]['deconvolved']
            if ch_name not in decon_group:
                decon_group.create_group(ch_name)

            ch_group = decon_group[ch_name]
            ch_group.create_dataset(
                f"Z{zplane:02d}",
                data=data,
                chunks=(min(1024, data.shape[0]), min(1024, data.shape[1])),
                compressor=DEFAULT_COMPRESSOR,
                overwrite=True
            )
            return str(self.zarr_store.path)

    def _write_tiff_deconvolved(self,
                                 cycle: int,
                                 channel: int,
                                 zplane: int,
                                 data: np.ndarray) -> str:
        """Write deconvolved image as TIFF."""
        try:
            import tifffile
            use_tifffile = True
        except ImportError:
            from skimage.io import imsave
            use_tifffile = False

        dest_dir = self.output_dir / f"cyc{cycle:02d}" / f"CH{channel}" / "deconvolved"
        dest_dir.mkdir(parents=True, exist_ok=True)

        output_path = dest_dir / f"{zplane:02d}.tif"

        if use_tifffile:
            tifffile.imwrite(str(output_path), data, compression=self.compression)
        else:
            imsave(str(output_path), data, check_contrast=False)

        return str(output_path)

    def write_edf(self,
                  cycle: int,
                  channel: int,
                  data: np.ndarray,
                  channel_name: Optional[str] = None) -> str:
        """Write EDF (extended depth of focus) result."""
        if self.format == 'tiff':
            return self._write_tiff_edf(cycle, channel, data, channel_name)
        else:
            self.zarr_store.write_edf(cycle, channel, data, channel_name)
            return str(self.zarr_store.path)

    def _write_tiff_edf(self,
                        cycle: int,
                        channel: int,
                        data: np.ndarray,
                        channel_name: Optional[str] = None) -> str:
        """Write EDF image as TIFF."""
        try:
            import tifffile
            use_tifffile = True
        except ImportError:
            from skimage.io import imsave
            use_tifffile = False

        dest_dir = self.output_dir / f"cyc{cycle:02d}"
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Use channel name if provided, otherwise use CH number
        filename = f"{channel_name}.tif" if channel_name else f"CH{channel}.tif"
        output_path = dest_dir / filename

        if use_tifffile:
            tifffile.imwrite(str(output_path), data, compression=self.compression)
        else:
            imsave(str(output_path), data, check_contrast=False)

        return str(output_path)

    def close(self):
        """Close any open resources."""
        if self.zarr_store is not None:
            self.zarr_store.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()
        return False


def convert_tiff_to_ome_zarr(tiff_base_dir: Union[str, Path],
                              zarr_path: Union[str, Path],
                              cycles: List[int],
                              channels: List[int],
                              stage: str = 'stitched',
                              channel_names: Optional[Dict[int, List[str]]] = None,
                              pixel_size_um: float = 0.377,
                              generate_pyramid: bool = True,
                              pyramid_levels: int = 4,
                              n_workers: int = 4) -> None:
    """
    Convert TIFF output directory to OME-Zarr format.

    This is useful for post-processing conversion after fast TIFF-based
    processing is complete.

    Parameters
    ----------
    tiff_base_dir : str or Path
        Base directory containing TIFF output (e.g., *_BaSiC_Stitched)
    zarr_path : str or Path
        Output OME-Zarr path
    cycles : list of int
        Cycle numbers to convert
    channels : list of int
        Channel numbers to convert
    stage : str
        Processing stage subdirectory name
    channel_names : dict, optional
        Mapping of cycle to list of channel names
    pixel_size_um : float
        Pixel size in micrometers
    generate_pyramid : bool
        Generate multi-resolution pyramid
    pyramid_levels : int
        Number of pyramid levels
    n_workers : int
        Number of parallel workers for loading TIFFs
    """
    from glob import glob
    from concurrent.futures import ThreadPoolExecutor

    try:
        import tifffile
        def read_tiff(path):
            return tifffile.imread(path)
    except ImportError:
        from skimage.io import imread
        def read_tiff(path):
            return imread(path)

    tiff_base_dir = Path(tiff_base_dir)
    store = KintsugiZarr(zarr_path, mode='w')

    print(f"Converting TIFF to OME-Zarr: {zarr_path}")

    for cycle in cycles:
        print(f"  Processing cycle {cycle}...")
        store.create_cycle(
            cycle,
            n_channels=len(channels),
            channel_names=channel_names.get(cycle) if channel_names else None
        )

        for channel in channels:
            ch_dir = tiff_base_dir / f"cyc{cycle:02d}" / f"CH{channel}"

            if stage == 'deconvolved':
                ch_dir = ch_dir / "deconvolved"

            tiff_files = sorted(glob(str(ch_dir / "*.tif")))

            if not tiff_files:
                print(f"    Warning: No TIFF files found in {ch_dir}")
                continue

            # Load TIFFs in parallel
            with ThreadPoolExecutor(max_workers=n_workers) as executor:
                images = list(executor.map(read_tiff, tiff_files))

            # Stack into 3D array
            stack = np.stack(images, axis=0)

            # Write to zarr
            if stage == 'stitched':
                for z, img in enumerate(images, start=1):
                    store.write_stitched(cycle, channel, z, img)
            elif stage == 'deconvolved':
                store.write_deconvolved(cycle, channel, stack)

            print(f"    CH{channel}: {len(images)} z-planes")

            del images, stack

    # Generate OME-NGFF compliant output with pyramids
    if generate_pyramid:
        print("  Generating multi-resolution pyramids...")
        for cycle in cycles:
            # Collect all channels for this cycle
            edf_group = store.root[f"cyc{cycle:02d}"].get('edf')
            if edf_group and len(edf_group.keys()) > 0:
                # Stack EDF channels
                ch_keys = sorted([k for k in edf_group.keys() if k.startswith('CH')])
                if ch_keys:
                    stack = np.stack([edf_group[k][:] for k in ch_keys], axis=0)
                    ch_names = channel_names.get(cycle) if channel_names else None
                    store.write_ome_zarr(
                        cycle, stack, axes="cyx",
                        channel_names=ch_names,
                        pixel_size_um=pixel_size_um,
                        generate_pyramid=True,
                        pyramid_levels=pyramid_levels
                    )

    store.close()
    print(f"Conversion complete: {zarr_path}")
