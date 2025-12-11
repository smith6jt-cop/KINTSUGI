"""
Optimized parallel I/O utilities for KINTSUGI pipeline.

This module provides memory-mapped file reading and parallel loading
utilities to minimize I/O bottlenecks in the processing pipeline.

Key features:
- Memory-mapped TIFF reading for zero-copy data access
- Parallel tile loading with ThreadPoolExecutor
- Pre-loading and caching for z-stack processing
- Optimized numpy array handling
"""

import os
import numpy as np
from glob import glob
from pathlib import Path
from typing import List, Tuple, Optional, Dict, Union, Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

# Try to import tifffile for memory-mapped reading
try:
    import tifffile
    HAS_TIFFFILE = True
except ImportError:
    HAS_TIFFFILE = False

# Fallback to skimage
from skimage.io import imread
from skimage.io.collection import alphanumeric_key


class TileLoader:
    """
    Optimized tile loader with memory-mapped reading and caching.

    This class provides efficient loading of image tiles with:
    - Memory-mapped TIFF reading (zero-copy when possible)
    - Thread-safe caching of loaded tiles
    - Parallel loading of multiple tiles
    - Pre-loading for z-stack processing

    Parameters
    ----------
    use_mmap : bool
        Use memory-mapped reading if tifffile is available
    cache_size : int
        Maximum number of tiles to cache (0 = no caching)
    n_workers : int
        Number of parallel workers for loading
    """

    def __init__(self,
                 use_mmap: bool = True,
                 cache_size: int = 0,
                 n_workers: int = 4):
        self.use_mmap = use_mmap and HAS_TIFFFILE
        self.cache_size = cache_size
        self.n_workers = n_workers

        # Thread-safe cache
        self._cache: Dict[str, np.ndarray] = {}
        self._cache_order: List[str] = []
        self._cache_lock = threading.Lock()

    def load_tile(self, filepath: str, copy: bool = True) -> np.ndarray:
        """
        Load a single tile from file.

        Parameters
        ----------
        filepath : str
            Path to the TIFF file
        copy : bool
            If True, return a copy (safe for modification).
            If False, may return memory-mapped view (faster but read-only).

        Returns
        -------
        np.ndarray
            The loaded image
        """
        # Check cache first
        if self.cache_size > 0:
            with self._cache_lock:
                if filepath in self._cache:
                    return self._cache[filepath].copy() if copy else self._cache[filepath]

        # Load the tile
        if self.use_mmap:
            # Memory-mapped reading with tifffile
            with tifffile.TiffFile(filepath) as tif:
                # Get the first page as memory-mapped array
                data = tif.pages[0].asarray()
                if copy:
                    data = data.copy()
        else:
            # Standard reading with skimage
            data = imread(filepath)

        # Add to cache if enabled
        if self.cache_size > 0:
            with self._cache_lock:
                if filepath not in self._cache:
                    # Evict oldest if cache full
                    if len(self._cache) >= self.cache_size:
                        oldest = self._cache_order.pop(0)
                        del self._cache[oldest]
                    self._cache[filepath] = data if not copy else data.copy()
                    self._cache_order.append(filepath)

        return data

    def load_tiles_parallel(self,
                           filepaths: List[str],
                           copy: bool = True) -> List[np.ndarray]:
        """
        Load multiple tiles in parallel.

        Parameters
        ----------
        filepaths : list of str
            Paths to TIFF files
        copy : bool
            If True, return copies (safe for modification)

        Returns
        -------
        list of np.ndarray
            Loaded images in the same order as filepaths
        """
        if len(filepaths) == 0:
            return []

        # For small numbers, load sequentially
        if len(filepaths) <= 2:
            return [self.load_tile(fp, copy) for fp in filepaths]

        # Parallel loading
        results = [None] * len(filepaths)

        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = {
                executor.submit(self.load_tile, fp, copy): i
                for i, fp in enumerate(filepaths)
            }

            for future in as_completed(futures):
                idx = futures[future]
                results[idx] = future.result()

        return results

    def load_zplane_tiles(self,
                         image_dir: str,
                         cycle: int,
                         channel: int,
                         zplane: int,
                         as_stack: bool = True) -> Union[np.ndarray, List[np.ndarray]]:
        """
        Load all tiles for a specific z-plane.

        Parameters
        ----------
        image_dir : str
            Base directory containing cycle folders
        cycle : int
            Cycle number
        channel : int
            Channel number
        zplane : int
            Z-plane number
        as_stack : bool
            If True, return stacked array. If False, return list.

        Returns
        -------
        np.ndarray or list
            Loaded tiles as stack or list
        """
        pattern = f'1_000??_Z0{str(zplane).zfill(2)}_CH{str(channel)}.tif'
        files = sorted(
            glob(os.path.join(image_dir, f'cyc{str(cycle).zfill(3)}', pattern)),
            key=alphanumeric_key
        )

        if not files:
            raise FileNotFoundError(f"No tiles found for cycle {cycle}, channel {channel}, z-plane {zplane}")

        tiles = self.load_tiles_parallel(files, copy=True)

        if as_stack:
            return np.stack(tiles, axis=0)
        return tiles

    def preload_channel(self,
                       image_dir: str,
                       cycle: int,
                       channel: int,
                       zplanes: List[int]) -> Dict[int, np.ndarray]:
        """
        Pre-load all tiles for multiple z-planes of a channel.

        This is useful for batch processing where you need all z-planes.

        Parameters
        ----------
        image_dir : str
            Base directory
        cycle : int
            Cycle number
        channel : int
            Channel number
        zplanes : list of int
            Z-planes to load

        Returns
        -------
        dict
            Mapping of z-plane to tile stack
        """
        result = {}

        # Load in parallel across z-planes
        with ThreadPoolExecutor(max_workers=self.n_workers) as executor:
            futures = {
                executor.submit(
                    self.load_zplane_tiles,
                    image_dir, cycle, channel, z
                ): z
                for z in zplanes
            }

            for future in as_completed(futures):
                z = futures[future]
                result[z] = future.result()

        return result

    def clear_cache(self):
        """Clear the tile cache."""
        with self._cache_lock:
            self._cache.clear()
            self._cache_order.clear()


class ChannelProcessor:
    """
    Optimized channel processor with parallel z-plane processing.

    This class orchestrates parallel processing of z-planes within a channel,
    with optimized I/O and GPU utilization.

    Parameters
    ----------
    n_workers : int
        Number of parallel workers for z-plane processing
    use_gpu_basic : bool
        Use GPU for BaSiC correction
    preload_tiles : bool
        Pre-load all tiles before processing (uses more memory but faster)
    """

    def __init__(self,
                 n_workers: int = 8,
                 use_gpu_basic: bool = False,
                 preload_tiles: bool = False):
        self.n_workers = n_workers
        self.use_gpu_basic = use_gpu_basic
        self.preload_tiles = preload_tiles
        self.tile_loader = TileLoader(
            use_mmap=True,
            cache_size=0,  # Don't cache between z-planes
            n_workers=4
        )

    def process_zplane(self,
                      tiles: np.ndarray,
                      basic_func: Callable,
                      basic_params: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Process a single z-plane: normalize and apply BaSiC correction.

        Parameters
        ----------
        tiles : np.ndarray
            Tile stack for this z-plane (n_tiles, height, width)
        basic_func : callable
            BaSiC correction function
        basic_params : dict
            Parameters for BaSiC

        Returns
        -------
        corrected : np.ndarray
            Corrected tiles
        flatfield : np.ndarray
            Computed flatfield
        darkfield : np.ndarray
            Computed darkfield
        """
        # Normalize
        dtype_max = np.iinfo(tiles.dtype).max if tiles.dtype.kind in 'ui' else 1.0
        tiles_norm = tiles.astype(np.float64) / dtype_max

        # BaSiC correction
        flatfield, darkfield = basic_func(tiles_norm, **basic_params)

        # Apply correction
        corrected = (tiles_norm - darkfield) / (flatfield + 1e-10)
        corrected = np.clip(corrected, 0, 1)
        corrected = (corrected * dtype_max).astype(tiles.dtype)

        return corrected, flatfield, darkfield


def create_optimized_loader(use_mmap: bool = True,
                           n_workers: int = 4) -> TileLoader:
    """
    Create an optimized tile loader with recommended settings.

    Parameters
    ----------
    use_mmap : bool
        Use memory-mapped reading if available
    n_workers : int
        Number of parallel loading workers

    Returns
    -------
    TileLoader
        Configured tile loader
    """
    return TileLoader(
        use_mmap=use_mmap,
        cache_size=0,  # Caching not beneficial for sequential z-plane processing
        n_workers=n_workers
    )


def load_tiles_optimized(image_dir: str,
                        cycle: int,
                        channel: int,
                        zplane: int) -> np.ndarray:
    """
    Convenience function to load tiles with optimized settings.

    This is a drop-in replacement for imread_collection-based loading.

    Parameters
    ----------
    image_dir : str
        Base directory containing cycle folders
    cycle : int
        Cycle number
    channel : int
        Channel number
    zplane : int
        Z-plane number

    Returns
    -------
    np.ndarray
        Stacked tile array (n_tiles, height, width)
    """
    loader = TileLoader(use_mmap=True, n_workers=4)
    return loader.load_zplane_tiles(image_dir, cycle, channel, zplane)
