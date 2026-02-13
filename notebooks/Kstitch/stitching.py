"""This module provides microscope image stitching with the algorithm by MIST."""
import itertools
import os
import sys
import warnings
from dataclasses import dataclass
import concurrent.futures
from typing import Any
from typing import Optional
from typing import Sequence
from typing import Tuple
from typing import Union
import gc
from multiprocessing import cpu_count

import numpy as np
import pandas as pd
from sklearn.covariance import EllipticEnvelope
from tqdm import tqdm


def _setup_cuda_path():
    """
    Ensure CUDA libraries and headers are accessible.

    On Windows: Conda installs CUDA libraries to Library/bin but doesn't add it to PATH.
    On Linux: CUDA headers are at targets/x86_64-linux/include but CuPy expects them
    at $CUDA_PATH/include. This function sets up the environment automatically.

    This function searches multiple locations for CUDA files:
    1. Conda environment paths
    2. NVIDIA CUDA Toolkit installation
    3. Common CUDA installation paths
    """
    conda_prefix = os.environ.get('CONDA_PREFIX', '')

    # Linux-specific setup
    if sys.platform.startswith('linux') and conda_prefix:
        # Set CUDA_PATH and CUDA_HOME for CuPy
        if not os.environ.get('CUDA_PATH'):
            os.environ['CUDA_PATH'] = conda_prefix
        if not os.environ.get('CUDA_HOME'):
            os.environ['CUDA_HOME'] = conda_prefix

        # Add library paths to LD_LIBRARY_PATH
        lib_paths = [
            os.path.join(conda_prefix, 'lib'),
            os.path.join(conda_prefix, 'targets', 'x86_64-linux', 'lib'),
        ]
        current_ld_path = os.environ.get('LD_LIBRARY_PATH', '')
        paths_to_add = [p for p in lib_paths if os.path.exists(p) and p not in current_ld_path]
        if paths_to_add:
            new_ld_path = os.pathsep.join(paths_to_add)
            if current_ld_path:
                new_ld_path += os.pathsep + current_ld_path
            os.environ['LD_LIBRARY_PATH'] = new_ld_path

        # Copy CUDA headers to include/ if needed (for CuPy JIT compilation)
        targets_include = os.path.join(conda_prefix, 'targets', 'x86_64-linux', 'include')
        conda_include = os.path.join(conda_prefix, 'include')
        cuda_header = os.path.join(conda_include, 'cuda_fp16.h')

        if os.path.exists(targets_include) and not os.path.exists(cuda_header):
            # Headers exist in targets but not in include - copy them
            import shutil
            try:
                for item in os.listdir(targets_include):
                    src = os.path.join(targets_include, item)
                    dst = os.path.join(conda_include, item)
                    if not os.path.exists(dst):
                        if os.path.isdir(src):
                            shutil.copytree(src, dst)
                        else:
                            shutil.copy2(src, dst)
            except (PermissionError, OSError):
                pass  # Can't copy, will rely on existing headers

        return

    # Windows-specific setup
    if sys.platform != 'win32':
        return

    import glob

    paths_to_add = []
    current_path = os.environ.get('PATH', '')

    # 1. Conda environment paths
    conda_prefix = os.environ.get('CONDA_PREFIX', '')
    if conda_prefix:
        conda_paths = [
            os.path.join(conda_prefix, 'Library', 'bin'),
            os.path.join(conda_prefix, 'bin'),
            # Some conda CUDA packages install to these locations
            os.path.join(conda_prefix, 'Library', 'lib'),
        ]
        for p in conda_paths:
            if os.path.exists(p) and p not in current_path:
                paths_to_add.append(p)

    # 2. CUDA_PATH environment variable (set by NVIDIA installer)
    cuda_path = os.environ.get('CUDA_PATH', '')
    if cuda_path:
        cuda_bin = os.path.join(cuda_path, 'bin')
        if os.path.exists(cuda_bin) and cuda_bin not in current_path:
            paths_to_add.append(cuda_bin)

    # 3. Common NVIDIA CUDA Toolkit locations
    program_files = os.environ.get('ProgramFiles', 'C:\\Program Files')
    cuda_base = os.path.join(program_files, 'NVIDIA GPU Computing Toolkit', 'CUDA')
    if os.path.exists(cuda_base):
        # Find all installed CUDA versions and prefer the highest
        cuda_versions = glob.glob(os.path.join(cuda_base, 'v*', 'bin'))
        cuda_versions.sort(reverse=True)  # Highest version first
        for cuda_bin in cuda_versions:
            if os.path.exists(cuda_bin) and cuda_bin not in current_path:
                paths_to_add.append(cuda_bin)
                break  # Only add the highest version

    # 4. Check for nvrtc specifically and find its location
    # This helps when CUDA is installed but not in standard locations
    if conda_prefix:
        nvrtc_files = glob.glob(os.path.join(conda_prefix, '**', 'nvrtc*.dll'), recursive=True)
        for nvrtc_file in nvrtc_files:
            nvrtc_dir = os.path.dirname(nvrtc_file)
            if nvrtc_dir not in current_path and nvrtc_dir not in paths_to_add:
                paths_to_add.append(nvrtc_dir)

    # Add all found paths to PATH
    if paths_to_add:
        new_path = os.pathsep.join(paths_to_add) + os.pathsep + current_path
        os.environ['PATH'] = new_path


# Fix CUDA PATH before any CUDA imports
_setup_cuda_path()

# GPU availability detection with actual functionality test
HAS_CUPY = False
GPU_ERROR_MSG = None
GPU_DEVICE_INFO = None

def _test_gpu_availability(verbose=False):
    """
    Test if CuPy and cuFFT are actually functional, not just importable.

    This performs a real cuFFT operation to verify the CUDA libraries are properly
    installed and accessible. Common failure modes:
    - cuFFT DLL not found (Windows): CUDA toolkit not installed or PATH issues
    - CUDA driver mismatch: CuPy version doesn't match installed CUDA version
    - No GPU available: System has no NVIDIA GPU or driver not installed

    Returns:
        bool: True if GPU is functional, False otherwise
    """
    global HAS_CUPY, GPU_ERROR_MSG, GPU_DEVICE_INFO
    try:
        import cupy as cp

        # Get device info first
        device_id = cp.cuda.runtime.getDevice()
        device_props = cp.cuda.runtime.getDeviceProperties(device_id)
        GPU_DEVICE_INFO = {
            'name': device_props['name'].decode() if isinstance(device_props['name'], bytes) else device_props['name'],
            'compute_capability': (device_props['major'], device_props['minor']),
            'total_memory_gb': device_props['totalGlobalMem'] / (1024**3),
        }

        # Test that cuFFT actually works by running a small FFT
        test_array = cp.zeros((4, 4), dtype=cp.float32)
        _ = cp.fft.fft2(test_array)
        cp.get_default_memory_pool().free_all_blocks()

        HAS_CUPY = True
        GPU_ERROR_MSG = None

        if verbose:
            print(f"GPU available: {GPU_DEVICE_INFO['name']}")
            print(f"  Compute capability: {GPU_DEVICE_INFO['compute_capability']}")
            print(f"  Memory: {GPU_DEVICE_INFO['total_memory_gb']:.1f} GB")

        return True

    except ImportError as e:
        HAS_CUPY = False
        GPU_ERROR_MSG = f"CuPy not installed: {e}"
        GPU_DEVICE_INFO = None
        return False

    except Exception as e:
        HAS_CUPY = False
        GPU_DEVICE_INFO = None
        error_str = str(e)

        # Provide specific guidance based on error type
        if "cufft" in error_str.lower() or "DLL load failed" in error_str:
            GPU_ERROR_MSG = (
                f"cuFFT library not found: {e}\n"
                "TROUBLESHOOTING:\n"
                "  1. Ensure CUDA Toolkit is installed (not just the driver)\n"
                "  2. Verify CUDA bin directory is in PATH:\n"
                "     Windows: C:\\Program Files\\NVIDIA GPU Computing Toolkit\\CUDA\\v1x.x\\bin\n"
                "  3. Ensure CuPy version matches CUDA version:\n"
                "     - CUDA 11.x: pip install cupy-cuda11x\n"
                "     - CUDA 12.x: pip install cupy-cuda12x\n"
                "  4. Restart Python/Jupyter after PATH changes"
            )
        elif "out of memory" in error_str.lower():
            GPU_ERROR_MSG = f"GPU out of memory: {e}"
        elif "driver" in error_str.lower():
            GPU_ERROR_MSG = (
                f"CUDA driver issue: {e}\n"
                "TROUBLESHOOTING:\n"
                "  1. Update NVIDIA GPU driver\n"
                "  2. Ensure driver version supports your CUDA toolkit version"
            )
        else:
            GPU_ERROR_MSG = f"GPU/cuFFT not functional: {e}"

        return False


def check_gpu_status(verbose=True):
    """Print detailed GPU status information for debugging.

    Parameters
    ----------
    verbose : bool
        If True, print additional diagnostic information about CUDA paths.
    """
    print("=" * 60)
    print("GPU STATUS CHECK")
    print("=" * 60)

    if HAS_CUPY:
        print(f"Status: AVAILABLE")
        if GPU_DEVICE_INFO:
            print(f"Device: {GPU_DEVICE_INFO['name']}")
            print(f"Compute Capability: {GPU_DEVICE_INFO['compute_capability']}")
            print(f"Total Memory: {GPU_DEVICE_INFO['total_memory_gb']:.1f} GB")
    else:
        print(f"Status: NOT AVAILABLE")
        print(f"\nError Details:\n{GPU_ERROR_MSG}")

    if verbose:
        print("\n" + "-" * 60)
        print("CUDA PATH DIAGNOSTICS")
        print("-" * 60)

        # Show environment variables
        conda_prefix = os.environ.get('CONDA_PREFIX', 'Not set')
        cuda_path = os.environ.get('CUDA_PATH', 'Not set')
        print(f"CONDA_PREFIX: {conda_prefix}")
        print(f"CUDA_PATH: {cuda_path}")

        # Search for nvrtc DLLs
        if sys.platform == 'win32':
            import glob
            print("\nSearching for NVRTC DLLs...")

            search_paths = []
            if conda_prefix != 'Not set':
                search_paths.append(conda_prefix)
            if cuda_path != 'Not set':
                search_paths.append(cuda_path)

            # Also check common CUDA locations
            program_files = os.environ.get('ProgramFiles', 'C:\\Program Files')
            cuda_base = os.path.join(program_files, 'NVIDIA GPU Computing Toolkit', 'CUDA')
            if os.path.exists(cuda_base):
                search_paths.append(cuda_base)

            nvrtc_found = []
            for base_path in search_paths:
                if os.path.exists(base_path):
                    files = glob.glob(os.path.join(base_path, '**', 'nvrtc*.dll'), recursive=True)
                    nvrtc_found.extend(files)

            if nvrtc_found:
                print(f"Found {len(nvrtc_found)} NVRTC DLL(s):")
                for f in nvrtc_found[:5]:  # Show first 5
                    print(f"  {f}")
                if len(nvrtc_found) > 5:
                    print(f"  ... and {len(nvrtc_found) - 5} more")
            else:
                print("No NVRTC DLLs found in searched paths!")
                print("\nTo fix this, try one of:")
                print("  1. Install CUDA Toolkit from nvidia.com")
                print("  2. Run: conda install -c nvidia cuda-nvrtc")
                print("  3. Reinstall cupy: conda install cupy")

    print("=" * 60)


# Test GPU on module load
_test_gpu_availability(verbose=False)

from ._constrained_refinement import refine_translations
from ._global_optimization import compute_final_position
from ._global_optimization import compute_maximum_spanning_tree
from ._stage_model import compute_image_overlap2
from ._stage_model import filter_by_overlap_and_correlation
from ._stage_model import filter_by_repeatability
from ._stage_model import filter_outliers
from ._stage_model import replace_invalid_translations
from ._translation_computation import interpret_translation
from ._translation_computation import multi_peak_max
from ._translation_computation import pcm
from ._translation_computation import get_gpu_accelerator
from ._translation_computation import get_multi_gpu_accelerator
from ._translation_computation import reset_gpu_accelerator
from ._typing_utils import BoolArray
from ._typing_utils import Float
from ._typing_utils import NumArray

# Get number of available threads to limit CPU thrashing
# From preadator: https://pypi.org/project/preadator/
if hasattr(os, "sched_getaffinity"):
    # On Linux, we can detect how many cores are assigned to this process.
    # This is especially useful when running in a Docker container, when the
    # number of cores is intentionally limited.
    NUM_THREADS = len(os.sched_getaffinity(0))  # type: ignore
else:
    # Default back to multiprocessing cpu_count, which is always going to count
    # the total number of cpus
    NUM_THREADS = cpu_count()

@dataclass
class ElipticEnvelopPredictor:
    contamination: float
    epsilon: float
    random_seed: int

    def __call__(self, X: NumArray) -> BoolArray:
        if len(X) < 2:
            return np.ones(len(X), dtype=bool)
        ee = EllipticEnvelope(contamination=self.contamination)
        rng = np.random.default_rng(self.random_seed)
        X = rng.normal(size=X.shape) * self.epsilon + X
        return ee.fit_predict(X) > 0

def _compute_fft_gpu(image1, image2):
    """Compute FFT on GPU using CuPy. Must be called from main process or thread pool."""
    import cupy as cp
    # Move data to GPU and compute FFT
    F1 = cp.fft.fft2(cp.asarray(image1))
    F2 = cp.fft.fft2(cp.asarray(image2))
    # Move back to CPU for PCM computation
    F1_cpu = cp.asnumpy(F1)
    F2_cpu = cp.asnumpy(F2)
    # Free GPU memory immediately
    del F1, F2
    cp.get_default_memory_pool().free_all_blocks()
    return F1_cpu, F2_cpu


def _compute_fft_cpu(image1, image2):
    """Compute FFT on CPU using NumPy."""
    F1 = np.fft.fft2(image1)
    F2 = np.fft.fft2(image2)
    return F1, F2


def process_image_pair_cpu(i2, g, direction, images, position_initial_guess, overlap_diff_threshold, sizeY, sizeX):
    """Process image pair on CPU. Safe for ProcessPoolExecutor."""
    i1 = g[direction]
    if pd.isna(i1):
        return None

    image1 = images[i1]
    image2 = images[i2]

    # Compute FFT on CPU
    F1, F2 = _compute_fft_cpu(image1, image2)
    PCM = pcm(F1, F2)

    if position_initial_guess is not None:
        def get_lims(dimension, size):
            val = g[f"{direction}_{dimension}_init_guess"]
            r = size * overlap_diff_threshold / 100.0
            return np.round([val - r, val + r]).astype(np.int64)

        lims = np.array(
            [
                get_lims(dimension, size)
                for dimension, size in zip("yx", [sizeY, sizeX])
            ]
        )
    else:
        lims = np.array([[-sizeY, sizeY], [-sizeX, sizeX]])

    yins, xins, _ = multi_peak_max(PCM)
    max_peak = interpret_translation(
        image1, image2, yins, xins, *lims[0], *lims[1]
    )
    return i2, direction, max_peak


def process_image_pair_gpu(i2, g, direction, images, position_initial_guess, overlap_diff_threshold, sizeY, sizeX, accelerator=None):
    """Process image pair on GPU. Must use ThreadPoolExecutor (not ProcessPoolExecutor).

    Args:
        accelerator: Optional GPUStitchingAccelerator for optimized processing.
                    If None, falls back to legacy GPU processing.
    """
    i1 = g[direction]
    if pd.isna(i1):
        return None

    image1 = images[i1]
    image2 = images[i2]

    # Use accelerator if available (keeps entire computation on GPU)
    if accelerator is not None:
        PCM = accelerator.compute_pcm_single(image1, image2)
    else:
        # Legacy path: FFT on GPU, PCM on CPU
        F1, F2 = _compute_fft_gpu(image1, image2)
        PCM = pcm(F1, F2)

    if position_initial_guess is not None:
        def get_lims(dimension, size):
            val = g[f"{direction}_{dimension}_init_guess"]
            r = size * overlap_diff_threshold / 100.0
            return np.round([val - r, val + r]).astype(np.int64)

        lims = np.array(
            [
                get_lims(dimension, size)
                for dimension, size in zip("yx", [sizeY, sizeX])
            ]
        )
    else:
        lims = np.array([[-sizeY, sizeY], [-sizeX, sizeX]])

    yins, xins, _ = multi_peak_max(PCM)
    max_peak = interpret_translation(
        image1, image2, yins, xins, *lims[0], *lims[1]
    )
    return i2, direction, max_peak


def process_pairs_batch_gpu(pairs_info, images, position_initial_guess, overlap_diff_threshold, sizeY, sizeX, accelerator):
    """Process multiple image pairs in a batch using GPU acceleration.

    This function batches FFT and PCM computation for better GPU utilization.

    Args:
        pairs_info: List of (i2, g, direction) tuples
        images: Array of all images
        position_initial_guess: Initial position guess or None
        overlap_diff_threshold: Overlap difference threshold
        sizeY, sizeX: Image dimensions
        accelerator: GPUStitchingAccelerator instance

    Returns:
        List of (i2, direction, max_peak) results
    """
    # Filter out invalid pairs and collect image data
    valid_pairs = []
    image_pairs = []

    for i2, g, direction in pairs_info:
        i1 = g[direction]
        if pd.isna(i1):
            continue
        valid_pairs.append((i2, g, direction, i1))
        image_pairs.append((images[i1], images[i2]))

    if not image_pairs:
        return []

    # Batch compute PCMs on GPU
    pcm_results = accelerator.compute_pcm_batch(image_pairs)

    # Process results
    results = []
    for (i2, g, direction, i1), PCM in zip(valid_pairs, pcm_results):
        image1 = images[i1]
        image2 = images[i2]

        if position_initial_guess is not None:
            def get_lims(dimension, size):
                val = g[f"{direction}_{dimension}_init_guess"]
                r = size * overlap_diff_threshold / 100.0
                return np.round([val - r, val + r]).astype(np.int64)

            lims = np.array(
                [
                    get_lims(dimension, size)
                    for dimension, size in zip("yx", [sizeY, sizeX])
                ]
            )
        else:
            lims = np.array([[-sizeY, sizeY], [-sizeX, sizeX]])

        yins, xins, _ = multi_peak_max(PCM)
        max_peak = interpret_translation(
            image1, image2, yins, xins, *lims[0], *lims[1]
        )
        results.append((i2, direction, max_peak))

    return results


def stitch_images(
    images: Union[Sequence[NumArray], NumArray],
    rows: Optional[Sequence[Any]] = None,
    cols: Optional[Sequence[Any]] = None,
    position_indices: Optional[NumArray] = None,
    position_initial_guess: Optional[NumArray] = None,
    overlap_diff_threshold: Float = 10,
    pou: Float = 0.5,
    full_output: bool = False,
    row_col_transpose: bool = False,
    initial_ncc_threshold: Float = 0.9,
    max_ncc_threshold: Float = -0.9,
    decrement_step: Float = -0.1,
    max_cores: int = NUM_THREADS//2,
    overlap_percentage: Optional[Float] = None,
    use_gpu: bool = False,
) -> Tuple[pd.DataFrame, dict]:
    """Compute image positions for stitching.

    Parameters
    ---------
    images : np.ndarray
        the images to stitch.

    rows : list, optional
        the row indices (tile position in the second last dimension) of the images.

    cols : list, optional
        the column indices (tile position in the last dimension) of the images

    position_indices : np.ndarray, optional
        the tile position indices in each dimension.
        the dimensions corresponds to (image, index)
        ignored if rows and cols are not None.

    position_initial_guess : np.ndarray, optional
        the initial guess for the positions of the images, in the unit of pixels.

    overlap_diff_threshold : 10
        the allowed difference from the initial guess, in percentage of the image size.
        ignored if position_initial_guess is None

    pou : Float, default 0.5
        the "percent overlap uncertainty" parameter

    full_output : bool, default False
        if True, returns the full comptutation result in the pd.DataFrame

    row_col_transpose : bool, default True
        if True, row and col indices are switched.
        only for compatibility and the default value will be False in the future.

    overlap_percentage : Float, optional
        the expected overlap percentage between adjacent tiles (0-100).
        If provided, this value will be used instead of automatic overlap detection.
        Should be the actual overlap percentage, e.g., 30 for 30% overlap.

    use_gpu : bool, default False
        if True, use GPU acceleration for FFT computations via CuPy.
        Requires CuPy and CUDA toolkit to be properly installed.
        Raises RuntimeError if GPU is not available (no CPU fallback).

    ncc_threshold : Float, default 0.5
        the threshold of the normalized cross correlation used to select the initial
        stitched pairs.

    Returns
    -------
    grid : pd.DataFrame
        the result dataframe with the rows "x_pos" and "y_pos" whose values are
        the absolute positions.

    prop_dict : dict
        the dict of estimated parameters. (to be documented)
    """
    ncc_threshold = initial_ncc_threshold
    
    images = np.array(images)
    assert (position_indices is not None) or (rows is not None and cols is not None)
    if position_indices is None:
        if row_col_transpose:
            warnings.warn(
                "row_col_transpose is True. The default value will be changed to False in the major release."
            )
            position_indices = np.array([cols, rows]).T
        else:
            position_indices = np.array([rows, cols]).T
    position_indices = np.array(position_indices)
    assert images.shape[0] == position_indices.shape[0]
    assert position_indices.shape[1] == images.ndim - 1
    if position_initial_guess is not None:
        position_initial_guess = np.array(position_initial_guess)
        assert images.shape[0] == position_indices.shape[0]
        assert position_initial_guess.shape[1] == images.ndim - 1
    assert 0 <= overlap_diff_threshold and overlap_diff_threshold <= 100
    _rows, _cols = position_indices.T

    sizeY, sizeX = images.shape[1:]

    grid = pd.DataFrame(
        {
            "col": _cols,
            "row": _rows,
        },
        index=np.arange(len(_cols)),
    )

    def get_index(col, row):
        df = grid[(grid["col"] == col) & (grid["row"] == row)]
        assert len(df) < 2
        if len(df) == 1:
            return df.index[0]
        else:
            return None

    grid["top"] = grid.apply(
        lambda g: get_index(g["col"], g["row"] - 1), axis=1
    ).astype(pd.Int32Dtype())
    grid["left"] = grid.apply(
        lambda g: get_index(g["col"] - 1, g["row"]), axis=1
    ).astype(pd.Int32Dtype())

    ### dimension order ... m.y.x
    if position_initial_guess is not None:
        for j, dimension in enumerate(["y", "x"]):
            grid[f"{dimension}_pos_init_guess"] = position_initial_guess[:, j]
        for direction, dimension in itertools.product(["left", "top"], ["y", "x"]):
            for ind, g in grid.iterrows():
                i1 = g[direction]
                if pd.isna(i1):
                    continue
                g2 = grid.loc[i1]
                grid.loc[ind, f"{direction}_{dimension}_init_guess"] = (
                    g[f"{dimension}_pos_init_guess"] - g2[f"{dimension}_pos_init_guess"]
                )

    # Determine if GPU can actually be used
    actual_use_gpu = use_gpu
    if use_gpu:
        if not HAS_CUPY:
            error_msg = (
                f"GPU requested but not available.\n"
                f"Error: {GPU_ERROR_MSG}\n\n"
                f"Run 'from Kstitch.stitching import check_gpu_status; check_gpu_status()' "
                f"for detailed diagnostics."
            )
            raise RuntimeError(error_msg)
        else:
            try:
                import cupy as cp
                # Verify GPU is still functional (not just at import time)
                device_id = cp.cuda.runtime.getDevice()
                device_props = cp.cuda.runtime.getDeviceProperties(device_id)
                gpu_name = device_props['name']
                if isinstance(gpu_name, bytes):
                    gpu_name = gpu_name.decode()
                print(f"Using GPU: {gpu_name}")
                cp.cuda.Device(device_id).use()
            except Exception as e:
                error_msg = (
                    f"GPU initialization failed: {e}\n\n"
                    f"Run 'from Kstitch.stitching import check_gpu_status; check_gpu_status()' "
                    f"for detailed diagnostics."
                )
                raise RuntimeError(error_msg)

    # Process all image pairs once to compute correlations
    print("Computing phase correlations for all image pairs...")

    # CRITICAL: GPU operations CANNOT use ProcessPoolExecutor because CUDA contexts
    # don't transfer across process boundaries. Use ThreadPoolExecutor for GPU
    # (threads share the same CUDA context) or ProcessPoolExecutor for CPU
    # (processes give better parallelism for CPU-bound work).

    if actual_use_gpu:
        # GPU mode with optimized batch processing
        # Use multi-GPU accelerator if available, otherwise single GPU
        try:
            accelerator = get_multi_gpu_accelerator((sizeY, sizeX), max_batch_size=16)
            n_gpus = len(accelerator.device_ids)
            print(f"GPU mode: using {n_gpus} GPU(s) with batch processing")
        except Exception:
            accelerator = get_gpu_accelerator((sizeY, sizeX), max_batch_size=16)
            print(f"GPU mode: using single GPU with batch processing")

        for direction in ["left", "top"]:
            # Collect all pairs for this direction
            pairs_info = [(i2, g, direction) for i2, g in grid.iterrows()]

            # Process in batches for better GPU utilization
            batch_size = 32  # Process 32 pairs at a time
            all_results = []

            for batch_start in tqdm(range(0, len(pairs_info), batch_size),
                                    desc=f"Processing {direction} pairs (batch)"):
                batch_end = min(batch_start + batch_size, len(pairs_info))
                batch = pairs_info[batch_start:batch_end]

                # Batch process on GPU
                batch_results = process_pairs_batch_gpu(
                    batch, images, position_initial_guess,
                    overlap_diff_threshold, sizeY, sizeX, accelerator
                )
                all_results.extend(batch_results)

            # Store results
            for result in all_results:
                if result is not None:
                    i2, dir_result, max_peak = result
                    for j, key in enumerate(["ncc", "y", "x"]):
                        grid.loc[i2, f"{dir_result}_{key}_first"] = max_peak[j]

    else:
        # CPU mode: Use ProcessPoolExecutor (better parallelism)
        ExecutorClass = concurrent.futures.ProcessPoolExecutor
        process_fn = process_image_pair_cpu
        print(f"CPU mode: using {max_cores} processes")

        with ExecutorClass(max_workers=max_cores) as executor:
            for direction in ["left", "top"]:
                futures = [
                    executor.submit(
                        process_fn, i2, g, direction, images,
                        position_initial_guess, overlap_diff_threshold, sizeY, sizeX
                    )
                    for i2, g in grid.iterrows()
                ]

                for future in tqdm(concurrent.futures.as_completed(futures), desc=f"Processing {direction} pairs", total=len(futures)):
                    result = future.result()
                    if result is not None:
                        i2, direction, max_peak = result
                        for j, key in enumerate(["ncc", "y", "x"]):
                            grid.loc[i2, f"{direction}_{key}_first"] = max_peak[j]

    # Analyze actual NCC values to select appropriate threshold
    top_ncc_values = grid["top_ncc_first"].dropna()
    left_ncc_values = grid["left_ncc_first"].dropna()
    all_ncc_values = pd.concat([top_ncc_values, left_ncc_values])
    
    if len(all_ncc_values) > 0:
        min_ncc = all_ncc_values.min()
        mean_ncc = all_ncc_values.mean()
        max_ncc = all_ncc_values.max()
        
        print(f"NCC Statistics: min={min_ncc:.3f}, mean={mean_ncc:.3f}, max={max_ncc:.3f}")
        
        # Select threshold based on actual data
        if initial_ncc_threshold > max_ncc:
            # Requested threshold is higher than any actual NCC - use a more conservative value
            ncc_threshold = max(min_ncc, mean_ncc * 0.5)  # Use half of mean or minimum, whichever is higher
            print(f"Initial threshold {initial_ncc_threshold:.1f} > max NCC {max_ncc:.3f}")
            print(f"Using conservative adaptive threshold: {ncc_threshold:.3f}")
        elif initial_ncc_threshold >= min_ncc:
            # Requested threshold is achievable
            ncc_threshold = initial_ncc_threshold
            print(f"Using requested threshold: {ncc_threshold:.3f}")
        else:
            # Even the initial threshold is below minimum - use minimum
            ncc_threshold = min_ncc
            print(f"Using minimum possible threshold: {ncc_threshold:.3f}")
    else:
        ncc_threshold = 0.0
        print("No NCC values found, using threshold 0.0")
    predictor = ElipticEnvelopPredictor(contamination=0.1, epsilon=0.01, random_seed=0)
    
    # Use provided overlap_percentage if available, otherwise compute automatically
    if overlap_percentage is not None:
        # User provided overlap percentage - don't clip by POU here, let filter function handle it
        overlap_top = overlap_percentage
        overlap_left = overlap_percentage
        print(f"Using provided overlap_percentage: {overlap_percentage}% with POU: {pou}%")
    else:
        # Compute overlap automatically from phase correlation analysis
        print("Computing overlap from phase correlation analysis...")
        left_displacement = compute_image_overlap2(
            grid[grid["left_ncc_first"] >= ncc_threshold], "left", sizeY, sizeX, predictor
        )
        top_displacement = compute_image_overlap2(
            grid[grid["top_ncc_first"] >= ncc_threshold], "top", sizeY, sizeX, predictor
        )
        overlap_top = np.clip(100 - top_displacement[0] * 100, pou, 100 - pou)
        overlap_left = np.clip(100 - left_displacement[1] * 100, pou, 100 - pou)
        print(f"Computed overlap automatically - top: {overlap_top:.1f}%, left: {overlap_left:.1f}%")

    ### compute_repeatability ###
    # Progressively relax POU if the initial filter rejects ALL translations
    # for a direction. This prevents tight POU values from collapsing an axis.
    pou_schedule = [pou]
    for candidate in [0.5, 1, 1.5, 2, 3, 5, 10, 20]:
        if candidate > pou and candidate not in pou_schedule:
            pou_schedule.append(candidate)

    for current_pou in pou_schedule:
        grid["top_valid1"] = filter_by_overlap_and_correlation(
            grid["top_y_first"],
            grid["top_ncc_first"],
            overlap_top,
            sizeY,
            current_pou,
            ncc_threshold,
            overlap_percentage
        )
        grid["top_valid2"] = filter_outliers(grid["top_y_first"], grid["top_valid1"])
        grid["left_valid1"] = filter_by_overlap_and_correlation(
            grid["left_x_first"],
            grid["left_ncc_first"],
            overlap_left,
            sizeX,
            current_pou,
            ncc_threshold,
            overlap_percentage
        )
        grid["left_valid2"] = filter_outliers(grid["left_x_first"], grid["left_valid1"])

        n_top = grid["top_valid2"].sum()
        n_left = grid["left_valid2"].sum()

        if n_top > 0 and n_left > 0:
            if current_pou != pou:
                print(f"POU relaxed from {pou} to {current_pou} "
                      f"(top_valid={n_top}, left_valid={n_left})")
            break
        else:
            failed = []
            if n_top == 0:
                failed.append("top")
            if n_left == 0:
                failed.append("left")
            if current_pou == pou_schedule[-1]:
                print(f"WARNING: {'/'.join(failed)} translations all invalid "
                      f"even at POU={current_pou}. "
                      f"replace_invalid_translations will use raw medians.")
            else:
                print(f"POU={current_pou}: {'/'.join(failed)} has 0 valid "
                      f"translations, relaxing...")

    rs = []
    for direction, dims, rowcol in zip(["top", "left"], ["yx", "xy"], ["col", "row"]):
        valid_key = f"{direction}_valid2"
        valid_grid = grid[grid[valid_key]]
        if len(valid_grid) > 0:
            w1s = valid_grid[f"{direction}_{dims[0]}_first"]
            r1 = np.ceil((w1s.max() - w1s.min()) / 2)
            _, w2s = zip(*valid_grid.groupby(rowcol)[f"{direction}_{dims[1]}_first"])
            r2 = np.ceil(np.max([np.max(w2) - np.min(w2) for w2 in w2s]) / 2)
            rs.append(max(r1, r2))
        rs.append(0)
    r = np.max(rs)
    
    
    grid = filter_by_repeatability(grid, r, ncc_threshold)
    grid = replace_invalid_translations(grid)

    grid = refine_translations(images, grid, r)

    tree = compute_maximum_spanning_tree(grid)
    grid = compute_final_position(grid, tree)

    # Clean up GPU memory if it was used
    if actual_use_gpu:
        try:
            # Reset the GPU accelerator to free pre-allocated buffers
            reset_gpu_accelerator()
            import cupy as cp
            cp.get_default_memory_pool().free_all_blocks()
            cp.get_default_pinned_memory_pool().free_all_blocks()
        except Exception:
            pass  # Ignore cleanup errors

    prop_dict = {
        "W": sizeY,
        "H": sizeX,
        "overlap_left": overlap_left,
        "overlap_top": overlap_top,
        "repeatability": r,
    }
    if row_col_transpose:
        grid = grid.rename(columns={"x_pos": "y_pos", "y_pos": "x_pos"})
    if full_output:
        return grid, prop_dict
    else:
        return grid[["row", "col", "y_pos", "x_pos"]], prop_dict
