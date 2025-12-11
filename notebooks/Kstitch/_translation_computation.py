import itertools
from typing import Tuple, List, Optional

import numpy as np
import numpy.typing as npt

from ._typing_utils import Float
from ._typing_utils import FloatArray
from ._typing_utils import Int
from ._typing_utils import IntArray
from ._typing_utils import NumArray

# GPU support
try:
    import cupy as cp
    HAS_CUPY = True
except ImportError:
    HAS_CUPY = False
    cp = None


def pcm(F1: np.ndarray, F2: np.ndarray) -> np.ndarray:
    """Phase correlation matrix computation (CPU version)."""
    FC = F1 * np.conjugate(F2)
    PCM = np.fft.ifft2(FC / np.abs(FC)).real.astype(np.float32)
    return PCM


def pcm_gpu(F1, F2):
    """Phase correlation matrix computation on GPU.

    Args:
        F1, F2: CuPy arrays (FFT results already on GPU)

    Returns:
        PCM as CuPy array on GPU
    """
    FC = F1 * cp.conjugate(F2)
    # Add small epsilon to avoid division by zero
    PCM = cp.fft.ifft2(FC / (cp.abs(FC) + 1e-10)).real.astype(cp.float32)
    return PCM


class GPUStitchingAccelerator:
    """GPU-accelerated stitching with pipelining and batch processing.

    This class provides optimized GPU operations for image stitching:
    - Batch FFT computation for multiple image pairs
    - PCM computation entirely on GPU (no CPU roundtrip)
    - Pre-allocated memory buffers to reduce allocation overhead
    - CUDA stream pipelining for overlapping compute and transfer
    """

    def __init__(self, image_shape: Tuple[int, int], max_batch_size: int = 16):
        """Initialize the GPU accelerator.

        Args:
            image_shape: (height, width) of images to process
            max_batch_size: Maximum number of image pairs to process in one batch
        """
        if not HAS_CUPY:
            raise RuntimeError("CuPy is required for GPU acceleration")

        self.image_shape = image_shape
        self.max_batch_size = max_batch_size
        self.height, self.width = image_shape

        # Create CUDA streams for pipelining
        self.compute_stream = cp.cuda.Stream(non_blocking=True)
        self.transfer_stream = cp.cuda.Stream(non_blocking=True)

        # Pre-allocate GPU buffers for batch processing
        self._allocate_buffers()

    def _allocate_buffers(self):
        """Pre-allocate GPU memory buffers."""
        h, w = self.image_shape
        batch = self.max_batch_size

        # Input image buffers (batch of pairs)
        self.img1_buffer = cp.empty((batch, h, w), dtype=cp.float32)
        self.img2_buffer = cp.empty((batch, h, w), dtype=cp.float32)

        # FFT result buffers
        self.fft1_buffer = cp.empty((batch, h, w), dtype=cp.complex64)
        self.fft2_buffer = cp.empty((batch, h, w), dtype=cp.complex64)

        # PCM result buffer
        self.pcm_buffer = cp.empty((batch, h, w), dtype=cp.float32)

    def compute_pcm_batch(self, image_pairs: List[Tuple[np.ndarray, np.ndarray]]) -> List[np.ndarray]:
        """Compute PCM for a batch of image pairs efficiently on GPU.

        Args:
            image_pairs: List of (image1, image2) tuples

        Returns:
            List of PCM arrays (on CPU)
        """
        n_pairs = len(image_pairs)
        if n_pairs == 0:
            return []

        # Process in chunks of max_batch_size
        results = []
        for start_idx in range(0, n_pairs, self.max_batch_size):
            end_idx = min(start_idx + self.max_batch_size, n_pairs)
            batch = image_pairs[start_idx:end_idx]
            batch_results = self._process_batch(batch)
            results.extend(batch_results)

        return results

    def _process_batch(self, batch: List[Tuple[np.ndarray, np.ndarray]]) -> List[np.ndarray]:
        """Process a single batch of image pairs."""
        n = len(batch)

        with self.transfer_stream:
            # Transfer images to GPU
            for i, (img1, img2) in enumerate(batch):
                self.img1_buffer[i] = cp.asarray(img1, dtype=cp.float32)
                self.img2_buffer[i] = cp.asarray(img2, dtype=cp.float32)

        # Wait for transfer to complete before computing
        self.transfer_stream.synchronize()

        with self.compute_stream:
            # Batch FFT - process all images at once
            self.fft1_buffer[:n] = cp.fft.fft2(self.img1_buffer[:n], axes=(1, 2))
            self.fft2_buffer[:n] = cp.fft.fft2(self.img2_buffer[:n], axes=(1, 2))

            # Batch PCM computation
            FC = self.fft1_buffer[:n] * cp.conjugate(self.fft2_buffer[:n])
            self.pcm_buffer[:n] = cp.fft.ifft2(FC / (cp.abs(FC) + 1e-10), axes=(1, 2)).real

        # Wait for compute to complete
        self.compute_stream.synchronize()

        # Transfer results back to CPU
        results = [cp.asnumpy(self.pcm_buffer[i]) for i in range(n)]

        return results

    def compute_fft_pair(self, image1: np.ndarray, image2: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """Compute FFT for a single image pair (legacy compatibility).

        Returns FFT results on CPU for compatibility with existing code.
        """
        with self.compute_stream:
            img1_gpu = cp.asarray(image1, dtype=cp.float32)
            img2_gpu = cp.asarray(image2, dtype=cp.float32)

            F1 = cp.fft.fft2(img1_gpu)
            F2 = cp.fft.fft2(img2_gpu)

            F1_cpu = cp.asnumpy(F1)
            F2_cpu = cp.asnumpy(F2)

        self.compute_stream.synchronize()
        return F1_cpu, F2_cpu

    def compute_pcm_single(self, image1: np.ndarray, image2: np.ndarray) -> np.ndarray:
        """Compute PCM for a single image pair entirely on GPU.

        This is more efficient than compute_fft_pair + pcm because
        the entire computation stays on GPU.
        """
        with self.compute_stream:
            img1_gpu = cp.asarray(image1, dtype=cp.float32)
            img2_gpu = cp.asarray(image2, dtype=cp.float32)

            F1 = cp.fft.fft2(img1_gpu)
            F2 = cp.fft.fft2(img2_gpu)

            # PCM entirely on GPU
            FC = F1 * cp.conjugate(F2)
            pcm_gpu = cp.fft.ifft2(FC / (cp.abs(FC) + 1e-10)).real.astype(cp.float32)

            pcm_cpu = cp.asnumpy(pcm_gpu)

        self.compute_stream.synchronize()
        return pcm_cpu

    def free_buffers(self):
        """Free pre-allocated GPU buffers."""
        del self.img1_buffer, self.img2_buffer
        del self.fft1_buffer, self.fft2_buffer
        del self.pcm_buffer
        cp.get_default_memory_pool().free_all_blocks()

    def __del__(self):
        """Cleanup on deletion."""
        try:
            self.free_buffers()
        except Exception:
            pass


# Global accelerator instance (lazy initialization)
_gpu_accelerator: Optional[GPUStitchingAccelerator] = None


def get_gpu_accelerator(image_shape: Tuple[int, int], max_batch_size: int = 16) -> GPUStitchingAccelerator:
    """Get or create the global GPU accelerator instance.

    Args:
        image_shape: (height, width) of images
        max_batch_size: Maximum batch size for processing

    Returns:
        GPUStitchingAccelerator instance
    """
    global _gpu_accelerator

    if _gpu_accelerator is None or _gpu_accelerator.image_shape != image_shape:
        if _gpu_accelerator is not None:
            _gpu_accelerator.free_buffers()
        _gpu_accelerator = GPUStitchingAccelerator(image_shape, max_batch_size)

    return _gpu_accelerator


def reset_gpu_accelerator():
    """Reset the global GPU accelerator (free memory)."""
    global _gpu_accelerator
    if _gpu_accelerator is not None:
        _gpu_accelerator.free_buffers()
        _gpu_accelerator = None


def multi_peak_max(PCM: FloatArray) -> Tuple[IntArray, IntArray, FloatArray]:
    """Find the first to n th largest peaks in PCM.

    Parameters
    ---------
    PCM : np.ndarray
        the peak correlation matrix

    Returns
    -------
    rows : np.ndarray
        the row indices for the peaks
    cols : np.ndarray
        the column indices for the peaks
    vals : np.ndarray
        the values of the peaks
    """
    row, col = np.unravel_index(np.argsort(PCM.ravel()), PCM.shape)
    vals: FloatArray = PCM[row[::-1], col[::-1]]
    return row[::-1], col[::-1], vals


def ncc(image1: NumArray, image2: NumArray) -> Float:
    """Compute the normalized cross correlation for two images.

    Parameters
    ---------
    image1 : np.ndarray
        the first image (the dimension must be 2)

    image2 : np.ndarray
        the second image (the dimension must be 2)

    Returns
    -------
    ncc : Float
        the normalized cross correlation
    """
    assert image1.ndim == 2
    assert image2.ndim == 2
    assert np.array_equal(image1.shape, image2.shape)
    image1 = image1.flatten()
    image2 = image2.flatten()
    n = np.dot(image1 - np.mean(image1), image2 - np.mean(image2))
    d = np.linalg.norm(image1) * np.linalg.norm(image2)
    return n / d



def extract_overlap_subregion(image: NumArray, y: Int, x: Int) -> NumArray:
    """Extract the overlapping subregion of the image.

    Parameters
    ---------
    image : np.ndarray
        the image (the dimension must be 2)
    y : Int
        the y (second last dim.) position
    x : Int
        the x (last dim.) position
    Returns
    -------
    subimage : np.ndarray
        the extracted subimage
    """
    sizeY = image.shape[0]
    sizeX = image.shape[1]
    assert (np.abs(y) < sizeY) and (np.abs(x) < sizeX)
    # clip x to (0, size_Y)
    xstart = int(max(0, min(y, sizeY, key=int), key=int))
    # clip x+sizeY to (0, size_Y)
    xend = int(max(0, min(y + sizeY, sizeY, key=int), key=int))
    ystart = int(max(0, min(x, sizeX, key=int), key=int))
    yend = int(max(0, min(x + sizeX, sizeX, key=int), key=int))
    return image[xstart:xend, ystart:yend]


def interpret_translation(
    image1: NumArray,
    image2: npt.NDArray,
    yins: IntArray,
    xins: IntArray,
    ymin: Int,
    ymax: Int,
    xmin: Int,
    xmax: Int,
    n: Int = 2,
) -> Tuple[float, int, int]:
    """Interpret the translation to find the translation with heighest ncc.

    Parameters
    ---------
    image1 : np.ndarray
        the first image (the dimension must be 2)
    image2 : np.ndarray
        the second image (the dimension must be 2)
    yins : IntArray
        the y positions estimated by PCM
    xins : IntArray
        the x positions estimated by PCM
    ymin : Int
        the minimum value of y (second last dim.)
    ymax : Int
        the maximum value of y (second last dim.)
    xmin : Int
        the minimum value of x (last dim.)
    xmax : Int
        the maximum value of x (last dim.)
    n : Int
        the number of peaks to check, default is 2.

    Returns
    -------
    _ncc : Float
        the highest ncc
    x : Int
        the selected x position
    y : Int
        the selected y position
    """
    assert image1.ndim == 2
    assert image2.ndim == 2
    assert np.array_equal(image1.shape, image2.shape)
    sizeY = image1.shape[0]
    sizeX = image1.shape[1]
    assert np.all(0 <= yins) and np.all(yins < sizeY)
    assert np.all(0 <= xins) and np.all(xins < sizeX)

    _ncc = -np.infty
    y = 0
    x = 0

    ymagss = [yins, sizeY - yins]
    ymagss[1][ymagss[0] == 0] = 0
    xmagss = [xins, sizeX - xins]
    xmagss[1][xmagss[0] == 0] = 0

    # concatenate all the candidates
    _poss = []
    for ymags, xmags, ysign, xsign in itertools.product(
        ymagss, xmagss, [-1, +1], [-1, +1]
    ):
        yvals = ymags * ysign
        xvals = xmags * xsign
        _poss.append([yvals, xvals])
    poss = np.array(_poss)
    valid_ind = (
        (ymin <= poss[:, 0, :])
        & (poss[:, 0, :] <= ymax)
        & (xmin <= poss[:, 1, :])
        & (poss[:, 1, :] <= xmax)
    )
    assert np.any(valid_ind)
    valid_ind = np.any(valid_ind, axis=0)
    for pos in np.moveaxis(poss[:, :, valid_ind], -1, 0)[: int(n)]:
        for yval, xval in pos:
            if (ymin <= yval) and (yval <= ymax) and (xmin <= xval) and (xval <= xmax):
                subI1 = extract_overlap_subregion(image1, yval, xval)
                subI2 = extract_overlap_subregion(image2, -yval, -xval)
                ncc_val = ncc(subI1, subI2)
                if ncc_val > _ncc:
                    _ncc = float(ncc_val)
                    y = int(yval)
                    x = int(xval)
    return _ncc, y, x
