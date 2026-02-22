"""
3D Vessel Segmentation Module

Provides vessel segmentation from deconvolved z-stacks using Frangi vesselness
filtering, morphological cleanup, skeletonization, and graph-based morphometry.

Pipeline position: Runs after deconvolution (step 3) on the full 3D z-stack,
before EDF collapses the z-dimension. The deconvolved stack retains full
volumetric information needed for skeletonization and radius estimation.

Supports:
- GPU acceleration via CuPy (recommended for Frangi filter)
- CPU processing via SciPy/scikit-image (fallback)
- Tiled processing for memory-constrained systems
- Anisotropy-aware processing (typical Z/XY ratio ~4:1)

Usage:
    from kintsugi.vessel3d import (
        compute_vesselness_frangi,
        binarize_vessel_mask,
        skeletonize_vessels,
        analyze_vessel_graph,
        export_vessel_results,
    )

    # Full pipeline
    vesselness = compute_vesselness_frangi(stack, sigmas=[1, 2, 4, 8], spacing=spacing)
    mask = binarize_vessel_mask(vesselness, min_size=500, closing_radius=1)
    skeleton = skeletonize_vessels(mask)
    features = analyze_vessel_graph(skeleton, mask, spacing=spacing)
    export_vessel_results(results, output_dir, marker_name="CD31")
"""

from __future__ import annotations

import gc
import logging
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np
from scipy import ndimage
from scipy.ndimage import binary_fill_holes, distance_transform_edt

if TYPE_CHECKING:
    import cupy as cp
    import pandas as pd

logger = logging.getLogger(__name__)

# Type aliases
ArrayLike = Union[np.ndarray, "cp.ndarray"]


# =============================================================================
# Data Classes
# =============================================================================


@dataclass
class VesselSpacing:
    """Physical voxel spacing for anisotropic volumes.

    Parameters
    ----------
    xy : float
        XY pixel size in micrometers (e.g., 0.377 for 377 nm).
    z : float
        Z step size in micrometers (e.g., 1.5 for 1500 nm).
    """

    xy: float = 0.377
    z: float = 1.5

    @property
    def ratio(self) -> float:
        """Z-to-XY anisotropy ratio."""
        return self.z / self.xy

    @property
    def as_zyx(self) -> tuple[float, float, float]:
        """Spacing as (z, y, x) tuple in micrometers."""
        return (self.z, self.xy, self.xy)

    @property
    def isotropic_zyx(self) -> tuple[float, float, float]:
        """Isotropic spacing after z-resampling (all equal to xy)."""
        return (self.xy, self.xy, self.xy)

    @classmethod
    def from_experiment(cls, experiment_config: dict) -> VesselSpacing:
        """Create from experiment.json config dict."""
        xy = experiment_config.get("xy_pixel_size", 0.377)
        z = experiment_config.get("z_step_size", 1.5)
        return cls(xy=xy, z=z)


@dataclass
class VesselSegmentationResult:
    """Container for all vessel segmentation outputs.

    Attributes
    ----------
    binary_mask : np.ndarray
        3D binary vessel mask (isotropic resolution).
    binary_mask_native : np.ndarray | None
        3D binary mask at native (anisotropic) resolution.
    skeleton : np.ndarray | None
        3D skeleton image (1-voxel-wide centerlines).
    labeled_segments : np.ndarray | None
        3D labeled vessel segments.
    features : pd.DataFrame | None
        Per-segment morphometric features.
    vesselness : np.ndarray | None
        Frangi vesselness map (probability-like).
    spacing : VesselSpacing
        Physical voxel spacing used.
    marker_name : str
        Name of the vessel marker (e.g., "CD31").
    """

    binary_mask: np.ndarray
    binary_mask_native: np.ndarray | None = None
    skeleton: np.ndarray | None = None
    labeled_segments: np.ndarray | None = None
    features: pd.DataFrame | None = None
    vesselness: np.ndarray | None = None
    spacing: VesselSpacing = field(default_factory=VesselSpacing)
    marker_name: str = "vessel"


# =============================================================================
# Device Detection (follows edf.py pattern)
# =============================================================================


def _check_cupy_available() -> bool:
    """Check if CuPy can execute GPU operations."""
    try:
        import cupy as cp

        _ = cp.array([1, 2, 3])
        return True
    except (ImportError, Exception):
        return False


def _detect_device(device: str = "auto") -> tuple[bool, int, str]:
    """Detect best device for processing.

    Parameters
    ----------
    device : str
        'gpu', 'cpu', or 'auto'.

    Returns
    -------
    tuple
        (use_gpu, device_id, device_name)
    """
    if device.lower() == "cpu":
        return False, -1, "CPU"

    try:
        from kintsugi.gpu import get_gpu_manager

        gpu = get_gpu_manager()
        if gpu.has_gpu and gpu.cupy_available:
            import cupy as cp

            device_id = gpu.primary_device
            cp.cuda.Device(device_id).use()
            device_name = gpu.gpu_info[0].name if gpu.gpu_info else "GPU"
            return True, device_id, f"GPU ({device_name})"
    except ImportError:
        pass
    except Exception as e:
        logger.warning(f"GPUManager not available: {e}")

    if device.lower() in ("gpu", "auto") and _check_cupy_available():
        return True, 0, "GPU 0"

    return False, -1, "CPU"


# =============================================================================
# Preprocessing
# =============================================================================


def resample_isotropic(
    volume: np.ndarray,
    spacing: VesselSpacing,
    order: int = 3,
    device: str = "auto",
) -> np.ndarray:
    """Resample a volume to isotropic voxel spacing by upsampling the z-axis.

    Parameters
    ----------
    volume : np.ndarray
        3D array with shape (z, y, x).
    spacing : VesselSpacing
        Physical voxel spacing.
    order : int
        Interpolation order (3 = cubic).
    device : str
        'gpu', 'cpu', or 'auto'.

    Returns
    -------
    np.ndarray
        Resampled volume with isotropic voxels (z upsampled by spacing.ratio).
    """
    zoom_factor = (spacing.ratio, 1.0, 1.0)
    logger.info(
        f"Resampling z-axis by {spacing.ratio:.2f}x "
        f"(z={spacing.z} um -> {spacing.xy} um, shape {volume.shape})"
    )

    use_gpu, device_id, device_name = _detect_device(device)

    if use_gpu:
        try:
            import cupy as cp
            from cupyx.scipy.ndimage import zoom as gpu_zoom

            cp.cuda.Device(device_id).use()
            vol_gpu = cp.asarray(volume.astype(np.float32))
            result = gpu_zoom(vol_gpu, zoom_factor, order=order)
            result = cp.asnumpy(result)
            logger.info(f"Resampled on {device_name}: {volume.shape} -> {result.shape}")

            del vol_gpu
            cp.get_default_memory_pool().free_all_blocks()

            if np.issubdtype(volume.dtype, np.integer):
                result = np.clip(result, np.iinfo(volume.dtype).min, np.iinfo(volume.dtype).max)
            return result.astype(volume.dtype)
        except Exception as e:
            logger.warning(f"GPU resampling failed ({e}), falling back to CPU")

    result = ndimage.zoom(volume.astype(np.float32), zoom_factor, order=order)
    logger.info(f"Resampled on CPU: {volume.shape} -> {result.shape}")

    if np.issubdtype(volume.dtype, np.integer):
        result = np.clip(result, np.iinfo(volume.dtype).min, np.iinfo(volume.dtype).max)
    return result.astype(volume.dtype)


def preprocess_volume(
    volume: np.ndarray,
    spacing: VesselSpacing,
    denoise_sigma: float = 0.5,
    make_isotropic: bool = True,
    device: str = "auto",
) -> np.ndarray:
    """Preprocess a deconvolved z-stack for vessel segmentation.

    Steps:
    1. Optionally resample to isotropic voxels (z upsampling).
    2. Light Gaussian denoising to suppress shot noise.

    Parameters
    ----------
    volume : np.ndarray
        3D array (z, y, x), typically uint16 from deconvolution.
    spacing : VesselSpacing
        Physical voxel spacing.
    denoise_sigma : float
        Gaussian sigma for denoising (in isotropic voxels). 0 to skip.
    make_isotropic : bool
        If True, resample z-axis to match xy spacing.
    device : str
        'gpu', 'cpu', or 'auto'.

    Returns
    -------
    np.ndarray
        Preprocessed volume (float32 if denoised, original dtype otherwise).
    """
    if volume.ndim != 3:
        raise ValueError(f"Expected 3D volume, got shape {volume.shape}")

    result = volume

    if make_isotropic:
        result = resample_isotropic(result, spacing, device=device)

    if denoise_sigma > 0:
        use_gpu, device_id, device_name = _detect_device(device)
        logger.info(f"Denoising with sigma={denoise_sigma}")

        if use_gpu:
            try:
                import cupy as cp
                from cupyx.scipy.ndimage import gaussian_filter as gpu_gaussian

                cp.cuda.Device(device_id).use()
                vol_gpu = cp.asarray(result.astype(np.float32))
                vol_gpu = gpu_gaussian(vol_gpu, sigma=denoise_sigma)
                result = cp.asnumpy(vol_gpu)
                del vol_gpu
                cp.get_default_memory_pool().free_all_blocks()
                return result
            except Exception as e:
                logger.warning(f"GPU denoising failed ({e}), falling back to CPU")

        result = ndimage.gaussian_filter(result.astype(np.float32), sigma=denoise_sigma)

    return result


# =============================================================================
# Frangi Vesselness Filter
# =============================================================================


def _hessian_eigenvalues_3d(
    volume: np.ndarray,
    sigma: float,
    use_gpu: bool = False,
    device_id: int = 0,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Compute sorted Hessian eigenvalues at a given scale.

    Parameters
    ----------
    volume : np.ndarray
        3D float32 volume (should be isotropic).
    sigma : float
        Gaussian scale for derivative computation.
    use_gpu : bool
        Use CuPy GPU acceleration.
    device_id : int
        GPU device ID.

    Returns
    -------
    tuple of np.ndarray
        (lambda1, lambda2, lambda3) sorted by absolute value,
        |lambda1| <= |lambda2| <= |lambda3|.
    """
    if use_gpu:
        import cupy as cp

        # VRAM guard: Cardano eigenvalue path needs ~15 float32 arrays at peak
        # (6 Hessian + input + ~8 working arrays during Cardano + sorting).
        # Fall back to CPU if GPU has insufficient VRAM.
        volume_bytes = volume.nbytes
        estimated_vram = volume_bytes * 15
        try:
            cp.cuda.Device(device_id).use()
            free_vram = cp.cuda.Device(device_id).mem_info[0]
            if estimated_vram > free_vram:
                logger.warning(
                    f"GPU VRAM insufficient for Hessian eigenvalues: "
                    f"need ~{estimated_vram / 1e9:.1f} GB, "
                    f"available {free_vram / 1e9:.1f} GB. "
                    f"Falling back to CPU."
                )
                use_gpu = False
        except Exception as e:
            logger.warning(f"Could not query GPU VRAM ({e}), falling back to CPU")
            use_gpu = False

    if use_gpu:
        import cupy as cp
        from cupyx.scipy.ndimage import gaussian_filter as gf

        cp.cuda.Device(device_id).use()
        vol = cp.asarray(volume)
        xp = cp
    else:
        from scipy.ndimage import gaussian_filter as gf

        vol = volume
        xp = np

    # Compute second-order Gaussian derivatives (Hessian elements)
    # H = [[Hzz, Hzy, Hzx],
    #       [Hyz, Hyy, Hyx],
    #       [Hxz, Hxy, Hxx]]
    # Using order parameter of gaussian_filter for derivatives

    # Scale-normalized: multiply by sigma^2
    s2 = sigma**2

    Hzz = gf(vol, sigma=sigma, order=[2, 0, 0]) * s2
    Hyy = gf(vol, sigma=sigma, order=[0, 2, 0]) * s2
    Hxx = gf(vol, sigma=sigma, order=[0, 0, 2]) * s2
    Hzy = gf(vol, sigma=sigma, order=[1, 1, 0]) * s2
    Hzx = gf(vol, sigma=sigma, order=[1, 0, 1]) * s2
    Hyx = gf(vol, sigma=sigma, order=[0, 1, 1]) * s2

    # Analytical eigenvalues for symmetric 3x3 Hessian via Cardano's method.
    # Operates element-wise on the 6 Hessian arrays — avoids allocating the
    # massive (N, 3, 3) matrix that caused GPU OOM on large stitched volumes.
    #
    # For symmetric 3x3 matrix [[a,b,c],[b,d,e],[c,e,f]]:
    #   a=Hzz, b=Hzy, c=Hzx, d=Hyy, e=Hyx, f=Hxx

    del vol  # free GPU copy of input volume

    # Trace / 3
    q = (Hzz + Hyy + Hxx) / 3.0

    # Shift diagonal: a-q, d-q, f-q (free diagonals immediately after)
    a_q = Hzz - q
    del Hzz
    d_q = Hyy - q
    del Hyy
    f_q = Hxx - q
    del Hxx

    # p² = ( (a-q)² + (d-q)² + (f-q)² + 2*(b² + c² + e²) ) / 6
    p2 = (a_q * a_q + d_q * d_q + f_q * f_q + 2.0 * (Hzy * Hzy + Hzx * Hzx + Hyx * Hyx)) / 6.0

    p = xp.sqrt(xp.maximum(p2, xp.float32(1e-30)))
    del p2

    # Determinant of B = (A - q*I) / p
    # det(B) = (1/p³) * det(A - q*I)
    # det([[a-q, b, c],[b, d-q, e],[c, e, f-q]])
    inv_p = 1.0 / p
    b_ = Hzy * inv_p
    c_ = Hzx * inv_p
    e_ = Hyx * inv_p
    aq = a_q * inv_p
    dq = d_q * inv_p
    fq = f_q * inv_p
    del a_q, d_q, f_q, inv_p

    det_B = aq * (dq * fq - e_ * e_) - b_ * (b_ * fq - e_ * c_) + c_ * (b_ * e_ - dq * c_)
    del aq, dq, fq, b_, c_, e_

    # r = det_B / 2, clamp to [-1, 1] for numerical stability
    r = xp.clip(det_B * xp.float32(0.5), xp.float32(-1.0), xp.float32(1.0))
    del det_B

    # phi = arccos(r) / 3
    phi = xp.arccos(r) / 3.0
    del r

    # Eigenvalues (sorted: eig1 >= eig2 >= eig3)
    TWO_PI_3 = xp.float32(2.0 * np.pi / 3.0)
    two_p = 2.0 * p
    del p

    eig1 = q + two_p * xp.cos(phi)
    eig3 = q + two_p * xp.cos(phi + TWO_PI_3)
    eig2 = 3.0 * q - eig1 - eig3  # trace identity
    del phi, two_p, q

    # Free off-diagonal Hessian components (diagonals already freed above)
    del Hzy, Hzx, Hyx

    # Sort by absolute value using a 3-element sorting network (Knuth).
    # Three compare-and-swap steps guarantee |l1| <= |l2| <= |l3|.
    # Each step creates 2 new arrays and frees the 2 old ones.

    # Step 1: compare-swap positions 0 and 2
    swap = xp.abs(eig1) > xp.abs(eig3)
    eig1, eig3 = xp.where(swap, eig3, eig1), xp.where(swap, eig1, eig3)
    del swap

    # Step 2: compare-swap positions 0 and 1
    swap = xp.abs(eig1) > xp.abs(eig2)
    eig1, eig2 = xp.where(swap, eig2, eig1), xp.where(swap, eig1, eig2)
    del swap

    # Step 3: compare-swap positions 1 and 2
    swap = xp.abs(eig2) > xp.abs(eig3)
    eig2, eig3 = xp.where(swap, eig3, eig2), xp.where(swap, eig2, eig3)
    del swap

    l1, l2, l3 = eig1, eig2, eig3

    if use_gpu:
        l1 = cp.asnumpy(l1)
        l2 = cp.asnumpy(l2)
        l3 = cp.asnumpy(l3)
        cp.get_default_memory_pool().free_all_blocks()

    return l1, l2, l3


def compute_vesselness_frangi(
    volume: np.ndarray,
    sigmas: list[float] | None = None,
    alpha: float = 0.5,
    beta: float = 0.5,
    gamma: float | None = None,
    spacing: VesselSpacing | None = None,
    black_ridges: bool = False,
    device: str = "auto",
) -> np.ndarray:
    """Compute 3D Frangi vesselness filter at multiple scales.

    The Frangi filter enhances tubular structures by analyzing the eigenvalues
    of the Hessian matrix at multiple Gaussian scales. It returns a vesselness
    response map where high values indicate tube-like structures.

    Parameters
    ----------
    volume : np.ndarray
        3D array (z, y, x). Should be isotropic or anisotropy-corrected.
    sigmas : list of float, optional
        Gaussian scales for multi-scale analysis. Default [1, 2, 4, 8].
        Should span the range of expected vessel radii in voxels.
    alpha : float
        Sensitivity to plate-like structures (Ra ratio). Default 0.5.
    beta : float
        Sensitivity to blob-like structures (Rb ratio). Default 0.5.
    gamma : float, optional
        Sensitivity to background noise (Frobenius norm). If None,
        automatically set to half the maximum Hessian norm.
    spacing : VesselSpacing, optional
        Physical spacing (for logging only; volume should already be isotropic).
    black_ridges : bool
        If True, detect dark vessels on bright background. Default False
        (bright vessels on dark background, typical for fluorescence).
    device : str
        'gpu', 'cpu', or 'auto'.

    Returns
    -------
    np.ndarray
        Vesselness response map (float32, range [0, 1]).
        Maximum response across all scales at each voxel.
    """
    if volume.ndim != 3:
        raise ValueError(f"Expected 3D volume, got shape {volume.shape}")

    if sigmas is None:
        sigmas = [1.0, 2.0, 4.0, 8.0]

    use_gpu, device_id, device_name = _detect_device(device)
    logger.info(
        f"Computing Frangi vesselness: sigmas={sigmas}, "
        f"alpha={alpha}, beta={beta}, device={device_name}"
    )

    vol = volume.astype(np.float32)
    vesselness = np.zeros_like(vol)

    for i, sigma in enumerate(sigmas):
        logger.info(f"  Scale {i + 1}/{len(sigmas)}: sigma={sigma}")

        l1, l2, l3 = _hessian_eigenvalues_3d(vol, sigma, use_gpu=use_gpu, device_id=device_id)

        # Frangi criteria for bright tubular structures:
        # Vessels: l1 ~ 0, l2 << 0, l3 << 0
        # (for bright ridges; flip sign for dark ridges)
        if not black_ridges:
            # Reject where l2 > 0 or l3 > 0 (not tube-like)
            tube_mask = (l2 < 0) & (l3 < 0)
        else:
            tube_mask = (l2 > 0) & (l3 > 0)
            l1, l2, l3 = -l1, -l2, -l3

        # Ratios (Frangi 1998, eq. 11-13; eigenvalues sorted |l1| <= |l2| <= |l3|)
        # Ra = |l2| / |l3| — plate vs tube (near 1 for tubes, near 0 for plates)
        # Rb = |l1| / sqrt(|l2 * l3|) — deviation from blob (small for tubes)
        eps = 1e-10
        abs_l1 = np.abs(l1)
        abs_l2 = np.abs(l2)
        abs_l3 = np.abs(l3)

        Ra = abs_l2 / (abs_l3 + eps)
        Rb = abs_l1 / (np.sqrt(abs_l2 * abs_l3) + eps)

        # Frobenius norm (S) — structuredness
        S = np.sqrt(l1**2 + l2**2 + l3**2)

        if gamma is None:
            # Auto-set gamma to half max S for this scale
            gamma_scale = 0.5 * S.max() if S.max() > 0 else 1.0
        else:
            gamma_scale = gamma

        # Frangi vesselness function
        v = (
            (1.0 - np.exp(-(Ra**2) / (2.0 * alpha**2)))
            * np.exp(-(Rb**2) / (2.0 * beta**2))
            * (1.0 - np.exp(-(S**2) / (2.0 * gamma_scale**2)))
        )

        v[~tube_mask] = 0.0

        # Take maximum across scales
        vesselness = np.maximum(vesselness, v)

        del l1, l2, l3, Ra, Rb, S, v, tube_mask
        gc.collect()

    # Normalize to [0, 1]
    v_max = vesselness.max()
    if v_max > 0:
        vesselness /= v_max

    logger.info(f"Vesselness complete: range [{vesselness.min():.4f}, {vesselness.max():.4f}]")
    return vesselness


# =============================================================================
# Binarization and Morphological Cleanup
# =============================================================================


def binarize_vessel_mask(
    probability_map: np.ndarray,
    threshold: float | None = None,
    min_size: int = 500,
    closing_radius: int = 1,
    fill_holes_2d: bool = True,
) -> np.ndarray:
    """Convert a vesselness probability map to a clean binary vessel mask.

    Steps:
    1. Threshold (Otsu if not specified)
    2. Remove small connected components
    3. Morphological closing to bridge small gaps
    4. Fill internal holes per 2D plane

    Parameters
    ----------
    probability_map : np.ndarray
        3D vesselness map (float, range [0, 1]).
    threshold : float, optional
        Binarization threshold. If None, uses Otsu's method.
    min_size : int
        Minimum connected component volume in voxels. Default 500.
    closing_radius : int
        Ball radius for morphological closing. Default 1. Set 0 to skip.
    fill_holes_2d : bool
        Fill holes in each 2D plane. Default True.

    Returns
    -------
    np.ndarray
        3D binary mask (bool).
    """
    from skimage.filters import threshold_otsu
    from skimage.morphology import ball, closing, remove_small_objects

    if threshold is None:
        # Otsu on non-zero values for better threshold estimation
        nonzero = probability_map[probability_map > 0]
        if len(nonzero) > 0:
            threshold = threshold_otsu(nonzero)
        else:
            threshold = 0.1
        logger.info(f"Otsu threshold: {threshold:.4f}")
    else:
        logger.info(f"Using threshold: {threshold:.4f}")

    # Step 1: Threshold
    mask = probability_map > threshold
    logger.info(f"After threshold: {mask.sum():,} voxels ({100.0 * mask.sum() / mask.size:.2f}%)")

    # Step 2: Remove small objects
    # max_size removes objects with size <= threshold (replaces deprecated min_size)
    if min_size > 0:
        mask = remove_small_objects(mask, max_size=min_size - 1)
        logger.info(f"After small object removal (min={min_size}): {mask.sum():,} voxels")

    # Step 3: Morphological closing
    if closing_radius > 0:
        selem = ball(closing_radius)
        mask = closing(mask, footprint=selem)
        logger.info(f"After closing (radius={closing_radius}): {mask.sum():,} voxels")

    # Step 4: Fill holes per 2D plane (handles tubular cross-sections)
    if fill_holes_2d:
        for z in range(mask.shape[0]):
            mask[z] = binary_fill_holes(mask[z])
        logger.info(f"After hole filling: {mask.sum():,} voxels")

    return mask.astype(bool)


# =============================================================================
# Skeletonization
# =============================================================================


def skeletonize_vessels(
    binary_mask: np.ndarray,
) -> np.ndarray:
    """Reduce binary vessel mask to 1-voxel-wide centerlines.

    Uses the Lee (1994) 3D thinning algorithm via scikit-image.

    Parameters
    ----------
    binary_mask : np.ndarray
        3D binary mask (bool or uint8).

    Returns
    -------
    np.ndarray
        3D skeleton image (bool).
    """
    from skimage.morphology import skeletonize

    logger.info(f"Skeletonizing volume {binary_mask.shape}...")
    skeleton = skeletonize(binary_mask.astype(bool))
    n_skeleton = skeleton.sum()
    logger.info(f"Skeleton: {n_skeleton:,} voxels")
    return skeleton


def prune_skeleton(
    skeleton: np.ndarray,
    min_branch_length_um: float = 5.0,
    spacing: VesselSpacing | None = None,
) -> np.ndarray:
    """Remove short terminal branches (spurs) from a skeleton.

    Uses skan to identify endpoint-to-junction branches shorter than
    a threshold and removes them.

    Parameters
    ----------
    skeleton : np.ndarray
        3D skeleton image (bool).
    min_branch_length_um : float
        Minimum branch length in micrometers. Branches shorter than
        this are pruned. Default 5.0 um.
    spacing : VesselSpacing, optional
        Physical spacing. If None, assumes isotropic 1 um voxels.

    Returns
    -------
    np.ndarray
        Pruned skeleton (bool).
    """
    try:
        import skan
    except ImportError:
        logger.warning("skan not installed; skipping pruning. Install with: pip install skan")
        return skeleton

    if spacing is None:
        pixel_spacing = (1.0, 1.0, 1.0)
    else:
        pixel_spacing = spacing.isotropic_zyx

    skel_obj = skan.Skeleton(skeleton, spacing=pixel_spacing)
    summary = skan.summarize(skel_obj, find_main_branch=False)

    # Identify short terminal branches (endpoint-to-junction or endpoint-to-endpoint)
    short_mask = summary["branch-distance"] < min_branch_length_um
    terminal_mask = summary["branch-type"] != 2  # type 2 = junction-to-junction
    to_remove = short_mask & terminal_mask

    n_remove = to_remove.sum()
    logger.info(
        f"Pruning {n_remove}/{len(summary)} branches (< {min_branch_length_um} um, terminal)"
    )

    if n_remove == 0:
        return skeleton

    # Build pruned skeleton by removing short branch voxels
    pruned = skeleton.copy()
    for idx in summary.index[to_remove]:
        coords = skel_obj.path_coordinates(idx)
        for c in coords:
            c = tuple(int(x) for x in c)
            pruned[c] = False

    return pruned


# =============================================================================
# Graph Construction and Vessel Segment Analysis
# =============================================================================


def analyze_vessel_graph(
    skeleton: np.ndarray,
    binary_mask: np.ndarray,
    spacing: VesselSpacing | None = None,
) -> pd.DataFrame:
    """Convert skeleton to graph and extract per-segment morphometric features.

    Uses skan for graph construction and distance transform for radius
    estimation at each skeleton voxel.

    Parameters
    ----------
    skeleton : np.ndarray
        3D skeleton image (bool).
    binary_mask : np.ndarray
        3D binary vessel mask (same shape as skeleton).
    spacing : VesselSpacing, optional
        Physical spacing. If None, assumes isotropic 1 um voxels.

    Returns
    -------
    pd.DataFrame
        Per-segment features including:
        - branch_length_um: Physical branch length
        - euclidean_distance_um: Straight-line distance between endpoints
        - tortuosity: branch_length / euclidean_distance
        - branch_type: 0=endpoint-endpoint, 1=junction-endpoint, 2=junction-junction
        - mean_radius_um: Mean vessel radius along segment
        - median_radius_um: Median vessel radius
        - min_radius_um: Minimum vessel radius
        - max_radius_um: Maximum vessel radius
        - endpoint_0_z/y/x: Coordinates of first endpoint
        - endpoint_1_z/y/x: Coordinates of second endpoint
    """
    import pandas as pd

    try:
        import skan
    except ImportError:
        raise ImportError(
            "skan is required for vessel graph analysis. Install with: pip install skan"
        )

    if spacing is None:
        pixel_spacing = (1.0, 1.0, 1.0)
    else:
        pixel_spacing = spacing.isotropic_zyx

    logger.info("Building vessel graph with skan...")
    skel_obj = skan.Skeleton(skeleton, spacing=pixel_spacing)
    summary = skan.summarize(skel_obj, find_main_branch=False)

    # Rename columns for clarity
    df = summary.rename(
        columns={
            "branch-distance": "branch_length_um",
            "branch-type": "branch_type",
            "euclidean-distance": "euclidean_distance_um",
        }
    )

    # Compute tortuosity
    df["tortuosity"] = df["branch_length_um"] / df["euclidean_distance_um"].clip(lower=1e-6)

    # Compute distance transform for radius estimation
    logger.info("Computing distance transform for radius estimation...")
    dt = distance_transform_edt(binary_mask, sampling=pixel_spacing)

    # Sample radius along each branch
    radii_stats = []
    for idx in range(len(df)):
        coords = skel_obj.path_coordinates(idx)
        if len(coords) == 0:
            radii_stats.append((0, 0, 0, 0))
            continue

        # Sample distance transform at skeleton coordinates
        radii = []
        for c in coords:
            ci = tuple(int(round(x)) for x in c)
            # Bounds check
            if all(0 <= ci[d] < binary_mask.shape[d] for d in range(3)):
                radii.append(dt[ci])

        if radii:
            radii = np.array(radii)
            radii_stats.append(
                (
                    float(np.mean(radii)),
                    float(np.median(radii)),
                    float(np.min(radii)),
                    float(np.max(radii)),
                )
            )
        else:
            radii_stats.append((0, 0, 0, 0))

    radii_df = pd.DataFrame(
        radii_stats,
        columns=["mean_radius_um", "median_radius_um", "min_radius_um", "max_radius_um"],
        index=df.index,
    )
    df = pd.concat([df, radii_df], axis=1)

    # Add diameter columns
    df["mean_diameter_um"] = 2.0 * df["mean_radius_um"]

    logger.info(
        f"Graph analysis complete: {len(df)} segments, "
        f"mean length={df['branch_length_um'].mean():.1f} um, "
        f"mean diameter={df['mean_diameter_um'].mean():.1f} um"
    )

    return df


# =============================================================================
# Export
# =============================================================================


def export_vessel_results(
    result: VesselSegmentationResult,
    output_dir: str | Path,
    marker_name: str | None = None,
    save_vesselness: bool = False,
) -> dict[str, Path]:
    """Save all vessel segmentation outputs to disk.

    Parameters
    ----------
    result : VesselSegmentationResult
        Complete segmentation results.
    output_dir : str or Path
        Output directory (e.g., data/processed/vessel_3d/).
    marker_name : str, optional
        Marker name override. If None, uses result.marker_name.
    save_vesselness : bool
        Also save the vesselness map (large file). Default False.

    Returns
    -------
    dict[str, Path]
        Mapping of output names to file paths.
    """
    from skimage.io import imsave

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if marker_name is None:
        marker_name = result.marker_name

    outputs = {}

    # 1. Binary mask (isotropic)
    mask_path = output_dir / f"binary_mask_{marker_name}.tif"
    imsave(str(mask_path), result.binary_mask.astype(np.uint8) * 255, check_contrast=False)
    outputs["binary_mask"] = mask_path
    logger.info(f"Saved: {mask_path.name}")

    # 2. Binary mask (native resolution)
    if result.binary_mask_native is not None:
        native_path = output_dir / f"binary_mask_{marker_name}_native.tif"
        imsave(
            str(native_path),
            result.binary_mask_native.astype(np.uint8) * 255,
            check_contrast=False,
        )
        outputs["binary_mask_native"] = native_path
        logger.info(f"Saved: {native_path.name}")

    # 3. Skeleton
    if result.skeleton is not None:
        skel_path = output_dir / f"skeleton_{marker_name}.tif"
        imsave(str(skel_path), result.skeleton.astype(np.uint8) * 255, check_contrast=False)
        outputs["skeleton"] = skel_path
        logger.info(f"Saved: {skel_path.name}")

    # 4. Labeled segments
    if result.labeled_segments is not None:
        label_path = output_dir / f"labeled_segments_{marker_name}.tif"
        imsave(str(label_path), result.labeled_segments.astype(np.uint32), check_contrast=False)
        outputs["labeled_segments"] = label_path
        logger.info(f"Saved: {label_path.name}")

    # 5. Features CSV
    if result.features is not None:
        csv_path = output_dir / f"vessel_features_{marker_name}.csv"
        result.features.to_csv(csv_path, index=False)
        outputs["features"] = csv_path
        logger.info(f"Saved: {csv_path.name} ({len(result.features)} segments)")

    # 6. Graph (GraphML)
    if result.features is not None:
        try:
            graph_path = output_dir / f"vessel_graph_{marker_name}.graphml"
            _export_graph(result, graph_path)
            outputs["graph"] = graph_path
        except ImportError:
            logger.warning("networkx not installed; skipping graph export")

    # 7. 2D MIP of vessel mask (for Notebook 4 spatial analysis)
    mip_path = output_dir / f"vessel_2d_projection_{marker_name}.tif"
    mip = result.binary_mask.max(axis=0).astype(np.uint8) * 255
    imsave(str(mip_path), mip, check_contrast=False)
    outputs["mip_2d"] = mip_path
    logger.info(f"Saved: {mip_path.name}")

    # 8. Vesselness map (optional, large)
    if save_vesselness and result.vesselness is not None:
        vess_path = output_dir / f"vesselness_{marker_name}.tif"
        # Scale to uint16 for storage efficiency
        v_scaled = (result.vesselness * 65535).astype(np.uint16)
        imsave(str(vess_path), v_scaled, check_contrast=False)
        outputs["vesselness"] = vess_path
        logger.info(f"Saved: {vess_path.name}")

    return outputs


def _export_graph(result: VesselSegmentationResult, path: Path) -> None:
    """Export vessel graph as GraphML."""
    import networkx as nx

    G = nx.Graph()

    if result.features is None:
        return

    df = result.features
    for idx, row in df.iterrows():
        # Each branch is an edge in the graph
        node_a = f"node_{idx}_a"
        node_b = f"node_{idx}_b"
        G.add_node(node_a)
        G.add_node(node_b)
        G.add_edge(
            node_a,
            node_b,
            length_um=float(row.get("branch_length_um", 0)),
            tortuosity=float(row.get("tortuosity", 1)),
            mean_radius_um=float(row.get("mean_radius_um", 0)),
            mean_diameter_um=float(row.get("mean_diameter_um", 0)),
            branch_type=int(row.get("branch_type", 0)),
        )

    nx.write_graphml(G, str(path))
    logger.info(f"Saved: {path.name} ({G.number_of_nodes()} nodes, {G.number_of_edges()} edges)")


# =============================================================================
# High-Level Pipeline
# =============================================================================


def segment_vessels_3d(
    volume: np.ndarray,
    spacing: VesselSpacing | None = None,
    sigmas: list[float] | None = None,
    threshold: float | None = None,
    min_size: int = 500,
    closing_radius: int = 1,
    denoise_sigma: float = 0.5,
    make_isotropic: bool = True,
    prune_min_length_um: float = 5.0,
    skip_skeleton: bool = False,
    skip_graph: bool = False,
    device: str = "auto",
    marker_name: str = "vessel",
) -> VesselSegmentationResult:
    """Run the complete 3D vessel segmentation pipeline.

    This is the main entry point combining all steps:
    1. Preprocess (isotropic resampling + denoising)
    2. Frangi vesselness filter
    3. Binarization and morphological cleanup
    4. Skeletonization and pruning
    5. Graph construction and morphometry

    Parameters
    ----------
    volume : np.ndarray
        3D deconvolved z-stack (z, y, x), typically uint16.
    spacing : VesselSpacing, optional
        Physical voxel spacing. Default 377 nm XY, 1500 nm Z.
    sigmas : list of float, optional
        Frangi filter scales. Default [1, 2, 4, 8].
    threshold : float, optional
        Binarization threshold. None for Otsu auto-detection.
    min_size : int
        Minimum vessel volume in voxels. Default 500.
    closing_radius : int
        Morphological closing radius. Default 1.
    denoise_sigma : float
        Gaussian denoising sigma. Default 0.5. Set 0 to skip.
    make_isotropic : bool
        Resample z to match xy spacing. Default True.
    prune_min_length_um : float
        Minimum branch length for pruning. Default 5.0 um.
    skip_skeleton : bool
        Skip skeletonization and graph analysis. Default False.
    skip_graph : bool
        Skip graph analysis (still skeletonize). Default False.
    device : str
        'gpu', 'cpu', or 'auto'.
    marker_name : str
        Vessel marker name (e.g., "CD31"). Default "vessel".

    Returns
    -------
    VesselSegmentationResult
        Complete results including mask, skeleton, features.
    """
    if spacing is None:
        spacing = VesselSpacing()

    logger.info(
        f"Starting 3D vessel segmentation: {volume.shape}, "
        f"spacing=({spacing.z}, {spacing.xy}, {spacing.xy}) um, "
        f"marker={marker_name}"
    )

    # Step 1: Preprocess
    logger.info("Step 1/5: Preprocessing...")
    preprocessed = preprocess_volume(
        volume,
        spacing,
        denoise_sigma=denoise_sigma,
        make_isotropic=make_isotropic,
        device=device,
    )

    # Step 2: Frangi vesselness
    logger.info("Step 2/5: Computing Frangi vesselness...")
    vesselness = compute_vesselness_frangi(
        preprocessed,
        sigmas=sigmas,
        spacing=spacing,
        device=device,
    )

    # Step 3: Binarize and cleanup
    logger.info("Step 3/5: Binarization and cleanup...")
    mask = binarize_vessel_mask(
        vesselness,
        threshold=threshold,
        min_size=min_size,
        closing_radius=closing_radius,
    )

    # Create native-resolution mask by downsampling z back
    native_mask = None
    if make_isotropic:
        native_zoom = (1.0 / spacing.ratio, 1.0, 1.0)
        native_mask = ndimage.zoom(mask.astype(np.float32), native_zoom, order=0) > 0.5

    result = VesselSegmentationResult(
        binary_mask=mask,
        binary_mask_native=native_mask,
        vesselness=vesselness,
        spacing=spacing,
        marker_name=marker_name,
    )

    if skip_skeleton:
        logger.info("Skipping skeletonization (skip_skeleton=True)")
        return result

    # Step 4: Skeletonize
    logger.info("Step 4/5: Skeletonization...")
    skeleton = skeletonize_vessels(mask)
    skeleton = prune_skeleton(skeleton, min_branch_length_um=prune_min_length_um, spacing=spacing)
    result.skeleton = skeleton

    if skip_graph:
        logger.info("Skipping graph analysis (skip_graph=True)")
        return result

    # Step 5: Graph analysis
    logger.info("Step 5/5: Graph construction and morphometry...")
    features = analyze_vessel_graph(skeleton, mask, spacing=spacing)
    result.features = features

    logger.info("3D vessel segmentation complete.")
    return result
