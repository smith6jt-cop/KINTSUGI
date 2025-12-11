"""
Classical Segmentation Methods

Traditional segmentation approaches (watershed, thresholding)
for fluorescence microscopy when deep learning is not available.
"""

from __future__ import annotations

import logging
from typing import Optional, Tuple, Union

import numpy as np
from scipy import ndimage
from scipy.ndimage import distance_transform_edt

logger = logging.getLogger("kintsugi.segment.classical")


def threshold_otsu(image: np.ndarray) -> float:
    """
    Compute Otsu's threshold.

    Parameters
    ----------
    image : np.ndarray
        Input image

    Returns
    -------
    float
        Optimal threshold value
    """
    # Histogram
    if image.dtype == np.uint16:
        hist, bin_edges = np.histogram(image.ravel(), bins=256, range=(0, 65535))
    else:
        hist, bin_edges = np.histogram(image.ravel(), bins=256, range=(0, 255))

    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    hist = hist.astype(float)

    # Class probabilities
    weight1 = np.cumsum(hist)
    weight2 = np.cumsum(hist[::-1])[::-1]

    # Class means
    mean1 = np.cumsum(hist * bin_centers) / (weight1 + 1e-8)
    mean2 = (np.cumsum((hist * bin_centers)[::-1]) / (weight2[::-1] + 1e-8))[::-1]

    # Between-class variance
    variance = weight1[:-1] * weight2[1:] * (mean1[:-1] - mean2[1:]) ** 2

    idx = np.argmax(variance)
    threshold = bin_centers[idx]

    return float(threshold)


def threshold_adaptive(
    image: np.ndarray,
    block_size: int = 51,
    offset: float = 0,
    method: str = "mean",
) -> np.ndarray:
    """
    Adaptive thresholding.

    Parameters
    ----------
    image : np.ndarray
        Input image
    block_size : int
        Size of local neighborhood
    offset : float
        Constant subtracted from local mean/median
    method : str
        "mean" or "median"

    Returns
    -------
    np.ndarray
        Binary mask
    """
    img = image.astype(np.float64)

    if method == "mean":
        local_thresh = ndimage.uniform_filter(img, size=block_size)
    elif method == "median":
        local_thresh = ndimage.median_filter(img, size=block_size)
    else:
        raise ValueError(f"Unknown method: {method}")

    return img > (local_thresh - offset)


def segment_nuclei_watershed(
    image: np.ndarray,
    min_distance: int = 10,
    threshold: Optional[float] = None,
    min_size: int = 50,
    max_size: int = 5000,
) -> np.ndarray:
    """
    Segment nuclei using watershed algorithm.

    Parameters
    ----------
    image : np.ndarray
        Nuclear stain image (e.g., DAPI)
    min_distance : int
        Minimum distance between nuclei centers
    threshold : float, optional
        Threshold for binarization (Otsu if None)
    min_size : int
        Minimum nucleus size in pixels
    max_size : int
        Maximum nucleus size in pixels

    Returns
    -------
    np.ndarray
        Labeled segmentation mask
    """
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed
    from skimage import morphology

    img = image.astype(np.float64)

    # Normalize
    p1, p99 = np.percentile(img, [1, 99])
    if p99 > p1:
        img = (img - p1) / (p99 - p1)
        img = np.clip(img, 0, 1)

    # Threshold
    if threshold is None:
        threshold = threshold_otsu(img)
    binary = img > threshold

    # Clean up binary mask
    binary = morphology.remove_small_objects(binary, min_size=min_size // 4)
    binary = morphology.binary_closing(binary, morphology.disk(2))
    binary = morphology.binary_opening(binary, morphology.disk(1))

    # Distance transform
    distance = distance_transform_edt(binary)

    # Find peaks (nuclei centers)
    local_max_coords = peak_local_max(
        distance,
        min_distance=min_distance,
        labels=binary,
    )

    # Create markers from peaks
    markers = np.zeros_like(image, dtype=np.int32)
    for i, (y, x) in enumerate(local_max_coords, start=1):
        markers[y, x] = i

    markers = ndimage.label(morphology.binary_dilation(markers > 0, morphology.disk(2)))[0]

    # Watershed
    labels = watershed(-distance, markers, mask=binary)

    # Filter by size
    labels = _filter_labels_by_size(labels, min_size, max_size)

    return labels


def segment_cells_watershed(
    membrane_image: np.ndarray,
    nuclei_mask: Optional[np.ndarray] = None,
    expansion: int = 20,
    min_size: int = 100,
    max_size: int = 10000,
) -> np.ndarray:
    """
    Segment cells using watershed with nuclei as seeds.

    Parameters
    ----------
    membrane_image : np.ndarray
        Membrane marker image
    nuclei_mask : np.ndarray, optional
        Pre-computed nuclei mask to use as seeds
    expansion : int
        Maximum expansion distance from nuclei
    min_size : int
        Minimum cell size
    max_size : int
        Maximum cell size

    Returns
    -------
    np.ndarray
        Labeled cell segmentation mask
    """
    from skimage.segmentation import watershed
    from skimage import morphology

    membrane = membrane_image.astype(np.float64)

    # Normalize membrane image
    p1, p99 = np.percentile(membrane, [1, 99])
    if p99 > p1:
        membrane = (membrane - p1) / (p99 - p1)
        membrane = np.clip(membrane, 0, 1)

    # Create gradient image for watershed
    gradient = ndimage.sobel(membrane)
    gradient = (gradient - gradient.min()) / (gradient.max() - gradient.min() + 1e-8)

    # Get or compute nuclei markers
    if nuclei_mask is not None:
        markers = nuclei_mask
    else:
        # Simple thresholding for nuclei
        markers = segment_nuclei_watershed(membrane_image)

    # Ensure markers are labeled
    if markers.max() == 1:
        markers = ndimage.label(markers)[0]

    # Create expansion mask
    nuclei_binary = markers > 0
    expanded = morphology.binary_dilation(nuclei_binary, morphology.disk(expansion))

    # Watershed
    labels = watershed(gradient, markers, mask=expanded)

    # Filter by size
    labels = _filter_labels_by_size(labels, min_size, max_size)

    return labels


def _filter_labels_by_size(
    labels: np.ndarray,
    min_size: int,
    max_size: int,
) -> np.ndarray:
    """Filter labeled objects by size."""
    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    new_label = 1
    for label in unique_labels:
        mask = labels == label
        size = np.sum(mask)
        if min_size <= size <= max_size:
            output[mask] = new_label
            new_label += 1

    return output


def segment_by_intensity(
    image: np.ndarray,
    low_threshold: Optional[float] = None,
    high_threshold: Optional[float] = None,
    min_size: int = 50,
) -> np.ndarray:
    """
    Simple intensity-based segmentation.

    Parameters
    ----------
    image : np.ndarray
        Input image
    low_threshold : float, optional
        Lower intensity threshold (percentile if None)
    high_threshold : float, optional
        Upper intensity threshold
    min_size : int
        Minimum object size

    Returns
    -------
    np.ndarray
        Labeled segmentation
    """
    from skimage import morphology

    if low_threshold is None:
        low_threshold = np.percentile(image, 50)

    binary = image > low_threshold

    if high_threshold is not None:
        binary &= image < high_threshold

    # Clean up
    binary = morphology.remove_small_objects(binary, min_size=min_size)
    binary = morphology.binary_closing(binary, morphology.disk(2))

    # Label
    labels = ndimage.label(binary)[0]

    return labels


def voronoi_expand(
    seeds: np.ndarray,
    boundary_image: Optional[np.ndarray] = None,
    max_distance: int = 50,
) -> np.ndarray:
    """
    Expand seeds using Voronoi tessellation.

    Parameters
    ----------
    seeds : np.ndarray
        Labeled seed mask (e.g., nuclei)
    boundary_image : np.ndarray, optional
        Image to use for boundary detection
    max_distance : int
        Maximum expansion distance

    Returns
    -------
    np.ndarray
        Expanded labels
    """
    from skimage.segmentation import watershed
    from skimage import morphology

    if seeds.max() == 0:
        return seeds

    # Ensure seeds are labeled
    if seeds.max() == 1:
        seeds = ndimage.label(seeds)[0]

    # Create expansion mask
    seed_binary = seeds > 0
    expansion_mask = morphology.binary_dilation(seed_binary, morphology.disk(max_distance))

    # Use boundary image for watershed if provided
    if boundary_image is not None:
        gradient = ndimage.sobel(boundary_image.astype(float))
        expanded = watershed(gradient, seeds, mask=expansion_mask)
    else:
        # Simple distance-based expansion
        distance = distance_transform_edt(~seed_binary)
        expanded = watershed(distance, seeds, mask=expansion_mask)

    return expanded
