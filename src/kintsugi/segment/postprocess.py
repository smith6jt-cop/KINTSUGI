"""
Segmentation Post-Processing

Utilities for refining and cleaning up segmentation masks.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy import ndimage

logger = logging.getLogger("kintsugi.segment.postprocess")


def refine_masks(
    labels: np.ndarray,
    min_size: int = 50,
    max_size: int = 5000,
    fill_holes: bool = True,
    smooth: bool = True,
    split_touching: bool = False,
) -> np.ndarray:
    """
    Apply comprehensive refinement to segmentation masks.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    min_size : int
        Minimum object size to keep
    max_size : int
        Maximum object size to keep
    fill_holes : bool
        Fill holes in objects
    smooth : bool
        Smooth object boundaries
    split_touching : bool
        Try to split touching objects

    Returns
    -------
    np.ndarray
        Refined labeled mask
    """
    # Filter by size
    refined = filter_masks_by_size(labels, min_size, max_size)

    # Fill holes
    if fill_holes:
        refined = fill_holes_in_labels(refined)

    # Smooth boundaries
    if smooth:
        refined = smooth_boundaries(refined)

    # Split touching objects
    if split_touching:
        refined = split_touching_objects(refined)

    # Relabel consecutively
    refined = relabel_sequential(refined)

    return refined


def filter_masks_by_size(
    labels: np.ndarray,
    min_size: int = 50,
    max_size: int = 5000,
) -> np.ndarray:
    """
    Filter labeled objects by size.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    min_size : int
        Minimum object size
    max_size : int
        Maximum object size

    Returns
    -------
    np.ndarray
        Filtered labeled mask
    """
    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    for label in unique_labels:
        mask = labels == label
        size = np.sum(mask)
        if min_size <= size <= max_size:
            output[mask] = label

    return output


def split_touching_objects(
    labels: np.ndarray,
    min_distance: int = 5,
) -> np.ndarray:
    """
    Split touching objects using distance transform watershed.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    min_distance : int
        Minimum distance between object centers

    Returns
    -------
    np.ndarray
        Labels with touching objects split
    """
    from scipy.ndimage import distance_transform_edt
    from skimage.feature import peak_local_max
    from skimage.segmentation import watershed

    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    new_label = 1

    for label in unique_labels:
        mask = labels == label

        # Distance transform
        distance = distance_transform_edt(mask)

        # Find local maxima
        coords = peak_local_max(
            distance,
            min_distance=min_distance,
            labels=mask,
        )

        if len(coords) <= 1:
            # Single object, no splitting needed
            output[mask] = new_label
            new_label += 1
            continue

        # Create markers
        markers = np.zeros_like(mask, dtype=np.int32)
        for i, (y, x) in enumerate(coords, start=1):
            markers[y, x] = i

        # Watershed to split
        split_labels = watershed(-distance, markers, mask=mask)

        # Add to output
        for split_label in np.unique(split_labels):
            if split_label > 0:
                output[split_labels == split_label] = new_label
                new_label += 1

    return output


def fill_holes(mask: np.ndarray) -> np.ndarray:
    """
    Fill holes in a binary mask.

    Parameters
    ----------
    mask : np.ndarray
        Binary mask

    Returns
    -------
    np.ndarray
        Mask with holes filled
    """
    return ndimage.binary_fill_holes(mask)


def fill_holes_in_labels(labels: np.ndarray) -> np.ndarray:
    """
    Fill holes in each labeled object.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask

    Returns
    -------
    np.ndarray
        Labels with holes filled
    """
    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    for label in unique_labels:
        mask = labels == label
        filled = ndimage.binary_fill_holes(mask)
        output[filled] = label

    return output


def smooth_boundaries(
    labels: np.ndarray,
    iterations: int = 2,
) -> np.ndarray:
    """
    Smooth object boundaries using morphological operations.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    iterations : int
        Number of smoothing iterations

    Returns
    -------
    np.ndarray
        Labels with smoothed boundaries
    """
    from skimage import morphology

    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    for label in unique_labels:
        mask = labels == label

        # Morphological closing then opening
        for _ in range(iterations):
            mask = morphology.binary_closing(mask, morphology.disk(1))
            mask = morphology.binary_opening(mask, morphology.disk(1))

        output[mask] = label

    return output


def relabel_sequential(labels: np.ndarray) -> np.ndarray:
    """
    Relabel objects with sequential integers starting from 1.

    Parameters
    ----------
    labels : np.ndarray
        Labeled mask with potentially non-sequential labels

    Returns
    -------
    np.ndarray
        Sequentially labeled mask
    """
    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    for new_label, old_label in enumerate(unique_labels, start=1):
        output[labels == old_label] = new_label

    return output


def remove_edge_objects(
    labels: np.ndarray,
    margin: int = 1,
) -> np.ndarray:
    """
    Remove objects touching the image edge.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    margin : int
        Margin from edge

    Returns
    -------
    np.ndarray
        Labels with edge objects removed
    """
    h, w = labels.shape

    # Find labels touching edges
    edge_labels = set()
    edge_labels.update(np.unique(labels[:margin, :]))
    edge_labels.update(np.unique(labels[-margin:, :]))
    edge_labels.update(np.unique(labels[:, :margin]))
    edge_labels.update(np.unique(labels[:, -margin:]))
    edge_labels.discard(0)

    # Remove edge labels
    output = labels.copy()
    for label in edge_labels:
        output[labels == label] = 0

    return output


def expand_labels(
    labels: np.ndarray,
    distance: int = 5,
) -> np.ndarray:
    """
    Expand labels by a fixed distance.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    distance : int
        Expansion distance in pixels

    Returns
    -------
    np.ndarray
        Expanded labels
    """
    from skimage.segmentation import expand_labels as skimage_expand

    try:
        return skimage_expand(labels, distance=distance)
    except ImportError:
        # Fallback implementation
        from skimage import morphology

        output = np.zeros_like(labels)
        unique_labels = np.unique(labels)
        unique_labels = unique_labels[unique_labels > 0]

        for label in unique_labels:
            mask = labels == label
            expanded = morphology.binary_dilation(mask, morphology.disk(distance))
            # Only expand into unlabeled regions
            output[(expanded) & (output == 0)] = label

        return output


def erode_labels(
    labels: np.ndarray,
    distance: int = 2,
) -> np.ndarray:
    """
    Erode labels by a fixed distance.

    Parameters
    ----------
    labels : np.ndarray
        Labeled segmentation mask
    distance : int
        Erosion distance in pixels

    Returns
    -------
    np.ndarray
        Eroded labels
    """
    from skimage import morphology

    output = np.zeros_like(labels)
    unique_labels = np.unique(labels)
    unique_labels = unique_labels[unique_labels > 0]

    for label in unique_labels:
        mask = labels == label
        eroded = morphology.binary_erosion(mask, morphology.disk(distance))
        if np.any(eroded):
            output[eroded] = label

    return output


def compute_overlap(
    labels1: np.ndarray,
    labels2: np.ndarray,
) -> np.ndarray:
    """
    Compute overlap matrix between two labeled images.

    Parameters
    ----------
    labels1 : np.ndarray
        First labeled mask
    labels2 : np.ndarray
        Second labeled mask

    Returns
    -------
    np.ndarray
        Overlap matrix where entry (i,j) is the overlap between
        label i from labels1 and label j from labels2
    """
    max_label1 = labels1.max()
    max_label2 = labels2.max()

    overlap = np.zeros((max_label1 + 1, max_label2 + 1), dtype=np.int32)

    for i in range(labels1.shape[0]):
        for j in range(labels1.shape[1]):
            l1 = labels1[i, j]
            l2 = labels2[i, j]
            overlap[l1, l2] += 1

    return overlap


def compute_iou_matrix(
    labels1: np.ndarray,
    labels2: np.ndarray,
) -> np.ndarray:
    """
    Compute IoU (Intersection over Union) matrix between two labeled images.

    Parameters
    ----------
    labels1 : np.ndarray
        First labeled mask
    labels2 : np.ndarray
        Second labeled mask

    Returns
    -------
    np.ndarray
        IoU matrix
    """
    overlap = compute_overlap(labels1, labels2)

    # Compute areas
    areas1 = np.bincount(labels1.ravel())
    areas2 = np.bincount(labels2.ravel())

    # Pad if necessary
    if len(areas1) < overlap.shape[0]:
        areas1 = np.pad(areas1, (0, overlap.shape[0] - len(areas1)))
    if len(areas2) < overlap.shape[1]:
        areas2 = np.pad(areas2, (0, overlap.shape[1] - len(areas2)))

    # IoU = intersection / union = intersection / (area1 + area2 - intersection)
    union = areas1[:, np.newaxis] + areas2[np.newaxis, :] - overlap
    iou = np.zeros_like(overlap, dtype=np.float64)
    valid = union > 0
    iou[valid] = overlap[valid] / union[valid]

    return iou
