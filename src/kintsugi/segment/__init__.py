"""
KINTSUGI Segmentation Module

Pure Python wrappers for segmentation models optimized for
fluorescence microscopy images.

Includes:
- SAM (Segment Anything Model) wrapper
- InstanSeg integration
- Classical watershed segmentation
- Post-processing utilities

Usage:
    from kintsugi.segment import (
        SAMSegmenter,
        segment_nuclei,
        segment_cells,
        refine_masks,
    )

    # SAM-based segmentation
    segmenter = SAMSegmenter()
    masks = segmenter.segment(image, prompt_type="auto")

    # Quick nuclei segmentation
    nuclei_masks = segment_nuclei(dapi_image)
"""

from kintsugi.segment.classical import (
    segment_cells_watershed,
    segment_nuclei_watershed,
    threshold_adaptive,
    threshold_otsu,
)
from kintsugi.segment.postprocess import (
    fill_holes,
    filter_masks_by_size,
    refine_masks,
    smooth_boundaries,
    split_touching_objects,
)
from kintsugi.segment.sam_wrapper import (
    SAMConfig,
    SAMSegmenter,
    segment_with_sam,
)

__all__ = [
    # SAM
    "SAMSegmenter",
    "SAMConfig",
    "segment_with_sam",
    # Classical
    "segment_nuclei_watershed",
    "segment_cells_watershed",
    "threshold_otsu",
    "threshold_adaptive",
    # Post-processing
    "refine_masks",
    "filter_masks_by_size",
    "split_touching_objects",
    "fill_holes",
    "smooth_boundaries",
]


# Convenience functions
def segment_nuclei(
    image,
    method: str = "auto",
    **kwargs,
):
    """
    Segment nuclei from a nuclear stain image.

    Parameters
    ----------
    image : array-like
        Nuclear stain image (e.g., DAPI)
    method : str
        Segmentation method: "auto", "sam", "watershed", "instanseg"
    **kwargs
        Additional arguments passed to segmenter

    Returns
    -------
    np.ndarray
        Labeled segmentation mask
    """
    if method == "sam":
        segmenter = SAMSegmenter()
        return segmenter.segment(image, **kwargs)

    elif method == "watershed":
        return segment_nuclei_watershed(image, **kwargs)

    elif method == "instanseg":
        try:
            from kintsugi.segment.instanseg_wrapper import segment_instanseg

            return segment_instanseg(image, model="nuclei", **kwargs)
        except ImportError:
            raise ImportError("InstanSeg not available. Install with: pip install instanseg")

    else:  # auto
        # Try SAM first, fall back to watershed
        try:
            segmenter = SAMSegmenter()
            return segmenter.segment(image, **kwargs)
        except Exception:
            return segment_nuclei_watershed(image, **kwargs)


def segment_cells(
    membrane_image,
    nuclei_mask=None,
    method: str = "auto",
    **kwargs,
):
    """
    Segment cells from a membrane marker image.

    Parameters
    ----------
    membrane_image : array-like
        Membrane marker image
    nuclei_mask : array-like, optional
        Pre-computed nuclei mask as seeds
    method : str
        Segmentation method: "auto", "sam", "watershed"
    **kwargs
        Additional arguments

    Returns
    -------
    np.ndarray
        Labeled cell segmentation mask
    """
    if method == "watershed":
        return segment_cells_watershed(membrane_image, nuclei_mask, **kwargs)

    elif method == "sam":
        segmenter = SAMSegmenter()
        return segmenter.segment(membrane_image, **kwargs)

    else:  # auto
        try:
            segmenter = SAMSegmenter()
            return segmenter.segment(membrane_image, **kwargs)
        except Exception:
            return segment_cells_watershed(membrane_image, nuclei_mask, **kwargs)
