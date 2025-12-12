"""
SAM (Segment Anything Model) Wrapper

Pure Python wrapper for Meta's Segment Anything Model optimized
for fluorescence microscopy cell/nuclei segmentation.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.request import urlretrieve

import numpy as np

logger = logging.getLogger("kintsugi.segment.sam")

# Model URLs and checksums
SAM_MODELS = {
    "vit_h": {
        "url": "https://dl.fbaipublicfiles.com/segment_anything/sam_vit_h_4b8939.pth",
        "checkpoint": "sam_vit_h_4b8939.pth",
    },
    "vit_l": {
        "url": "https://dl.fbaipublicfiles.com/segment_anything/sam_vit_l_0b3195.pth",
        "checkpoint": "sam_vit_l_0b3195.pth",
    },
    "vit_b": {
        "url": "https://dl.fbaipublicfiles.com/segment_anything/sam_vit_b_01ec64.pth",
        "checkpoint": "sam_vit_b_01ec64.pth",
    },
}


@dataclass
class SAMConfig:
    """Configuration for SAM segmentation."""

    # Model selection
    model_type: str = "vit_b"  # vit_h, vit_l, vit_b

    # Processing
    points_per_side: int = 32  # For automatic mask generation
    pred_iou_thresh: float = 0.88
    stability_score_thresh: float = 0.95
    crop_n_layers: int = 0
    min_mask_region_area: int = 100

    # Image preprocessing for microscopy
    normalize_input: bool = True
    enhance_contrast: bool = True
    invert_image: bool = False  # Invert for dark objects on light background

    # Device
    device: str = "auto"

    # Tiling for large images
    tile_size: int = 1024
    tile_overlap: int = 256


class SAMSegmenter:
    """
    SAM-based segmentation for fluorescence microscopy.

    Wraps Meta's Segment Anything Model with preprocessing
    optimized for fluorescence microscopy images.

    Example
    -------
    >>> segmenter = SAMSegmenter()
    >>> masks = segmenter.segment(dapi_image, prompt_type="auto")
    >>> labeled = segmenter.masks_to_labels(masks)
    """

    def __init__(
        self,
        config: SAMConfig | None = None,
        model_path: str | Path | None = None,
    ):
        """
        Initialize SAM segmenter.

        Parameters
        ----------
        config : SAMConfig, optional
            Configuration for segmentation
        model_path : str or Path, optional
            Path to model checkpoint (downloaded if not present)
        """
        self.config = config or SAMConfig()
        self.model_path = model_path
        self.model = None
        self.predictor = None
        self.mask_generator = None
        self.device = None

    def _setup_device(self):
        """Set up PyTorch device with multi-GPU support."""
        try:
            import torch

            from kintsugi.gpu import get_gpu_manager

            self._gpu_manager = get_gpu_manager()

            if self.config.device == "auto":
                self.device = self._gpu_manager.get_torch_device()
            else:
                self.device = torch.device(self.config.device)

            if self._gpu_manager.device_count > 1:
                logger.info(
                    f"SAM using {self._gpu_manager.device_count} GPUs: "
                    f"{self._gpu_manager.device_ids}"
                )
            else:
                logger.info(f"SAM using device: {self.device}")

        except ImportError:
            raise ImportError("PyTorch is required for SAM. Install with: pip install torch")

    def _load_model(self):
        """Load SAM model."""
        try:
            from segment_anything import SamAutomaticMaskGenerator, SamPredictor, sam_model_registry
        except ImportError:
            raise ImportError(
                "segment_anything is required. Install with: pip install segment-anything"
            )

        self._setup_device()

        # Get model checkpoint
        if self.model_path is None:
            self.model_path = self._download_model()

        model_path = Path(self.model_path)
        if not model_path.exists():
            raise FileNotFoundError(f"Model not found: {model_path}")

        # Load model
        logger.info(f"Loading SAM model: {self.config.model_type}")
        self.model = sam_model_registry[self.config.model_type](checkpoint=str(model_path))

        # Use GPUManager for multi-GPU wrapping
        if hasattr(self, "_gpu_manager"):
            self.model = self._gpu_manager.wrap_model(self.model)
        else:
            self.model.to(self.device)

        # Create predictor and mask generator
        self.predictor = SamPredictor(self.model)
        self.mask_generator = SamAutomaticMaskGenerator(
            self.model,
            points_per_side=self.config.points_per_side,
            pred_iou_thresh=self.config.pred_iou_thresh,
            stability_score_thresh=self.config.stability_score_thresh,
            crop_n_layers=self.config.crop_n_layers,
            min_mask_region_area=self.config.min_mask_region_area,
        )

        logger.info("SAM model loaded successfully")

    def _download_model(self) -> Path:
        """Download SAM model checkpoint."""
        cache_dir = Path.home() / ".kintsugi" / "models" / "sam"
        cache_dir.mkdir(parents=True, exist_ok=True)

        model_info = SAM_MODELS.get(self.config.model_type)
        if model_info is None:
            raise ValueError(f"Unknown model type: {self.config.model_type}")

        checkpoint_path = cache_dir / model_info["checkpoint"]

        if not checkpoint_path.exists():
            logger.info(f"Downloading SAM {self.config.model_type} model...")
            urlretrieve(model_info["url"], checkpoint_path)
            logger.info("Download complete")

        return checkpoint_path

    def _preprocess_image(self, image: np.ndarray) -> np.ndarray:
        """Preprocess microscopy image for SAM."""
        img = image.astype(np.float64)

        # Normalize
        if self.config.normalize_input:
            p1, p99 = np.percentile(img, [1, 99])
            if p99 > p1:
                img = (img - p1) / (p99 - p1)
                img = np.clip(img, 0, 1)

        # Enhance contrast if needed
        if self.config.enhance_contrast:
            # Simple contrast stretching
            img = img**0.8  # Gamma correction

        # Invert if needed
        if self.config.invert_image:
            img = 1.0 - img

        # Convert to 8-bit RGB (SAM expects RGB images)
        img_uint8 = (img * 255).astype(np.uint8)
        if img_uint8.ndim == 2:
            img_rgb = np.stack([img_uint8] * 3, axis=-1)
        else:
            img_rgb = img_uint8

        return img_rgb

    def segment(
        self,
        image: np.ndarray,
        prompt_type: str = "auto",
        point_coords: np.ndarray | None = None,
        point_labels: np.ndarray | None = None,
        box: np.ndarray | None = None,
    ) -> list[dict[str, Any]]:
        """
        Segment image using SAM.

        Parameters
        ----------
        image : np.ndarray
            Input image (2D grayscale or RGB)
        prompt_type : str
            Type of prompting:
            - "auto": Automatic mask generation
            - "points": Use point prompts
            - "box": Use bounding box prompt
        point_coords : np.ndarray, optional
            Point coordinates for point prompts (N x 2)
        point_labels : np.ndarray, optional
            Point labels (1=foreground, 0=background)
        box : np.ndarray, optional
            Bounding box [x1, y1, x2, y2]

        Returns
        -------
        list
            List of mask dictionaries with keys:
            - 'segmentation': boolean mask
            - 'area': mask area
            - 'predicted_iou': confidence score
            - 'bbox': bounding box
        """
        if self.model is None:
            self._load_model()

        # Preprocess
        img_rgb = self._preprocess_image(image)

        # Handle large images with tiling
        if image.shape[0] > self.config.tile_size or image.shape[1] > self.config.tile_size:
            return self._segment_tiled(img_rgb, prompt_type)

        if prompt_type == "auto":
            masks = self.mask_generator.generate(img_rgb)
            return masks

        elif prompt_type == "points":
            if point_coords is None:
                raise ValueError("point_coords required for point prompts")
            if point_labels is None:
                point_labels = np.ones(len(point_coords), dtype=int)

            self.predictor.set_image(img_rgb)
            masks, scores, _ = self.predictor.predict(
                point_coords=point_coords,
                point_labels=point_labels,
                multimask_output=True,
            )

            # Convert to list of dicts
            return self._masks_to_dict_list(masks, scores)

        elif prompt_type == "box":
            if box is None:
                raise ValueError("box required for box prompts")

            self.predictor.set_image(img_rgb)
            masks, scores, _ = self.predictor.predict(
                box=box,
                multimask_output=True,
            )

            return self._masks_to_dict_list(masks, scores)

        else:
            raise ValueError(f"Unknown prompt_type: {prompt_type}")

    def _segment_tiled(
        self,
        image: np.ndarray,
        prompt_type: str,
    ) -> list[dict[str, Any]]:
        """Segment large image using tiling."""
        h, w = image.shape[:2]
        tile_size = self.config.tile_size
        overlap = self.config.tile_overlap
        step = tile_size - overlap

        all_masks = []

        for y in range(0, h, step):
            for x in range(0, w, step):
                # Extract tile
                y_end = min(y + tile_size, h)
                x_end = min(x + tile_size, w)
                tile = image[y:y_end, x:x_end]

                # Pad if needed
                if tile.shape[0] < tile_size or tile.shape[1] < tile_size:
                    padded = np.zeros((tile_size, tile_size, 3), dtype=tile.dtype)
                    padded[: tile.shape[0], : tile.shape[1]] = tile
                    tile = padded

                # Segment tile
                if prompt_type == "auto":
                    tile_masks = self.mask_generator.generate(tile)
                else:
                    continue  # Skip other prompt types for tiling

                # Offset masks to global coordinates
                for mask in tile_masks:
                    # Create full-size mask
                    full_mask = np.zeros((h, w), dtype=bool)
                    seg = mask["segmentation"][: y_end - y, : x_end - x]
                    full_mask[y:y_end, x:x_end] = seg

                    mask["segmentation"] = full_mask
                    mask["bbox"] = [
                        mask["bbox"][0] + x,
                        mask["bbox"][1] + y,
                        mask["bbox"][2],
                        mask["bbox"][3],
                    ]

                all_masks.extend(tile_masks)

        # Merge overlapping masks
        all_masks = self._merge_overlapping_masks(all_masks)

        return all_masks

    def _merge_overlapping_masks(
        self,
        masks: list[dict[str, Any]],
        iou_threshold: float = 0.5,
    ) -> list[dict[str, Any]]:
        """Merge overlapping masks from different tiles."""
        if len(masks) <= 1:
            return masks

        # Sort by area (descending)
        masks = sorted(masks, key=lambda x: x["area"], reverse=True)

        merged = []
        used = set()

        for i, mask1 in enumerate(masks):
            if i in used:
                continue

            current_mask = mask1["segmentation"].copy()
            current_area = mask1["area"]

            # Find overlapping masks
            for j, mask2 in enumerate(masks[i + 1 :], start=i + 1):
                if j in used:
                    continue

                # Check IoU
                intersection = np.sum(current_mask & mask2["segmentation"])
                union = np.sum(current_mask | mask2["segmentation"])

                if union > 0 and intersection / union > iou_threshold:
                    # Merge masks
                    current_mask = current_mask | mask2["segmentation"]
                    current_area = np.sum(current_mask)
                    used.add(j)

            merged.append(
                {
                    "segmentation": current_mask,
                    "area": int(current_area),
                    "predicted_iou": mask1["predicted_iou"],
                    "bbox": self._compute_bbox(current_mask),
                }
            )

        return merged

    def _masks_to_dict_list(
        self,
        masks: np.ndarray,
        scores: np.ndarray,
    ) -> list[dict[str, Any]]:
        """Convert mask array to list of dicts."""
        results = []
        for i in range(len(masks)):
            mask = masks[i]
            results.append(
                {
                    "segmentation": mask,
                    "area": int(np.sum(mask)),
                    "predicted_iou": float(scores[i]),
                    "bbox": self._compute_bbox(mask),
                }
            )
        return results

    def _compute_bbox(self, mask: np.ndarray) -> list[int]:
        """Compute bounding box from mask."""
        rows = np.any(mask, axis=1)
        cols = np.any(mask, axis=0)

        if not np.any(rows) or not np.any(cols):
            return [0, 0, 0, 0]

        y1, y2 = np.where(rows)[0][[0, -1]]
        x1, x2 = np.where(cols)[0][[0, -1]]

        return [int(x1), int(y1), int(x2 - x1 + 1), int(y2 - y1 + 1)]

    def masks_to_labels(
        self,
        masks: list[dict[str, Any]],
        image_shape: tuple[int, int] | None = None,
    ) -> np.ndarray:
        """
        Convert list of masks to labeled image.

        Parameters
        ----------
        masks : list
            List of mask dictionaries
        image_shape : tuple, optional
            Output shape (inferred from masks if not provided)

        Returns
        -------
        np.ndarray
            Labeled image where each object has a unique integer
        """
        if not masks:
            return np.zeros(image_shape or (1, 1), dtype=np.int32)

        if image_shape is None:
            # Infer from masks
            h = max(m["segmentation"].shape[0] for m in masks)
            w = max(m["segmentation"].shape[1] for m in masks)
            image_shape = (h, w)

        labels = np.zeros(image_shape, dtype=np.int32)

        # Sort by area (largest first) to handle overlaps
        masks = sorted(masks, key=lambda x: x["area"], reverse=True)

        for i, mask in enumerate(masks, start=1):
            seg = mask["segmentation"]
            # Handle size mismatch
            if seg.shape != image_shape:
                full_seg = np.zeros(image_shape, dtype=bool)
                h, w = min(seg.shape[0], image_shape[0]), min(seg.shape[1], image_shape[1])
                full_seg[:h, :w] = seg[:h, :w]
                seg = full_seg

            # Only label unlabeled pixels
            labels[seg & (labels == 0)] = i

        return labels


def segment_with_sam(
    image: np.ndarray,
    model_type: str = "vit_b",
    return_labels: bool = True,
    **kwargs,
) -> np.ndarray | list[dict[str, Any]]:
    """
    Convenience function for SAM segmentation.

    Parameters
    ----------
    image : np.ndarray
        Input image
    model_type : str
        SAM model type: "vit_h", "vit_l", "vit_b"
    return_labels : bool
        If True, return labeled image; if False, return mask list
    **kwargs
        Additional arguments passed to segment()

    Returns
    -------
    np.ndarray or list
        Labeled image or list of masks
    """
    config = SAMConfig(model_type=model_type)
    segmenter = SAMSegmenter(config)
    masks = segmenter.segment(image, **kwargs)

    if return_labels:
        return segmenter.masks_to_labels(masks, image.shape[:2])
    return masks
