"""
Channel Quality Assessment using Deep Learning

This module provides automated channel quality assessment for multiplex imaging,
reducing manual review workload by automatically identifying problematic channels.
"""

import numpy as np
import torch
import zarr
from typing import Dict, List, Tuple, Optional, Union
from dataclasses import dataclass, asdict
from pathlib import Path
import warnings

try:
    import dask.array as da
    DASK_AVAILABLE = True
except ImportError:
    DASK_AVAILABLE = False
    warnings.warn("Dask not available. Large array support will be limited.")


@dataclass
class ChannelQualityResult:
    """Results from channel quality assessment.

    Attributes:
        channel_name: Identifier for the channel
        confidence_score: Overall quality confidence (0-1)
        quality_metrics: Dictionary of computed quality metrics
        requires_review: Whether manual review is needed
        problem_regions: List of (x1, y1, x2, y2) bounding boxes for problem areas
        metadata: Additional metadata about the assessment
    """
    channel_name: str
    confidence_score: float
    quality_metrics: Dict[str, float]
    requires_review: bool
    problem_regions: List[Tuple[int, int, int, int]]
    metadata: Optional[Dict] = None

    def to_dict(self) -> Dict:
        """Convert to dictionary for serialization."""
        return asdict(self)


class ChannelAssessor:
    """Assess channel quality using deep learning models.

    This class provides automated quality assessment for microscopy channels,
    processing large images efficiently using tiling strategies and GPU acceleration.

    Example:
        >>> assessor = ChannelAssessor(
        ...     model_path='models/quality_net.pth',
        ...     device='cuda',
        ...     tile_size=1024
        ... )
        >>> result = assessor.process_image(channel_data, 'CD3_cycle1')
        >>> if result.requires_review:
        ...     print(f"Manual review needed: {result.confidence_score}")
    """

    def __init__(
        self,
        model_path: Optional[str] = None,
        device: Optional[str] = None,
        tile_size: int = 512,
        overlap: int = 64,
        confidence_threshold: float = 0.85,
        batch_size: int = 4,
        use_heuristics: bool = True
    ):
        """Initialize the channel quality assessor.

        Args:
            model_path: Path to trained PyTorch model. If None, uses heuristic-only mode.
            device: Device for PyTorch ('cuda', 'cpu', or None for auto-detect)
            tile_size: Size of tiles for processing large images
            overlap: Overlap between tiles to avoid edge artifacts
            confidence_threshold: Threshold below which manual review is needed
            batch_size: Number of tiles to process in parallel
            use_heuristics: Whether to use heuristic metrics in addition to DL model
        """
        self.tile_size = tile_size
        self.overlap = overlap
        self.confidence_threshold = confidence_threshold
        self.batch_size = batch_size
        self.use_heuristics = use_heuristics

        # Setup device
        if device is None:
            self.device = 'cuda' if torch.cuda.is_available() else 'cpu'
        else:
            self.device = device

        # Load model if provided
        self.model = None
        if model_path is not None:
            self.model = self._load_model(model_path)
        elif not use_heuristics:
            raise ValueError("Either model_path must be provided or use_heuristics must be True")

    def _load_model(self, model_path: str) -> torch.nn.Module:
        """Load PyTorch model from path.

        Args:
            model_path: Path to saved model checkpoint

        Returns:
            Loaded PyTorch model in eval mode
        """
        if not Path(model_path).exists():
            raise FileNotFoundError(f"Model not found: {model_path}")

        try:
            # Try loading as full model
            model = torch.load(model_path, map_location=self.device)
        except Exception:
            # Try loading as state dict
            try:
                model = self._create_default_model()
                state_dict = torch.load(model_path, map_location=self.device)
                model.load_state_dict(state_dict)
            except Exception as e:
                raise RuntimeError(f"Failed to load model from {model_path}: {e}")

        model.eval()
        model.to(self.device)
        return model

    def _create_default_model(self) -> torch.nn.Module:
        """Create a default quality assessment model architecture.

        This is a simple CNN that can be replaced with your custom architecture.
        """
        return torch.nn.Sequential(
            torch.nn.Conv2d(1, 32, 3, padding=1),
            torch.nn.ReLU(),
            torch.nn.MaxPool2d(2),
            torch.nn.Conv2d(32, 64, 3, padding=1),
            torch.nn.ReLU(),
            torch.nn.AdaptiveAvgPool2d(1),
            torch.nn.Flatten(),
            torch.nn.Linear(64, 1),
            torch.nn.Sigmoid()
        )

    def process_image(
        self,
        image: Union[np.ndarray, da.Array],
        channel_name: str,
        metadata: Optional[Dict] = None
    ) -> ChannelQualityResult:
        """Process a full image and assess quality.

        Args:
            image: 2D numpy or dask array
            channel_name: Identifier for this channel
            metadata: Optional metadata to include in results

        Returns:
            ChannelQualityResult with assessment details
        """
        # Convert dask to numpy if needed (compute chunks lazily during tiling)
        if DASK_AVAILABLE and isinstance(image, da.Array):
            # Keep as dask array, we'll compute tiles as needed
            h, w = image.shape
            is_dask = True
        else:
            h, w = image.shape[:2]
            is_dask = False

        # Create confidence map
        confidence_map = np.zeros((h, w), dtype=np.float32)
        count_map = np.zeros((h, w), dtype=np.int32)
        quality_scores = []
        problem_regions = []

        # Calculate tiling parameters
        step_size = self.tile_size - self.overlap
        n_tiles_y = (h - self.overlap) // step_size + 1
        n_tiles_x = (w - self.overlap) // step_size + 1

        # Process in tiles
        tiles_batch = []
        tile_coords = []

        for y_idx in range(n_tiles_y):
            y = min(y_idx * step_size, h - self.tile_size)
            if y < 0:
                y = 0

            for x_idx in range(n_tiles_x):
                x = min(x_idx * step_size, w - self.tile_size)
                if x < 0:
                    x = 0

                # Extract tile
                if is_dask:
                    tile = image[y:y+self.tile_size, x:x+self.tile_size].compute()
                else:
                    tile = image[y:y+self.tile_size, x:x+self.tile_size]

                # Ensure tile is correct size (handle edges)
                if tile.shape[0] < self.tile_size or tile.shape[1] < self.tile_size:
                    tile_padded = np.zeros((self.tile_size, self.tile_size), dtype=tile.dtype)
                    tile_padded[:tile.shape[0], :tile.shape[1]] = tile
                    tile = tile_padded

                tiles_batch.append(tile)
                tile_coords.append((x, y))

                # Process batch when full
                if len(tiles_batch) >= self.batch_size:
                    self._process_tile_batch(
                        tiles_batch, tile_coords, confidence_map,
                        count_map, quality_scores, problem_regions
                    )
                    tiles_batch = []
                    tile_coords = []

        # Process remaining tiles
        if tiles_batch:
            self._process_tile_batch(
                tiles_batch, tile_coords, confidence_map,
                count_map, quality_scores, problem_regions
            )

        # Average overlapping regions
        confidence_map = np.divide(
            confidence_map, count_map,
            out=np.zeros_like(confidence_map),
            where=count_map > 0
        )

        # Calculate overall metrics
        avg_confidence = np.mean(quality_scores) if quality_scores else 0.0
        quality_metrics = self._calculate_quality_metrics(image, confidence_map)

        # Merge nearby problem regions
        problem_regions = self._merge_nearby_regions(problem_regions)

        return ChannelQualityResult(
            channel_name=channel_name,
            confidence_score=float(avg_confidence),
            quality_metrics=quality_metrics,
            requires_review=avg_confidence < self.confidence_threshold,
            problem_regions=problem_regions,
            metadata=metadata
        )

    def _process_tile_batch(
        self,
        tiles: List[np.ndarray],
        coords: List[Tuple[int, int]],
        confidence_map: np.ndarray,
        count_map: np.ndarray,
        quality_scores: List[float],
        problem_regions: List[Tuple[int, int, int, int]]
    ):
        """Process a batch of tiles."""
        for tile, (x, y) in zip(tiles, coords):
            tile_quality = self._assess_tile(tile)
            quality_scores.append(tile_quality['score'])

            # Update confidence map
            tile_conf = tile_quality['confidence']
            h, w = tile_conf.shape
            confidence_map[y:y+h, x:x+w] += tile_conf
            count_map[y:y+h, x:x+w] += 1

            # Track problematic regions
            if tile_quality['score'] < self.confidence_threshold:
                problem_regions.append((x, y, x+self.tile_size, y+self.tile_size))

    def _assess_tile(self, tile: np.ndarray) -> Dict:
        """Assess quality of a single tile.

        Args:
            tile: 2D numpy array

        Returns:
            Dictionary with 'score' and 'confidence' map
        """
        # Normalize tile
        tile_normalized = self._normalize_tile(tile)

        # DL model assessment
        if self.model is not None:
            with torch.no_grad():
                tile_tensor = torch.from_numpy(tile_normalized).float().to(self.device)
                if len(tile_tensor.shape) == 2:
                    tile_tensor = tile_tensor.unsqueeze(0).unsqueeze(0)

                # Get model prediction
                output = self.model(tile_tensor)

                # Handle different output formats
                if isinstance(output, dict):
                    dl_score = output.get('quality_score', output.get('score', 0.5)).item()
                    conf_map = output.get('confidence_map', None)
                elif isinstance(output, tuple):
                    dl_score = output[0].item()
                    conf_map = output[1] if len(output) > 1 else None
                else:
                    dl_score = output.item() if output.numel() == 1 else output.mean().item()
                    conf_map = None

                if conf_map is not None:
                    if isinstance(conf_map, torch.Tensor):
                        conf_map = conf_map.cpu().numpy().squeeze()
                else:
                    conf_map = np.full(tile.shape, dl_score, dtype=np.float32)
        else:
            dl_score = 0.5
            conf_map = np.full(tile.shape, 0.5, dtype=np.float32)

        # Heuristic assessment
        if self.use_heuristics:
            heuristic_score = self._assess_tile_heuristics(tile)
            # Combine DL and heuristic scores
            final_score = 0.7 * dl_score + 0.3 * heuristic_score
        else:
            final_score = dl_score

        return {
            'score': final_score,
            'confidence': conf_map,
            'features': {
                'dl_score': dl_score,
                'heuristic_score': heuristic_score if self.use_heuristics else None
            }
        }

    def _normalize_tile(self, tile: np.ndarray) -> np.ndarray:
        """Normalize tile intensities to [0, 1]."""
        tile = tile.astype(np.float32)
        tile_min = tile.min()
        tile_max = tile.max()
        if tile_max > tile_min:
            return (tile - tile_min) / (tile_max - tile_min)
        return tile

    def _assess_tile_heuristics(self, tile: np.ndarray) -> float:
        """Assess tile quality using heuristic metrics.

        Args:
            tile: 2D numpy array

        Returns:
            Quality score between 0 and 1
        """
        scores = []

        # SNR score
        signal = np.percentile(tile, 90)
        noise = np.std(tile[tile < np.percentile(tile, 20)])
        if noise > 0:
            snr = signal / noise
            snr_score = min(1.0, snr / 10.0)  # Normalize to ~10 as good SNR
            scores.append(snr_score)

        # Contrast score
        contrast = tile.max() - tile.min()
        max_possible = 65535 if tile.dtype == np.uint16 else 255
        contrast_score = min(1.0, contrast / (0.3 * max_possible))
        scores.append(contrast_score)

        # Non-uniformity score (penalize flat regions)
        gradient_magnitude = np.abs(np.gradient(tile.astype(float))).mean()
        uniformity_score = min(1.0, gradient_magnitude / 100.0)
        scores.append(uniformity_score)

        # Saturation penalty
        if tile.dtype == np.uint16:
            saturated_pixels = np.sum((tile == 0) | (tile >= 65535))
        else:
            saturated_pixels = np.sum((tile == 0) | (tile >= 255))
        saturation_ratio = saturated_pixels / tile.size
        saturation_score = 1.0 - min(1.0, saturation_ratio * 10)
        scores.append(saturation_score)

        return np.mean(scores)

    def _calculate_quality_metrics(
        self,
        image: Union[np.ndarray, da.Array],
        confidence_map: np.ndarray
    ) -> Dict[str, float]:
        """Calculate comprehensive quality metrics for the full image.

        Args:
            image: Full image array
            confidence_map: Confidence scores per pixel

        Returns:
            Dictionary of quality metrics
        """
        # Compute image if it's a dask array
        if DASK_AVAILABLE and isinstance(image, da.Array):
            # Sample for efficiency
            sample_size = min(1000, image.shape[0])
            sample_indices = np.linspace(0, image.shape[0]-1, sample_size, dtype=int)
            image_sample = image[sample_indices].compute()
        else:
            image_sample = image

        metrics = {}

        # Basic intensity metrics
        metrics['mean_intensity'] = float(np.mean(image_sample))
        metrics['std_intensity'] = float(np.std(image_sample))
        metrics['min_intensity'] = float(np.min(image_sample))
        metrics['max_intensity'] = float(np.max(image_sample))

        # SNR
        signal = np.percentile(image_sample, 90)
        background = np.percentile(image_sample, 10)
        noise = np.std(image_sample[image_sample < np.percentile(image_sample, 20)])
        if noise > 0:
            metrics['snr'] = float(signal / noise)
        else:
            metrics['snr'] = 0.0

        # Confidence metrics
        metrics['mean_confidence'] = float(np.mean(confidence_map))
        metrics['std_confidence'] = float(np.std(confidence_map))
        metrics['low_confidence_ratio'] = float(
            np.sum(confidence_map < self.confidence_threshold) / confidence_map.size
        )

        # Autofluorescence estimate
        if signal > 0:
            metrics['autofluorescence_score'] = float(background / signal)
        else:
            metrics['autofluorescence_score'] = 1.0

        # Dynamic range
        metrics['dynamic_range'] = float(signal - background)

        return metrics

    def _merge_nearby_regions(
        self,
        regions: List[Tuple[int, int, int, int]],
        merge_distance: int = 50
    ) -> List[Tuple[int, int, int, int]]:
        """Merge nearby problem regions to reduce fragmentation.

        Args:
            regions: List of (x1, y1, x2, y2) bounding boxes
            merge_distance: Maximum distance to merge regions

        Returns:
            Merged list of regions
        """
        if not regions:
            return []

        # Simple greedy merging
        merged = [list(regions[0])]

        for x1, y1, x2, y2 in regions[1:]:
            merged_with_existing = False

            for i, (mx1, my1, mx2, my2) in enumerate(merged):
                # Check if regions are close
                if (abs(x1 - mx2) < merge_distance or abs(mx1 - x2) < merge_distance) and \
                   (abs(y1 - my2) < merge_distance or abs(my1 - y2) < merge_distance):
                    # Merge by expanding bounding box
                    merged[i] = [
                        min(mx1, x1),
                        min(my1, y1),
                        max(mx2, x2),
                        max(my2, y2)
                    ]
                    merged_with_existing = True
                    break

            if not merged_with_existing:
                merged.append([x1, y1, x2, y2])

        return [tuple(r) for r in merged]


class HeuristicChannelAssessor(ChannelAssessor):
    """Channel assessor that only uses heuristic metrics (no DL model required).

    This is useful for quick assessment without trained models or as a baseline.

    Example:
        >>> assessor = HeuristicChannelAssessor(confidence_threshold=0.7)
        >>> result = assessor.process_image(channel_data, 'DAPI')
    """

    def __init__(
        self,
        tile_size: int = 512,
        overlap: int = 64,
        confidence_threshold: float = 0.75,
        **kwargs
    ):
        """Initialize heuristic-only assessor."""
        super().__init__(
            model_path=None,
            tile_size=tile_size,
            overlap=overlap,
            confidence_threshold=confidence_threshold,
            use_heuristics=True,
            **kwargs
        )
