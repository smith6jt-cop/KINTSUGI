"""
AutofluorescenceSubtractor Class

High-level class-based interface for autofluorescence subtraction with
integrated parameter learning and quality assessment.
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np

from .autofluorescence import (
    analyze_for_subtraction,
    compute_subtraction_quality,
    subtract_autofluorescence,
    subtract_autofluorescence_dask,
)

logger = logging.getLogger("kintsugi.signal.subtractor")


@dataclass
class SubtractionParameters:
    """Parameters for autofluorescence subtraction."""

    blank_clip_factor: int = 0
    blank_scale_factor: float = 1.0
    smooth_low: bool = False
    low_size: int = 2
    low_percentile: int = 60
    smooth_high: bool = False
    high_size: int = 2
    high_percentile: int = 90
    erosion: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> SubtractionParameters:
        """Create from dictionary."""
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


@dataclass
class SubtractionResult:
    """Result of autofluorescence subtraction."""

    image: np.ndarray
    parameters: SubtractionParameters
    quality_metrics: dict = field(default_factory=dict)
    analysis: dict = field(default_factory=dict)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def quality_passed(self, threshold: float = 0.5) -> bool:
        """Check if quality score passes threshold."""
        return self.quality_metrics.get('quality_score', 0) >= threshold


class AutofluorescenceSubtractor:
    """
    High-level interface for autofluorescence subtraction with learning.

    This class provides:
    - Intelligent parameter suggestion based on image analysis
    - Integration with parameter learning database
    - Quality assessment of results
    - Batch processing capabilities

    Parameters
    ----------
    project_dir : str or Path, optional
        Project directory for learning database
    tissue_type : str, optional
        Default tissue type for recommendations
    auto_learn : bool
        Automatically record successful parameters (default: True)
    quality_threshold : float
        Minimum quality score to consider successful (default: 0.5)

    Examples
    --------
    >>> # Basic usage
    >>> subtractor = AutofluorescenceSubtractor(project_dir='./my_project')
    >>> result = subtractor.process(signal_image, blank_image, marker='CD3')
    >>> print(f"Quality: {result.quality_metrics['quality_score']:.3f}")

    >>> # With parameter suggestions
    >>> params = subtractor.suggest_parameters(signal_image, blank_image, marker='CD3')
    >>> result = subtractor.process(signal_image, blank_image, params=params)

    >>> # Batch processing
    >>> results = subtractor.process_batch(channels, blank_channels)
    """

    def __init__(
        self,
        project_dir: str | Path | None = None,
        tissue_type: str | None = None,
        auto_learn: bool = True,
        quality_threshold: float = 0.5,
    ):
        self.project_dir = Path(project_dir) if project_dir else None
        self.tissue_type = tissue_type
        self.auto_learn = auto_learn
        self.quality_threshold = quality_threshold

        # Learning engine (lazy initialization)
        self._learning_engine = None

    @property
    def learning_engine(self):
        """Lazy initialization of learning engine."""
        if self._learning_engine is None and self.project_dir:
            try:
                from kintsugi.mcp.tools.learning import ParameterLearningEngine
                self._learning_engine = ParameterLearningEngine(str(self.project_dir))
            except ImportError:
                logger.warning("Parameter learning not available. Install kintsugi[claude]")
        return self._learning_engine

    def suggest_parameters(
        self,
        signal: np.ndarray,
        blank: np.ndarray,
        marker: str | None = None,
        tissue_type: str | None = None,
        use_learning: bool = True,
    ) -> SubtractionParameters:
        """
        Suggest optimal parameters for subtraction.

        Combines image analysis with learned parameters from previous
        successful runs.

        Parameters
        ----------
        signal : np.ndarray
            Signal channel image
        blank : np.ndarray
            Blank channel image
        marker : str, optional
            Marker name for learned recommendations
        tissue_type : str, optional
            Tissue type (uses default if not specified)
        use_learning : bool
            Use learned parameters if available (default: True)

        Returns
        -------
        SubtractionParameters
            Suggested parameters
        """
        tissue = tissue_type or self.tissue_type

        # Get analysis-based suggestions
        analysis = analyze_for_subtraction(signal, blank, tissue, marker)

        # Try to get learned parameters
        learned_params = None
        if use_learning and self.learning_engine and marker:
            try:
                recommendation = self.learning_engine.recommend_parameters(
                    operation='blank_subtraction',
                    tissue_type=tissue or 'unknown',
                    marker_name=marker,
                )
                if recommendation.get('found') and recommendation.get('confidence', 0) > 0.5:
                    learned_params = recommendation.get('recommended_parameters', {})
                    logger.info(f"Using learned parameters for {marker} "
                              f"(confidence: {recommendation['confidence']:.2f})")
            except Exception as e:
                logger.warning(f"Failed to get learned parameters: {e}")

        # Combine analysis and learned parameters
        if learned_params:
            # Weight learned params more heavily if confidence is high
            params = self._merge_parameters(analysis, learned_params)
        else:
            params = SubtractionParameters(
                blank_clip_factor=analysis['blank_clip_factor'],
                blank_scale_factor=analysis['blank_scale_factor'],
                smooth_low=analysis['smooth_low'],
                low_size=analysis['low_size'],
                low_percentile=analysis['low_percentile'],
                smooth_high=analysis['smooth_high'],
                high_size=analysis['high_size'],
                high_percentile=analysis['high_percentile'],
                erosion=analysis['erosion'],
            )

        return params

    def process(
        self,
        signal: np.ndarray,
        blank: np.ndarray,
        params: SubtractionParameters | None = None,
        marker: str | None = None,
        tissue_type: str | None = None,
        use_dask: bool = False,
    ) -> SubtractionResult:
        """
        Process signal channel to remove autofluorescence.

        Parameters
        ----------
        signal : np.ndarray or dask.array.Array
            Signal channel image
        blank : np.ndarray or dask.array.Array
            Blank channel image
        params : SubtractionParameters, optional
            Subtraction parameters (auto-suggested if not provided)
        marker : str, optional
            Marker name (for learning)
        tissue_type : str, optional
            Tissue type (for learning)
        use_dask : bool
            Use dask-based processing for large images

        Returns
        -------
        SubtractionResult
            Result containing processed image, parameters, and quality metrics
        """
        tissue = tissue_type or self.tissue_type

        # Auto-suggest parameters if not provided
        if params is None:
            # Convert dask to numpy for analysis
            signal_np = signal.compute() if hasattr(signal, 'compute') else signal
            blank_np = blank.compute() if hasattr(blank, 'compute') else blank
            params = self.suggest_parameters(signal_np, blank_np, marker, tissue)

        # Perform subtraction
        if use_dask or hasattr(signal, 'compute'):
            result_image = subtract_autofluorescence_dask(
                signal, blank, **params.to_dict()
            )
            # Compute for quality assessment
            if hasattr(result_image, 'compute'):
                result_image_np = result_image.compute()
            else:
                result_image_np = result_image
        else:
            result_image = subtract_autofluorescence(
                signal, blank, **params.to_dict()
            )
            result_image_np = result_image

        # Compute quality metrics
        signal_np = signal.compute() if hasattr(signal, 'compute') else signal
        blank_np = blank.compute() if hasattr(blank, 'compute') else blank

        quality = compute_subtraction_quality(signal_np, result_image_np, blank_np)

        # Create result
        result = SubtractionResult(
            image=result_image,
            parameters=params,
            quality_metrics=quality,
            analysis={
                'marker': marker,
                'tissue_type': tissue,
            }
        )

        # Auto-learn if quality passed
        if (self.auto_learn and result.quality_passed(self.quality_threshold)
                and self.learning_engine and marker):
            try:
                self.learning_engine.record_parameters(
                    operation='blank_subtraction',
                    tissue_type=tissue or 'unknown',
                    marker_name=marker,
                    parameters=params.to_dict(),
                    quality_score=quality['quality_score'],
                )
                logger.info(f"Recorded parameters for {marker} (quality: {quality['quality_score']:.3f})")
            except Exception as e:
                logger.warning(f"Failed to record parameters: {e}")

        return result

    def process_batch(
        self,
        channels: dict[str, np.ndarray],
        blank_channels: dict[str, np.ndarray],
        channel_blank_mapping: dict[str, str] | None = None,
        tissue_type: str | None = None,
        progress_callback: callable | None = None,
    ) -> dict[str, SubtractionResult]:
        """
        Process multiple channels in batch.

        Parameters
        ----------
        channels : dict
            Dictionary of {marker_name: signal_image}
        blank_channels : dict
            Dictionary of {blank_name: blank_image}
        channel_blank_mapping : dict, optional
            Mapping of marker names to blank names.
            If not provided, best blank is auto-selected.
        tissue_type : str, optional
            Tissue type for all channels
        progress_callback : callable, optional
            Function called with (current, total) for progress

        Returns
        -------
        dict
            Dictionary of {marker_name: SubtractionResult}
        """
        from .autofluorescence import suggest_blank_channel

        results = {}
        total = len(channels)

        for i, (marker, signal) in enumerate(channels.items()):
            if progress_callback:
                progress_callback(i, total)

            # Determine blank to use
            if channel_blank_mapping and marker in channel_blank_mapping:
                blank_name = channel_blank_mapping[marker]
                blank = blank_channels.get(blank_name)
            else:
                # Auto-select best blank
                blank_name, _ = suggest_blank_channel(signal, blank_channels)
                blank = blank_channels.get(blank_name)

            if blank is None:
                logger.warning(f"No blank found for {marker}, skipping")
                continue

            try:
                result = self.process(
                    signal, blank,
                    marker=marker,
                    tissue_type=tissue_type,
                )
                results[marker] = result
                logger.info(f"Processed {marker}: quality={result.quality_metrics['quality_score']:.3f}")
            except Exception as e:
                logger.error(f"Failed to process {marker}: {e}")

        if progress_callback:
            progress_callback(total, total)

        return results

    def _merge_parameters(
        self,
        analysis: dict,
        learned: dict,
        analysis_weight: float = 0.4,
    ) -> SubtractionParameters:
        """Merge analysis-based and learned parameters."""
        # For numeric parameters, weighted average
        clip = int(
            analysis['blank_clip_factor'] * analysis_weight +
            learned.get('blank_clip_factor', analysis['blank_clip_factor']) * (1 - analysis_weight)
        )
        scale = round(
            analysis['blank_scale_factor'] * analysis_weight +
            learned.get('blank_scale_factor', analysis['blank_scale_factor']) * (1 - analysis_weight),
            2
        )

        # For boolean parameters, use learned if available
        smooth_low = learned.get('smooth_low', analysis['smooth_low'])
        smooth_high = learned.get('smooth_high', analysis['smooth_high'])

        # Other parameters from learned if available, else analysis
        return SubtractionParameters(
            blank_clip_factor=clip,
            blank_scale_factor=scale,
            smooth_low=smooth_low,
            low_size=learned.get('low_size', analysis['low_size']),
            low_percentile=learned.get('low_percentile', analysis['low_percentile']),
            smooth_high=smooth_high,
            high_size=learned.get('high_size', analysis['high_size']),
            high_percentile=learned.get('high_percentile', analysis['high_percentile']),
            erosion=learned.get('erosion', analysis['erosion']),
        )

    def save_results(
        self,
        results: dict[str, SubtractionResult],
        output_dir: str | Path,
        save_images: bool = True,
    ) -> None:
        """
        Save batch results to disk.

        Parameters
        ----------
        results : dict
            Results from process_batch
        output_dir : str or Path
            Output directory
        save_images : bool
            Whether to save processed images as TIFFs
        """
        import tifffile

        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Save summary
        summary = {
            'timestamp': datetime.now().isoformat(),
            'channels': {}
        }

        for marker, result in results.items():
            summary['channels'][marker] = {
                'parameters': result.parameters.to_dict(),
                'quality_metrics': result.quality_metrics,
                'quality_passed': result.quality_passed(self.quality_threshold),
            }

            if save_images:
                img = result.image
                if hasattr(img, 'compute'):
                    img = img.compute()
                tifffile.imwrite(
                    str(output_dir / f"{marker}.tif"),
                    img.astype(np.uint16)
                )

        # Save summary JSON
        with open(output_dir / 'subtraction_summary.json', 'w') as f:
            json.dump(summary, f, indent=2)

        # Save parameters for each channel
        params_dir = output_dir / 'parameters'
        params_dir.mkdir(exist_ok=True)

        for marker, result in results.items():
            with open(params_dir / f"{marker}_params.json", 'w') as f:
                json.dump({
                    'marker': marker,
                    'parameters': result.parameters.to_dict(),
                    'quality_metrics': result.quality_metrics,
                    'timestamp': result.timestamp,
                }, f, indent=2)

        logger.info(f"Saved {len(results)} results to {output_dir}")
