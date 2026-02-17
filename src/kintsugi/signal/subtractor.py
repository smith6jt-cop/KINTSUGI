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
    analyze_for_weighted_subtraction,
    compute_subtraction_quality,
    compute_weighted_subtraction_quality,
    subtract_autofluorescence,
    subtract_autofluorescence_dask,
    subtract_autofluorescence_weighted,
)

logger = logging.getLogger("kintsugi.signal.subtractor")


@dataclass
class IntensityRange:
    """A single intensity range for weighted subtraction."""

    lower_bound: float = 0.0
    upper_bound: float = 65535.0
    weight: float = 1.0
    label: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> IntensityRange:
        return cls(**{k: v for k, v in d.items() if k in cls.__dataclass_fields__})


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
class WeightedSubtractionParameters:
    """Parameters for weighted autofluorescence subtraction."""

    blank_clip_factor: int = 0
    base_scale_factor: float = 1.0
    n_ranges: int = 5
    range_method: str = "percentile"
    ranges: list[IntensityRange] = field(default_factory=list)
    transition_width: float = 0.1
    smooth_low: bool = False
    low_size: int = 2
    low_percentile: int = 60
    smooth_high: bool = False
    high_size: int = 2
    high_percentile: int = 90
    erosion: int = 0

    def to_dict(self) -> dict:
        """Convert to dictionary with nested ranges."""
        d = asdict(self)
        d["ranges"] = [r.to_dict() if isinstance(r, IntensityRange) else r for r in self.ranges]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> WeightedSubtractionParameters:
        """Create from dictionary with nested range deserialization."""
        d = dict(d)
        if "ranges" in d:
            d["ranges"] = [
                IntensityRange.from_dict(r) if isinstance(r, dict) else r for r in d["ranges"]
            ]
        valid_keys = cls.__dataclass_fields__
        return cls(**{k: v for k, v in d.items() if k in valid_keys})

    def to_subtraction_kwargs(self) -> dict:
        """Convert to kwargs for subtract_autofluorescence_weighted()."""
        ranges_dicts = [r.to_dict() if isinstance(r, IntensityRange) else r for r in self.ranges]
        return {
            "blank_clip_factor": self.blank_clip_factor,
            "base_scale_factor": self.base_scale_factor,
            "n_ranges": self.n_ranges,
            "range_method": self.range_method,
            "ranges": ranges_dicts if ranges_dicts else None,
            "transition_width": self.transition_width,
            "smooth_low": self.smooth_low,
            "low_size": self.low_size,
            "low_percentile": self.low_percentile,
            "smooth_high": self.smooth_high,
            "high_size": self.high_size,
            "high_percentile": self.high_percentile,
            "erosion": self.erosion,
        }


@dataclass
class SubtractionResult:
    """Result of autofluorescence subtraction."""

    image: np.ndarray
    parameters: SubtractionParameters | WeightedSubtractionParameters
    quality_metrics: dict = field(default_factory=dict)
    analysis: dict = field(default_factory=dict)
    range_metrics: list[dict] | None = None
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    def quality_passed(self, threshold: float = 0.5) -> bool:
        """Check if quality score passes threshold."""
        # Support both flat and nested quality_metrics
        score = self.quality_metrics.get("quality_score", 0)
        if score == 0 and "global" in self.quality_metrics:
            score = self.quality_metrics["global"].get("quality_score", 0)
        return score >= threshold


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
        method: str = "global",
    ):
        self.project_dir = Path(project_dir) if project_dir else None
        self.tissue_type = tissue_type
        self.auto_learn = auto_learn
        self.quality_threshold = quality_threshold
        self.method = method  # "global" or "weighted"

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
                    operation="blank_subtraction",
                    tissue_type=tissue or "unknown",
                    marker_name=marker,
                )
                if recommendation.get("found") and recommendation.get("confidence", 0) > 0.5:
                    learned_params = recommendation.get("recommended_parameters", {})
                    logger.info(
                        f"Using learned parameters for {marker} "
                        f"(confidence: {recommendation['confidence']:.2f})"
                    )
            except Exception as e:
                logger.warning(f"Failed to get learned parameters: {e}")

        # Combine analysis and learned parameters
        if learned_params:
            # Weight learned params more heavily if confidence is high
            params = self._merge_parameters(analysis, learned_params)
        else:
            params = SubtractionParameters(
                blank_clip_factor=analysis["blank_clip_factor"],
                blank_scale_factor=analysis["blank_scale_factor"],
                smooth_low=analysis["smooth_low"],
                low_size=analysis["low_size"],
                low_percentile=analysis["low_percentile"],
                smooth_high=analysis["smooth_high"],
                high_size=analysis["high_size"],
                high_percentile=analysis["high_percentile"],
                erosion=analysis["erosion"],
            )

        return params

    def suggest_weighted_parameters(
        self,
        signal: np.ndarray,
        blank: np.ndarray,
        marker: str | None = None,
        tissue_type: str | None = None,
        n_ranges: int = 5,
        range_method: str = "percentile",
    ) -> WeightedSubtractionParameters:
        """
        Suggest optimal parameters for weighted subtraction.

        Parameters
        ----------
        signal, blank : np.ndarray
            Signal and blank images
        marker : str, optional
            Marker name
        tissue_type : str, optional
            Tissue type
        n_ranges : int
            Number of intensity ranges
        range_method : str
            Method for range boundaries

        Returns
        -------
        WeightedSubtractionParameters
        """
        tissue = tissue_type or self.tissue_type
        analysis = analyze_for_weighted_subtraction(
            signal, blank, tissue, marker, n_ranges=n_ranges, range_method=range_method
        )

        intensity_ranges = [
            IntensityRange(
                lower_bound=r["lower_bound"],
                upper_bound=r["upper_bound"],
                weight=r["weight"],
                label=r.get("label", ""),
            )
            for r in analysis.get("ranges", [])
        ]

        return WeightedSubtractionParameters(
            blank_clip_factor=analysis["blank_clip_factor"],
            base_scale_factor=analysis["base_scale_factor"],
            n_ranges=analysis.get("n_ranges", n_ranges),
            range_method=analysis.get("range_method", range_method),
            ranges=intensity_ranges,
            transition_width=analysis.get("transition_width", 0.1),
            smooth_low=analysis["smooth_low"],
            low_size=analysis["low_size"],
            low_percentile=analysis["low_percentile"],
            smooth_high=analysis["smooth_high"],
            high_size=analysis["high_size"],
            high_percentile=analysis["high_percentile"],
            erosion=analysis["erosion"],
        )

    def process(
        self,
        signal: np.ndarray,
        blank: np.ndarray,
        params: SubtractionParameters | WeightedSubtractionParameters | None = None,
        marker: str | None = None,
        tissue_type: str | None = None,
        use_dask: bool = False,
        method: str | None = None,
    ) -> SubtractionResult:
        """
        Process signal channel to remove autofluorescence.

        Parameters
        ----------
        signal : np.ndarray or dask.array.Array
            Signal channel image
        blank : np.ndarray or dask.array.Array
            Blank channel image
        params : SubtractionParameters or WeightedSubtractionParameters, optional
            Subtraction parameters (auto-suggested if not provided)
        marker : str, optional
            Marker name (for learning)
        tissue_type : str, optional
            Tissue type (for learning)
        use_dask : bool
            Use dask-based processing for large images
        method : str, optional
            Override instance method ("global" or "weighted")

        Returns
        -------
        SubtractionResult
            Result containing processed image, parameters, and quality metrics
        """
        tissue = tissue_type or self.tissue_type
        active_method = method or self.method

        # Convert dask to numpy for analysis if needed
        signal_np = signal.compute() if hasattr(signal, "compute") else signal
        blank_np = blank.compute() if hasattr(blank, "compute") else blank

        # Dispatch based on method
        if active_method == "weighted":
            return self._process_weighted(
                signal, blank, signal_np, blank_np, params, marker, tissue, use_dask
            )
        else:
            return self._process_global(
                signal, blank, signal_np, blank_np, params, marker, tissue, use_dask
            )

    def _process_global(
        self, signal, blank, signal_np, blank_np, params, marker, tissue, use_dask
    ) -> SubtractionResult:
        """Process using global (single scale factor) method."""
        if params is None:
            params = self.suggest_parameters(signal_np, blank_np, marker, tissue)

        if use_dask or hasattr(signal, "compute"):
            result_image = subtract_autofluorescence_dask(signal, blank, **params.to_dict())
            result_image_np = (
                result_image.compute() if hasattr(result_image, "compute") else result_image
            )
        else:
            result_image = subtract_autofluorescence(signal, blank, **params.to_dict())
            result_image_np = result_image

        quality = compute_subtraction_quality(signal_np, result_image_np, blank_np)

        result = SubtractionResult(
            image=result_image,
            parameters=params,
            quality_metrics=quality,
            analysis={"marker": marker, "tissue_type": tissue, "method": "global"},
        )

        self._auto_learn(result, marker, tissue, quality)
        return result

    def _process_weighted(
        self, signal, blank, signal_np, blank_np, params, marker, tissue, use_dask
    ) -> SubtractionResult:
        """Process using per-intensity-range weighted method."""
        if params is None:
            params = self.suggest_weighted_parameters(signal_np, blank_np, marker, tissue)

        if isinstance(params, WeightedSubtractionParameters):
            kwargs = params.to_subtraction_kwargs()
        else:
            # Fallback: convert SubtractionParameters to weighted kwargs
            kwargs = {
                "blank_clip_factor": params.blank_clip_factor,
                "base_scale_factor": getattr(params, "blank_scale_factor", 1.0),
            }

        result_image = subtract_autofluorescence_weighted(signal_np, blank_np, **kwargs)
        result_image_np = result_image

        # Compute weighted quality metrics if we have ranges
        ranges = kwargs.get("ranges")
        if ranges:
            from .autofluorescence import compute_intensity_ranges

            quality = compute_weighted_subtraction_quality(
                signal_np, result_image_np, blank_np, ranges
            )
            range_metrics = quality.get("per_range")
            quality_flat = quality.get("global", quality)
        else:
            # Compute ranges for quality assessment
            signal_f = signal_np.astype(np.float32)
            blank_f = blank_np.astype(np.float32)
            clip = kwargs.get("blank_clip_factor", 0)
            blank_clipped = np.clip(blank_f, clip, blank_f.max())
            blank_clipped[blank_clipped <= clip] = 0

            from .autofluorescence import compute_intensity_ranges

            computed_ranges = compute_intensity_ranges(signal_f, blank_clipped)
            quality = compute_weighted_subtraction_quality(
                signal_np, result_image_np, blank_np, computed_ranges
            )
            range_metrics = quality.get("per_range")
            quality_flat = quality.get("global", quality)

        result = SubtractionResult(
            image=result_image,
            parameters=params,
            quality_metrics=quality_flat,
            analysis={"marker": marker, "tissue_type": tissue, "method": "weighted"},
            range_metrics=range_metrics,
        )

        self._auto_learn(result, marker, tissue, quality_flat, algorithm_version="weighted_v1")
        return result

    def _auto_learn(
        self,
        result: SubtractionResult,
        marker: str | None,
        tissue: str | None,
        quality: dict,
        algorithm_version: str = "v1",
    ):
        """Record parameters if quality passed."""
        quality_score = quality.get("quality_score", 0)
        if (
            self.auto_learn
            and result.quality_passed(self.quality_threshold)
            and self.learning_engine
            and marker
        ):
            try:
                params_dict = result.parameters.to_dict()
                params_dict["algorithm_version"] = algorithm_version
                self.learning_engine.record_parameters(
                    operation="blank_subtraction",
                    tissue_type=tissue or "unknown",
                    marker_name=marker,
                    parameters=params_dict,
                    quality_score=quality_score,
                )
                logger.info(
                    f"Recorded {algorithm_version} parameters for {marker} "
                    f"(quality: {quality_score:.3f})"
                )
            except Exception as e:
                logger.warning(f"Failed to record parameters: {e}")

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
                    signal,
                    blank,
                    marker=marker,
                    tissue_type=tissue_type,
                )
                results[marker] = result
                logger.info(
                    f"Processed {marker}: quality={result.quality_metrics['quality_score']:.3f}"
                )
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
            analysis["blank_clip_factor"] * analysis_weight
            + learned.get("blank_clip_factor", analysis["blank_clip_factor"])
            * (1 - analysis_weight)
        )
        scale = round(
            analysis["blank_scale_factor"] * analysis_weight
            + learned.get("blank_scale_factor", analysis["blank_scale_factor"])
            * (1 - analysis_weight),
            2,
        )

        # For boolean parameters, use learned if available
        smooth_low = learned.get("smooth_low", analysis["smooth_low"])
        smooth_high = learned.get("smooth_high", analysis["smooth_high"])

        # Other parameters from learned if available, else analysis
        return SubtractionParameters(
            blank_clip_factor=clip,
            blank_scale_factor=scale,
            smooth_low=smooth_low,
            low_size=learned.get("low_size", analysis["low_size"]),
            low_percentile=learned.get("low_percentile", analysis["low_percentile"]),
            smooth_high=smooth_high,
            high_size=learned.get("high_size", analysis["high_size"]),
            high_percentile=learned.get("high_percentile", analysis["high_percentile"]),
            erosion=learned.get("erosion", analysis["erosion"]),
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
        summary = {"timestamp": datetime.now().isoformat(), "channels": {}}

        for marker, result in results.items():
            summary["channels"][marker] = {
                "parameters": result.parameters.to_dict(),
                "quality_metrics": result.quality_metrics,
                "quality_passed": result.quality_passed(self.quality_threshold),
            }

            if save_images:
                img = result.image
                if hasattr(img, "compute"):
                    img = img.compute()
                tifffile.imwrite(str(output_dir / f"{marker}.tif"), img.astype(np.uint16))

        # Save summary JSON
        with open(output_dir / "subtraction_summary.json", "w") as f:
            json.dump(summary, f, indent=2)

        # Save parameters for each channel
        params_dir = output_dir / "parameters"
        params_dir.mkdir(exist_ok=True)

        for marker, result in results.items():
            with open(params_dir / f"{marker}_params.json", "w") as f:
                json.dump(
                    {
                        "marker": marker,
                        "parameters": result.parameters.to_dict(),
                        "quality_metrics": result.quality_metrics,
                        "timestamp": result.timestamp,
                    },
                    f,
                    indent=2,
                )

        logger.info(f"Saved {len(results)} results to {output_dir}")
