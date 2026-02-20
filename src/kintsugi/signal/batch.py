"""
Batch Signal Isolation Processing

Automated per-marker parameter optimization and batch autofluorescence subtraction
for CODEX multiplex immunofluorescence datasets.

Key functions:
- discover_channels: Parse config.yaml to find signal/blank pairs
- select_method: Choose global vs weighted subtraction per marker
- process_channel: Process a single channel with optimized parameters
- process_batch: Process all channels in a project
- smooth_blank_for_subtraction: Remove tile-grid intensity artifacts from blank
"""

from __future__ import annotations

import json
import logging
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path

import numpy as np
import tifffile

from .autofluorescence import (
    analyze_for_subtraction,
    analyze_for_weighted_subtraction,
    compute_subtraction_quality,
    compute_weighted_subtraction_quality,
    subtract_autofluorescence,
    subtract_autofluorescence_weighted,
)

logger = logging.getLogger("kintsugi.signal.batch")

# Channel names that are not signal markers
_SKIP_PREFIXES = ("DAPI", "Blank", "Empty")


@dataclass
class ChannelSpec:
    """Specification for a signal channel to process."""

    marker_name: str
    cycle: int
    channel_position: int  # 1-based index within cycle (1=CH1, 2=CH2, 3=CH3, 4=CH4)
    signal_path: Path
    blank_name: str
    blank_path: Path

    def to_dict(self) -> dict:
        d = asdict(self)
        d["signal_path"] = str(self.signal_path)
        d["blank_path"] = str(self.blank_path)
        return d


@dataclass
class ChannelResult:
    """Result of processing a single channel."""

    marker_name: str
    cycle: int
    blank_used: str
    method: str  # "global" or "weighted"
    parameters: dict = field(default_factory=dict)
    quality_metrics: dict = field(default_factory=dict)
    analysis: dict = field(default_factory=dict)
    zero_percent: float = 0.0
    status: str = "success"  # "success", "skipped", "error"
    warning: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class BatchResult:
    """Result of processing all channels in a project."""

    timestamp: str = ""
    method_override: str = "auto"
    tissue_type: str = ""
    tile_smooth_sigma: float = 500.0
    summary: dict = field(default_factory=dict)
    channels: dict[str, ChannelResult] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = {
            "timestamp": self.timestamp,
            "method_override": self.method_override,
            "tissue_type": self.tissue_type,
            "tile_smooth_sigma": self.tile_smooth_sigma,
            "summary": self.summary,
            "channels": {k: v.to_dict() for k, v in self.channels.items()},
        }
        return d


def smooth_blank_for_subtraction(
    blank: np.ndarray,
    sigma: float = 500.0,
    preserve_mean: bool = True,
) -> np.ndarray:
    """Smooth blank image to remove tile-grid intensity artifacts.

    Applies large-sigma Gaussian blur that averages out tile-to-tile
    intensity differences while preserving the gradual AF spatial gradient.

    Parameters
    ----------
    blank : np.ndarray
        Blank/AF image (uint16, from registered stitched output)
    sigma : float
        Gaussian sigma in pixels. Should be ~1/3 to 1/2 of tile pitch.
        For 9x7 grid on 12662x7479 image: 400-700 px recommended.
        Default: 500.0. Set to 0 to disable smoothing.
    preserve_mean : bool
        Scale smoothed image to match original mean intensity.
        Prevents systematic over/under subtraction.

    Returns
    -------
    np.ndarray
        Smoothed blank image (uint16)
    """
    if sigma <= 0:
        return blank.copy()

    from scipy.ndimage import gaussian_filter

    blank_f = blank.astype(np.float64)
    smoothed = gaussian_filter(blank_f, sigma=sigma)

    if preserve_mean:
        # Only consider non-zero pixels for mean calculation
        orig_nonzero = blank_f[blank_f > 0]
        smooth_nonzero = smoothed[smoothed > 0]
        if len(orig_nonzero) > 0 and len(smooth_nonzero) > 0:
            orig_mean = np.mean(orig_nonzero)
            smooth_mean = np.mean(smooth_nonzero)
            if smooth_mean > 0:
                smoothed *= orig_mean / smooth_mean

    return np.clip(smoothed, 0, 65535).astype(np.uint16)


def discover_channels(project_dir: str | Path) -> list[ChannelSpec]:
    """Discover signal channels and their blank pairs from project config.

    Parses workflow/config.yaml for channel_names dict and maps each signal
    marker to its corresponding cycle-1 blank based on channel position:
    - CH2 (index 1) → Blank_1a
    - CH3 (index 2) → Blank_1b
    - CH4 (index 3) → Blank_1c

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory.

    Returns
    -------
    list[ChannelSpec]
        List of channel specifications with signal/blank paths.

    Raises
    ------
    FileNotFoundError
        If config.yaml or required image files are missing.
    ValueError
        If config.yaml lacks channel_names.
    """
    import yaml

    project_dir = Path(project_dir).resolve()
    config_path = project_dir / "workflow" / "config.yaml"

    if not config_path.exists():
        raise FileNotFoundError(f"Config not found: {config_path}")

    with open(config_path) as f:
        config = yaml.safe_load(f)

    channel_names = config.get("channel_names")
    if not channel_names:
        raise ValueError(f"No channel_names in {config_path}")

    registered_dir = project_dir / "data" / "processed" / "registered"

    # Map channel position (0-indexed) to blank name
    blank_map = {1: "Blank_1a", 2: "Blank_1b", 3: "Blank_1c"}

    specs = []
    for cycle_num, names in sorted(channel_names.items(), key=lambda x: int(x[0])):
        cycle_num = int(cycle_num)
        cycle_dir = registered_dir / f"cyc{cycle_num:02d}"

        for idx, name in enumerate(names):
            # Skip DAPI, Blank, Empty channels
            if any(name.startswith(prefix) for prefix in _SKIP_PREFIXES):
                continue

            # Channel position is 0-indexed in the list, but CH positions are 1-based
            # idx=0 is CH1 (DAPI), idx=1 is CH2, idx=2 is CH3, idx=3 is CH4
            if idx not in blank_map:
                logger.warning(f"No blank mapping for index {idx} ({name}), skipping")
                continue

            blank_name = blank_map[idx]
            signal_path = cycle_dir / f"{name}.tif"
            blank_path = registered_dir / "cyc01" / f"{blank_name}.tif"

            if not signal_path.exists():
                logger.warning(f"Signal not found: {signal_path}")
                continue
            if not blank_path.exists():
                logger.warning(f"Blank not found: {blank_path}")
                continue

            specs.append(
                ChannelSpec(
                    marker_name=name,
                    cycle=cycle_num,
                    channel_position=idx + 1,
                    signal_path=signal_path,
                    blank_name=blank_name,
                    blank_path=blank_path,
                )
            )

    logger.info(f"Discovered {len(specs)} signal channels")
    return specs


def select_method(
    signal: np.ndarray,
    blank: np.ndarray,
    analysis: dict,
) -> str:
    """Select optimal subtraction method (global vs weighted) for a marker.

    Decision logic:
    - blank dominates signal (p99 ratio > 1.2) → weighted
    - structured AF (contribution > 0.3 and correlation > 0.4) → weighted
    - correlated AF across wide range (corr > 0.5 and range > 5000) → weighted
    - Otherwise → global

    Parameters
    ----------
    signal : np.ndarray
        Signal channel image.
    blank : np.ndarray
        Blank channel image.
    analysis : dict
        Output from analyze_for_subtraction().

    Returns
    -------
    str
        "global" or "weighted"
    """
    stats = analysis.get("analysis", {})
    signal_stats = stats.get("signal_stats", {})
    blank_stats = stats.get("blank_stats", {})
    correlation = stats.get("correlation", 0)
    af_contribution = stats.get("af_contribution", 0)

    signal_p99 = signal_stats.get("p99", signal_stats.get("p95", 1))
    blank_p99 = blank_stats.get("p99", blank_stats.get("p95", 0))

    # Blank dominates signal
    if signal_p99 > 0 and blank_p99 / signal_p99 > 1.2:
        logger.debug(f"Blank dominates (ratio={blank_p99/signal_p99:.2f}), selecting weighted")
        return "weighted"

    # Structured AF
    if af_contribution > 0.3 and correlation > 0.4:
        logger.debug(f"Structured AF (contrib={af_contribution:.2f}, corr={correlation:.2f})")
        return "weighted"

    # Correlated AF across wide range
    signal_range = signal_stats.get("max", 0) - signal_stats.get("min", 0)
    if correlation > 0.5 and signal_range > 5000:
        logger.debug(f"Correlated AF (corr={correlation:.2f}, range={signal_range:.0f})")
        return "weighted"

    return "global"


def process_channel(
    spec: ChannelSpec,
    method: str = "auto",
    tissue_type: str | None = None,
    tile_smooth_sigma: float = 500.0,
    dry_run: bool = False,
    output_dir: Path | None = None,
) -> ChannelResult:
    """Process a single channel with optimized parameters.

    Parameters
    ----------
    spec : ChannelSpec
        Channel specification with paths.
    method : str
        "auto", "global", or "weighted".
    tissue_type : str, optional
        Tissue type for parameter suggestions.
    tile_smooth_sigma : float
        Gaussian sigma for blank smoothing. 0 to disable.
    dry_run : bool
        If True, analyze but don't save output.
    output_dir : Path, optional
        Output directory. Defaults to project's signal_isolated dir.

    Returns
    -------
    ChannelResult
        Processing result with parameters and quality metrics.
    """
    logger.info(f"Processing {spec.marker_name} (cycle {spec.cycle}, {spec.blank_name})")

    try:
        # Load images
        signal = tifffile.imread(str(spec.signal_path))
        blank = tifffile.imread(str(spec.blank_path))

        # Smooth blank to remove tile-grid artifacts
        if tile_smooth_sigma > 0:
            blank_processed = smooth_blank_for_subtraction(blank, sigma=tile_smooth_sigma)
        else:
            blank_processed = blank

        # Analyze for parameter suggestions
        analysis = analyze_for_subtraction(
            signal, blank_processed, tissue_type=tissue_type, marker_name=spec.marker_name
        )

        # Select method
        if method == "auto":
            chosen_method = select_method(signal, blank_processed, analysis)
        else:
            chosen_method = method

        # Build result skeleton
        result = ChannelResult(
            marker_name=spec.marker_name,
            cycle=spec.cycle,
            blank_used=spec.blank_name,
            method=chosen_method,
        )

        # Store analysis info
        result.analysis = {
            "correlation": analysis["analysis"]["correlation"],
            "af_contribution": analysis["analysis"]["af_contribution"],
            "signal_p99": analysis["analysis"]["signal_stats"].get("p99", 0),
            "blank_p99": analysis["analysis"]["blank_stats"].get("p99", 0),
        }

        if dry_run:
            result.parameters = {
                "blank_clip_factor": analysis["blank_clip_factor"],
                "blank_scale_factor": analysis["blank_scale_factor"],
            }
            result.status = "dry_run"
            return result

        # Apply subtraction
        if chosen_method == "weighted":
            weighted_analysis = analyze_for_weighted_subtraction(
                signal, blank_processed, tissue_type=tissue_type, marker_name=spec.marker_name
            )
            subtracted = subtract_autofluorescence_weighted(
                signal,
                blank_processed,
                blank_clip_factor=weighted_analysis["blank_clip_factor"],
                base_scale_factor=weighted_analysis["base_scale_factor"],
                ranges=[
                    r if isinstance(r, dict) else r
                    for r in weighted_analysis.get("ranges", [])
                ],
                transition_width=weighted_analysis.get("transition_width", 0.1),
                smooth_low=weighted_analysis.get("smooth_low", False),
                low_size=weighted_analysis.get("low_size", 2),
                smooth_high=weighted_analysis.get("smooth_high", False),
                high_size=weighted_analysis.get("high_size", 2),
                erosion=weighted_analysis.get("erosion", 0),
            )
            result.parameters = {
                "blank_clip_factor": weighted_analysis["blank_clip_factor"],
                "base_scale_factor": weighted_analysis["base_scale_factor"],
                "n_ranges": weighted_analysis.get("n_ranges", 5),
                "transition_width": weighted_analysis.get("transition_width", 0.1),
            }
            # Quality metrics — extract global dict for flat structure
            weighted_quality = compute_weighted_subtraction_quality(
                signal, subtracted, blank_processed,
                ranges=weighted_analysis.get("ranges", []),
            )
            quality = weighted_quality.get("global", weighted_quality)
        else:
            subtracted = subtract_autofluorescence(
                signal,
                blank_processed,
                blank_clip_factor=analysis["blank_clip_factor"],
                blank_scale_factor=analysis["blank_scale_factor"],
                smooth_low=analysis.get("smooth_low", False),
                low_size=analysis.get("low_size", 2),
                low_percentile=analysis.get("low_percentile", 60),
                smooth_high=analysis.get("smooth_high", False),
                high_size=analysis.get("high_size", 2),
                high_percentile=analysis.get("high_percentile", 90),
                erosion=analysis.get("erosion", 0),
            )
            result.parameters = {
                "blank_clip_factor": analysis["blank_clip_factor"],
                "blank_scale_factor": analysis["blank_scale_factor"],
            }
            quality = compute_subtraction_quality(signal, subtracted, blank_processed)

        result.quality_metrics = quality

        # Compute zero percentage
        total_pixels = subtracted.size
        zero_pixels = np.sum(subtracted == 0)
        result.zero_percent = round(100.0 * zero_pixels / total_pixels, 2)

        # Flag warnings
        warnings = []
        if result.zero_percent > 70:
            warnings.append(f"high_zero={result.zero_percent:.1f}%")
        if quality.get("quality_score", 1) < 0.5:
            warnings.append(f"low_quality={quality['quality_score']:.2f}")
        p99 = float(np.percentile(subtracted, 99))
        p1 = float(np.percentile(subtracted, 1))
        if p99 - p1 < 100:
            warnings.append(f"low_range={p99-p1:.0f}")
        result.warning = "; ".join(warnings)

        # Save output
        if output_dir is not None:
            output_dir = Path(output_dir)
            output_dir.mkdir(parents=True, exist_ok=True)
            out_path = output_dir / f"{spec.marker_name}.tif"
            tifffile.imwrite(str(out_path), subtracted, compression="zlib")
            logger.info(f"Saved {out_path}")

        return result

    except Exception as e:
        logger.error(f"Error processing {spec.marker_name}: {e}")
        return ChannelResult(
            marker_name=spec.marker_name,
            cycle=spec.cycle,
            blank_used=spec.blank_name,
            method=method,
            status="error",
            warning=str(e),
        )


def process_batch(
    project_dir: str | Path,
    method: str = "auto",
    tissue_type: str | None = None,
    tile_smooth_sigma: float = 500.0,
    dry_run: bool = False,
    output_dir: str | Path | None = None,
    force: bool = False,
    channels: list[str] | None = None,
) -> BatchResult:
    """Process all signal channels in a project.

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory.
    method : str
        "auto", "global", or "weighted". Default: "auto".
    tissue_type : str, optional
        Tissue type for parameter suggestions.
    tile_smooth_sigma : float
        Gaussian sigma for blank smoothing. 0 to disable.
    dry_run : bool
        If True, analyze but don't save outputs.
    output_dir : str or Path, optional
        Output directory. Defaults to data/processed/signal_isolated.
    force : bool
        If True, overwrite existing outputs.
    channels : list[str], optional
        Process only these marker names. None = all.

    Returns
    -------
    BatchResult
        Batch processing result with per-channel details.
    """
    project_dir = Path(project_dir).resolve()

    if output_dir is None:
        output_dir = project_dir / "data" / "processed" / "signal_isolated"
    else:
        output_dir = Path(output_dir).resolve()

    specs = discover_channels(project_dir)

    # Filter to requested channels
    if channels:
        channel_set = set(channels)
        specs = [s for s in specs if s.marker_name in channel_set]

    batch = BatchResult(
        timestamp=datetime.now().isoformat(),
        method_override=method,
        tissue_type=tissue_type or "",
        tile_smooth_sigma=tile_smooth_sigma,
    )

    counts = {"total": 0, "global": 0, "weighted": 0, "skipped": 0, "error": 0}
    quality_scores = []

    for spec in specs:
        counts["total"] += 1

        # Skip existing unless force
        if not dry_run and not force:
            existing = output_dir / f"{spec.marker_name}.tif"
            if existing.exists():
                logger.info(f"Skipping {spec.marker_name} (exists, use --force to overwrite)")
                batch.channels[spec.marker_name] = ChannelResult(
                    marker_name=spec.marker_name,
                    cycle=spec.cycle,
                    blank_used=spec.blank_name,
                    method="",
                    status="skipped",
                )
                counts["skipped"] += 1
                continue

        result = process_channel(
            spec,
            method=method,
            tissue_type=tissue_type,
            tile_smooth_sigma=tile_smooth_sigma,
            dry_run=dry_run,
            output_dir=None if dry_run else output_dir,
        )

        batch.channels[spec.marker_name] = result

        if result.status == "error":
            counts["error"] += 1
        elif result.method == "global":
            counts["global"] += 1
        elif result.method == "weighted":
            counts["weighted"] += 1

        qs = result.quality_metrics.get("quality_score")
        if qs is not None:
            quality_scores.append(qs)

    counts["mean_quality"] = round(np.mean(quality_scores), 4) if quality_scores else 0.0
    batch.summary = counts

    # Write manifest
    if not dry_run:
        output_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = output_dir / "signal_isolation_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(batch.to_dict(), f, indent=2, default=str)
        logger.info(f"Manifest written to {manifest_path}")

    return batch
