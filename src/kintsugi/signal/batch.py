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
- load_legacy_recipes: Parse Notebook 3 parameter files into MarkerRecipe objects
- clean_background: Background removal (threshold + median smooth + small object removal)
"""

from __future__ import annotations

import ast
import json
import logging
import re
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

# Quality gate thresholds for over-subtraction detection
_QUALITY_GATE_THRESHOLD = 0.6
_OVER_SUBTRACTION_ZERO_PCT = 95.0

# Default percentile for post-subtraction outlier clipping
_DEFAULT_CLIP_PERCENTILE = 99.5


# =============================================================================
# Data Structures
# =============================================================================


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
    method: str  # "global", "weighted", or "recipe"
    parameters: dict = field(default_factory=dict)
    quality_metrics: dict = field(default_factory=dict)
    analysis: dict = field(default_factory=dict)
    zero_percent: float = 0.0
    status: str = "success"  # "success", "skipped", "error"
    warning: str = ""
    recipe_source: str = "auto"  # "own_recipe", "learned", "template", "auto"

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class BatchResult:
    """Result of processing all channels in a project."""

    timestamp: str = ""
    method_override: str = "auto"
    tissue_type: str = ""
    tile_smooth_sigma: float = 0.0
    recipe_dir: str = ""
    summary: dict = field(default_factory=dict)
    channels: dict[str, ChannelResult] = field(default_factory=dict)

    def to_dict(self) -> dict:
        d = {
            "timestamp": self.timestamp,
            "method_override": self.method_override,
            "tissue_type": self.tissue_type,
            "tile_smooth_sigma": self.tile_smooth_sigma,
            "recipe_dir": self.recipe_dir,
            "summary": self.summary,
            "channels": {k: v.to_dict() for k, v in self.channels.items()},
        }
        return d


# =============================================================================
# Recipe Data Structures
# =============================================================================


@dataclass
class SubtractionParams:
    """Parameters for a single subtraction step (primary or secondary)."""

    blank_name: str
    blank_clip_factor: int = 0
    blank_scale_factor: float = 1.0
    smooth_low: bool = False
    low_size: int = 2
    low_percentile: int = 60
    smooth_high: bool = False
    high_size: int = 2
    high_percentile: int = 90
    erosion: int = 0


@dataclass
class CleanParams:
    """Parameters for background cleaning step."""

    background_threshold: int = 0
    smooth: bool = False
    smooth_threshold: int = 0
    footprint: int = 1
    remove_small: bool = False
    small_size: int = 30


@dataclass
class MarkerRecipe:
    """Complete processing recipe for a marker, parsed from legacy param files."""

    marker_name: str
    primary: SubtractionParams
    second: SubtractionParams | None = None
    clean: CleanParams | None = None


# =============================================================================
# Recipe Parser
# =============================================================================


def _parse_param_value(raw: str) -> object:
    """Parse a single value from a param file line.

    Handles Python literals (int, float, bool, str, tuple) via ast.literal_eval,
    with fallback to raw string for unparseable values like dask arrays and datetime.
    """
    raw = raw.strip()
    # Skip known unparseable runtime artifacts
    if raw.startswith("dask.array<") or raw.startswith("datetime.datetime("):
        return None
    try:
        return ast.literal_eval(raw)
    except (ValueError, SyntaxError):
        return raw


def load_legacy_recipes(recipe_dir: str | Path) -> dict[str, MarkerRecipe]:
    """Parse legacy Notebook 3 parameter files into MarkerRecipe objects.

    Parameters
    ----------
    recipe_dir : str or Path
        Directory containing *_param.txt files.

    Returns
    -------
    dict[str, MarkerRecipe]
        Mapping from marker name to its processing recipe.
    """
    recipe_dir = Path(recipe_dir)
    recipes = {}

    for param_file in sorted(recipe_dir.glob("*_param.txt")):
        try:
            recipe = _parse_param_file(param_file)
            if recipe is not None:
                recipes[recipe.marker_name] = recipe
        except Exception as e:
            logger.warning(f"Failed to parse {param_file.name}: {e}")

    logger.info(f"Loaded {len(recipes)} recipes from {recipe_dir}")
    return recipes


def _parse_param_file(path: Path) -> MarkerRecipe | None:
    """Parse a single *_param.txt file into a MarkerRecipe."""
    params = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line or ":" not in line:
                continue
            # Split on first colon only
            key, _, value = line.partition(":")
            key = key.strip()
            parsed = _parse_param_value(value)
            if parsed is not None:
                params[key] = parsed

    marker_name = params.get("signal_channel")
    if not marker_name:
        logger.warning(f"No signal_channel in {path.name}")
        return None

    # Strip quotes if present
    if isinstance(marker_name, str):
        marker_name = marker_name.strip("'\"")

    blank_name = params.get("blankID", "Blank_1a")
    if isinstance(blank_name, str):
        blank_name = blank_name.strip("'\"")

    primary = SubtractionParams(
        blank_name=blank_name,
        blank_clip_factor=int(params.get("blank_clip_factor", 0)),
        blank_scale_factor=float(params.get("blank_scale_factor", 1.0)),
        smooth_low=bool(params.get("smooth_low", False)),
        low_size=int(params.get("low_size", 2)),
        low_percentile=int(params.get("low_percentile", 60)),
        smooth_high=bool(params.get("smooth_high", False)),
        high_size=int(params.get("high_size", 2)),
        high_percentile=int(params.get("high_percentile", 90)),
        erosion=int(params.get("erosion", 0)),
    )

    # Second subtraction
    second = None
    if params.get("second_subtract"):
        blank_name_2 = params.get("blankID2", "")
        if isinstance(blank_name_2, str):
            blank_name_2 = blank_name_2.strip("'\"")
        if blank_name_2:
            second = SubtractionParams(
                blank_name=blank_name_2,
                blank_clip_factor=int(params.get("blank_clip_factor_2", 0)),
                blank_scale_factor=float(params.get("blank_scale_factor_2", 1.0)),
                smooth_low=bool(params.get("smooth_low_2", False)),
                low_size=int(params.get("low_size_2", 2)),
                low_percentile=int(params.get("low_percentile_2", 60)),
                smooth_high=bool(params.get("smooth_high_2", False)),
                high_size=int(params.get("high_size_2", 2)),
                high_percentile=int(params.get("high_percentile_2", 90)),
                erosion=int(params.get("erosion_2", 0)),
            )

    # Clean parameters
    clean = None
    if params.get("clean_used"):
        clean = CleanParams(
            background_threshold=int(params.get("background_threshold", 0)),
            smooth=bool(params.get("smooth", False)),
            smooth_threshold=int(params.get("smooth_threshold", 0)),
            footprint=int(params.get("footprint", 1)),
            remove_small=bool(params.get("remove_small", False)),
            small_size=int(params.get("small_size", 30)),
        )

    return MarkerRecipe(
        marker_name=marker_name,
        primary=primary,
        second=second,
        clean=clean,
    )


# =============================================================================
# Blank Name Resolution
# =============================================================================


def _normalize_blank_name(name: str) -> str:
    """Normalize blank names from legacy param files.

    Transforms: Blank1b -> Blank_1b, Blank13c -> Blank_13c
    Passes through non-blank names unchanged (CD3e, DAPI, etc.).
    """
    # Match patterns like Blank1b, Blank13c — missing underscore
    m = re.match(r"^Blank(\d+)([a-c])$", name, re.IGNORECASE)
    if m:
        return f"Blank_{m.group(1)}{m.group(2)}"
    return name


def _strip_for_comparison(name: str) -> str:
    """Strip hyphens, underscores, and lowercase for fuzzy comparison."""
    return re.sub(r"[-_\s]", "", name).lower()


def _extract_blank_position(name: str) -> str | None:
    """Extract the position suffix (a/b/c) from a blank channel name.

    Works on both normalized and unnormalized forms:
    Blank_1b -> b, Blank_13c -> c, Blank1b -> b

    Returns None for non-blank names or blanks without a position suffix.
    """
    m = re.match(r"^Blank_?\d+([a-c])$", name, re.IGNORECASE)
    return m.group(1).lower() if m else None


def _build_channel_location_map(channel_names: dict, registered_dir: Path) -> dict[str, Path]:
    """Build a mapping from channel name to its registered TIFF path.

    Scans all cycles and all channel names in config to create a
    comprehensive lookup for blank resolution.
    """
    location_map = {}
    for cycle_num, names in channel_names.items():
        cycle_num = int(cycle_num)
        cycle_dir = registered_dir / f"cyc{cycle_num:02d}"
        for name in names:
            path = cycle_dir / f"{name}.tif"
            if path.exists():
                location_map[name] = path
    return location_map


def _resolve_blank_path(blank_name: str, location_map: dict[str, Path]) -> tuple[str, Path] | None:
    """Resolve a blank name from a recipe to an actual file path.

    Resolution chain:
    1. Normalize (Blank1b -> Blank_1b)
    2. Exact match in location_map
    3. Fuzzy match (strip hyphens/underscores, case-insensitive)

    Returns
    -------
    tuple[str, Path] or None
        (resolved_name, path) or None if not found.
    """
    normalized = _normalize_blank_name(blank_name)

    # Exact match
    if normalized in location_map:
        return normalized, location_map[normalized]

    # Fuzzy match
    stripped = _strip_for_comparison(normalized)
    for name, path in location_map.items():
        if _strip_for_comparison(name) == stripped:
            return name, path

    # Positional fallback: Blank_13b -> Blank_1b (same position, cycle 1)
    position = _extract_blank_position(normalized)
    if position:
        fallback_name = f"Blank_1{position}"
        if fallback_name in location_map:
            logger.warning(
                f"Positional remap: '{blank_name}' -> '{fallback_name}' "
                f"(original not found, using cycle 1 blank with same position)"
            )
            return fallback_name, location_map[fallback_name]
        # Try fuzzy match on fallback
        fallback_stripped = _strip_for_comparison(fallback_name)
        for name, path in location_map.items():
            if _strip_for_comparison(name) == fallback_stripped:
                logger.warning(
                    f"Positional remap: '{blank_name}' -> '{name}' "
                    f"(original not found, using cycle 1 blank with same position)"
                )
                return name, path

    return None


# =============================================================================
# Background Cleaning
# =============================================================================


def clean_background(image: np.ndarray, params: CleanParams) -> np.ndarray:
    """Remove background from processed image.

    Pure numpy reimplementation of Kutils.clean(). Three operations:
    1. Zero pixels below background_threshold
    2. Median filter in transition zone (between 0 and smooth_threshold)
    3. Remove small objects + morphological closing

    Parameters
    ----------
    image : np.ndarray
        Input image (uint16).
    params : CleanParams
        Cleaning parameters.

    Returns
    -------
    np.ndarray
        Cleaned image (uint16).
    """
    from scipy.ndimage import median_filter
    from skimage.morphology import closing, disk, remove_small_objects

    result = image.copy().astype(np.float64)

    # 1. Zero pixels below threshold
    bg_thresh = max(1, params.background_threshold)
    result[result <= bg_thresh] = 0

    # 2. Median filter in transition zone
    if params.smooth and params.smooth_threshold > 0:
        smooth_thresh = max(1, params.smooth_threshold)
        fp = max(1, params.footprint)
        transition_mask = (result > 0) & (result < smooth_thresh)
        if np.any(transition_mask):
            filtered = median_filter(result, size=fp)
            result[transition_mask] = filtered[transition_mask]

    # 3. Remove small objects + closing
    if params.remove_small:
        min_size = max(1, params.small_size)
        binary_mask = result > 0
        if np.any(binary_mask):
            cleaned_mask = remove_small_objects(binary_mask, min_size=min_size, connectivity=2)
            result[~cleaned_mask] = 0
            # Morphological closing to smooth boundaries
            closed_mask = closing(cleaned_mask, disk(1))
            # Restore pixels exposed by closing but zeroed by small object removal
            result[closed_mask & ~cleaned_mask] = image.astype(np.float64)[
                closed_mask & ~cleaned_mask
            ]

    return np.clip(result, 0, 65535).astype(np.uint16)


# =============================================================================
# Tile Smoothing
# =============================================================================


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


# =============================================================================
# Channel Discovery
# =============================================================================


def discover_channels(
    project_dir: str | Path,
    recipes: dict[str, MarkerRecipe] | None = None,
) -> list[ChannelSpec]:
    """Discover signal channels and their blank pairs from project config.

    Parses workflow/config.yaml for channel_names dict and maps each signal
    marker to its corresponding blank. When a recipe exists for a marker,
    uses the recipe's blank_name; otherwise uses positional mapping:
    - CH2 (index 1) -> Blank_1a
    - CH3 (index 2) -> Blank_1b
    - CH4 (index 3) -> Blank_1c

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory.
    recipes : dict[str, MarkerRecipe], optional
        Recipe lookup. When present, overrides positional blank mapping.

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

    # Build location map for recipe blank resolution
    location_map = _build_channel_location_map(channel_names, registered_dir)

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
            signal_path = cycle_dir / f"{name}.tif"
            if not signal_path.exists():
                logger.warning(f"Signal not found: {signal_path}")
                continue

            # Resolve blank: recipe takes priority over positional mapping
            blank_name = None
            blank_path = None

            if recipes and name in recipes:
                resolved = _resolve_blank_path(recipes[name].primary.blank_name, location_map)
                if resolved:
                    blank_name, blank_path = resolved
                else:
                    logger.warning(
                        f"Recipe blank '{recipes[name].primary.blank_name}' "
                        f"not found for {name}, falling back to positional"
                    )

            # Positional fallback
            if blank_name is None:
                if idx not in blank_map:
                    logger.warning(f"No blank mapping for index {idx} ({name}), skipping")
                    continue
                blank_name = blank_map[idx]
                blank_path = registered_dir / "cyc01" / f"{blank_name}.tif"

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


# =============================================================================
# Method Selection
# =============================================================================


def select_method(
    signal: np.ndarray,
    blank: np.ndarray,
    analysis: dict,
) -> str:
    """Select optimal subtraction method (global vs weighted) for a marker.

    Decision logic:
    - blank dominates signal (p99 ratio > 1.2) -> weighted
    - structured AF (contribution > 0.3 and correlation > 0.4) -> weighted
    - correlated AF across wide range (corr > 0.5 and range > 5000) -> weighted
    - Otherwise -> global

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
        logger.debug(f"Blank dominates (ratio={blank_p99 / signal_p99:.2f}), selecting weighted")
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


# =============================================================================
# Channel Processing
# =============================================================================


def clip_outliers(
    image: np.ndarray,
    percentile: float = _DEFAULT_CLIP_PERCENTILE,
) -> np.ndarray:
    """Clip outlier pixel values above a percentile threshold.

    Hot pixels and detector artifacts survive subtraction because the blank
    at those locations is typically normal. This clips values above the
    given percentile (computed from non-zero pixels only) to compress the
    dynamic range into the real signal distribution.

    Parameters
    ----------
    image : np.ndarray
        Subtracted image (uint16).
    percentile : float
        Percentile threshold for clipping (0-100). Values above this
        percentile are clipped to the threshold value. 0 disables clipping.
        Default: 99.5.

    Returns
    -------
    np.ndarray
        Clipped image (uint16).
    """
    if percentile <= 0 or percentile >= 100:
        return image

    # Compute percentile from non-zero pixels only to avoid
    # the zero-floor from subtraction biasing the threshold
    nonzero = image[image > 0]
    if nonzero.size == 0:
        return image

    threshold = np.percentile(nonzero, percentile)
    if threshold <= 0:
        return image

    clipped = np.minimum(image, np.uint16(min(threshold, 65535)))
    return clipped.astype(np.uint16)


def _check_over_subtraction(
    zero_percent: float,
    quality_score: float,
    quality_gate: float = _QUALITY_GATE_THRESHOLD,
    zero_threshold: float = _OVER_SUBTRACTION_ZERO_PCT,
) -> tuple[bool, str]:
    """Check if subtraction result is over-subtracted.

    Over-subtraction is detected when EITHER:
    - Zero percentage exceeds threshold (>= 95%) — image is essentially destroyed
    - Quality score is below gate AND zero percentage is high (> 70%)

    Set quality_gate <= 0 to disable the check entirely.

    Returns
    -------
    tuple[bool, str]
        (is_over_subtracted, reason)
    """
    if quality_gate <= 0:
        return False, ""
    if zero_percent >= zero_threshold:
        return True, f"zero={zero_percent:.1f}%"
    if quality_score < quality_gate and zero_percent > 70:
        return True, f"quality={quality_score:.2f}; zero={zero_percent:.1f}%"
    return False, ""


def _validate_secondary_blank(marker_name: str, blank_name: str) -> tuple[bool, str]:
    """Validate that a secondary blank is not the marker itself.

    Returns
    -------
    tuple[bool, str]
        (is_valid, reason) — False with reason if self-referential.
    """
    if _strip_for_comparison(blank_name) == _strip_for_comparison(marker_name):
        return False, "self-referential"
    return True, ""


def _validate_duplicate_blank(
    primary_blank: str,
    secondary_blank: str,
    primary_resolved: str | None = None,
    secondary_resolved: str | None = None,
) -> tuple[bool, str]:
    """Validate that secondary blank is not the same as primary blank.

    Checks both raw recipe names and resolved file names to catch cases
    where different recipe names (e.g. Blank_1b vs Blank_13b) resolve to
    the same underlying file.

    Returns
    -------
    tuple[bool, str]
        (is_valid, reason) — False with reason if duplicate.
    """
    if _strip_for_comparison(primary_blank) == _strip_for_comparison(
        secondary_blank
    ):
        return False, "duplicate_blank"
    if (
        primary_resolved
        and secondary_resolved
        and _strip_for_comparison(primary_resolved)
        == _strip_for_comparison(secondary_resolved)
    ):
        return False, "duplicate_blank_resolved"
    return True, ""


def _auto_fallback(
    spec: ChannelSpec,
    output_dir: Path | None = None,
) -> ChannelResult:
    """Run auto-analysis fallback for a channel (no recipe).

    Used when recipe-based processing produces over-subtracted results.
    """
    signal = tifffile.imread(str(spec.signal_path))
    blank = tifffile.imread(str(spec.blank_path))

    analysis = analyze_for_subtraction(signal, blank)
    chosen_method = select_method(signal, blank, analysis)

    subtracted = subtract_autofluorescence(
        signal,
        blank,
        blank_clip_factor=analysis["blank_clip_factor"],
        blank_scale_factor=analysis["blank_scale_factor"],
        smooth_low=analysis.get("smooth_low", False),
        low_size=analysis.get("low_size", 2),
        smooth_high=analysis.get("smooth_high", False),
        high_size=analysis.get("high_size", 2),
        erosion=analysis.get("erosion", 0),
    )

    quality = compute_subtraction_quality(signal, subtracted, blank)
    total_pixels = subtracted.size
    zero_pixels = np.sum(subtracted == 0)
    zero_pct = round(100.0 * zero_pixels / total_pixels, 2)

    result = ChannelResult(
        marker_name=spec.marker_name,
        cycle=spec.cycle,
        blank_used=spec.blank_name,
        method=chosen_method,
        parameters={
            "blank_clip_factor": analysis["blank_clip_factor"],
            "blank_scale_factor": analysis["blank_scale_factor"],
        },
        quality_metrics=quality,
        zero_percent=zero_pct,
        recipe_source="auto_fallback",
    )

    if output_dir is not None:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        out_path = output_dir / f"{spec.marker_name}.tif"
        tifffile.imwrite(str(out_path), subtracted, compression="zlib")

    return result


def _process_channel_recipe(
    spec: ChannelSpec,
    recipe: MarkerRecipe,
    location_map: dict[str, Path],
    output_dir: Path | None = None,
    dry_run: bool = False,
    quality_gate: float = _QUALITY_GATE_THRESHOLD,
    clip_percentile: float = _DEFAULT_CLIP_PERCENTILE,
) -> ChannelResult:
    """Process a single channel using a legacy recipe.

    Applies primary subtraction, optional second subtraction, optional
    background cleaning, and outlier clipping — matching the legacy
    Notebook 3 pipeline with hot pixel correction.
    """
    logger.info(
        f"Processing {spec.marker_name} with recipe "
        f"(blank={recipe.primary.blank_name}, "
        f"scale={recipe.primary.blank_scale_factor})"
    )

    result = ChannelResult(
        marker_name=spec.marker_name,
        cycle=spec.cycle,
        blank_used=spec.blank_name,
        method="recipe",
    )

    # Store recipe parameters in result
    params_dict = {
        "primary_blank": recipe.primary.blank_name,
        "primary_clip": recipe.primary.blank_clip_factor,
        "primary_scale": recipe.primary.blank_scale_factor,
        "primary_erosion": recipe.primary.erosion,
    }
    if recipe.second:
        params_dict["second_blank"] = recipe.second.blank_name
        params_dict["second_clip"] = recipe.second.blank_clip_factor
        params_dict["second_scale"] = recipe.second.blank_scale_factor
        params_dict["second_erosion"] = recipe.second.erosion
    if recipe.clean:
        params_dict["clean_threshold"] = recipe.clean.background_threshold
        params_dict["clean_smooth"] = recipe.clean.smooth
        params_dict["clean_remove_small"] = recipe.clean.remove_small

    result.parameters = params_dict

    if dry_run:
        result.status = "dry_run"
        return result

    try:
        # Load signal
        signal = tifffile.imread(str(spec.signal_path))

        # Load primary blank (already resolved in spec)
        blank = tifffile.imread(str(spec.blank_path))

        warnings = []

        # 1. Primary subtraction — no blank smoothing with recipes
        subtracted = subtract_autofluorescence(
            signal,
            blank,
            blank_clip_factor=recipe.primary.blank_clip_factor,
            blank_scale_factor=recipe.primary.blank_scale_factor,
            smooth_low=recipe.primary.smooth_low,
            low_size=recipe.primary.low_size,
            low_percentile=recipe.primary.low_percentile,
            smooth_high=recipe.primary.smooth_high,
            high_size=recipe.primary.high_size,
            high_percentile=recipe.primary.high_percentile,
            erosion=recipe.primary.erosion,
        )

        # 2. Second subtraction
        if recipe.second:
            # Validate secondary blank is not the marker itself
            valid, reason = _validate_secondary_blank(
                spec.marker_name, recipe.second.blank_name
            )
            if not valid:
                logger.error(
                    f"Secondary blank '{recipe.second.blank_name}' for "
                    f"{spec.marker_name} is {reason} — skipping second subtraction"
                )
                warnings.append(f"secondary_blank_{reason}")
            else:
                # Check for duplicate blank (same as primary)
                dup_valid, dup_reason = _validate_duplicate_blank(
                    recipe.primary.blank_name, recipe.second.blank_name
                )
                if not dup_valid:
                    logger.warning(
                        f"Secondary blank '{recipe.second.blank_name}' for "
                        f"{spec.marker_name} is same as primary blank "
                        f"'{recipe.primary.blank_name}' — skipping second "
                        f"subtraction"
                    )
                    warnings.append(f"secondary_blank_{dup_reason}")
                else:
                    resolved = _resolve_blank_path(
                        recipe.second.blank_name, location_map
                    )
                    if resolved:
                        resolved_name, blank2_path = resolved
                        # Re-validate resolved name against marker
                        valid2, reason2 = _validate_secondary_blank(
                            spec.marker_name, resolved_name
                        )
                        if not valid2:
                            logger.error(
                                f"Secondary blank '{recipe.second.blank_name}' "
                                f"resolved to '{resolved_name}' which is "
                                f"{reason2} for {spec.marker_name} — skipping "
                                f"second subtraction"
                            )
                            warnings.append(f"secondary_blank_{reason2}")
                        else:
                            # Check resolved names for duplicate
                            dup2_valid, dup2_reason = _validate_duplicate_blank(
                                recipe.primary.blank_name,
                                recipe.second.blank_name,
                                primary_resolved=spec.blank_name,
                                secondary_resolved=resolved_name,
                            )
                            if not dup2_valid:
                                logger.warning(
                                    f"Secondary blank "
                                    f"'{recipe.second.blank_name}' resolves "
                                    f"to '{resolved_name}' which is same as "
                                    f"primary blank '{spec.blank_name}' — "
                                    f"skipping second subtraction"
                                )
                                warnings.append(
                                    f"secondary_blank_{dup2_reason}"
                                )
                            else:
                                blank2 = tifffile.imread(str(blank2_path))
                                subtracted = subtract_autofluorescence(
                                    subtracted,
                                    blank2,
                                    blank_clip_factor=recipe.second.blank_clip_factor,
                                    blank_scale_factor=recipe.second.blank_scale_factor,
                                    smooth_low=recipe.second.smooth_low,
                                    low_size=recipe.second.low_size,
                                    low_percentile=recipe.second.low_percentile,
                                    smooth_high=recipe.second.smooth_high,
                                    high_size=recipe.second.high_size,
                                    high_percentile=recipe.second.high_percentile,
                                    erosion=recipe.second.erosion,
                                )
                                result.blank_used = (
                                    f"{spec.blank_name}+{recipe.second.blank_name}"
                                )
                    else:
                        logger.error(
                            f"Second blank '{recipe.second.blank_name}' not "
                            f"found for {spec.marker_name}, skipping second "
                            f"subtraction"
                        )

        # 3. Background cleaning
        if recipe.clean:
            subtracted = clean_background(subtracted, recipe.clean)

        # 4. Outlier clipping — hot pixels survive subtraction unchanged
        if clip_percentile > 0:
            subtracted = clip_outliers(subtracted, percentile=clip_percentile)

        # Quality metrics
        quality = compute_subtraction_quality(signal, subtracted, blank)
        result.quality_metrics = quality

        # Compute zero percentage
        total_pixels = subtracted.size
        zero_pixels = np.sum(subtracted == 0)
        result.zero_percent = round(100.0 * zero_pixels / total_pixels, 2)

        # Flag warnings
        if result.zero_percent > 70:
            warnings.append(f"high_zero={result.zero_percent:.1f}%")
        if quality.get("quality_score", 1) < 0.5:
            warnings.append(f"low_quality={quality['quality_score']:.2f}")
        p99 = float(np.percentile(subtracted, 99))
        p1 = float(np.percentile(subtracted, 1))
        if p99 - p1 < 100:
            warnings.append(f"low_range={p99 - p1:.0f}")

        # Quality gate: detect over-subtraction and attempt auto fallback
        qs = quality.get("quality_score", 1.0)
        over_sub, over_reason = _check_over_subtraction(
            result.zero_percent, qs, quality_gate=quality_gate
        )
        if over_sub:
            logger.warning(
                f"Over-subtraction detected for {spec.marker_name} "
                f"(recipe): {over_reason}"
            )
            # Attempt auto-analysis fallback
            try:
                fallback_result = _auto_fallback(
                    spec,
                    output_dir=None,  # Don't save yet
                )
                if (
                    fallback_result.zero_percent < result.zero_percent
                    and fallback_result.quality_metrics.get("quality_score", 0)
                    > qs
                ):
                    # Re-check fallback against quality gate
                    fb_qs = fallback_result.quality_metrics.get(
                        "quality_score", 0
                    )
                    fb_over, fb_reason = _check_over_subtraction(
                        fallback_result.zero_percent,
                        fb_qs,
                        quality_gate=quality_gate,
                    )
                    if fb_over:
                        logger.warning(
                            f"Auto fallback for {spec.marker_name} also "
                            f"over-subtracted: {fb_reason}"
                        )
                        fallback_result.status = "failed"
                        warnings.append(
                            f"auto_fallback_also_failed ({fb_reason})"
                        )
                    else:
                        logger.info(
                            f"Auto fallback improved {spec.marker_name}: "
                            f"zero {result.zero_percent:.1f}% -> "
                            f"{fallback_result.zero_percent:.1f}%, "
                            f"quality {qs:.2f} -> {fb_qs:.2f}"
                        )
                        warnings.append(
                            f"auto_fallback (recipe: {over_reason})"
                        )

                    fallback_result.recipe_source = "auto_fallback"
                    fallback_result.warning = "; ".join(warnings)

                    # Save fallback output
                    if output_dir is not None:
                        output_dir = Path(output_dir)
                        output_dir.mkdir(parents=True, exist_ok=True)
                        # Re-process with output to save
                        fallback_saved = _auto_fallback(
                            spec, output_dir=output_dir
                        )
                        fallback_saved.recipe_source = "auto_fallback"
                        fallback_saved.warning = fallback_result.warning
                    return fallback_result
                else:
                    logger.warning(
                        f"Auto fallback did not improve {spec.marker_name}, "
                        f"keeping recipe result with status=failed"
                    )
                    result.status = "failed"
                    warnings.append(f"over_subtracted ({over_reason})")
            except Exception as fb_err:
                logger.warning(
                    f"Auto fallback failed for {spec.marker_name}: {fb_err}"
                )
                result.status = "failed"
                warnings.append(f"over_subtracted ({over_reason})")

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
        logger.error(f"Error processing {spec.marker_name} with recipe: {e}")
        return ChannelResult(
            marker_name=spec.marker_name,
            cycle=spec.cycle,
            blank_used=spec.blank_name,
            method="recipe",
            status="error",
            warning=str(e),
        )


def process_channel(
    spec: ChannelSpec,
    method: str = "auto",
    tissue_type: str | None = None,
    tile_smooth_sigma: float = 0.0,
    dry_run: bool = False,
    output_dir: Path | None = None,
    recipe: MarkerRecipe | None = None,
    location_map: dict[str, Path] | None = None,
    quality_gate: float = _QUALITY_GATE_THRESHOLD,
    clip_percentile: float = _DEFAULT_CLIP_PERCENTILE,
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
    recipe : MarkerRecipe, optional
        Legacy recipe for this marker. When provided, uses recipe pipeline.
    location_map : dict[str, Path], optional
        Channel name to path map for resolving recipe blank names.
    quality_gate : float
        Minimum quality score threshold. Channels below this with >95%
        zeros are marked as "failed". Default: 0.6.
    clip_percentile : float
        Percentile for post-subtraction outlier clipping. Hot pixels that
        survive subtraction are clipped to this percentile of the non-zero
        distribution. 0 disables clipping. Default: 99.5.

    Returns
    -------
    ChannelResult
        Processing result with parameters and quality metrics.
    """
    # Recipe path: use legacy multi-step pipeline
    if recipe is not None:
        return _process_channel_recipe(
            spec, recipe, location_map or {}, output_dir, dry_run,
            quality_gate=quality_gate,
            clip_percentile=clip_percentile,
        )

    # Auto-analysis path (original behavior)
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
                    r if isinstance(r, dict) else r for r in weighted_analysis.get("ranges", [])
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
                signal,
                subtracted,
                blank_processed,
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

        # Outlier clipping — hot pixels survive subtraction unchanged
        if clip_percentile > 0:
            subtracted = clip_outliers(subtracted, percentile=clip_percentile)

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
            warnings.append(f"low_range={p99 - p1:.0f}")

        # Quality gate: detect over-subtraction (no fallback in auto path)
        qs = quality.get("quality_score", 1.0)
        over_sub, over_reason = _check_over_subtraction(
            result.zero_percent, qs, quality_gate=quality_gate
        )
        if over_sub:
            logger.warning(
                f"Over-subtraction detected for {spec.marker_name} "
                f"(auto/{chosen_method}): {over_reason}"
            )
            result.status = "failed"
            warnings.append(f"over_subtracted ({over_reason})")

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


# =============================================================================
# Batch Processing
# =============================================================================


def _record_to_learning_db(
    engine,
    tissue_type: str,
    marker_name: str,
    recipe: MarkerRecipe,
    quality_score: float,
):
    """Record recipe outcomes to parameter learning DB."""
    # Record primary subtraction params
    sub_params = {
        "blank_name": recipe.primary.blank_name,
        "blank_clip_factor": recipe.primary.blank_clip_factor,
        "blank_scale_factor": recipe.primary.blank_scale_factor,
        "smooth_low": recipe.primary.smooth_low,
        "low_size": recipe.primary.low_size,
        "smooth_high": recipe.primary.smooth_high,
        "high_size": recipe.primary.high_size,
        "erosion": recipe.primary.erosion,
        "algorithm_version": "recipe_v1",
    }
    if recipe.second:
        sub_params["second_blank"] = recipe.second.blank_name
        sub_params["second_clip_factor"] = recipe.second.blank_clip_factor
        sub_params["second_scale_factor"] = recipe.second.blank_scale_factor
        sub_params["second_erosion"] = recipe.second.erosion

    engine.record_parameters(
        tissue_type=tissue_type,
        marker_name=marker_name,
        operation="blank_subtraction",
        parameters=sub_params,
        quality_after=quality_score,
        user_approved=True,
        algorithm_version="recipe_v1",
    )

    # Record clean params as separate operation
    if recipe.clean:
        clean_params = {
            "background_threshold": recipe.clean.background_threshold,
            "smooth": recipe.clean.smooth,
            "smooth_threshold": recipe.clean.smooth_threshold,
            "footprint": recipe.clean.footprint,
            "remove_small": recipe.clean.remove_small,
            "small_size": recipe.clean.small_size,
            "algorithm_version": "recipe_v1",
        }
        engine.record_parameters(
            tissue_type=tissue_type,
            marker_name=marker_name,
            operation="clean_background",
            parameters=clean_params,
            quality_after=quality_score,
            user_approved=True,
            algorithm_version="recipe_v1",
        )


def process_batch(
    project_dir: str | Path,
    method: str = "auto",
    tissue_type: str | None = None,
    tile_smooth_sigma: float = 0.0,
    dry_run: bool = False,
    output_dir: str | Path | None = None,
    force: bool = False,
    channels: list[str] | None = None,
    recipe_dir: str | Path | None = None,
    learn: bool = True,
    quality_gate: float = _QUALITY_GATE_THRESHOLD,
    clip_percentile: float = _DEFAULT_CLIP_PERCENTILE,
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
        Gaussian sigma for blank smoothing. 0 to disable. Default: 0.
    dry_run : bool
        If True, analyze but don't save outputs.
    output_dir : str or Path, optional
        Output directory. Defaults to data/processed/signal_isolated.
    force : bool
        If True, overwrite existing outputs.
    channels : list[str], optional
        Process only these marker names. None = all.
    recipe_dir : str or Path, optional
        Directory containing *_param.txt legacy recipe files.
    learn : bool
        If True, record outcomes to parameter learning DB. Default: True.
    quality_gate : float
        Minimum quality score threshold. Channels below this with >95%
        zeros are marked as "failed". Default: 0.6.
    clip_percentile : float
        Percentile for post-subtraction outlier clipping. 0 to disable.
        Default: 99.5.

    Returns
    -------
    BatchResult
        Batch processing result with per-channel details.
    """
    import yaml

    project_dir = Path(project_dir).resolve()

    if output_dir is None:
        output_dir = project_dir / "data" / "processed" / "signal_isolated"
    else:
        output_dir = Path(output_dir).resolve()

    # Load recipes if provided
    recipes = None
    if recipe_dir is not None:
        recipes = load_legacy_recipes(recipe_dir)
        logger.info(f"Loaded {len(recipes)} recipes from {recipe_dir}")

    # Build location map for blank resolution
    config_path = project_dir / "workflow" / "config.yaml"
    location_map = {}
    if config_path.exists():
        with open(config_path) as f:
            config = yaml.safe_load(f)
        channel_names = config.get("channel_names", {})
        registered_dir = project_dir / "data" / "processed" / "registered"
        location_map = _build_channel_location_map(channel_names, registered_dir)

    specs = discover_channels(project_dir, recipes=recipes)

    # Filter to requested channels
    if channels:
        channel_set = set(channels)
        specs = [s for s in specs if s.marker_name in channel_set]

    batch = BatchResult(
        timestamp=datetime.now().isoformat(),
        method_override=method,
        tissue_type=tissue_type or "",
        tile_smooth_sigma=tile_smooth_sigma,
        recipe_dir=str(recipe_dir) if recipe_dir else "",
    )

    # Initialize learning engine if enabled
    engine = None
    if learn and not dry_run and tissue_type:
        try:
            from kintsugi.claude.parameter_learning import ParameterLearningEngine

            engine = ParameterLearningEngine(project_path=project_dir)
        except Exception as e:
            logger.warning(f"Could not initialize learning engine: {e}")

    counts = {
        "total": 0,
        "global": 0,
        "weighted": 0,
        "recipe": 0,
        "skipped": 0,
        "error": 0,
        "failed": 0,
    }
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

        # Determine if we have a recipe for this marker
        marker_recipe = recipes.get(spec.marker_name) if recipes else None

        result = process_channel(
            spec,
            method=method,
            tissue_type=tissue_type,
            tile_smooth_sigma=tile_smooth_sigma,
            dry_run=dry_run,
            output_dir=None if dry_run else output_dir,
            recipe=marker_recipe,
            location_map=location_map,
            quality_gate=quality_gate,
            clip_percentile=clip_percentile,
        )

        batch.channels[spec.marker_name] = result

        if result.status == "error":
            counts["error"] += 1
        elif result.status == "failed":
            counts["failed"] += 1
            # Still count the method used
            if result.method == "recipe":
                counts["recipe"] += 1
            elif result.method == "global":
                counts["global"] += 1
            elif result.method == "weighted":
                counts["weighted"] += 1
        elif result.method == "recipe":
            counts["recipe"] += 1
        elif result.method == "global":
            counts["global"] += 1
        elif result.method == "weighted":
            counts["weighted"] += 1

        qs = result.quality_metrics.get("quality_score")
        if qs is not None:
            quality_scores.append(qs)

        # Record to learning DB
        if engine and marker_recipe and result.status == "success" and qs is not None:
            try:
                _record_to_learning_db(engine, tissue_type, spec.marker_name, marker_recipe, qs)
            except Exception as e:
                logger.warning(f"Failed to record {spec.marker_name} to learning DB: {e}")

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
