"""
Multi-Project Batch Signal Isolation

Orchestrates autofluorescence subtraction across multiple KINTSUGI projects
with cascading recipe resolution:
  1. Project's own recipes (Processing_parameters/*.txt)
  2. Learned DB parameters (from previous successful runs)
  3. Template recipes (from a reference project, e.g. CX_19-001)
  4. Auto-analysis (global vs weighted selection per marker)

Key functions:
- discover_projects: Find projects with registered images ready for isolation
- parse_tissue_type: Infer tissue type from project name
- resolve_recipes_for_project: Cascading recipe resolution per marker
- learned_params_to_recipe: Convert learning DB params to MarkerRecipe
- process_all_projects: Sequential multi-project orchestrator
"""

from __future__ import annotations

import json
import logging
import re
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .batch import BatchResult, MarkerRecipe

logger = logging.getLogger("kintsugi.signal.batch_multi")

# Default template recipe project (has 26 validated recipes)
DEFAULT_TEMPLATE_PROJECT = "CX_19-001_SP_CC2-A28"

# Directories to search for legacy recipe files, in priority order (relative to project)
_RECIPE_SEARCH_DIRS = [
    "data/processed/Processing_parameters",
    "Processing_parameters",
    "notebooks/Processing_parameters",
]

# External base directories containing per-project recipe folders
# Structure: {base}/{dataset_name}/Processing_parameters/*_param.txt
_EXTERNAL_RECIPE_BASES = [
    "/mnt/ahc_share/PERSONAL/smith6jt/Processed_CODEX",
]


# =============================================================================
# Data Structures
# =============================================================================


@dataclass
class ProjectInfo:
    """Metadata for a discovered project."""

    name: str
    path: Path
    tissue_type: str
    n_registered_channels: int = 0
    has_manifest: bool = False
    has_own_recipes: bool = False
    own_recipe_dir: str = ""
    status: str = "pending"  # "pending", "completed", "skipped", "error"

    def to_dict(self) -> dict:
        d = asdict(self)
        d["path"] = str(self.path)
        return d


@dataclass
class ResolvedRecipes:
    """Result of cascading recipe resolution for a project."""

    recipes: dict  # marker_name -> MarkerRecipe
    sources: dict[str, str]  # marker_name -> source tier name

    def summary(self) -> dict[str, int]:
        """Count markers per source tier."""
        counts: dict[str, int] = {}
        for source in self.sources.values():
            counts[source] = counts.get(source, 0) + 1
        return counts


@dataclass
class MultiProjectResult:
    """Result of processing multiple projects."""

    timestamp: str = ""
    projects_dir: str = ""
    template_recipe_dir: str = ""
    total_projects: int = 0
    completed: int = 0
    skipped: int = 0
    errors: int = 0
    project_results: dict[str, dict] = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


# =============================================================================
# Tissue Type Parsing
# =============================================================================

# Patterns ordered from most specific to least specific
_TISSUE_PATTERNS: list[tuple[str, str]] = [
    (r"(?:_SP_|spleen|splenic)", "spleen"),
    (r"(?:_LN_|lymph.?node|lymph)", "lymph_node"),
    (r"(?:_TH_|thymus|thymic)", "thymus"),
    (r"(?:pancrea)", "pancreas"),
]


def parse_tissue_type(project_name: str, project_path: Path | None = None) -> str:
    """Infer tissue type from project name using regex patterns.

    Falls back to experiment.json 'name' field, then 'unknown'.

    Parameters
    ----------
    project_name : str
        Project directory name.
    project_path : Path, optional
        Project directory path for experiment.json fallback.

    Returns
    -------
    str
        Normalized tissue type: 'spleen', 'lymph_node', 'thymus', etc.
    """
    for pattern, tissue in _TISSUE_PATTERNS:
        if re.search(pattern, project_name, re.IGNORECASE):
            return tissue

    # Fallback: check experiment.json
    if project_path is not None:
        for meta_dir in ["meta", "."]:
            exp_path = project_path / meta_dir / "experiment.json"
            if exp_path.exists():
                try:
                    with open(exp_path) as f:
                        exp = json.load(f)
                    name = exp.get("name", "")
                    if name:
                        for pattern, tissue in _TISSUE_PATTERNS:
                            if re.search(pattern, name, re.IGNORECASE):
                                return tissue
                except (json.JSONDecodeError, OSError):
                    pass

    logger.warning(f"Could not determine tissue type for '{project_name}', using 'unknown'")
    return "unknown"


# =============================================================================
# External Recipe Lookup
# =============================================================================


def _find_external_recipes(project_name: str) -> str:
    """Search external recipe bases for matching recipes.

    Tries exact match first, then fuzzy substring matching (min 5 chars)
    to handle mixed naming conventions (e.g. KINTSUGI '..._CC2-A28' vs
    Processed_CODEX '1904CC2-A28').

    Returns recipe directory path or empty string.
    """
    for base in _EXTERNAL_RECIPE_BASES:
        base_path = Path(base)
        if not base_path.is_dir():
            continue

        # Exact match
        candidate = base_path / project_name / "Processing_parameters"
        if candidate.is_dir() and list(candidate.glob("*_param.txt")):
            return str(candidate)

        # Fuzzy: check if any external subdir name overlaps with project name
        # Use bidirectional substring match with min length 5 to avoid false positives
        project_lower = project_name.lower()
        for entry in sorted(base_path.iterdir()):
            if not entry.is_dir() or entry.name == project_name:
                continue
            entry_lower = entry.name.lower()
            if len(entry_lower) >= 5 and (
                entry_lower in project_lower or project_lower in entry_lower
            ):
                candidate = entry / "Processing_parameters"
                if candidate.is_dir() and list(candidate.glob("*_param.txt")):
                    logger.info(f"Fuzzy recipe match: '{project_name}' -> '{entry.name}' in {base}")
                    return str(candidate)

    return ""


# =============================================================================
# Project Discovery
# =============================================================================


def discover_projects(
    projects_dir: str | Path,
    force: bool = False,
    dataset: str | None = None,
) -> list[ProjectInfo]:
    """Find projects with registered images ready for signal isolation.

    Parameters
    ----------
    projects_dir : str or Path
        Parent directory containing project subdirectories.
    force : bool
        If True, include projects that already have signal_isolation_manifest.json.
    dataset : str, optional
        If provided, only return the named project.

    Returns
    -------
    list[ProjectInfo]
        Projects sorted by name, filtered to those with config + registered images.
    """
    import yaml

    projects_dir = Path(projects_dir).resolve()
    projects = []

    for project_path in sorted(projects_dir.iterdir()):
        if not project_path.is_dir():
            continue

        name = project_path.name

        # Filter to single dataset if requested
        if dataset and name != dataset:
            continue

        # Skip src_ prefixed directories that lack config (raw staging dirs)
        config_path = project_path / "workflow" / "config.yaml"
        if not config_path.exists():
            continue

        # Must have channel_names in config
        try:
            with open(config_path) as f:
                config = yaml.safe_load(f)
            if not config or not config.get("channel_names"):
                continue
        except (yaml.YAMLError, OSError):
            continue

        # Must have registered images
        registered_dir = project_path / "data" / "processed" / "registered"
        if not registered_dir.exists():
            continue

        # Count registered TIFs
        tifs = list(registered_dir.glob("cyc*/*.tif"))
        if not tifs:
            continue

        # Check manifest
        manifest_path = (
            project_path
            / "data"
            / "processed"
            / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        has_manifest = manifest_path.exists()

        if has_manifest and not force:
            continue

        # Check for own recipes (local first, then external)
        own_recipe_dir = ""
        has_own_recipes = False
        for rel_dir in _RECIPE_SEARCH_DIRS:
            candidate = project_path / rel_dir
            if candidate.exists() and list(candidate.glob("*_param.txt")):
                has_own_recipes = True
                own_recipe_dir = str(candidate)
                break

        if not has_own_recipes:
            ext_dir = _find_external_recipes(name)
            if ext_dir:
                has_own_recipes = True
                own_recipe_dir = ext_dir

        tissue_type = parse_tissue_type(name, project_path)

        projects.append(
            ProjectInfo(
                name=name,
                path=project_path,
                tissue_type=tissue_type,
                n_registered_channels=len(tifs),
                has_manifest=has_manifest,
                has_own_recipes=has_own_recipes,
                own_recipe_dir=own_recipe_dir,
            )
        )

    logger.info(f"Discovered {len(projects)} projects ready for signal isolation")
    return projects


# =============================================================================
# Recipe Resolution
# =============================================================================


def _get_marker_names_from_config(project_path: Path) -> list[str]:
    """Extract signal marker names from project config.yaml."""
    import yaml

    config_path = project_path / "workflow" / "config.yaml"
    if not config_path.exists():
        return []

    with open(config_path) as f:
        config = yaml.safe_load(f)

    channel_names = config.get("channel_names", {})
    skip_prefixes = ("DAPI", "Blank", "Empty")
    markers = []
    for _cycle, names in channel_names.items():
        for name in names:
            if not any(name.startswith(p) for p in skip_prefixes):
                if name not in markers:
                    markers.append(name)
    return markers


def learned_params_to_recipe(params: dict, marker_name: str) -> MarkerRecipe:
    """Convert learned DB parameters back to a MarkerRecipe.

    Parameters
    ----------
    params : dict
        Recommended parameters from ParameterLearningEngine.recommend_parameters().
        Expected keys: blank_name, blank_clip_factor, blank_scale_factor, etc.
    marker_name : str
        Marker name for the recipe.

    Returns
    -------
    MarkerRecipe
        Recipe constructed from learned parameters.
    """
    from .batch import MarkerRecipe, SubtractionParams

    rec = params.get("recommended_parameters", params)
    if rec is None:
        rec = {}

    primary = SubtractionParams(
        blank_name=rec.get("blank_name", "Blank_1a"),
        blank_clip_factor=int(rec.get("blank_clip_factor", 0)),
        blank_scale_factor=float(rec.get("blank_scale_factor", 1.0)),
        smooth_low=bool(rec.get("smooth_low", False)),
        low_size=int(rec.get("low_size", 2)),
        low_percentile=int(rec.get("low_percentile", 60)),
        smooth_high=bool(rec.get("smooth_high", False)),
        high_size=int(rec.get("high_size", 2)),
        high_percentile=int(rec.get("high_percentile", 90)),
        erosion=int(rec.get("erosion", 0)),
    )

    # Second subtraction from learned params
    second = None
    if rec.get("second_blank"):
        second = SubtractionParams(
            blank_name=rec["second_blank"],
            blank_clip_factor=int(rec.get("second_clip_factor", 0)),
            blank_scale_factor=float(rec.get("second_scale_factor", 1.0)),
            erosion=int(rec.get("second_erosion", 0)),
        )

    return MarkerRecipe(
        marker_name=marker_name,
        primary=primary,
        second=second,
        clean=None,  # Clean params stored as separate operation in learning DB
    )


def resolve_recipes_for_project(
    project: ProjectInfo,
    template_recipe_dir: str | Path | None = None,
    min_confidence: float = 0.6,
    projects_dir: Path | None = None,
) -> ResolvedRecipes:
    """Resolve recipes for each marker using cascading priority.

    Tier 1: Project's own recipe files
    Tier 2: Parameter learning DB (tissue + marker match, confidence >= threshold)
    Tier 3: Template recipes from reference project
    Tier 4: Auto-analysis (no recipe — process_batch will select method automatically)

    Parameters
    ----------
    project : ProjectInfo
        Target project.
    template_recipe_dir : str or Path, optional
        Path to template recipe directory. Defaults to CX_19-001's recipes.
    min_confidence : float
        Minimum confidence for learned DB parameters. Default: 0.6.
    projects_dir : Path, optional
        Parent projects directory (for resolving default template path).

    Returns
    -------
    ResolvedRecipes
        Per-marker recipes and their source tier.
    """
    from .batch import load_legacy_recipes

    markers = _get_marker_names_from_config(project.path)
    recipes = {}
    sources = {}

    # --- Tier 1: Own recipes ---
    own_recipes = {}
    if project.has_own_recipes and project.own_recipe_dir:
        own_recipes = load_legacy_recipes(project.own_recipe_dir)
        for marker in markers:
            if marker in own_recipes:
                recipes[marker] = own_recipes[marker]
                sources[marker] = "own_recipe"

    # --- Tier 2: Learned DB ---
    remaining = [m for m in markers if m not in recipes]
    if remaining and project.tissue_type != "unknown":
        try:
            from kintsugi.claude.parameter_learning import ParameterLearningEngine

            engine = ParameterLearningEngine(project_path=project.path)
            for marker in remaining:
                rec = engine.recommend_parameters(
                    tissue_type=project.tissue_type,
                    marker_name=marker,
                    operation="blank_subtraction",
                )
                confidence = rec.get("confidence", 0.0)
                if confidence >= min_confidence and rec.get("recommended_parameters"):
                    try:
                        recipe = learned_params_to_recipe(rec, marker)
                        recipes[marker] = recipe
                        sources[marker] = "learned"
                    except Exception as e:
                        logger.debug(f"Failed to convert learned params for {marker}: {e}")
        except Exception as e:
            logger.warning(f"Could not query learning DB: {e}")

    # --- Tier 3: Template recipes ---
    remaining = [m for m in markers if m not in recipes]
    if remaining:
        template_dir = _resolve_template_dir(template_recipe_dir, projects_dir)
        if template_dir:
            template_recipes = load_legacy_recipes(template_dir)
            for marker in remaining:
                if marker in template_recipes:
                    recipes[marker] = template_recipes[marker]
                    sources[marker] = "template"

    # --- Tier 4: Auto-analysis (no recipe needed) ---
    for marker in markers:
        if marker not in sources:
            sources[marker] = "auto"

    return ResolvedRecipes(recipes=recipes, sources=sources)


def _resolve_template_dir(
    template_recipe_dir: str | Path | None,
    projects_dir: Path | None,
) -> Path | None:
    """Resolve template recipe directory path."""
    if template_recipe_dir:
        p = Path(template_recipe_dir)
        if p.exists() and list(p.glob("*_param.txt")):
            return p
        logger.warning(f"Template recipe dir not found or empty: {p}")
        return None

    # Default: CX_19-001's recipes
    if projects_dir:
        default = (
            projects_dir / DEFAULT_TEMPLATE_PROJECT / "data" / "processed" / "Processing_parameters"
        )
        if default.exists() and list(default.glob("*_param.txt")):
            return default

    return None


# =============================================================================
# Processing Orchestration
# =============================================================================


def process_all_projects(
    projects_dir: str | Path,
    dry_run: bool = False,
    dataset: str | None = None,
    force: bool = False,
    template_recipe_dir: str | Path | None = None,
    skip_qc: bool = False,
    learn: bool = True,
    method: str = "auto",
    min_confidence: float = 0.6,
    tile_smooth_sigma: float = 0.0,
) -> MultiProjectResult:
    """Process signal isolation across multiple projects.

    Sequential processing with error isolation per project.

    Parameters
    ----------
    projects_dir : str or Path
        Parent directory containing project subdirectories.
    dry_run : bool
        Preview without processing.
    dataset : str, optional
        Process only this named project.
    force : bool
        Reprocess projects that already have manifests.
    template_recipe_dir : str or Path, optional
        Override template recipe directory.
    skip_qc : bool
        Skip QC page generation.
    learn : bool
        Record outcomes to parameter learning DB.
    method : str
        Override method for auto-analysis channels: 'auto', 'global', 'weighted'.
    min_confidence : float
        Minimum confidence for learned DB parameters.
    tile_smooth_sigma : float
        Gaussian sigma for blank smoothing.

    Returns
    -------
    MultiProjectResult
        Cross-project results summary.
    """
    projects_dir = Path(projects_dir).resolve()
    projects = discover_projects(projects_dir, force=force, dataset=dataset)

    result = MultiProjectResult(
        timestamp=datetime.now().isoformat(),
        projects_dir=str(projects_dir),
        template_recipe_dir=str(template_recipe_dir or ""),
        total_projects=len(projects),
    )

    if not projects:
        logger.info("No projects to process")
        return result

    for i, project in enumerate(projects, 1):
        logger.info(
            f"\n{'=' * 60}\n"
            f"[{i}/{len(projects)}] {project.name} ({project.tissue_type})\n"
            f"{'=' * 60}"
        )

        try:
            # Resolve recipes with cascade
            resolved = resolve_recipes_for_project(
                project,
                template_recipe_dir=template_recipe_dir,
                min_confidence=min_confidence,
                projects_dir=projects_dir,
            )

            source_summary = resolved.summary()
            logger.info(f"Recipe sources: {source_summary}")

            # Build per-channel recipe source tracking
            recipe_source_map = resolved.sources

            # Call process_batch with resolved recipes
            batch_result = _process_project_with_resolved_recipes(
                project=project,
                resolved=resolved,
                method=method,
                tile_smooth_sigma=tile_smooth_sigma,
                dry_run=dry_run,
                force=force,
                learn=learn,
            )

            # Tag recipe sources on channel results
            for marker, ch_result in batch_result.channels.items():
                ch_result.recipe_source = recipe_source_map.get(marker, "auto")

            # Generate QC
            if not dry_run and not skip_qc and batch_result.summary.get("total", 0) > 0:
                try:
                    from .isolation_qc import generate_qc_pages

                    generate_qc_pages(project.path)
                except Exception as e:
                    logger.warning(f"QC generation failed for {project.name}: {e}")

            project.status = "completed" if not dry_run else "dry_run"
            result.completed += 1

            result.project_results[project.name] = {
                "tissue_type": project.tissue_type,
                "status": project.status,
                "recipe_sources": source_summary,
                "summary": batch_result.summary,
                "warnings": [
                    f"{m}: {ch.warning}" for m, ch in batch_result.channels.items() if ch.warning
                ],
            }

        except Exception as e:
            logger.error(f"Error processing {project.name}: {e}")
            project.status = "error"
            result.errors += 1
            result.project_results[project.name] = {
                "tissue_type": project.tissue_type,
                "status": "error",
                "error": str(e),
            }

    result.skipped = result.total_projects - result.completed - result.errors

    # Write cross-project report
    if not dry_run:
        report_path = projects_dir / f"batch_isolation_report_{_timestamp_slug()}.json"
        with open(report_path, "w") as f:
            json.dump(result.to_dict(), f, indent=2, default=str)
        logger.info(f"Cross-project report: {report_path}")

    return result


def _process_project_with_resolved_recipes(
    project: ProjectInfo,
    resolved: ResolvedRecipes,
    method: str = "auto",
    tile_smooth_sigma: float = 0.0,
    dry_run: bool = False,
    force: bool = False,
    learn: bool = True,
) -> BatchResult:
    """Process a single project using pre-resolved recipes.

    Markers with recipes (tiers 1-3) use recipe path; markers without (tier 4)
    use auto-analysis path. Both go through process_batch, but we split them
    into two calls to correctly route recipes vs auto.
    """
    from .batch import (
        BatchResult,
        process_batch,
    )

    # Markers with recipes vs auto
    recipe_markers = {m for m, s in resolved.sources.items() if s != "auto"}
    auto_markers = {m for m, s in resolved.sources.items() if s == "auto"}

    # We need a temporary recipe dir if we have resolved recipes from mixed sources
    # Instead, call process_batch twice: once with recipes, once without
    combined_result = BatchResult(
        timestamp=datetime.now().isoformat(),
        method_override=method,
        tissue_type=project.tissue_type,
        tile_smooth_sigma=tile_smooth_sigma,
    )

    # Part 1: Process recipe-sourced markers
    if recipe_markers and resolved.recipes:
        # Write resolved recipes to a temp dir so process_batch can load them
        import tempfile

        with tempfile.TemporaryDirectory(prefix="kintsugi_recipes_") as tmpdir:
            _write_recipes_as_param_files(resolved.recipes, Path(tmpdir))

            batch_r = process_batch(
                project.path,
                method=method,
                tissue_type=project.tissue_type,
                tile_smooth_sigma=tile_smooth_sigma,
                dry_run=dry_run,
                force=force,
                channels=list(recipe_markers),
                recipe_dir=tmpdir,
                learn=learn,
            )
            combined_result.channels.update(batch_r.channels)

    # Part 2: Process auto-analysis markers
    if auto_markers:
        batch_a = process_batch(
            project.path,
            method=method,
            tissue_type=project.tissue_type,
            tile_smooth_sigma=tile_smooth_sigma,
            dry_run=dry_run,
            force=force,
            channels=list(auto_markers),
            learn=learn,
        )
        combined_result.channels.update(batch_a.channels)

    # Rebuild summary
    counts = {
        "total": 0,
        "global": 0,
        "weighted": 0,
        "recipe": 0,
        "skipped": 0,
        "error": 0,
    }
    quality_scores = []
    for ch in combined_result.channels.values():
        counts["total"] += 1
        if ch.status == "skipped":
            counts["skipped"] += 1
        elif ch.status == "error":
            counts["error"] += 1
        elif ch.method == "recipe":
            counts["recipe"] += 1
        elif ch.method == "global":
            counts["global"] += 1
        elif ch.method == "weighted":
            counts["weighted"] += 1
        qs = ch.quality_metrics.get("quality_score")
        if qs is not None:
            quality_scores.append(qs)

    import numpy as np

    counts["mean_quality"] = round(float(np.mean(quality_scores)), 4) if quality_scores else 0.0
    combined_result.summary = counts

    # Write manifest (process_batch writes per-call, but we need a combined one)
    if not dry_run:
        output_dir = project.path / "data" / "processed" / "signal_isolated"
        output_dir.mkdir(parents=True, exist_ok=True)
        manifest_path = output_dir / "signal_isolation_manifest.json"
        with open(manifest_path, "w") as f:
            json.dump(combined_result.to_dict(), f, indent=2, default=str)
        logger.info(f"Manifest written to {manifest_path}")

    return combined_result


def _write_recipes_as_param_files(recipes: dict, output_dir: Path) -> None:
    """Write MarkerRecipe objects as legacy *_param.txt files.

    This allows passing resolved recipes (from learned DB or template)
    to process_batch via its existing recipe_dir interface.
    """
    for marker_name, recipe in recipes.items():
        lines = [
            f"signal_channel: '{marker_name}'",
            f"blankID: '{recipe.primary.blank_name}'",
            f"blank_clip_factor: {recipe.primary.blank_clip_factor}",
            f"blank_scale_factor: {recipe.primary.blank_scale_factor}",
            f"smooth_low: {recipe.primary.smooth_low}",
            f"low_size: {recipe.primary.low_size}",
            f"low_percentile: {recipe.primary.low_percentile}",
            f"smooth_high: {recipe.primary.smooth_high}",
            f"high_size: {recipe.primary.high_size}",
            f"high_percentile: {recipe.primary.high_percentile}",
            f"erosion: {recipe.primary.erosion}",
        ]

        if recipe.second:
            lines.extend(
                [
                    "second_subtract: True",
                    f"blankID2: '{recipe.second.blank_name}'",
                    f"blank_clip_factor_2: {recipe.second.blank_clip_factor}",
                    f"blank_scale_factor_2: {recipe.second.blank_scale_factor}",
                    f"smooth_low_2: {recipe.second.smooth_low}",
                    f"low_size_2: {recipe.second.low_size}",
                    f"low_percentile_2: {recipe.second.low_percentile}",
                    f"smooth_high_2: {recipe.second.smooth_high}",
                    f"high_size_2: {recipe.second.high_size}",
                    f"high_percentile_2: {recipe.second.high_percentile}",
                    f"erosion_2: {recipe.second.erosion}",
                ]
            )

        if recipe.clean:
            lines.extend(
                [
                    "clean_used: True",
                    f"background_threshold: {recipe.clean.background_threshold}",
                    f"smooth: {recipe.clean.smooth}",
                    f"smooth_threshold: {recipe.clean.smooth_threshold}",
                    f"footprint: {recipe.clean.footprint}",
                    f"remove_small: {recipe.clean.remove_small}",
                    f"small_size: {recipe.clean.small_size}",
                ]
            )

        param_file = output_dir / f"{marker_name}_param.txt"
        param_file.write_text("\n".join(lines) + "\n")


def _timestamp_slug() -> str:
    """Generate a short timestamp for filenames."""
    return datetime.now().strftime("%Y%m%d_%H%M%S")
