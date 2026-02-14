"""
TissUUmaps export pipeline for KINTSUGI projects.

Converts processed outputs (registered images, segmentation masks, phenotyping
data) to web-optimized DZI tile pyramids and generates .tmap project files for
TissUUmaps visualization.

Usage:
    from kintsugi.export import prepare_export, deploy_export, get_export_status

    # Prepare export (DZI conversion + .tmap generation)
    result = prepare_export(project_dir)

    # Deploy to remote server
    deploy_export(export_dir, "user@host:/var/www/tissuumaps/projects/")

Requires: pyvips (already a KINTSUGI dependency)
"""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any


# ---------------------------------------------------------------------------
# Marker → color mapping for common immune markers
# Colors as hex strings for TissUUmaps layerFilters
# ---------------------------------------------------------------------------

_MARKER_COLORS: dict[str, str] = {
    "dapi": "#0000ff",
    "cd3": "#00ff00",
    "cd3e": "#00ff00",
    "cd8": "#ff0000",
    "cd4": "#00ffff",
    "cd45": "#ffff00",
    "cd20": "#ff00ff",
    "ki67": "#ff8000",
    "cd68": "#80ff00",
    "panck": "#ffffff",
    "cd31": "#ff80c0",
    "podoplan": "#00c0c0",
    "actin": "#c0c0c0",
    "cd21": "#c000ff",
    "cd11c": "#ff4040",
    "cd163": "#40ff40",
    "foxp3": "#ffff80",
    "pd1": "#ff80ff",
    "pdl1": "#80ffff",
    "collagen": "#c0a060",
    "vimentin": "#a0a0ff",
    "sma": "#ffa0a0",
}

# Fallback color cycle for unknown markers
_FALLBACK_COLORS = [
    "#00ff00", "#ff0000", "#ffff00", "#00ffff", "#ff00ff",
    "#ff8000", "#80ff00", "#0080ff", "#ff0080", "#80ff80",
    "#ff8080", "#8080ff", "#c0ff00", "#00ffc0", "#c000ff",
]


def _marker_color(name: str, index: int) -> str:
    """Get color hex for a marker name, falling back to cycle colors."""
    key = re.sub(r"[-_]\d+$", "", name.lower()).strip()
    key = re.sub(r"_cyc\d+$", "", key)
    if key in _MARKER_COLORS:
        return _MARKER_COLORS[key]
    return _FALLBACK_COLORS[index % len(_FALLBACK_COLORS)]


def _is_blank(name: str) -> bool:
    """Check if a marker name indicates a blank/autofluorescence channel."""
    lower = name.lower()
    return lower.startswith("blank") or lower.startswith("empty")


def _is_dapi(name: str) -> bool:
    """Check if a marker name is a DAPI channel."""
    return name.lower().startswith("dapi")


# ---------------------------------------------------------------------------
# Export manifest for skip-existing logic
# ---------------------------------------------------------------------------


@dataclass
class _ManifestEntry:
    """Record for a single exported file."""

    source: str
    mtime: float
    size: int
    dzi_path: str
    converted: str  # ISO timestamp


@dataclass
class ExportManifest:
    """Tracks exported files for skip-existing logic.

    Uses source file mtime + size (not MD5) for fast change detection.
    """

    created: str = field(default_factory=lambda: datetime.now().isoformat())
    updated: str = field(default_factory=lambda: datetime.now().isoformat())
    entries: dict[str, _ManifestEntry] = field(default_factory=dict)

    def save(self, path: Path) -> None:
        self.updated = datetime.now().isoformat()
        data = {
            "created": self.created,
            "updated": self.updated,
            "entries": {
                k: {
                    "source": v.source,
                    "mtime": v.mtime,
                    "size": v.size,
                    "dzi_path": v.dzi_path,
                    "converted": v.converted,
                }
                for k, v in self.entries.items()
            },
        }
        path.write_text(json.dumps(data, indent=2))

    @classmethod
    def load(cls, path: Path) -> ExportManifest:
        if not path.exists():
            return cls()
        data = json.loads(path.read_text())
        manifest = cls(
            created=data.get("created", ""),
            updated=data.get("updated", ""),
        )
        for key, entry_data in data.get("entries", {}).items():
            manifest.entries[key] = _ManifestEntry(**entry_data)
        return manifest

    def needs_update(self, source_path: Path, key: str) -> bool:
        """Check if source file has changed since last export."""
        if key not in self.entries:
            return True
        entry = self.entries[key]
        stat = source_path.stat()
        return stat.st_mtime != entry.mtime or stat.st_size != entry.size


# ---------------------------------------------------------------------------
# Discovery functions
# ---------------------------------------------------------------------------


def discover_registered_images(
    registered_dir: Path,
) -> dict[str, list[tuple[str, Path]]]:
    """Find registered marker images organized by cycle.

    Returns ``{cycle_name: [(marker_name, path), ...]}`` sorted by cycle.
    Skips ``registration_data/`` and sentinel files.
    """
    results: dict[str, list[tuple[str, Path]]] = {}
    if not registered_dir.exists():
        return results

    cycle_pattern = re.compile(r"^cyc\d+$", re.IGNORECASE)

    for cycle_dir in sorted(registered_dir.iterdir()):
        if not cycle_dir.is_dir():
            continue
        if not cycle_pattern.match(cycle_dir.name):
            continue

        markers = []
        for tif in sorted(cycle_dir.glob("*.tif")):
            markers.append((tif.stem, tif))
        # Also check .tiff extension
        for tif in sorted(cycle_dir.glob("*.tiff")):
            if tif.stem not in {m[0] for m in markers}:
                markers.append((tif.stem, tif))

        if markers:
            results[cycle_dir.name] = markers

    return results


def discover_segmentation_masks(
    segmented_dir: Path,
) -> list[tuple[str, Path]]:
    """Find segmentation mask TIFFs.

    Looks for ``*label*.tif``, ``*mask*.tif``, or ``*segmentation*.tif``.
    """
    results: list[tuple[str, Path]] = []
    if not segmented_dir.exists():
        return results

    seen: set[Path] = set()
    for pattern in ["*label*.tif", "*mask*.tif", "*segmentation*.tif",
                     "*label*.tiff", "*mask*.tiff", "*segmentation*.tiff"]:
        for tif in sorted(segmented_dir.rglob(pattern)):
            if tif not in seen:
                seen.add(tif)
                results.append((tif.stem, tif))

    return results


def discover_overlays(project_dir: Path) -> dict[str, Path]:
    """Find CSV overlays and h5ad files for TissUUmaps markers.

    Searches ``segmented/`` and ``analysis/`` directories for:

    - ``cell_stats.csv`` — cell features + coordinates
    - ``*_som.csv``, ``*phenotype*.csv`` — clustering results
    - ``*.h5ad`` — AnnData objects
    """
    results: dict[str, Path] = {}
    processed = project_dir / "data" / "processed"

    for search_dir in [processed / "segmented", processed / "analysis"]:
        if not search_dir.exists():
            continue
        for csv in search_dir.rglob("*.csv"):
            results[csv.stem] = csv
        for h5ad in search_dir.rglob("*.h5ad"):
            results[h5ad.stem] = h5ad

    return results


# ---------------------------------------------------------------------------
# DZI conversion
# ---------------------------------------------------------------------------


def convert_to_dzi(
    source_path: Path,
    output_dir: Path,
    name: str,
    *,
    tile_size: int = 256,
    overlap: int = 1,
    suffix: str = ".png",
    is_label: bool = False,
) -> Path:
    """Convert a TIFF to DZI tile pyramid using pyvips.

    Args:
        source_path: Path to source TIFF file.
        output_dir: Directory to write DZI files into.
        name: Base name for the ``.dzi`` file (without extension).
        tile_size: Tile size in pixels.
        overlap: Tile overlap in pixels.
        suffix: Tile image format suffix.
        is_label: Use nearest-neighbor resampling (preserves label IDs).

    Returns:
        Path to the generated ``.dzi`` descriptor file.
    """
    import pyvips

    img = pyvips.Image.new_from_file(str(source_path), access="sequential")

    kwargs: dict[str, Any] = {
        "tile_size": tile_size,
        "overlap": overlap,
        "suffix": suffix,
    }
    if is_label:
        kwargs["region_shrink"] = "nearest"

    output_dir.mkdir(parents=True, exist_ok=True)
    img.dzsave(str(output_dir / name), **kwargs)

    return output_dir / f"{name}.dzi"


# ---------------------------------------------------------------------------
# .tmap project file generation
# ---------------------------------------------------------------------------


def generate_tmap(
    project_name: str,
    export_dir: Path,
    registered_images: dict[str, list[tuple[str, Path]]],
    segmentation_masks: list[tuple[str, Path]],
    overlays: dict[str, Path],
) -> Path:
    """Generate a ``.tmap`` project file for TissUUmaps.

    Creates a JSON project file that references DZI tile layers with
    per-channel color filters and CSV/h5ad marker overlays.

    Returns:
        Path to the generated ``.tmap`` file.
    """
    layers: list[dict[str, Any]] = []
    layer_filters: list[list[dict[str, str]]] = []
    color_index = 0

    # --- Registered marker layers ---
    for cycle_name, markers in registered_images.items():
        for marker_name, _ in markers:
            tile_source = f"images/registered/{cycle_name}/{marker_name}.dzi"

            # First DAPI visible by default, everything else hidden
            visible = _is_dapi(marker_name) and cycle_name == "cyc01"

            layers.append({
                "tileSource": tile_source,
                "name": f"{marker_name} ({cycle_name})",
                "visible": visible,
            })

            color = _marker_color(marker_name, color_index)
            layer_filters.append([{"name": "Color", "value": color}])

            if not _is_blank(marker_name) and not _is_dapi(marker_name):
                color_index += 1

    # --- Segmentation mask layers ---
    for mask_name, _ in segmentation_masks:
        tile_source = f"images/segmented/{mask_name}.dzi"
        layers.append({
            "tileSource": tile_source,
            "name": mask_name,
            "visible": False,
        })
        layer_filters.append([{"name": "Color", "value": "#ffffff"}])

    # --- Marker files (CSV overlays) ---
    marker_files: list[dict[str, Any]] = []
    for overlay_name, overlay_path in overlays.items():
        if overlay_path.suffix == ".csv":
            marker_files.append({
                "path": f"data/{overlay_path.name}",
                "title": overlay_name.replace("_", " ").title(),
                "expectedHeader": {
                    "X": "centroid_1",
                    "Y": "centroid_0",
                },
            })

    # --- Regions (GeoJSON) ---
    regions: dict[str, str] = {}
    regions_dir = export_dir / "regions"
    if regions_dir.exists():
        for geojson in sorted(regions_dir.glob("*.geojson")):
            regions[geojson.stem] = f"regions/{geojson.name}"

    # --- Assemble .tmap ---
    tmap: dict[str, Any] = {
        "filename": f"{project_name}.tmap",
        "layers": layers,
        "layerFilters": layer_filters,
        "compositeMode": "lighter",
        "markerFiles": marker_files,
        "plugins": ["Feature_Space"] if marker_files else [],
        "regions": regions,
    }

    tmap_path = export_dir / f"{project_name}.tmap"
    tmap_path.write_text(json.dumps(tmap, indent=2))
    return tmap_path


# ---------------------------------------------------------------------------
# Overlay preparation
# ---------------------------------------------------------------------------


def _copy_csv_overlay(source: Path, data_dir: Path) -> None:
    """Copy a CSV overlay file to the export data/ directory."""
    data_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(data_dir / source.name))


def _extract_h5ad_phenotypes(source: Path, data_dir: Path) -> Path | None:
    """Copy h5ad and attempt to extract a phenotypes CSV.

    Returns the phenotypes CSV path if extraction succeeded, else None.
    """
    data_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(str(source), str(data_dir / source.name))

    try:
        import anndata

        adata = anndata.read_h5ad(source)
        obs = adata.obs

        # Find spatial coordinate columns
        x_col = next(
            (c for c in ["centroid_1", "X_centroid", "x"] if c in obs.columns), None
        )
        y_col = next(
            (c for c in ["centroid_0", "Y_centroid", "y"] if c in obs.columns), None
        )
        if not x_col or not y_col:
            return None

        # Collect phenotype columns (clusters, leiden, etc.)
        export_cols = [x_col, y_col]
        for col in obs.columns:
            if (
                col.startswith("leiden")
                or col.startswith("merged")
                or "cluster" in col.lower()
                or "phenotype" in col.lower()
            ):
                export_cols.append(col)

        if len(export_cols) <= 2:
            return None

        phenotype_path = data_dir / "phenotypes.csv"
        obs[export_cols].to_csv(phenotype_path)
        return phenotype_path
    except ImportError:
        return None
    except Exception:
        return None


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------


def prepare_export(
    project_dir: Path,
    *,
    force: bool = False,
    dry_run: bool = False,
    exclude_blanks: bool = True,
    progress_callback: Any | None = None,
) -> dict[str, Any]:
    """Prepare TissUUmaps export from a KINTSUGI project.

    Converts registered images and segmentation masks to DZI tile pyramids,
    copies overlay CSVs and h5ad files, and generates a ``.tmap`` project file.

    Args:
        project_dir: Path to KINTSUGI project root.
        force: Reconvert all files regardless of manifest.
        dry_run: Preview what would be exported without writing.
        exclude_blanks: Skip blank/autofluorescence channels (default True).
        progress_callback: Optional ``callable(current, total, message)`` for
            progress reporting.

    Returns:
        Dict with keys: ``status``, ``export_dir``, ``images_converted``, etc.
    """
    project_dir = Path(project_dir).resolve()
    registered_dir = project_dir / "data" / "processed" / "registered"
    segmented_dir = project_dir / "data" / "processed" / "segmented"
    export_dir = project_dir / "data" / "exported"

    # Validate
    if not registered_dir.exists():
        return {
            "status": "error",
            "message": f"No registered images found: {registered_dir}",
        }

    # Discover
    registered = discover_registered_images(registered_dir)
    if not registered:
        return {
            "status": "error",
            "message": "No cycle directories found in registered/",
        }

    masks = discover_segmentation_masks(segmented_dir)
    overlays = discover_overlays(project_dir)

    # Filter blanks
    if exclude_blanks:
        for cycle in registered:
            registered[cycle] = [
                (name, path) for name, path in registered[cycle] if not _is_blank(name)
            ]
        registered = {k: v for k, v in registered.items() if v}

    total_images = sum(len(m) for m in registered.values()) + len(masks)
    project_name = project_dir.name

    # --- Dry run ---
    if dry_run:
        summary: dict[str, Any] = {
            "status": "dry_run",
            "export_dir": str(export_dir),
            "project_name": project_name,
            "registered_cycles": len(registered),
            "registered_images": sum(len(m) for m in registered.values()),
            "segmentation_masks": len(masks),
            "overlays": len(overlays),
            "total_dzi_conversions": total_images,
            "cycles": {},
            "overlay_files": list(overlays.keys()),
        }
        for cycle, markers in registered.items():
            summary["cycles"][cycle] = [name for name, _ in markers]
        if masks:
            summary["masks"] = [name for name, _ in masks]
        return summary

    # --- Real export ---
    manifest_path = export_dir / "export_manifest.json"
    manifest = ExportManifest() if force else ExportManifest.load(manifest_path)

    # Create directory structure
    for subdir in [
        "images/registered",
        "images/segmented",
        "data",
        "regions",
    ]:
        (export_dir / subdir).mkdir(parents=True, exist_ok=True)

    converted = 0
    skipped = 0
    current = 0

    # Convert registered images
    for cycle_name, markers in registered.items():
        cycle_out = export_dir / "images" / "registered" / cycle_name
        for marker_name, source_path in markers:
            current += 1
            key = f"registered/{cycle_name}/{marker_name}"

            if not force and not manifest.needs_update(source_path, key):
                skipped += 1
                if progress_callback:
                    progress_callback(current, total_images, f"Skip {marker_name}")
                continue

            if progress_callback:
                progress_callback(
                    current, total_images, f"{marker_name} ({cycle_name})"
                )

            dzi_path = convert_to_dzi(source_path, cycle_out, marker_name)
            stat = source_path.stat()
            manifest.entries[key] = _ManifestEntry(
                source=str(source_path),
                mtime=stat.st_mtime,
                size=stat.st_size,
                dzi_path=str(dzi_path.relative_to(export_dir)),
                converted=datetime.now().isoformat(),
            )
            converted += 1

    # Convert segmentation masks
    seg_out = export_dir / "images" / "segmented"
    for mask_name, source_path in masks:
        current += 1
        key = f"segmented/{mask_name}"

        if not force and not manifest.needs_update(source_path, key):
            skipped += 1
            if progress_callback:
                progress_callback(current, total_images, f"Skip {mask_name}")
            continue

        if progress_callback:
            progress_callback(current, total_images, f"mask: {mask_name}")

        dzi_path = convert_to_dzi(
            source_path, seg_out, mask_name, is_label=True
        )
        stat = source_path.stat()
        manifest.entries[key] = _ManifestEntry(
            source=str(source_path),
            mtime=stat.st_mtime,
            size=stat.st_size,
            dzi_path=str(dzi_path.relative_to(export_dir)),
            converted=datetime.now().isoformat(),
        )
        converted += 1

    # Copy/prepare overlays
    data_dir = export_dir / "data"
    overlay_count = 0
    extra_overlays: dict[str, Path] = {}

    for overlay_name, overlay_path in overlays.items():
        if overlay_path.suffix == ".csv":
            _copy_csv_overlay(overlay_path, data_dir)
            overlay_count += 1
        elif overlay_path.suffix == ".h5ad":
            phenotype_csv = _extract_h5ad_phenotypes(overlay_path, data_dir)
            overlay_count += 1
            if phenotype_csv:
                extra_overlays[f"{overlay_name}_phenotypes"] = phenotype_csv

    # Merge extracted phenotype CSVs into overlays for .tmap generation
    all_overlays = {**overlays, **extra_overlays}

    # Copy GeoJSON regions
    region_count = 0
    for search_dir in [segmented_dir, project_dir / "data" / "processed" / "analysis"]:
        if not search_dir.exists():
            continue
        for geojson in search_dir.rglob("*.geojson"):
            shutil.copy2(str(geojson), str(export_dir / "regions" / geojson.name))
            region_count += 1

    # Generate .tmap
    tmap_path = generate_tmap(
        project_name, export_dir, registered, masks, all_overlays
    )

    # Save manifest
    manifest.save(manifest_path)

    return {
        "status": "success",
        "export_dir": str(export_dir),
        "tmap_path": str(tmap_path),
        "project_name": project_name,
        "images_converted": converted,
        "images_skipped": skipped,
        "segmentation_masks": len(masks),
        "overlays_copied": overlay_count,
        "regions_copied": region_count,
        "total_images": total_images,
    }


# ---------------------------------------------------------------------------
# Deploy
# ---------------------------------------------------------------------------


def deploy_export(
    export_dir: Path,
    remote_spec: str,
    *,
    dry_run: bool = False,
    bwlimit: int | None = None,
    ssh_key: str | None = None,
) -> int:
    """Deploy exported files to a remote server via rsync.

    Args:
        export_dir: Local export directory (``data/exported/``).
        remote_spec: Remote destination (``user@host:/path``).
        dry_run: Preview what would be transferred without sending.
        bwlimit: Bandwidth limit in KB/s (passed to ``rsync --bwlimit``).
        ssh_key: Path to SSH private key file.

    Returns:
        rsync exit code (0 = success).
    """
    if not export_dir.exists():
        raise FileNotFoundError(f"Export directory not found: {export_dir}")

    cmd = ["rsync", "-avz", "--progress"]
    if dry_run:
        cmd.append("--dry-run")
    if bwlimit:
        cmd.extend(["--bwlimit", str(bwlimit)])
    if ssh_key:
        cmd.extend(["-e", f"ssh -i {ssh_key}"])
    cmd.extend([str(export_dir) + "/", remote_spec])

    result = subprocess.run(cmd, check=False)
    return result.returncode


# ---------------------------------------------------------------------------
# Status
# ---------------------------------------------------------------------------


def get_export_status(project_dir: Path) -> dict[str, Any]:
    """Get export status for a project.

    Returns dict with ``status`` key (``"not_exported"`` or ``"exported"``)
    and summary statistics (DZI count, size, manifest info).
    """
    export_dir = Path(project_dir).resolve() / "data" / "exported"

    if not export_dir.exists():
        return {"status": "not_exported", "message": "No export directory found"}

    manifest_path = export_dir / "export_manifest.json"
    tmap_files = list(export_dir.glob("*.tmap"))
    dzi_files = list(export_dir.rglob("*.dzi"))

    # Total size of exported data
    total_size = sum(f.stat().st_size for f in export_dir.rglob("*") if f.is_file())

    status: dict[str, Any] = {
        "status": "exported",
        "export_dir": str(export_dir),
        "tmap_files": [f.name for f in tmap_files],
        "dzi_count": len(dzi_files),
        "total_size_mb": round(total_size / (1024 * 1024), 1),
    }

    if manifest_path.exists():
        manifest = ExportManifest.load(manifest_path)
        status["manifest_entries"] = len(manifest.entries)
        status["last_updated"] = manifest.updated

    data_dir = export_dir / "data"
    if data_dir.exists():
        status["overlay_files"] = [f.name for f in data_dir.iterdir() if f.is_file()]

    regions_dir = export_dir / "regions"
    if regions_dir.exists():
        status["geojson_regions"] = len(list(regions_dir.glob("*.geojson")))

    return status
