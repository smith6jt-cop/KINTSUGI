"""
Pipeline-aware cleanup safety system.

Encodes the full pipeline dependency graph as a single source of truth,
checks ALL consumers before deletion, and stages data to a recoverable
trash directory instead of permanent deletion.

The core problem this solves: intermediate data (e.g., deconvolved z-stacks)
may be consumed by multiple downstream stages (EDF *and* vessel3d). Deleting
intermediates after one consumer finishes but before another runs causes
irreversible data loss.

Usage:
    from kintsugi.cleanup import assess_cleanup_safety, execute_cleanup

    manifest = assess_cleanup_safety(project_dir)
    for entry in manifest.entries:
        print(f"{entry.directory}: {entry.status} — {entry.reason}")

    execute_cleanup(project_dir, manifest, use_trash=True)
"""

from __future__ import annotations

import json
import shutil
import time
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any

# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class PipelineStage:
    """A processing stage in the KINTSUGI pipeline."""

    name: str  # "stitch", "deconvolve", "edf", "vessel3d", "registration"
    directory: str  # relative path under data/processed/
    sentinel_pattern: str  # e.g. "cyc{cycle}/.snakemake_complete"
    is_optional: bool = False
    config_key: str | None = None  # key in config.yaml optional_stages


@dataclass
class DataDependency:
    """A directory consumed by one or more pipeline stages."""

    directory: str  # relative path under data/processed/ (e.g. "deconvolved")
    consumed_by: list[str]  # stage names that need this data


@dataclass
class CleanupEntry:
    """Assessment result for a single directory."""

    directory: str  # relative path under data/processed/
    abs_path: Path
    status: str  # "safe", "blocked", "missing", "already_deleted"
    reason: str  # human-readable explanation
    consumers: list[str] = field(default_factory=list)
    blocking_consumers: list[str] = field(default_factory=list)
    size_bytes: int = 0

    @property
    def size_mb(self) -> float:
        return self.size_bytes / (1024 * 1024)

    @property
    def size_gb(self) -> float:
        return self.size_bytes / (1024 * 1024 * 1024)


@dataclass
class CleanupManifest:
    """Complete assessment of cleanup safety for a project."""

    project_dir: Path
    entries: list[CleanupEntry] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    timestamp: str = field(default_factory=lambda: datetime.now().isoformat())

    @property
    def safe_entries(self) -> list[CleanupEntry]:
        return [e for e in self.entries if e.status == "safe"]

    @property
    def blocked_entries(self) -> list[CleanupEntry]:
        return [e for e in self.entries if e.status == "blocked"]

    @property
    def total_reclaimable_bytes(self) -> int:
        return sum(e.size_bytes for e in self.safe_entries)

    @property
    def total_reclaimable_gb(self) -> float:
        return self.total_reclaimable_bytes / (1024 * 1024 * 1024)

    def to_dict(self) -> dict[str, Any]:
        return {
            "project_dir": str(self.project_dir),
            "timestamp": self.timestamp,
            "entries": [
                {
                    "directory": e.directory,
                    "status": e.status,
                    "reason": e.reason,
                    "consumers": e.consumers,
                    "blocking_consumers": e.blocking_consumers,
                    "size_bytes": e.size_bytes,
                }
                for e in self.entries
            ],
            "warnings": self.warnings,
        }


@dataclass
class TrashReceipt:
    """Record of a trashed directory for recovery."""

    original_path: str
    trash_path: str
    timestamp: str
    size_bytes: int
    manifest_snapshot: dict[str, Any]


# ---------------------------------------------------------------------------
# Pipeline dependency graph (single source of truth)
# ---------------------------------------------------------------------------

PIPELINE_STAGES: list[PipelineStage] = [
    PipelineStage(
        name="stitch",
        directory="stitched",
        sentinel_pattern="cyc{cycle}/.snakemake_complete",
    ),
    PipelineStage(
        name="deconvolve",
        directory="deconvolved",
        sentinel_pattern="cyc{cycle}/.snakemake_complete",
    ),
    PipelineStage(
        name="edf",
        directory="edf",
        sentinel_pattern="cyc{cycle}/.snakemake_complete",
    ),
    PipelineStage(
        name="vessel3d",
        directory="vessel_3d",
        sentinel_pattern="cyc{cycle}/.snakemake_complete",
        is_optional=True,
        config_key="vessel3d",
    ),
    PipelineStage(
        name="registration",
        directory="registered",
        sentinel_pattern=".snakemake_complete",
    ),
]

# Note: raw/ is at data/raw/ (not data/processed/raw/) so it's handled
# separately in execute_cleanup() rather than through this list.
DATA_DEPENDENCIES: list[DataDependency] = [
    DataDependency(directory="stitched", consumed_by=["deconvolve"]),
    DataDependency(directory="deconvolved", consumed_by=["edf", "vessel3d"]),
]

# QC sentinel files that must exist before cleanup
QC_SENTINELS = ["stitch", "decon", "edf"]


def _get_stage(name: str) -> PipelineStage | None:
    for s in PIPELINE_STAGES:
        if s.name == name:
            return s
    return None


# ---------------------------------------------------------------------------
# Config helpers
# ---------------------------------------------------------------------------


def _load_workflow_config(project_dir: Path) -> dict[str, Any]:
    """Load workflow/config.yaml if it exists."""
    config_path = project_dir / "workflow" / "config.yaml"
    if not config_path.exists():
        return {}
    try:
        import yaml

        with open(config_path) as f:
            return yaml.safe_load(f) or {}
    except Exception:
        return {}


def _get_optional_stages(config: dict[str, Any]) -> dict[str, dict[str, Any]]:
    """Extract optional_stages section from config."""
    return config.get("optional_stages", {})


def _get_cycles(config: dict[str, Any], project_dir: Path) -> list[int]:
    """Get cycle numbers from config or discover from EDF outputs."""
    if config.get("cycles"):
        return [int(c) for c in config["cycles"]]

    # Discover from EDF directory
    edf_dir = project_dir / "data" / "processed" / "edf"
    if edf_dir.exists():
        cycles = []
        for d in sorted(edf_dir.iterdir()):
            if d.is_dir() and d.name.startswith("cyc"):
                try:
                    num = int(d.name[3:])
                    cycles.append(num)
                except ValueError:
                    pass
        if cycles:
            return cycles

    return []


# ---------------------------------------------------------------------------
# Sentinel checking
# ---------------------------------------------------------------------------


def _check_stage_sentinels(
    project_dir: Path,
    stage: PipelineStage,
    cycles: list[int],
    required_cycles: list[int] | None = None,
) -> tuple[bool, list[str]]:
    """Check if all sentinel files exist for a stage.

    Returns (all_complete, list_of_missing_sentinels).
    """
    proc_dir = project_dir / "data" / "processed"
    check_cycles = required_cycles if required_cycles is not None else cycles
    missing = []

    if "{cycle}" in stage.sentinel_pattern:
        for cyc in check_cycles:
            sentinel = proc_dir / stage.directory / stage.sentinel_pattern.format(
                cycle=f"{cyc:02d}"
            )
            if not sentinel.exists():
                missing.append(str(sentinel.relative_to(project_dir)))
    else:
        sentinel = proc_dir / stage.directory / stage.sentinel_pattern
        if not sentinel.exists():
            missing.append(str(sentinel.relative_to(project_dir)))

    return len(missing) == 0, missing


def _check_qc_sentinels(project_dir: Path) -> tuple[bool, list[str]]:
    """Check QC sentinel files (stitch, decon, edf).

    Registration QC is NOT required for intermediate cleanup.
    """
    qc_dir = project_dir / "qc_plots"
    missing = []
    for stage in QC_SENTINELS:
        sentinel = qc_dir / f".snakemake_complete_{stage}"
        if not sentinel.exists():
            missing.append(stage)
    return len(missing) == 0, missing


# ---------------------------------------------------------------------------
# Directory size calculation
# ---------------------------------------------------------------------------


def _dir_size(path: Path) -> int:
    """Calculate total size of a directory in bytes."""
    if not path.exists():
        return 0
    total = 0
    try:
        for f in path.rglob("*"):
            if f.is_file():
                total += f.stat().st_size
    except (PermissionError, OSError):
        pass
    return total


# ---------------------------------------------------------------------------
# Core assessment function
# ---------------------------------------------------------------------------


def assess_cleanup_safety(
    project_dir: Path | str,
    config: dict[str, Any] | None = None,
    skip_size: bool = False,
) -> CleanupManifest:
    """Assess which intermediate directories can be safely deleted.

    For each deletable directory, checks ALL consumers (not just the
    obvious one). Returns a manifest with per-directory safe/blocked
    status and human-readable reasons.

    Args:
        project_dir: Path to the KINTSUGI project.
        config: Pre-loaded workflow config (loaded from config.yaml if None).
        skip_size: Skip directory size calculation (faster for dry runs).

    Returns:
        CleanupManifest with per-directory assessment.
    """
    project_dir = Path(project_dir).resolve()
    proc_dir = project_dir / "data" / "processed"

    if config is None:
        config = _load_workflow_config(project_dir)

    optional_stages = _get_optional_stages(config)
    cycles = _get_cycles(config, project_dir)
    manifest = CleanupManifest(project_dir=project_dir)

    if not cycles:
        manifest.warnings.append(
            "Could not determine cycle list. Check workflow/config.yaml."
        )
        return manifest

    # Check QC sentinels first
    qc_ok, qc_missing = _check_qc_sentinels(project_dir)
    if not qc_ok:
        manifest.warnings.append(
            f"QC sentinels missing: {', '.join(qc_missing)}. "
            "Run QC rules before cleanup."
        )

    # Assess each deletable directory
    for dep in DATA_DEPENDENCIES:
        dir_path = proc_dir / dep.directory
        entry = CleanupEntry(
            directory=dep.directory,
            abs_path=dir_path,
            status="safe",
            reason="",
            consumers=list(dep.consumed_by),
        )

        # Check if directory exists
        if not dir_path.exists():
            entry.status = "already_deleted"
            entry.reason = f"{dep.directory}/ does not exist"
            manifest.entries.append(entry)
            continue

        # Calculate size
        if not skip_size:
            entry.size_bytes = _dir_size(dir_path)

        # Check QC requirement
        if not qc_ok:
            entry.status = "blocked"
            entry.reason = f"QC incomplete (missing: {', '.join(qc_missing)})"
            entry.blocking_consumers = ["qc"]
            manifest.entries.append(entry)
            continue

        # Check each consumer
        blocking = []
        reasons = []

        for consumer_name in dep.consumed_by:
            stage = _get_stage(consumer_name)
            if stage is None:
                continue

            if stage.is_optional:
                # Check optional stage config
                stage_config = optional_stages.get(stage.config_key or stage.name, {})

                if not optional_stages:
                    # No optional_stages section at all — conservative default
                    blocking.append(consumer_name)
                    reasons.append(
                        f"{consumer_name}: optional_stages not declared in config.yaml. "
                        f"Set optional_stages.{stage.config_key or stage.name}.enabled: false "
                        "to allow cleanup, or re-run 'kintsugi workflow config .'"
                    )
                    continue

                if not isinstance(stage_config, dict):
                    # Key exists but not a dict — treat as not configured
                    continue

                enabled = stage_config.get("enabled", False)
                if not enabled:
                    # Explicitly disabled — not blocking
                    continue

                # Stage is enabled — check sentinels
                declared_cycles = stage_config.get("cycles", [])
                if declared_cycles:
                    # Only specific cycles need this stage
                    required = [int(c) for c in declared_cycles]
                else:
                    # Empty list or omitted = all cycles
                    required = cycles

                complete, missing = _check_stage_sentinels(
                    project_dir, stage, cycles, required_cycles=required
                )
                if not complete:
                    blocking.append(consumer_name)
                    n_missing = len(missing)
                    reasons.append(
                        f"{consumer_name}: {n_missing} sentinel(s) missing "
                        f"(enabled for {len(required)} cycles)"
                    )
            else:
                # Required stage — just check sentinels
                complete, missing = _check_stage_sentinels(
                    project_dir, stage, cycles
                )
                if not complete:
                    blocking.append(consumer_name)
                    n_missing = len(missing)
                    reasons.append(f"{consumer_name}: {n_missing} sentinel(s) missing")

        if blocking:
            entry.status = "blocked"
            entry.blocking_consumers = blocking
            entry.reason = "; ".join(reasons)
        else:
            entry.status = "safe"
            consumer_names = ", ".join(dep.consumed_by)
            entry.reason = f"All consumers complete ({consumer_names})"

        manifest.entries.append(entry)

    return manifest


# ---------------------------------------------------------------------------
# Staged deletion (trash directory)
# ---------------------------------------------------------------------------

TRASH_DIR_NAME = ".trash"
TRASH_RECEIPT_FILE = "receipt.json"
TRASH_MAX_AGE_DAYS = 7


def _trash_dir(project_dir: Path) -> Path:
    """Get the trash directory for a project."""
    return project_dir / "data" / TRASH_DIR_NAME


def _move_to_trash(
    dir_path: Path,
    project_dir: Path,
    manifest_snapshot: dict[str, Any] | None = None,
) -> TrashReceipt:
    """Move a directory to the project's trash directory.

    Uses shutil.move() which is atomic on the same filesystem (POSIX rename).
    """
    trash = _trash_dir(project_dir)
    trash.mkdir(parents=True, exist_ok=True)

    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    dir_name = dir_path.name
    trash_dest = trash / f"{dir_name}_{timestamp}"

    # Ensure unique name
    counter = 0
    while trash_dest.exists():
        counter += 1
        trash_dest = trash / f"{dir_name}_{timestamp}_{counter}"

    size = _dir_size(dir_path)
    shutil.move(str(dir_path), str(trash_dest))

    receipt = TrashReceipt(
        original_path=str(dir_path),
        trash_path=str(trash_dest),
        timestamp=timestamp,
        size_bytes=size,
        manifest_snapshot=manifest_snapshot or {},
    )

    # Write receipt
    receipt_file = trash_dest / TRASH_RECEIPT_FILE
    receipt_file.write_text(
        json.dumps(
            {
                "original_path": receipt.original_path,
                "trash_path": receipt.trash_path,
                "timestamp": receipt.timestamp,
                "size_bytes": receipt.size_bytes,
                "manifest_snapshot": receipt.manifest_snapshot,
            },
            indent=2,
        )
    )

    return receipt


def _permanent_delete(dir_path: Path) -> None:
    """Permanently delete a directory."""
    if dir_path.exists():
        shutil.rmtree(str(dir_path))


# ---------------------------------------------------------------------------
# Execute cleanup
# ---------------------------------------------------------------------------


def execute_cleanup(
    project_dir: Path | str,
    manifest: CleanupManifest | None = None,
    use_trash: bool = True,
    skip_vessel3d: bool = False,
    config: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Execute cleanup based on manifest assessment.

    Args:
        project_dir: Path to the KINTSUGI project.
        manifest: Pre-computed manifest. If None, assess_cleanup_safety() is called.
        use_trash: Move to trash instead of permanent delete.
        skip_vessel3d: Bypass vessel3d check (escape hatch).
        config: Pre-loaded workflow config.

    Returns:
        Summary dict with deleted/skipped/errors counts.
    """
    project_dir = Path(project_dir).resolve()

    if manifest is None:
        manifest = assess_cleanup_safety(project_dir, config=config)

    result = {
        "deleted": [],
        "skipped": [],
        "errors": [],
        "trash_dir": str(_trash_dir(project_dir)) if use_trash else None,
    }

    # Pre-check raw/ eligibility BEFORE deleting anything (stitched/ contains
    # the stitch sentinels — deleting it first would make this check fail).
    raw_dir = project_dir / "data" / "raw"
    raw_eligible = False
    if raw_dir.exists():
        stitch_stage = _get_stage("stitch")
        cycles_list = _get_cycles(config or {}, project_dir)
        if stitch_stage and cycles_list:
            raw_eligible, _ = _check_stage_sentinels(
                project_dir, stitch_stage, cycles_list
            )
            # Also require QC to be complete
            qc_ok, _ = _check_qc_sentinels(project_dir)
            raw_eligible = raw_eligible and qc_ok

    for entry in manifest.entries:
        if entry.status == "already_deleted":
            continue

        if entry.status == "blocked":
            # Check if we should override vessel3d block
            if skip_vessel3d and entry.blocking_consumers == ["vessel3d"]:
                pass  # Allow deletion
            elif skip_vessel3d and "vessel3d" in entry.blocking_consumers:
                # vessel3d is one of multiple blockers — check if others resolved
                other_blockers = [
                    b for b in entry.blocking_consumers if b != "vessel3d"
                ]
                if other_blockers:
                    result["skipped"].append(
                        {
                            "directory": entry.directory,
                            "reason": entry.reason,
                        }
                    )
                    continue
            else:
                result["skipped"].append(
                    {
                        "directory": entry.directory,
                        "reason": entry.reason,
                    }
                )
                continue

        # Safe to delete
        try:
            if use_trash:
                receipt = _move_to_trash(
                    entry.abs_path,
                    project_dir,
                    manifest_snapshot=manifest.to_dict(),
                )
                result["deleted"].append(
                    {
                        "directory": entry.directory,
                        "trash_path": receipt.trash_path,
                        "size_bytes": entry.size_bytes,
                    }
                )
            else:
                _permanent_delete(entry.abs_path)
                result["deleted"].append(
                    {
                        "directory": entry.directory,
                        "size_bytes": entry.size_bytes,
                    }
                )
        except Exception as e:
            result["errors"].append(
                {
                    "directory": entry.directory,
                    "error": str(e),
                }
            )

    # Delete raw/ if pre-check determined it's eligible
    if raw_eligible and raw_dir.exists():
        try:
            if use_trash:
                receipt = _move_to_trash(
                    raw_dir,
                    project_dir,
                    manifest_snapshot=manifest.to_dict(),
                )
                result["deleted"].append(
                    {
                        "directory": "raw",
                        "trash_path": receipt.trash_path,
                        "size_bytes": _dir_size(raw_dir)
                        if raw_dir.exists()
                        else 0,
                    }
                )
            else:
                size = _dir_size(raw_dir)
                _permanent_delete(raw_dir)
                # Recreate empty dir for project structure
                raw_dir.mkdir(parents=True, exist_ok=True)
                result["deleted"].append(
                    {
                        "directory": "raw",
                        "size_bytes": size,
                    }
                )
        except Exception as e:
            result["errors"].append(
                {"directory": "raw", "error": str(e)}
            )

    return result


# ---------------------------------------------------------------------------
# Trash management (recover, purge, list)
# ---------------------------------------------------------------------------


def list_trash(project_dir: Path | str) -> list[dict[str, Any]]:
    """List all entries in the project's trash directory."""
    project_dir = Path(project_dir).resolve()
    trash = _trash_dir(project_dir)
    entries = []

    if not trash.exists():
        return entries

    for item in sorted(trash.iterdir()):
        if not item.is_dir():
            continue

        receipt_file = item / TRASH_RECEIPT_FILE
        if receipt_file.exists():
            try:
                receipt = json.loads(receipt_file.read_text())
                receipt["current_path"] = str(item)
                entries.append(receipt)
            except (json.JSONDecodeError, KeyError):
                entries.append(
                    {
                        "current_path": str(item),
                        "original_path": "unknown",
                        "timestamp": "unknown",
                        "size_bytes": _dir_size(item),
                    }
                )
        else:
            entries.append(
                {
                    "current_path": str(item),
                    "original_path": "unknown",
                    "timestamp": "unknown",
                    "size_bytes": _dir_size(item),
                }
            )

    return entries


def recover_trash(
    project_dir: Path | str,
    trash_entry_name: str | None = None,
) -> dict[str, Any]:
    """Recover trashed data back to its original location.

    Args:
        project_dir: Path to the KINTSUGI project.
        trash_entry_name: Name of the trash entry directory to recover.
            If None, lists available entries.

    Returns:
        Result dict with recovered path or list of available entries.
    """
    project_dir = Path(project_dir).resolve()
    trash = _trash_dir(project_dir)

    if trash_entry_name is None:
        return {"entries": list_trash(project_dir)}

    entry_path = trash / trash_entry_name
    if not entry_path.exists():
        return {"error": f"Trash entry not found: {trash_entry_name}"}

    receipt_file = entry_path / TRASH_RECEIPT_FILE
    if not receipt_file.exists():
        return {"error": f"No receipt found in {trash_entry_name}"}

    try:
        receipt = json.loads(receipt_file.read_text())
    except (json.JSONDecodeError, KeyError):
        return {"error": f"Corrupt receipt in {trash_entry_name}"}

    original = Path(receipt["original_path"])
    if original.exists():
        return {
            "error": f"Original path already exists: {original}. "
            "Remove it first or recover manually."
        }

    # Remove receipt before moving (it's not part of original data)
    receipt_file.unlink()

    # Move back
    shutil.move(str(entry_path), str(original))

    return {
        "recovered": str(original),
        "from": str(entry_path),
        "size_bytes": receipt.get("size_bytes", 0),
    }


def purge_trash(
    project_dir: Path | str,
    max_age_days: int = TRASH_MAX_AGE_DAYS,
    purge_all: bool = False,
) -> dict[str, Any]:
    """Permanently delete old trash entries.

    Args:
        project_dir: Path to the KINTSUGI project.
        max_age_days: Delete entries older than this many days.
        purge_all: Delete all entries regardless of age.

    Returns:
        Summary dict with purged entries.
    """
    project_dir = Path(project_dir).resolve()
    trash = _trash_dir(project_dir)
    purged = []

    if not trash.exists():
        return {"purged": purged, "total_bytes": 0}

    now = time.time()
    cutoff_seconds = max_age_days * 86400

    for item in sorted(trash.iterdir()):
        if not item.is_dir():
            continue

        should_purge = purge_all

        if not should_purge:
            # Check age from receipt timestamp or directory mtime
            receipt_file = item / TRASH_RECEIPT_FILE
            if receipt_file.exists():
                try:
                    receipt = json.loads(receipt_file.read_text())
                    ts = receipt.get("timestamp", "")
                    # Parse YYYYMMDD_HHMMSS
                    if ts and len(ts) >= 8:
                        entry_time = datetime.strptime(
                            ts[:15], "%Y%m%d_%H%M%S"
                        ).timestamp()
                        if (now - entry_time) > cutoff_seconds:
                            should_purge = True
                except (json.JSONDecodeError, ValueError):
                    pass

            if not should_purge:
                # Fall back to directory mtime
                try:
                    mtime = item.stat().st_mtime
                    if (now - mtime) > cutoff_seconds:
                        should_purge = True
                except OSError:
                    pass

        if should_purge:
            size = _dir_size(item)
            shutil.rmtree(str(item))
            purged.append({"path": str(item), "size_bytes": size})

    # Remove trash dir if empty
    if trash.exists() and not any(trash.iterdir()):
        trash.rmdir()

    total_bytes = sum(e["size_bytes"] for e in purged)
    return {"purged": purged, "total_bytes": total_bytes}
