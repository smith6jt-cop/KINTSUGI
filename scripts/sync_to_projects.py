#!/usr/bin/env python3
"""
Sync notebook modules from main KINTSUGI repo to project folders.

This script copies the latest versions of notebook modules (Kdecon, Kstitch,
Kreg, Kview, etc.) and key Python files from the main repository to all
registered project folders.

Uses MD5 checksums for comparison to ensure content changes are detected even
when destination files have newer timestamps (e.g., from being saved with
notebook output during processing runs).

Usage:
    python scripts/sync_to_projects.py           # Sync all projects
    python scripts/sync_to_projects.py --dry-run # Preview changes
    python scripts/sync_to_projects.py --verbose # Show detailed output
    python scripts/sync_to_projects.py --force   # Force sync all files

Project folders are defined in KINTSUGI_PROJECT_DIRS environment variable
or default to known project locations.
"""

import argparse
import glob
import hashlib
import os
import shutil
import sys
from pathlib import Path
from datetime import datetime


# Static project folders (always synced)
_STATIC_PROJECT_FOLDERS = [
    # Test project (mini_project)
    "/blue/maigan/smith6jt/KINTSUGI/test_data/mini_project/notebooks",
    # Full project (1904CC1-1L)
    "/blue/maigan/smith6jt/KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L/notebooks",
]

# Auto-discover batch project directories under KINTSUGI_Projects/
_BATCH_PROJECT_GLOB = "/blue/maigan/smith6jt/KINTSUGI_Projects/*/notebooks"


def _discover_default_project_folders() -> list[str]:
    """Build the default project folder list with auto-discovery of batch projects."""
    folders = list(_STATIC_PROJECT_FOLDERS)
    for path in sorted(glob.glob(_BATCH_PROJECT_GLOB)):
        if path not in folders:
            folders.append(path)
    return folders


DEFAULT_PROJECT_FOLDERS = _discover_default_project_folders()

# Directories to sync (relative to notebooks/)
SYNC_DIRECTORIES = [
    "Kdecon",
    "Kstitch",
    "Kreg",
    "Kview",
    "Kview2",
    "Kseg",
]

# Individual files to sync (relative to notebooks/)
SYNC_FILES = [
    "Kio.py",
    "Kprocess.py",
    "Kutils.py",
    "Kview_qc.py",
    "Kpipeline.py",
    "Kvis.py",
    "1_Single_Channel_Eval.ipynb",
    "2_Cycle_Processing.ipynb",
    "3_Signal_Isolation_QC.ipynb",
    "4_Segmentation_Analysis.ipynb",
    "config_example.json",
    "MIGRATION_GUIDE.md",
]


def compute_file_checksum(filepath: Path) -> str:
    """Compute MD5 checksum of a file."""
    hasher = hashlib.md5()
    with open(filepath, 'rb') as f:
        # Read in chunks for large files
        for chunk in iter(lambda: f.read(65536), b''):
            hasher.update(chunk)
    return hasher.hexdigest()


def files_are_identical(src: Path, dst: Path) -> bool:
    """
    Check if two files have identical content using MD5 checksum.

    This is more reliable than timestamp comparison because:
    - Notebooks saved with output get newer timestamps
    - Git operations can change timestamps
    - Network file systems may have timestamp issues
    """
    if not dst.exists():
        return False

    # Quick check: different sizes means definitely different
    if src.stat().st_size != dst.stat().st_size:
        return False

    # Compare checksums for definitive answer
    return compute_file_checksum(src) == compute_file_checksum(dst)


def get_repo_root() -> Path:
    """Get the KINTSUGI repository root directory."""
    # Try to find repo root from this script's location
    script_path = Path(__file__).resolve()
    repo_root = script_path.parent.parent

    # Verify it's the right repo
    if (repo_root / "src" / "kintsugi").exists():
        return repo_root

    # Fall back to environment variable
    if "KINTSUGI_REPO" in os.environ:
        return Path(os.environ["KINTSUGI_REPO"])

    # Fall back to hardcoded path
    return Path("/blue/maigan/smith6jt/KINTSUGI")


def get_project_folders() -> list[Path]:
    """Get list of project folders to sync to."""
    # Check environment variable first
    if "KINTSUGI_PROJECT_DIRS" in os.environ:
        paths = os.environ["KINTSUGI_PROJECT_DIRS"].split(":")
        return [Path(p) for p in paths if p]

    # Use defaults
    return [Path(p) for p in DEFAULT_PROJECT_FOLDERS]


def sync_directory(
    src: Path,
    dst: Path,
    dry_run: bool = False,
    verbose: bool = False,
    force: bool = False,
) -> tuple[int, int]:
    """
    Sync a directory from source to destination.

    Returns:
        Tuple of (files_copied, files_skipped)
    """
    files_copied = 0
    files_skipped = 0

    if not src.exists():
        if verbose:
            print(f"  SKIP: Source not found: {src}")
        return 0, 0

    # Create destination if needed
    if not dry_run:
        dst.mkdir(parents=True, exist_ok=True)

    # Sync all Python files and __init__.py
    for src_file in src.glob("**/*.py"):
        rel_path = src_file.relative_to(src)
        dst_file = dst / rel_path

        # Check if file needs updating using checksum comparison
        needs_update = force or not files_are_identical(src_file, dst_file)

        if needs_update:
            if dry_run:
                print(f"  Would copy: {rel_path}")
            else:
                dst_file.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(src_file, dst_file)
                if verbose:
                    print(f"  Copied: {rel_path}")
            files_copied += 1
        else:
            if verbose:
                print(f"  SKIP: Up to date: {rel_path}")
            files_skipped += 1

    return files_copied, files_skipped


def sync_file(
    src: Path,
    dst: Path,
    dry_run: bool = False,
    verbose: bool = False,
    force: bool = False,
) -> bool:
    """
    Sync a single file from source to destination.

    Uses checksum comparison to detect content changes, which is more reliable
    than timestamp comparison (notebooks saved with output get newer timestamps).

    Returns:
        True if file was copied, False if skipped
    """
    if not src.exists():
        if verbose:
            print(f"  SKIP: Source not found: {src.name}")
        return False

    # Check if file needs updating using checksum comparison
    if not force and files_are_identical(src, dst):
        if verbose:
            print(f"  SKIP: Up to date: {src.name}")
        return False

    if dry_run:
        print(f"  Would copy: {src.name}")
    else:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        if verbose:
            print(f"  Copied: {src.name}")

    return True


def sync_to_project(
    repo_root: Path,
    project_notebooks: Path,
    dry_run: bool = False,
    verbose: bool = False,
    force: bool = False,
) -> dict:
    """
    Sync all modules from repo to a single project folder.

    Returns:
        Dictionary with sync statistics
    """
    stats = {
        "directories": {"copied": 0, "skipped": 0},
        "files": {"copied": 0, "skipped": 0},
    }

    src_notebooks = repo_root / "notebooks"

    # Sync directories
    for dir_name in SYNC_DIRECTORIES:
        src_dir = src_notebooks / dir_name
        dst_dir = project_notebooks / dir_name

        if src_dir.exists():
            copied, skipped = sync_directory(src_dir, dst_dir, dry_run, verbose, force)
            stats["directories"]["copied"] += copied
            stats["directories"]["skipped"] += skipped

    # Sync individual files
    for file_name in SYNC_FILES:
        src_file = src_notebooks / file_name
        dst_file = project_notebooks / file_name

        if sync_file(src_file, dst_file, dry_run, verbose, force):
            stats["files"]["copied"] += 1
        else:
            stats["files"]["skipped"] += 1

    return stats


def main():
    parser = argparse.ArgumentParser(
        description="Sync KINTSUGI notebook modules to project folders",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument("--dry-run", "-n", action="store_true",
                       help="Preview changes without copying")
    parser.add_argument("--verbose", "-v", action="store_true",
                       help="Show detailed output")
    parser.add_argument("--force", "-f", action="store_true",
                       help="Force sync all files regardless of checksum")
    parser.add_argument("--project", "-p", type=str, action="append",
                       help="Specific project folder to sync (can be repeated)")

    args = parser.parse_args()

    repo_root = get_repo_root()

    if args.project:
        project_folders = [Path(p) for p in args.project]
    else:
        project_folders = get_project_folders()

    print(f"KINTSUGI Sync - {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Source: {repo_root}/notebooks/")
    if args.dry_run:
        print("DRY RUN - no changes will be made")
    if args.force:
        print("FORCE MODE - syncing all files regardless of content")
    print()

    total_copied = 0
    total_skipped = 0

    for project_path in project_folders:
        if not project_path.exists():
            print(f"SKIP: Project folder not found: {project_path}")
            continue

        print(f"Syncing to: {project_path}")
        stats = sync_to_project(repo_root, project_path, args.dry_run, args.verbose, args.force)

        dir_copied = stats["directories"]["copied"]
        dir_skipped = stats["directories"]["skipped"]
        file_copied = stats["files"]["copied"]
        file_skipped = stats["files"]["skipped"]

        copied = dir_copied + file_copied
        skipped = dir_skipped + file_skipped

        if copied > 0:
            print(f"  -> {copied} files {'would be ' if args.dry_run else ''}copied, {skipped} up to date")
        else:
            print(f"  -> All {skipped} files up to date")

        total_copied += copied
        total_skipped += skipped

    print()
    if total_copied > 0:
        action = "would be synced" if args.dry_run else "synced"
        print(f"Total: {total_copied} files {action}, {total_skipped} already up to date")
    else:
        print(f"All {total_skipped} files already up to date across all projects")

    return 0


if __name__ == "__main__":
    sys.exit(main())
