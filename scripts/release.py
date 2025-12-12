#!/usr/bin/env python3
"""
KINTSUGI Release Script

Automates version bumping and release preparation using commitizen.

Usage:
    python scripts/release.py --dry-run      # Preview changes
    python scripts/release.py --bump patch   # Bump patch version
    python scripts/release.py --bump minor   # Bump minor version
    python scripts/release.py --bump major   # Bump major version
    python scripts/release.py --auto         # Auto-detect from commits
"""

import argparse
import re
import subprocess
import sys
from pathlib import Path


def run_command(cmd: list[str], check: bool = True, capture: bool = True) -> subprocess.CompletedProcess:
    """Run a shell command and return the result."""
    print(f"  Running: {' '.join(cmd)}")
    return subprocess.run(cmd, capture_output=capture, text=True, check=check)


def get_current_version() -> str:
    """Get the current version from pyproject.toml."""
    pyproject = Path("pyproject.toml")
    content = pyproject.read_text()
    match = re.search(r'^version\s*=\s*"([^"]+)"', content, re.MULTILINE)
    if match:
        return match.group(1)
    raise ValueError("Could not find version in pyproject.toml")


def get_commits_since_tag(tag: str | None = None) -> list[str]:
    """Get list of commits since the last tag."""
    if tag is None:
        # Get the most recent tag
        result = run_command(["git", "describe", "--tags", "--abbrev=0"], check=False)
        if result.returncode != 0:
            # No tags exist, get all commits
            result = run_command(["git", "log", "--oneline"])
            return result.stdout.strip().split("\n") if result.stdout else []
        tag = result.stdout.strip()

    result = run_command(["git", "log", f"{tag}..HEAD", "--oneline"])
    return result.stdout.strip().split("\n") if result.stdout else []


def analyze_commits(commits: list[str]) -> dict:
    """Analyze commits to determine bump type."""
    analysis = {
        "breaking": [],
        "feat": [],
        "fix": [],
        "refactor": [],
        "perf": [],
        "docs": [],
        "style": [],
        "test": [],
        "build": [],
        "ci": [],
        "chore": [],
        "other": [],
    }

    for commit in commits:
        if not commit:
            continue

        # Check for breaking changes (! or BREAKING CHANGE)
        if "!" in commit.split()[1] if len(commit.split()) > 1 else False:
            analysis["breaking"].append(commit)
        elif commit.lower().find("breaking change") >= 0:
            analysis["breaking"].append(commit)

        # Categorize by type
        for type_name in ["feat", "fix", "refactor", "perf", "docs", "style", "test", "build", "ci", "chore"]:
            if re.match(rf"^[a-f0-9]+\s+{type_name}(\(.+\))?:", commit):
                analysis[type_name].append(commit)
                break
        else:
            analysis["other"].append(commit)

    return analysis


def determine_bump_type(analysis: dict) -> str:
    """Determine the version bump type based on commit analysis."""
    if analysis["breaking"]:
        return "major"
    elif analysis["feat"]:
        return "minor"
    elif any(analysis[t] for t in ["fix", "refactor", "perf", "docs", "style", "test", "build", "ci", "chore"]):
        return "patch"
    else:
        return "none"


def bump_version(bump_type: str, dry_run: bool = False) -> str | None:
    """Bump the version using commitizen."""
    if bump_type == "none":
        print("No commits require a version bump.")
        return None

    cmd = ["cz", "bump", "--increment", bump_type.upper()]
    if dry_run:
        cmd.append("--dry-run")

    result = run_command(cmd, check=False, capture=not dry_run)

    if result.returncode != 0:
        if "No commits found" in (result.stderr or ""):
            print("No commits found since last tag.")
            return None
        print(f"Error: {result.stderr}")
        return None

    # Get new version
    return get_current_version()


def generate_changelog(dry_run: bool = False) -> None:
    """Generate changelog using commitizen."""
    cmd = ["cz", "changelog"]
    if dry_run:
        cmd.append("--dry-run")

    run_command(cmd, check=False)


def main():
    parser = argparse.ArgumentParser(description="KINTSUGI Release Script")
    parser.add_argument("--dry-run", action="store_true", help="Preview changes without making them")
    parser.add_argument("--bump", choices=["patch", "minor", "major"], help="Manually specify bump type")
    parser.add_argument("--auto", action="store_true", help="Auto-detect bump type from commits")
    parser.add_argument("--changelog-only", action="store_true", help="Only update changelog")

    args = parser.parse_args()

    print("=" * 60)
    print("KINTSUGI Release Script")
    print("=" * 60)

    # Get current version
    current_version = get_current_version()
    print(f"\nCurrent version: {current_version}")

    # Get commits since last tag
    commits = get_commits_since_tag()
    print(f"Commits since last release: {len(commits)}")

    if not commits:
        print("\nNo new commits to release.")
        return 0

    # Analyze commits
    analysis = analyze_commits(commits)
    print("\nCommit Analysis:")
    for type_name, type_commits in analysis.items():
        if type_commits:
            print(f"  {type_name}: {len(type_commits)}")

    # Determine or use specified bump type
    if args.bump:
        bump_type = args.bump
        print(f"\nManual bump type: {bump_type}")
    else:
        bump_type = determine_bump_type(analysis)
        print(f"\nAuto-detected bump type: {bump_type}")

    if args.changelog_only:
        print("\nGenerating changelog only...")
        generate_changelog(dry_run=args.dry_run)
        return 0

    if args.auto or args.bump:
        print(f"\n{'[DRY RUN] ' if args.dry_run else ''}Bumping version...")
        new_version = bump_version(bump_type, dry_run=args.dry_run)
        if new_version and not args.dry_run:
            print(f"\nSuccessfully bumped to version: {new_version}")
            print("\nNext steps:")
            print(f"  1. Review changes: git diff")
            print(f"  2. Push changes: git push && git push --tags")
            print(f"  3. Create GitHub release for v{new_version}")
    else:
        print("\nTo bump version, use --auto or --bump <type>")

    return 0


if __name__ == "__main__":
    sys.exit(main())
