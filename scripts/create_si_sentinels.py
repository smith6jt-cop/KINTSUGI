#!/usr/bin/env python3
"""Validate signal isolation output and create missing .snakemake_complete sentinels.

Usage:
    python scripts/create_si_sentinels.py ../KINTSUGI_Projects --dry-run
    python scripts/create_si_sentinels.py ../KINTSUGI_Projects
"""

import argparse
import json
import sys
from pathlib import Path


def validate_and_create_sentinel(project_dir: Path, dry_run: bool = True) -> dict:
    """Validate SI output for a project and optionally create sentinel.

    Returns dict with validation results.
    """
    result = {
        "project": project_dir.name,
        "valid": False,
        "sentinel_created": False,
        "errors": [],
    }

    si_dir = project_dir / "data" / "processed" / "signal_isolated"
    manifest_path = si_dir / "signal_isolation_manifest.json"
    sentinel_path = si_dir / ".snakemake_complete"

    # Skip if sentinel already exists
    if sentinel_path.exists():
        result["errors"].append("sentinel already exists")
        return result

    # Check manifest exists and is valid JSON
    if not manifest_path.exists():
        result["errors"].append("manifest not found")
        return result

    try:
        manifest = json.loads(manifest_path.read_text())
    except (json.JSONDecodeError, OSError) as e:
        result["errors"].append(f"manifest invalid: {e}")
        return result

    # Extract manifest metadata
    channels = manifest.get("channels", {})
    summary = manifest.get("summary", {})
    timestamp = manifest.get("timestamp", "")
    tissue_type = manifest.get("tissue_type", "unknown")
    method_override = manifest.get("method_override", "auto")
    recipe_dir = manifest.get("recipe_dir", "")

    result["total_channels"] = len(channels)
    result["timestamp"] = timestamp
    result["tissue_type"] = tissue_type

    if not channels:
        result["errors"].append("manifest has no channels")
        return result

    # Verify every channel in manifest has a corresponding TIF file
    missing_tifs = []
    zero_size_tifs = []
    for ch_key, ch_data in channels.items():
        marker = ch_data.get("marker_name", ch_key)
        cycle = ch_data.get("cycle", "")

        # Find matching TIF: pattern is typically "{marker}_cyc{NN}.tif" or similar
        # Check for any TIF that matches this channel
        tif_candidates = list(si_dir.glob(f"*{marker}*.tif"))
        if not tif_candidates:
            # Also try by channel key directly
            tif_candidates = list(si_dir.glob(f"*{ch_key}*.tif"))

        if not tif_candidates:
            missing_tifs.append(f"{marker} (cycle {cycle})")
        else:
            for tif in tif_candidates:
                if tif.stat().st_size == 0:
                    zero_size_tifs.append(tif.name)

    if missing_tifs:
        result["errors"].append(f"missing TIFs: {', '.join(missing_tifs)}")

    if zero_size_tifs:
        result["errors"].append(f"zero-size TIFs: {', '.join(zero_size_tifs)}")

    # Count actual TIF files
    all_tifs = list(si_dir.glob("*.tif"))
    result["tif_count"] = len(all_tifs)

    # Verify TIF count matches manifest channel count
    expected = len(channels)
    actual = len(all_tifs)
    if actual < expected:
        result["errors"].append(f"TIF count mismatch: {actual} files vs {expected} channels")

    # Extract method counts from summary
    n_global = summary.get("global", 0)
    n_weighted = summary.get("weighted", 0)
    n_recipe = summary.get("recipe", 0)
    n_skipped = summary.get("skipped", 0)
    n_errors = summary.get("error", 0)
    mean_quality = summary.get("mean_quality", 0.0)
    total = summary.get("total", len(channels))

    result["summary"] = {
        "total": total,
        "global": n_global,
        "weighted": n_weighted,
        "recipe": n_recipe,
        "skipped": n_skipped,
        "errors": n_errors,
        "mean_quality": mean_quality,
    }

    # If there are blocking errors, don't create sentinel
    if result["errors"]:
        return result

    result["valid"] = True

    # Create sentinel
    if not dry_run:
        sentinel_content = (
            f"stage=signal_isolation\n"
            f"completed={timestamp}\n"
            f"method={method_override}\n"
            f"tissue_type={tissue_type}\n"
            f"recipe_dir={recipe_dir}\n"
            f"total={total}\n"
            f"global={n_global}\n"
            f"weighted={n_weighted}\n"
            f"recipe={n_recipe}\n"
            f"skipped={n_skipped}\n"
            f"errors={n_errors}\n"
            f"mean_quality={mean_quality:.3f}\n"
            f"files={actual}\n"
            f"duration_minutes=0.0\n"
        )
        sentinel_path.write_text(sentinel_content)
        result["sentinel_created"] = True

    return result


def main():
    parser = argparse.ArgumentParser(
        description="Validate SI output and create missing sentinels"
    )
    parser.add_argument("projects_dir", type=Path, help="Directory containing projects")
    parser.add_argument(
        "--dry-run", action="store_true", help="Preview without writing sentinels"
    )
    args = parser.parse_args()

    projects_dir = args.projects_dir.resolve()
    if not projects_dir.is_dir():
        print(f"Error: {projects_dir} is not a directory")
        sys.exit(1)

    # Find projects with SI manifest but no sentinel
    candidates = []
    for proj_dir in sorted(projects_dir.iterdir()):
        if not proj_dir.is_dir():
            continue
        si_dir = proj_dir / "data" / "processed" / "signal_isolated"
        manifest = si_dir / "signal_isolation_manifest.json"
        sentinel = si_dir / ".snakemake_complete"
        if manifest.exists() and not sentinel.exists():
            candidates.append(proj_dir)

    if not candidates:
        print("No projects found with SI manifest but missing sentinel.")
        return

    print(f"Found {len(candidates)} projects with manifest but no sentinel")
    if args.dry_run:
        print("=== DRY RUN (no files will be written) ===")
    print()

    valid_count = 0
    error_count = 0

    for proj_dir in candidates:
        result = validate_and_create_sentinel(proj_dir, dry_run=args.dry_run)
        name = result["project"]

        if result["valid"]:
            valid_count += 1
            s = result["summary"]
            status = "WOULD CREATE" if args.dry_run else "CREATED"
            print(
                f"  [OK] {name}: {result['tif_count']} TIFs, "
                f"quality={s['mean_quality']:.3f}, "
                f"global={s['global']}/weighted={s['weighted']}/recipe={s['recipe']} "
                f"-> {status} sentinel"
            )
        else:
            error_count += 1
            errors = "; ".join(result["errors"])
            print(f"  [SKIP] {name}: {errors}")

    print()
    print(f"Summary: {valid_count} valid, {error_count} skipped/errors")
    if args.dry_run and valid_count > 0:
        print(f"\nRe-run without --dry-run to create {valid_count} sentinels.")


if __name__ == "__main__":
    main()
