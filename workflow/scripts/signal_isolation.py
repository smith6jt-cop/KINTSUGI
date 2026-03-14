"""
Snakemake wrapper for KINTSUGI signal isolation (autofluorescence subtraction).

Processes all registered channels via batch signal isolation with recipe-driven
or auto-analyzed parameter optimization. Runs as an aggregate rule (no {cycle}
wildcard) after registration completes — like registration itself.

CPU-only processing (numpy/scipy, no GPU needed).

Snakemake guarantees:
  - All registered outputs exist before this script runs
  - No need for wait_for_input() polling

Skip-existing: If manifest exists and all marker TIFs are present, the entire
job is skipped. Re-run with force=True in config to reprocess.
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CHANNELS = list(snakemake.params.channels)
CYCLES = list(snakemake.params.cycles)
CHANNEL_NAMES = dict(snakemake.params.channel_names)

# Setup Python path (same 3-tier as other wrapper scripts)
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_footer, log_header, log_info, log_warn, log_error
from log_utils import summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration from config.yaml
# ---------------------------------------------------------------------------
wf_config = snakemake.config
si_cfg = wf_config.get("signal_isolation", {})

METHOD = str(si_cfg.get("method", "auto"))
TISSUE_TYPE = str(si_cfg.get("tissue_type", "unknown"))
TILE_SMOOTH_SIGMA = float(si_cfg.get("tile_smooth_sigma", 0.0))
RECIPE_DIR = si_cfg.get("recipe_dir", "")
LEARN = bool(si_cfg.get("learn", True))
FORCE = bool(si_cfg.get("force", False))
QUALITY_GATE = float(si_cfg.get("quality_gate", 0.6))
CLIP_PERCENTILE = float(si_cfg.get("clip_percentile", 99.5))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data" / "processed"
REGISTERED_DIR = DATA_DIR / "registered"
OUTPUT_DIR = DATA_DIR / "signal_isolated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Structured logging header
# ---------------------------------------------------------------------------
log_header("signal_isolation", 0, PROJECT_DIR)
summary_before("signal_isolation", PROJECT_DIR)

print(f"\n{'='*60}")
print(f"Batch Signal Isolation")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Cycles: {CYCLES}")
print(f"Channels: {CHANNELS}")
print(f"Method: {METHOD}")
print(f"Tissue type: {TISSUE_TYPE}")
print(f"Tile smooth sigma: {TILE_SMOOTH_SIGMA}")
print(f"Recipe dir: {RECIPE_DIR or '(none)'}")
print(f"Learn: {LEARN}")
print(f"Force: {FORCE}")
print(f"Input: {REGISTERED_DIR}")
print(f"Output: {OUTPUT_DIR}")

# ---------------------------------------------------------------------------
# Skip-existing check
# ---------------------------------------------------------------------------
manifest_path = OUTPUT_DIR / "signal_isolation_manifest.json"


def check_all_isolated():
    """Check if manifest exists and all marker TIFs are present."""
    if not manifest_path.exists():
        return False
    try:
        with open(manifest_path) as f:
            manifest = json.load(f)
        channels_data = manifest.get("channels", {})
        if not channels_data:
            return False
        for marker, info in channels_data.items():
            if info.get("status") == "error":
                continue
            tif_path = OUTPUT_DIR / f"{marker}.tif"
            if not tif_path.exists():
                return False
        return True
    except (json.JSONDecodeError, KeyError):
        return False


# ---------------------------------------------------------------------------
# Main processing
# ---------------------------------------------------------------------------
start_time = time.time()

if not FORCE and check_all_isolated():
    print("\nAll signal-isolated files already exist — nothing to do")
    elapsed = (time.time() - start_time) / 60
    n_files = len(list(OUTPUT_DIR.glob("*.tif")))

    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=signal_isolation\n"
        f"completed={datetime.now().isoformat()}\n"
        f"skipped=all\n"
        f"files={n_files}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    summary_after("signal_isolation", PROJECT_DIR, elapsed * 60, exit_code=0)
    log_footer(0)
    sys.exit(0)

# ---------------------------------------------------------------------------
# Recipe auto-discovery
# ---------------------------------------------------------------------------
recipe_dir_resolved = None
if RECIPE_DIR:
    recipe_dir_resolved = Path(RECIPE_DIR)
    if not recipe_dir_resolved.exists():
        log_warn(f"Configured recipe_dir not found: {RECIPE_DIR}")
        recipe_dir_resolved = None

# Search standard locations if not explicitly configured
if recipe_dir_resolved is None:
    search_paths = [
        DATA_DIR / "Processing_parameters",
        PROJECT_DIR / "Processing_parameters",
        PROJECT_DIR / "meta" / "Processing_parameters",
    ]
    for sp in search_paths:
        if sp.exists() and list(sp.glob("*_param.txt")):
            recipe_dir_resolved = sp
            log_info(f"Auto-discovered recipes: {sp}")
            break

if recipe_dir_resolved:
    print(f"Using recipes from: {recipe_dir_resolved}")
else:
    print("No recipes found — using auto-analysis mode")

# ---------------------------------------------------------------------------
# Run batch signal isolation
# ---------------------------------------------------------------------------
try:
    from kintsugi.signal.batch import process_batch

    result = process_batch(
        project_dir=str(PROJECT_DIR),
        method=METHOD,
        tissue_type=TISSUE_TYPE if TISSUE_TYPE != "unknown" else None,
        tile_smooth_sigma=TILE_SMOOTH_SIGMA,
        output_dir=str(OUTPUT_DIR),
        force=FORCE,
        recipe_dir=str(recipe_dir_resolved) if recipe_dir_resolved else None,
        learn=LEARN,
        quality_gate=QUALITY_GATE,
        clip_percentile=CLIP_PERCENTILE,
    )

    # Print summary
    summary = result.summary if hasattr(result, "summary") else {}
    total = summary.get("total", 0)
    n_global = summary.get("global", 0)
    n_weighted = summary.get("weighted", 0)
    n_recipe = summary.get("recipe", 0)
    n_skipped = summary.get("skipped", 0)
    n_errors = summary.get("errors", 0)
    n_failed = summary.get("failed", 0)
    mean_quality = summary.get("mean_quality", 0.0)

    elapsed = (time.time() - start_time) / 60

    print(f"\n{'='*60}")
    print(f"Signal Isolation Complete")
    print(f"{'='*60}")
    print(f"Total channels: {total}")
    print(f"  Global: {n_global}")
    print(f"  Weighted: {n_weighted}")
    print(f"  Recipe: {n_recipe}")
    print(f"  Skipped: {n_skipped}")
    print(f"  Errors: {n_errors}")
    print(f"  Failed (over-subtracted): {n_failed}")
    print(f"Mean quality: {mean_quality:.3f}")
    print(f"Time: {elapsed:.1f} minutes")

    n_files = len(list(OUTPUT_DIR.glob("*.tif")))

    summary_after("signal_isolation", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Write sentinel with metadata
    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=signal_isolation\n"
        f"completed={datetime.now().isoformat()}\n"
        f"method={METHOD}\n"
        f"tissue_type={TISSUE_TYPE}\n"
        f"recipe_dir={recipe_dir_resolved or ''}\n"
        f"total={total}\n"
        f"global={n_global}\n"
        f"weighted={n_weighted}\n"
        f"recipe={n_recipe}\n"
        f"skipped={n_skipped}\n"
        f"errors={n_errors}\n"
        f"failed={n_failed}\n"
        f"mean_quality={mean_quality:.3f}\n"
        f"files={n_files}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    log_footer(0)

except Exception as e:
    log_error(f"Signal isolation failed: {e}")
    import traceback

    traceback.print_exc()

    elapsed = (time.time() - start_time) / 60
    summary_after("signal_isolation", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    raise
