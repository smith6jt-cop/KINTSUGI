"""
Snakemake wrapper for Signal Isolation Phase 2 (HITL workflow).

Reads the approved_parameters.json (or pending_approval.json with approved
clusters), then applies the approved parameters to all cluster members using
propagate_cluster_parameters(). Generates QC montages per cluster.

This is the second half of the two-phase HITL signal isolation workflow:
  Phase 1: Feature extraction + clustering + parameter suggestion → checkpoint
  (human reviews QC montages, approves/edits parameters)
  Phase 2: Applies approved parameters to all cluster members (this script)

CPU-only processing (numpy/scipy — no GPU needed).
"""

import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Snakemake interface
# ---------------------------------------------------------------------------
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CHANNELS = list(snakemake.params.channels)
CYCLES = list(snakemake.params.cycles)
CHANNEL_NAMES = dict(snakemake.params.channel_names)

# Setup Python path
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_footer, log_header, log_info, log_warn, log_error
from log_utils import summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
wf_config = snakemake.config
si_cfg = wf_config.get("signal_isolation", {})
TISSUE_TYPE = str(si_cfg.get("tissue_type", "unknown"))
LEARN = bool(si_cfg.get("learn", True))

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
DATA_DIR = PROJECT_DIR / "data" / "processed"
REGISTERED_DIR = DATA_DIR / "registered"
OUTPUT_DIR = DATA_DIR / "signal_isolated"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

QC_DIR = OUTPUT_DIR / "cluster_qc"
QC_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
log_header("signal_isolation_phase2", 0, PROJECT_DIR)
summary_before("signal_isolation_phase2", PROJECT_DIR)

print(f"\n{'='*60}")
print("Signal Isolation Phase 2: Apply Approved Parameters")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Tissue type: {TISSUE_TYPE}")
print(f"Input: {REGISTERED_DIR}")
print(f"Output: {OUTPUT_DIR}")

# ---------------------------------------------------------------------------
# Load approval file
# ---------------------------------------------------------------------------
start_time = time.time()

# Check for approved_parameters.json first, fall back to pending_approval.json
approved_file = OUTPUT_DIR / "approved_parameters.json"
pending_file = OUTPUT_DIR / "pending_approval.json"

if approved_file.exists():
    approval_path = approved_file
elif pending_file.exists():
    approval_path = pending_file
else:
    log_error("No approval file found. Run signal_isolation_phase1 first.")
    log_error(f"Expected: {approved_file} or {pending_file}")
    log_footer(1)
    sys.exit(1)

print(f"Loading approval file: {approval_path}")

with open(approval_path) as f:
    approval_data = json.load(f)

clusters = approval_data.get("clusters", [])
approved_clusters = [c for c in clusters if c.get("approved", False)]

if not approved_clusters:
    log_error("No clusters have been approved. Edit the approval file first.")
    log_error(f"Set 'approved': true in {approval_path}")
    log_footer(1)
    sys.exit(1)

total_clusters = len(clusters)
n_approved = len(approved_clusters)
print(f"Approved clusters: {n_approved}/{total_clusters}")

# ---------------------------------------------------------------------------
# Load channel images
# ---------------------------------------------------------------------------
try:
    import tifffile
    import matplotlib
    matplotlib.use("Agg")

    from kintsugi.signal.clustering import propagate_cluster_parameters, generate_cluster_qc

    # Build marker_dict from registered images
    marker_dict = {}
    for cycle_dir in sorted(REGISTERED_DIR.glob("cyc*")):
        if not cycle_dir.is_dir():
            continue
        for tif_path in sorted(cycle_dir.glob("*.tif")):
            ch_name = tif_path.stem
            if ch_name not in marker_dict:
                marker_dict[ch_name] = tifffile.imread(str(tif_path))

    print(f"Loaded {len(marker_dict)} channels from {REGISTERED_DIR}")

    # Build blank map from channel names
    blank_channels = [ch for ch in marker_dict if "blank" in ch.lower()]
    default_blank = blank_channels[0] if blank_channels else None
    blank_map = {}
    if default_blank:
        for ch in marker_dict:
            if "blank" not in ch.lower() and "dapi" not in ch.lower():
                blank_map[ch] = default_blank
        print(f"Default blank: {default_blank}")

    # ---------------------------------------------------------------------------
    # Process each approved cluster
    # ---------------------------------------------------------------------------
    all_results = {}
    total_processed = 0
    total_failed = 0

    for cluster in approved_clusters:
        cid = cluster["cluster_id"]
        params = cluster.get("approved_params") or cluster["suggested_params"]
        members = cluster["members"]

        # Filter to signal channels only
        signal_members = [
            ch for ch in members
            if "dapi" not in ch.lower() and "blank" not in ch.lower()
            and ch in marker_dict
        ]

        if not signal_members:
            print(f"\nCluster {cid}: no signal channels to process")
            continue

        print(f"\nCluster {cid}: processing {len(signal_members)} channels")

        results = propagate_cluster_parameters(
            marker_dict=marker_dict,
            member_channels=signal_members,
            params=params,
            blank_map=blank_map,
            output_dir=OUTPUT_DIR,
        )
        all_results.update(results)

        # Generate QC montage
        fig = generate_cluster_qc(
            marker_dict=marker_dict,
            processing_results=results,
            member_channels=signal_members,
            cluster_id=cid,
            save_path=QC_DIR / f"cluster_{cid}_processed.png",
        )
        import matplotlib.pyplot as plt
        plt.close(fig)

        for ch, r in results.items():
            status = r["status"]
            score = r.get("quality_score", 0)
            print(f"  {ch}: {status} (quality={score:.3f})")
            if status == "success":
                total_processed += 1
            else:
                total_failed += 1

    # ---------------------------------------------------------------------------
    # Parameter learning (optional)
    # ---------------------------------------------------------------------------
    if LEARN and TISSUE_TYPE != "unknown":
        try:
            from kintsugi.claude.parameter_learning import ParameterLearningEngine

            engine = ParameterLearningEngine(str(PROJECT_DIR))
            recorded = 0

            for cluster in approved_clusters:
                params = cluster.get("approved_params") or cluster["suggested_params"]
                bp = params.get("blank_params", params)

                for ch in cluster["members"]:
                    r = all_results.get(ch)
                    if r and r["status"] == "success" and r["quality_score"] > 0.5:
                        engine.record_parameters(
                            operation="blank_subtraction",
                            tissue_type=TISSUE_TYPE,
                            marker_name=ch,
                            parameters=bp,
                            quality_score=r["quality_score"],
                        )
                        recorded += 1

            log_info(f"Recorded {recorded} parameter sets to learning database")
        except Exception as e:
            log_warn(f"Parameter learning failed: {e}")

    # ---------------------------------------------------------------------------
    # Summary
    # ---------------------------------------------------------------------------
    elapsed = (time.time() - start_time) / 60
    scores = [r["quality_score"] for r in all_results.values() if r["status"] == "success"]
    mean_quality = np.mean(scores) if scores else 0

    n_files = len(list(OUTPUT_DIR.glob("*.tif")))

    print(f"\n{'='*60}")
    print("Phase 2 Complete")
    print(f"{'='*60}")
    print(f"Clusters processed: {n_approved}")
    print(f"Channels processed: {total_processed}")
    print(f"Channels failed: {total_failed}")
    print(f"Mean quality: {mean_quality:.3f}")
    print(f"Output files: {n_files}")
    print(f"Time: {elapsed:.1f} minutes")

    summary_after("signal_isolation_phase2", PROJECT_DIR, elapsed * 60, exit_code=0)

    # Write sentinel
    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=signal_isolation_phase2\n"
        f"completed={datetime.now().isoformat()}\n"
        f"mode=hitl\n"
        f"tissue_type={TISSUE_TYPE}\n"
        f"clusters_approved={n_approved}\n"
        f"channels_processed={total_processed}\n"
        f"channels_failed={total_failed}\n"
        f"mean_quality={mean_quality:.3f}\n"
        f"files={n_files}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )
    print(f"Sentinel written: {sentinel}")
    log_footer(0)

except Exception as e:
    log_error(f"Phase 2 failed: {e}")
    import traceback
    traceback.print_exc()

    elapsed = (time.time() - start_time) / 60
    summary_after("signal_isolation_phase2", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    raise
