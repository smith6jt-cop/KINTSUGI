"""
Snakemake wrapper for Signal Isolation Phase 1 (HITL workflow).

Extracts image features from all registered channels, clusters them by
similarity, suggests parameters for each cluster representative, and writes
a pending_approval.json checkpoint file for human review.

This is the first half of the two-phase HITL signal isolation workflow:
  Phase 1: Feature extraction + clustering + parameter suggestion → checkpoint
  (human reviews QC montages, approves/edits parameters)
  Phase 2: Applies approved parameters to all cluster members

CPU-only processing (sklearn, numpy — no GPU needed).
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

# Setup Python path
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Logging utilities
sys.path.insert(0, snakemake.scriptdir)
from log_utils import log_footer, log_header, log_info, log_warn
from log_utils import summary_before, summary_after

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
wf_config = snakemake.config
si_cfg = wf_config.get("signal_isolation", {})
TISSUE_TYPE = str(si_cfg.get("tissue_type", "unknown"))
PREDICTOR_MODEL = si_cfg.get("predictor_model_path", "")

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
log_header("signal_isolation_phase1", 0, PROJECT_DIR)
summary_before("signal_isolation_phase1", PROJECT_DIR)

print(f"\n{'='*60}")
print("Signal Isolation Phase 1: Clustering & Parameter Suggestion")
print(f"{'='*60}")
print(f"Project: {PROJECT_DIR.name}")
print(f"Tissue type: {TISSUE_TYPE}")
print(f"Input: {REGISTERED_DIR}")
print(f"Output: {OUTPUT_DIR}")

# ---------------------------------------------------------------------------
# Feature extraction
# ---------------------------------------------------------------------------
start_time = time.time()

try:
    import tifffile
    from kintsugi.signal.features import extract_channel_features, features_to_vector
    from kintsugi.signal.clustering import (
        cluster_channels,
        get_cluster_representatives,
        plot_cluster_summary,
    )

    # Discover channels from registered directory
    channel_features = {}
    channel_paths = {}

    for cycle_dir in sorted(REGISTERED_DIR.glob("cyc*")):
        if not cycle_dir.is_dir():
            continue
        for tif_path in sorted(cycle_dir.glob("*.tif")):
            ch_name = tif_path.stem
            # Skip DAPI and blank channels
            if ch_name.lower().startswith("dapi") or ch_name.lower().startswith("blank"):
                continue
            if ch_name in channel_features:
                continue  # Already processed from earlier cycle

            try:
                img = tifffile.imread(str(tif_path))
                if img.ndim > 2:
                    img = img.max(axis=0) if img.ndim == 3 else img[0]
                channel_features[ch_name] = extract_channel_features(
                    img, channel_name=ch_name
                )
                channel_paths[ch_name] = str(tif_path)
            except Exception as e:
                log_warn(f"Failed to extract features from {tif_path}: {e}")

    n_channels = len(channel_features)
    print(f"\nExtracted features from {n_channels} channels")

    if n_channels == 0:
        raise RuntimeError(f"No channels found in {REGISTERED_DIR}")

    # ---------------------------------------------------------------------------
    # Clustering
    # ---------------------------------------------------------------------------
    cluster_assignments = cluster_channels(
        channel_features,
        wavelength_aware=True,
    )

    representatives = get_cluster_representatives(channel_features, cluster_assignments)

    clusters_by_id = {}
    for ch, cid in cluster_assignments.items():
        clusters_by_id.setdefault(cid, []).append(ch)

    n_clusters = len(clusters_by_id)
    print(f"Clustered into {n_clusters} groups:")
    for cid in sorted(clusters_by_id):
        members = clusters_by_id[cid]
        rep_name = representatives[cid][0]
        print(f"  Cluster {cid}: {len(members)} channels, representative={rep_name}")

    # Save PCA scatter plot
    import matplotlib
    matplotlib.use("Agg")

    fig = plot_cluster_summary(
        channel_features,
        cluster_assignments,
        save_path=QC_DIR / "cluster_pca.png",
    )

    # ---------------------------------------------------------------------------
    # Parameter suggestion
    # ---------------------------------------------------------------------------
    # Try predictor first, fall back to defaults
    predictor = None
    if PREDICTOR_MODEL:
        try:
            from kintsugi.signal.predictor import ParameterPredictor
            model_path = Path(PREDICTOR_MODEL)
            if model_path.exists():
                predictor = ParameterPredictor(model_path=model_path)
                if predictor.is_trained:
                    log_info(
                        f"Loaded predictor ({predictor.n_training_examples} examples, "
                        f"R²={predictor.validation_r2:.2f})"
                    )
                else:
                    predictor = None
        except Exception as e:
            log_warn(f"Failed to load predictor: {e}")

    default_params = {
        "blank_params": {
            "blank_clip_factor": 5000,
            "blank_scale_factor": 1.0,
            "smooth_low": False,
            "low_size": 2,
            "low_percentile": 50,
            "smooth_high": False,
            "high_size": 2,
            "high_percentile": 50,
            "erosion": 1,
        },
        "clean_params": {
            "backgrnd_thresh": 200,
            "smooth_thresh": 0,
            "smooth": False,
            "remove_small": False,
            "small_size": 64,
            "footprint": 3,
        },
    }

    pending_clusters = []
    for cid in sorted(representatives):
        rep_name, member_count = representatives[cid]
        members = clusters_by_id[cid]

        suggested = None
        confidence = 0.0
        method = "default"

        if predictor:
            fv = features_to_vector(channel_features[rep_name])
            pred = predictor.predict(fv)
            if pred["status"] == "success":
                suggested = pred["predicted_params"]
                confidence = pred["confidence"]
                method = "predictor"

        if suggested is None:
            suggested = default_params
            confidence = 0.0
            method = "default"

        pending_clusters.append({
            "cluster_id": cid,
            "representative": rep_name,
            "member_count": member_count,
            "members": members,
            "suggested_params": suggested,
            "prediction_confidence": confidence,
            "suggestion_method": method,
            "qc_montage_path": str(QC_DIR / f"cluster_{cid}.png"),
            "approved": False,
            "approved_params": None,
        })

    # ---------------------------------------------------------------------------
    # Write checkpoint
    # ---------------------------------------------------------------------------
    pending = {
        "project": str(PROJECT_DIR),
        "tissue_type": TISSUE_TYPE,
        "created_at": datetime.now().isoformat(),
        "n_channels": n_channels,
        "n_clusters": n_clusters,
        "channel_paths": channel_paths,
        "clusters": pending_clusters,
        "instructions": (
            "Review each cluster's suggested_params and QC montage.\n"
            "Set 'approved': true for clusters you accept.\n"
            "Optionally edit 'approved_params' to override suggestions.\n"
            "Save this file, then run signal_isolation_phase2."
        ),
    }

    approval_file = OUTPUT_DIR / "pending_approval.json"
    with open(approval_file, "w") as f:
        json.dump(pending, f, indent=2, default=str)

    elapsed = (time.time() - start_time) / 60

    print(f"\n{'='*60}")
    print("Phase 1 Complete")
    print(f"{'='*60}")
    print(f"Channels: {n_channels}")
    print(f"Clusters: {n_clusters}")
    print(f"Approval file: {approval_file}")
    print(f"QC montages: {QC_DIR}/")
    print(f"Time: {elapsed:.1f} minutes")
    print("\nNext steps:")
    print(f"  1. Review {approval_file}")
    print("  2. Set 'approved': true for each cluster")
    print("  3. Run signal_isolation_phase2")

    # Write sentinel
    sentinel = Path(snakemake.output.sentinel)
    sentinel.parent.mkdir(parents=True, exist_ok=True)
    sentinel.write_text(
        f"stage=signal_isolation_phase1\n"
        f"completed={datetime.now().isoformat()}\n"
        f"n_channels={n_channels}\n"
        f"n_clusters={n_clusters}\n"
        f"approval_file={approval_file}\n"
        f"duration_minutes={elapsed:.1f}\n"
    )

    summary_after("signal_isolation_phase1", PROJECT_DIR, elapsed * 60, exit_code=0)
    log_footer(0)

except Exception as e:
    from log_utils import log_error
    log_error(f"Phase 1 failed: {e}")
    import traceback
    traceback.print_exc()

    elapsed = (time.time() - start_time) / 60
    summary_after("signal_isolation_phase1", PROJECT_DIR, elapsed * 60, exit_code=1)
    log_footer(1)
    raise
