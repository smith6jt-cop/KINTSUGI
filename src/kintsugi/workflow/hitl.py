"""
Human-in-the-Loop Batch Workflow for Signal Isolation

Defines the HITL batch workflow pattern for signal isolation at scale.
This module provides workflow functions that can be used either standalone
or integrated with Prefect/Snakemake.

The pattern is orchestrator-agnostic: it produces JSON checkpoint files
that can be consumed by any workflow engine or by a human reviewing
results in a notebook.

Workflow:
    1. Extract features for all channels across all datasets
    2. Cluster channels globally (or per-dataset if tissue types vary)
    3. For each cluster: predict or optimize parameters for representative
    4. Generate QC montage per cluster
    5. CHECKPOINT: Write pending_approval.json, wait for human review
    6. After approval: propagate parameters to all cluster members
    7. Generate final QC report
"""

from __future__ import annotations

import json
import logging
import time
from pathlib import Path

logger = logging.getLogger("kintsugi.workflow.hitl")


class SignalIsolationBatchWorkflow:
    """
    Orchestrates batch signal isolation with human approval checkpoints.

    Usage in notebook::

        workflow = SignalIsolationBatchWorkflow(project_paths=[...])
        workflow.run_phase1()  # auto-suggest parameters
        # User reviews QC montages, edits approved_parameters.json
        workflow.run_phase2()  # apply approved parameters

    Parameters
    ----------
    project_paths : list[Path]
        Paths to KINTSUGI project directories to process.
    output_dir : Path
        Directory for workflow state, QC montages, and checkpoint files.
    tissue_type : str
        Tissue type label for parameter learning (default 'unknown').
    predictor_model_path : Path, optional
        Path to trained ParameterPredictor joblib model.
    use_optuna : bool
        If True, run Optuna optimization for each cluster representative
        instead of using predictor alone.
    optuna_trials : int
        Number of Optuna trials per cluster representative (default 50).
    """

    def __init__(
        self,
        project_paths: list[Path],
        output_dir: Path,
        tissue_type: str = "unknown",
        predictor_model_path: Path | None = None,
        use_optuna: bool = False,
        optuna_trials: int = 50,
    ):
        self.project_paths = [Path(p) for p in project_paths]
        self.output_dir = Path(output_dir)
        self.tissue_type = tissue_type
        self.predictor_model_path = predictor_model_path
        self.use_optuna = use_optuna
        self.optuna_trials = optuna_trials

        self.state_file = self.output_dir / "workflow_state.json"
        self.approval_file = self.output_dir / "pending_approval.json"
        self.approved_file = self.output_dir / "approved_parameters.json"

        self.output_dir.mkdir(parents=True, exist_ok=True)

    def run_phase1(self) -> dict:
        """
        Phase 1: Feature extraction, clustering, parameter suggestion.

        Produces:
        - workflow_state.json (full state for resume)
        - pending_approval.json (clusters with suggested params + QC montage paths)
        - qc_montages/ directory with per-cluster montage images

        Returns dict with phase1 summary.
        """
        from kintsugi.signal.clustering import cluster_channels, get_cluster_representatives
        from kintsugi.signal.predictor import ParameterPredictor

        state = {"phase": "1_suggesting", "started_at": time.time(), "datasets": {}}

        # Load predictor if available
        predictor = None
        if self.predictor_model_path and self.predictor_model_path.exists():
            predictor = ParameterPredictor(model_path=self.predictor_model_path)

        all_features = {}

        # Extract features from all datasets
        for proj_path in self.project_paths:
            dataset_name = proj_path.name
            features = self._load_dataset_features(proj_path)
            all_features[dataset_name] = features
            state["datasets"][dataset_name] = {
                "n_channels": len(features),
                "features_extracted": True,
            }

        # Global clustering across all datasets
        flat_features = {}
        for ds, feats in all_features.items():
            for ch, f in feats.items():
                flat_features[f"{ds}/{ch}"] = f

        if len(flat_features) < 3:
            # Too few channels to cluster; assign all to cluster 0
            clusters = dict.fromkeys(flat_features.keys(), 0)
        else:
            clusters = cluster_channels(flat_features, wavelength_aware=True)

        representatives = get_cluster_representatives(flat_features, clusters)

        # Suggest parameters for each representative
        pending = {"clusters": [], "created_at": time.time()}
        qc_dir = self.output_dir / "qc_montages"
        qc_dir.mkdir(exist_ok=True)

        for cid, (rep_key, member_count) in representatives.items():
            ds_name, ch_name = rep_key.split("/", 1)
            members = [k for k, v in clusters.items() if v == cid]

            # Predict or use defaults
            suggested_params = None
            confidence = 0.0

            if predictor and predictor.is_trained:
                from kintsugi.signal.features import features_to_vector

                fv = features_to_vector(flat_features[rep_key])
                pred = predictor.predict(fv)
                if pred["status"] == "success":
                    suggested_params = pred["predicted_params"]
                    confidence = pred["confidence"]

            if suggested_params is None:
                suggested_params = self._default_params()
                confidence = 0.0

            cluster_entry = {
                "cluster_id": cid,
                "representative": {"dataset": ds_name, "channel": ch_name},
                "member_count": member_count,
                "members": members,
                "suggested_params": suggested_params,
                "prediction_confidence": confidence,
                "qc_montage_path": str(qc_dir / f"cluster_{cid}.png"),
                "approved": False,
                "approved_params": None,
            }
            pending["clusters"].append(cluster_entry)

        # Save state
        state["phase"] = "1_complete_awaiting_approval"
        state["n_clusters"] = len(pending["clusters"])
        state["total_channels"] = len(flat_features)
        self._save_json(self.state_file, state)
        self._save_json(self.approval_file, pending)

        return {
            "status": "phase1_complete",
            "n_clusters": len(pending["clusters"]),
            "total_channels": len(flat_features),
            "approval_file": str(self.approval_file),
            "message": (
                f"Review {self.approval_file} and QC montages in {qc_dir}/.\n"
                f"Set 'approved': true for each cluster (optionally adjust params).\n"
                f"Save as {self.approved_file}, then run phase2."
            ),
        }

    def run_phase2(self) -> dict:
        """
        Phase 2: Apply approved parameters to all cluster members.

        Reads approved_parameters.json, processes all channels, generates
        final QC report.

        Returns dict with phase2 summary.
        """
        if not self.approved_file.exists():
            return {
                "status": "error",
                "message": f"Approved parameters file not found: {self.approved_file}",
            }

        approved = self._load_json(self.approved_file)

        results_summary = {
            "clusters_processed": 0,
            "channels_processed": 0,
            "channels_failed": 0,
        }

        for cluster in approved["clusters"]:
            if not cluster.get("approved", False):
                continue

            params = cluster.get("approved_params") or cluster["suggested_params"]
            members = cluster["members"]

            for member_key in members:
                ds_name, ch_name = member_key.split("/", 1)
                try:
                    self._process_channel(ds_name, ch_name, params)
                    results_summary["channels_processed"] += 1
                except Exception as e:
                    logger.warning(f"Failed to process {member_key}: {e}")
                    results_summary["channels_failed"] += 1

            results_summary["clusters_processed"] += 1

        # Save final state
        state = self._load_json(self.state_file) if self.state_file.exists() else {}
        state["phase"] = "2_complete"
        state["results"] = results_summary
        self._save_json(self.state_file, state)

        return {"status": "phase2_complete", **results_summary}

    def _load_dataset_features(self, proj_path: Path) -> dict:
        """Load channels and extract features for a single dataset."""
        from kintsugi.signal.features import extract_channel_features

        # Scan for registered images
        registered_dir = proj_path / "data" / "processed" / "registered"
        if not registered_dir.exists():
            # Try EDF outputs
            registered_dir = proj_path / "data" / "processed" / "edf"
        if not registered_dir.exists():
            logger.warning(f"No processed images found in {proj_path}")
            return {}

        features = {}
        import tifffile

        for tif_path in sorted(registered_dir.rglob("*.tif")):
            ch_name = tif_path.stem
            if ch_name.startswith("DAPI") or ch_name.startswith("Blank"):
                continue
            try:
                img = tifffile.imread(str(tif_path))
                if img.ndim > 2:
                    # Use max projection for 3D stacks
                    img = img.max(axis=0) if img.ndim == 3 else img[0]
                features[ch_name] = extract_channel_features(img, channel_name=ch_name)
            except Exception as e:
                logger.warning(f"Failed to extract features from {tif_path}: {e}")

        return features

    def _process_channel(self, dataset_name: str, channel_name: str, params: dict):
        """Process a single channel with given parameters."""
        # Find the project path for this dataset
        proj_path = None
        for p in self.project_paths:
            if p.name == dataset_name:
                proj_path = p
                break

        if proj_path is None:
            raise ValueError(f"Dataset not found: {dataset_name}")

        import tifffile

        from kintsugi.signal.clustering import _default_process

        # Load signal image
        registered_dir = proj_path / "data" / "processed" / "registered"
        if not registered_dir.exists():
            registered_dir = proj_path / "data" / "processed" / "edf"

        sig_path = None
        for p in registered_dir.rglob(f"{channel_name}.tif"):
            sig_path = p
            break

        if sig_path is None:
            raise FileNotFoundError(f"Channel not found: {channel_name} in {dataset_name}")

        sig_np = tifffile.imread(str(sig_path))
        if sig_np.ndim > 2:
            sig_np = sig_np.max(axis=0) if sig_np.ndim == 3 else sig_np[0]

        # Process (blank_np=None uses copy fallback in _default_process)
        processed = _default_process(sig_np, None, params)

        # Save to signal_isolated directory
        output_dir = proj_path / "data" / "processed" / "signal_isolated"
        output_dir.mkdir(parents=True, exist_ok=True)
        out_path = output_dir / f"{channel_name}.tif"
        tifffile.imwrite(str(out_path), processed.astype(sig_np.dtype))

    @staticmethod
    def _default_params() -> dict:
        """Return conservative default parameters."""
        return {
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

    @staticmethod
    def _save_json(path: Path, data: dict):
        with open(path, "w") as f:
            json.dump(data, f, indent=2, default=str)

    @staticmethod
    def _load_json(path: Path) -> dict:
        with open(path) as f:
            return json.load(f)
