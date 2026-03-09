"""
Channel Clustering and Parameter Propagation

Groups channels by image feature similarity so that parameters tuned on
one representative channel can be applied to all members of the same cluster.
Supports wavelength-aware partitioning, automatic cluster count selection,
and visual QC of propagation results.
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("kintsugi.signal.clustering")


def cluster_channels(
    channel_features: dict,
    n_clusters: int | None = None,
    method: str = "agglomerative",
    wavelength_aware: bool = True,
    min_cluster_size: int = 2,
    max_clusters: int = 20,
) -> dict:
    """
    Cluster channels by image feature similarity.

    Parameters
    ----------
    channel_features : dict
        Channel name -> feature dict (from extract_channel_features).
    n_clusters : int, optional
        Number of clusters. If None, auto-select via silhouette score
        testing range [2, min(max_clusters, n_channels // 2)].
    method : str
        'agglomerative' (default) or 'kmeans'.
    wavelength_aware : bool
        If True, enforce that channels from different wavelength groups
        cannot be in the same cluster by partitioning by wavelength group
        first, then clustering within each partition.
    min_cluster_size : int
        Minimum channels per cluster (currently unused, reserved for future).
    max_clusters : int
        Upper bound for auto-selection.

    Returns
    -------
    dict
        channel_name -> cluster_id (int)
    """
    from .features import features_to_vector

    names = list(channel_features.keys())
    if len(names) < 3:
        return {name: 0 for name in names}

    X = np.array([features_to_vector(channel_features[n]) for n in names])

    if wavelength_aware:
        return _cluster_wavelength_aware(
            names, X, channel_features, n_clusters, method, max_clusters
        )
    else:
        return _cluster_global(names, X, n_clusters, method, max_clusters)


def _cluster_wavelength_aware(
    names: list[str],
    X: np.ndarray,
    channel_features: dict,
    n_clusters: int | None,
    method: str,
    max_clusters: int,
) -> dict:
    """Partition by wavelength group, cluster within each partition."""
    from sklearn.preprocessing import StandardScaler

    wl_groups: dict[str, list[tuple[int, str]]] = {}
    for i, name in enumerate(names):
        wl = channel_features[name]["wavelength_group"]
        wl_groups.setdefault(wl, []).append((i, name))

    result = {}
    cluster_offset = 0
    for wl, members in wl_groups.items():
        if len(members) < 3:
            for _, name in members:
                result[name] = cluster_offset
            cluster_offset += 1
            continue

        indices = [m[0] for m in members]
        member_names = [m[1] for m in members]
        X_sub = X[indices]

        # Scale numeric features only (exclude one-hot wavelength suffix)
        scaler = StandardScaler()
        X_scaled = scaler.fit_transform(X_sub[:, :15])

        sub_n = _auto_n_clusters(X_scaled, n_clusters, method, max_clusters)
        labels = _fit_cluster_model(X_scaled, sub_n, method)

        for name, label in zip(member_names, labels):
            result[name] = int(label) + cluster_offset
        cluster_offset += sub_n

    return result


def _cluster_global(
    names: list[str],
    X: np.ndarray,
    n_clusters: int | None,
    method: str,
    max_clusters: int,
) -> dict:
    """Global clustering without wavelength partitioning."""
    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X[:, :15])

    auto_n = _auto_n_clusters(X_scaled, n_clusters, method, max_clusters)
    labels = _fit_cluster_model(X_scaled, auto_n, method)

    return {name: int(label) for name, label in zip(names, labels)}


def _fit_cluster_model(X: np.ndarray, n_clusters: int, method: str) -> np.ndarray:
    """Fit a clustering model and return labels."""
    from sklearn.cluster import AgglomerativeClustering, KMeans

    if method == "agglomerative":
        model = AgglomerativeClustering(n_clusters=n_clusters)
    else:
        model = KMeans(n_clusters=n_clusters, n_init=10, random_state=42)

    return model.fit_predict(X)


def _auto_n_clusters(
    X: np.ndarray,
    n_clusters: int | None,
    method: str,
    max_clusters: int,
) -> int:
    """Select optimal n_clusters via silhouette score if not specified."""
    from sklearn.metrics import silhouette_score

    if n_clusters is not None:
        return n_clusters

    max_k = min(max_clusters, X.shape[0] // 2)
    if max_k < 2:
        return 2

    best_k = 2
    best_score = -1.0
    for k in range(2, max_k + 1):
        labels = _fit_cluster_model(X, k, method)
        if len(set(labels)) < 2:
            continue
        score = silhouette_score(X, labels)
        if score > best_score:
            best_score = score
            best_k = k

    return best_k


def get_cluster_representatives(
    channel_features: dict,
    cluster_assignments: dict,
) -> dict:
    """
    For each cluster, select the channel closest to the cluster centroid.

    Parameters
    ----------
    channel_features : dict
        Channel name -> feature dict.
    cluster_assignments : dict
        Channel name -> cluster_id.

    Returns
    -------
    dict
        cluster_id -> (representative_channel_name, member_count)
    """
    from .features import features_to_vector

    clusters: dict[int, list[str]] = {}
    for name, cid in cluster_assignments.items():
        clusters.setdefault(cid, []).append(name)

    representatives = {}
    for cid, members in clusters.items():
        vectors = np.array([features_to_vector(channel_features[m]) for m in members])
        centroid = vectors.mean(axis=0)
        distances = np.linalg.norm(vectors - centroid, axis=1)
        rep_idx = int(np.argmin(distances))
        representatives[cid] = (members[rep_idx], len(members))

    return representatives


def propagate_cluster_parameters(
    marker_dict: dict,
    member_channels: list[str],
    params: dict,
    blank_map: dict,
    output_dir: Path | str,
    process_fn=None,
) -> dict:
    """
    Apply tuned parameters from a representative to all cluster members.

    Parameters
    ----------
    marker_dict : dict
        Channel name -> image data (numpy or dask).
    member_channels : list[str]
        Channel names to process with these parameters.
    params : dict
        Processing parameters with 'blank_params' and/or 'clean_params' dicts.
    blank_map : dict
        Channel name -> blank channel name.
    output_dir : Path or str
        Directory to save processed channel TIFFs.
    process_fn : callable, optional
        Custom processing function(signal, blank, params) -> np.ndarray.
        Default uses subtract_autofluorescence + clean_background.

    Returns
    -------
    dict
        channel_name -> {output_path, quality_score, quality_details, status, error}
    """
    import tifffile as tiff

    from kintsugi.qc.signal_quality import compute_signal_isolation_quality

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    results = {}

    for ch_name in member_channels:
        try:
            sig = marker_dict[ch_name]
            sig_np = sig.compute() if hasattr(sig, "compute") else np.array(sig)

            blank_name = blank_map.get(ch_name)
            blank_np = None
            if blank_name and blank_name in marker_dict:
                b = marker_dict[blank_name]
                blank_np = b.compute() if hasattr(b, "compute") else np.array(b)

            if process_fn is not None:
                processed = process_fn(sig_np, blank_np, params)
            else:
                processed = _default_process(sig_np, blank_np, params)

            qc = compute_signal_isolation_quality(sig_np, processed, blank=blank_np)

            out_path = output_dir / f"{ch_name}.tif"
            tiff.imwrite(str(out_path), processed.astype(np.uint16))

            results[ch_name] = {
                "output_path": out_path,
                "quality_score": qc["quality_score"],
                "quality_details": qc,
                "status": "success",
            }
        except Exception as e:
            logger.warning(f"Failed to process {ch_name}: {e}")
            results[ch_name] = {
                "output_path": None,
                "quality_score": 0.0,
                "status": "error",
                "error": str(e),
            }

    return results


def _default_process(
    sig_np: np.ndarray,
    blank_np: np.ndarray | None,
    params: dict,
) -> np.ndarray:
    """Apply standard subtraction + optional cleaning using the signal module API."""
    from kintsugi.signal.autofluorescence import subtract_autofluorescence

    bp = params.get("blank_params", params)

    if blank_np is None:
        return sig_np.copy()

    result = subtract_autofluorescence(
        sig_np,
        blank_np,
        blank_clip_factor=bp.get("blank_clip_factor", 0),
        blank_scale_factor=bp.get("blank_scale_factor", 1.0),
        smooth_low=bp.get("smooth_low", False),
        low_size=bp.get("low_size", 2),
        low_percentile=bp.get("low_percentile", 60),
        smooth_high=bp.get("smooth_high", False),
        high_size=bp.get("high_size", 2),
        high_percentile=bp.get("high_percentile", 90),
        erosion=bp.get("erosion", 0),
    )

    cp = params.get("clean_params")
    if cp:
        from kintsugi.signal.batch import CleanParams, clean_background

        clean_p = CleanParams(
            background_threshold=cp.get("backgrnd_thresh", cp.get("background_threshold", 0)),
            smooth=cp.get("smooth", False),
            smooth_threshold=cp.get("smooth_thresh", cp.get("smooth_threshold", 0)),
            footprint=cp.get("footprint", 1),
            remove_small=cp.get("remove_small", False),
            small_size=cp.get("small_size", 30),
        )
        result = clean_background(result, clean_p)

    return result


# --- Visualization ---


def plot_cluster_summary(
    channel_features: dict,
    cluster_assignments: dict,
    save_path: Path | str | None = None,
):
    """
    2D PCA scatter of channels colored by cluster assignment.

    Parameters
    ----------
    channel_features : dict
        Channel name -> feature dict.
    cluster_assignments : dict
        Channel name -> cluster_id.
    save_path : Path, optional
        If provided, save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler

    from .features import features_to_vector

    names = list(channel_features.keys())
    X = np.array([features_to_vector(channel_features[n]) for n in names])
    labels = [cluster_assignments[n] for n in names]

    X_scaled = StandardScaler().fit_transform(X[:, :15])
    pca = PCA(n_components=2)
    coords = pca.fit_transform(X_scaled)

    fig, ax = plt.subplots(1, 1, figsize=(10, 8))
    scatter = ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=labels,
        cmap="tab20",
        s=60,
        alpha=0.8,
        edgecolors="k",
        linewidth=0.5,
    )
    for i, name in enumerate(names):
        ax.annotate(
            name,
            (coords[i, 0], coords[i, 1]),
            fontsize=6,
            alpha=0.7,
            ha="center",
            va="bottom",
        )
    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]:.1%} var)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]:.1%} var)")
    ax.set_title(f"Channel Clusters ({len(set(labels))} clusters from {len(names)} channels)")
    plt.colorbar(scatter, ax=ax, label="Cluster ID")
    plt.tight_layout()

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches="tight")

    return fig


def generate_cluster_qc(
    marker_dict: dict,
    processing_results: dict,
    member_channels: list[str],
    cluster_id: int,
    save_path: Path | str | None = None,
    max_thumbnails: int = 8,
):
    """
    Generate a before/after QC montage for visual verification of a cluster.

    Shows original (top row) and processed (bottom row) thumbnails,
    annotated with channel name and quality score.

    Parameters
    ----------
    marker_dict : dict
        Channel name -> image data.
    processing_results : dict
        Channel name -> result dict from propagate_cluster_parameters.
    member_channels : list[str]
        Channel names in this cluster.
    cluster_id : int
        Cluster identifier for the title.
    save_path : Path, optional
        If provided, save figure to this path.
    max_thumbnails : int
        Maximum channels to show (default 8).

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    n = min(len(member_channels), max_thumbnails)
    channels = member_channels[:n]

    fig, axes = plt.subplots(2, n, figsize=(3 * n, 6))
    if n == 1:
        axes = axes.reshape(2, 1)

    for i, ch in enumerate(channels):
        # Original
        orig = marker_dict[ch]
        orig_np = orig.compute() if hasattr(orig, "compute") else np.array(orig)
        ds = max(1, orig_np.shape[0] // 256)
        thumb_orig = orig_np[::ds, ::ds]

        axes[0, i].imshow(
            thumb_orig,
            cmap="gray",
            vmin=np.percentile(thumb_orig, 1),
            vmax=np.percentile(thumb_orig, 99),
        )
        axes[0, i].set_title(ch, fontsize=8)
        axes[0, i].axis("off")

        # Processed
        result = processing_results.get(ch, {})
        if result.get("output_path") and Path(result["output_path"]).exists():
            import tifffile as tiff

            proc = tiff.imread(str(result["output_path"]))
            thumb_proc = proc[::ds, ::ds]
            score = result.get("quality_score", 0)
            color = "green" if score > 0.5 else "orange" if score > 0.3 else "red"
            axes[1, i].imshow(
                thumb_proc,
                cmap="gray",
                vmin=np.percentile(thumb_proc, 1),
                vmax=np.percentile(thumb_proc, 99),
            )
            axes[1, i].set_title(f"Q={score:.2f}", fontsize=8, color=color)
        else:
            axes[1, i].text(
                0.5,
                0.5,
                "FAILED",
                ha="center",
                va="center",
                transform=axes[1, i].transAxes,
                color="red",
            )
        axes[1, i].axis("off")

    axes[0, 0].set_ylabel("Original", fontsize=10)
    axes[1, 0].set_ylabel("Processed", fontsize=10)
    fig.suptitle(f"Cluster {cluster_id} — {len(member_channels)} channels", fontsize=12)
    plt.tight_layout()

    if save_path:
        fig.savefig(str(save_path), dpi=150, bbox_inches="tight")

    return fig
