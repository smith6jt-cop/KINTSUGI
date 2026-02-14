"""
KRONOS Embedding Analysis Utilities

Downstream analysis of KRONOS embeddings extracted from KINTSUGI projects:
- Dimensionality reduction (UMAP, PCA)
- Clustering (Leiden, KMeans)
- Spatial visualization of embedding clusters
- Spatial search (find similar regions across datasets)
- Cross-dataset patient stratification
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import numpy as np

from .model import EmbeddingResult

logger = logging.getLogger(__name__)


def cluster_embeddings(
    embeddings: np.ndarray | EmbeddingResult,
    method: str = "leiden",
    n_clusters: int | None = None,
    resolution: float = 1.0,
    n_neighbors: int = 15,
    random_state: int = 42,
) -> np.ndarray:
    """Cluster patch embeddings.

    Parameters
    ----------
    embeddings : np.ndarray or EmbeddingResult
        Patch embeddings, shape (N, D). If EmbeddingResult, uses patch_embeddings.
    method : str
        Clustering method: "leiden", "kmeans", "spectral".
    n_clusters : int, optional
        Number of clusters (required for kmeans/spectral, ignored for leiden).
    resolution : float
        Resolution parameter for Leiden clustering.
    n_neighbors : int
        Number of neighbors for graph construction (Leiden).
    random_state : int
        Random seed.

    Returns
    -------
    np.ndarray
        Cluster labels, shape (N,).
    """
    if isinstance(embeddings, EmbeddingResult):
        embeddings = embeddings.patch_embeddings

    if method == "leiden":
        return _cluster_leiden(embeddings, resolution, n_neighbors, random_state)
    elif method == "kmeans":
        return _cluster_kmeans(embeddings, n_clusters or 10, random_state)
    elif method == "spectral":
        return _cluster_spectral(embeddings, n_clusters or 10, n_neighbors, random_state)
    else:
        raise ValueError(f"Unknown clustering method: {method}")


def _cluster_leiden(
    embeddings: np.ndarray,
    resolution: float,
    n_neighbors: int,
    random_state: int,
) -> np.ndarray:
    """Leiden clustering via scanpy."""
    import anndata as ad
    import scanpy as sc

    adata = ad.AnnData(embeddings)
    sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep="X", random_state=random_state)
    sc.tl.leiden(adata, resolution=resolution, random_state=random_state)
    return adata.obs["leiden"].values.astype(int)


def _cluster_kmeans(
    embeddings: np.ndarray,
    n_clusters: int,
    random_state: int,
) -> np.ndarray:
    """KMeans clustering via sklearn."""
    from sklearn.cluster import KMeans

    km = KMeans(n_clusters=n_clusters, random_state=random_state, n_init=10)
    return km.fit_predict(embeddings)


def _cluster_spectral(
    embeddings: np.ndarray,
    n_clusters: int,
    n_neighbors: int,
    random_state: int,
) -> np.ndarray:
    """Spectral clustering via sklearn."""
    from sklearn.cluster import SpectralClustering

    sc = SpectralClustering(
        n_clusters=n_clusters,
        n_neighbors=n_neighbors,
        random_state=random_state,
        affinity="nearest_neighbors",
    )
    return sc.fit_predict(embeddings)


def reduce_dimensions(
    embeddings: np.ndarray | EmbeddingResult,
    method: str = "umap",
    n_components: int = 2,
    n_neighbors: int = 15,
    min_dist: float = 0.1,
    random_state: int = 42,
) -> np.ndarray:
    """Reduce embedding dimensionality for visualization.

    Parameters
    ----------
    embeddings : np.ndarray or EmbeddingResult
        Patch embeddings, shape (N, D).
    method : str
        Reduction method: "umap", "pca", "tsne".
    n_components : int
        Number of output dimensions.
    n_neighbors : int
        UMAP neighbors parameter.
    min_dist : float
        UMAP minimum distance parameter.
    random_state : int
        Random seed.

    Returns
    -------
    np.ndarray
        Reduced embeddings, shape (N, n_components).
    """
    if isinstance(embeddings, EmbeddingResult):
        embeddings = embeddings.patch_embeddings

    if method == "umap":
        import umap

        reducer = umap.UMAP(
            n_components=n_components,
            n_neighbors=n_neighbors,
            min_dist=min_dist,
            random_state=random_state,
        )
        return reducer.fit_transform(embeddings)

    elif method == "pca":
        from sklearn.decomposition import PCA

        pca = PCA(n_components=n_components, random_state=random_state)
        return pca.fit_transform(embeddings)

    elif method == "tsne":
        from sklearn.manifold import TSNE

        tsne = TSNE(
            n_components=n_components,
            random_state=random_state,
            perplexity=min(30, len(embeddings) - 1),
        )
        return tsne.fit_transform(embeddings)

    else:
        raise ValueError(f"Unknown reduction method: {method}")


def visualize_embeddings(
    result: EmbeddingResult,
    labels: np.ndarray | None = None,
    method: str = "umap",
    show_spatial: bool = True,
    figsize: tuple[int, int] = (16, 6),
    save_path: str | Path | None = None,
    **kwargs: Any,
) -> Any:
    """Visualize KRONOS embeddings with optional spatial overlay.

    Creates a figure with:
    - Left: UMAP/PCA scatter colored by cluster labels
    - Right: Spatial map showing cluster assignments at patch locations

    Parameters
    ----------
    result : EmbeddingResult
        Embedding extraction results.
    labels : np.ndarray, optional
        Cluster labels. If None, clusters with Leiden.
    method : str
        Dimensionality reduction method.
    show_spatial : bool
        Show spatial map alongside embedding plot.
    figsize : tuple
        Figure size.
    save_path : Path, optional
        Save figure to this path.
    **kwargs
        Passed to reduce_dimensions().

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if labels is None:
        labels = cluster_embeddings(result)

    reduced = reduce_dimensions(result, method=method, **kwargs)
    n_clusters = len(np.unique(labels))
    cmap = plt.cm.get_cmap("tab20", n_clusters)

    if show_spatial and result.patch_coords is not None:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=figsize)
    else:
        fig, ax1 = plt.subplots(1, 1, figsize=(figsize[0] // 2, figsize[1]))
        ax2 = None

    # Embedding scatter
    scatter = ax1.scatter(
        reduced[:, 0],
        reduced[:, 1],
        c=labels,
        cmap=cmap,
        s=5,
        alpha=0.7,
    )
    ax1.set_title(f"KRONOS Embeddings ({method.upper()})")
    ax1.set_xlabel(f"{method.upper()} 1")
    ax1.set_ylabel(f"{method.upper()} 2")
    plt.colorbar(scatter, ax=ax1, label="Cluster")

    # Spatial map
    if ax2 is not None:
        coords = result.patch_coords
        ax2.scatter(
            coords[:, 1],  # col = x
            coords[:, 0],  # row = y (invert for image coords)
            c=labels,
            cmap=cmap,
            s=20,
            marker="s",
            alpha=0.8,
        )
        ax2.set_title("Spatial Cluster Map")
        ax2.set_xlabel("X (pixels)")
        ax2.set_ylabel("Y (pixels)")
        ax2.invert_yaxis()
        ax2.set_aspect("equal")

    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info("Saved visualization to %s", save_path)

    return fig


def spatial_search(
    query_embedding: np.ndarray,
    database: EmbeddingResult | list[EmbeddingResult],
    top_k: int = 10,
    metric: str = "cosine",
) -> list[dict[str, Any]]:
    """Find patches most similar to a query embedding.

    Parameters
    ----------
    query_embedding : np.ndarray
        Query embedding vector, shape (D,) or (1, D).
    database : EmbeddingResult or list of EmbeddingResult
        Database of embeddings to search.
    top_k : int
        Number of top matches to return.
    metric : str
        Distance metric: "cosine", "euclidean".

    Returns
    -------
    list of dict
        Top-k matches with keys: index, distance, coords, dataset_idx.
    """
    if query_embedding.ndim == 1:
        query_embedding = query_embedding[None, :]

    # Build database matrix
    if isinstance(database, EmbeddingResult):
        database = [database]

    all_embeddings = []
    all_coords = []
    all_dataset_idx = []

    for ds_idx, result in enumerate(database):
        all_embeddings.append(result.patch_embeddings)
        all_coords.append(result.patch_coords)
        all_dataset_idx.extend([ds_idx] * result.n_patches)

    db_matrix = np.concatenate(all_embeddings, axis=0)
    all_coords = np.concatenate(all_coords, axis=0)

    # Compute distances
    if metric == "cosine":
        # Normalize
        query_norm = query_embedding / (np.linalg.norm(query_embedding, axis=1, keepdims=True) + 1e-8)
        db_norm = db_matrix / (np.linalg.norm(db_matrix, axis=1, keepdims=True) + 1e-8)
        similarities = (query_norm @ db_norm.T).squeeze()
        distances = 1.0 - similarities
    elif metric == "euclidean":
        diff = db_matrix - query_embedding
        distances = np.linalg.norm(diff, axis=1)
    else:
        raise ValueError(f"Unknown metric: {metric}")

    # Get top-k
    top_indices = np.argsort(distances)[:top_k]

    results = []
    for idx in top_indices:
        results.append({
            "index": int(idx),
            "distance": float(distances[idx]),
            "coords": all_coords[idx].tolist(),
            "dataset_idx": all_dataset_idx[idx],
        })

    return results


def compare_datasets(
    results: list[EmbeddingResult],
    dataset_names: list[str] | None = None,
    method: str = "umap",
    figsize: tuple[int, int] = (12, 8),
    save_path: str | Path | None = None,
) -> Any:
    """Compare embeddings across multiple datasets.

    Creates a UMAP visualization colored by dataset origin,
    useful for cross-dataset patient stratification.

    Parameters
    ----------
    results : list of EmbeddingResult
        Embeddings from multiple datasets.
    dataset_names : list of str, optional
        Names for each dataset.
    method : str
        Reduction method.
    figsize : tuple
        Figure size.
    save_path : Path, optional
        Save figure to this path.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if dataset_names is None:
        dataset_names = [f"Dataset {i}" for i in range(len(results))]

    # Concatenate all embeddings
    all_embeddings = np.concatenate([r.patch_embeddings for r in results], axis=0)
    dataset_labels = np.concatenate(
        [np.full(r.n_patches, i) for i, r in enumerate(results)]
    )

    # Reduce dimensions
    reduced = reduce_dimensions(all_embeddings, method=method)

    fig, ax = plt.subplots(figsize=figsize)
    cmap = plt.cm.get_cmap("tab10", len(results))

    for i, name in enumerate(dataset_names):
        mask = dataset_labels == i
        ax.scatter(
            reduced[mask, 0],
            reduced[mask, 1],
            c=[cmap(i)],
            label=name,
            s=3,
            alpha=0.5,
        )

    ax.legend(markerscale=5, fontsize=8)
    ax.set_title(f"Cross-Dataset KRONOS Embeddings ({method.upper()})")
    ax.set_xlabel(f"{method.upper()} 1")
    ax.set_ylabel(f"{method.upper()} 2")
    plt.tight_layout()

    if save_path:
        fig.savefig(save_path, dpi=150, bbox_inches="tight")
        logger.info("Saved comparison to %s", save_path)

    return fig
