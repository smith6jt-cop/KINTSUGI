"""
KINTSUGI Analysis Utilities Module

Extracted from 4_Segmentation_Analysis.ipynb to eliminate cell execution-order
fragility. Functions defined inline in notebook cells are moved here for
importable, testable access.

Location: notebooks/Kanalysis.py

Functions:
---------
Clustering:
    - visualize_som: SOM clustering with progressive visualization
    - load_marker_names: Load marker names from CHANNELNAMES.txt or project config

Marker Configuration:
    - get_cell_markers: Get markers suitable for cell/cytoplasm compartment
    - get_nuclear_markers: Get markers suitable for nuclear compartment
    - get_ecm_markers: Get markers suitable for ECM compartment
    - load_channel_names: Load channel names from project metadata

Checkpointing:
    - save_checkpoint: Save intermediate DataFrame to Parquet
    - load_checkpoint: Load intermediate DataFrame from Parquet if exists
    - checkpoint_exists: Check if a checkpoint file exists
"""

from __future__ import annotations

import os
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np
import pandas as pd

if TYPE_CHECKING:
    from kintsugi.project import KintsugiProject

# Lazy imports for optional heavy dependencies
_pyFlowSOM = None
_StandardScaler = None


def _get_pyFlowSOM():
    global _pyFlowSOM
    if _pyFlowSOM is None:
        import pyFlowSOM

        _pyFlowSOM = pyFlowSOM
    return _pyFlowSOM


def _get_StandardScaler():
    global _StandardScaler
    if _StandardScaler is None:
        from sklearn.preprocessing import StandardScaler

        _StandardScaler = StandardScaler
    return _StandardScaler


# =============================================================================
# Marker Configuration
# =============================================================================

# Default markers excluded from each compartment analysis.
# These are guidance defaults; actual panels vary by experiment.
_NUCLEAR_EXCLUDE = {"CollagenIV", "ECAD"}
_CELL_EXCLUDE = {"DAPI", "FoxP3", "Ki67", "CollagenIV", "ECAD"}
_ECM_EXCLUDE = {"DAPI", "CD3e", "CD15", "CD20", "CD45", "CD45RO", "FoxP3", "Ki67"}


def load_channel_names(project_dir: str | Path) -> list[str]:
    """
    Load channel/marker names from CHANNELNAMES.txt in the project metadata.

    Supports two formats:
    - Simple list (one name per line, CODEX format)
    - Cycle-prefixed format ("1: DAPI, Blank, Blank, Blank")

    Blank and DAPI-## entries (cycle markers) are included but can be
    filtered downstream.

    Parameters
    ----------
    project_dir : str or Path
        Path to the KINTSUGI project directory.

    Returns
    -------
    list[str]
        Ordered list of unique marker names (excluding blanks and
        cycle-specific DAPI duplicates like DAPI-02, DAPI-03, etc.).
    """
    project_dir = Path(project_dir)
    channel_file = project_dir / "meta" / "CHANNELNAMES.txt"

    if not channel_file.exists():
        raise FileNotFoundError(
            f"CHANNELNAMES.txt not found at {channel_file}. "
            "Create this file with one marker name per line, or use "
            "cycle-prefixed format (e.g., '1: DAPI, CD3, CD20, Blank')."
        )

    all_names = []
    with open(channel_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            # Cycle-prefixed format: "1: DAPI, Blank, Blank, Blank"
            if ":" in line and line.split(":")[0].strip().isdigit():
                names = [n.strip() for n in line.split(":", 1)[1].split(",")]
                all_names.extend(names)
            else:
                all_names.append(line)

    # Deduplicate while preserving order, skip blanks and cycle DAPI duplicates
    seen = set()
    unique_markers = []
    for name in all_names:
        # Skip blanks
        if name.lower() == "blank" or not name:
            continue
        # Normalize cycle-specific DAPI (DAPI-01, DAPI-02) to just "DAPI"
        if name.upper().startswith("DAPI-") or name.upper().startswith("DAPI_"):
            name = "DAPI"
        if name not in seen:
            seen.add(name)
            unique_markers.append(name)

    return unique_markers


def get_cell_markers(
    all_markers: list[str],
    exclude: set[str] | None = None,
) -> list[str]:
    """
    Get markers suitable for cell/cytoplasm compartment analysis.

    Excludes nuclear-specific (DAPI, FoxP3, Ki67) and structural
    (CollagenIV, ECAD) markers by default.

    Parameters
    ----------
    all_markers : list[str]
        Full list of marker names from the panel.
    exclude : set[str], optional
        Custom set of markers to exclude. If None, uses default exclusions.

    Returns
    -------
    list[str]
        Filtered marker list for cell compartment.
    """
    if exclude is None:
        exclude = _CELL_EXCLUDE
    return [m for m in all_markers if m not in exclude]


def get_nuclear_markers(
    all_markers: list[str],
    exclude: set[str] | None = None,
) -> list[str]:
    """
    Get markers suitable for nuclear compartment analysis.

    Excludes extracellular-only markers (CollagenIV, ECAD) by default.
    Includes DAPI, FoxP3, Ki67.

    Parameters
    ----------
    all_markers : list[str]
        Full list of marker names from the panel.
    exclude : set[str], optional
        Custom set of markers to exclude. If None, uses default exclusions.

    Returns
    -------
    list[str]
        Filtered marker list for nuclear compartment.
    """
    if exclude is None:
        exclude = _NUCLEAR_EXCLUDE
    return [m for m in all_markers if m not in exclude]


def get_ecm_markers(
    all_markers: list[str],
    exclude: set[str] | None = None,
) -> list[str]:
    """
    Get markers suitable for ECM compartment analysis.

    Excludes nuclear/intracellular markers by default.

    Parameters
    ----------
    all_markers : list[str]
        Full list of marker names from the panel.
    exclude : set[str], optional
        Custom set of markers to exclude. If None, uses default exclusions.

    Returns
    -------
    list[str]
        Filtered marker list for ECM compartment.
    """
    if exclude is None:
        exclude = _ECM_EXCLUDE
    return [m for m in all_markers if m not in exclude]


# =============================================================================
# Checkpointing
# =============================================================================


def save_checkpoint(
    data: pd.DataFrame,
    checkpoint_dir: str | Path,
    name: str,
    format: str = "parquet",
) -> Path:
    """
    Save intermediate DataFrame as a checkpoint for crash recovery.

    Parameters
    ----------
    data : pd.DataFrame
        Data to checkpoint.
    checkpoint_dir : str or Path
        Directory to save checkpoints.
    name : str
        Checkpoint name (e.g., 'cell_features', 'som_clusters', 'qc_filtered').
    format : str
        Save format: 'parquet' (default, recommended) or 'csv'.

    Returns
    -------
    Path
        Path to the saved checkpoint file.
    """
    checkpoint_dir = Path(checkpoint_dir)
    checkpoint_dir.mkdir(parents=True, exist_ok=True)

    if format == "parquet":
        path = checkpoint_dir / f"{name}.parquet"
        data.to_parquet(path, index=False)
    elif format == "csv":
        path = checkpoint_dir / f"{name}.csv"
        data.to_csv(path, index=False)
    else:
        raise ValueError(f"Unsupported format: {format}. Use 'parquet' or 'csv'.")

    return path


def load_checkpoint(
    checkpoint_dir: str | Path,
    name: str,
    format: str = "parquet",
) -> pd.DataFrame | None:
    """
    Load a checkpoint if it exists, otherwise return None.

    Parameters
    ----------
    checkpoint_dir : str or Path
        Directory containing checkpoints.
    name : str
        Checkpoint name.
    format : str
        Expected format: 'parquet' or 'csv'.

    Returns
    -------
    pd.DataFrame or None
        Loaded data, or None if checkpoint doesn't exist.
    """
    checkpoint_dir = Path(checkpoint_dir)

    if format == "parquet":
        path = checkpoint_dir / f"{name}.parquet"
        if path.exists():
            return pd.read_parquet(path)
    elif format == "csv":
        path = checkpoint_dir / f"{name}.csv"
        if path.exists():
            return pd.read_csv(path)

    return None


def checkpoint_exists(
    checkpoint_dir: str | Path,
    name: str,
    format: str = "parquet",
) -> bool:
    """
    Check if a checkpoint file exists.

    Parameters
    ----------
    checkpoint_dir : str or Path
        Directory containing checkpoints.
    name : str
        Checkpoint name.
    format : str
        Expected format.

    Returns
    -------
    bool
        True if checkpoint exists.
    """
    checkpoint_dir = Path(checkpoint_dir)
    ext = "parquet" if format == "parquet" else "csv"
    return (checkpoint_dir / f"{name}.{ext}").exists()


# =============================================================================
# SOM Clustering Visualization
# =============================================================================


def visualize_som(
    pixel_data: pd.DataFrame,
    channels: list[str],
    image: np.ndarray,
    channel_vis: str,
    num_nodes: int = 10,
    num_passes: int = 5,
    lr_start: float = 0.05,
    lr_end: float = 0.01,
    live_plots: bool = False,
    seed: int = 627,
) -> tuple[np.ndarray, pd.DataFrame, pd.DataFrame]:
    """
    Perform SOM clustering with progressive visualization.

    Runs pyFlowSOM self-organizing map on the selected marker channels
    and generates composite plots showing spatial clusters + marker
    expression heatmaps for each training iteration.

    Parameters
    ----------
    pixel_data : pd.DataFrame
        Cell feature table. Must contain the columns listed in ``channels``
        plus 'centroid_0' and 'centroid_1'.
    channels : list[str]
        Marker column names to use for clustering.
    image : np.ndarray
        2D reference image for spatial context (e.g., one marker channel).
    channel_vis : str
        Name of the channel being visualized (for plot titles).
    num_nodes : int
        Total number of SOM nodes (grid_size = sqrt(num_nodes)).
    num_passes : int
        Number of SOM training iterations.
    lr_start : float
        Initial learning rate.
    lr_end : float
        Final learning rate.
    live_plots : bool
        If True, display plots during each iteration.
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    scatter_stack : np.ndarray
        Stack of rendered plot frames (n_passes, height, width, 3).
    train_data : pd.DataFrame
        Input data with added 'cluster' and 'distances' columns.
    df_mean : pd.DataFrame
        Mean marker expression per cluster.
    """
    import matplotlib.pyplot as plt
    import seaborn as sns
    from matplotlib.colors import ListedColormap
    from tqdm.notebook import tqdm

    pyFlowSOM = _get_pyFlowSOM()
    StandardScaler = _get_StandardScaler()

    scaler = StandardScaler()
    grid_size = int(np.sqrt(num_nodes))

    train_data = pixel_data[channels].copy()
    train_data["centroid_0"] = pixel_data["centroid_0"]
    train_data["centroid_1"] = pixel_data["centroid_1"]

    train_array_pre = pixel_data[channels].values.astype(np.float64)
    train_array = scaler.fit_transform(train_array_pre)

    colors = sns.color_palette("tab20", n_colors=num_nodes)
    node_cmap = ListedColormap(colors)

    pbar = tqdm(
        total=100,
        unit="Percent",
        bar_format="{l_bar}{bar}| {n_fmt}/{total_fmt} [{elapsed}<{remaining}]",
        colour="green",
        position=0,
        leave=True,
    )

    scatter_frames = []
    for rlen_iter in range(1, num_passes + 1):
        som = pyFlowSOM.som(
            train_array,
            grid_size,
            grid_size,
            rlen_iter,
            alpha_range=(lr_start, lr_end),
            seed=seed,
        )
        clusters, dist = pyFlowSOM.map_data_to_nodes(som, train_array)

        unique_clusters = np.unique(clusters)
        cluster_map = {old: new for new, old in enumerate(unique_clusters)}
        clusters = np.array([cluster_map[c] for c in clusters])
        train_data["cluster"] = clusters

        df_mean = train_data.groupby(["cluster"]).mean()
        df_mean = df_mean.drop(columns=["centroid_0", "centroid_1"])
        df_mean.index = [f"Cluster_{i + 1}" for i in range(len(df_mean))]

        fig = plt.figure(figsize=(18, 8))
        gs = fig.add_gridspec(3, 5)

        ax_top_left = fig.add_subplot(gs[0, 0])
        ax_top_right = fig.add_subplot(gs[1:, 0:2])
        ax_bottom = fig.add_subplot(gs[:, 2:])

        im = ax_top_left.imshow(image, vmax=np.percentile(image, 99.9), cmap="turbo")
        fig.colorbar(im, shrink=0.5, ax=ax_top_left, cmap=node_cmap, aspect=30, pad=0.02)
        ax_top_left.set_title(f"Original Image: {channel_vis}")

        im = ax_top_right.scatter(
            pixel_data["centroid_0"],
            pixel_data["centroid_1"],
            c=train_data["cluster"],
            cmap=node_cmap,
            s=1,
        )
        ax_top_right.invert_yaxis()
        ax_top_right.set_aspect(image.shape[1] / image.shape[0])
        ax_top_right.set_title(f"SOM Iteration {rlen_iter}")
        cbar = fig.colorbar(
            im,
            shrink=0.8,
            ax=ax_top_right,
            cmap=node_cmap,
            aspect=30,
            pad=0.02,
            boundaries=np.arange(-0.5, len(unique_clusters) + 0.5, 1),
            ticks=np.arange(0, len(unique_clusters)),
        )
        cbar.set_ticklabels(np.arange(1, len(unique_clusters) + 1))

        g = sns.clustermap(df_mean.T, z_score=1, cmap="vlag", center=0, yticklabels=True)
        heatmap_data = g.data2d
        plt.close(g.figure)
        sns.heatmap(
            data=heatmap_data,
            ax=ax_bottom,
            linewidths=0.05,
            linecolor="grey",
            cmap="coolwarm",
            robust=True,
            xticklabels=g.data2d.columns,
            yticklabels=g.data2d.index,
        )
        ax_bottom.set_title(f"Cluster heatmap for SOM iteration {rlen_iter}")
        plt.setp(ax_bottom.get_xticklabels(), rotation=90)
        plt.setp(ax_bottom.get_yticklabels(), rotation=0)
        fig.tight_layout()

        canvas = fig.canvas
        canvas.draw()
        width, height = fig.get_size_inches() * fig.dpi
        buf = np.frombuffer(canvas.buffer_rgba(), dtype=np.uint8)
        buf.shape = (int(height), int(width), 4)
        scatter_frame = buf[:, :, :3]

        scatter_frames.append(scatter_frame)
        if live_plots:
            plt.show()
        plt.close(fig)
        pbar.update(100 / num_passes)

    train_data["distances"] = dist
    pbar.close()
    scatter_stack = np.stack(scatter_frames, axis=0)

    return scatter_stack, train_data, df_mean
