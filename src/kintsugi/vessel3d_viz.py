"""
3D Vessel Segmentation Visualization Utilities

Separated from vessel3d.py to keep the core module importable without
napari or heavy visualization dependencies.

Usage:
    from kintsugi.vessel3d_viz import (
        ortho_view,
        overlay_mask_on_raw,
        render_skeleton_mip,
        plot_vessel_features,
    )
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:
    import matplotlib.figure
    import pandas as pd

logger = logging.getLogger(__name__)


def ortho_view(
    volume: np.ndarray,
    z: int | None = None,
    y: int | None = None,
    x: int | None = None,
    title: str = "",
    cmap: str = "gray",
    figsize: tuple[int, int] = (14, 5),
    vmin: float | None = None,
    vmax: float | None = None,
) -> "matplotlib.figure.Figure":
    """Three-panel orthogonal slice viewer (XY, XZ, YZ).

    Parameters
    ----------
    volume : np.ndarray
        3D array (z, y, x).
    z, y, x : int, optional
        Slice indices. Default: center of each axis.
    title : str
        Figure title.
    cmap : str
        Matplotlib colormap.
    figsize : tuple
        Figure size.
    vmin, vmax : float, optional
        Display range. Default: 1st-99th percentile.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object.
    """
    import matplotlib.pyplot as plt

    nz, ny, nx = volume.shape
    if z is None:
        z = nz // 2
    if y is None:
        y = ny // 2
    if x is None:
        x = nx // 2

    if vmin is None or vmax is None:
        p1, p99 = np.percentile(volume, [1, 99])
        vmin = vmin if vmin is not None else p1
        vmax = vmax if vmax is not None else p99

    fig, axes = plt.subplots(1, 3, figsize=figsize)

    # XY plane (axial)
    axes[0].imshow(volume[z], cmap=cmap, vmin=vmin, vmax=vmax, aspect="equal")
    axes[0].set_title(f"XY (z={z})")
    axes[0].axhline(y, color="cyan", linewidth=0.5, alpha=0.7)
    axes[0].axvline(x, color="cyan", linewidth=0.5, alpha=0.7)

    # XZ plane (coronal)
    axes[1].imshow(volume[:, y, :], cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    axes[1].set_title(f"XZ (y={y})")
    axes[1].axhline(z, color="cyan", linewidth=0.5, alpha=0.7)
    axes[1].axvline(x, color="cyan", linewidth=0.5, alpha=0.7)

    # YZ plane (sagittal)
    axes[2].imshow(volume[:, :, x], cmap=cmap, vmin=vmin, vmax=vmax, aspect="auto")
    axes[2].set_title(f"YZ (x={x})")
    axes[2].axhline(z, color="cyan", linewidth=0.5, alpha=0.7)
    axes[2].axvline(y, color="cyan", linewidth=0.5, alpha=0.7)

    for ax in axes:
        ax.axis("off")

    if title:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    return fig


def overlay_mask_on_raw(
    raw: np.ndarray,
    mask: np.ndarray,
    z: int | None = None,
    alpha: float = 0.3,
    mask_color: tuple[float, float, float] = (1.0, 0.2, 0.2),
    figsize: tuple[int, int] = (12, 6),
    title: str = "",
) -> "matplotlib.figure.Figure":
    """Side-by-side comparison of raw data and mask overlay.

    Parameters
    ----------
    raw : np.ndarray
        3D raw/vesselness volume.
    mask : np.ndarray
        3D binary mask.
    z : int, optional
        Z-plane to display. Default: center.
    alpha : float
        Mask overlay transparency.
    mask_color : tuple
        RGB color for mask overlay.
    figsize : tuple
        Figure size.
    title : str
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    if z is None:
        z = raw.shape[0] // 2

    raw_slice = raw[z].astype(np.float32)
    mask_slice = mask[z].astype(bool)

    # Normalize raw to [0, 1]
    p1, p99 = np.percentile(raw_slice, [1, 99])
    raw_norm = np.clip((raw_slice - p1) / max(p99 - p1, 1e-6), 0, 1)

    # Create RGB overlay
    overlay = np.stack([raw_norm] * 3, axis=-1)
    for c in range(3):
        overlay[mask_slice, c] = (
            overlay[mask_slice, c] * (1 - alpha) + mask_color[c] * alpha
        )

    fig, axes = plt.subplots(1, 2, figsize=figsize)

    axes[0].imshow(raw_norm, cmap="gray")
    axes[0].set_title(f"Raw (z={z})")
    axes[0].axis("off")

    axes[1].imshow(overlay)
    axes[1].set_title(f"Mask overlay (z={z})")
    axes[1].axis("off")

    if title:
        fig.suptitle(title, fontsize=14)
    fig.tight_layout()
    return fig


def render_skeleton_mip(
    skeleton: np.ndarray,
    binary_mask: np.ndarray | None = None,
    raw: np.ndarray | None = None,
    figsize: tuple[int, int] = (12, 6),
    title: str = "",
) -> "matplotlib.figure.Figure":
    """Max-intensity projection with skeleton overlay.

    Parameters
    ----------
    skeleton : np.ndarray
        3D skeleton image (bool).
    binary_mask : np.ndarray, optional
        3D binary mask for background context.
    raw : np.ndarray, optional
        3D raw volume for background. If provided, overrides binary_mask.
    figsize : tuple
        Figure size.
    title : str
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    skel_mip = skeleton.max(axis=0).astype(bool)

    if raw is not None:
        bg = raw.max(axis=0).astype(np.float32)
        p1, p99 = np.percentile(bg, [1, 99])
        bg = np.clip((bg - p1) / max(p99 - p1, 1e-6), 0, 1)
    elif binary_mask is not None:
        bg = binary_mask.max(axis=0).astype(np.float32) * 0.3
    else:
        bg = np.zeros(skel_mip.shape, dtype=np.float32)

    # RGB: background in gray, skeleton in red
    rgb = np.stack([bg, bg, bg], axis=-1)
    rgb[skel_mip, 0] = 1.0
    rgb[skel_mip, 1] = 0.1
    rgb[skel_mip, 2] = 0.1

    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.imshow(rgb)
    ax.set_title(title or "Skeleton MIP (red) on mask (gray)")
    ax.axis("off")
    fig.tight_layout()
    return fig


def plot_vessel_features(
    df: "pd.DataFrame",
    figsize: tuple[int, int] = (16, 10),
    title: str = "",
) -> "matplotlib.figure.Figure":
    """Histogram panel of vessel morphometric features.

    Parameters
    ----------
    df : pd.DataFrame
        Vessel features from analyze_vessel_graph().
    figsize : tuple
        Figure size.
    title : str
        Figure title.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 3, figsize=figsize)

    # 1. Diameter distribution
    if "mean_diameter_um" in df.columns:
        ax = axes[0, 0]
        diameters = df["mean_diameter_um"].dropna()
        ax.hist(diameters, bins=50, color="steelblue", edgecolor="white", alpha=0.8)
        ax.set_xlabel("Diameter (um)")
        ax.set_ylabel("Count")
        ax.set_title(f"Vessel Diameter (n={len(diameters)})")
        ax.axvline(diameters.median(), color="red", linestyle="--", label=f"median={diameters.median():.1f}")
        ax.legend(fontsize=8)

    # 2. Branch length distribution
    if "branch_length_um" in df.columns:
        ax = axes[0, 1]
        lengths = df["branch_length_um"].dropna()
        ax.hist(lengths, bins=50, color="forestgreen", edgecolor="white", alpha=0.8)
        ax.set_xlabel("Branch Length (um)")
        ax.set_ylabel("Count")
        ax.set_title(f"Branch Length (n={len(lengths)})")
        ax.axvline(lengths.median(), color="red", linestyle="--", label=f"median={lengths.median():.1f}")
        ax.legend(fontsize=8)

    # 3. Tortuosity distribution
    if "tortuosity" in df.columns:
        ax = axes[0, 2]
        tort = df["tortuosity"].dropna()
        tort = tort[tort < tort.quantile(0.99)]  # Clip outliers for visualization
        ax.hist(tort, bins=50, color="darkorange", edgecolor="white", alpha=0.8)
        ax.set_xlabel("Tortuosity")
        ax.set_ylabel("Count")
        ax.set_title(f"Tortuosity (n={len(tort)})")
        ax.axvline(tort.median(), color="red", linestyle="--", label=f"median={tort.median():.2f}")
        ax.legend(fontsize=8)

    # 4. Branch type distribution
    if "branch_type" in df.columns:
        ax = axes[1, 0]
        types = df["branch_type"].value_counts().sort_index()
        type_labels = {0: "EP-EP", 1: "JP-EP", 2: "JP-JP", 3: "Isolated"}
        labels = [type_labels.get(t, str(t)) for t in types.index]
        ax.bar(labels, types.values, color="mediumpurple", edgecolor="white")
        ax.set_xlabel("Branch Type")
        ax.set_ylabel("Count")
        ax.set_title("Branch Types")

    # 5. Diameter vs length scatter
    if "mean_diameter_um" in df.columns and "branch_length_um" in df.columns:
        ax = axes[1, 1]
        ax.scatter(
            df["branch_length_um"], df["mean_diameter_um"],
            s=3, alpha=0.3, color="steelblue",
        )
        ax.set_xlabel("Branch Length (um)")
        ax.set_ylabel("Diameter (um)")
        ax.set_title("Diameter vs Length")

    # 6. Summary statistics text
    ax = axes[1, 2]
    ax.axis("off")
    stats_text = (
        f"Total segments: {len(df)}\n"
        f"Total network length: {df.get('branch_length_um', [0]).sum():.0f} um\n"
    )
    if "mean_diameter_um" in df.columns:
        stats_text += f"Mean diameter: {df['mean_diameter_um'].mean():.1f} um\n"
        stats_text += f"Median diameter: {df['mean_diameter_um'].median():.1f} um\n"
    if "branch_length_um" in df.columns:
        stats_text += f"Mean branch length: {df['branch_length_um'].mean():.1f} um\n"
    if "tortuosity" in df.columns:
        stats_text += f"Mean tortuosity: {df['tortuosity'].mean():.2f}\n"
    ax.text(
        0.1, 0.5, stats_text, transform=ax.transAxes,
        fontsize=12, verticalalignment="center", fontfamily="monospace",
    )
    ax.set_title("Summary Statistics")

    if title:
        fig.suptitle(title, fontsize=14, fontweight="bold")
    fig.tight_layout()
    return fig
