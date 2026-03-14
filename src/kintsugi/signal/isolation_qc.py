"""
Signal Isolation QC Visualization

Generates multi-page QC reports for batch signal isolation results.
Three-column layout: Before (self-normalized), After (self-normalized),
Difference (inferno colormap).
"""

from __future__ import annotations

import json
import logging
from pathlib import Path

import numpy as np

logger = logging.getLogger("kintsugi.signal.isolation_qc")


def _self_normalize(image: np.ndarray, p_low: float = 1.0, p_high: float = 99.0) -> np.ndarray:
    """Normalize image to [0, 1] using its own percentile range.

    Parameters
    ----------
    image : np.ndarray
        Input image (any dtype).
    p_low, p_high : float
        Percentiles for clipping range.

    Returns
    -------
    np.ndarray
        Float64 image in [0, 1].
    """
    img_f = image.astype(np.float64)
    vmin = np.percentile(img_f, p_low)
    vmax = np.percentile(img_f, p_high)
    if vmax <= vmin:
        vmax = vmin + 1
    normalized = np.clip((img_f - vmin) / (vmax - vmin), 0, 1)
    return normalized


def _matched_normalize(
    before: np.ndarray,
    after: np.ndarray,
    p_low: float = 1.0,
    p_high: float = 99.0,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize both images to the Before image's percentile range.

    Prevents artificial amplification of subtle tile artifacts in the After
    image that occurs when self-normalization stretches 1-5% intensity
    differences across the full contrast range.

    Parameters
    ----------
    before, after : np.ndarray
        Before and after images (any dtype).
    p_low, p_high : float
        Percentiles for clipping range, computed from ``before``.

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        (before_norm, after_norm), both float64 in [0, 1].
    """
    before_f = before.astype(np.float64)
    after_f = after.astype(np.float64)

    vmin = np.percentile(before_f, p_low)
    vmax = np.percentile(before_f, p_high)
    if vmax <= vmin:
        vmax = vmin + 1

    before_norm = np.clip((before_f - vmin) / (vmax - vmin), 0, 1)
    after_norm = np.clip((after_f - vmin) / (vmax - vmin), 0, 1)
    return before_norm, after_norm


def _downsample(image: np.ndarray, factor: int = 16) -> np.ndarray:
    """Downsample image by block averaging.

    Parameters
    ----------
    image : np.ndarray
        2D input image.
    factor : int
        Downsample factor.

    Returns
    -------
    np.ndarray
        Downsampled image.
    """
    if factor <= 1:
        return image
    h, w = image.shape[:2]
    # Trim to multiple of factor
    h_trim = (h // factor) * factor
    w_trim = (w // factor) * factor
    trimmed = image[:h_trim, :w_trim].astype(np.float64)
    # Block average
    reshaped = trimmed.reshape(h_trim // factor, factor, w_trim // factor, factor)
    return reshaped.mean(axis=(1, 3))


def generate_qc_pages(
    project_dir: str | Path,
    page_size: int = 6,
    dpi: int = 120,
    downsample: int = 16,
    normalize_mode: str = "independent",
) -> list[Path]:
    """Generate QC pages for signal isolation results.

    Three-column layout per channel row:
    - Before (registered): normalized (gray)
    - After (signal_isolated): normalized (gray)
    - Difference (before - after): self-normalized (inferno)

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory.
    page_size : int
        Number of channels per page. Default: 6.
    dpi : int
        Output DPI. Default: 120.
    downsample : int
        Downsample factor for thumbnails. Default: 16.
    normalize_mode : str
        ``"independent"`` (default): each image uses its own p1-p99 range.
        ``"matched"``: both Before and After use the Before image's range,
        preventing artificial amplification of tile-boundary artifacts.

    Returns
    -------
    list[Path]
        Paths to generated QC page images.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import tifffile

    project_dir = Path(project_dir).resolve()
    registered_dir = project_dir / "data" / "processed" / "registered"
    isolated_dir = project_dir / "data" / "processed" / "signal_isolated"
    qc_dir = project_dir / "qc_plots"
    qc_dir.mkdir(parents=True, exist_ok=True)

    # Load manifest for per-channel metadata
    manifest_path = isolated_dir / "signal_isolation_manifest.json"
    manifest = {}
    if manifest_path.exists():
        with open(manifest_path) as f:
            manifest = json.load(f)

    channel_info = manifest.get("channels", {})

    # Discover processed channels
    isolated_files = sorted(isolated_dir.glob("*.tif"))
    if not isolated_files:
        logger.warning("No signal-isolated TIFFs found")
        return []

    # Build channel list with before/after paths
    channels = []
    for iso_path in isolated_files:
        marker = iso_path.stem
        info = channel_info.get(marker, {})
        cycle = info.get("cycle")

        # Find the before image in registered dir
        before_path = None
        if cycle:
            before_path = registered_dir / f"cyc{int(cycle):02d}" / f"{marker}.tif"
        else:
            # Search all cycle dirs
            for cyc_dir in sorted(registered_dir.glob("cyc*")):
                candidate = cyc_dir / f"{marker}.tif"
                if candidate.exists():
                    before_path = candidate
                    break

        if before_path is None or not before_path.exists():
            logger.warning(f"Before image not found for {marker}, skipping QC")
            continue

        channels.append(
            {
                "marker": marker,
                "before_path": before_path,
                "after_path": iso_path,
                "info": info,
            }
        )

    if not channels:
        return []

    # Generate pages
    pages = []
    n_pages = (len(channels) + page_size - 1) // page_size

    for page_idx in range(n_pages):
        start = page_idx * page_size
        end = min(start + page_size, len(channels))
        page_channels = channels[start:end]
        n_rows = len(page_channels)

        fig, axes = plt.subplots(n_rows, 3, figsize=(15, 3.5 * n_rows))
        if n_rows == 1:
            axes = axes[np.newaxis, :]

        for row, ch in enumerate(page_channels):
            before = tifffile.imread(str(ch["before_path"]))
            after = tifffile.imread(str(ch["after_path"]))

            # Downsample for display
            before_ds = _downsample(before, downsample)
            after_ds = _downsample(after, downsample)

            # Normalize
            if normalize_mode == "matched":
                before_norm, after_norm = _matched_normalize(
                    before_ds, after_ds
                )
            else:
                before_norm = _self_normalize(before_ds)
                after_norm = _self_normalize(after_ds)

            # Difference (what was removed)
            diff = before_ds.astype(np.float64) - after_ds.astype(np.float64)
            diff = np.clip(diff, 0, None)
            diff_norm = _self_normalize(diff) if np.max(diff) > 0 else diff

            # Plot
            axes[row, 0].imshow(before_norm, cmap="gray", vmin=0, vmax=1)
            axes[row, 0].set_title("Before" if row == 0 else "")
            axes[row, 0].axis("off")

            axes[row, 1].imshow(after_norm, cmap="gray", vmin=0, vmax=1)
            axes[row, 1].set_title("After" if row == 0 else "")
            axes[row, 1].axis("off")

            axes[row, 2].imshow(diff_norm, cmap="inferno", vmin=0, vmax=1)
            axes[row, 2].set_title("Removed" if row == 0 else "")
            axes[row, 2].axis("off")

            # Annotation
            info = ch["info"]
            method = info.get("method", "?")
            scale = info.get("parameters", {}).get(
                "blank_scale_factor",
                info.get("parameters", {}).get("base_scale_factor", "?"),
            )
            zero_pct = info.get("zero_percent", "?")
            q_score = info.get("quality_metrics", {}).get("quality_score", "?")
            resid = info.get("quality_metrics", {}).get("residual_correlation", "?")
            warning = info.get("warning", "")
            flag = f" [{warning}]" if warning else ""

            label = (
                f"{ch['marker']} | {method} | scale={scale} | "
                f"zero={zero_pct}% | Q={q_score} | resid={resid}{flag}"
            )
            axes[row, 0].set_ylabel(label, fontsize=7, rotation=0, labelpad=220, va="center")

        plt.tight_layout()
        page_path = qc_dir / f"signal_isolation_qc_p{page_idx + 1}.png"
        fig.savefig(page_path, dpi=dpi, bbox_inches="tight")
        plt.close(fig)
        pages.append(page_path)
        logger.info(f"QC page {page_idx + 1}/{n_pages}: {page_path}")

    return pages


def generate_summary_table(
    project_dir: str | Path,
) -> str:
    """Generate a Rich-formatted summary table from the manifest.

    Parameters
    ----------
    project_dir : str or Path
        Path to KINTSUGI project directory.

    Returns
    -------
    str
        Rendered summary table string.
    """
    from rich.console import Console
    from rich.table import Table

    project_dir = Path(project_dir).resolve()
    manifest_path = (
        project_dir / "data" / "processed" / "signal_isolated" / "signal_isolation_manifest.json"
    )

    if not manifest_path.exists():
        return "No manifest found. Run 'kintsugi workflow isolate run' first."

    with open(manifest_path) as f:
        manifest = json.load(f)

    table = Table(title="Signal Isolation Summary")
    table.add_column("Marker", style="cyan")
    table.add_column("Cycle", justify="right")
    table.add_column("Blank", style="dim")
    table.add_column("Method", style="bold")
    table.add_column("Scale", justify="right")
    table.add_column("B/S Ratio", justify="right", style="dim")
    table.add_column("Corr", justify="right")
    table.add_column("Q Score", justify="right")
    table.add_column("Zero%", justify="right")
    table.add_column("Flags", style="yellow", max_width=30)

    for marker, info in sorted(manifest.get("channels", {}).items()):
        if info.get("status") == "skipped":
            continue

        params = info.get("parameters", {})
        analysis = info.get("analysis", {})
        quality = info.get("quality_metrics", {})

        scale = params.get("blank_scale_factor", params.get("base_scale_factor", ""))
        signal_p99 = analysis.get("signal_p99", 0)
        blank_p99 = analysis.get("blank_p99", 0)
        ratio = f"{blank_p99 / signal_p99:.2f}" if signal_p99 > 0 else "-"

        corr = analysis.get("correlation", "")
        q_score = quality.get("quality_score", "")
        zero_pct = info.get("zero_percent", "")
        warning = info.get("warning", "")

        # Color-code quality score
        if isinstance(q_score, (int, float)):
            if q_score >= 0.7:
                q_str = f"[green]{q_score:.3f}[/green]"
            elif q_score >= 0.5:
                q_str = f"[yellow]{q_score:.3f}[/yellow]"
            else:
                q_str = f"[red]{q_score:.3f}[/red]"
        else:
            q_str = str(q_score)

        # Color-code method
        method = info.get("method", "")
        if method == "weighted":
            method_str = f"[magenta]{method}[/magenta]"
        else:
            method_str = method

        table.add_row(
            marker,
            str(info.get("cycle", "")),
            info.get("blank_used", ""),
            method_str,
            f"{scale}" if isinstance(scale, str) else f"{scale:.2f}",
            ratio,
            f"{corr:.3f}" if isinstance(corr, (int, float)) else str(corr),
            q_str,
            f"{zero_pct:.1f}" if isinstance(zero_pct, (int, float)) else str(zero_pct),
            warning,
        )

    # Summary row
    summary = manifest.get("summary", {})
    if summary:
        table.add_section()
        table.add_row(
            f"[bold]Total: {summary.get('total', 0)}[/bold]",
            "",
            "",
            f"G:{summary.get('global', 0)} W:{summary.get('weighted', 0)}",
            "",
            "",
            "",
            f"avg={summary.get('mean_quality', 0):.3f}",
            "",
            f"skip={summary.get('skipped', 0)} err={summary.get('error', 0)}",
        )

    console = Console(record=True)
    console.print(table)
    return console.export_text()
