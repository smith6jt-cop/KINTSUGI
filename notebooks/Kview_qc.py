"""
Plotly-based visualization panel for KINTSUGI pipeline QC.

This module provides interactive viewers with sliders for inspecting images
at each stage of the processing pipeline:
- Raw tiles
- Stitched images
- Deconvolved images
- EDF output

Usage:
    from Kview_qc import (
        view_raw_tiles,
        view_zstack,
        view_edf,
        view_pipeline_comparison
    )

    # View raw tiles with tile/z-plane sliders
    view_raw_tiles(raw_dir, cycle=1, channel=3)

    # View stitched or deconvolved z-stack
    view_zstack(stitched_dir / "cyc01" / "CH3", title="Stitched cyc01 CH3")

    # View EDF output
    view_edf(edf_path)

    # Compare all stages side-by-side
    view_pipeline_comparison(project_dir, cycle=1, channel=3)
"""

from __future__ import annotations

import numpy as np
from pathlib import Path
from typing import Callable, Literal
import warnings

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    HAS_PLOTLY = True
except ImportError:
    HAS_PLOTLY = False
    warnings.warn("Plotly not installed. Install with: pip install plotly")

try:
    import tifffile
    HAS_TIFFFILE = True
except ImportError:
    HAS_TIFFFILE = False


def _check_deps():
    """Check required dependencies."""
    if not HAS_PLOTLY:
        raise ImportError("Plotly required. Install with: pip install plotly")
    if not HAS_TIFFFILE:
        raise ImportError("tifffile required. Install with: pip install tifffile")


def _downsample(img: np.ndarray, max_size: int = 1000) -> np.ndarray:
    """Downsample image for display performance."""
    h, w = img.shape[:2]
    scale = min(1.0, max_size / max(h, w))
    if scale < 1.0:
        step = int(1 / scale)
        return img[::step, ::step]
    return img


def _normalize_for_display(img: np.ndarray, percentile: tuple = (1, 99)) -> np.ndarray:
    """Normalize image to 0-255 for display."""
    img = img.astype(np.float32)
    p_low, p_high = np.percentile(img, percentile)
    img = np.clip((img - p_low) / (p_high - p_low + 1e-10), 0, 1)
    return (img * 255).astype(np.uint8)


def _compute_stats(img: np.ndarray) -> dict:
    """Compute image statistics."""
    return {
        'min': int(img.min()),
        'max': int(img.max()),
        'mean': float(img.mean()),
        'std': float(img.std()),
        'shape': img.shape,
    }


def _stats_text(stats: dict) -> str:
    """Format stats as text for display."""
    return (
        f"Shape: {stats['shape']}<br>"
        f"Min: {stats['min']}, Max: {stats['max']}<br>"
        f"Mean: {stats['mean']:.1f}, Std: {stats['std']:.1f}"
    )


# =============================================================================
# RAW TILE VIEWER
# =============================================================================

def view_raw_tiles(
    raw_dir: str | Path,
    cycle: int,
    channel: int,
    n_rows: int = 13,
    n_cols: int = 9,
    max_display_size: int = 800,
    height: int = 700,
) -> go.Figure:
    """
    Interactive viewer for raw tiles with tile and z-plane sliders.

    Parameters
    ----------
    raw_dir : Path
        Path to raw data directory containing cyc001/, cyc002/, etc.
    cycle : int
        Cycle number (1-based)
    channel : int
        Channel number (1-based)
    n_rows, n_cols : int
        Grid dimensions for tile layout
    max_display_size : int
        Maximum dimension for downsampled display
    height : int
        Figure height in pixels

    Returns
    -------
    go.Figure
        Plotly figure with sliders
    """
    _check_deps()

    raw_dir = Path(raw_dir)
    cycle_dir = raw_dir / f"cyc{cycle:03d}"

    if not cycle_dir.exists():
        raise FileNotFoundError(f"Cycle directory not found: {cycle_dir}")

    # Find all tiles for this channel
    pattern = f"*_CH{channel}.tif"
    all_files = sorted(cycle_dir.glob(pattern))

    if not all_files:
        raise FileNotFoundError(f"No files found matching {pattern} in {cycle_dir}")

    # Parse z-planes and tile indices
    z_planes = sorted(set(int(f.stem.split('_')[2][1:]) for f in all_files))
    n_tiles = len(all_files) // len(z_planes)

    print(f"Found {len(all_files)} files: {n_tiles} tiles × {len(z_planes)} z-planes")

    # Load first image to get dimensions
    sample = tifffile.imread(all_files[0])
    tile_h, tile_w = sample.shape

    # Create figure
    fig = go.Figure()

    # Initial image (first z-plane, first tile)
    def get_tile_path(z_idx: int, tile_idx: int) -> Path:
        z = z_planes[z_idx]
        # Tile index is 1-based in filename
        pattern = f"*_{tile_idx+1:05d}_Z{z:03d}_CH{channel}.tif"
        matches = list(cycle_dir.glob(pattern))
        return matches[0] if matches else None

    def load_tile(z_idx: int, tile_idx: int) -> np.ndarray:
        path = get_tile_path(z_idx, tile_idx)
        if path and path.exists():
            return tifffile.imread(path)
        return np.zeros((tile_h, tile_w), dtype=np.uint16)

    # Load initial tile
    initial_img = load_tile(0, 0)
    img_display = _normalize_for_display(_downsample(initial_img, max_display_size))
    stats = _compute_stats(initial_img)

    # Add heatmap trace
    fig.add_trace(go.Heatmap(
        z=img_display,
        colorscale='gray',
        showscale=False,
        hovertemplate='x: %{x}<br>y: %{y}<br>value: %{z}<extra></extra>'
    ))

    # Create sliders
    z_steps = []
    for i, z in enumerate(z_planes):
        z_steps.append({
            'args': [{'visible': [True]}, {'title': f'Z-plane {z}'}],
            'label': str(z),
            'method': 'update'
        })

    tile_steps = []
    for i in range(n_tiles):
        row = i // n_cols
        col = i % n_cols
        tile_steps.append({
            'args': [{'visible': [True]}, {'title': f'Tile {i+1} (row {row}, col {col})'}],
            'label': str(i + 1),
            'method': 'update'
        })

    fig.update_layout(
        title=dict(
            text=f"Raw Tiles: cyc{cycle:02d} CH{channel}<br><sub>{_stats_text(stats)}</sub>",
            x=0.5
        ),
        height=height,
        yaxis=dict(scaleanchor='x', scaleratio=1, autorange='reversed'),
        sliders=[
            dict(
                active=0,
                currentvalue={'prefix': 'Z-plane: '},
                pad={'t': 50},
                steps=z_steps,
                y=0,
                yanchor='top'
            ),
            dict(
                active=0,
                currentvalue={'prefix': 'Tile: '},
                pad={'t': 100},
                steps=tile_steps,
                y=-0.15,
                yanchor='top'
            )
        ]
    )

    # Store state for callbacks
    state = {'z_idx': 0, 'tile_idx': 0}

    def update_image(z_idx: int, tile_idx: int):
        """Update displayed image."""
        img = load_tile(z_idx, tile_idx)
        img_display = _normalize_for_display(_downsample(img, max_display_size))
        stats = _compute_stats(img)

        with fig.batch_update():
            fig.data[0].z = img_display
            row = tile_idx // n_cols
            col = tile_idx % n_cols
            z = z_planes[z_idx]
            fig.layout.title.text = (
                f"Raw Tiles: cyc{cycle:02d} CH{channel} - Z{z:03d} Tile {tile_idx+1} "
                f"(row {row}, col {col})<br><sub>{_stats_text(stats)}</sub>"
            )

    # Note: Plotly sliders in FigureWidget require JavaScript callbacks
    # For full interactivity, use in Jupyter with ipywidgets integration
    # This provides the visual framework - manual slider interaction updates via JS

    return fig


# =============================================================================
# Z-STACK VIEWER (Stitched / Deconvolved)
# =============================================================================

def view_zstack(
    stack_dir: str | Path,
    title: str = "Z-Stack",
    max_display_size: int = 1000,
    height: int = 700,
    initial_z: int | None = None,
) -> go.Figure:
    """
    Interactive viewer for z-stack (stitched or deconvolved) with z-slider.

    Uses go.Figure with animation frames for slider and play/pause controls.
    For ipywidgets-based interactivity in Jupyter, use create_qc_panel() instead.

    Parameters
    ----------
    stack_dir : Path
        Directory containing z-plane TIFF files (01.tif, 02.tif, etc.)
    title : str
        Title for the figure
    max_display_size : int
        Maximum dimension for downsampled display
    height : int
        Figure height in pixels
    initial_z : int
        Initial z-plane to display (0-indexed). Default: middle plane.

    Returns
    -------
    go.Figure
        Plotly figure with animation slider
    """
    _check_deps()

    stack_dir = Path(stack_dir)

    # Find all z-plane files
    tiff_files = sorted(stack_dir.glob("*.tif"))

    if not tiff_files:
        raise FileNotFoundError(f"No TIFF files found in {stack_dir}")

    n_planes = len(tiff_files)

    if initial_z is None:
        initial_z = n_planes // 2
    initial_z = max(0, min(initial_z, n_planes - 1))

    print(f"Found {n_planes} z-planes in {stack_dir.name}")

    # Preload all images (downsampled) for smooth slider interaction
    print("Loading and downsampling images...")
    images = []
    stats_list = []

    for f in tiff_files:
        img = tifffile.imread(f)
        stats_list.append(_compute_stats(img))
        img_small = _downsample(img, max_display_size)
        img_norm = _normalize_for_display(img_small)
        images.append(img_norm)

    print(f"Loaded {len(images)} images")

    # Create frames for animation
    frames = []
    for i, (img, stats) in enumerate(zip(images, stats_list)):
        frames.append(go.Frame(
            data=[go.Heatmap(z=img, colorscale='gray', showscale=False)],
            name=str(i),
            layout=go.Layout(
                title=dict(
                    text=f"{title} - Z{i+1:02d}/{n_planes}<br><sub>{_stats_text(stats)}</sub>",
                    x=0.5
                )
            )
        ))

    # Initial figure - use go.Figure (not FigureWidget) to support frames
    fig = go.Figure(
        data=[go.Heatmap(
            z=images[initial_z],
            colorscale='gray',
            showscale=False,
            hovertemplate='x: %{x}<br>y: %{y}<br>value: %{z}<extra></extra>'
        )],
        frames=frames
    )

    # Slider steps
    slider_steps = []
    for i in range(n_planes):
        slider_steps.append({
            'args': [
                [str(i)],
                {'frame': {'duration': 0, 'redraw': True}, 'mode': 'immediate'}
            ],
            'label': f'{i+1:02d}',
            'method': 'animate'
        })

    fig.update_layout(
        title=dict(
            text=f"{title} - Z{initial_z+1:02d}/{n_planes}<br><sub>{_stats_text(stats_list[initial_z])}</sub>",
            x=0.5
        ),
        height=height,
        yaxis=dict(scaleanchor='x', scaleratio=1, autorange='reversed'),
        sliders=[{
            'active': initial_z,
            'currentvalue': {'prefix': 'Z-plane: ', 'visible': True},
            'pad': {'t': 50, 'b': 10},
            'steps': slider_steps,
            'transition': {'duration': 0}
        }],
        updatemenus=[{
            'type': 'buttons',
            'showactive': False,
            'y': 0,
            'x': 0.1,
            'xanchor': 'right',
            'yanchor': 'top',
            'buttons': [
                {
                    'label': '> Play',
                    'method': 'animate',
                    'args': [None, {
                        'frame': {'duration': 200, 'redraw': True},
                        'fromcurrent': True,
                        'mode': 'immediate'
                    }]
                },
                {
                    'label': '|| Pause',
                    'method': 'animate',
                    'args': [[None], {
                        'frame': {'duration': 0, 'redraw': False},
                        'mode': 'immediate'
                    }]
                }
            ]
        }]
    )

    return fig


# =============================================================================
# EDF VIEWER
# =============================================================================

def view_edf(
    edf_path: str | Path,
    title: str = "Extended Depth of Focus",
    max_display_size: int = 1000,
    height: int = 700,
) -> go.Figure:
    """
    Viewer for EDF (single 2D image) with zoom/pan and stats.

    Parameters
    ----------
    edf_path : Path
        Path to EDF TIFF file
    title : str
        Title for the figure
    max_display_size : int
        Maximum dimension for downsampled display
    height : int
        Figure height in pixels

    Returns
    -------
    go.Figure
        Plotly figure
    """
    _check_deps()

    edf_path = Path(edf_path)

    if not edf_path.exists():
        raise FileNotFoundError(f"EDF file not found: {edf_path}")

    img = tifffile.imread(edf_path)
    stats = _compute_stats(img)

    print(f"Loaded EDF: {img.shape}, range [{img.min()}, {img.max()}]")

    img_display = _normalize_for_display(_downsample(img, max_display_size))

    fig = go.Figure(
        data=[go.Heatmap(
            z=img_display,
            colorscale='gray',
            showscale=True,
            colorbar=dict(title='Intensity'),
            hovertemplate='x: %{x}<br>y: %{y}<br>value: %{z}<extra></extra>'
        )]
    )

    fig.update_layout(
        title=dict(
            text=f"{title}<br><sub>{_stats_text(stats)}</sub>",
            x=0.5
        ),
        height=height,
        yaxis=dict(scaleanchor='x', scaleratio=1, autorange='reversed'),
    )

    return fig


# =============================================================================
# PIPELINE COMPARISON VIEWER
# =============================================================================

def view_pipeline_comparison(
    project_dir: str | Path,
    cycle: int,
    channel: int,
    z_plane: int | None = None,
    max_display_size: int = 600,
    height: int = 500,
) -> go.Figure:
    """
    Side-by-side comparison of all pipeline stages.

    Parameters
    ----------
    project_dir : Path
        Project directory containing data/processed/
    cycle : int
        Cycle number (1-based)
    channel : int
        Channel number (1-based)
    z_plane : int
        Z-plane to display (1-based). Default: middle plane.
    max_display_size : int
        Maximum dimension for each image
    height : int
        Figure height in pixels

    Returns
    -------
    go.Figure
        Plotly figure with subplots
    """
    _check_deps()

    project_dir = Path(project_dir)
    proc_dir = project_dir / "data" / "processed"

    # Find paths
    stitch_dir = proc_dir / "stitched" / f"cyc{cycle:02d}" / f"CH{channel}"
    decon_dir = proc_dir / "deconvolved" / f"cyc{cycle:02d}" / f"CH{channel}"
    edf_dir = proc_dir / "edf" / f"cyc{cycle:02d}" / f"CH{channel}"

    # Determine z-plane
    if stitch_dir.exists():
        z_files = sorted(stitch_dir.glob("*.tif"))
        n_planes = len(z_files)
        if z_plane is None:
            z_plane = n_planes // 2 + 1
        z_plane = max(1, min(z_plane, n_planes))
    else:
        z_plane = z_plane or 7
        n_planes = 13

    # Load images
    images = {}
    titles = []

    # Stitched
    stitch_path = stitch_dir / f"{z_plane:02d}.tif"
    if stitch_path.exists():
        img = tifffile.imread(stitch_path)
        images['stitched'] = _normalize_for_display(_downsample(img, max_display_size))
        stats = _compute_stats(img)
        titles.append(f"Stitched Z{z_plane:02d}<br><sub>Mean: {stats['mean']:.0f}</sub>")
    else:
        titles.append("Stitched (not found)")

    # Deconvolved
    decon_path = decon_dir / f"deconv_{z_plane:03d}.tif"
    if decon_path.exists():
        img = tifffile.imread(decon_path)
        images['decon'] = _normalize_for_display(_downsample(img, max_display_size))
        stats = _compute_stats(img)
        titles.append(f"Deconvolved Z{z_plane:02d}<br><sub>Mean: {stats['mean']:.0f}</sub>")
    else:
        titles.append("Deconvolved (not found)")

    # EDF
    edf_path = edf_dir / "edf.tif"
    if edf_path.exists():
        img = tifffile.imread(edf_path)
        images['edf'] = _normalize_for_display(_downsample(img, max_display_size))
        stats = _compute_stats(img)
        titles.append(f"EDF<br><sub>Mean: {stats['mean']:.0f}</sub>")
    else:
        titles.append("EDF (not found)")

    # Create subplots
    n_cols = len(titles)
    fig = make_subplots(
        rows=1, cols=n_cols,
        subplot_titles=titles,
        horizontal_spacing=0.02
    )

    # Add images
    col = 1
    for key in ['stitched', 'decon', 'edf']:
        if key in images:
            fig.add_trace(
                go.Heatmap(
                    z=images[key],
                    colorscale='gray',
                    showscale=False
                ),
                row=1, col=col
            )
        col += 1

    fig.update_layout(
        title=dict(
            text=f"Pipeline Comparison: cyc{cycle:02d} CH{channel}",
            x=0.5
        ),
        height=height,
    )

    # Match aspect ratios
    for i in range(1, n_cols + 1):
        fig.update_yaxes(scaleanchor=f'x{i}' if i > 1 else 'x', scaleratio=1,
                        autorange='reversed', row=1, col=i)

    return fig


# =============================================================================
# MULTI-CHANNEL VIEWER
# =============================================================================

def view_multichannel(
    stack_dir: str | Path,
    channels: list[int] = [1, 2, 3, 4],
    z_plane: int | None = None,
    max_display_size: int = 500,
    height: int = 400,
) -> go.Figure:
    """
    View multiple channels side-by-side for a single z-plane.

    Parameters
    ----------
    stack_dir : Path
        Base directory containing CH1/, CH2/, etc. subdirectories
    channels : list
        Channel numbers to display
    z_plane : int
        Z-plane to display (1-based). Default: middle plane.
    max_display_size : int
        Maximum dimension for each image
    height : int
        Figure height in pixels

    Returns
    -------
    go.Figure
        Plotly figure with subplots
    """
    _check_deps()

    stack_dir = Path(stack_dir)

    # Find z-planes
    ch1_dir = stack_dir / f"CH{channels[0]}"
    if ch1_dir.exists():
        z_files = sorted(ch1_dir.glob("*.tif"))
        n_planes = len(z_files)
        if z_plane is None:
            z_plane = n_planes // 2 + 1
    else:
        z_plane = z_plane or 7

    # Load images
    images = {}
    titles = []

    for ch in channels:
        ch_dir = stack_dir / f"CH{ch}"
        z_path = ch_dir / f"{z_plane:02d}.tif"

        if z_path.exists():
            img = tifffile.imread(z_path)
            images[ch] = _normalize_for_display(_downsample(img, max_display_size))
            stats = _compute_stats(img)
            titles.append(f"CH{ch}<br><sub>Mean: {stats['mean']:.0f}</sub>")
        else:
            titles.append(f"CH{ch} (not found)")

    # Create subplots
    n_cols = len(channels)
    fig = make_subplots(
        rows=1, cols=n_cols,
        subplot_titles=titles,
        horizontal_spacing=0.02
    )

    # Add images
    for i, ch in enumerate(channels, 1):
        if ch in images:
            fig.add_trace(
                go.Heatmap(
                    z=images[ch],
                    colorscale='gray',
                    showscale=False
                ),
                row=1, col=i
            )

    fig.update_layout(
        title=dict(text=f"Multi-Channel View - Z{z_plane:02d}", x=0.5),
        height=height,
    )

    for i in range(1, n_cols + 1):
        fig.update_yaxes(scaleanchor=f'x{i}' if i > 1 else 'x', scaleratio=1,
                        autorange='reversed', row=1, col=i)

    return fig


# =============================================================================
# INTERACTIVE QC PANEL (ipywidgets integration)
# =============================================================================

def create_qc_panel(
    project_dir: str | Path,
    cycle: int,
    channel: int,
    stage: Literal['raw', 'stitched', 'deconvolved', 'edf'] = 'stitched',
) -> "QCPanel":
    """
    Create an interactive QC panel with ipywidgets for full interactivity.

    Parameters
    ----------
    project_dir : Path
        Project directory
    cycle : int
        Cycle number
    channel : int
        Channel number
    stage : str
        Pipeline stage to view

    Returns
    -------
    QCPanel
        Interactive panel object with .show() method
    """
    return QCPanel(project_dir, cycle, channel, stage)


class QCPanel:
    """
    Interactive QC panel using ipywidgets for slider callbacks.

    This provides full interactivity in Jupyter notebooks.
    """

    def __init__(
        self,
        project_dir: str | Path,
        cycle: int,
        channel: int,
        stage: str = 'stitched',
        max_display_size: int = 800,
    ):
        self.project_dir = Path(project_dir)
        self.cycle = cycle
        self.channel = channel
        self.stage = stage
        self.max_display_size = max_display_size

        self._setup_paths()
        self._load_data()

    def _setup_paths(self):
        """Setup directory paths."""
        proc_dir = self.project_dir / "data" / "processed"
        raw_dir = self.project_dir / "data" / "raw"

        self.paths = {
            'raw': raw_dir / f"cyc{self.cycle:03d}",
            'stitched': proc_dir / "stitched" / f"cyc{self.cycle:02d}" / f"CH{self.channel}",
            'deconvolved': proc_dir / "deconvolved" / f"cyc{self.cycle:02d}" / f"CH{self.channel}",
            'edf': proc_dir / "edf" / f"cyc{self.cycle:02d}" / f"CH{self.channel}",
        }

    def _load_data(self):
        """Load image data for current stage."""
        self.images = []
        self.stats = []

        path = self.paths[self.stage]

        if self.stage == 'edf':
            edf_file = path / "edf.tif"
            if edf_file.exists():
                img = tifffile.imread(edf_file)
                self.images.append(_downsample(img, self.max_display_size))
                self.stats.append(_compute_stats(img))
        else:
            if self.stage == 'deconvolved':
                pattern = "deconv_*.tif"
            else:
                pattern = "*.tif"

            files = sorted(path.glob(pattern))
            for f in files:
                img = tifffile.imread(f)
                self.images.append(_downsample(img, self.max_display_size))
                self.stats.append(_compute_stats(img))

        self.n_images = len(self.images)
        print(f"Loaded {self.n_images} images from {self.stage}")

    def show(self):
        """Display the interactive panel."""
        try:
            import ipywidgets as widgets
            from IPython.display import display
        except ImportError:
            print("ipywidgets required for interactive panel. Install with: pip install ipywidgets")
            # Fall back to static view
            if self.stage == 'edf':
                return view_edf(self.paths['edf'] / "edf.tif")
            else:
                return view_zstack(self.paths[self.stage], title=f"{self.stage.title()} cyc{self.cycle:02d} CH{self.channel}")

        # Create figure widget
        initial_idx = min(self.n_images // 2, self.n_images - 1) if self.n_images > 0 else 0

        if self.n_images == 0:
            print(f"No images found for {self.stage}")
            return None

        fig = go.FigureWidget(
            data=[go.Heatmap(
                z=_normalize_for_display(self.images[initial_idx]),
                colorscale='gray',
                showscale=False
            )]
        )

        fig.update_layout(
            title=f"{self.stage.title()} cyc{self.cycle:02d} CH{self.channel}",
            height=600,
            yaxis=dict(scaleanchor='x', scaleratio=1, autorange='reversed')
        )

        # Create slider
        if self.n_images > 1:
            slider = widgets.IntSlider(
                value=initial_idx,
                min=0,
                max=self.n_images - 1,
                step=1,
                description='Z-plane:',
                continuous_update=True
            )

            # Stats display
            stats_html = widgets.HTML(value=_stats_text(self.stats[initial_idx]).replace('<br>', ' | '))

            def on_slider_change(change):
                idx = change['new']
                img_display = _normalize_for_display(self.images[idx])
                with fig.batch_update():
                    fig.data[0].z = img_display
                    fig.layout.title = f"{self.stage.title()} cyc{self.cycle:02d} CH{self.channel} - Z{idx+1:02d}"
                stats_html.value = _stats_text(self.stats[idx]).replace('<br>', ' | ')

            slider.observe(on_slider_change, names='value')

            # Play button
            play = widgets.Play(
                value=initial_idx,
                min=0,
                max=self.n_images - 1,
                step=1,
                interval=200,
                description="Play"
            )
            widgets.jslink((play, 'value'), (slider, 'value'))

            display(widgets.VBox([
                widgets.HBox([play, slider]),
                stats_html,
                fig
            ]))
        else:
            display(fig)

        return fig
