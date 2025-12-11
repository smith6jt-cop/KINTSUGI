"""
Interactive Parameter Tuning for Claude-Guided Workflow

Provides Jupyter widgets for real-time parameter adjustment
with visual feedback using Kview2 visualization tools.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, Optional, Callable
import numpy as np

try:
    import ipywidgets as widgets
    from IPython.display import display, clear_output
    WIDGETS_AVAILABLE = True
except ImportError:
    WIDGETS_AVAILABLE = False

try:
    import stackview
    STACKVIEW_AVAILABLE = True
except ImportError:
    STACKVIEW_AVAILABLE = False


@dataclass
class TunerResult:
    """Result from parameter tuning session."""

    operation: str
    channel: str
    params: Dict[str, Any] = field(default_factory=dict)
    processed_data: Optional[np.ndarray] = None
    approved: bool = False

    def get_params(self) -> Dict[str, Any]:
        """Get the tuned parameters."""
        return self.params.copy()

    def get_data(self) -> Optional[np.ndarray]:
        """Get the processed data array."""
        return self.processed_data


def _check_widgets():
    """Check if Jupyter widgets are available."""
    if not WIDGETS_AVAILABLE:
        raise ImportError(
            "ipywidgets not installed. Install with: pip install ipywidgets"
        )


def _load_image(channel_name: str) -> np.ndarray:
    """Load image from MCP session state or file."""
    try:
        from kintsugi.mcp.tools.signal_isolation import _loaded_images

        if channel_name in _loaded_images:
            data = _loaded_images[channel_name]["data"]
            if hasattr(data, "compute"):
                return data.compute()
            return data
    except ImportError:
        pass

    raise ValueError(f"Channel not loaded: {channel_name}")


def blank_subtraction_tuner(
    signal_channel: str,
    blank_channel: str,
    initial_params: Optional[Dict[str, Any]] = None,
    zoom_factor: float = 0.25,
) -> TunerResult:
    """
    Interactive blank subtraction parameter tuner.

    Creates sliders for blank_clip_factor, blank_scale_factor, etc.
    with real-time visual feedback.
    """
    _check_widgets()

    initial_params = initial_params or {}
    result = TunerResult(
        operation="blank_subtraction",
        channel=signal_channel,
        params=initial_params.copy(),
    )

    # Load images
    try:
        signal_data = _load_image(signal_channel)
        blank_data = _load_image(blank_channel)
    except ValueError as e:
        print(f"Error: {e}")
        print("Make sure channels are loaded via MCP tools first.")
        return result

    # Downsample for interactive display
    step = max(1, int(1 / zoom_factor))
    signal_display = signal_data[::step, ::step].astype(np.float32)
    blank_display = blank_data[::step, ::step].astype(np.float32)

    # Create widgets
    clip_slider = widgets.IntSlider(
        value=initial_params.get("blank_clip_factor", 0),
        min=0,
        max=5000,
        step=100,
        description="Clip Factor:",
        continuous_update=False,
    )

    scale_slider = widgets.FloatSlider(
        value=initial_params.get("blank_scale_factor", 0.8),
        min=0.0,
        max=2.0,
        step=0.05,
        description="Scale Factor:",
        continuous_update=False,
    )

    smooth_low_check = widgets.Checkbox(
        value=initial_params.get("smooth_low", False),
        description="Smooth Low Intensities",
    )

    smooth_high_check = widgets.Checkbox(
        value=initial_params.get("smooth_high", False),
        description="Smooth High Intensities",
    )

    erosion_slider = widgets.IntSlider(
        value=initial_params.get("erosion", 0),
        min=0,
        max=5,
        step=1,
        description="Erosion:",
    )

    approve_button = widgets.Button(
        description="Approve Parameters",
        button_style="success",
    )

    output = widgets.Output()

    def update_preview(_=None):
        """Update the preview with current parameters."""
        with output:
            clear_output(wait=True)

            # Apply processing
            blank_clipped = np.clip(blank_display, clip_slider.value, blank_display.max())
            blank_clipped[blank_clipped <= clip_slider.value] = 0

            processed = signal_display - np.minimum(
                signal_display, blank_clipped * scale_slider.value
            )
            processed = np.maximum(processed, 0)

            # Update result params
            result.params = {
                "blank_clip_factor": clip_slider.value,
                "blank_scale_factor": scale_slider.value,
                "smooth_low": smooth_low_check.value,
                "smooth_high": smooth_high_check.value,
                "erosion": erosion_slider.value,
            }

            # Display comparison
            if STACKVIEW_AVAILABLE:
                print("Original vs Processed (drag slider to compare)")
                stackview.curtain(
                    signal_display,
                    processed,
                    colormap=["gray", "turbo"],
                )
            else:
                # Fallback to matplotlib
                import matplotlib.pyplot as plt

                fig, axes = plt.subplots(1, 3, figsize=(15, 5))
                axes[0].imshow(signal_display, cmap="gray")
                axes[0].set_title("Signal")
                axes[1].imshow(blank_display, cmap="gray")
                axes[1].set_title("Blank")
                axes[2].imshow(processed, cmap="turbo")
                axes[2].set_title("Processed")
                plt.tight_layout()
                plt.show()

            # Stats
            print(f"\nSignal range: {signal_display.min():.0f} - {signal_display.max():.0f}")
            print(f"Processed range: {processed.min():.0f} - {processed.max():.0f}")
            print(f"Reduction: {100 * (1 - processed.mean() / signal_display.mean()):.1f}%")

    def on_approve(_):
        """Handle approve button click."""
        result.approved = True
        print("✓ Parameters approved!")
        print(f"Parameters: {result.params}")

    # Connect callbacks
    clip_slider.observe(update_preview, names="value")
    scale_slider.observe(update_preview, names="value")
    smooth_low_check.observe(update_preview, names="value")
    smooth_high_check.observe(update_preview, names="value")
    erosion_slider.observe(update_preview, names="value")
    approve_button.on_click(on_approve)

    # Layout
    controls = widgets.VBox([
        widgets.HTML("<h3>Blank Subtraction Parameters</h3>"),
        clip_slider,
        scale_slider,
        smooth_low_check,
        smooth_high_check,
        erosion_slider,
        approve_button,
    ])

    display(widgets.HBox([controls, output]))
    update_preview()

    return result


def denoise_tuner(
    channel: str,
    initial_params: Optional[Dict[str, Any]] = None,
    zoom_factor: float = 0.25,
) -> TunerResult:
    """Interactive denoising parameter tuner."""
    _check_widgets()

    initial_params = initial_params or {}
    result = TunerResult(
        operation="denoise",
        channel=channel,
        params=initial_params.copy(),
    )

    try:
        data = _load_image(channel)
    except ValueError as e:
        print(f"Error: {e}")
        return result

    step = max(1, int(1 / zoom_factor))
    display_data = data[::step, ::step].astype(np.float32)

    # Widgets
    method_dropdown = widgets.Dropdown(
        options=["median", "uniform", "percentile"],
        value=initial_params.get("method", "median"),
        description="Method:",
    )

    size_slider = widgets.IntSlider(
        value=initial_params.get("filter_size", 3),
        min=1,
        max=9,
        step=2,
        description="Filter Size:",
    )

    percentile_slider = widgets.IntSlider(
        value=initial_params.get("percentile", 10),
        min=1,
        max=50,
        description="Percentile:",
    )

    approve_button = widgets.Button(description="Approve", button_style="success")
    output = widgets.Output()

    def update_preview(_=None):
        with output:
            clear_output(wait=True)
            from scipy.ndimage import median_filter, uniform_filter, percentile_filter

            method = method_dropdown.value
            size = size_slider.value

            if method == "median":
                processed = median_filter(display_data, size=size)
            elif method == "uniform":
                processed = uniform_filter(display_data, size=size)
            else:
                processed = percentile_filter(display_data, percentile_slider.value, size=size)

            result.params = {
                "method": method,
                "filter_size": size,
                "percentile": percentile_slider.value,
            }

            if STACKVIEW_AVAILABLE:
                stackview.curtain(display_data, processed, colormap=["gray", "turbo"])
            else:
                import matplotlib.pyplot as plt
                fig, axes = plt.subplots(1, 2, figsize=(10, 5))
                axes[0].imshow(display_data, cmap="gray")
                axes[0].set_title("Original")
                axes[1].imshow(processed, cmap="turbo")
                axes[1].set_title("Denoised")
                plt.show()

    def on_approve(_):
        result.approved = True
        print("✓ Parameters approved:", result.params)

    method_dropdown.observe(update_preview, names="value")
    size_slider.observe(update_preview, names="value")
    percentile_slider.observe(update_preview, names="value")
    approve_button.on_click(on_approve)

    controls = widgets.VBox([
        widgets.HTML("<h3>Denoising Parameters</h3>"),
        method_dropdown,
        size_slider,
        percentile_slider,
        approve_button,
    ])

    display(widgets.HBox([controls, output]))
    update_preview()

    return result


def clahe_tuner(
    channel: str,
    initial_params: Optional[Dict[str, Any]] = None,
    zoom_factor: float = 0.25,
) -> TunerResult:
    """Interactive CLAHE parameter tuner."""
    _check_widgets()

    initial_params = initial_params or {}
    result = TunerResult(operation="clahe", channel=channel, params=initial_params.copy())

    try:
        data = _load_image(channel)
    except ValueError as e:
        print(f"Error: {e}")
        return result

    step = max(1, int(1 / zoom_factor))
    display_data = data[::step, ::step]

    clip_slider = widgets.FloatSlider(
        value=initial_params.get("clip_limit", 0.01),
        min=0.001,
        max=0.1,
        step=0.005,
        description="Clip Limit:",
        readout_format=".3f",
    )

    tile_slider = widgets.IntSlider(
        value=initial_params.get("tile_grid_size", 70),
        min=10,
        max=200,
        step=10,
        description="Tile Grid:",
    )

    approve_button = widgets.Button(description="Approve", button_style="success")
    output = widgets.Output()

    def update_preview(_=None):
        with output:
            clear_output(wait=True)
            from skimage.exposure import equalize_adapthist, rescale_intensity

            data_uint16 = display_data.astype(np.uint16)
            kernel_size = max(1, data_uint16.shape[0] // tile_slider.value)

            processed = equalize_adapthist(
                data_uint16,
                clip_limit=clip_slider.value,
                kernel_size=kernel_size,
            )
            processed = rescale_intensity(
                processed,
                out_range=(float(data_uint16.min()), float(data_uint16.max())),
            )

            result.params = {
                "clip_limit": clip_slider.value,
                "tile_grid_size": tile_slider.value,
            }

            if STACKVIEW_AVAILABLE:
                stackview.curtain(display_data, processed, colormap=["gray", "turbo"])
            else:
                import matplotlib.pyplot as plt
                fig, axes = plt.subplots(1, 2, figsize=(10, 5))
                axes[0].imshow(display_data, cmap="gray")
                axes[0].set_title("Original")
                axes[1].imshow(processed, cmap="turbo")
                axes[1].set_title("CLAHE")
                plt.show()

    def on_approve(_):
        result.approved = True
        print("✓ Parameters approved:", result.params)

    clip_slider.observe(update_preview, names="value")
    tile_slider.observe(update_preview, names="value")
    approve_button.on_click(on_approve)

    controls = widgets.VBox([
        widgets.HTML("<h3>CLAHE Parameters</h3>"),
        clip_slider,
        tile_slider,
        approve_button,
    ])

    display(widgets.HBox([controls, output]))
    update_preview()

    return result


def clean_tuner(
    channel: str,
    initial_params: Optional[Dict[str, Any]] = None,
    zoom_factor: float = 0.25,
) -> TunerResult:
    """Interactive background cleaning parameter tuner."""
    _check_widgets()

    initial_params = initial_params or {}
    result = TunerResult(operation="clean", channel=channel, params=initial_params.copy())

    try:
        data = _load_image(channel)
    except ValueError as e:
        print(f"Error: {e}")
        return result

    step = max(1, int(1 / zoom_factor))
    display_data = data[::step, ::step].astype(np.float32)

    threshold_slider = widgets.IntSlider(
        value=initial_params.get("background_threshold", 100),
        min=0,
        max=2000,
        step=50,
        description="BG Threshold:",
    )

    remove_small_check = widgets.Checkbox(
        value=initial_params.get("remove_small_objects", False),
        description="Remove Small Objects",
    )

    min_size_slider = widgets.IntSlider(
        value=initial_params.get("min_object_size", 30),
        min=5,
        max=200,
        step=5,
        description="Min Size:",
    )

    approve_button = widgets.Button(description="Approve", button_style="success")
    output = widgets.Output()

    def update_preview(_=None):
        with output:
            clear_output(wait=True)

            processed = display_data.copy()
            processed[processed <= threshold_slider.value] = 0

            if remove_small_check.value:
                from skimage.morphology import remove_small_objects
                binary = processed > 0
                binary = remove_small_objects(binary, min_size=min_size_slider.value)
                processed = np.where(binary, processed, 0)

            result.params = {
                "background_threshold": threshold_slider.value,
                "remove_small_objects": remove_small_check.value,
                "min_object_size": min_size_slider.value,
            }

            if STACKVIEW_AVAILABLE:
                stackview.curtain(display_data, processed, colormap=["gray", "turbo"])
            else:
                import matplotlib.pyplot as plt
                fig, axes = plt.subplots(1, 2, figsize=(10, 5))
                axes[0].imshow(display_data, cmap="gray")
                axes[1].imshow(processed, cmap="turbo")
                plt.show()

    def on_approve(_):
        result.approved = True
        print("✓ Parameters approved:", result.params)

    threshold_slider.observe(update_preview, names="value")
    remove_small_check.observe(update_preview, names="value")
    min_size_slider.observe(update_preview, names="value")
    approve_button.on_click(on_approve)

    controls = widgets.VBox([
        widgets.HTML("<h3>Background Cleaning Parameters</h3>"),
        threshold_slider,
        remove_small_check,
        min_size_slider,
        approve_button,
    ])

    display(widgets.HBox([controls, output]))
    update_preview()

    return result


def gaussian_subtract_tuner(
    channel: str,
    initial_params: Optional[Dict[str, Any]] = None,
    zoom_factor: float = 0.25,
) -> TunerResult:
    """Interactive Gaussian subtraction parameter tuner."""
    _check_widgets()

    initial_params = initial_params or {}
    result = TunerResult(
        operation="gaussian_subtract", channel=channel, params=initial_params.copy()
    )

    try:
        data = _load_image(channel)
    except ValueError as e:
        print(f"Error: {e}")
        return result

    step = max(1, int(1 / zoom_factor))
    display_data = data[::step, ::step].astype(np.float32)

    sigma_slider = widgets.IntSlider(
        value=initial_params.get("sigma", 10),
        min=1,
        max=50,
        step=1,
        description="Sigma:",
    )

    scale_slider = widgets.FloatSlider(
        value=initial_params.get("scale_factor", 0.1),
        min=0.0,
        max=1.0,
        step=0.05,
        description="Scale:",
    )

    approve_button = widgets.Button(description="Approve", button_style="success")
    output = widgets.Output()

    def update_preview(_=None):
        with output:
            clear_output(wait=True)
            from scipy.ndimage import gaussian_filter

            blurred = gaussian_filter(display_data, sigma=sigma_slider.value)
            processed = display_data - np.minimum(display_data, blurred * scale_slider.value)
            processed = np.maximum(processed, 0)

            result.params = {
                "sigma": sigma_slider.value,
                "scale_factor": scale_slider.value,
            }

            if STACKVIEW_AVAILABLE:
                stackview.curtain(display_data, processed, colormap=["gray", "turbo"])
            else:
                import matplotlib.pyplot as plt
                fig, axes = plt.subplots(1, 3, figsize=(15, 5))
                axes[0].imshow(display_data, cmap="gray")
                axes[0].set_title("Original")
                axes[1].imshow(blurred, cmap="gray")
                axes[1].set_title("Gaussian Blur")
                axes[2].imshow(processed, cmap="turbo")
                axes[2].set_title("Subtracted")
                plt.show()

    def on_approve(_):
        result.approved = True
        print("✓ Parameters approved:", result.params)

    sigma_slider.observe(update_preview, names="value")
    scale_slider.observe(update_preview, names="value")
    approve_button.on_click(on_approve)

    controls = widgets.VBox([
        widgets.HTML("<h3>Gaussian Subtraction Parameters</h3>"),
        sigma_slider,
        scale_slider,
        approve_button,
    ])

    display(widgets.HBox([controls, output]))
    update_preview()

    return result
