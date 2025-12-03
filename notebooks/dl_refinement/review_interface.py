"""
Interactive Channel Review Interface

Napari-based GUI for reviewing channels flagged by quality assessment.
"""

import numpy as np
import zarr
from pathlib import Path
from typing import Dict, List, Optional, Callable, Tuple
import json
import logging

try:
    import napari
    from magicgui import magicgui
    NAPARI_AVAILABLE = True
except ImportError:
    NAPARI_AVAILABLE = False
    logging.warning("Napari not available. Install with: pip install napari[all]")

logger = logging.getLogger(__name__)


class ChannelReviewInterface:
    """Interactive interface for reviewing flagged channels.

    Provides a Napari-based GUI with tools for inspecting problem regions,
    making decisions, and optionally applying adjustments.

    Example:
        >>> reviewer = ChannelReviewInterface(
        ...     results_summary_path='results/processing_summary.json',
        ...     image_data_path='data/processed.zarr'
        ... )
        >>> reviewer.start_review()  # Opens interactive GUI
        >>> decisions = reviewer.get_review_decisions()
    """

    def __init__(
        self,
        results_summary_path: str,
        image_data_path: str,
        adjustment_callback: Optional[Callable] = None,
        autosave: bool = True
    ):
        """Initialize review interface.

        Args:
            results_summary_path: Path to JSON summary from BatchChannelProcessor
            image_data_path: Path to Zarr array or directory with images
            adjustment_callback: Optional function to apply adjustments
            autosave: Automatically save decisions after each channel
        """
        if not NAPARI_AVAILABLE:
            raise ImportError("Napari is required for the review interface. "
                            "Install with: pip install napari[all]")

        self.results = self._load_results(results_summary_path)
        self.results_path = Path(results_summary_path)
        self.image_data_path = Path(image_data_path)
        self.adjustment_callback = adjustment_callback
        self.autosave = autosave

        # Open image data
        if self.image_data_path.suffix == '.zarr' or self.image_data_path.is_dir():
            try:
                self.image_data = zarr.open(str(self.image_data_path), mode='r')
                self.data_type = 'zarr'
            except Exception:
                self.data_type = 'directory'
                self.image_data = None
        else:
            self.data_type = 'file'
            self.image_data = None

        # State
        self.viewer = None
        self.current_channel_idx = 0
        self.channels_to_review = self.results['channels_requiring_review']
        self.review_decisions = {}
        self.adjustment_params = {}

        # Load existing decisions if available
        self._load_decisions()

        logger.info(f"Initialized review interface for {len(self.channels_to_review)} channels")

    def _load_results(self, path: str) -> Dict:
        """Load processing results from JSON."""
        with open(path, 'r') as f:
            return json.load(f)

    def _load_decisions(self):
        """Load existing review decisions if available."""
        decisions_path = self.results_path.parent / 'review_decisions.json'
        if decisions_path.exists():
            with open(decisions_path, 'r') as f:
                saved = json.load(f)
                self.review_decisions = saved.get('decisions', {})
                self.adjustment_params = saved.get('adjustments', {})
                logger.info(f"Loaded {len(self.review_decisions)} existing decisions")

    def _save_decisions(self):
        """Save review decisions to JSON."""
        decisions_path = self.results_path.parent / 'review_decisions.json'
        output = {
            'timestamp': self.results['timestamp'],
            'decisions': self.review_decisions,
            'adjustments': self.adjustment_params,
            'summary': {
                'total_reviewed': len(self.review_decisions),
                'approved': sum(1 for d in self.review_decisions.values() if d == 'approved'),
                'manual': sum(1 for d in self.review_decisions.values() if d == 'manual'),
                'adjusted': sum(1 for d in self.review_decisions.values() if d == 'adjusted'),
            }
        }

        with open(decisions_path, 'w') as f:
            json.dump(output, f, indent=2)

        logger.info(f"Saved decisions to {decisions_path}")

    def start_review(self):
        """Launch interactive review interface."""
        if not self.channels_to_review:
            logger.info("No channels require review!")
            print("✓ All channels passed quality assessment!")
            return

        # Create viewer
        self.viewer = napari.Viewer(title="KINTSUGI Channel Review")

        # Load first channel
        self._load_channel(0)

        # Add widgets
        self._setup_widgets()

        # Start napari event loop
        napari.run()

        # Save final decisions
        self._save_decisions()

        logger.info("Review session complete")

    def _setup_widgets(self):
        """Setup GUI widgets."""
        # Navigation slider
        @magicgui(
            auto_call=True,
            channel_idx={
                "widget_type": "Slider",
                "min": 0,
                "max": max(0, len(self.channels_to_review) - 1),
                "label": "Channel Index"
            }
        )
        def navigate_channels(channel_idx: int = 0):
            """Navigate to specific channel."""
            self._load_channel(channel_idx)

        # Info display
        @magicgui(
            call_button=False,
            channel_info={"widget_type": "TextEdit", "label": "Channel Info"}
        )
        def info_widget(channel_info: str = ""):
            """Display channel information."""
            pass

        # Update info
        self._update_info_widget(info_widget)

        # Decision buttons
        @magicgui(call_button="✓ Approve Channel")
        def approve_channel():
            """Approve current channel."""
            self._save_decision('approved')
            self._next_channel()

        @magicgui(call_button="⚠ Flag for Manual Processing")
        def flag_channel():
            """Flag channel for manual processing."""
            self._save_decision('manual')
            self._next_channel()

        @magicgui(call_button="← Previous Channel")
        def previous_channel():
            """Go to previous channel."""
            self._previous_channel()

        @magicgui(call_button="Next Channel →")
        def next_channel():
            """Go to next channel."""
            self._next_channel()

        @magicgui(call_button="Skip (Decide Later)")
        def skip_channel():
            """Skip this channel for now."""
            self._next_channel()

        # Adjustment widget
        @magicgui(
            blank_scale={"widget_type": "FloatSpinBox", "min": 0.0, "max": 2.0, "step": 0.1, "value": 1.0},
            blank_clip={"widget_type": "FloatSpinBox", "min": 0, "max": 65535, "step": 100, "value": 1000},
            call_button="Apply Adjustment & Approve"
        )
        def apply_adjustment(blank_scale: float = 1.0, blank_clip: float = 1000):
            """Apply adjustment parameters."""
            params = {
                'blank_scale_factor': blank_scale,
                'blank_clip_factor': blank_clip
            }
            self._save_decision('adjusted', params)
            self._next_channel()

        # Add all widgets to viewer
        self.viewer.window.add_dock_widget(navigate_channels, area='left', name='Navigation')
        self.viewer.window.add_dock_widget(info_widget, area='left', name='Info')
        self.viewer.window.add_dock_widget(approve_channel, area='right', name='Approve')
        self.viewer.window.add_dock_widget(flag_channel, area='right', name='Flag')
        self.viewer.window.add_dock_widget(previous_channel, area='bottom', name='Previous')
        self.viewer.window.add_dock_widget(next_channel, area='bottom', name='Next')
        self.viewer.window.add_dock_widget(skip_channel, area='bottom', name='Skip')
        self.viewer.window.add_dock_widget(apply_adjustment, area='right', name='Adjust')

        # Store widget references
        self.info_widget = info_widget
        self.navigate_widget = navigate_channels

    def _update_info_widget(self, widget):
        """Update info widget with current channel information."""
        if self.current_channel_idx >= len(self.channels_to_review):
            widget.channel_info.value = "Review complete!"
            return

        channel_name = self.channels_to_review[self.current_channel_idx]
        details = self.results['channel_details'].get(channel_name, {})

        info_text = f"""Channel: {channel_name}
Index: {self.current_channel_idx + 1} / {len(self.channels_to_review)}

Confidence Score: {details.get('confidence_score', 0):.4f}
Threshold: {0.85}

Quality Metrics:
  SNR: {details.get('quality_metrics', {}).get('snr', 0):.2f}
  Mean Intensity: {details.get('quality_metrics', {}).get('mean_intensity', 0):.2f}
  Low Confidence Ratio: {details.get('quality_metrics', {}).get('low_confidence_ratio', 0):.4f}

Problem Regions: {details.get('n_problem_regions', 0)}

Previous Decision: {self.review_decisions.get(channel_name, 'None')}
"""
        widget.channel_info.value = info_text

    def _load_channel(self, idx: int):
        """Load channel for review.

        Args:
            idx: Index in channels_to_review list
        """
        if idx < 0 or idx >= len(self.channels_to_review):
            logger.warning(f"Invalid channel index: {idx}")
            return

        self.current_channel_idx = idx
        channel_name = self.channels_to_review[idx]

        logger.info(f"Loading channel {idx + 1}/{len(self.channels_to_review)}: {channel_name}")

        # Parse channel name to get indices
        channel_data = self._load_channel_data(channel_name)

        if channel_data is None:
            logger.error(f"Failed to load data for {channel_name}")
            return

        # Clear existing layers
        self.viewer.layers.clear()

        # Add channel data to viewer
        self.viewer.add_image(
            channel_data,
            name=channel_name,
            colormap='gray',
            contrast_limits=self._get_contrast_limits(channel_data)
        )

        # Highlight problem regions
        self._highlight_problem_regions(channel_name)

        # Update title and info
        self.viewer.title = f"Review: {channel_name} ({idx + 1}/{len(self.channels_to_review)})"
        if hasattr(self, 'info_widget'):
            self._update_info_widget(self.info_widget)

    def _load_channel_data(self, channel_name: str) -> Optional[np.ndarray]:
        """Load image data for a channel.

        Args:
            channel_name: Channel identifier

        Returns:
            2D numpy array or None if failed
        """
        if self.data_type == 'zarr':
            # Parse channel name: "cycle_X_channel_Y" or similar
            parts = channel_name.split('_')
            try:
                # Find cycle and channel indices
                cycle_idx = None
                chan_idx = None

                for i, part in enumerate(parts):
                    if part == 'cycle' and i + 1 < len(parts):
                        cycle_idx = int(parts[i + 1])
                    elif part in ['channel', 'chan'] and i + 1 < len(parts):
                        # Handle case where channel name might be the last part
                        try:
                            chan_idx = int(parts[i + 1])
                        except ValueError:
                            # Channel name is not an integer, try to find it in metadata
                            pass

                # Load from zarr
                if len(self.image_data.shape) == 4 and cycle_idx is not None and chan_idx is not None:
                    return np.array(self.image_data[cycle_idx, chan_idx, :, :])
                elif len(self.image_data.shape) == 3 and chan_idx is not None:
                    return np.array(self.image_data[chan_idx, :, :])
                elif len(self.image_data.shape) == 2:
                    return np.array(self.image_data[:, :])
                else:
                    logger.error(f"Unable to parse channel indices from: {channel_name}")
                    return None

            except Exception as e:
                logger.error(f"Error loading from zarr: {e}")
                return None

        elif self.data_type == 'directory':
            # Try to find matching file
            from imageio import imread
            pattern = f"*{channel_name}*.tif*"
            matching_files = list(self.image_data_path.glob(pattern))

            if matching_files:
                return imread(matching_files[0])
            else:
                logger.error(f"No matching file found for {channel_name}")
                return None

        return None

    def _get_contrast_limits(self, data: np.ndarray) -> Tuple[float, float]:
        """Calculate appropriate contrast limits for display.

        Args:
            data: Image data

        Returns:
            Tuple of (min, max) contrast values
        """
        p2, p98 = np.percentile(data, [2, 98])
        return (p2, p98)

    def _highlight_problem_regions(self, channel_name: str):
        """Highlight problematic regions identified by assessment.

        Args:
            channel_name: Channel identifier
        """
        details = self.results['channel_details'].get(channel_name, {})
        problem_regions = details.get('problem_regions', [])

        if not problem_regions:
            return

        # Create rectangles for each problem region
        shapes_data = []
        for x1, y1, x2, y2 in problem_regions:
            # Napari uses (y, x) coordinates
            rect = np.array([[y1, x1], [y1, x2], [y2, x2], [y2, x1]])
            shapes_data.append(rect)

        # Add shapes layer
        self.viewer.add_shapes(
            shapes_data,
            shape_type='rectangle',
            edge_color='red',
            edge_width=3,
            face_color='transparent',
            name='Problem Regions',
            opacity=0.7
        )

    def _save_decision(self, decision: str, params: Optional[Dict] = None):
        """Save decision for current channel.

        Args:
            decision: 'approved', 'manual', or 'adjusted'
            params: Optional adjustment parameters
        """
        channel_name = self.channels_to_review[self.current_channel_idx]
        self.review_decisions[channel_name] = decision

        if params:
            self.adjustment_params[channel_name] = params

        logger.info(f"Decision for {channel_name}: {decision}")

        if self.autosave:
            self._save_decisions()

    def _next_channel(self):
        """Move to next channel."""
        if self.current_channel_idx < len(self.channels_to_review) - 1:
            self._load_channel(self.current_channel_idx + 1)
            if hasattr(self, 'navigate_widget'):
                self.navigate_widget.channel_idx.value = self.current_channel_idx
        else:
            logger.info("Reached last channel")
            if hasattr(self, 'info_widget'):
                self.info_widget.channel_info.value = "✓ Review Complete!\nAll channels reviewed."

    def _previous_channel(self):
        """Move to previous channel."""
        if self.current_channel_idx > 0:
            self._load_channel(self.current_channel_idx - 1)
            if hasattr(self, 'navigate_widget'):
                self.navigate_widget.channel_idx.value = self.current_channel_idx

    def get_review_decisions(self) -> Dict:
        """Get all review decisions.

        Returns:
            Dictionary mapping channel names to decisions
        """
        return {
            'decisions': self.review_decisions,
            'adjustments': self.adjustment_params,
            'summary': {
                'total_reviewed': len(self.review_decisions),
                'approved': sum(1 for d in self.review_decisions.values() if d == 'approved'),
                'manual': sum(1 for d in self.review_decisions.values() if d == 'manual'),
                'adjusted': sum(1 for d in self.review_decisions.values() if d == 'adjusted'),
            }
        }

    def export_approved_channels(self, output_path: str):
        """Export list of approved channels for downstream processing.

        Args:
            output_path: Path to save approved channels list
        """
        approved = [
            ch for ch, decision in self.review_decisions.items()
            if decision in ['approved', 'adjusted']
        ]

        with open(output_path, 'w') as f:
            json.dump({
                'approved_channels': approved,
                'adjustments': {
                    ch: self.adjustment_params[ch]
                    for ch in approved
                    if ch in self.adjustment_params
                }
            }, f, indent=2)

        logger.info(f"Exported {len(approved)} approved channels to {output_path}")


class BatchReviewInterface:
    """Review multiple datasets in sequence.

    Useful for reviewing results from multiple experiments or batches.
    """

    def __init__(self, summary_paths: List[str], image_data_paths: List[str]):
        """Initialize batch review interface.

        Args:
            summary_paths: List of processing summary JSON paths
            image_data_paths: List of corresponding image data paths
        """
        if len(summary_paths) != len(image_data_paths):
            raise ValueError("Number of summary paths must match image data paths")

        self.datasets = list(zip(summary_paths, image_data_paths))
        self.current_dataset_idx = 0

        logger.info(f"Initialized batch review for {len(self.datasets)} datasets")

    def start_review(self):
        """Start reviewing datasets in sequence."""
        for idx, (summary_path, data_path) in enumerate(self.datasets):
            logger.info(f"=== Reviewing dataset {idx + 1}/{len(self.datasets)} ===")
            print(f"\n{'='*60}")
            print(f"Dataset {idx + 1}/{len(self.datasets)}: {Path(summary_path).parent.name}")
            print(f"{'='*60}\n")

            reviewer = ChannelReviewInterface(
                results_summary_path=summary_path,
                image_data_path=data_path
            )

            reviewer.start_review()

        logger.info("All datasets reviewed!")
        print("\n✓ All datasets have been reviewed!")
