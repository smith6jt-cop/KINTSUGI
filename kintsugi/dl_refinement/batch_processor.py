"""
Batch Channel Processing

Parallel processing of multiple channels with automatic quality assessment
and result aggregation.
"""

import zarr
import numpy as np
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, ProcessPoolExecutor, as_completed
from typing import Dict, List, Optional, Union, Callable
from datetime import datetime
import json
import logging
from tqdm import tqdm
import warnings

from .channel_assessor import ChannelAssessor, ChannelQualityResult

# Setup logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


class BatchChannelProcessor:
    """Process multiple channels in parallel with quality assessment.

    This class orchestrates batch processing of imaging datasets, applying
    quality assessment to each channel and aggregating results for review.

    Example:
        >>> from kintsugi.dl_refinement import ChannelAssessor, BatchChannelProcessor
        >>>
        >>> assessor = ChannelAssessor(model_path='models/quality_net.pth')
        >>> processor = BatchChannelProcessor(
        ...     assessor=assessor,
        ...     output_dir='results/quality_assessment'
        ... )
        >>> results = processor.process_zarr_dataset('data/processed.zarr')
        >>> print(f"Channels needing review: {len(results['channels_requiring_review'])}")
    """

    def __init__(
        self,
        assessor: ChannelAssessor,
        output_dir: Union[str, Path],
        n_workers: int = 4,
        use_processes: bool = False,
        save_confidence_maps: bool = False,
        progress_bar: bool = True
    ):
        """Initialize batch processor.

        Args:
            assessor: ChannelAssessor instance to use for quality assessment
            output_dir: Directory to save results
            n_workers: Number of parallel workers
            use_processes: Use ProcessPoolExecutor instead of ThreadPoolExecutor
            save_confidence_maps: Whether to save per-channel confidence maps
            progress_bar: Show progress bar during processing
        """
        self.assessor = assessor
        self.output_dir = Path(output_dir)
        self.n_workers = n_workers
        self.use_processes = use_processes
        self.save_confidence_maps = save_confidence_maps
        self.progress_bar = progress_bar

        # Create output directory
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.confidence_maps_dir = self.output_dir / 'confidence_maps'
        if save_confidence_maps:
            self.confidence_maps_dir.mkdir(exist_ok=True)

        logger.info(f"Initialized BatchChannelProcessor with {n_workers} workers")
        logger.info(f"Output directory: {self.output_dir}")

    def process_zarr_dataset(
        self,
        dataset_path: Union[str, Path],
        channel_names: Optional[List[str]] = None,
        cycle_indices: Optional[List[int]] = None,
        channel_indices: Optional[List[int]] = None
    ) -> Dict:
        """Process a Zarr dataset with multiple cycles and channels.

        Args:
            dataset_path: Path to Zarr array
            channel_names: Optional list of channel names (must match dataset dimensions)
            cycle_indices: Optional list of specific cycle indices to process
            channel_indices: Optional list of specific channel indices to process

        Returns:
            Dictionary with processing summary and results
        """
        logger.info(f"Opening Zarr dataset: {dataset_path}")
        z_arr = zarr.open(str(dataset_path), mode='r')

        # Get dataset dimensions
        shape = z_arr.shape
        logger.info(f"Dataset shape: {shape}")

        # Handle different dataset structures
        if len(shape) == 4:
            # Assume (cycles, channels, height, width)
            n_cycles, n_channels, height, width = shape
            logger.info(f"Detected 4D dataset: {n_cycles} cycles, {n_channels} channels")
        elif len(shape) == 3:
            # Assume (channels, height, width) - single cycle
            n_cycles = 1
            n_channels, height, width = shape
            logger.info(f"Detected 3D dataset: {n_channels} channels, single cycle")
        elif len(shape) == 2:
            # Single image
            n_cycles = 1
            n_channels = 1
            height, width = shape
            logger.info(f"Detected 2D dataset: single channel")
        else:
            raise ValueError(f"Unsupported dataset shape: {shape}")

        # Determine which cycles and channels to process
        if cycle_indices is None:
            cycle_indices = list(range(n_cycles))
        if channel_indices is None:
            channel_indices = list(range(n_channels))

        # Generate channel names if not provided
        if channel_names is None:
            channel_names = [f"channel_{i}" for i in range(n_channels)]
        elif len(channel_names) != n_channels:
            warnings.warn(f"Number of channel names ({len(channel_names)}) doesn't match "
                         f"number of channels ({n_channels}). Using default names.")
            channel_names = [f"channel_{i}" for i in range(n_channels)]

        # Prepare tasks
        tasks = []
        for cycle_idx in cycle_indices:
            for chan_idx in channel_indices:
                task_name = f"cycle_{cycle_idx}_{channel_names[chan_idx]}"
                tasks.append({
                    'cycle_idx': cycle_idx,
                    'channel_idx': chan_idx,
                    'channel_name': task_name,
                    'zarr_path': str(dataset_path)
                })

        logger.info(f"Processing {len(tasks)} channels...")

        # Process in parallel
        results = self._process_tasks_parallel(tasks)

        # Aggregate results
        summary = self._aggregate_results(results, n_cycles, n_channels)

        # Save summary
        self._save_summary(summary)

        logger.info("Processing complete!")
        logger.info(f"Auto-approved: {summary['auto_approved_ratio']*100:.1f}%")
        logger.info(f"Requiring review: {len(summary['channels_requiring_review'])}")

        return summary

    def process_directory(
        self,
        image_dir: Union[str, Path],
        file_pattern: str = "*.tif*",
        channel_parser: Optional[Callable[[str], str]] = None
    ) -> Dict:
        """Process all images in a directory.

        Args:
            image_dir: Directory containing images
            file_pattern: Glob pattern to match image files
            channel_parser: Function to extract channel name from filename

        Returns:
            Dictionary with processing summary and results
        """
        image_dir = Path(image_dir)
        image_files = sorted(image_dir.glob(file_pattern))

        if not image_files:
            raise ValueError(f"No files found matching {file_pattern} in {image_dir}")

        logger.info(f"Found {len(image_files)} image files")

        # Prepare tasks
        tasks = []
        for img_path in image_files:
            if channel_parser:
                channel_name = channel_parser(img_path.name)
            else:
                channel_name = img_path.stem

            tasks.append({
                'file_path': str(img_path),
                'channel_name': channel_name
            })

        # Process in parallel
        results = self._process_tasks_parallel(tasks)

        # Aggregate results
        summary = self._aggregate_results(results, 1, len(tasks))

        # Save summary
        self._save_summary(summary)

        logger.info("Processing complete!")
        logger.info(f"Auto-approved: {summary['auto_approved_ratio']*100:.1f}%")
        logger.info(f"Requiring review: {len(summary['channels_requiring_review'])}")

        return summary

    def _process_tasks_parallel(self, tasks: List[Dict]) -> Dict[str, ChannelQualityResult]:
        """Process tasks in parallel.

        Args:
            tasks: List of task dictionaries

        Returns:
            Dictionary mapping channel names to results
        """
        results = {}
        executor_class = ProcessPoolExecutor if self.use_processes else ThreadPoolExecutor

        with executor_class(max_workers=self.n_workers) as executor:
            # Submit all tasks
            future_to_task = {}
            for task in tasks:
                future = executor.submit(self._process_single_task, task)
                future_to_task[future] = task['channel_name']

            # Collect results with progress bar
            iterator = as_completed(future_to_task)
            if self.progress_bar:
                iterator = tqdm(iterator, total=len(tasks), desc="Processing channels")

            for future in iterator:
                channel_name = future_to_task[future]
                try:
                    result = future.result()
                    results[channel_name] = result

                    # Log if review needed
                    if result.requires_review:
                        logger.debug(f"{channel_name}: Review needed "
                                   f"(confidence={result.confidence_score:.3f})")

                except Exception as e:
                    logger.error(f"Error processing {channel_name}: {e}")
                    results[channel_name] = None

        return results

    def _process_single_task(self, task: Dict) -> ChannelQualityResult:
        """Process a single channel assessment task.

        Args:
            task: Task dictionary with channel information

        Returns:
            ChannelQualityResult
        """
        channel_name = task['channel_name']

        # Load channel data
        if 'zarr_path' in task:
            # Load from Zarr
            z_arr = zarr.open(task['zarr_path'], mode='r')
            cycle_idx = task['cycle_idx']
            channel_idx = task['channel_idx']

            # Handle different array structures
            if len(z_arr.shape) == 4:
                channel_data = z_arr[cycle_idx, channel_idx, :, :]
            elif len(z_arr.shape) == 3:
                channel_data = z_arr[channel_idx, :, :]
            else:
                channel_data = z_arr[:, :]

        elif 'file_path' in task:
            # Load from file
            from imageio import imread
            channel_data = imread(task['file_path'])

        else:
            raise ValueError(f"Task must contain either 'zarr_path' or 'file_path': {task}")

        # Process channel
        result = self.assessor.process_image(
            channel_data,
            channel_name,
            metadata=task
        )

        return result

    def _aggregate_results(
        self,
        results: Dict[str, Optional[ChannelQualityResult]],
        n_cycles: int,
        n_channels: int
    ) -> Dict:
        """Aggregate processing results into summary.

        Args:
            results: Dictionary of processing results
            n_cycles: Number of cycles processed
            n_channels: Number of channels per cycle

        Returns:
            Summary dictionary
        """
        channels_requiring_review = []
        auto_approved = []
        failed = []

        for channel_name, result in results.items():
            if result is None:
                failed.append(channel_name)
            elif result.requires_review:
                channels_requiring_review.append(channel_name)
            else:
                auto_approved.append(channel_name)

        total_processed = len([r for r in results.values() if r is not None])
        auto_approved_ratio = len(auto_approved) / total_processed if total_processed > 0 else 0

        summary = {
            'timestamp': datetime.now().isoformat(),
            'total_channels': len(results),
            'n_cycles': n_cycles,
            'n_channels_per_cycle': n_channels,
            'successfully_processed': total_processed,
            'failed': failed,
            'auto_approved': auto_approved,
            'channels_requiring_review': channels_requiring_review,
            'auto_approved_ratio': auto_approved_ratio,
            'channel_details': {}
        }

        # Add detailed results
        for channel_name, result in results.items():
            if result:
                summary['channel_details'][channel_name] = {
                    'confidence_score': result.confidence_score,
                    'requires_review': result.requires_review,
                    'quality_metrics': result.quality_metrics,
                    'n_problem_regions': len(result.problem_regions),
                    'problem_regions': result.problem_regions
                }
            else:
                summary['channel_details'][channel_name] = None

        return summary

    def _save_summary(self, summary: Dict):
        """Save processing summary to JSON file.

        Args:
            summary: Summary dictionary
        """
        output_path = self.output_dir / 'processing_summary.json'
        logger.info(f"Saving summary to {output_path}")

        with open(output_path, 'w') as f:
            json.dump(summary, f, indent=2)

        # Also save a simplified CSV for quick review
        self._save_summary_csv(summary)

    def _save_summary_csv(self, summary: Dict):
        """Save a simplified CSV summary.

        Args:
            summary: Summary dictionary
        """
        import csv

        output_path = self.output_dir / 'processing_summary.csv'

        with open(output_path, 'w', newline='') as f:
            writer = csv.writer(f)
            writer.writerow([
                'Channel', 'Confidence', 'Requires Review',
                'Mean Intensity', 'SNR', 'Low Confidence Ratio', 'N Problem Regions'
            ])

            for channel_name, details in summary['channel_details'].items():
                if details:
                    metrics = details['quality_metrics']
                    writer.writerow([
                        channel_name,
                        f"{details['confidence_score']:.4f}",
                        details['requires_review'],
                        f"{metrics.get('mean_intensity', 0):.2f}",
                        f"{metrics.get('snr', 0):.2f}",
                        f"{metrics.get('low_confidence_ratio', 0):.4f}",
                        details['n_problem_regions']
                    ])
                else:
                    writer.writerow([channel_name, 'FAILED', 'TRUE', '', '', '', ''])

        logger.info(f"Saved CSV summary to {output_path}")


class StreamingBatchProcessor(BatchChannelProcessor):
    """Batch processor that processes data in streaming fashion for memory efficiency.

    This is useful for very large datasets that don't fit in memory.
    """

    def __init__(self, *args, chunk_size: int = 10, **kwargs):
        """Initialize streaming processor.

        Args:
            chunk_size: Number of channels to process before writing intermediate results
            *args, **kwargs: Arguments for BatchChannelProcessor
        """
        super().__init__(*args, **kwargs)
        self.chunk_size = chunk_size

    def process_zarr_dataset(self, *args, **kwargs) -> Dict:
        """Process Zarr dataset in chunks."""
        # This would implement chunked processing
        # For now, delegate to parent implementation
        logger.info(f"Streaming processing with chunk_size={self.chunk_size}")
        return super().process_zarr_dataset(*args, **kwargs)
