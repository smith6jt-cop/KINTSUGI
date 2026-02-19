"""
Unified artifact scanning for KINTSUGI projects.

This module provides a single, integrated workflow for detecting artifacts
in microscopy image stacks. It works directly with the project structure
and supports both TIFF and Zarr data sources.

Usage:
    from kintsugi.project import KintsugiProject
    from kintsugi.qc import ArtifactScanner

    project = KintsugiProject.load("/path/to/project")
    scanner = ArtifactScanner(project)

    # Scan all raw data
    report = scanner.scan_raw_data()

    # Get affected images
    for item in report.affected_items:
        print(f"{item.cycle}/{item.channel}/z{item.z_plane}: {item.severity}")

    # Save report for later use
    report.save(project.paths.meta / "artifact_report.json")
"""

from __future__ import annotations

import json
import re
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import TYPE_CHECKING, Any, Literal

import numpy as np

if TYPE_CHECKING:
    from kintsugi.project import KintsugiProject

# Lazy imports for optional dependencies
_tifffile = None
_zarr = None


def _get_tifffile():
    global _tifffile
    if _tifffile is None:
        import tifffile

        _tifffile = tifffile
    return _tifffile


def _get_zarr():
    global _zarr
    if _zarr is None:
        import zarr

        _zarr = zarr
    return _zarr


@dataclass
class ArtifactItem:
    """A single detected artifact."""

    cycle: str
    channel: int
    z_plane: int
    tile: str | None = None
    severity: Literal["mild", "moderate", "severe"] = "mild"
    score: float = 0.0
    metrics: dict[str, float] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {
            "cycle": self.cycle,
            "channel": self.channel,
            "z_plane": self.z_plane,
            "tile": self.tile,
            "severity": self.severity,
            "score": self.score,
            "metrics": self.metrics,
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ArtifactItem:
        return cls(**data)


@dataclass
class ChannelScanResult:
    """Scan result for a single channel across all z-planes."""

    cycle: str
    channel: int
    n_planes: int
    n_tiles: int
    affected_planes: list[int] = field(default_factory=list)
    worst_severity: Literal["none", "mild", "moderate", "severe"] = "none"
    baseline_metric: float = 0.0
    threshold: float = 0.0
    plane_metrics: dict[int, float] = field(default_factory=dict)

    @property
    def has_artifacts(self) -> bool:
        return len(self.affected_planes) > 0

    @property
    def artifact_ratio(self) -> float:
        """Fraction of planes with artifacts."""
        return len(self.affected_planes) / self.n_planes if self.n_planes > 0 else 0.0


@dataclass
class ArtifactReport:
    """Complete artifact report for a project or dataset."""

    created: str = field(default_factory=lambda: datetime.now().isoformat())
    project_name: str = ""
    data_source: str = ""  # "raw", "stitched", etc.
    scan_parameters: dict[str, Any] = field(default_factory=dict)

    # Summary statistics
    total_cycles: int = 0
    total_channels: int = 0
    total_planes_scanned: int = 0
    total_artifacts_found: int = 0

    # Detailed results
    channel_results: list[ChannelScanResult] = field(default_factory=list)
    affected_items: list[ArtifactItem] = field(default_factory=list)

    @property
    def has_artifacts(self) -> bool:
        return self.total_artifacts_found > 0

    @property
    def worst_severity(self) -> Literal["none", "mild", "moderate", "severe"]:
        if not self.affected_items:
            return "none"
        severity_order = {"mild": 1, "moderate": 2, "severe": 3}
        worst = max(self.affected_items, key=lambda x: severity_order.get(x.severity, 0))
        return worst.severity

    def summary(self) -> str:
        """Generate human-readable summary."""
        lines = [
            f"Artifact Scan Report: {self.project_name}",
            f"{'=' * 50}",
            f"Source: {self.data_source}",
            f"Scanned: {self.total_cycles} cycles, {self.total_channels} channels",
            f"Total planes: {self.total_planes_scanned}",
            "",
        ]

        if not self.has_artifacts:
            lines.append("No artifacts detected.")
        else:
            lines.append(f"Artifacts found: {self.total_artifacts_found}")
            lines.append(f"Worst severity: {self.worst_severity}")
            lines.append("")

            # Group by cycle/channel
            by_cycle = {}
            for item in self.affected_items:
                key = f"{item.cycle}/CH{item.channel}"
                if key not in by_cycle:
                    by_cycle[key] = []
                by_cycle[key].append(item.z_plane)

            lines.append("Affected locations:")
            for key, planes in sorted(by_cycle.items()):
                lines.append(f"  {key}: z-planes {planes}")

        return "\n".join(lines)

    def save(self, path: str | Path) -> None:
        """Save report to JSON file."""
        path = Path(path)
        data = {
            "created": self.created,
            "project_name": self.project_name,
            "data_source": self.data_source,
            "scan_parameters": self.scan_parameters,
            "summary": {
                "total_cycles": self.total_cycles,
                "total_channels": self.total_channels,
                "total_planes_scanned": self.total_planes_scanned,
                "total_artifacts_found": self.total_artifacts_found,
                "worst_severity": self.worst_severity,
            },
            "affected_items": [item.to_dict() for item in self.affected_items],
        }
        with open(path, "w") as f:
            json.dump(data, f, indent=2)

    @classmethod
    def load(cls, path: str | Path) -> ArtifactReport:
        """Load report from JSON file."""
        path = Path(path)
        with open(path) as f:
            data = json.load(f)

        report = cls(
            created=data["created"],
            project_name=data["project_name"],
            data_source=data["data_source"],
            scan_parameters=data.get("scan_parameters", {}),
            total_cycles=data["summary"]["total_cycles"],
            total_channels=data["summary"]["total_channels"],
            total_planes_scanned=data["summary"]["total_planes_scanned"],
            total_artifacts_found=data["summary"]["total_artifacts_found"],
        )
        report.affected_items = [ArtifactItem.from_dict(item) for item in data["affected_items"]]
        return report


class ArtifactScanner:
    """
    Unified artifact scanner for KINTSUGI projects.

    This scanner uses statistical comparison across z-planes to detect
    artifacts, which is robust to dataset-specific baselines.

    Parameters
    ----------
    project : KintsugiProject
        The project to scan
    threshold_sigma : float
        Number of standard deviations above baseline to flag as artifact.
        Default 2.0 catches significant outliers without too many false positives.

    Example
    -------
    >>> from kintsugi.project import KintsugiProject
    >>> from kintsugi.qc import ArtifactScanner
    >>>
    >>> project = KintsugiProject.load("/path/to/project")
    >>> scanner = ArtifactScanner(project)
    >>> report = scanner.scan_raw_data()
    >>> print(report.summary())
    """

    def __init__(
        self,
        project: KintsugiProject,
        threshold_sigma: float = 2.0,
    ):
        self.project = project
        self.threshold_sigma = threshold_sigma

    def scan_raw_data(
        self,
        cycles: list[str] | None = None,
        channels: list[int] | None = None,
        progress_callback: callable | None = None,
    ) -> ArtifactReport:
        """
        Scan raw data for artifacts.

        Parameters
        ----------
        cycles : list[str], optional
            Specific cycles to scan (default: all)
        channels : list[int], optional
            Specific channels to scan (default: all)
        progress_callback : callable, optional
            Function called with (current, total, message) for progress updates

        Returns
        -------
        ArtifactReport
            Complete scan report
        """
        report = ArtifactReport(
            project_name=self.project.config.name,
            data_source="raw",
            scan_parameters={"threshold_sigma": self.threshold_sigma},
        )

        # Discover cycles
        if cycles is None:
            cycle_dirs = self._find_cycle_dirs(self.project.paths.raw)
        else:
            cycle_dirs = [
                self.project.paths.raw / c for c in cycles if (self.project.paths.raw / c).exists()
            ]

        if not cycle_dirs:
            return report

        report.total_cycles = len(cycle_dirs)

        # Scan each cycle
        for cycle_dir in cycle_dirs:
            cycle_name = cycle_dir.name
            channel_groups = self._group_by_channel(cycle_dir)

            for channel, files_by_z in channel_groups.items():
                if channels is not None and channel not in channels:
                    continue

                report.total_channels += 1

                # Load z-stack (sample one tile per z-plane)
                z_planes = sorted(files_by_z.keys())
                if not z_planes:
                    continue

                # Load representative images for each z-plane
                z_stack = self._load_zstack_sample(files_by_z, z_planes)
                if z_stack is None or len(z_stack) < 3:
                    continue

                report.total_planes_scanned += len(z_planes)

                # Scan this z-stack
                result = self._scan_zstack(z_stack, cycle_name, channel, z_planes)
                report.channel_results.append(result)

                # Add affected items
                for z_idx in result.affected_planes:
                    z_plane = z_planes[z_idx]
                    deviation = result.plane_metrics[z_idx] - result.baseline_metric
                    deviation /= max(0.0001, result.threshold - result.baseline_metric)

                    if deviation > 2:
                        severity = "severe"
                    elif deviation > 1:
                        severity = "moderate"
                    else:
                        severity = "mild"

                    item = ArtifactItem(
                        cycle=cycle_name,
                        channel=channel,
                        z_plane=z_plane,
                        severity=severity,
                        score=min(1.0, deviation / 3.0),
                        metrics={"normalized_row_diff": result.plane_metrics[z_idx]},
                    )
                    report.affected_items.append(item)
                    report.total_artifacts_found += 1

                if progress_callback:
                    progress_callback(
                        report.total_channels,
                        report.total_cycles * 4,  # Estimate 4 channels per cycle
                        f"Scanned {cycle_name}/CH{channel}",
                    )

        return report

    def scan_zarr(
        self,
        stage: str = "raw",
        cycles: list[str] | None = None,
        progress_callback: callable | None = None,
    ) -> ArtifactReport:
        """
        Scan Zarr store for artifacts.

        Parameters
        ----------
        stage : str
            Processing stage to scan ("raw", "stitched", "deconvolved", etc.)
        cycles : list[str], optional
            Specific cycles to scan
        progress_callback : callable, optional
            Progress callback function

        Returns
        -------
        ArtifactReport
            Scan report
        """
        zarr_path = self.project.paths.processed / f"{self.project.root.name}.zarr"

        if not zarr_path.exists():
            return ArtifactReport(
                project_name=self.project.config.name,
                data_source=f"zarr/{stage}",
            )

        zarr = _get_zarr()
        store = zarr.open(str(zarr_path), mode="r")

        report = ArtifactReport(
            project_name=self.project.config.name,
            data_source=f"zarr/{stage}",
            scan_parameters={"threshold_sigma": self.threshold_sigma},
        )

        # Iterate over cycles in zarr
        for cycle_name in store.keys():
            if cycles is not None and cycle_name not in cycles:
                continue

            cycle_group = store[cycle_name]
            if stage not in cycle_group:
                continue

            stage_group = cycle_group[stage]
            report.total_cycles += 1

            # Scan each channel array
            for channel_name in stage_group.keys():
                channel = int(re.search(r"\d+", channel_name).group())
                arr = stage_group[channel_name]

                report.total_channels += 1

                # Load z-stack
                if arr.ndim >= 3:
                    z_stack = np.array(arr[:])  # Load all z-planes
                    z_planes = list(range(z_stack.shape[0]))

                    report.total_planes_scanned += len(z_planes)

                    # Scan
                    result = self._scan_zstack(z_stack, cycle_name, channel, z_planes)
                    report.channel_results.append(result)

                    for z_idx in result.affected_planes:
                        deviation = result.plane_metrics[z_idx] - result.baseline_metric
                        deviation /= max(0.0001, result.threshold - result.baseline_metric)

                        severity = (
                            "severe" if deviation > 2 else ("moderate" if deviation > 1 else "mild")
                        )

                        report.affected_items.append(
                            ArtifactItem(
                                cycle=cycle_name,
                                channel=channel,
                                z_plane=z_planes[z_idx],
                                severity=severity,
                                score=min(1.0, deviation / 3.0),
                            )
                        )
                        report.total_artifacts_found += 1

        return report

    def _scan_zstack(
        self,
        z_stack: np.ndarray,
        cycle: str,
        channel: int,
        z_indices: list[int],
    ) -> ChannelScanResult:
        """
        Scan a z-stack for artifacts using comparative analysis.

        Uses normalized row-difference metric compared across z-planes
        to find statistical outliers.
        """
        nz = z_stack.shape[0]
        stack = z_stack.astype(np.float32)

        # Compute row means for each z-plane
        # Shape: (nz, height) after mean over width
        row_means = np.mean(stack, axis=2) if stack.ndim == 3 else np.mean(stack, axis=-1)

        # Compute row-diff std for each plane (stripe metric)
        row_diffs = np.diff(row_means, axis=1)
        row_diff_stds = np.std(row_diffs, axis=1)

        # Normalize by signal std
        signal_stds = np.std(stack.reshape(nz, -1), axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            metrics = row_diff_stds / (signal_stds + 1e-10)

        # Compute baseline statistics
        baseline_mean = np.mean(metrics)
        baseline_std = np.std(metrics)

        # Flag outliers
        threshold = baseline_mean + self.threshold_sigma * baseline_std
        affected_indices = np.where(metrics > threshold)[0].tolist()

        # Determine worst severity
        if not affected_indices:
            worst_severity = "none"
        else:
            max_deviation = max(
                (metrics[i] - baseline_mean) / (baseline_std + 1e-10) for i in affected_indices
            )
            if max_deviation > 4:
                worst_severity = "severe"
            elif max_deviation > 3:
                worst_severity = "moderate"
            else:
                worst_severity = "mild"

        return ChannelScanResult(
            cycle=cycle,
            channel=channel,
            n_planes=nz,
            n_tiles=1,  # Sampled
            affected_planes=affected_indices,
            worst_severity=worst_severity,
            baseline_metric=float(baseline_mean),
            threshold=float(threshold),
            plane_metrics={i: float(metrics[i]) for i in range(nz)},
        )

    def _find_cycle_dirs(self, raw_dir: Path) -> list[Path]:
        """Find cycle directories in raw data folder."""
        if not raw_dir.exists():
            return []

        cycle_pattern = re.compile(r"^cyc\d+$|^cycle\d+$", re.IGNORECASE)
        cycles = []

        for item in raw_dir.iterdir():
            if item.is_dir() and cycle_pattern.match(item.name):
                cycles.append(item)

        return sorted(cycles, key=lambda p: p.name)

    def _group_by_channel(self, cycle_dir: Path) -> dict[int, dict[int, list[Path]]]:
        """
        Group TIFF files by channel and z-plane.

        Returns dict[channel, dict[z_plane, list[files]]]
        """
        # Pattern: *_Z###_CH#.tif
        pattern = re.compile(r"_Z(\d+)_CH(\d+)\.tif", re.IGNORECASE)

        groups: dict[int, dict[int, list[Path]]] = {}

        for f in cycle_dir.glob("*.tif"):
            match = pattern.search(f.name)
            if match:
                z_plane = int(match.group(1))
                channel = int(match.group(2))

                if channel not in groups:
                    groups[channel] = {}
                if z_plane not in groups[channel]:
                    groups[channel][z_plane] = []

                groups[channel][z_plane].append(f)

        return groups

    def _load_zstack_sample(
        self,
        files_by_z: dict[int, list[Path]],
        z_planes: list[int],
    ) -> np.ndarray | None:
        """Load one sample tile per z-plane to create a z-stack."""
        tifffile = _get_tifffile()

        images = []
        for z in z_planes:
            files = files_by_z.get(z, [])
            if not files:
                continue
            # Load first tile for this z-plane
            try:
                img = tifffile.imread(str(files[0]))
                images.append(img)
            except Exception:
                continue

        if not images:
            return None

        return np.array(images)


def detect_general_artifacts(
    image: np.ndarray,
    downscale_factor: int = 4,
    erosion_size: int = 3,
    h_tolerance: float | None = None,
    flood_tolerance: float | None = None,
    min_artifact_size: int = 100,
) -> np.ndarray:
    """
    Detect general imaging artifacts using flood-fill from intensity peaks.

    Inspired by CyLinter's ``artifact_detector_v3`` algorithm. This approach
    is complementary to KINTSUGI's existing stripe artifact detection — it
    detects arbitrary high-contrast regions that may correspond to tissue
    folding, air bubbles, bright debris, or other imaging artifacts.

    Algorithm:
    1. Downscale image for computational efficiency
    2. Apply erosion + local contrast enhancement to highlight artifact boundaries
    3. Use h-maxima transform to find local intensity peaks (artifact centers)
    4. Apply flood-fill from each peak with automatically determined tolerance
    5. Return a binary artifact mask

    Parameters
    ----------
    image : np.ndarray
        2D grayscale image (uint8 or uint16).
    downscale_factor : int
        Factor to downscale by for detection. Results are upscaled to
        original resolution. Default 4.
    erosion_size : int
        Size of erosion structuring element. Default 3.
    h_tolerance : float, optional
        Height parameter for h-maxima transform. If None, auto-computed
        as 0.2 * (max - median) of the processed image.
    flood_tolerance : float, optional
        Tolerance for flood-fill. If None, auto-computed as
        0.15 * (max - median) of the processed image.
    min_artifact_size : int
        Minimum size (in pixels at downscaled resolution) for a region
        to be considered an artifact. Default 100.

    Returns
    -------
    np.ndarray
        Binary mask (same shape as input) where True = artifact region.
        Can be used for cell exclusion or as a ``Kview2`` overlay.
    """
    from scipy import ndimage as ndi
    from skimage.morphology import reconstruction
    from skimage.transform import resize

    img = np.asarray(image, dtype=np.float32)
    original_shape = img.shape

    # Step 1: Downscale
    if downscale_factor > 1:
        new_shape = (
            img.shape[0] // downscale_factor,
            img.shape[1] // downscale_factor,
        )
        small = resize(img, new_shape, preserve_range=True, anti_aliasing=True)
    else:
        small = img.copy()

    # Step 2: Erosion to suppress noise + highlight artifact edges
    struct = np.ones((erosion_size, erosion_size))
    eroded = ndi.grey_erosion(small, footprint=struct)

    # Check if image has sufficient dynamic range for artifact detection
    img_global_range = small.max() - small.min()
    if img_global_range < 50:  # Nearly uniform image
        return np.zeros(original_shape, dtype=bool)

    # Local contrast enhancement (difference of Gaussians)
    smooth_low = ndi.gaussian_filter(eroded, sigma=2)
    smooth_high = ndi.gaussian_filter(eroded, sigma=10)
    enhanced = np.clip(smooth_low - smooth_high, 0, None)

    # Step 3: H-maxima transform to find peaks
    img_max = enhanced.max()
    img_median = np.median(enhanced)
    img_range = img_max - img_median

    if img_range < 1e-6:
        # Uniform image — no artifacts
        return np.zeros(original_shape, dtype=bool)

    if h_tolerance is None:
        h_tolerance = 0.2 * img_range
    if flood_tolerance is None:
        flood_tolerance = 0.15 * img_range

    # H-maxima: suppress peaks lower than h_tolerance
    seed = enhanced - h_tolerance
    seed = np.clip(seed, enhanced.min(), None)
    dilated = reconstruction(seed, enhanced, method="dilation")
    h_maxima = enhanced - dilated

    # Find peak locations
    peak_mask = h_maxima > (h_tolerance * 0.5)
    labeled_peaks, n_peaks = ndi.label(peak_mask)

    if n_peaks == 0:
        return np.zeros(original_shape, dtype=bool)

    # Step 4: Flood-fill from each peak
    artifact_mask = np.zeros(small.shape, dtype=bool)

    for peak_id in range(1, n_peaks + 1):
        peak_coords = np.where(labeled_peaks == peak_id)
        # Use the brightest pixel in this peak region as seed
        peak_values = enhanced[peak_coords]
        brightest_idx = np.argmax(peak_values)
        seed_y = peak_coords[0][brightest_idx]
        seed_x = peak_coords[1][brightest_idx]
        seed_value = enhanced[seed_y, seed_x]

        # Simple tolerance-based flood fill
        flood_low = seed_value - flood_tolerance
        flood_high = seed_value + flood_tolerance

        # Use connected component analysis on thresholded region
        region = (enhanced >= flood_low) & (enhanced <= flood_high)
        labeled_region, _ = ndi.label(region)
        seed_label = labeled_region[seed_y, seed_x]

        if seed_label > 0:
            component = labeled_region == seed_label
            if np.sum(component) >= min_artifact_size:
                artifact_mask |= component

    # Step 5: Upscale mask to original resolution
    if downscale_factor > 1:
        full_mask = resize(
            artifact_mask.astype(np.float32),
            original_shape,
            preserve_range=True,
            order=0,  # Nearest-neighbor to preserve binary boundaries
        )
        return full_mask > 0.5
    else:
        return artifact_mask


def scan_project_artifacts(
    project: KintsugiProject,
    source: str = "raw",
    threshold_sigma: float = 2.0,
    save_report: bool = True,
) -> ArtifactReport:
    """
    Convenience function to scan a project for artifacts.

    Parameters
    ----------
    project : KintsugiProject
        Project to scan
    source : str
        Data source to scan ("raw" or stage name for zarr)
    threshold_sigma : float
        Detection threshold
    save_report : bool
        If True, save report to project meta directory

    Returns
    -------
    ArtifactReport
        Scan results
    """
    scanner = ArtifactScanner(project, threshold_sigma=threshold_sigma)

    if source == "raw":
        report = scanner.scan_raw_data()
    else:
        report = scanner.scan_zarr(stage=source)

    if save_report:
        report_path = project.paths.meta / f"artifact_report_{source}.json"
        report_path.parent.mkdir(parents=True, exist_ok=True)
        report.save(report_path)
        print(f"Report saved to: {report_path}")

    return report
