"""
KINTSUGI Project Management Framework

Provides a clean separation between:
- Repository code (KINTSUGI installation)
- Project data (raw images, processed outputs)
- Working notebooks (copied from templates)

Usage:
    from kintsugi.project import KintsugiProject

    # Create or load a project
    project = KintsugiProject.create("C:/Projects/MyExperiment")
    # or
    project = KintsugiProject.load("C:/Projects/MyExperiment")

    # Access directories
    project.raw_dir      # Raw input images
    project.processed_dir # Processed outputs
    project.notebooks_dir # Working notebooks
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import sys
import warnings
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, ClassVar, Literal

# Project configuration filename
PROJECT_CONFIG_FILE = "kintsugi_project.json"
PROJECT_VERSION = "1.0.0"


def _resolve_kintsugi_executable() -> str:
    """Resolve absolute path to the kintsugi CLI executable.

    Tries shutil.which first (works if env is activated), then falls back
    to the sibling of the current Python interpreter (works from the same
    conda env even when PATH doesn't include the env's bin/).
    """
    path = shutil.which("kintsugi")
    if path:
        return str(Path(path).resolve())

    # Fall back to a sibling of the current Python interpreter. Layouts differ
    # by platform: POSIX venv/conda put console scripts in the same bin/ as
    # python, while Windows uses a Scripts/ dir and an .exe suffix.
    interp_dir = Path(sys.executable).parent
    candidates = [
        interp_dir / "kintsugi",
        interp_dir / "kintsugi.exe",
        interp_dir / "Scripts" / "kintsugi.exe",
        interp_dir.parent / "Scripts" / "kintsugi.exe",
    ]
    for candidate in candidates:
        if candidate.exists():
            return str(candidate.resolve())

    # Last resort: bare command name (resolved from PATH at launch time).
    return "kintsugi"


# Common image file extensions
IMAGE_EXTENSIONS = {
    ".tif",
    ".tiff",
    ".ome.tif",
    ".ome.tiff",
    ".png",
    ".jpg",
    ".jpeg",
    ".nd2",
    ".czi",
    ".lif",
    ".ims",
    ".vsi",
    ".svs",
    ".qptiff",
}


@dataclass
class ExistingDataReport:
    """Report of existing data found in a directory."""

    has_data: bool = False
    image_files: list[Path] = field(default_factory=list)
    image_count: int = 0
    total_size_mb: float = 0.0
    cycle_folders: list[str] = field(default_factory=list)
    other_files: list[Path] = field(default_factory=list)
    metadata_samples: list[dict[str, Any]] = field(default_factory=list)
    filename_patterns: list[str] = field(default_factory=list)

    # Separate tracking for raw vs processed data
    has_raw_data: bool = False
    has_processed_data: bool = False
    raw_image_count: int = 0
    raw_size_mb: float = 0.0
    raw_cycle_folders: list[str] = field(default_factory=list)
    processed_stages: dict[str, int] = field(default_factory=dict)  # stage_name -> file_count
    processed_size_mb: float = 0.0

    def summary(self) -> str:
        """Generate a human-readable summary."""
        lines = []
        if not self.has_data:
            lines.append("No existing data found.")
            return "\n".join(lines)

        # Report raw data
        if self.has_raw_data:
            lines.append(f"Raw data: {self.raw_image_count} images ({self.raw_size_mb:.1f} MB)")
            if self.raw_cycle_folders:
                lines.append(f"  Cycles: {', '.join(self.raw_cycle_folders)}")

        # Report processed data
        if self.has_processed_data:
            lines.append(f"Processed data: {self.processed_size_mb:.1f} MB")
            for stage, count in self.processed_stages.items():
                lines.append(f"  {stage}: {count} files")

        # Legacy summary for backwards compatibility
        if not self.has_raw_data and not self.has_processed_data:
            lines.append(
                f"Found {self.image_count} image files ({self.total_size_mb:.1f} MB total)"
            )
            if self.cycle_folders:
                lines.append(f"Cycle folders detected: {', '.join(self.cycle_folders)}")

        if self.filename_patterns:
            lines.append(f"Filename patterns: {', '.join(self.filename_patterns[:5])}")

        if self.other_files:
            lines.append(f"Other files: {len(self.other_files)}")

        return "\n".join(lines)


@dataclass
class ImageMetadata:
    """Metadata extracted from an image file."""

    path: Path
    filename: str
    size_mb: float
    dimensions: tuple[int, ...] | None = None
    dtype: str | None = None
    channels: int | None = None
    channel_names: list[str] = field(default_factory=list)
    pixel_size: float | None = None
    pixel_unit: str | None = None
    software: str | None = None
    datetime: str | None = None
    description: str | None = None
    is_ome: bool = False
    extra: dict[str, Any] = field(default_factory=dict)

    def summary(self) -> str:
        """Generate a human-readable summary."""
        lines = [f"File: {self.filename}"]
        if self.dimensions:
            lines.append(f"  Dimensions: {self.dimensions}")
        if self.dtype:
            lines.append(f"  Data type: {self.dtype}")
        if self.channels:
            lines.append(f"  Channels: {self.channels}")
        if self.channel_names:
            lines.append(f"  Channel names: {', '.join(self.channel_names[:5])}")
        if self.pixel_size:
            lines.append(f"  Pixel size: {self.pixel_size} {self.pixel_unit or 'units'}")
        if self.software:
            lines.append(f"  Software: {self.software}")
        if self.datetime:
            lines.append(f"  Acquired: {self.datetime}")
        return "\n".join(lines)


# Processing stage names in order
PROCESSING_STAGES = [
    "stitched",
    "deconvolved",
    "edf",
    "registered",
    "signal_isolated",
    "segmented",
    "analysis",
]

# Experiment metadata filename
EXPERIMENT_CONFIG_FILE = "experiment.json"

# Known CODEX filter sets: excitation → (excitation, emission) in nm
# Used to convert CODEX flat wavelength lists to proper (ex, em) pairs
_CODEX_FILTER_SETS: dict[int, tuple[float, float]] = {
    358: (358.0, 461.0),  # DAPI
    488: (488.0, 525.0),  # FITC / Alexa 488
    550: (560.0, 575.0),  # TRITC / Cy3
    650: (648.0, 668.0),  # Cy5
    750: (753.0, 775.0),  # Cy7
}


def _excitation_to_filter_pair(ex_nm: float) -> tuple[float, float]:
    """Map a CODEX excitation wavelength to (excitation, emission) pair.

    Uses known filter sets with a ±10 nm tolerance for matching.
    Falls back to the exact value for both if no known filter matches.
    """
    ex_int = round(ex_nm)
    # Exact match first
    if ex_int in _CODEX_FILTER_SETS:
        return _CODEX_FILTER_SETS[ex_int]
    # Fuzzy match within ±10 nm
    for ref_ex, pair in _CODEX_FILTER_SETS.items():
        if abs(ex_int - ref_ex) <= 10:
            return pair
    # Unknown filter - return excitation only with warning
    import warnings

    warnings.warn(
        f"Unknown CODEX excitation wavelength {ex_nm} nm; "
        f"emission wavelength unknown. Add to _CODEX_FILTER_SETS in project.py.",
        stacklevel=3,
    )
    return (ex_nm, ex_nm)


@dataclass
class ExperimentConfig:
    """
    Experiment metadata and microscope parameters.

    Stored in /meta/experiment.json and used by processing scripts.
    """

    # Tile grid configuration
    tile_rows: int = 5
    tile_cols: int = 5
    tile_overlap: float = 0.3  # Fraction (0.3 = 30%)

    # Pixel dimensions (nanometers)
    xy_pixel_size: float = 377.0  # nm per pixel XY
    z_step_size: float = 1500.0  # nm per z-slice

    # Optical parameters
    numerical_aperture: float = 0.75
    tissue_refractive_index: float = 1.44
    immersion_medium: str = "air"  # air, water, oil

    # Wavelengths per channel: {channel_num: (excitation_nm, emission_nm)}
    wavelengths: dict[int, tuple[float, float]] = field(
        default_factory=lambda: {
            1: (358.0, 461.0),  # DAPI
            2: (753.0, 775.0),  # Cy7
            3: (560.0, 575.0),  # TRITC
            4: (648.0, 668.0),  # Cy5
        }
    )

    # Data organization
    channels_per_cycle: int = 4
    n_cycles: int | None = None  # Auto-detected if None
    n_zplanes: int | None = None  # Auto-detected if None

    # Stitching parameters
    ncc_threshold: float = 0.078
    pou: float = 0.5  # Percentage of overlap uncertainty

    # Deconvolution parameters
    decon_iterations: int = 25
    decon_damping: float = 0.0
    decon_stop_criterion: float = 5.0

    # Metadata
    microscope: str = ""
    acquisition_date: str = ""
    notes: str = ""

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary for JSON serialization."""
        data = asdict(self)
        # Convert wavelengths dict keys to strings for JSON
        data["wavelengths"] = {str(k): v for k, v in self.wavelengths.items()}
        return data

    # CODEX field name -> KINTSUGI field name
    _CODEX_FIELD_MAP: ClassVar[dict[str, str]] = {
        "numCycles": "n_cycles",
        "numZPlanes": "n_zplanes",
        "regionHeight": "tile_rows",
        "regionWidth": "tile_cols",
        "numChannels": "channels_per_cycle",
        "xyResolution": "xy_pixel_size",
        "zPitch": "z_step_size",
        "aperture": "numerical_aperture",
        "tileOverlapX": "tile_overlap",
    }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ExperimentConfig:
        """Create from dictionary (handles both KINTSUGI and CODEX field names)."""
        # Translate CODEX field names to KINTSUGI equivalents
        for codex_key, kintsugi_key in cls._CODEX_FIELD_MAP.items():
            if codex_key in data and kintsugi_key not in data:
                data[kintsugi_key] = data[codex_key]

        # Warn if CODEX tileOverlapX != tileOverlapY (asymmetric overlap)
        if "tileOverlapX" in data and "tileOverlapY" in data:
            ox, oy = data["tileOverlapX"], data["tileOverlapY"]
            if ox != oy:
                warnings.warn(
                    f"Asymmetric tile overlap: tileOverlapX={ox}, tileOverlapY={oy}. "
                    f"Using tileOverlapX={ox} for tile_overlap.",
                    stacklevel=2,
                )

        # Convert wavelengths keys back to integers
        if "wavelengths" in data:
            wl = data["wavelengths"]
            if isinstance(wl, dict):
                data["wavelengths"] = {int(k): tuple(v) for k, v in wl.items()}
            elif isinstance(wl, list):
                if wl and isinstance(wl[0], (list, tuple)):
                    # List-of-pairs format: [[1, [358.0, 461.0]], ...]
                    data["wavelengths"] = {int(pair[0]): tuple(pair[1]) for pair in wl}
                elif wl and isinstance(wl[0], (int, float)):
                    # CODEX flat list of excitation wavelengths: [358, 750, 550, 650]
                    # Map to (excitation, emission) pairs using known filter sets
                    data["wavelengths"] = {
                        i + 1: _excitation_to_filter_pair(float(w)) for i, w in enumerate(wl)
                    }
                else:
                    del data["wavelengths"]
            else:
                del data["wavelengths"]  # Invalid format, use default
        # Filter to valid fields only
        valid_fields = {f.name for f in cls.__dataclass_fields__.values()}
        filtered_data = {k: v for k, v in data.items() if k in valid_fields}
        return cls(**filtered_data)

    def to_wavelengths_string(self) -> str:
        """Convert wavelengths to config.sh format: '1:358:461,2:753:775,...'"""
        parts = []
        for ch in sorted(self.wavelengths.keys()):
            ex, em = self.wavelengths[ch]
            parts.append(f"{ch}:{int(ex)}:{int(em)}")
        return ",".join(parts)

    @classmethod
    def from_wavelengths_string(cls, wl_str: str) -> dict[int, tuple[float, float]]:
        """Parse wavelengths from config.sh format."""
        wavelengths = {}
        for item in wl_str.split(","):
            parts = item.split(":")
            if len(parts) == 3:
                ch, ex, em = int(parts[0]), float(parts[1]), float(parts[2])
                wavelengths[ch] = (ex, em)
        return wavelengths


@dataclass
class ProcessedStageInfo:
    """Information about data in a single processing stage."""

    name: str
    path: Path
    exists: bool = False
    file_count: int = 0
    total_size_mb: float = 0.0
    cycle_folders: list[str] = field(default_factory=list)
    sample_files: list[str] = field(default_factory=list)

    def summary(self) -> str:
        """Generate a human-readable summary."""
        if not self.exists or self.file_count == 0:
            return f"  {self.name}: (empty)"
        cycles_str = f", cycles: {', '.join(self.cycle_folders)}" if self.cycle_folders else ""
        return f"  {self.name}: {self.file_count} files ({self.total_size_mb:.1f} MB){cycles_str}"


@dataclass
class ProcessedDataReport:
    """Report of existing data in processed subdirectories."""

    stages: dict[str, ProcessedStageInfo] = field(default_factory=dict)
    total_files: int = 0
    total_size_mb: float = 0.0

    @property
    def has_data(self) -> bool:
        """Check if any processing stage has data."""
        return self.total_files > 0

    @property
    def stages_with_data(self) -> list[str]:
        """Get list of stage names that have data."""
        return [name for name, info in self.stages.items() if info.file_count > 0]

    def summary(self) -> str:
        """Generate a human-readable summary."""
        lines = []
        if not self.has_data:
            lines.append("No processed data found.")
            return "\n".join(lines)

        lines.append(
            f"Found {self.total_files} files ({self.total_size_mb:.1f} MB) in processed directories:"
        )
        for name in PROCESSING_STAGES:
            if name in self.stages:
                lines.append(self.stages[name].summary())
        return "\n".join(lines)


def scan_processed_data(project_paths: ProjectPaths) -> ProcessedDataReport:
    """
    Scan processed subdirectories for existing data (TIFF and Zarr).

    Parameters
    ----------
    project_paths : ProjectPaths
        Project paths object with processed directory paths

    Returns
    -------
    ProcessedDataReport
        Report of existing data in each processing stage
    """
    import re

    report = ProcessedDataReport()

    stage_paths = {
        "stitched": project_paths.stitched,
        "deconvolved": project_paths.deconvolved,
        "edf": project_paths.edf,
        "registered": project_paths.registered,
        "signal_isolated": project_paths.signal_isolated,
        "segmented": project_paths.segmented,
        "analysis": project_paths.analysis,
    }

    cycle_pattern = re.compile(r"cyc\d+", re.IGNORECASE)

    # Find zarr file if it exists
    zarr_path = project_paths.processed / f"{project_paths.root.name}.zarr"

    for stage_name, stage_path in stage_paths.items():
        info = ProcessedStageInfo(name=stage_name, path=stage_path)

        tiff_exists = stage_path.exists()
        zarr_stage_exists = False
        files = []
        cycles = set()

        # Scan TIFF directory for files
        if tiff_exists:
            for item in stage_path.rglob("*"):
                if item.is_file():
                    files.append(item)
                    info.total_size_mb += item.stat().st_size / (1024 * 1024)

                    # Check for cycle folders
                    rel_path = item.relative_to(stage_path)
                    for part in rel_path.parts:
                        if cycle_pattern.match(part):
                            cycles.add(part)

        # Scan Zarr for this stage (in each cycle)
        if zarr_path.exists():
            for cycle_dir in zarr_path.iterdir():
                if cycle_dir.is_dir() and not cycle_dir.name.startswith("."):
                    zarr_stage_path = cycle_dir / stage_name
                    if zarr_stage_path.exists():
                        zarr_stage_exists = True
                        cycles.add(cycle_dir.name)
                        # Count zarr chunk files and size
                        for item in zarr_stage_path.rglob("*"):
                            if item.is_file():
                                files.append(item)
                                info.total_size_mb += item.stat().st_size / (1024 * 1024)

        info.exists = tiff_exists or zarr_stage_exists
        info.file_count = len(files)
        info.cycle_folders = sorted(cycles)
        info.sample_files = [f.name for f in files[:5]]

        report.stages[stage_name] = info
        report.total_files += info.file_count
        report.total_size_mb += info.total_size_mb

    return report


def prompt_for_processed_data(
    report: ProcessedDataReport,
    interactive: bool = True,
) -> dict[str, Literal["keep", "delete", "cancel"]]:
    """
    Prompt user for what to do with existing processed data.

    Parameters
    ----------
    report : ProcessedDataReport
        Report of existing processed data
    interactive : bool
        If True, prompt user. If False, default to 'keep'.

    Returns
    -------
    dict
        Dictionary mapping stage names to actions ('keep', 'delete', 'cancel')
    """
    if not report.has_data:
        return {}

    stages_with_data = report.stages_with_data

    print("\n" + "=" * 60)
    print("EXISTING PROCESSED DATA DETECTED")
    print("=" * 60)
    print(report.summary())
    print()

    if not interactive:
        print("Non-interactive mode: keeping all existing data.")
        return dict.fromkeys(stages_with_data, "keep")

    print("Options:")
    print("  1. Keep all - Continue without deleting any data")
    print("  2. Delete all - Remove all processed data and start fresh")
    print("  3. Select per stage - Choose what to do for each stage")
    print("  4. Cancel - Exit without changes")
    print()

    while True:
        try:
            choice = input("Enter choice [1-4] (default: 1): ").strip() or "1"
            if choice == "1":
                print("Keeping all existing data.")
                return dict.fromkeys(stages_with_data, "keep")
            elif choice == "2":
                confirm = (
                    input("Are you sure you want to delete ALL processed data? [y/N]: ")
                    .strip()
                    .lower()
                )
                if confirm == "y":
                    print("Will delete all processed data.")
                    return dict.fromkeys(stages_with_data, "delete")
                else:
                    print("Cancelled. Keeping all data.")
                    return dict.fromkeys(stages_with_data, "keep")
            elif choice == "3":
                decisions = {}
                print("\nFor each stage, enter: k=keep, d=delete")
                for stage in stages_with_data:
                    info = report.stages[stage]
                    while True:
                        action = (
                            input(
                                f"  {stage} ({info.file_count} files, {info.total_size_mb:.1f} MB) [k/d]: "
                            )
                            .strip()
                            .lower()
                            or "k"
                        )
                        if action == "k":
                            decisions[stage] = "keep"
                            break
                        elif action == "d":
                            decisions[stage] = "delete"
                            break
                        else:
                            print("    Invalid choice. Enter k or d.")
                return decisions
            elif choice == "4":
                return dict.fromkeys(stages_with_data, "cancel")
            else:
                print("Invalid choice. Enter 1, 2, 3, or 4.")
        except (EOFError, KeyboardInterrupt):
            print("\nCancelled.")
            return dict.fromkeys(stages_with_data, "cancel")


def apply_processed_data_decisions(
    project_paths: ProjectPaths,
    decisions: dict[str, Literal["keep", "delete", "cancel"]],
) -> bool:
    """
    Apply user decisions for processed data (delete stages as requested).

    Deletes both TIFF directories and corresponding Zarr groups for each stage.

    Parameters
    ----------
    project_paths : ProjectPaths
        Project paths object
    decisions : dict
        Dictionary mapping stage names to actions

    Returns
    -------
    bool
        True if successful, False if cancelled
    """
    if any(action == "cancel" for action in decisions.values()):
        return False

    stage_paths = {
        "stitched": project_paths.stitched,
        "deconvolved": project_paths.deconvolved,
        "edf": project_paths.edf,
        "registered": project_paths.registered,
        "signal_isolated": project_paths.signal_isolated,
        "segmented": project_paths.segmented,
        "analysis": project_paths.analysis,
    }

    # Find the zarr file if it exists
    zarr_path = project_paths.processed / f"{project_paths.root.name}.zarr"

    for stage_name, action in decisions.items():
        if action == "delete" and stage_name in stage_paths:
            # Delete TIFF directory
            stage_path = stage_paths[stage_name]
            if stage_path.exists():
                print(f"  Deleting {stage_name} (TIFF)...")
                shutil.rmtree(stage_path)
                stage_path.mkdir(parents=True, exist_ok=True)

            # Delete corresponding Zarr groups (in each cycle)
            if zarr_path.exists():
                deleted_zarr = False
                for cycle_dir in zarr_path.iterdir():
                    if cycle_dir.is_dir() and not cycle_dir.name.startswith("."):
                        zarr_stage_path = cycle_dir / stage_name
                        if zarr_stage_path.exists():
                            if not deleted_zarr:
                                print(f"  Deleting {stage_name} (Zarr)...")
                                deleted_zarr = True
                            shutil.rmtree(zarr_stage_path)

            # Delete associated checkpoint files to prevent stale progress tracking
            cache_dir = project_paths.root / "cache"
            if cache_dir.exists():
                # Map stage names to their checkpoint file names
                checkpoint_map = {
                    "stitched": "stitching_checkpoint.json",
                    "deconvolved": "deconvolution_checkpoint.json",
                    "edf": "edf_checkpoint.json",
                    "registered": "registration_checkpoint.json",
                }
                if stage_name in checkpoint_map:
                    checkpoint_file = cache_dir / checkpoint_map[stage_name]
                    if checkpoint_file.exists():
                        print(f"  Clearing {stage_name} checkpoint...")
                        checkpoint_file.unlink()

    return True


def scan_existing_data(
    directory: str | Path,
    max_depth: int = 3,
    sample_count: int = 5,
) -> ExistingDataReport:
    """
    Scan a directory for existing data files.

    Differentiates between raw data (in data/raw/) and processed data
    (in data/processed/).

    Parameters
    ----------
    directory : str or Path
        Directory to scan
    max_depth : int
        Maximum depth to scan for files
    sample_count : int
        Number of images to sample for metadata extraction

    Returns
    -------
    ExistingDataReport
        Report of existing data found
    """
    import re

    directory = Path(directory)
    report = ExistingDataReport()

    if not directory.exists():
        return report

    image_files = []
    other_files = []
    raw_images = []
    processed_images: dict[str, list[Path]] = {}  # stage -> files
    cycle_pattern = re.compile(r"^cyc\d+", re.IGNORECASE)  # Matches cyc001, cyc001_DAPI_..., etc.
    filename_bases = set()

    # Define expected directory structure
    raw_dir = directory / "data" / "raw"
    processed_dir = directory / "data" / "processed"

    def is_image_file(path: Path) -> bool:
        """Check if file is an image based on extension."""
        suffix = path.suffix.lower()
        if path.name.lower().endswith(".ome.tif") or path.name.lower().endswith(".ome.tiff"):
            suffix = ".ome.tif"
        return suffix in IMAGE_EXTENSIONS

    def scan_dir(path: Path, depth: int = 0, location: str = "other"):
        """Scan directory, tracking location (raw, processed stage, or other)."""
        if depth > max_depth:
            return
        try:
            for item in path.iterdir():
                if item.is_dir():
                    # Check for cycle folder pattern
                    if cycle_pattern.match(item.name):
                        if location == "raw":
                            report.raw_cycle_folders.append(item.name)
                        report.cycle_folders.append(item.name)

                    # Determine location for subdirectory
                    sub_location = location
                    if location == "processed_root":
                        # Inside data/processed/, subdirs are processing stages
                        if item.name in PROCESSING_STAGES:
                            sub_location = f"processed_{item.name}"
                        else:
                            sub_location = "processed_other"

                    scan_dir(item, depth + 1, sub_location)
                elif item.is_file():
                    if is_image_file(item):
                        image_files.append(item)
                        # Extract filename pattern (remove numbers)
                        base = re.sub(r"\d+", "#", item.stem)
                        filename_bases.add(base)

                        # Track by location
                        if location == "raw":
                            raw_images.append(item)
                        elif location.startswith("processed_"):
                            stage = location.replace("processed_", "")
                            if stage not in processed_images:
                                processed_images[stage] = []
                            processed_images[stage].append(item)
                    elif item.suffix.lower() not in {
                        ".pyc",
                        ".log",
                        ".json",
                        ".md",
                        ".txt",
                        ".gitignore",
                    }:
                        other_files.append(item)
        except PermissionError:
            pass

    # Scan raw data directory
    if raw_dir.exists():
        scan_dir(raw_dir, depth=0, location="raw")

    # Scan processed data directory
    if processed_dir.exists():
        scan_dir(processed_dir, depth=0, location="processed_root")

    # Scan root directory for loose files (not in data/ structure)
    # This handles legacy layouts or files that haven't been organized
    for item in directory.iterdir():
        if item.is_dir() and item.name not in {
            "data",
            "meta",
            "notebooks",
            "slurm",
            ".claude",
            ".vscode",
            ".git",
            "__pycache__",
            "configs",
        }:
            # Check for cycle folders at root level (legacy layout)
            if cycle_pattern.match(item.name):
                report.cycle_folders.append(item.name)
            scan_dir(item, depth=1, location="other")
        elif item.is_file() and is_image_file(item):
            image_files.append(item)

    # Populate report
    report.image_files = image_files
    report.image_count = len(image_files)
    report.other_files = other_files
    report.has_data = len(image_files) > 0 or len(other_files) > 0
    report.filename_patterns = sorted(filename_bases)[:10]

    # Raw data stats
    report.has_raw_data = len(raw_images) > 0
    report.raw_image_count = len(raw_images)
    try:
        report.raw_size_mb = sum(f.stat().st_size for f in raw_images) / (1024 * 1024)
    except (OSError, PermissionError):
        pass

    # Processed data stats
    report.has_processed_data = len(processed_images) > 0
    report.processed_stages = {stage: len(files) for stage, files in processed_images.items()}
    try:
        all_processed = [f for files in processed_images.values() for f in files]
        report.processed_size_mb = sum(f.stat().st_size for f in all_processed) / (1024 * 1024)
    except (OSError, PermissionError):
        pass

    # Calculate total size
    try:
        report.total_size_mb = sum(f.stat().st_size for f in image_files) / (1024 * 1024)
    except (OSError, PermissionError):
        pass

    # Sort cycle folders naturally
    report.cycle_folders = sorted(set(report.cycle_folders))
    report.raw_cycle_folders = sorted(set(report.raw_cycle_folders))

    # Extract metadata from sample images (prefer raw images)
    sample_source = raw_images if raw_images else image_files
    if sample_source and sample_count > 0:
        sample_files = sample_source[: min(sample_count, len(sample_source))]
        for img_path in sample_files:
            try:
                metadata = extract_image_metadata(img_path)
                if metadata:
                    report.metadata_samples.append(
                        {
                            "file": str(img_path.name),
                            "dimensions": metadata.dimensions,
                            "channels": metadata.channels,
                            "channel_names": metadata.channel_names,
                            "pixel_size": metadata.pixel_size,
                            "pixel_unit": metadata.pixel_unit,
                            "is_ome": metadata.is_ome,
                        }
                    )
            except Exception:
                pass

    return report


def extract_image_metadata(path: str | Path) -> ImageMetadata | None:
    """
    Extract metadata from an image file.

    Parameters
    ----------
    path : str or Path
        Path to image file

    Returns
    -------
    ImageMetadata or None
        Extracted metadata, or None if extraction failed
    """
    path = Path(path)

    if not path.exists():
        return None

    metadata = ImageMetadata(
        path=path,
        filename=path.name,
        size_mb=path.stat().st_size / (1024 * 1024),
    )

    suffix = path.suffix.lower()

    # Try tifffile for TIFF images
    if suffix in {".tif", ".tiff"} or path.name.lower().endswith((".ome.tif", ".ome.tiff")):
        try:
            import tifffile

            with tifffile.TiffFile(str(path)) as tif:
                # Basic dimensions
                series = tif.series[0] if tif.series else None
                if series:
                    metadata.dimensions = series.shape
                    metadata.dtype = str(series.dtype)
                elif tif.pages:
                    page = tif.pages[0]
                    metadata.dimensions = (len(tif.pages), page.shape[0], page.shape[1])
                    metadata.dtype = str(page.dtype)

                # Check for OME metadata
                if tif.ome_metadata:
                    metadata.is_ome = True
                    _parse_ome_metadata(tif.ome_metadata, metadata)

                # Check TIFF tags
                if tif.pages:
                    page = tif.pages[0]
                    tags = page.tags

                    # Software
                    if "Software" in tags:
                        metadata.software = str(tags["Software"].value)[:100]

                    # DateTime
                    if "DateTime" in tags:
                        metadata.datetime = str(tags["DateTime"].value)

                    # ImageDescription (may contain metadata)
                    if "ImageDescription" in tags:
                        desc = str(tags["ImageDescription"].value)
                        metadata.description = desc[:500] if len(desc) > 500 else desc

                    # Resolution
                    if "XResolution" in tags and "ResolutionUnit" in tags:
                        try:
                            xres = tags["XResolution"].value
                            if isinstance(xres, tuple):
                                xres = xres[0] / xres[1] if xres[1] else xres[0]
                            unit = tags["ResolutionUnit"].value
                            if unit == 3:  # centimeter
                                metadata.pixel_size = 10000 / xres  # convert to microns
                                metadata.pixel_unit = "um"
                            elif unit == 2:  # inch
                                metadata.pixel_size = 25400 / xres
                                metadata.pixel_unit = "um"
                        except (TypeError, ZeroDivisionError):
                            pass

        except ImportError:
            pass
        except Exception:
            pass

    # Try to determine channels from shape
    if metadata.dimensions:
        shape = metadata.dimensions
        if len(shape) >= 3:
            # Assume smallest dim < 10 is channels
            for _i, dim in enumerate(shape):
                if dim < 10 and dim > 1:
                    metadata.channels = dim
                    break

    return metadata


def _parse_ome_metadata(ome_xml: str, metadata: ImageMetadata) -> None:
    """Parse OME-XML metadata and populate ImageMetadata."""
    try:
        import xml.etree.ElementTree as ET

        root = ET.fromstring(ome_xml)
        ns = {"ome": "http://www.openmicroscopy.org/Schemas/OME/2016-06"}

        # Try to find Image element
        image = root.find(".//ome:Image", ns)
        if image is None:
            # Try without namespace
            image = root.find(".//Image")

        if image is not None:
            # Get acquisition date
            acq_date = image.get("AcquisitionDate")
            if acq_date:
                metadata.datetime = acq_date

        # Find Pixels element
        pixels = root.find(".//ome:Pixels", ns)
        if pixels is None:
            pixels = root.find(".//Pixels")

        if pixels is not None:
            # Get physical sizes
            px_size = pixels.get("PhysicalSizeX")
            if px_size:
                try:
                    metadata.pixel_size = float(px_size)
                    metadata.pixel_unit = pixels.get("PhysicalSizeXUnit", "um")
                except ValueError:
                    pass

            # Get channel count and names
            channels = root.findall(".//ome:Channel", ns)
            if not channels:
                channels = root.findall(".//Channel")

            if channels:
                metadata.channels = len(channels)
                metadata.channel_names = []
                for ch in channels:
                    name = ch.get("Name") or ch.get("ID", "")
                    if name:
                        metadata.channel_names.append(name)

    except Exception:
        pass


@dataclass
class ProjectPaths:
    """Structured paths for a KINTSUGI project."""

    # Root project directory
    root: Path

    # Input data
    raw: Path  # Raw images from microscope

    # Processing outputs
    processed: Path  # Base processed directory
    stitched: Path  # Stitched images (illumination corrected)
    deconvolved: Path  # Deconvolved images
    edf: Path  # Extended depth of focus projections
    registered: Path  # Registered multi-cycle images
    signal_isolated: Path  # Signal-isolated images (autofluorescence removed)
    segmented: Path  # Segmentation masks
    analysis: Path  # Analysis outputs (tables, plots)

    # Working files
    notebooks: Path  # Working copies of notebooks
    configs: Path  # Project-specific configurations
    logs: Path  # Processing logs
    cache: Path  # Temporary/cache files

    # Metadata
    meta: Path  # Metadata files

    @classmethod
    def from_root(cls, root: str | Path) -> ProjectPaths:
        """Create ProjectPaths from a root directory."""
        root = Path(root).resolve()

        return cls(
            root=root,
            raw=root / "data" / "raw",
            processed=root / "data" / "processed",
            stitched=root / "data" / "processed" / "stitched",
            deconvolved=root / "data" / "processed" / "deconvolved",
            edf=root / "data" / "processed" / "edf",
            registered=root / "data" / "processed" / "registered",
            signal_isolated=root / "data" / "processed" / "signal_isolated",
            segmented=root / "data" / "processed" / "segmented",
            analysis=root / "data" / "processed" / "analysis",
            notebooks=root / "notebooks",
            configs=root / "configs",
            logs=root / "logs",
            cache=root / ".cache",
            meta=root / "meta",
        )

    def create_all(self) -> None:
        """Create all directories."""
        for name, path in self.__dict__.items():
            if isinstance(path, Path) and name != "root":
                path.mkdir(parents=True, exist_ok=True)

    def to_dict(self) -> dict[str, str]:
        """Convert to dictionary of string paths."""
        return {k: str(v) for k, v in self.__dict__.items()}


@dataclass
class ProjectConfig:
    """Project configuration and metadata."""

    # Project identification
    name: str
    description: str = ""
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    modified: str = field(default_factory=lambda: datetime.now().isoformat())
    version: str = PROJECT_VERSION

    # KINTSUGI info
    kintsugi_version: str = ""
    kintsugi_path: str = ""

    # Environment info
    python_version: str = field(default_factory=lambda: sys.version)
    platform: str = field(default_factory=platform.platform)

    # Processing parameters (defaults)
    parameters: dict[str, Any] = field(default_factory=dict)

    # Cycle information
    cycles: list[dict[str, Any]] = field(default_factory=list)

    # Channel metadata
    channels: list[dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ProjectConfig:
        """Create from dictionary."""
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})

    def update_modified(self) -> None:
        """Update the modified timestamp."""
        self.modified = datetime.now().isoformat()


class KintsugiProject:
    """
    KINTSUGI Project Manager

    Manages the separation of:
    - Repository code (read-only, version controlled)
    - Project data (user's experimental data)
    - Working notebooks (copies that can be modified)

    Example
    -------
    >>> project = KintsugiProject.create(
    ...     "C:/Projects/Experiment_2024",
    ...     name="My Experiment",
    ...     description="CODEX data from pancreas samples"
    ... )
    >>>
    >>> # Get paths
    >>> print(project.paths.raw)  # Where to put raw data
    >>> print(project.paths.stitched)  # Where stitched images go
    >>>
    >>> # Copy notebooks to work on
    >>> project.setup_notebooks()
    >>>
    >>> # Save/load project
    >>> project.save()
    >>> project = KintsugiProject.load("C:/Projects/Experiment_2024")
    """

    def __init__(
        self,
        root: str | Path,
        config: ProjectConfig | None = None,
        kintsugi_path: str | Path | None = None,
    ):
        """
        Initialize a project.

        Parameters
        ----------
        root : str or Path
            Root directory for the project
        config : ProjectConfig, optional
            Project configuration
        kintsugi_path : str or Path, optional
            Path to KINTSUGI repository. Auto-detected if not provided.
        """
        self.root = Path(root).resolve()
        self.paths = ProjectPaths.from_root(self.root)
        self.config = config or ProjectConfig(name=self.root.name)

        # Detect or use provided KINTSUGI path
        if kintsugi_path:
            self._kintsugi_path = Path(kintsugi_path).resolve()
        else:
            self._kintsugi_path = self._detect_kintsugi_path()

        self.config.kintsugi_path = str(self._kintsugi_path)

        # Get KINTSUGI version
        try:
            import kintsugi

            self.config.kintsugi_version = getattr(kintsugi, "__version__", "unknown")
        except ImportError:
            self.config.kintsugi_version = "not installed"

    @staticmethod
    def _detect_kintsugi_path() -> Path:
        """Auto-detect KINTSUGI repository path."""
        # Try to find from package location
        try:
            import kintsugi

            pkg_path = Path(kintsugi.__file__).parent
            # Go up from src/kintsugi to repo root
            repo_path = pkg_path.parent.parent
            if (repo_path / "pyproject.toml").exists():
                return repo_path
        except ImportError:
            pass

        # Try common locations
        common_paths = [
            Path.home() / "KINTSUGI",
            Path("C:/Users") / os.getlogin() / "KINTSUGI",
            Path.cwd() / "KINTSUGI",
        ]

        for path in common_paths:
            if path.exists() and (path / "pyproject.toml").exists():
                return path

        # Return current directory as fallback
        return Path.cwd()

    @classmethod
    def scan_directory(cls, root: str | Path) -> ExistingDataReport:
        """
        Scan a directory for existing data before project creation.

        Use this to check for existing files and extract metadata
        before calling create().

        Parameters
        ----------
        root : str or Path
            Directory to scan

        Returns
        -------
        ExistingDataReport
            Report of existing data found
        """
        return scan_existing_data(root)

    @classmethod
    def create(
        cls,
        root: str | Path,
        name: str | None = None,
        description: str = "",
        kintsugi_path: str | Path | None = None,
        setup_notebooks: bool = True,
        existing_data_report: ExistingDataReport | None = None,
        slurm: bool = False,
        slurm_account: str | None = None,
        slurm_partition: str | None = None,
        slurm_qos: str | None = None,
        slurm_gpu_type: str | None = None,
        # Microscope parameters (stored in /meta/experiment.json)
        tile_rows: int | None = None,
        tile_cols: int | None = None,
        xy_pixel_size: float = 377.0,
        z_step_size: float = 1500.0,
        numerical_aperture: float = 0.75,
        tissue_refractive_index: float = 1.44,
    ) -> KintsugiProject:
        """
        Create a new KINTSUGI project.

        IMPORTANT: This method will NEVER delete existing files. It only creates
        new directories and files. Existing data is preserved.

        Parameters
        ----------
        root : str or Path
            Root directory for the project
        name : str, optional
            Project name (defaults to directory name)
        description : str
            Project description
        kintsugi_path : str or Path, optional
            Path to KINTSUGI repository
        setup_notebooks : bool
            Whether to copy notebook templates
        existing_data_report : ExistingDataReport, optional
            Pre-scanned data report (avoids re-scanning)
        slurm : bool
            If True, set up SLURM job submission infrastructure
        slurm_account : str, optional
            HPC account name for SLURM. Auto-detected if not provided.
        slurm_partition : str, optional
            SLURM partition. Defaults to 'gpu'.
        slurm_qos : str, optional
            Quality of Service for SLURM.
        slurm_gpu_type : str, optional
            GPU type for SLURM constraints. Auto-detected from nvidia-smi.
        tile_rows : int
            Number of tile rows in the acquisition grid (default: 5)
        tile_cols : int
            Number of tile columns in the acquisition grid (default: 5)
        xy_pixel_size : float
            XY pixel size in nanometers (default: 377)
        z_step_size : float
            Z step size in nanometers (default: 1500)
        numerical_aperture : float
            Objective numerical aperture (default: 0.75)
        tissue_refractive_index : float
            Tissue refractive index (default: 1.44)

        Returns
        -------
        KintsugiProject
            Initialized project
        """
        root = Path(root).resolve()

        # Check if project already exists
        config_file = root / PROJECT_CONFIG_FILE
        if config_file.exists():
            warnings.warn(
                f"Project already exists at {root}. Use load() to open it.",
                UserWarning,
                stacklevel=2,
            )
            return cls.load(root)

        # Scan for existing data if not provided
        if existing_data_report is None:
            existing_data_report = scan_existing_data(root)

        # Create project
        config = ProjectConfig(
            name=name or root.name,
            description=description,
        )

        project = cls(root, config, kintsugi_path)

        # Store discovered metadata in project config
        if existing_data_report.has_data:
            config.parameters["discovered_data"] = {
                "image_count": existing_data_report.image_count,
                "total_size_mb": existing_data_report.total_size_mb,
                "cycle_folders": existing_data_report.cycle_folders,
                "filename_patterns": existing_data_report.filename_patterns,
                "metadata_samples": existing_data_report.metadata_samples,
            }

        # Store microscope parameters for experiment config creation
        config.parameters["microscope"] = {
            "tile_rows": tile_rows,
            "tile_cols": tile_cols,
            "xy_pixel_size": xy_pixel_size,
            "z_step_size": z_step_size,
            "numerical_aperture": numerical_aperture,
            "tissue_refractive_index": tissue_refractive_index,
        }

        # Create directory structure (uses exist_ok=True, never overwrites)
        project.paths.create_all()

        # Copy notebooks if requested
        if setup_notebooks:
            project.setup_notebooks()

        # Create default config files (including experiment.json)
        project._create_default_configs()

        # Set up SLURM if requested
        if slurm:
            project.setup_slurm(
                account=slurm_account,
                partition=slurm_partition,
                qos=slurm_qos,
                gpu_type=slurm_gpu_type,
            )

        # Save project
        project.save()

        print(f"[OK] Created KINTSUGI project: {project.config.name}")
        print(f"  Location: {root}")
        print(f"  KINTSUGI: {project._kintsugi_path}")
        print()
        print("Directory structure:")
        print(f"  {root}/")
        print("  +-- data/")
        print("  |   +-- raw/          <- Put raw images here")
        print("  |   +-- processed/    <- Outputs go here")
        print("  +-- notebooks/        <- Working notebooks")
        print("  +-- configs/          <- Processing configs")
        print("  +-- meta/             <- Metadata files")
        print("  +-- logs/             <- Processing logs")
        if slurm:
            print("  +-- slurm/            <- SLURM job submission")
        print("  +-- .claude/          <- Claude Code config")
        print("  +-- .vscode/          <- VS Code config")
        print()

        if existing_data_report.has_data:
            print("Existing data detected:")
            if existing_data_report.has_raw_data:
                print(f"  Raw: {existing_data_report.raw_image_count} images")
                if existing_data_report.raw_cycle_folders:
                    print(f"  Cycles: {', '.join(existing_data_report.raw_cycle_folders)}")
            if existing_data_report.has_processed_data:
                print(f"  Processed: {existing_data_report.processed_size_mb:.1f} MB")
                for stage, count in existing_data_report.processed_stages.items():
                    print(f"    {stage}: {count} files")
            print()

        print("Next steps:")
        if existing_data_report.has_raw_data:
            print("  1. Raw data found - ready to process")
        else:
            print("  1. Copy raw images to data/raw/")
        print("  2. Open folder in VS Code: code .")
        print("  3. Start with notebooks/1_Single_Channel_Eval.ipynb")

        return project

    def _adopt_existing_data(self, report: ExistingDataReport) -> None:
        """
        Organize existing data into the project structure.

        Moves cycle folders and images to data/raw/ while preserving structure.
        NEVER deletes data - only moves it within the project.
        """
        import shutil

        # If cycle folders exist at root level, move them to raw
        for cycle_name in report.cycle_folders:
            src = self.root / cycle_name
            if src.exists() and src.is_dir():
                dest = self.paths.raw / cycle_name
                if not dest.exists():
                    print(f"  Moving {cycle_name}/ -> data/raw/{cycle_name}/")
                    shutil.move(str(src), str(dest))

        # Check for images directly in root (not in cycle folders)
        for img_path in report.image_files:
            # Only move images that are directly in root (not in subdirs)
            if img_path.parent == self.root:
                dest = self.paths.raw / img_path.name
                if not dest.exists():
                    print(f"  Moving {img_path.name} -> data/raw/{img_path.name}")
                    shutil.move(str(img_path), str(dest))

    @classmethod
    def load(cls, root: str | Path) -> KintsugiProject:
        """
        Load an existing KINTSUGI project.

        Parameters
        ----------
        root : str or Path
            Root directory of the project

        Returns
        -------
        KintsugiProject
            Loaded project
        """
        root = Path(root).resolve()
        config_file = root / PROJECT_CONFIG_FILE

        if not config_file.exists():
            raise FileNotFoundError(
                f"No KINTSUGI project found at {root}. Use create() to create a new project."
            )

        with open(config_file) as f:
            data = json.load(f)

        config = ProjectConfig.from_dict(data.get("config", {}))
        kintsugi_path = data.get("kintsugi_path")

        project = cls(root, config, kintsugi_path)

        print(f"[OK] Loaded project: {project.config.name}")
        print(f"  Created: {project.config.created[:10]}")
        print(f"  Modified: {project.config.modified[:10]}")

        return project

    def save(self) -> None:
        """Save project configuration."""
        self.config.update_modified()

        data = {
            "version": PROJECT_VERSION,
            "kintsugi_path": str(self._kintsugi_path),
            "paths": self.paths.to_dict(),
            "config": self.config.to_dict(),
        }

        config_file = self.root / PROJECT_CONFIG_FILE
        with open(config_file, "w") as f:
            json.dump(data, f, indent=2)

    # =========================================================================
    # EXPERIMENT CONFIGURATION
    # =========================================================================

    def get_experiment_config_path(self) -> Path:
        """Get path to experiment.json file."""
        return self.paths.meta / EXPERIMENT_CONFIG_FILE

    def load_experiment_config(self) -> ExperimentConfig | None:
        """
        Load experiment configuration from /meta/experiment.json.

        Returns
        -------
        ExperimentConfig or None
            Experiment configuration, or None if file doesn't exist
        """
        config_path = self.get_experiment_config_path()
        if not config_path.exists():
            return None

        try:
            with open(config_path) as f:
                data = json.load(f)
            return ExperimentConfig.from_dict(data)
        except (json.JSONDecodeError, TypeError, KeyError, AttributeError, ValueError) as e:
            print(f"Warning: Could not load experiment config: {e}")
            return None

    def save_experiment_config(self, config: ExperimentConfig) -> Path:
        """
        Save experiment configuration to /meta/experiment.json.

        Parameters
        ----------
        config : ExperimentConfig
            Configuration to save

        Returns
        -------
        Path
            Path to saved file
        """
        self.paths.meta.mkdir(parents=True, exist_ok=True)
        config_path = self.get_experiment_config_path()

        with open(config_path, "w") as f:
            json.dump(config.to_dict(), f, indent=2)

        return config_path

    def create_experiment_config(
        self,
        tile_rows: int | None = None,
        tile_cols: int | None = None,
        xy_pixel_size: float = 377.0,
        z_step_size: float = 1500.0,
        numerical_aperture: float = 0.75,
        tissue_refractive_index: float = 1.44,
        wavelengths: dict[int, tuple[float, float]] | None = None,
        auto_detect: bool = True,
        overwrite: bool = False,
        **kwargs,
    ) -> ExperimentConfig:
        """
        Create experiment configuration with optional auto-detection.

        Parameters
        ----------
        tile_rows, tile_cols : int or None
            Tile grid dimensions. If None, auto-detected from raw data.
        xy_pixel_size : float
            XY pixel size in nanometers
        z_step_size : float
            Z step size in nanometers
        numerical_aperture : float
            Objective numerical aperture
        tissue_refractive_index : float
            Tissue refractive index
        wavelengths : dict, optional
            Channel wavelengths {channel: (excitation, emission)}
        auto_detect : bool
            If True, auto-detect n_cycles, n_zplanes, tile grid from raw data
        overwrite : bool
            If True, overwrite existing config
        **kwargs
            Additional ExperimentConfig fields

        Returns
        -------
        ExperimentConfig
            Created configuration
        """
        config_path = self.get_experiment_config_path()
        if config_path.exists() and not overwrite:
            existing = self.load_experiment_config()
            if existing:
                print(f"  Experiment config already exists: {config_path}")
                return existing

        # Track whether user explicitly provided tile dimensions
        user_set_tiles = tile_rows is not None or tile_cols is not None

        # Create config with provided parameters (use defaults for None)
        # Use default wavelengths if not provided
        default_wavelengths = {
            1: (358.0, 461.0),  # DAPI
            2: (753.0, 775.0),  # Cy7
            3: (560.0, 575.0),  # TRITC
            4: (648.0, 668.0),  # Cy5
        }
        config = ExperimentConfig(
            tile_rows=tile_rows if tile_rows is not None else 5,
            tile_cols=tile_cols if tile_cols is not None else 5,
            xy_pixel_size=xy_pixel_size,
            z_step_size=z_step_size,
            numerical_aperture=numerical_aperture,
            tissue_refractive_index=tissue_refractive_index,
            wavelengths=wavelengths if wavelengths is not None else default_wavelengths,
            **kwargs,
        )

        # Auto-detect from raw data if requested
        if auto_detect:
            import re

            cycles = find_raw_cycles(self.paths.raw)
            if cycles:
                config.n_cycles = len(cycles)
                print(f"  Auto-detected {config.n_cycles} cycles")

                # Try to detect z-planes from first cycle
                first_cycle = cycles[0]
                z_planes = set()
                for f in first_cycle.glob("*_Z*_CH1.tif"):
                    match = re.search(r"_Z(\d+)_", f.name)
                    if match:
                        z_planes.add(int(match.group(1)))
                if z_planes:
                    config.n_zplanes = len(z_planes)
                    print(f"  Auto-detected {config.n_zplanes} z-planes")

                # Auto-detect tile grid and channels_per_cycle
                # Count tiles at Z001 CH1 to get tile count
                z001_ch1_files = list(first_cycle.glob("*_Z001_CH1.tif"))
                if z001_ch1_files:
                    n_tiles = len(z001_ch1_files)
                    if not user_set_tiles:
                        rows, cols = _tile_count_to_grid(n_tiles)
                        config.tile_rows = rows
                        config.tile_cols = cols
                        print(f"  Auto-detected tile grid: {rows}x{cols} ({n_tiles} tiles)")

                    # Detect channels_per_cycle from distinct CH values
                    channels = set()
                    for f in first_cycle.glob("*_Z001_CH*.tif"):
                        ch_match = re.search(r"_CH(\d+)", f.name)
                        if ch_match:
                            channels.add(int(ch_match.group(1)))
                    if channels:
                        config.channels_per_cycle = len(channels)
                        print(f"  Auto-detected {config.channels_per_cycle} channels per cycle")

        # Save config
        self.save_experiment_config(config)
        print(f"  Created experiment config: {config_path}")

        return config

    def load_channel_names(self) -> dict[int, list[str]] | None:
        """
        Load channel names from /meta/CHANNELNAMES.txt.

        Returns
        -------
        dict or None
            Dictionary mapping cycle number to list of channel names,
            or None if file not found
        """
        # Import from Kio if available, otherwise use inline implementation
        try:
            from notebooks.Kio import load_channel_names as kio_load_channel_names

            return kio_load_channel_names(self.paths.meta)
        except ImportError:
            pass

        # Inline implementation (simplified)
        import re

        meta_dir = self.paths.meta
        channel_file = None

        # Try to find channel names file
        for name in ["CHANNELNAMES.txt", "channelnames.txt", "channel_names.txt", "channels.txt"]:
            candidate = meta_dir / name
            if candidate.exists():
                channel_file = candidate
                break

        if channel_file is None:
            return None

        # Read and parse file
        lines = []
        with open(channel_file) as f:
            for line in f:
                line = line.strip()
                if line and not line.startswith("#"):
                    lines.append(line)

        if not lines:
            return None

        channel_dict = {}
        first_line = lines[0]

        # Check if cycle-prefixed format
        if (
            ":" in first_line
            or "\t" in first_line
            or (first_line.split(",")[0].strip().isdigit() and len(first_line.split(",")) > 2)
        ):
            # Cycle-prefixed format
            for line in lines:
                try:
                    if ":" in line:
                        cycle_str, names_str = line.split(":", 1)
                        cycle = int(cycle_str.strip())
                        names = [n.strip() for n in names_str.split(",")]
                    elif "\t" in line:
                        parts = line.split("\t")
                        cycle = int(parts[0].strip())
                        names = [n.strip() for n in parts[1:]]
                    else:
                        parts = line.split(",")
                        cycle = int(parts[0].strip())
                        names = [n.strip() for n in parts[1:]]
                    channel_dict[cycle] = names
                except (ValueError, IndexError):
                    continue
        else:
            # Simple list format - detect cycles from DAPI-XX pattern
            current_cycle = 0
            cycle_channels = []
            channels_per_cycle = 4

            for line in lines:
                dapi_match = re.match(r"DAPI[-_]?(\d+)", line, re.IGNORECASE)
                dapi_bare = re.match(r"DAPI\s*$", line, re.IGNORECASE)
                if dapi_match:
                    if cycle_channels and current_cycle > 0:
                        channel_dict[current_cycle] = cycle_channels
                    current_cycle = int(dapi_match.group(1))
                    cycle_channels = [line]
                elif dapi_bare:
                    # Bare "DAPI" without cycle suffix — auto-increment cycle
                    if cycle_channels and current_cycle > 0:
                        channel_dict[current_cycle] = cycle_channels
                    current_cycle += 1
                    cycle_channels = [line]
                elif current_cycle > 0:
                    cycle_channels.append(line)
                    if len(cycle_channels) == channels_per_cycle:
                        channel_dict[current_cycle] = cycle_channels
                        cycle_channels = []

            if cycle_channels and current_cycle > 0:
                channel_dict[current_cycle] = cycle_channels

        return channel_dict if channel_dict else None

    def setup_notebooks(self, overwrite: bool = False) -> list[Path]:
        """
        Copy notebook templates to project directory.

        Parameters
        ----------
        overwrite : bool
            Whether to overwrite existing notebooks

        Returns
        -------
        List[Path]
            Paths to copied notebooks
        """
        source_dir = self._kintsugi_path / "notebooks"
        dest_dir = self.paths.notebooks

        if not source_dir.exists():
            raise FileNotFoundError(f"Notebook templates not found at {source_dir}")
        dest_dir.mkdir(parents=True, exist_ok=True)

        # Notebooks to copy (main workflow)
        workflow_notebooks = [
            "1_Single_Channel_Eval.ipynb",
            "2_Cycle_Processing.ipynb",
            "3_Signal_Isolation_QC.ipynb",
            "4_Segmentation_Analysis.ipynb",
        ]

        copied = []
        for nb_name in workflow_notebooks:
            src = source_dir / nb_name
            dst = dest_dir / nb_name

            if not src.exists():
                continue

            existed = dst.exists()
            if dst.exists() and not overwrite:
                print(f"  Skipping {nb_name} (exists)")
                continue

            shutil.copy2(src, dst)
            copied.append(dst)
            action = "Updated" if existed else "Copied"
            print(f"  {action} {nb_name}")

        # Copy supporting modules
        support_dirs = ["Kreg", "Kstitch", "Kview", "Kdecon", "Kseg"]
        for dir_name in support_dirs:
            src_dir = source_dir / dir_name
            dst_dir = dest_dir / dir_name

            if not src_dir.exists():
                continue

            existed = dst_dir.exists()
            if dst_dir.exists() and not overwrite:
                continue

            # dirs_exist_ok updates existing files without removing user data
            shutil.copytree(src_dir, dst_dir, dirs_exist_ok=True)
            action = "Updated" if existed else "Copied"
            print(f"  {action} {dir_name}/")

        # Copy supporting files
        support_files = ["Kutils.py", "config_example.json", "MIGRATION_GUIDE.md"]
        for file_name in support_files:
            src_file = source_dir / file_name
            dst_file = dest_dir / file_name

            if not src_file.exists():
                continue

            existed = dst_file.exists()
            if dst_file.exists() and not overwrite:
                continue

            shutil.copy2(src_file, dst_file)
            action = "Updated" if existed else "Copied"
            print(f"  {action} {file_name}")

        return copied

    def _create_default_params_file(self) -> Path:
        """Create or regenerate default_parameters.json from experiment.json.

        Returns
        -------
        Path
            Path to the written default_parameters.json file.
        """
        # Read overlap and reference cycle from experiment.json if available
        overlap_pct = 30  # CODEX standard default
        ref_cycle = 1
        exp_path = self.paths.meta / "experiment.json"
        if exp_path.exists():
            try:
                with open(exp_path) as f:
                    exp_data = json.load(f)
                overlap_raw = exp_data.get("tileOverlapX", exp_data.get("tile_overlap", None))
                if overlap_raw is not None:
                    val = float(overlap_raw)
                    # Handle both fraction (0.3) and percent (30) formats
                    overlap_pct = int(val * 100) if val < 1 else int(val)
                ref_raw = exp_data.get("referenceCycle", exp_data.get("reference_cycle"))
                if ref_raw is not None:
                    ref_cycle = int(ref_raw)
            except Exception:
                pass  # Use defaults on read error

        default_params = {
            "illumination_correction": {
                "method": "BaSiC",
                "working_size": 128,
                "max_iterations": 500,
            },
            "stitching": {
                "overlap_percent": overlap_pct,
                "blend_mode": "linear",
            },
            "registration": {
                "method": "rigid",
                "reference_cycle": ref_cycle,
            },
            "segmentation": {
                "method": "instanseg",
                "cell_diameter": 30,
            },
        }

        self.paths.configs.mkdir(parents=True, exist_ok=True)
        params_file = self.paths.configs / "default_parameters.json"
        with open(params_file, "w") as f:
            json.dump(default_params, f, indent=2)

        return params_file

    def _create_default_configs(self) -> None:
        """Create default configuration files."""
        self._create_default_params_file()

        # Create experiment configuration with auto-detection
        # Use microscope parameters from project config if available
        microscope_params = self.config.parameters.get("microscope", {})
        self.create_experiment_config(
            tile_rows=microscope_params.get("tile_rows"),
            tile_cols=microscope_params.get("tile_cols"),
            xy_pixel_size=microscope_params.get("xy_pixel_size", 377.0),
            z_step_size=microscope_params.get("z_step_size", 1500.0),
            numerical_aperture=microscope_params.get("numerical_aperture", 0.75),
            tissue_refractive_index=microscope_params.get("tissue_refractive_index", 1.44),
            auto_detect=True,
            overwrite=False,
        )

        # Create Claude Code configuration
        self._create_claude_config()

        # Create VS Code configuration
        self._create_vscode_config()

    def _create_claude_config(self) -> None:
        """Create Claude Code MCP configuration.

        Creates .mcp.json at project root (MCP server definition) and
        .claude/settings.local.json (enables the server).
        """
        # Create .mcp.json with the MCP server definition
        mcp_file = self.paths.root / ".mcp.json"
        if not mcp_file.exists():
            mcp_config = {
                "mcpServers": {
                    "kintsugi": {
                        "command": _resolve_kintsugi_executable(),
                        "args": ["mcp", "start"],
                        "cwd": str(self.paths.root),
                    }
                }
            }
            with open(mcp_file, "w") as f:
                json.dump(mcp_config, f, indent=2)
                f.write("\n")
            print("  Created .mcp.json")

        # Create .claude/settings.local.json to enable the MCP server
        claude_dir = self.paths.root / ".claude"
        claude_dir.mkdir(exist_ok=True)
        settings_file = claude_dir / "settings.local.json"
        if not settings_file.exists():
            settings_config = {
                "enableAllProjectMcpServers": True,
                "enabledMcpjsonServers": ["kintsugi"],
            }
            with open(settings_file, "w") as f:
                json.dump(settings_config, f, indent=2)
                f.write("\n")
            print("  Created .claude/settings.local.json")

    def _create_vscode_config(self) -> None:
        """Create VS Code configuration."""
        vscode_dir = self.paths.root / ".vscode"
        vscode_dir.mkdir(exist_ok=True)

        settings_file = vscode_dir / "settings.json"
        if not settings_file.exists():
            vscode_config = {
                "python.defaultInterpreterPath": "${env:CONDA_PREFIX}/python",
                "jupyter.notebookFileRoot": "${workspaceFolder}/notebooks",
                "files.exclude": {
                    "**/__pycache__": True,
                    "**/.ipynb_checkpoints": True,
                    "**/*.pyc": True,
                },
                "python.analysis.extraPaths": [
                    "${workspaceFolder}/notebooks",
                ],
            }
            with open(settings_file, "w") as f:
                json.dump(vscode_config, f, indent=2)
            print("  Created .vscode/settings.json")

    def setup_slurm(
        self,
        account: str | None = None,
        partition: str | None = None,
        qos: str | None = None,
        gpu_type: str | None = None,
        auto_detect: bool = True,
    ) -> Path:
        """
        Set up SLURM job submission infrastructure for the project.

        Creates the slurm/ directory with:
        - config.sh: Pre-populated configuration file
        - jobs/: Symlink to KINTSUGI job scripts
        - README.md: Quick-start guide

        Parameters
        ----------
        account : str, optional
            HPC account name. Auto-detected if not provided.
        partition : str, optional
            SLURM partition. Auto-detected or defaults to 'gpu'.
        qos : str, optional
            Quality of Service. Auto-detected if not provided.
        gpu_type : str, optional
            GPU type for constraints. Auto-detected from nvidia-smi.
        auto_detect : bool
            Whether to auto-detect settings from environment (default: True)

        Returns
        -------
        Path
            Path to the created slurm directory
        """
        from kintsugi.hpc import (
            detect_hpc_settings,
            generate_slurm_config,
            generate_slurm_readme,
        )

        slurm_dir = self.root / "slurm"
        slurm_dir.mkdir(exist_ok=True)

        # Auto-detect HPC settings if requested
        detected: dict[str, Any] = {}
        if auto_detect:
            detected = detect_hpc_settings()

        # Build user settings from provided arguments
        user_settings: dict[str, Any] = {}
        if account is not None:
            user_settings["account"] = account
        if partition is not None:
            user_settings["partition"] = partition
        if qos is not None:
            user_settings["qos"] = qos
        if gpu_type is not None:
            user_settings["gpu_type"] = gpu_type

        # Load experiment config to populate processing parameters
        exp_config = self.load_experiment_config()
        experiment_dict = exp_config.to_dict() if exp_config else None

        # Generate config.sh
        config_content = generate_slurm_config(
            project_name=self.config.name,
            project_dir=str(self.root),
            kintsugi_dir=str(self._kintsugi_path),
            detected_settings=detected,
            user_settings=user_settings,
            experiment_config=experiment_dict,
        )

        config_file = slurm_dir / "config.sh"
        with open(config_file, "w") as f:
            f.write(config_content)
        print("  Created slurm/config.sh")

        # Create symlink to job scripts
        jobs_link = slurm_dir / "jobs"
        jobs_source = self._kintsugi_path / "slurm" / "jobs"

        if jobs_source.exists():
            if jobs_link.exists() or jobs_link.is_symlink():
                jobs_link.unlink()
            try:
                jobs_link.symlink_to(jobs_source)
                print(f"  Created slurm/jobs -> {jobs_source}")
            except OSError as e:
                # Symlink may fail on Windows without admin rights
                print(f"  Warning: Could not create symlink for jobs: {e}")
                print(f"  Job scripts are at: {jobs_source}")

        # Generate README
        readme_content = generate_slurm_readme(
            project_name=self.config.name,
            kintsugi_dir=str(self._kintsugi_path),
        )

        readme_file = slurm_dir / "README.md"
        with open(readme_file, "w") as f:
            f.write(readme_content)
        print("  Created slurm/README.md")

        # Log detected settings
        if detected:
            detected_items = [f"{k}={v}" for k, v in detected.items() if v]
            if detected_items:
                print(f"  Auto-detected: {', '.join(detected_items)}")

        return slurm_dir

    # -------------------------------------------------------------------------
    # Convenience properties for common paths
    # -------------------------------------------------------------------------

    @property
    def raw_dir(self) -> Path:
        """Raw data directory."""
        return self.paths.raw

    @property
    def processed_dir(self) -> Path:
        """Processed data directory."""
        return self.paths.processed

    @property
    def stitched_dir(self) -> Path:
        """Stitched images directory."""
        return self.paths.stitched

    @property
    def deconvolved_dir(self) -> Path:
        """Deconvolved images directory."""
        return self.paths.deconvolved

    @property
    def edf_dir(self) -> Path:
        """Extended depth of focus directory."""
        return self.paths.edf

    @property
    def registered_dir(self) -> Path:
        """Registered images directory."""
        return self.paths.registered

    @property
    def signal_isolated_dir(self) -> Path:
        """Signal-isolated images directory."""
        return self.paths.signal_isolated

    @property
    def notebooks_dir(self) -> Path:
        """Working notebooks directory."""
        return self.paths.notebooks

    @property
    def config_dir(self) -> Path:
        """Configuration directory."""
        return self.paths.configs

    # -------------------------------------------------------------------------
    # Cycle management
    # -------------------------------------------------------------------------

    def add_cycle(
        self,
        name: str,
        path: str | Path | None = None,
        channels: list[str] | None = None,
        metadata: dict[str, Any] | None = None,
    ) -> None:
        """
        Add a cycle to the project.

        Parameters
        ----------
        name : str
            Cycle name (e.g., "cyc001")
        path : str or Path, optional
            Path to cycle data (relative to raw_dir)
        channels : List[str], optional
            Channel names for this cycle
        metadata : dict, optional
            Additional metadata
        """
        cycle_info = {
            "name": name,
            "path": str(path) if path else name,
            "channels": channels or [],
            "metadata": metadata or {},
        }

        # Check if cycle already exists
        for i, cyc in enumerate(self.config.cycles):
            if cyc["name"] == name:
                self.config.cycles[i] = cycle_info
                return

        self.config.cycles.append(cycle_info)
        self.save()

    def list_cycles(self) -> list[str]:
        """List all cycles in the project."""
        return [c["name"] for c in self.config.cycles]

    def get_cycle_path(self, name: str) -> Path:
        """Get the full path to a cycle's data."""
        for cyc in self.config.cycles:
            if cyc["name"] == name:
                return self.paths.raw / cyc["path"]
        raise KeyError(f"Cycle not found: {name}")

    def refresh_from_disk(self, update_cycles: bool = True) -> None:
        """Refresh project configuration to match on-disk state."""

        if update_cycles:
            disk_cycles = find_raw_cycles(self.paths.raw)
            existing_cycles = {c["name"]: c for c in self.config.cycles}
            refreshed_cycles = []

            for cycle_path in disk_cycles:
                existing = existing_cycles.get(cycle_path.name, {})
                refreshed_cycles.append(
                    {
                        "name": cycle_path.name,
                        "path": cycle_path.name,
                        "channels": existing.get("channels", []),
                        "metadata": existing.get("metadata", {}),
                    }
                )

            if refreshed_cycles != self.config.cycles:
                self.config.cycles = refreshed_cycles
                self.save()

    # -------------------------------------------------------------------------
    # Utility methods
    # -------------------------------------------------------------------------

    def setup_cuda_path(self) -> bool:
        """
        Setup CUDA PATH for CuPy on Windows.

        Returns
        -------
        bool
            True if CUDA was found and PATH updated
        """
        if platform.system() != "Windows":
            return True

        cuda_paths = [
            Path(r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.1\bin"),
            Path(r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.6\bin"),
            Path(r"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin"),
        ]

        for cuda_bin in cuda_paths:
            nvrtc = cuda_bin / "nvrtc64_121.dll"
            if cuda_bin.exists() and (nvrtc.exists() or list(cuda_bin.glob("nvrtc*.dll"))):
                os.environ["PATH"] = str(cuda_bin) + os.pathsep + os.environ.get("PATH", "")
                print(f"[OK] Added CUDA to PATH: {cuda_bin}")
                return True

        print("[WARN] CUDA not found - GPU acceleration may not work")
        return False

    def get_project_summary(self) -> str:
        """Get a summary of the project."""
        lines = [
            f"KINTSUGI Project: {self.config.name}",
            f"{'=' * 50}",
            f"Location: {self.root}",
            f"Created: {self.config.created[:10]}",
            f"Modified: {self.config.modified[:10]}",
            "",
            f"Cycles: {len(self.config.cycles)}",
        ]

        for cyc in self.config.cycles:
            lines.append(f"  - {cyc['name']}: {len(cyc.get('channels', []))} channels")

        lines.extend(
            [
                "",
                "Directories:",
                f"  Raw data: {self.paths.raw}",
                f"  Processed: {self.paths.processed}",
                f"  Notebooks: {self.paths.notebooks}",
            ]
        )

        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"KintsugiProject('{self.config.name}', root='{self.root}')"


# -----------------------------------------------------------------------------
# Convenience functions for notebook use
# -----------------------------------------------------------------------------


def _handle_processed_data_check(project: KintsugiProject, interactive: bool) -> bool:
    """
    Check for existing processed data and handle user interaction.

    Parameters
    ----------
    project : KintsugiProject
        The loaded project
    interactive : bool
        Whether to prompt the user

    Returns
    -------
    bool
        True if should continue, False if cancelled
    """
    report = scan_processed_data(project.paths)
    if not report.has_data:
        return True  # No data, nothing to do

    decisions = prompt_for_processed_data(report, interactive=interactive)
    if not decisions:
        return True  # No data or empty decisions

    return apply_processed_data_decisions(project.paths, decisions)


def init_project(
    project_dir: str | Path,
    name: str | None = None,
    description: str = "",
    mode: Literal["auto", "create", "load"] = "auto",
    refresh: bool = True,
    interactive: bool = True,
    on_existing: Literal["prompt", "keep", "clear", "cancel"] = "prompt",
    check_processed: bool = True,
) -> KintsugiProject | None:
    """
    Initialize or load a KINTSUGI project.

    This is the main entry point for notebooks. It will:
    1. Create a new project if it doesn't exist
    2. Load an existing project if it does
    3. Check for existing processed data and prompt user for action
    4. Setup CUDA paths for GPU acceleration
    5. Return the project object with all paths configured

    Parameters
    ----------
    project_dir : str or Path
        Directory for the project
    name : str, optional
        Project name (for new projects)
    description : str
        Project description (for new projects)
    mode : {"auto", "create", "load"}
        Choose whether to automatically create or load (default), force creation
        of a new project, or require loading an existing project.
    refresh : bool
        If True, synchronize the loaded project with the on-disk folder state
        (e.g., renamed cycle folders) before returning.
    interactive : bool
        If True (default), prompt the user when the directory is non-empty.
        If False, use the `on_existing` default behavior.
    on_existing : {"prompt", "keep", "clear", "cancel"}
        What to do when the directory is non-empty:
        - "prompt": Ask the user (default, requires interactive=True)
        - "keep": Keep existing files and create project around them
        - "clear": Delete existing files and start fresh
        - "cancel": Cancel project creation
    check_processed : bool
        If True (default), check for existing processed data when loading
        a project and prompt user for what to do (keep/delete per stage).

    Returns
    -------
    KintsugiProject or None
        Initialized project, or None if creation was cancelled

    Example
    -------
    >>> from kintsugi.project import init_project
    >>> project = init_project("C:/Projects/MyExperiment")
    >>>
    >>> # Use project paths
    >>> raw_images = project.raw_dir / "cyc001"
    >>> output = project.stitched_dir / "cyc001_stitched.tiff"
    """
    project_dir = Path(project_dir).resolve()
    config_file = project_dir / PROJECT_CONFIG_FILE

    if mode not in {"auto", "create", "load"}:
        raise ValueError("mode must be one of 'auto', 'create', or 'load'")

    if on_existing not in {"prompt", "keep", "clear", "cancel"}:
        raise ValueError("on_existing must be one of 'prompt', 'keep', 'clear', or 'cancel'")

    # Check if loading existing project
    if mode == "load":
        project = KintsugiProject.load(project_dir)
        project.setup_cuda_path()
        if name is not None and name != project.config.name:
            old_name = project.config.name
            project.config.name = name
            project.save()
            print(f"  Updated project name: '{old_name}' -> '{name}'")
        if refresh:
            project.refresh_from_disk()
        # Check for existing processed data
        if check_processed:
            if not _handle_processed_data_check(project, interactive):
                return None
        return project

    # Check if project already exists (auto mode)
    if mode == "auto" and config_file.exists():
        project = KintsugiProject.load(project_dir)
        project.setup_cuda_path()
        if name is not None and name != project.config.name:
            old_name = project.config.name
            project.config.name = name
            project.save()
            print(f"  Updated project name: '{old_name}' -> '{name}'")
        if refresh:
            project.refresh_from_disk()
        # Check for existing processed data
        if check_processed:
            if not _handle_processed_data_check(project, interactive):
                return None
        return project

    # Creating a new project - check for existing contents
    contents = check_directory_contents(project_dir)

    if contents["exists"] and not contents["is_empty"] and not contents["has_project_config"]:
        # Directory has files but no project config - need to decide what to do
        if on_existing == "prompt":
            action = prompt_for_existing_directory(project_dir, contents, interactive=interactive)
        else:
            action = on_existing

        if action == "cancel":
            return None
        elif action == "clear":
            deleted = clear_directory_contents(project_dir)
            print(f"Deleted {deleted} items from {project_dir}")

    # Create the project
    project = KintsugiProject.create(
        project_dir,
        name=name,
        description=description,
    )

    # Setup CUDA for GPU acceleration
    project.setup_cuda_path()

    if refresh:
        project.refresh_from_disk()

    return project


def check_directory_contents(directory: str | Path) -> dict[str, Any]:
    """
    Check if a directory has existing contents.

    Parameters
    ----------
    directory : str or Path
        Directory to check

    Returns
    -------
    dict
        Dictionary with:
        - 'exists': bool - whether directory exists
        - 'is_empty': bool - whether directory is empty
        - 'file_count': int - number of files (recursive)
        - 'dir_count': int - number of subdirectories
        - 'total_size_mb': float - total size in MB
        - 'has_project_config': bool - whether kintsugi_project.json exists
        - 'sample_files': list - sample of file names
    """
    directory = Path(directory)

    result = {
        "exists": False,
        "is_empty": True,
        "file_count": 0,
        "dir_count": 0,
        "total_size_mb": 0.0,
        "has_project_config": False,
        "sample_files": [],
    }

    if not directory.exists():
        return result

    result["exists"] = True
    result["has_project_config"] = (directory / PROJECT_CONFIG_FILE).exists()

    files = []
    dirs = []
    total_size = 0

    try:
        for item in directory.rglob("*"):
            if item.is_file():
                files.append(item)
                try:
                    total_size += item.stat().st_size
                except (OSError, PermissionError):
                    pass
            elif item.is_dir():
                dirs.append(item)
    except (OSError, PermissionError):
        pass

    result["file_count"] = len(files)
    result["dir_count"] = len(dirs)
    result["total_size_mb"] = total_size / (1024 * 1024)
    result["is_empty"] = len(files) == 0 and len(dirs) == 0
    result["sample_files"] = [f.name for f in files[:10]]

    return result


def prompt_for_existing_directory(
    directory: str | Path,
    contents: dict[str, Any],
    interactive: bool = True,
) -> str:
    """
    Prompt the user about what to do with an existing non-empty directory.

    Parameters
    ----------
    directory : str or Path
        Directory path
    contents : dict
        Result from check_directory_contents()
    interactive : bool
        If True, prompt for user input. If False, return 'keep' by default.

    Returns
    -------
    str
        One of: 'keep', 'clear', 'cancel'
    """
    directory = Path(directory)

    if contents["is_empty"]:
        return "keep"

    if contents["has_project_config"]:
        # Directory already has a KINTSUGI project - just load it
        return "keep"

    # Build information message
    print("\n" + "=" * 60)
    print("WARNING: Directory is not empty!")
    print("=" * 60)
    print(f"Location: {directory}")
    print(f"Files: {contents['file_count']}")
    print(f"Subdirectories: {contents['dir_count']}")
    print(f"Total size: {contents['total_size_mb']:.1f} MB")
    if contents["sample_files"]:
        print(f"Sample files: {', '.join(contents['sample_files'][:5])}")
    print()

    if not interactive:
        print("Non-interactive mode: keeping existing files.")
        return "keep"

    # Check if we're in a Jupyter notebook
    try:
        from IPython import get_ipython

        shell = get_ipython()
        in_notebook = shell is not None and "IPKernelApp" in shell.config
    except (ImportError, AttributeError):
        in_notebook = False

    # Prompt user for choice
    print("Options:")
    print("  [K]eep   - Keep existing files and create project around them")
    print("  [C]lear  - Delete all existing files and start fresh")
    print("  [A]bort  - Cancel project creation")
    print()

    while True:
        try:
            if in_notebook:
                # In Jupyter, use input() which shows a text field
                response = input("Choose an option [K/C/A] (default: K): ").strip().lower()
            else:
                # CLI mode
                response = input("Choose an option [K/C/A] (default: K): ").strip().lower()

            if response == "" or response == "k" or response == "keep":
                print("Keeping existing files.")
                return "keep"
            elif response == "c" or response == "clear":
                # Double-check for clear operation
                confirm = (
                    input(
                        f"Are you sure you want to DELETE all {contents['file_count']} files? [y/N]: "
                    )
                    .strip()
                    .lower()
                )
                if confirm == "y" or confirm == "yes":
                    print("Clearing directory...")
                    return "clear"
                else:
                    print("Clear cancelled. Keeping existing files.")
                    return "keep"
            elif response == "a" or response == "abort":
                print("Project creation cancelled.")
                return "cancel"
            else:
                print("Invalid option. Please enter K, C, or A.")
        except (EOFError, KeyboardInterrupt):
            print("\nProject creation cancelled.")
            return "cancel"


def clear_directory_contents(directory: str | Path, preserve_hidden: bool = True) -> int:
    """
    Clear all contents from a directory.

    Parameters
    ----------
    directory : str or Path
        Directory to clear
    preserve_hidden : bool
        If True, preserve hidden files/directories (starting with '.')

    Returns
    -------
    int
        Number of items deleted
    """
    directory = Path(directory)
    deleted = 0

    if not directory.exists():
        return 0

    for item in list(directory.iterdir()):
        if preserve_hidden and item.name.startswith("."):
            continue

        try:
            if item.is_file():
                item.unlink()
                deleted += 1
            elif item.is_dir():
                shutil.rmtree(item)
                deleted += 1
        except (OSError, PermissionError) as e:
            print(f"Warning: Could not delete {item}: {e}")

    return deleted


def find_raw_cycles(raw_dir: str | Path) -> list[Path]:
    """
    Find cycle directories in raw data folder.

    Looks for directories matching common naming patterns:
    - cyc001, cyc002, ...
    - Cycle_1, Cycle_2, ...
    - 001, 002, ...

    Parameters
    ----------
    raw_dir : str or Path
        Raw data directory to search

    Returns
    -------
    List[Path]
        List of cycle directories, sorted
    """
    raw_dir = Path(raw_dir)

    if not raw_dir.exists():
        return []

    cycles = []
    for item in raw_dir.iterdir():
        if item.is_dir():
            name = item.name.lower()
            # Match common cycle naming patterns
            if (
                name.startswith("cyc")
                or name.startswith("cycle")
                or name.isdigit()
                or (len(name) >= 3 and name[:3].isdigit())
            ):
                cycles.append(item)

    # Sort naturally
    try:
        from natsort import natsorted

        return natsorted(cycles, key=lambda p: p.name)
    except ImportError:
        return sorted(cycles, key=lambda p: p.name)


def _tile_count_to_grid(n_tiles: int) -> tuple[int, int]:
    """
    Convert a tile count to (rows, cols) by finding the factor pair closest to square.

    Parameters
    ----------
    n_tiles : int
        Total number of tiles

    Returns
    -------
    tuple[int, int]
        (rows, cols) where rows <= cols
    """
    if n_tiles <= 0:
        return (1, 1)

    import math

    sqrt = int(math.isqrt(n_tiles))
    # Search downward from sqrt for the largest factor
    for i in range(sqrt, 0, -1):
        if n_tiles % i == 0:
            return (i, n_tiles // i)
    return (1, n_tiles)
