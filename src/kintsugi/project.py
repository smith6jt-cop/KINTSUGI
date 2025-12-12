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
from typing import Any

# Project configuration filename
PROJECT_CONFIG_FILE = "kintsugi_project.json"
PROJECT_VERSION = "1.0.0"

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

    def summary(self) -> str:
        """Generate a human-readable summary."""
        lines = []
        if not self.has_data:
            lines.append("No existing data found.")
            return "\n".join(lines)

        lines.append(f"Found {self.image_count} image files ({self.total_size_mb:.1f} MB total)")

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


def scan_existing_data(
    directory: str | Path,
    max_depth: int = 3,
    sample_count: int = 5,
) -> ExistingDataReport:
    """
    Scan a directory for existing data files.

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
    cycle_pattern = re.compile(r"^cyc\d+$|^cycle\d+$|^Cycle\d+$", re.IGNORECASE)
    filename_bases = set()

    def scan_dir(path: Path, depth: int = 0):
        if depth > max_depth:
            return
        try:
            for item in path.iterdir():
                if item.is_dir():
                    # Check for cycle folder pattern
                    if cycle_pattern.match(item.name):
                        report.cycle_folders.append(item.name)
                    scan_dir(item, depth + 1)
                elif item.is_file():
                    # Check file extension
                    suffix = item.suffix.lower()
                    # Handle .ome.tif specially
                    if item.name.lower().endswith(".ome.tif") or item.name.lower().endswith(
                        ".ome.tiff"
                    ):
                        suffix = ".ome.tif"

                    if suffix in IMAGE_EXTENSIONS:
                        image_files.append(item)
                        # Extract filename pattern (remove numbers)
                        base = re.sub(r"\d+", "#", item.stem)
                        filename_bases.add(base)
                    elif suffix not in {".pyc", ".log", ".json", ".md", ".txt", ".gitignore"}:
                        other_files.append(item)
        except PermissionError:
            pass

    scan_dir(directory)

    report.image_files = image_files
    report.image_count = len(image_files)
    report.other_files = other_files
    report.has_data = len(image_files) > 0 or len(other_files) > 0
    report.filename_patterns = sorted(filename_bases)[:10]

    # Calculate total size
    try:
        report.total_size_mb = sum(f.stat().st_size for f in image_files) / (1024 * 1024)
    except (OSError, PermissionError):
        pass

    # Sort cycle folders naturally
    report.cycle_folders = sorted(set(report.cycle_folders))

    # Extract metadata from sample images
    if image_files and sample_count > 0:
        sample_files = image_files[: min(sample_count, len(image_files))]
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
    stitched: Path  # Stitched images
    corrected: Path  # Illumination-corrected images
    registered: Path  # Registered multi-cycle images
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
            corrected=root / "data" / "processed" / "corrected",
            registered=root / "data" / "processed" / "registered",
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
        adopt_existing_data: bool = False,
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
        adopt_existing_data : bool
            If True and existing data is found, automatically organize it
            into the project structure (move to data/raw/)

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

        # Create directory structure (uses exist_ok=True, never overwrites)
        project.paths.create_all()

        # Handle existing data adoption
        if adopt_existing_data and existing_data_report.has_data:
            project._adopt_existing_data(existing_data_report)

        # Copy notebooks if requested
        if setup_notebooks:
            project.setup_notebooks()

        # Create default config files
        project._create_default_configs()

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
        print("  +-- .claude/          <- Claude Code config")
        print("  +-- .vscode/          <- VS Code config")
        print()

        if existing_data_report.has_data:
            print("Existing data detected:")
            print(f"  {existing_data_report.image_count} image files")
            if existing_data_report.cycle_folders:
                print(f"  Cycles: {', '.join(existing_data_report.cycle_folders)}")
            print()

        print("Next steps:")
        if existing_data_report.has_data and not adopt_existing_data:
            print("  1. Review existing data locations")
            print("  2. Move/copy raw images to data/raw/ if needed")
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
                f"No KINTSUGI project found at {root}. " f"Use create() to create a new project."
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
        support_dirs = ["Kreg", "Kstitch", "Kview2", "KDecon", "instanseg"]
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

    def _create_default_configs(self) -> None:
        """Create default configuration files."""
        # Default processing parameters
        default_params = {
            "illumination_correction": {
                "method": "BaSiC",
                "working_size": 128,
                "max_iterations": 500,
            },
            "stitching": {
                "overlap_percent": 10,
                "blend_mode": "linear",
            },
            "registration": {
                "method": "rigid",
                "reference_cycle": 1,
            },
            "segmentation": {
                "method": "instanseg",
                "cell_diameter": 30,
            },
        }

        params_file = self.paths.configs / "default_parameters.json"
        with open(params_file, "w") as f:
            json.dump(default_params, f, indent=2)

        # Create Claude Code configuration
        self._create_claude_config()

        # Create VS Code configuration
        self._create_vscode_config()

    def _create_claude_config(self) -> None:
        """Create Claude Code MCP configuration."""
        claude_dir = self.paths.root / ".claude"
        claude_dir.mkdir(exist_ok=True)

        settings_file = claude_dir / "settings.local.json"
        if not settings_file.exists():
            claude_config = {
                "mcpServers": {
                    "kintsugi": {
                        "command": "kintsugi",
                        "args": ["mcp", "start"],
                        "cwd": str(self.paths.root),
                    }
                }
            }
            with open(settings_file, "w") as f:
                json.dump(claude_config, f, indent=2)
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
    def registered_dir(self) -> Path:
        """Registered images directory."""
        return self.paths.registered

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


def init_project(
    project_dir: str | Path,
    name: str | None = None,
    description: str = "",
) -> KintsugiProject:
    """
    Initialize or load a KINTSUGI project.

    This is the main entry point for notebooks. It will:
    1. Create a new project if it doesn't exist
    2. Load an existing project if it does
    3. Setup CUDA paths for GPU acceleration
    4. Return the project object with all paths configured

    Parameters
    ----------
    project_dir : str or Path
        Directory for the project
    name : str, optional
        Project name (for new projects)
    description : str
        Project description (for new projects)

    Returns
    -------
    KintsugiProject
        Initialized project

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

    if config_file.exists():
        project = KintsugiProject.load(project_dir)
    else:
        project = KintsugiProject.create(
            project_dir,
            name=name,
            description=description,
        )

    # Setup CUDA for GPU acceleration
    project.setup_cuda_path()

    return project


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
