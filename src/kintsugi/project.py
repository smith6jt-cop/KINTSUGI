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
from dataclasses import dataclass, field, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Any, Union
import warnings


# Project configuration filename
PROJECT_CONFIG_FILE = "kintsugi_project.json"
PROJECT_VERSION = "1.0.0"


@dataclass
class ProjectPaths:
    """Structured paths for a KINTSUGI project."""

    # Root project directory
    root: Path

    # Input data
    raw: Path              # Raw images from microscope

    # Processing outputs
    processed: Path        # Base processed directory
    stitched: Path         # Stitched images
    corrected: Path        # Illumination-corrected images
    registered: Path       # Registered multi-cycle images
    segmented: Path        # Segmentation masks
    analysis: Path         # Analysis outputs (tables, plots)

    # Working files
    notebooks: Path        # Working copies of notebooks
    configs: Path          # Project-specific configurations
    logs: Path             # Processing logs
    cache: Path            # Temporary/cache files

    # Metadata
    meta: Path             # Metadata files

    @classmethod
    def from_root(cls, root: Union[str, Path]) -> "ProjectPaths":
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

    def to_dict(self) -> Dict[str, str]:
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
    parameters: Dict[str, Any] = field(default_factory=dict)

    # Cycle information
    cycles: List[Dict[str, Any]] = field(default_factory=list)

    # Channel metadata
    channels: List[Dict[str, str]] = field(default_factory=list)

    def to_dict(self) -> Dict[str, Any]:
        """Convert to dictionary."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: Dict[str, Any]) -> "ProjectConfig":
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
        root: Union[str, Path],
        config: Optional[ProjectConfig] = None,
        kintsugi_path: Optional[Union[str, Path]] = None,
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
    def create(
        cls,
        root: Union[str, Path],
        name: Optional[str] = None,
        description: str = "",
        kintsugi_path: Optional[Union[str, Path]] = None,
        setup_notebooks: bool = True,
    ) -> "KintsugiProject":
        """
        Create a new KINTSUGI project.

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
                UserWarning
            )
            return cls.load(root)

        # Create project
        config = ProjectConfig(
            name=name or root.name,
            description=description,
        )

        project = cls(root, config, kintsugi_path)

        # Create directory structure
        project.paths.create_all()

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
        print(f"  +-- data/")
        print(f"  |   +-- raw/          <- Put raw images here")
        print(f"  |   +-- processed/    <- Outputs go here")
        print(f"  +-- notebooks/        <- Working notebooks")
        print(f"  +-- configs/          <- Processing configs")
        print(f"  +-- meta/             <- Metadata files")
        print(f"  +-- logs/             <- Processing logs")

        return project

    @classmethod
    def load(cls, root: Union[str, Path]) -> "KintsugiProject":
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
                f"No KINTSUGI project found at {root}. "
                f"Use create() to create a new project."
            )

        with open(config_file, "r") as f:
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

    def setup_notebooks(self, overwrite: bool = False) -> List[Path]:
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

        # Notebooks to copy (main workflow)
        workflow_notebooks = [
            "1_Single_Channel_Eval.ipynb",
            "2_Cycle_Processing.ipynb",
            "3_Signal_Isolation.ipynb",
            "4_Segmentation_Analysis.ipynb",
            "5_DL_Channel_Refinement.ipynb",
        ]

        copied = []
        for nb_name in workflow_notebooks:
            src = source_dir / nb_name
            dst = dest_dir / nb_name

            if not src.exists():
                continue

            if dst.exists() and not overwrite:
                print(f"  Skipping {nb_name} (exists)")
                continue

            shutil.copy2(src, dst)
            copied.append(dst)
            print(f"  Copied {nb_name}")

        # Copy supporting modules
        support_dirs = ["Kreg", "Kstitch", "Kview2", "KDecon", "instanseg"]
        for dir_name in support_dirs:
            src_dir = source_dir / dir_name
            dst_dir = dest_dir / dir_name

            if src_dir.exists() and not dst_dir.exists():
                shutil.copytree(src_dir, dst_dir)
                print(f"  Copied {dir_name}/")

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
        path: Optional[Union[str, Path]] = None,
        channels: Optional[List[str]] = None,
        metadata: Optional[Dict[str, Any]] = None,
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

    def list_cycles(self) -> List[str]:
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
            f"",
            f"Cycles: {len(self.config.cycles)}",
        ]

        for cyc in self.config.cycles:
            lines.append(f"  - {cyc['name']}: {len(cyc.get('channels', []))} channels")

        lines.extend([
            f"",
            f"Directories:",
            f"  Raw data: {self.paths.raw}",
            f"  Processed: {self.paths.processed}",
            f"  Notebooks: {self.paths.notebooks}",
        ])

        return "\n".join(lines)

    def __repr__(self) -> str:
        return f"KintsugiProject('{self.config.name}', root='{self.root}')"


# -----------------------------------------------------------------------------
# Convenience functions for notebook use
# -----------------------------------------------------------------------------

def init_project(
    project_dir: Union[str, Path],
    name: Optional[str] = None,
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


def find_raw_cycles(raw_dir: Union[str, Path]) -> List[Path]:
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
            if (name.startswith("cyc") or
                name.startswith("cycle") or
                name.isdigit() or
                (len(name) >= 3 and name[:3].isdigit())):
                cycles.append(item)

    # Sort naturally
    try:
        from natsort import natsorted
        return natsorted(cycles, key=lambda p: p.name)
    except ImportError:
        return sorted(cycles, key=lambda p: p.name)
