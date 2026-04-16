"""
Dependency validation module for KINTSUGI.

Provides runtime checking of all external dependencies including
native libraries (libvips), CUDA/GPU, and Python packages.

Key features:
- `require()`: Use at the top of notebooks to ensure dependencies are installed
- `check_dependencies()`: Full dependency check for CLI
- `install_optional()`: Install optional dependency groups via pip

Note: Java/Maven dependencies are deprecated and no longer required.
The pipeline now uses pure Python implementations for all processing.
"""

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Literal


class DependencyStatus(Enum):
    """Status of a dependency check."""

    OK = "ok"
    MISSING = "missing"
    VERSION_MISMATCH = "version_mismatch"
    ERROR = "error"
    OPTIONAL_MISSING = "optional_missing"
    # Transient: library will load fine once conda activation scripts take effect
    # (e.g., libstdc++ version mismatch resolved by reactivating the env).
    PENDING = "pending"


def _is_libstdcxx_load_error(message: str) -> bool:
    """Return True if `message` looks like a libstdc++ version mismatch at load time.

    Matches the signatures produced when a Python extension is linked against a
    newer libstdc++ than the one currently on LD_LIBRARY_PATH — typically
    right after conda activation scripts are deployed but before the env is
    reactivated. Those failures all clear on the next `conda activate`.
    """
    return any(tok in message for tok in ("CXXABI_", "GLIBCXX_", "libstdc++"))


@dataclass
class DependencyResult:
    """Result of checking a single dependency."""

    name: str
    status: DependencyStatus
    version: str | None = None
    required_version: str | None = None
    message: str = ""
    is_optional: bool = False
    details: dict = field(default_factory=dict)


# =============================================================================
# Optional Dependency Groups
# =============================================================================
# These map to pyproject.toml optional dependencies and define what each
# notebook/feature requires.

OPTIONAL_GROUPS = {
    "gpu": {
        "description": "GPU acceleration (PyTorch + CuPy for CUDA)",
        "packages": ["torch", "torchvision", "cupy"],
        # CUDA runtime libraries (libcufft, etc.) must be installed via conda;
        # cupy-cuda12x (pip) only provides Python bindings, not the native libs
        # Use CUDA 12.8 for broad GPU support including Blackwell (B200, sm_100)
        # cuda-cudart-dev provides headers needed for CuPy JIT compilation
        # Headers must be copied to $CONDA_PREFIX/include for CuPy to find them
        "install_cmd": "conda install cuda-libraries cuda-cudart-dev -c nvidia -y && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null; pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128 && pip install cupy-cuda12x",
        "conda_cmd": "conda install pytorch torchvision pytorch-cuda=12.8 cupy cuda-libraries cuda-cudart-dev -c pytorch -c nvidia -c conda-forge && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null",
    },
    "viz": {
        "description": "Napari interactive visualization",
        "packages": ["napari", "magicgui"],
        # Install packages directly to avoid PyPI 'kintsugi' name collision
        "install_cmd": "pip install napari magicgui napari-console napari-crop napari-simpleitk-image-processing napari-skimage-regionprops",
        "conda_cmd": "conda install napari pyqt -c conda-forge",
    },
    "dl": {
        "description": "Deep learning segmentation (InstanSeg)",
        "packages": ["torch", "instanseg", "kornia"],
        "install_cmd": "pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128 && pip install instanseg instanseg-torch kornia",
        "conda_cmd": "conda install pytorch torchvision pytorch-cuda=12.8 -c pytorch -c nvidia -c conda-forge && pip install instanseg instanseg-torch kornia",
    },
    "analysis": {
        "description": "Spatial analysis (scanpy, scimap)",
        "packages": ["scanpy", "anndata", "phenograph", "scimap", "skan", "networkx"],
        "install_cmd": "pip install scanpy anndata phenograph scimap umap-learn hdbscan skan networkx",
        "conda_cmd": "conda install scanpy anndata networkx -c conda-forge && pip install scimap phenograph skan",
    },
    "bio": {
        "description": "Bio formats I/O (OME-TIFF, LIF, etc.)",
        "packages": ["bioio", "bioio-ome-tiff", "ome-zarr", "slideio", "readlif"],
        "install_cmd": "pip install bioio bioio-ome-tiff ome-zarr slideio readlif",
    },
    "claude": {
        "description": "Claude Code MCP integration",
        "packages": ["mcp", "httpx"],
        "install_cmd": "pip install mcp httpx",
    },
    "dev": {
        "description": "Development tools (pytest, ruff, black, mypy)",
        "packages": ["pytest", "ruff", "black", "mypy", "pre_commit"],
        "install_cmd": "pip install pytest pytest-asyncio pytest-cov pytest-xdist ruff black mypy pre-commit commitizen",
    },
    "docs": {
        "description": "Documentation (Sphinx)",
        "packages": ["sphinx", "sphinx_rtd_theme", "myst_parser"],
        "install_cmd": "pip install sphinx sphinx-rtd-theme myst-parser",
    },
    "kronos": {
        "description": "KRONOS foundation model for spatial proteomics",
        "packages": ["torch", "h5py", "umap"],
        "install_cmd": "pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128 && pip install h5py umap-learn scikit-learn anndata scanpy tifffile",
        "conda_cmd": "conda install pytorch torchvision pytorch-cuda=12.8 -c pytorch -c nvidia -c conda-forge && pip install h5py umap-learn scikit-learn anndata scanpy tifffile",
        "note": "Also requires: git clone https://github.com/mahmoodlab/KRONOS.git && pip install -e KRONOS",
    },
    "denoise": {
        "description": "Advanced denoising (N2V, CARE)",
        "packages": ["torch"],
        "install_cmd": "pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128",
        "conda_cmd": "conda install pytorch torchvision pytorch-cuda=12.8 -c pytorch -c nvidia -c conda-forge",
    },
    "optimize": {
        "description": "Automated parameter optimization (Optuna + SMO)",
        "packages": ["optuna", "smo", "joblib"],
        "install_cmd": "pip install optuna smo joblib",
    },
    "rapids": {
        "description": "RAPIDS GPU-accelerated data science",
        "packages": ["cudf", "cuml"],
        "install_cmd": "pip install cudf-cu12 cuml-cu12 --extra-index-url=https://pypi.nvidia.com",
        "conda_cmd": "conda install cudf cuml cugraph rmm -c rapidsai -c conda-forge -c nvidia",
    },
    "workflow": {
        "description": "Snakemake workflow orchestration (HPC/SLURM)",
        "packages": ["snakemake", "snakemake_executor_plugin_slurm"],
        "install_cmd": "pip install snakemake snakemake-executor-plugin-slurm snakemake-executor-plugin-slurm-jobstep",
    },
    "full": {
        "description": "All optional features",
        "packages": [],  # Composite group
        "install_cmd": "conda install cuda-libraries cuda-cudart-dev -c nvidia -y && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null; pip install torch torchvision --index-url https://download.pytorch.org/whl/cu128 && pip install cupy-cuda12x napari magicgui instanseg instanseg-torch kornia scanpy anndata phenograph scimap bioio bioio-ome-tiff ome-zarr slideio readlif",
    },
}

# Mapping of notebooks to required optional groups
NOTEBOOK_REQUIREMENTS = {
    "1_Single_Channel_Eval": ["gpu"],
    "2_Cycle_Processing": ["gpu"],
    "3_Signal_Isolation": [],  # Core only
    "4_Segmentation_Analysis": ["dl", "viz", "analysis"],
    "5_Cluster_Analysis": ["analysis"],
    "5_DL_Channel_Refinement": [],  # Core only
    "Image_Registration_Workflow": [],  # Core only
    "Vessel_Analysis": ["viz", "analysis"],
    "2.5_Vessel_3D_Segmentation": ["analysis"],
}


def _find_constraints_file() -> Path | None:
    """Find the constraints.txt file in the project root."""
    root = _find_project_root()
    if root and (root / "constraints.txt").exists():
        return root / "constraints.txt"
    return None


def _find_project_root() -> Path | None:
    """Find KINTSUGI project root containing pyproject.toml.

    Walks up from the kintsugi package directory to find the repo root.
    Returns None if pyproject.toml cannot be found (e.g., pip-installed wheel).
    """
    try:
        import kintsugi

        pkg_dir = Path(kintsugi.__file__).parent  # src/kintsugi/
        for parent in [pkg_dir.parent.parent, pkg_dir.parent]:
            if (parent / "pyproject.toml").exists():
                return parent
    except Exception:
        pass
    return None


def detect_hpc() -> bool:
    """Detect if running on an HPC cluster (typically SLURM-based).

    Uses strong indicators to avoid false positives on non-HPC systems:
    - SLURM environment (``SLURM_CONF``)
    - Environment modules (``MODULESHOME``)
    - Presence of SLURM commands (``srun``, ``sbatch``) on PATH
    """
    if os.environ.get("SLURM_CONF") or os.environ.get("MODULESHOME"):
        return True
    if shutil.which("srun") or shutil.which("sbatch"):
        return True
    return False


def deploy_activation_scripts() -> bool:
    """Deploy conda activate/deactivate scripts for HPC LD_LIBRARY_PATH fix.

    Copies ``envs/activate.d/env_vars.sh`` and ``envs/deactivate.d/env_vars.sh``
    from the KINTSUGI project root into the active conda environment's
    activation directories. These scripts prepend ``$CONDA_PREFIX/lib`` to
    ``LD_LIBRARY_PATH`` so conda's ``libstdc++.so.6`` is loaded before the
    system's (which is often too old on HPC clusters).

    Returns True on success, False if unable to deploy.
    """
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return False

    project_root = _find_project_root()
    if not project_root:
        return False

    activate_src = project_root / "envs" / "activate.d" / "env_vars.sh"
    deactivate_src = project_root / "envs" / "deactivate.d" / "env_vars.sh"

    if not activate_src.exists() or not deactivate_src.exists():
        return False

    activate_dst = Path(conda_prefix) / "etc" / "conda" / "activate.d"
    deactivate_dst = Path(conda_prefix) / "etc" / "conda" / "deactivate.d"

    try:
        activate_dst.mkdir(parents=True, exist_ok=True)
        deactivate_dst.mkdir(parents=True, exist_ok=True)
        shutil.copy2(activate_src, activate_dst / "env_vars.sh")
        shutil.copy2(deactivate_src, deactivate_dst / "env_vars.sh")
        # Make executable
        (activate_dst / "env_vars.sh").chmod(0o755)
        (deactivate_dst / "env_vars.sh").chmod(0o755)
        return True
    except OSError:
        return False


def ensure_hpc_activation_scripts() -> None:
    """Auto-deploy HPC activation scripts if missing (idempotent, <1ms).

    Called at import time on HPC systems. Only acts when:
    1. Running on an HPC cluster (SLURM detected)
    2. A conda environment is active
    3. The activation scripts are NOT already deployed

    This is the permanent fix for the libstdc++ version mismatch on HiPerGator.
    Without these scripts, conda's ``libstdc++.so.6`` is not on
    ``LD_LIBRARY_PATH``, causing import failures for scikit-learn, matplotlib,
    pyvips, and stackview.
    """
    conda_prefix = os.environ.get("CONDA_PREFIX")
    if not conda_prefix:
        return

    # Only on HPC (SLURM or module system detected)
    if not (os.environ.get("SLURM_CONF") or os.environ.get("MODULESHOME")):
        return

    # Already deployed — fast path
    target = Path(conda_prefix) / "etc" / "conda" / "activate.d" / "env_vars.sh"
    if target.exists():
        return

    # Deploy and notify user to re-activate
    if deploy_activation_scripts():
        import sys

        print(
            "KINTSUGI: Deployed HPC activation scripts. "
            "Run 'conda deactivate && conda activate KINTSUGI' to apply.",
            file=sys.stderr,
        )


def _pre_install_pims() -> None:
    """Pre-install pims with --no-build-isolation to work around setuptools bug.

    pims uses a legacy setup.py that references the removed ``install_layout``
    attribute, causing build failures with setuptools >= 69. Building without
    isolation uses the already-installed setuptools (which may be patched or
    older), avoiding the error.  This is a no-op if pims is already installed.
    """
    try:
        import pims  # noqa: F401

        return  # Already installed
    except ImportError:
        pass

    import logging

    logger = logging.getLogger(__name__)
    constraints = _find_constraints_file()
    cmd = ["pip", "install", "pims", "--no-build-isolation"]
    if constraints:
        cmd = ["pip", "install", "-c", str(constraints), "pims", "--no-build-isolation"]
    try:
        subprocess.run(cmd, check=True, capture_output=True)
    except subprocess.CalledProcessError as exc:
        stderr = exc.stderr.decode() if exc.stderr else ""
        logger.warning("Pre-install of pims failed (non-fatal): %s", stderr[:200])


class MissingDependencyError(Exception):
    """Raised when required dependencies are not installed."""

    def __init__(self, missing: list[str], groups: list[str], install_hint: str):
        self.missing = missing
        self.groups = groups
        self.install_hint = install_hint
        super().__init__(self._format_message())

    def _format_message(self) -> str:
        lines = [
            "",
            "=" * 70,
            "MISSING DEPENDENCIES",
            "=" * 70,
            "",
            "The following packages are required but not installed:",
            "",
        ]
        for pkg in self.missing:
            lines.append(f"  - {pkg}")
        lines.extend(
            [
                "",
                "To install, run:",
                "",
                f"  {self.install_hint}",
                "",
                "Or install all optional dependencies:",
                "",
                "  pip install 'kintsugi[full]'",
                "",
                "=" * 70,
            ]
        )
        return "\n".join(lines)


def _check_package(package: str) -> bool:
    """Check if a package is importable."""
    import_name = package.replace("-", "_")
    # Handle special cases
    if package == "scikit-image":
        import_name = "skimage"
    elif package == "scikit-learn":
        import_name = "sklearn"
    elif package == "opencv-contrib-python-headless":
        import_name = "cv2"
    elif package == "cupy-cuda12x":
        import_name = "cupy"

    try:
        __import__(import_name)
        return True
    except ImportError:
        return False


def require(
    *groups: str,
    notebook: str | None = None,
    strict: bool = True,
) -> dict[str, bool]:
    """
    Check that required optional dependencies are installed.

    Use this at the top of notebooks to ensure all required packages
    are available before running. If dependencies are missing, raises
    a clear error with installation instructions.

    Parameters
    ----------
    *groups : str
        Optional dependency groups to check: 'gpu', 'viz', 'dl', 'analysis', 'bio'
    notebook : str, optional
        Notebook name to auto-detect required groups (e.g., '4_Segmentation_Analysis')
    strict : bool, default True
        If True, raise MissingDependencyError when dependencies are missing.
        If False, print a warning and return status dict.

    Returns
    -------
    dict[str, bool]
        Dictionary mapping package names to availability status.

    Raises
    ------
    MissingDependencyError
        If strict=True and any required packages are missing.

    Examples
    --------
    At the top of a notebook:

    >>> from kintsugi.deps import require
    >>> require('gpu', 'viz')  # Check specific groups

    Or auto-detect from notebook name:

    >>> require(notebook='4_Segmentation_Analysis')  # Checks dl, viz, analysis
    """
    # Auto-detect groups from notebook name
    if notebook:
        detected_groups = NOTEBOOK_REQUIREMENTS.get(notebook, [])
        groups = tuple(set(groups) | set(detected_groups))

    if not groups:
        return {}

    # Expand 'full' group
    if "full" in groups:
        groups = tuple(set(groups) | {"gpu", "viz", "dl", "analysis", "bio"} - {"full"})

    # Collect all packages to check
    packages_to_check = set()
    for group in groups:
        if group in OPTIONAL_GROUPS:
            packages_to_check.update(OPTIONAL_GROUPS[group]["packages"])

    # Check each package
    status = {}
    missing = []
    for pkg in packages_to_check:
        available = _check_package(pkg)
        status[pkg] = available
        if not available:
            missing.append(pkg)

    if missing:
        # Build install hint
        install_cmds = []
        for group in groups:
            if group in OPTIONAL_GROUPS:
                install_cmds.append(OPTIONAL_GROUPS[group]["install_cmd"])
        install_hint = " && ".join(install_cmds) if install_cmds else "pip install 'kintsugi[full]'"

        if strict:
            raise MissingDependencyError(missing, list(groups), install_hint)
        else:
            print(f"\n⚠️  Warning: Missing optional dependencies: {', '.join(missing)}")
            print(f"   Install with: {install_hint}\n")

    return status


def _inject_constraints(cmd: str) -> str:
    """Inject -c constraints.txt into pip install commands if the file exists."""
    constraints = _find_constraints_file()
    if constraints and "pip install" in cmd:
        # Replace each 'pip install' with 'pip install -c constraints.txt'
        return cmd.replace("pip install ", f"pip install -c {constraints} ")
    return cmd


# Groups that require CUDA-enabled PyTorch
_TORCH_GROUPS = frozenset({"gpu", "dl", "denoise", "kronos"})


def _pre_install_guard(group: str) -> list[str]:
    """Check for conditions that will cause install to break the environment.

    Returns list of warning messages. Empty list means safe to proceed.
    """
    warnings = []

    # Guard 1: If installing a torch-using group (not gpu), check if CUDA torch exists
    if group in _TORCH_GROUPS and group != "gpu":
        try:
            import torch

            if not hasattr(torch.version, "cuda") or torch.version.cuda is None:
                warnings.append(
                    f"CPU-only PyTorch {torch.__version__} detected. "
                    f"Run 'kintsugi install gpu' first for CUDA-enabled PyTorch."
                )
        except ImportError:
            pass  # No torch yet — install will add it

    # Guard 2: If installing analysis/bio, verify numpy constraint will hold
    if group in {"analysis", "bio", "kronos"}:
        try:
            import numpy as np
            from packaging.version import parse

            if parse(np.__version__) >= parse("2.0.0"):
                warnings.append(
                    f"numpy {np.__version__} >= 2.0 detected! "
                    f"This is incompatible with KINTSUGI core. "
                    f"Fix with: pip install 'numpy>=1.24,<2.0'"
                )
        except ImportError:
            pass

    return warnings


def install_optional(
    group: Literal[
        "gpu",
        "viz",
        "dl",
        "analysis",
        "bio",
        "claude",
        "dev",
        "docs",
        "kronos",
        "denoise",
        "optimize",
        "rapids",
        "workflow",
        "full",
    ],
    use_conda: bool = False,
) -> bool:
    """
    Install optional dependency group.

    Parameters
    ----------
    group : str
        Group to install: 'gpu', 'viz', 'dl', 'analysis', 'bio', or 'full'
    use_conda : bool, default False
        Use conda instead of pip where available.

    Returns
    -------
    bool
        True if installation succeeded.

    Examples
    --------
    >>> from kintsugi.deps import install_optional
    >>> install_optional('gpu')  # Install GPU support
    >>> install_optional('viz', use_conda=True)  # Install Napari via conda
    """
    if group not in OPTIONAL_GROUPS:
        print(f"Unknown group: {group}")
        print(f"Available groups: {', '.join(OPTIONAL_GROUPS.keys())}")
        return False

    # Pre-install guard checks
    guard_warnings = _pre_install_guard(group)
    for warning in guard_warnings:
        print(f"\n⚠️  Warning: {warning}")

    info = OPTIONAL_GROUPS[group]
    print(f"\nInstalling {group}: {info['description']}")
    print("-" * 50)

    if use_conda and "conda_cmd" in info:
        cmd = info["conda_cmd"]
    else:
        cmd = info["install_cmd"]

    # Inject constraints file for pip commands
    cmd = _inject_constraints(cmd)

    print(f"Running: {cmd}\n")

    try:
        subprocess.run(cmd, shell=True, check=True)
        print(f"\n✓ Successfully installed {group} dependencies")
        return True
    except subprocess.CalledProcessError as e:
        print(f"\n✗ Failed to install {group}: {e}")
        return False


# =============================================================================
# Full Dependency Checker (for CLI)
# =============================================================================


class DependencyChecker:
    """
    Validates all KINTSUGI dependencies at runtime.

    Checks:
    - Python packages (core and optional)
    - Native libraries (libvips)
    - CUDA/GPU availability
    """

    def __init__(self):
        self.results: list[DependencyResult] = []
        # Set by _detect_libstdcxx_pending() at the top of check_all(). When True,
        # library-loading failures that match libstdc++ version-mismatch signatures
        # are reported as PENDING instead of MISSING/ERROR.
        self._libstdcxx_pending: bool = False

    def check_all(self, verbose: bool = True) -> dict:
        """
        Run all dependency checks.

        Parameters
        ----------
        verbose : bool
            Print status messages during checking.

        Returns
        -------
        dict
            Summary of all dependency checks.
        """
        self.results = []
        self._libstdcxx_pending = False

        if verbose:
            print("=" * 60)
            print("KINTSUGI Dependency Check")
            print("=" * 60)

        # Detect stale LD_LIBRARY_PATH before running cascade-prone checks, so
        # transient libstdc++ load failures get reported as PENDING instead of
        # as a wall of ERROR lines.
        self._detect_libstdcxx_pending()

        # Core Python packages
        self._check_core_packages(verbose)

        # Native libraries
        self._check_libvips(verbose)

        # GPU/CUDA
        self._check_cuda(verbose)

        # Optional packages
        self._check_optional_packages(verbose)

        # Environment safety (deep validation)
        self._check_environment_safety(verbose)

        # Summary
        summary = self._generate_summary(verbose)

        return summary

    def _check_core_packages(self, verbose: bool):
        """Check core Python package dependencies."""
        if verbose:
            print("\n[Core Python Packages]")

        core_packages = [
            ("numpy", "1.24.0"),
            ("scipy", "1.10.0"),
            ("pandas", "2.0.0"),
            ("scikit-image", "0.21.0"),
            ("scikit-learn", "1.3.0"),
            ("opencv-python", None),
            ("tifffile", "2023.7.0"),
            ("matplotlib", "3.7.0"),
            ("zarr", "2.16.0"),
            ("dask", "2023.7.0"),
            ("tqdm", "4.65.0"),
            ("click", "8.1.0"),
            ("rich", "13.0.0"),
            ("pyvips", "2.2.0"),
            ("ome-types", "0.5.0"),
            ("stackview", "0.18.0"),
        ]

        for package, min_version in core_packages:
            result = self._check_python_package(package, min_version, optional=False)
            self.results.append(result)
            if verbose:
                self._print_result(result)

    def _check_optional_packages(self, verbose: bool):
        """Check optional Python package dependencies by group."""
        if verbose:
            print("\n[Optional Packages]")

        # Group packages with their group names
        # Note: torch/torchvision/cupy are checked in _check_cuda (GPU/CUDA section)
        # Note: pyvips/ome-types/stackview are checked in _check_core_packages (base env)
        optional_packages = [
            # Visualization
            ("napari", "0.4.19", "viz"),
            ("magicgui", "0.7.0", "viz"),
            # Deep learning
            ("instanseg", "0.0.2", "dl"),
            ("kornia", "0.7.0", "dl"),
            # Analysis
            ("scanpy", "1.9.0", "analysis"),
            ("anndata", "0.9.0", "analysis"),
            ("phenograph", "1.5.0", "analysis"),
            ("scimap", None, "analysis"),
            # Bio formats
            ("bioio", "1.0.0", "bio"),
            ("ome-zarr", None, "bio"),
            ("slideio", None, "bio"),
        ]

        for package, min_version, group in optional_packages:
            result = self._check_python_package(package, min_version, optional=True)
            result.details["group"] = group
            self.results.append(result)
            if verbose:
                self._print_result(result)

    def _check_python_package(
        self, package: str, min_version: str | None, optional: bool = False
    ) -> DependencyResult:
        """Check if a Python package is installed and meets version requirements."""
        import_name = package.replace("-", "_").replace("opencv_python", "cv2")
        if package == "scikit-image":
            import_name = "skimage"
        elif package == "scikit-learn":
            import_name = "sklearn"
        elif package == "jpype1":
            import_name = "jpype"

        try:
            module = __import__(import_name)
            version = getattr(module, "__version__", None)

            if version is None:
                try:
                    from importlib.metadata import version as get_version

                    version = get_version(package)
                except Exception:
                    version = "unknown"

            # Version comparison
            if min_version and version != "unknown":
                from packaging.version import parse

                if parse(version) < parse(min_version):
                    return DependencyResult(
                        name=package,
                        status=DependencyStatus.VERSION_MISMATCH,
                        version=version,
                        required_version=min_version,
                        message=f"Version {version} < required {min_version}",
                        is_optional=optional,
                    )

            return DependencyResult(
                name=package,
                status=DependencyStatus.OK,
                version=version,
                required_version=min_version,
                message="",
                is_optional=optional,
            )

        except ImportError as e:
            err_msg = str(e)
            if self._libstdcxx_pending and _is_libstdcxx_load_error(err_msg):
                return DependencyResult(
                    name=package,
                    status=DependencyStatus.PENDING,
                    required_version=min_version,
                    message="Will verify after env reactivation (stale LD_LIBRARY_PATH)",
                    is_optional=optional,
                )
            status = DependencyStatus.OPTIONAL_MISSING if optional else DependencyStatus.MISSING
            return DependencyResult(
                name=package,
                status=status,
                required_version=min_version,
                message=err_msg,
                is_optional=optional,
            )
        except Exception as e:
            err_msg = str(e)
            if self._libstdcxx_pending and _is_libstdcxx_load_error(err_msg):
                return DependencyResult(
                    name=package,
                    status=DependencyStatus.PENDING,
                    message="Will verify after env reactivation (stale LD_LIBRARY_PATH)",
                    is_optional=optional,
                )
            return DependencyResult(
                name=package,
                status=DependencyStatus.ERROR,
                message=err_msg,
                is_optional=optional,
            )

    def _check_libvips(self, verbose: bool):
        """Check libvips native library availability."""
        if verbose:
            print("\n[Native Libraries]")

        try:
            import pyvips

            version = pyvips.version(0)
            minor = pyvips.version(1)
            micro = pyvips.version(2)
            full_version = f"{version}.{minor}.{micro}"

            result = DependencyResult(
                name="libvips",
                status=DependencyStatus.OK,
                version=full_version,
                required_version="8.10.0",
                message="",
                details={"cache_max": pyvips.cache_get_max()},
            )

        except Exception as e:
            err_msg = str(e)
            if self._libstdcxx_pending and _is_libstdcxx_load_error(err_msg):
                result = DependencyResult(
                    name="libvips",
                    status=DependencyStatus.PENDING,
                    message="Will verify after env reactivation (stale LD_LIBRARY_PATH)",
                )
            else:
                result = DependencyResult(
                    name="libvips",
                    status=DependencyStatus.MISSING,
                    message=err_msg,
                    details={
                        "hint": "Install libvips: conda install -c conda-forge libvips "
                        "or download from Zenodo (Windows)"
                    },
                )

        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_cuda(self, verbose: bool):
        """Check GPU packages and CUDA hardware availability."""
        if verbose:
            print("\n[GPU/CUDA]")

        # Check torch package
        torch_imported = False
        try:
            import torch

            torch_imported = True
            result = DependencyResult(
                name="torch",
                status=DependencyStatus.OK,
                version=torch.__version__,
                is_optional=True,
                details={"group": "gpu"},
            )
        except ImportError:
            result = DependencyResult(
                name="torch",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="Not installed",
                is_optional=True,
                details={"group": "gpu", "hint": "Install with: kintsugi install gpu"},
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

        # Check torchvision package
        try:
            import torchvision

            result = DependencyResult(
                name="torchvision",
                status=DependencyStatus.OK,
                version=torchvision.__version__,
                is_optional=True,
                details={"group": "gpu"},
            )
        except ImportError:
            result = DependencyResult(
                name="torchvision",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="Not installed",
                is_optional=True,
                details={"group": "gpu", "hint": "Install with: kintsugi install gpu"},
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

        # Check PyTorch CUDA hardware availability
        if torch_imported:
            try:
                import torch

                cuda_available = torch.cuda.is_available()

                if cuda_available:
                    device_count = torch.cuda.device_count()
                    cuda_version = torch.version.cuda
                    devices = [torch.cuda.get_device_name(i) for i in range(device_count)]

                    result = DependencyResult(
                        name="torch-cuda",
                        status=DependencyStatus.OK,
                        version=cuda_version,
                        message=f"{device_count} GPU(s) available",
                        is_optional=True,
                        details={"devices": devices, "device_count": device_count},
                    )
                else:
                    result = DependencyResult(
                        name="torch-cuda",
                        status=DependencyStatus.OPTIONAL_MISSING,
                        message="CUDA not available (CPU mode)",
                        is_optional=True,
                    )

            except Exception as e:
                result = DependencyResult(
                    name="torch-cuda",
                    status=DependencyStatus.ERROR,
                    message=str(e),
                    is_optional=True,
                )

            self.results.append(result)
            if verbose:
                self._print_result(result)

        # Check CuPy package
        try:
            import cupy

            result = DependencyResult(
                name="cupy",
                status=DependencyStatus.OK,
                version=cupy.__version__,
                message="GPU acceleration available",
                is_optional=True,
                details={"group": "gpu"},
            )
        except ImportError:
            result = DependencyResult(
                name="cupy",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="Not installed (CPU fallback available)",
                is_optional=True,
                details={"group": "gpu", "hint": "Install with: kintsugi install gpu"},
            )
        except Exception as e:
            result = DependencyResult(
                name="cupy",
                status=DependencyStatus.ERROR,
                message=str(e),
                is_optional=True,
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_environment_safety(self, verbose: bool):
        """Run deep validation checks for known conflict conditions."""
        if verbose:
            print("\n[Environment Safety]")

        self._check_numpy_version_constraint(verbose)
        self._check_libstdcxx(verbose)
        self._check_torch_cuda_build(verbose)
        self._check_cupy_cuda_libs(verbose)
        self._check_snakemake_tres_patch(verbose)

    def _check_numpy_version_constraint(self, verbose: bool):
        """Verify numpy < 2.0 (critical for all processing)."""
        try:
            import numpy as np
            from packaging.version import parse

            if parse(np.__version__) >= parse("2.0.0"):
                result = DependencyResult(
                    name="numpy-constraint",
                    status=DependencyStatus.ERROR,
                    version=np.__version__,
                    message=(
                        "numpy >= 2.0 detected! This WILL break processing. "
                        "Fix: pip install 'numpy>=1.24,<2.0'"
                    ),
                )
            else:
                result = DependencyResult(
                    name="numpy-constraint",
                    status=DependencyStatus.OK,
                    version=np.__version__,
                    message="numpy < 2.0 constraint satisfied",
                )
        except ImportError:
            result = DependencyResult(
                name="numpy-constraint",
                status=DependencyStatus.MISSING,
                message="numpy not installed",
            )
        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _detect_libstdcxx_pending(self) -> None:
        """Set ``self._libstdcxx_pending`` if conda's libstdc++ is newer but
        ``LD_LIBRARY_PATH`` still points at the (older) system libstdc++.

        This is a silent pre-check that does NOT emit a result record — its only
        job is to flip the flag so downstream checks can downgrade cascade-prone
        load errors to PENDING. The full ``_check_libstdcxx()`` still runs later
        in the Environment Safety section and emits the user-facing record.
        """
        import sys

        if sys.platform != "linux":
            return
        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        if not conda_prefix:
            return
        ld_path = os.environ.get("LD_LIBRARY_PATH", "")
        conda_lib = os.path.join(conda_prefix, "lib")
        if conda_lib in ld_path.split(os.pathsep):
            return  # conda lib already on LD_LIBRARY_PATH — no cascade risk
        conda_libstdcxx = Path(conda_prefix) / "lib" / "libstdc++.so.6"
        if not conda_libstdcxx.exists():
            return  # no newer libstdc++ to shadow the system one
        self._libstdcxx_pending = True

    def _check_libstdcxx(self, verbose: bool):
        """Check that conda's libstdc++ is loaded instead of the system's.

        On HPC systems (e.g., HiPerGator), the system ``/lib64/libstdc++.so.6``
        may be too old (missing GLIBCXX_3.4.30, CXXABI_1.3.15). If conda's
        newer version isn't on LD_LIBRARY_PATH, packages like scikit-learn,
        matplotlib, pyvips, and stackview will fail to import despite being
        installed.
        """
        import sys

        if sys.platform != "linux":
            return  # Only relevant on Linux (HPC)

        conda_prefix = os.environ.get("CONDA_PREFIX", "")
        if not conda_prefix:
            return  # Not in a conda env

        # Check if LD_LIBRARY_PATH includes conda lib
        ld_path = os.environ.get("LD_LIBRARY_PATH", "")
        conda_lib = os.path.join(conda_prefix, "lib")
        conda_lib_in_path = conda_lib in ld_path.split(os.pathsep)

        # Check if conda's libstdc++ exists and is newer than system's
        conda_libstdcxx = Path(conda_prefix) / "lib" / "libstdc++.so.6"
        if not conda_libstdcxx.exists():
            return  # No conda libstdc++ to check against

        if conda_lib_in_path:
            result = DependencyResult(
                name="libstdcxx",
                status=DependencyStatus.OK,
                version=conda_lib,
                message="Conda lib directory on LD_LIBRARY_PATH",
            )
        else:
            # Check if this is actually causing problems by testing a
            # known-affected import
            try:
                import sklearn  # noqa: F401

                # It works — system lib might be new enough
                result = DependencyResult(
                    name="libstdcxx",
                    status=DependencyStatus.OK,
                    message="System libstdc++ is sufficient",
                )
            except (ImportError, OSError):
                # Try to auto-deploy activation scripts before reporting error
                deployed = deploy_activation_scripts()
                msg = f"Conda lib dir not on LD_LIBRARY_PATH ({conda_lib}). "
                if deployed:
                    msg += (
                        "Activation scripts deployed. "
                        "Run: conda deactivate && conda activate KINTSUGI"
                    )
                    # Transient: will clear on next activate. Report PENDING so
                    # the install flow doesn't exit non-zero on a self-healing
                    # condition, and so the cascade of CXXABI/GLIBCXX import
                    # failures below gets downgraded to PENDING too.
                    status = DependencyStatus.PENDING
                else:
                    msg += (
                        "System libstdc++ may be too old, causing import failures "
                        "for scikit-learn, matplotlib, pyvips, stackview. "
                        "Fix: kintsugi fix-hpc"
                    )
                    status = DependencyStatus.ERROR
                result = DependencyResult(
                    name="libstdcxx",
                    status=status,
                    message=msg,
                )
        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_torch_cuda_build(self, verbose: bool):
        """Verify torch is CUDA build, not CPU-only."""
        try:
            import torch

            if hasattr(torch.version, "cuda") and torch.version.cuda is not None:
                result = DependencyResult(
                    name="torch-build",
                    status=DependencyStatus.OK,
                    version=f"CUDA {torch.version.cuda}",
                    message="CUDA-enabled PyTorch build",
                    is_optional=True,
                    details={"group": "gpu"},
                )
            else:
                result = DependencyResult(
                    name="torch-build",
                    status=DependencyStatus.ERROR,
                    version=torch.__version__,
                    message=(
                        "CPU-only PyTorch detected! GPU processing will fail. "
                        "Fix: kintsugi install gpu"
                    ),
                    is_optional=True,
                    details={"group": "gpu"},
                )
        except ImportError:
            return  # torch not installed = fine (optional)
        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_cupy_cuda_libs(self, verbose: bool):
        """Verify CuPy can actually load CUDA libraries (not just installed)."""
        try:
            import cupy  # noqa: F401
        except ImportError:
            return  # cupy not installed = fine
        except Exception as e:
            err_str = str(e)
            if "libcufft" in err_str or "libcublas" in err_str or "libcusolver" in err_str:
                result = DependencyResult(
                    name="cupy-libs",
                    status=DependencyStatus.ERROR,
                    message=(
                        f"CuPy installed but CUDA libraries missing: {err_str}. "
                        "Fix: conda install cuda-libraries cuda-cudart-dev -c nvidia"
                    ),
                    is_optional=True,
                    details={"group": "gpu"},
                )
                self.results.append(result)
                if verbose:
                    self._print_result(result)

    def _check_snakemake_tres_patch(self, verbose: bool):
        """Check if the SLURM_TRES_PER_TASK patch is applied to snakemake jobstep plugin."""
        try:
            import snakemake_executor_plugin_slurm_jobstep

            init_file = Path(snakemake_executor_plugin_slurm_jobstep.__file__)
            source = init_file.read_text()
            if "SLURM_TRES_PER_TASK" in source:
                result = DependencyResult(
                    name="slurm-tres-patch",
                    status=DependencyStatus.OK,
                    message="SLURM_TRES_PER_TASK patch applied",
                    is_optional=True,
                    details={"group": "workflow"},
                )
            else:
                result = DependencyResult(
                    name="slurm-tres-patch",
                    status=DependencyStatus.ERROR,
                    message=(
                        "SLURM_TRES_PER_TASK patch NOT applied! GPU jobs will fail. "
                        "Fix: kintsugi patch slurm"
                    ),
                    is_optional=True,
                    details={"group": "workflow"},
                )
        except ImportError:
            return  # snakemake not installed = fine
        except Exception:
            return  # can't read file = skip
        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _print_result(self, result: DependencyResult):
        """Print a single dependency check result."""
        symbols = {
            DependencyStatus.OK: "[OK]",
            DependencyStatus.MISSING: "[MISSING]",
            DependencyStatus.VERSION_MISMATCH: "[VERSION]",
            DependencyStatus.ERROR: "[ERROR]",
            DependencyStatus.OPTIONAL_MISSING: "[SKIP]",
            DependencyStatus.PENDING: "[PENDING]",
        }

        symbol = symbols.get(result.status, "[?]")
        version_str = f" v{result.version}" if result.version else ""
        group = result.details.get("group", "")
        group_str = f" ({group})" if group else ""
        # Show (optional) only for items without a group tag (group already implies optional)
        optional_str = " (optional)" if result.is_optional and not group else ""

        print(f"  {symbol:10} {result.name}{version_str}{group_str}{optional_str}")

        if result.message and result.status not in (
            DependencyStatus.OK,
            DependencyStatus.OPTIONAL_MISSING,
        ):
            print(f"             {result.message}")

    def _generate_summary(self, verbose: bool) -> dict:
        """Generate summary of all checks."""
        required_ok = sum(
            1 for r in self.results if not r.is_optional and r.status == DependencyStatus.OK
        )
        required_total = sum(1 for r in self.results if not r.is_optional)
        # PENDING means "will verify after env reactivation" — don't count as
        # missing (the env isn't broken, it just needs activate).
        required_missing = [
            r.name
            for r in self.results
            if not r.is_optional
            and r.status not in (DependencyStatus.OK, DependencyStatus.PENDING)
        ]
        required_pending = [
            r.name
            for r in self.results
            if not r.is_optional and r.status == DependencyStatus.PENDING
        ]

        optional_ok = sum(
            1 for r in self.results if r.is_optional and r.status == DependencyStatus.OK
        )
        optional_total = sum(1 for r in self.results if r.is_optional)

        # Group optional packages by group
        groups_status = {}
        for r in self.results:
            if r.is_optional and "group" in r.details:
                group = r.details["group"]
                if group not in groups_status:
                    groups_status[group] = {"ok": 0, "total": 0}
                groups_status[group]["total"] += 1
                if r.status == DependencyStatus.OK:
                    groups_status[group]["ok"] += 1

        summary = {
            "required": {
                "passed": required_ok,
                "total": required_total,
                "missing": required_missing,
                "pending": required_pending,
            },
            "optional": {
                "passed": optional_ok,
                "total": optional_total,
            },
            "groups": groups_status,
            # PENDING counts as OK for exit-status purposes — the env isn't broken.
            "all_required_ok": (required_ok + len(required_pending)) == required_total,
            "libstdcxx_pending": self._libstdcxx_pending,
            "results": [
                {
                    "name": r.name,
                    "status": r.status.value,
                    "version": r.version,
                    "optional": r.is_optional,
                    "group": r.details.get("group"),
                }
                for r in self.results
            ],
        }

        if verbose:
            print("\n" + "=" * 60)
            print("SUMMARY")
            print("=" * 60)
            print(f"  Required: {required_ok}/{required_total} passed")
            print(f"  Optional: {optional_ok}/{optional_total} available")

            print("\n  Optional groups:")
            for group, status in groups_status.items():
                checkmark = "✓" if status["ok"] == status["total"] else "○"
                print(f"    {checkmark} {group}: {status['ok']}/{status['total']}")

            if required_missing:
                print("\n  Missing required dependencies:")
                for name in required_missing:
                    print(f"    - {name}")
                print("\n  Install with: pip install -e .")
            elif required_pending:
                print(
                    "\n  ⚠ Environment activation scripts were just deployed."
                )
                print("     Reactivate to verify the remaining checks:")
                print("       conda deactivate && conda activate KINTSUGI")
                print("       kintsugi check")
                print("\n  Pending (will clear after reactivate):")
                for name in required_pending:
                    print(f"    - {name}")
            else:
                print("\n  ✓ All required dependencies satisfied!")

            print("\n  To install optional features:")
            for name, info in OPTIONAL_GROUPS.items():
                if name == "full":
                    continue
                print(f"    kintsugi install {name:14s} # {info['description']}")
            print(f"    {'kintsugi install all':>30s}       # All optional features")

        return summary


def check_dependencies(verbose: bool = True, strict: bool = False) -> dict:
    """
    Convenience function to check all dependencies.

    Parameters
    ----------
    verbose : bool
        Print status messages.
    strict : bool
        If True, returns non-zero exit info when errors are found.

    Returns
    -------
    dict
        Summary of dependency checks.
    """
    checker = DependencyChecker()
    summary = checker.check_all(verbose=verbose)

    if strict:
        errors = [r for r in checker.results if r.status == DependencyStatus.ERROR]
        summary["has_errors"] = len(errors) > 0
        summary["errors"] = [{"name": r.name, "message": r.message} for r in errors]

    return summary


if __name__ == "__main__":
    check_dependencies()
