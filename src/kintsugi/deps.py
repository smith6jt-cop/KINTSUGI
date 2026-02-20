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

import subprocess
from dataclasses import dataclass, field
from enum import Enum
from typing import Literal


class DependencyStatus(Enum):
    """Status of a dependency check."""

    OK = "ok"
    MISSING = "missing"
    VERSION_MISMATCH = "version_mismatch"
    ERROR = "error"
    OPTIONAL_MISSING = "optional_missing"


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
        # Use CUDA 12.9 for broad GPU support including Blackwell (B200, compute 10.0)
        # cuda-cudart-dev provides headers needed for CuPy JIT compilation
        # Headers must be copied to $CONDA_PREFIX/include for CuPy to find them
        "install_cmd": "conda install cuda-libraries cuda-cudart-dev -c nvidia -y && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null; pip install torch torchvision --index-url https://download.pytorch.org/whl/cu124 && pip install cupy-cuda12x",
        "conda_cmd": "conda install pytorch torchvision pytorch-cuda=12.4 cupy cuda-libraries cuda-cudart-dev -c pytorch -c nvidia -c conda-forge && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null",
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
        "install_cmd": "pip install torch torchvision instanseg instanseg-torch kornia",
    },
    "analysis": {
        "description": "Spatial analysis (scanpy, scimap)",
        "packages": ["scanpy", "anndata", "phenograph", "scimap", "skan", "networkx"],
        "install_cmd": "pip install scanpy anndata phenograph scimap umap-learn hdbscan skan networkx",
        "conda_cmd": "conda install scanpy anndata networkx -c conda-forge && pip install scimap phenograph skan",
    },
    "bio": {
        "description": "Bio formats I/O (OME-TIFF, LIF, etc.)",
        "packages": ["aicsimageio", "bioio", "ome-zarr", "slideio"],
        "install_cmd": "pip install aicsimageio bioio bioio-ome-tiff ome-zarr slideio readlif",
    },
    "full": {
        "description": "All optional features",
        "packages": [],  # Composite group
        "install_cmd": "conda install cuda-libraries cuda-cudart-dev -c nvidia -y && cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null; pip install torch torchvision --index-url https://download.pytorch.org/whl/cu124 && pip install cupy-cuda12x napari magicgui instanseg instanseg-torch kornia scanpy anndata phenograph scimap aicsimageio bioio bioio-ome-tiff ome-zarr slideio readlif",
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


def install_optional(
    group: Literal["gpu", "viz", "dl", "analysis", "bio", "full"],
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

    info = OPTIONAL_GROUPS[group]
    print(f"\nInstalling {group}: {info['description']}")
    print("-" * 50)

    if use_conda and "conda_cmd" in info:
        cmd = info["conda_cmd"]
    else:
        cmd = info["install_cmd"]

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

        if verbose:
            print("=" * 60)
            print("KINTSUGI Dependency Check")
            print("=" * 60)

        # Core Python packages
        self._check_core_packages(verbose)

        # Native libraries
        self._check_libvips(verbose)

        # GPU/CUDA
        self._check_cuda(verbose)

        # Optional packages
        self._check_optional_packages(verbose)

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
        optional_packages = [
            # GPU
            ("torch", "2.0.0", "gpu"),
            ("torchvision", "0.15.0", "gpu"),
            ("cupy", None, "gpu"),
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
            ("aicsimageio", None, "bio"),
            ("ome-types", "0.5.0", "bio"),
            ("ome-zarr", None, "bio"),
            # Image processing
            ("valis-wsi", "1.2.0", "core"),
            ("pyvips", "2.2.0", "core"),
            ("stackview", "0.18.0", "viz"),
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
            status = DependencyStatus.OPTIONAL_MISSING if optional else DependencyStatus.MISSING
            return DependencyResult(
                name=package,
                status=status,
                required_version=min_version,
                message=str(e),
                is_optional=optional,
            )
        except Exception as e:
            return DependencyResult(
                name=package,
                status=DependencyStatus.ERROR,
                message=str(e),
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
            result = DependencyResult(
                name="libvips",
                status=DependencyStatus.MISSING,
                message=str(e),
                details={
                    "hint": "Install libvips: conda install -c conda-forge libvips "
                    "or download from Zenodo (Windows)"
                },
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_cuda(self, verbose: bool):
        """Check CUDA/GPU availability."""
        if verbose:
            print("\n[GPU/CUDA]")

        # Check PyTorch CUDA
        try:
            import torch

            cuda_available = torch.cuda.is_available()

            if cuda_available:
                device_count = torch.cuda.device_count()
                cuda_version = torch.version.cuda
                devices = [torch.cuda.get_device_name(i) for i in range(device_count)]

                result = DependencyResult(
                    name="cuda-pytorch",
                    status=DependencyStatus.OK,
                    version=cuda_version,
                    message=f"{device_count} GPU(s) available",
                    is_optional=True,
                    details={"devices": devices, "device_count": device_count},
                )
            else:
                result = DependencyResult(
                    name="cuda-pytorch",
                    status=DependencyStatus.OPTIONAL_MISSING,
                    message="CUDA not available (CPU mode)",
                    is_optional=True,
                )

        except ImportError:
            result = DependencyResult(
                name="cuda-pytorch",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="PyTorch not installed",
                is_optional=True,
                details={"hint": "Install with: kintsugi install gpu"},
            )
        except Exception as e:
            result = DependencyResult(
                name="cuda-pytorch",
                status=DependencyStatus.ERROR,
                message=str(e),
                is_optional=True,
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

        # Check CuPy
        try:
            import cupy

            result = DependencyResult(
                name="cupy",
                status=DependencyStatus.OK,
                version=cupy.__version__,
                message="GPU acceleration available",
                is_optional=True,
            )
        except ImportError:
            result = DependencyResult(
                name="cupy",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="CuPy not installed (CPU fallback available)",
                is_optional=True,
                details={"hint": "Install with: kintsugi install gpu"},
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

    def _print_result(self, result: DependencyResult):
        """Print a single dependency check result."""
        symbols = {
            DependencyStatus.OK: "[OK]",
            DependencyStatus.MISSING: "[MISSING]",
            DependencyStatus.VERSION_MISMATCH: "[VERSION]",
            DependencyStatus.ERROR: "[ERROR]",
            DependencyStatus.OPTIONAL_MISSING: "[SKIP]",
        }

        symbol = symbols.get(result.status, "[?]")
        version_str = f" v{result.version}" if result.version else ""
        group_str = f" ({result.details.get('group', '')})" if result.details.get("group") else ""
        optional_str = " (optional)" if result.is_optional else ""

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
        required_missing = [
            r.name for r in self.results if not r.is_optional and r.status != DependencyStatus.OK
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
            },
            "optional": {
                "passed": optional_ok,
                "total": optional_total,
            },
            "groups": groups_status,
            "all_required_ok": required_ok == required_total,
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
            else:
                print("\n  ✓ All required dependencies satisfied!")

            print("\n  To install optional features:")
            print("    kintsugi install gpu       # GPU acceleration (CuPy for CUDA)")
            print("    kintsugi install torch     # PyTorch for deep learning models")
            print(
                "    kintsugi install bio       # Spatial biology analysis (scanpy, scimap, squidpy)"
            )
            print("    kintsugi install viz       # Napari visualization")
            print("    kintsugi install all       # All optional features")

        return summary


def check_dependencies(verbose: bool = True) -> dict:
    """
    Convenience function to check all dependencies.

    Parameters
    ----------
    verbose : bool
        Print status messages.

    Returns
    -------
    dict
        Summary of dependency checks.
    """
    checker = DependencyChecker()
    return checker.check_all(verbose=verbose)


if __name__ == "__main__":
    check_dependencies()
