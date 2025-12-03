"""
Dependency validation module for KINTSUGI.

Provides runtime checking of all external dependencies including
native libraries (libvips), Java/Maven, and Python packages.
"""

import os
import shutil
import subprocess
from dataclasses import dataclass, field
from enum import Enum


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


class DependencyChecker:
    """
    Validates all KINTSUGI dependencies at runtime.

    Checks:
    - Python packages (core and optional)
    - Native libraries (libvips)
    - Java runtime and Maven
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

        # Java/Maven
        self._check_java(verbose)

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
            ("opencv-python", None),  # opencv-contrib-python-headless
            ("tifffile", "2023.7.0"),
            ("valis-wsi", "1.2.0"),
            ("pyvips", "2.2.0"),
            ("stackview", "0.18.0"),
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
        """Check optional Python package dependencies."""
        if verbose:
            print("\n[Optional Packages]")

        optional_packages = [
            # GPU
            ("torch", "2.0.0", "gpu"),
            ("torchvision", "0.15.0", "gpu"),
            ("cupy", None, "gpu"),
            # Java
            ("jpype1", "1.5.0", "java"),
            ("scyjava", "1.0.0", "java"),
            ("pyimagej", "1.4.0", "java"),
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
            # Bio formats
            ("aicsimageio", None, "bio"),
            ("ome-types", "0.5.0", "bio"),
            ("slideio", "2.6.0", "bio"),
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
        # Handle package name variations
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
                # Try importlib.metadata
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

            # Try to get version (this validates the native library)
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
                    "or download from Zenodo"
                },
            )

        self.results.append(result)
        if verbose:
            self._print_result(result)

    def _check_java(self, verbose: bool):
        """Check Java runtime and Maven availability."""
        if verbose:
            print("\n[Java/Maven]")

        # Check Java
        java_result = self._check_java_runtime()
        self.results.append(java_result)
        if verbose:
            self._print_result(java_result)

        # Check Maven
        maven_result = self._check_maven()
        self.results.append(maven_result)
        if verbose:
            self._print_result(maven_result)

        # Check JPype JVM
        jpype_result = self._check_jpype()
        self.results.append(jpype_result)
        if verbose:
            self._print_result(jpype_result)

    def _check_java_runtime(self) -> DependencyResult:
        """Check if Java runtime is available."""
        try:
            # Check JAVA_HOME
            java_home = os.environ.get("JAVA_HOME", "")

            # Try to run java -version
            java_cmd = shutil.which("java")
            if java_cmd:
                result = subprocess.run(
                    [java_cmd, "-version"], capture_output=True, text=True, timeout=10
                )
                # Java version info goes to stderr
                output = result.stderr or result.stdout
                version_line = output.split("\n")[0] if output else ""

                return DependencyResult(
                    name="java",
                    status=DependencyStatus.OK,
                    version=version_line,
                    message="",
                    is_optional=True,
                    details={"java_home": java_home, "path": java_cmd},
                )

            return DependencyResult(
                name="java",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="Java not found in PATH",
                is_optional=True,
                details={"java_home": java_home},
            )

        except Exception as e:
            return DependencyResult(
                name="java",
                status=DependencyStatus.ERROR,
                message=str(e),
                is_optional=True,
            )

    def _check_maven(self) -> DependencyResult:
        """Check if Maven is available."""
        try:
            mvn_cmd = shutil.which("mvn")
            if mvn_cmd:
                result = subprocess.run(
                    [mvn_cmd, "--version"], capture_output=True, text=True, timeout=10
                )
                output = result.stdout or ""
                version_line = output.split("\n")[0] if output else ""

                return DependencyResult(
                    name="maven",
                    status=DependencyStatus.OK,
                    version=version_line,
                    message="",
                    is_optional=True,
                    details={"path": mvn_cmd},
                )

            return DependencyResult(
                name="maven",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="Maven not found in PATH",
                is_optional=True,
            )

        except Exception as e:
            return DependencyResult(
                name="maven",
                status=DependencyStatus.ERROR,
                message=str(e),
                is_optional=True,
            )

    def _check_jpype(self) -> DependencyResult:
        """Check JPype JVM integration."""
        try:
            import jpype

            if jpype.isJVMStarted():
                jvm_path = jpype.getDefaultJVMPath()
                return DependencyResult(
                    name="jpype-jvm",
                    status=DependencyStatus.OK,
                    version=jpype.__version__,
                    message="JVM running",
                    is_optional=True,
                    details={"jvm_path": jvm_path, "jvm_started": True},
                )
            else:
                # Try to get default JVM path without starting
                try:
                    jvm_path = jpype.getDefaultJVMPath()
                    return DependencyResult(
                        name="jpype-jvm",
                        status=DependencyStatus.OK,
                        version=jpype.__version__,
                        message="JVM available (not started)",
                        is_optional=True,
                        details={"jvm_path": jvm_path, "jvm_started": False},
                    )
                except Exception as e:
                    return DependencyResult(
                        name="jpype-jvm",
                        status=DependencyStatus.OPTIONAL_MISSING,
                        version=jpype.__version__,
                        message=f"JVM not found: {e}",
                        is_optional=True,
                    )

        except ImportError:
            return DependencyResult(
                name="jpype-jvm",
                status=DependencyStatus.OPTIONAL_MISSING,
                message="JPype not installed",
                is_optional=True,
            )
        except Exception as e:
            return DependencyResult(
                name="jpype-jvm",
                status=DependencyStatus.ERROR,
                message=str(e),
                is_optional=True,
            )

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

    def _print_result(self, result: DependencyResult):
        """Print a single dependency check result."""
        # Status symbols
        symbols = {
            DependencyStatus.OK: "[OK]",
            DependencyStatus.MISSING: "[MISSING]",
            DependencyStatus.VERSION_MISMATCH: "[VERSION]",
            DependencyStatus.ERROR: "[ERROR]",
            DependencyStatus.OPTIONAL_MISSING: "[SKIP]",
        }

        symbol = symbols.get(result.status, "[?]")
        version_str = f" v{result.version}" if result.version else ""
        optional_str = " (optional)" if result.is_optional else ""

        print(f"  {symbol:10} {result.name}{version_str}{optional_str}")

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
            "all_required_ok": required_ok == required_total,
            "results": [
                {
                    "name": r.name,
                    "status": r.status.value,
                    "version": r.version,
                    "optional": r.is_optional,
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

            if required_missing:
                print("\n  Missing required dependencies:")
                for name in required_missing:
                    print(f"    - {name}")
                print("\n  Install with: pip install kintsugi[full]")
            else:
                print("\n  All required dependencies satisfied!")

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
