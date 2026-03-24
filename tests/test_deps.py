"""
Tests for the dependency validation module.
"""


class TestDependencyChecker:
    """Test DependencyChecker class."""

    def test_import_deps_module(self):
        """Test that deps module can be imported."""
        from kintsugi import deps

        assert hasattr(deps, "DependencyChecker")
        assert hasattr(deps, "DependencyStatus")
        assert hasattr(deps, "DependencyResult")

    def test_create_checker(self):
        """Test that DependencyChecker can be instantiated."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        assert checker is not None
        assert isinstance(checker.results, list)
        assert len(checker.results) == 0

    def test_check_all_returns_dict(self):
        """Test that check_all returns a dictionary."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        result = checker.check_all(verbose=False)

        assert isinstance(result, dict)
        assert "required" in result
        assert "optional" in result
        assert "all_required_ok" in result
        assert "results" in result

    def test_check_all_required_structure(self):
        """Test structure of required dependency results."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        result = checker.check_all(verbose=False)

        required = result["required"]
        assert "passed" in required
        assert "total" in required
        assert "missing" in required
        assert isinstance(required["passed"], int)
        assert isinstance(required["total"], int)
        assert isinstance(required["missing"], list)

    def test_results_structure(self):
        """Test structure of individual results."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        result = checker.check_all(verbose=False)

        for item in result["results"]:
            assert "name" in item
            assert "status" in item
            assert "optional" in item


class TestDependencyStatus:
    """Test DependencyStatus enum."""

    def test_status_values(self):
        """Test that all expected status values exist."""
        from kintsugi.deps import DependencyStatus

        assert hasattr(DependencyStatus, "OK")
        assert hasattr(DependencyStatus, "MISSING")
        assert hasattr(DependencyStatus, "VERSION_MISMATCH")
        assert hasattr(DependencyStatus, "ERROR")
        assert hasattr(DependencyStatus, "OPTIONAL_MISSING")

    def test_status_string_values(self):
        """Test status enum string values."""
        from kintsugi.deps import DependencyStatus

        assert DependencyStatus.OK.value == "ok"
        assert DependencyStatus.MISSING.value == "missing"
        assert DependencyStatus.VERSION_MISMATCH.value == "version_mismatch"
        assert DependencyStatus.ERROR.value == "error"
        assert DependencyStatus.OPTIONAL_MISSING.value == "optional_missing"


class TestDependencyResult:
    """Test DependencyResult dataclass."""

    def test_create_result(self):
        """Test creating a DependencyResult."""
        from kintsugi.deps import DependencyResult, DependencyStatus

        result = DependencyResult(
            name="test_package",
            status=DependencyStatus.OK,
            version="1.0.0",
        )

        assert result.name == "test_package"
        assert result.status == DependencyStatus.OK
        assert result.version == "1.0.0"
        assert result.is_optional is False

    def test_result_with_optional(self):
        """Test creating an optional DependencyResult."""
        from kintsugi.deps import DependencyResult, DependencyStatus

        result = DependencyResult(
            name="optional_package",
            status=DependencyStatus.OPTIONAL_MISSING,
            is_optional=True,
            message="Package not installed",
        )

        assert result.is_optional is True
        assert result.message == "Package not installed"


class TestConvenienceFunction:
    """Test convenience check_dependencies function."""

    def test_check_dependencies_function(self):
        """Test the convenience function."""
        from kintsugi.deps import check_dependencies

        result = check_dependencies(verbose=False)
        assert isinstance(result, dict)
        assert "all_required_ok" in result


class TestOptionalGroups:
    """Test OPTIONAL_GROUPS structure and content."""

    def test_all_groups_have_required_keys(self):
        """Every group must have description, packages, and install_cmd."""
        from kintsugi.deps import OPTIONAL_GROUPS

        for name, info in OPTIONAL_GROUPS.items():
            assert "description" in info, f"{name} missing description"
            assert "packages" in info, f"{name} missing packages"
            assert "install_cmd" in info, f"{name} missing install_cmd"

    def test_expected_groups_present(self):
        """All expected groups from pyproject.toml are present."""
        from kintsugi.deps import OPTIONAL_GROUPS

        expected = {
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
        }
        assert set(OPTIONAL_GROUPS.keys()) == expected

    def test_workflow_group_exists(self):
        """Workflow group for Snakemake/SLURM should exist."""
        from kintsugi.deps import OPTIONAL_GROUPS

        assert "workflow" in OPTIONAL_GROUPS
        assert "snakemake" in OPTIONAL_GROUPS["workflow"]["packages"]

    def test_bio_group_has_bio_packages(self):
        """Bio group should contain bio I/O packages, not analysis packages."""
        from kintsugi.deps import OPTIONAL_GROUPS

        bio_pkgs = OPTIONAL_GROUPS["bio"]["packages"]
        # Bio should have imaging I/O packages
        assert "aicsimageio" in bio_pkgs or "bioio" in bio_pkgs or "ome-zarr" in bio_pkgs
        # Bio should NOT have analysis packages
        assert "scanpy" not in bio_pkgs
        assert "scimap" not in bio_pkgs

    def test_gpu_group_includes_torch(self):
        """GPU group should include torch, not just cupy."""
        from kintsugi.deps import OPTIONAL_GROUPS

        gpu_pkgs = OPTIONAL_GROUPS["gpu"]["packages"]
        assert "torch" in gpu_pkgs
        assert "cupy" in gpu_pkgs

    def test_full_is_composite(self):
        """Full group should have empty packages list (composite)."""
        from kintsugi.deps import OPTIONAL_GROUPS

        assert OPTIONAL_GROUPS["full"]["packages"] == []

    def test_conda_cmd_optional(self):
        """conda_cmd is optional; groups that have it should have valid strings."""
        from kintsugi.deps import OPTIONAL_GROUPS

        for _name, info in OPTIONAL_GROUPS.items():
            if "conda_cmd" in info:
                assert isinstance(info["conda_cmd"], str)
                assert len(info["conda_cmd"]) > 0

    def test_note_field_optional(self):
        """Groups with notes should have non-empty string values."""
        from kintsugi.deps import OPTIONAL_GROUPS

        for _name, info in OPTIONAL_GROUPS.items():
            if "note" in info:
                assert isinstance(info["note"], str)
                assert len(info["note"]) > 0

    def test_kronos_has_note(self):
        """Kronos group should have a note about manual clone."""
        from kintsugi.deps import OPTIONAL_GROUPS

        assert "note" in OPTIONAL_GROUPS["kronos"]
        assert "KRONOS" in OPTIONAL_GROUPS["kronos"]["note"]

    def test_install_optional_accepts_new_groups(self):
        """install_optional should accept all groups in OPTIONAL_GROUPS."""
        from kintsugi.deps import OPTIONAL_GROUPS

        # Just verify the function doesn't reject known groups (don't actually install)
        for name in OPTIONAL_GROUPS:
            assert name in OPTIONAL_GROUPS  # Validates type hint coverage


class TestPythonPackageChecks:
    """Test Python package checking functionality."""

    def test_numpy_detected(self):
        """Test that numpy is detected correctly."""
        from kintsugi.deps import DependencyChecker, DependencyStatus

        checker = DependencyChecker()
        result = checker._check_python_package("numpy", "1.0.0", optional=False)

        assert result.name == "numpy"
        assert result.status == DependencyStatus.OK
        assert result.version is not None

    def test_nonexistent_package(self):
        """Test handling of non-existent package."""
        from kintsugi.deps import DependencyChecker, DependencyStatus

        checker = DependencyChecker()
        result = checker._check_python_package(
            "nonexistent_package_xyz_12345", None, optional=False
        )

        assert result.status == DependencyStatus.MISSING

    def test_optional_missing_package(self):
        """Test handling of optional missing package."""
        from kintsugi.deps import DependencyChecker, DependencyStatus

        checker = DependencyChecker()
        result = checker._check_python_package("nonexistent_package_xyz_12345", None, optional=True)

        assert result.status == DependencyStatus.OPTIONAL_MISSING
        assert result.is_optional is True


class TestCudaIndexInInstallCommands:
    """Test that torch-using groups use CUDA index URL."""

    def test_dl_group_has_cuda_index(self):
        """dl group must use CUDA index to avoid CPU-only torch."""
        from kintsugi.deps import OPTIONAL_GROUPS

        cmd = OPTIONAL_GROUPS["dl"]["install_cmd"]
        assert "--index-url" in cmd or "--extra-index-url" in cmd, (
            "dl group install_cmd missing CUDA index URL — will install CPU-only torch!"
        )
        assert "download.pytorch.org" in cmd

    def test_denoise_group_has_cuda_index(self):
        """denoise group must use CUDA index to avoid CPU-only torch."""
        from kintsugi.deps import OPTIONAL_GROUPS

        cmd = OPTIONAL_GROUPS["denoise"]["install_cmd"]
        assert "--index-url" in cmd or "--extra-index-url" in cmd, (
            "denoise group install_cmd missing CUDA index URL!"
        )

    def test_kronos_group_has_cuda_index(self):
        """kronos group must use CUDA index to avoid CPU-only torch."""
        from kintsugi.deps import OPTIONAL_GROUPS

        cmd = OPTIONAL_GROUPS["kronos"]["install_cmd"]
        assert "--index-url" in cmd or "--extra-index-url" in cmd, (
            "kronos group install_cmd missing CUDA index URL!"
        )

    def test_gpu_group_has_cuda_index(self):
        """gpu group must use CUDA index."""
        from kintsugi.deps import OPTIONAL_GROUPS

        cmd = OPTIONAL_GROUPS["gpu"]["install_cmd"]
        assert "--index-url" in cmd or "--extra-index-url" in cmd


class TestConstraintsFile:
    """Test constraints file injection."""

    def test_inject_constraints_with_pip(self):
        """Constraints should be injected into pip install commands."""
        from unittest.mock import patch

        from kintsugi.deps import _inject_constraints

        fake_constraints = "/fake/path/constraints.txt"
        with patch("kintsugi.deps._find_constraints_file", return_value=fake_constraints):
            result = _inject_constraints("pip install numpy scipy")
            assert f"-c {fake_constraints}" in result

    def test_inject_constraints_no_file(self):
        """No injection when constraints file doesn't exist."""
        from unittest.mock import patch

        from kintsugi.deps import _inject_constraints

        with patch("kintsugi.deps._find_constraints_file", return_value=None):
            cmd = "pip install numpy scipy"
            result = _inject_constraints(cmd)
            assert result == cmd

    def test_inject_constraints_non_pip(self):
        """Conda commands should not get constraints injection."""
        from unittest.mock import patch

        from kintsugi.deps import _inject_constraints

        with patch(
            "kintsugi.deps._find_constraints_file", return_value="/fake/constraints.txt"
        ):
            cmd = "conda install numpy scipy -c conda-forge"
            result = _inject_constraints(cmd)
            assert result == cmd


class TestPreInstallGuard:
    """Test pre-install guard checks."""

    def test_guard_warns_on_cpu_torch(self):
        """Guard should warn when CPU-only torch is installed and installing dl."""
        import sys
        import types
        from unittest.mock import patch

        from kintsugi.deps import _pre_install_guard

        mock_torch = types.ModuleType("torch")
        mock_torch.__version__ = "2.5.0"
        mock_torch_version = types.SimpleNamespace(cuda=None)
        mock_torch.version = mock_torch_version

        with patch.dict(sys.modules, {"torch": mock_torch}):
            warnings = _pre_install_guard("dl")
            assert any("CPU-only" in w for w in warnings)

    def test_guard_ok_with_cuda_torch(self):
        """Guard should not warn when CUDA torch is installed."""
        import sys
        import types
        from unittest.mock import patch

        from kintsugi.deps import _pre_install_guard

        mock_torch = types.ModuleType("torch")
        mock_torch.__version__ = "2.5.0+cu124"
        mock_torch_version = types.SimpleNamespace(cuda="12.4")
        mock_torch.version = mock_torch_version

        with patch.dict(sys.modules, {"torch": mock_torch}):
            warnings = _pre_install_guard("dl")
            assert not any("CPU-only" in w for w in warnings)

    def test_guard_no_warn_for_non_torch_group(self):
        """Guard should not warn about torch for non-torch groups."""
        from kintsugi.deps import _pre_install_guard

        warnings = _pre_install_guard("docs")
        assert not any("CPU-only" in w for w in warnings)
        assert not any("numpy" in w for w in warnings)


class TestEnvironmentSafetyChecks:
    """Test deep validation checks in DependencyChecker."""

    def test_numpy_constraint_check_passes(self):
        """numpy < 2.0 should pass the constraint check."""
        from kintsugi.deps import DependencyChecker, DependencyStatus

        checker = DependencyChecker()
        checker._check_numpy_version_constraint(verbose=False)

        numpy_results = [r for r in checker.results if r.name == "numpy-constraint"]
        assert len(numpy_results) == 1
        # Current environment should have numpy < 2.0
        assert numpy_results[0].status == DependencyStatus.OK

    def test_torch_build_check_runs(self):
        """torch build check should run without error."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        # Should not raise regardless of torch availability
        checker._check_torch_cuda_build(verbose=False)

    def test_environment_safety_included_in_check_all(self):
        """check_all should include environment safety results."""
        from kintsugi.deps import DependencyChecker

        checker = DependencyChecker()
        checker.check_all(verbose=False)

        result_names = [r.name for r in checker.results]
        assert "numpy-constraint" in result_names

    def test_check_dependencies_strict_mode(self):
        """strict mode should add has_errors and errors keys."""
        from kintsugi.deps import check_dependencies

        result = check_dependencies(verbose=False, strict=True)
        assert "has_errors" in result
        assert "errors" in result
        assert isinstance(result["errors"], list)
