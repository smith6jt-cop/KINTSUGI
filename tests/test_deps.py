"""
Tests for the dependency validation module.
"""

import pytest


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
            "nonexistent_package_xyz_12345",
            None,
            optional=False
        )

        assert result.status == DependencyStatus.MISSING

    def test_optional_missing_package(self):
        """Test handling of optional missing package."""
        from kintsugi.deps import DependencyChecker, DependencyStatus

        checker = DependencyChecker()
        result = checker._check_python_package(
            "nonexistent_package_xyz_12345",
            None,
            optional=True
        )

        assert result.status == DependencyStatus.OPTIONAL_MISSING
        assert result.is_optional is True
