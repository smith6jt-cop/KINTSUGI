"""
Tests for the main kintsugi package.
"""

import pytest


class TestPackageImport:
    """Test package import and basic functionality."""

    def test_import_kintsugi(self):
        """Test that kintsugi package can be imported."""
        import kintsugi

        assert hasattr(kintsugi, "__version__")
        assert isinstance(kintsugi.__version__, str)

    def test_version_format(self):
        """Test that version follows semantic versioning."""
        import kintsugi

        version = kintsugi.__version__
        parts = version.split(".")
        assert len(parts) >= 2, f"Version should have at least major.minor: {version}"

        # Check that major and minor are numeric
        assert parts[0].isdigit(), f"Major version should be numeric: {parts[0]}"
        assert parts[1].isdigit(), f"Minor version should be numeric: {parts[1]}"

    def test_author_info(self):
        """Test that author information is present."""
        import kintsugi

        assert hasattr(kintsugi, "__author__")
        assert hasattr(kintsugi, "__email__")

    def test_check_dependencies_function(self):
        """Test that check_dependencies function exists."""
        import kintsugi

        assert hasattr(kintsugi, "check_dependencies")
        assert callable(kintsugi.check_dependencies)

    def test_get_config_template(self):
        """Test that config template function works."""
        import kintsugi

        assert hasattr(kintsugi, "get_config_template")

        config = kintsugi.get_config_template()
        assert isinstance(config, dict)

        # Check required keys
        assert "src_dir" in config
        assert "dst_dir" in config
        assert "image_type" in config


class TestLazyImports:
    """Test lazy import functionality."""

    def test_lazy_import_deps(self):
        """Test lazy import of deps module."""
        import kintsugi

        deps = kintsugi.deps
        assert hasattr(deps, "DependencyChecker")
        assert hasattr(deps, "check_dependencies")

    def test_invalid_attribute_raises(self):
        """Test that invalid attributes raise AttributeError."""
        import kintsugi

        with pytest.raises(AttributeError):
            _ = kintsugi.nonexistent_module


class TestAllExports:
    """Test __all__ exports."""

    def test_all_defined(self):
        """Test that __all__ is defined."""
        import kintsugi

        assert hasattr(kintsugi, "__all__")
        assert isinstance(kintsugi.__all__, list)

    def test_all_exports_exist(self):
        """Test that all items in __all__ can be accessed."""
        import kintsugi

        for name in kintsugi.__all__:
            assert hasattr(kintsugi, name), f"Missing export: {name}"
