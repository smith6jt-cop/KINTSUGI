"""
Tests for the CLI module.
"""

import json
import pytest
from click.testing import CliRunner


class TestCLIMain:
    """Test main CLI entry point."""

    @pytest.fixture
    def runner(self):
        """Create a CLI test runner."""
        return CliRunner()

    def test_cli_help(self, runner):
        """Test that --help works."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["--help"])
        assert result.exit_code == 0
        assert "KINTSUGI" in result.output

    def test_cli_version(self, runner):
        """Test that --version works."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["--version"])
        assert result.exit_code == 0


class TestCLICheck:
    """Test the check subcommand."""

    @pytest.fixture
    def runner(self):
        """Create a CLI test runner."""
        return CliRunner()

    def test_check_command(self, runner):
        """Test that check command runs."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["check"])
        # Should complete (exit code 0) even if some deps missing
        assert result.exit_code == 0

    def test_check_quiet(self, runner):
        """Test check command with quiet flag."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["check", "--quiet"])
        assert result.exit_code == 0

    def test_check_json_output(self, runner):
        """Test check command with JSON output."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["check", "--json"])
        assert result.exit_code == 0

        # Output should be valid JSON
        data = json.loads(result.output)
        assert "required" in data
        assert "optional" in data


class TestCLITemplate:
    """Test the template subcommand."""

    @pytest.fixture
    def runner(self):
        """Create a CLI test runner."""
        return CliRunner()

    def test_template_stdout(self, runner):
        """Test template command outputs to stdout."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["template"])
        assert result.exit_code == 0

        # Output should be valid JSON
        data = json.loads(result.output)
        assert "src_dir" in data
        assert "dst_dir" in data

    def test_template_to_file(self, runner, tmp_path):
        """Test template command outputs to file."""
        from kintsugi.cli import main

        output_file = tmp_path / "config.json"
        result = runner.invoke(main, ["template", "-o", str(output_file)])
        assert result.exit_code == 0

        # File should exist and contain valid JSON
        assert output_file.exists()
        with open(output_file) as f:
            data = json.load(f)
        assert "src_dir" in data


class TestCLIInfo:
    """Test the info subcommand."""

    @pytest.fixture
    def runner(self):
        """Create a CLI test runner."""
        return CliRunner()

    def test_info_command(self, runner):
        """Test that info command runs."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["info"])
        assert result.exit_code == 0
        assert "Version" in result.output or "KINTSUGI" in result.output


class TestCLIRegister:
    """Test the register subcommand."""

    @pytest.fixture
    def runner(self):
        """Create a CLI test runner."""
        return CliRunner()

    def test_register_help(self, runner):
        """Test register command help."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["register", "--help"])
        assert result.exit_code == 0
        assert "CONFIG" in result.output

    def test_register_missing_config(self, runner):
        """Test register command with missing config."""
        from kintsugi.cli import main

        result = runner.invoke(main, ["register", "nonexistent.json"])
        assert result.exit_code != 0

    def test_register_dry_run(self, runner, tmp_path, sample_config):
        """Test register command with dry run."""
        from kintsugi.cli import main

        # Create a test config file
        config_file = tmp_path / "config.json"
        src_dir = tmp_path / "src"
        dst_dir = tmp_path / "dst"
        src_dir.mkdir()

        sample_config["src_dir"] = str(src_dir)
        sample_config["dst_dir"] = str(dst_dir)

        with open(config_file, "w") as f:
            json.dump(sample_config, f)

        result = runner.invoke(main, ["register", str(config_file), "--dry-run"])
        assert result.exit_code == 0
        assert "Dry Run" in result.output
