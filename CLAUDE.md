# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KINTSUGI (Knowledge Integration with New Technologies for Simplified User-Guided Image processing) is a multiplex immunofluorescence image processing pipeline for scientific research (Smith et al., 2025, STAR Protocols). Key capabilities include illumination correction, image stitching, deconvolution, extended depth of focus, multi-cycle registration, autofluorescence removal, segmentation, and spatial analysis.

## Claude Code Integration

KINTSUGI includes an MCP (Model Context Protocol) server that exposes image processing tools to Claude Code. This enables Claude to act as an AI image processing assistant.

### Starting the MCP Server

```bash
# Start the MCP server
kintsugi mcp start

# List available tools
kintsugi mcp tools

# Generate Claude Code configuration
kintsugi mcp config /path/to/project
```

### Configuring Claude Code

Add to `.claude/settings.local.json`:

```json
{
    "mcpServers": {
        "kintsugi": {
            "command": "kintsugi",
            "args": ["mcp", "start"],
            "cwd": "/path/to/your/project"
        }
    }
}
```

### Available MCP Tools

**Signal Isolation:**
- `load_channel` - Load channel image from project
- `subtract_blank` - Autofluorescence subtraction
- `denoise` - Basic denoising (median, Gaussian, bilateral)
- `denoise_advanced` - Advanced denoising (N2V, NLM, BM3D, adaptive)
- `apply_clahe` - Contrast enhancement
- `clean_background` - Background removal
- `gaussian_subtract` - Structured background removal

**Quality Assessment:**
- `assess_quality` - Automatic quality assessment
- `compute_snr` - Signal-to-noise ratio

**Visualization:**
- `get_image_stats` - Image statistics
- `get_thumbnail` - Downsampled preview

**Workflow:**
- `list_channels` - List available channels
- `save_processed` - Save processed output
- `suggest_parameters` - AI-powered parameter suggestions
- `generate_jupyter_cell` - Generate interactive tuning code

**Parameter Learning:**
- `get_learned_parameters` - Get recommendations from history
- `record_successful_parameters` - Record approved params
- `suggest_with_learning` - Combine analysis + learned history
- `approve_and_learn` - Approve and record for future learning
- `get_learning_statistics` - Database statistics

## Build and Development Commands

```bash
# Install in development mode
pip install -e ".[dev]"

# Install with Claude Code integration
pip install -e ".[claude]"

# Run tests
pytest tests/ -v
pytest tests/ --cov=src/kintsugi --cov-report=html

# Linting and formatting
ruff check src/ tests/
black src/ tests/

# Type checking
mypy src/ tests/

# Build documentation (from docs/ directory)
./make.bat html  # Windows
make html        # Linux/macOS

# Build package
python -m build
```

## Architecture

**Core Package** (`src/kintsugi/`):
- `cli.py` - Click-based CLI with subcommands: check, register, template, info, mcp, init
- `kreg.py`, `kstitch.py`, `kview2.py` - Bridge modules to notebook implementations
- `deps.py` - Runtime dependency validation (Python packages, libvips, CUDA)
- `edf.py` - Extended depth of focus processing
- `project.py` - Project structure management
- `deprecation.py` - Deprecation warnings and migration guidance

**MCP Server** (`src/kintsugi/mcp/`):
- `server.py` - FastMCP server implementation
- `tools/signal_isolation.py` - Signal processing tools
- `tools/quality_assessment.py` - QC tools
- `tools/visualization.py` - Image preview tools
- `tools/workflow.py` - Workflow management tools
- `tools/learning.py` - Parameter learning with SQLite storage

**Pure Python Modules** (`src/kintsugi/`):
- `signal/` - Autofluorescence subtraction with intelligent parameter suggestion
- `denoise/` - Denoising algorithms (median, Gaussian, bilateral, NLM, N2V, CARE, BM3D-lite)
- `qc/` - Quality control (ImageQC, CellQC, MarkerQC, BatchQC)
- `segment/` - Segmentation (SAM wrapper, watershed, postprocessing)

**Notebook Modules** (`notebooks/`):
- `Kreg/` - Registration system (14 modules) wrapping VALIS for rigid/non-rigid registration
- `Kstitch/` - GPU-accelerated image stitching with CuPy
- `Kview2/` - 20+ interactive visualization functions for Jupyter
- `KDecon/` - Deconvolution utilities
- `instanseg/` - Instance segmentation models (local customized version)
- `Kutils.py` - Signal isolation utilities (legacy, used by interactive tuners)

**Processing Notebooks** (sequential workflow):
1. `1_Single_Channel_Eval.ipynb` - Parameter tuning and setup
2. `2_Cycle_Processing.ipynb` - Batch cycle processing
3. `3_Signal_Isolation_QC.ipynb` - Signal isolation with integrated QC (Claude-guided or interactive)
4. `4_Segmentation_Analysis.ipynb` - Segmentation & spatial analysis

## Key Patterns

- **Lazy loading**: `__init__.py` uses lazy imports to avoid loading heavy dependencies at startup
- **Bridge pattern**: `src/kintsugi/*.py` modules wrap notebook implementations for package use
- **Configuration-driven**: Registration and processing pipelines use JSON config files (see `notebooks/config_example.json`)
- **Platform-specific**: Windows/Linux/macOS have different environment files (`envs/`) and dependency handling
- **Parameter learning**: Successful parameters are stored in SQLite databases, indexed by tissue type and marker name
- **MCP integration**: Claude Code can access image processing tools via the MCP server

## Dependencies

Core: numpy<2.0, scipy, pandas, scikit-image, opencv-contrib-python-headless, pyvips, valis-wsi

Optional groups in pyproject.toml:
- `[gpu]` - PyTorch + CuPy for GPU acceleration
- `[claude]` - MCP server for Claude Code integration
- `[denoise]` - PyTorch for N2V/CARE denoising
- `[viz]` - Napari visualization
- `[analysis]` - scanpy, scimap for spatial analysis
- `[bio]` - Bio formats I/O (OME-TIFF, LIF, etc.)
- `[dl]` - Deep learning segmentation (InstanSeg)
- `[full]` - All optional dependencies
- `[java]` - (DEPRECATED) JPype + PyImageJ for BioFormats

**External requirements**: libvips (native library). Java/Maven no longer required.

## Testing

Tests are in `tests/` with fixtures in `conftest.py`. Key fixtures: `sample_image`, `sample_multichannel_image`, `sample_stack`, `sample_tiff`, `sample_config`, `temp_dir`.

CI runs on Windows/Linux/macOS with Python 3.10-3.12.

## Code Style

- Line length: 100 characters
- Formatter: Black
- Linter: Ruff (with isort, pyupgrade, flake8-bugbear)
- Python 3.10+ required

## Release and Versioning

This project uses **Conventional Commits** and **automatic semantic versioning**.

### Commit Message Format

All commits should follow the [Conventional Commits](https://www.conventionalcommits.org/) format:

```
<type>(<scope>): <description>

[optional body]

[optional footer(s)]
```

**Types:**
- `feat`: New feature (bumps MINOR version)
- `fix`: Bug fix (bumps PATCH version)
- `docs`: Documentation changes
- `style`: Code style changes (formatting, etc.)
- `refactor`: Code refactoring
- `perf`: Performance improvements
- `test`: Adding/updating tests
- `build`: Build system changes
- `ci`: CI/CD changes
- `chore`: Maintenance tasks

**Breaking Changes:** Add `!` after type or include `BREAKING CHANGE:` in footer (bumps MAJOR version).

### Examples

```bash
feat(mcp): add new image processing tool
fix(segmentation): resolve edge detection issue
docs: update installation instructions
refactor(signal)!: change API for background subtraction
```

### Pre-commit Hooks

Install hooks for automatic commit validation:

```bash
pre-commit install
pre-commit install --hook-type commit-msg
```

### Creating Releases

Releases are automated via GitHub Actions:

1. **Automatic releases**: Push to `main` triggers version analysis and release creation
2. **Manual releases**: Use workflow dispatch in GitHub Actions
3. **Local release script**: `python scripts/release.py --auto`

```bash
# Preview release (dry run)
python scripts/release.py --dry-run --auto

# Manual bump
python scripts/release.py --bump minor
```

### Changelog

The `CHANGELOG.md` is automatically updated during releases based on commit messages.

## Migration Notes

The old Notebooks 3 (Signal Isolation) and 5 (DL Channel Refinement) have been replaced with a unified `3_Signal_Isolation_QC.ipynb` that supports both Claude-guided and interactive workflows. See `notebooks/MIGRATION_GUIDE.md` for detailed transition guidance.
