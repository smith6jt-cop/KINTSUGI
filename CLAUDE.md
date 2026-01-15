# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

KINTSUGI (Knowledge Integration with New Technologies for Simplified User-Guided Image processing) is a multiplex immunofluorescence image processing pipeline for scientific research (Smith et al., 2025, STAR Protocols). Key capabilities include illumination correction, image stitching, deconvolution, extended depth of focus, multi-cycle registration, autofluorescence removal, segmentation, and spatial analysis.

## Claude Code Integration

KINTSUGI includes an MCP (Model Context Protocol) server that exposes image processing tools to Claude Code. This enables Claude to act as an AI image processing assistant.

### Creating a Project

Use `kintsugi init` to create a new project with the standard directory structure:

```bash
kintsugi init /path/to/my_project --name "My Experiment"
```

This creates:
```
my_project/
├── data/
│   ├── raw/           ← Put your raw images here (cyc001/, cyc002/, etc.)
│   └── processed/     ← Outputs go here automatically
├── notebooks/         ← Working copies of processing notebooks
├── configs/           ← Processing configuration files
├── .claude/           ← Claude Code MCP config (auto-generated)
└── .vscode/           ← VS Code settings (auto-generated)
```

The `.claude/settings.local.json` file is **automatically created** with the MCP server configuration.

## Development Workspace

Use the VS Code multi-root workspace (`kintsugi-dev.code-workspace`) which combines:

| Workspace Folder | Path | Purpose |
|------------------|------|---------|
| KINTSUGI (repo) | `KINTSUGI/` | Source code |
| Test Project | `KINTSUGI/test_data/mini_project/` | 2x2 test dataset |
| Full Project | `KINTSUGI_Projects/.../1904CC1-1L/` | Real data |

Notebooks default to mini_project for testing. Change `PROJECT_DIR` for other projects.

### Adding MCP Support to Existing Projects

If you're not using `kintsugi init`, or need to add MCP support to an existing project, run:

```bash
kintsugi mcp config /path/to/project
```

This automatically creates `.claude/settings.local.json` with the proper configuration. If the file already exists, it will add the KINTSUGI MCP server to your existing configuration.

Use `--print-only` to display the configuration without creating files:
```bash
kintsugi mcp config /path/to/project --print-only
```

### MCP Server Commands

```bash
# Start the MCP server (usually called by Claude Code automatically)
kintsugi mcp start

# List available tools
kintsugi mcp tools
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
- `edf.py` - Extended depth of focus processing (CuPy GPU or NumPy CPU)
- `project.py` - Project structure management (`KintsugiProject` class)
- `deprecation.py` - Deprecation warnings and migration guidance
- `kcorrect_gpu.py` - GPU-accelerated BaSiC illumination correction (CuPy/SciPy)
- `zarr_io.py` - OME-Zarr format I/O with Dask lazy loading
- `parallel_io.py` - Parallel image loading/saving utilities
- `dl_refinement.py` - Deep learning channel refinement

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
- `qc/` - Quality control:
  - `artifact_scanner.py` - Unified `ArtifactScanner` with project/Zarr integration
  - `stripe_artifact.py` - Core detection algorithms (`scan_zstack`, `detect_stripe_artifact`)
  - `stripe_mitigation.py` - Correction methods (notch filter, directional filter, z-plane interpolation)
  - `ImageQC`, `CellQC`, `MarkerQC`, `BatchQC` - Standard QC classes
- `segment/` - Segmentation (SAM wrapper, watershed, postprocessing)
- `claude/` - Claude Code integration (parameter_learning, param_tuner, workflow_state)

**Notebook Modules** (`notebooks/`):
- `Kreg/` - Registration system (15 modules) wrapping VALIS for rigid/non-rigid registration
- `Kstitch/` - GPU-accelerated image stitching with CuPy
- `Kview2/` - 20+ interactive visualization functions for Jupyter
- `KDecon/` - Deconvolution utilities (pure Python, replaces MATLAB)
- `Kio.py` - I/O and parallel processing: channel name parsing, OME-TIFF extraction, parallel tile/stack loading, ProcessingConfig dataclass, ProgressCounter
- `instanseg/` - Instance segmentation models (local customized version)
- `Kutils.py` - Signal isolation utilities (legacy, used by interactive tuners)

**Processing Notebooks** (sequential workflow):
1. `1_Single_Channel_Eval.ipynb` - Parameter tuning and setup
2. `2_Cycle_Processing.ipynb` - Batch cycle processing
3. `3_Signal_Isolation_QC.ipynb` - Signal isolation with integrated QC (Claude-guided or interactive)
4. `4_Segmentation_Analysis.ipynb` - Segmentation, feature extraction, SOM clustering & spatial analysis

## Jupyter Notebook Workflow

**CRITICAL: Autoreload is enabled in all KINTSUGI notebooks:**

```python
%load_ext autoreload
%autoreload 2
```

This means:
- **NEVER tell the user to restart the kernel** after editing Python modules
- **NEVER tell the user to reload/reopen the notebook**
- After editing `.py` files (Kio.py, kintsugi/*.py, etc.), changes take effect automatically on next cell execution
- Just re-run the relevant cell - no restart needed

Only suggest kernel restart for: new package installation, environment variable changes, or C extension recompilation.

## Key Patterns

- **Lazy loading**: `__init__.py` uses lazy imports to avoid loading heavy dependencies at startup
- **Bridge pattern**: `src/kintsugi/*.py` modules wrap notebook implementations for package use
- **Configuration-driven**: Registration and processing pipelines use JSON config files (see `notebooks/config_example.json`)
- **Platform-specific**: Windows/Linux/macOS have different environment files (`envs/`) and dependency handling
- **Parameter learning**: Successful parameters are stored in SQLite databases, indexed by tissue type and marker name
- **MCP integration**: Claude Code can access image processing tools via the MCP server

## Recent Enhancements

**Multi-GPU Support:**
- Explicit device ID parameter for GPU functions (`device_id` parameter)
- Queue-based GPU allocation to prevent OOM during parallel processing
- Robust memory management with automatic OOM retry and chunk size reduction

**Deconvolution Optimizations:**
- Smart padding strategy for optimal FFT performance (power-of-2 friendly sizes)
- Configurable output directory (`decon_dir` parameter)

**Deconvolution Intensity Scaling (Critical Fix):**
- Output scaling now preserves original intensity relationships between channels
- Previously, histogram clipping always stretched output to full 16-bit range (0-65535)
- This caused artificial saturation and inverted relative intensities between channels
- Fix: Scale to `raw_max` (input maximum) instead of 65535
- Critical for quantitative marker expression analysis where relative intensities matter
- Affected code: `notebooks/Kdecon/main.py` lines 283-288

**BaSiC Illumination Correction:**
- Flatfield minimum threshold (`BASIC_FLATFIELD_MIN`) prevents division artifacts
- Default 0.1 limits maximum amplification to 10x, preventing extreme value clipping
- Diagnostic warning when significant flatfield clamping occurs
- Particularly important for sparse markers (e.g., CD8) and out-of-focus z-planes

**Artifact Detection & Mitigation:**
- Unified `ArtifactScanner` class integrated with project structure
- Comparative z-stack analysis for robust detection (avoids dataset-specific threshold issues)
- Supports both TIFF files and Zarr stores
- Generates saveable JSON reports
- Multiple mitigation methods: notch filter, directional filter, z-plane interpolation

Usage:
```python
from kintsugi.project import KintsugiProject
from kintsugi.qc import ArtifactScanner

project = KintsugiProject.load("/path/to/project")
scanner = ArtifactScanner(project)
report = scanner.scan_raw_data()
print(report.summary())
report.save(project.paths.meta / "artifacts.json")
```

**Skills Registry Integration:**
- `/advise` - Search for relevant learnings before starting new work
- `/retrospective` - Save session learnings as new skills
- `/skills` - List all available skills with trigger conditions

**Cache Validation (Kio.py):**
- `validate_stats_cache()` - Check if cached statistics reference files that still exist
- `invalidate_phase_cache()` - Delete cached statistics for a specific processing phase
- Prevents using stale cached data when images are deleted or reprocessed
- Usage:
```python
from Kio import validate_stats_cache, invalidate_phase_cache

# Validate before using cache
if CACHE_FILE.exists():
    with open(CACHE_FILE, 'rb') as f:
        cached_df = pickle.load(f)
    if validate_stats_cache(cached_df, source_dir, phase='stitched'):
        # Safe to use cache
    else:
        invalidate_phase_cache(cache_dir, phase='stitched')
        # Recompute
```

**Quantification Comparisons (Kutils.py):**
- `compare_raw_to_stitched()` - Compare raw tile to corrected/stitched region
- `compare_stitched_to_deconvolved()` - Compare stitched to deconvolved image
- `quantify_processing_pipeline()` - Comprehensive pipeline quantification
- Metrics: SNR improvement, uniformity improvement, sharpness improvement, contrast improvement
- Usage:
```python
from Kutils import quantify_processing_pipeline

result = quantify_processing_pipeline(
    raw_tile, stitched_region, deconvolved_region,
    channel_name='CD3e', cycle='cyc01', zplane=7
)
print(f"Pipeline SNR improvement: {result['pipeline_snr_improvement']:.1f}x")
```

**EDF Smooth Transitions:**
- `blend_depth` parameter - Number of adjacent z-slices to blend for smooth transitions
- `z_smooth_sigma` parameter - Gaussian smoothing for z-index map
- Fixes abrupt transitions in areas of changing contrast
- Usage:
```python
from kintsugi.edf import extended_depth_of_focus_variance

edf = extended_depth_of_focus_variance(
    stack,
    blend_depth=2,      # Blend 2 adjacent slices
    z_smooth_sigma=1.0  # Smooth z-index map
)
```

**Stitched Image Quality Control:**
- New skill `stitched-image-qc` detects saturation and tile grid patterns
- Reprocessing script: `notebooks/reprocess_problematic_images.py`
- Saturation detection: >50% pixels at maximum value
- Tile grid detection: std/mean ratio > 0.8
- See TROUBLESHOOTING.md for detailed diagnostic steps

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

## Git Submodules

This repository uses Git submodules for shared components. See [docs/SUBMODULES.md](docs/SUBMODULES.md) for detailed documentation.

### Quick Reference

```bash
# Clone with submodules
git clone --recurse-submodules https://github.com/smith6jt-cop/KINTSUGI.git

# Initialize submodules after regular clone
git submodule update --init --recursive

# Pull with submodule updates
git pull --recurse-submodules

# Update submodules to latest remote
git submodule update --remote
```

### Current Submodules

| Submodule | Path | Purpose |
|-----------|------|---------|
| Skills_Registry | `Skills_Registry/` | Shared skills and learnings for Claude Code |

### Recommended Aliases

```bash
git config --global alias.clone-all 'clone --recurse-submodules'
git config --global alias.pull-all 'pull --recurse-submodules'
```

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

## Notebook 4: Segmentation & Spatial Analysis

Notebook 4 provides a comprehensive workflow for cell segmentation, feature extraction, and spatial analysis:

### Current Capabilities

**Segmentation (InstanSeg):**
- Nuclear segmentation from DAPI
- Combined cell + nucleus segmentation
- ECM (extracellular matrix) region segmentation via watershed from cell boundaries
- Matched compartment masks (cytoplasm, nuclear, ECM) with shared labels

**Feature Extraction:**
- Uses `napari-simpleitk-image-processing.label_statistics` for:
  - Morphological features (size, perimeter, shape descriptors)
  - Position features (centroids)
  - Intensity statistics (mean, median, min, max, sigma, variance, sum)
  - Moments for texture analysis
- Separate feature DataFrames for cell, nuclear, and ECM compartments

**SOM Clustering (pyFlowSOM):**
- Self-organizing map clustering with configurable grid size
- Learning rate scheduling (start → end)
- Iterative refinement with visualization
- Separate clustering for: cells, nuclei, ECM, combined features

**Spatial Analysis (scanpy + scimap):**
- AnnData integration for single-cell analysis framework
- PCA, UMAP dimensionality reduction
- PAGA graph-based analysis
- Leiden community detection clustering
- Combined SOM + Leiden consensus phenotyping
- Spatial scatter plots with cluster overlays
- Hierarchical clustering for merged phenotypes
- Differential marker analysis (Wilcoxon rank-sum)

### Key Dependencies
- `instanseg` - Deep learning instance segmentation
- `pyFlowSOM` - Self-organizing maps for clustering
- `napari-simpleitk-image-processing` - Feature extraction
- `scanpy` - Single-cell analysis framework
- `scimap` - Spatial analysis toolkit
- `napari` - Interactive visualization

## Migration Notes

The old Notebooks 3 (Signal Isolation) and 5 (DL Channel Refinement) have been replaced with a unified `3_Signal_Isolation_QC.ipynb` that supports both Claude-guided and interactive workflows. See `notebooks/MIGRATION_GUIDE.md` for detailed transition guidance.
