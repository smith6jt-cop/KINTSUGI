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
├── meta/              ← Experiment metadata
│   ├── experiment.json    ← Microscope parameters (auto-generated)
│   └── CHANNELNAMES.txt   ← Channel/marker names (user-provided)
├── notebooks/         ← Working copies of processing notebooks
├── configs/           ← Processing configuration files
├── .mcp.json          ← MCP server config (auto-generated, absolute path)
├── .claude/           ← Claude Code settings (auto-generated)
└── .vscode/           ← VS Code settings (auto-generated)
```

The `.mcp.json` file is **automatically created** with the MCP server definition (using an absolute path to the `kintsugi` executable). The `.claude/settings.local.json` enables the server via `enabledMcpjsonServers`.

### Experiment Metadata

`/meta/experiment.json` stores microscope parameters (tile grid, pixel sizes, wavelengths, NA, RI). Auto-generated during `kintsugi init` with sensible defaults. CLI options: `--tile-rows`, `--tile-cols`, `--xy-pixel-size`, `--z-step-size`, `--numerical-aperture`, `--tissue-ri`.

**CODEX/Akoya compatibility**: `ExperimentConfig.from_dict()` auto-translates CODEX fields (`numCycles`→`n_cycles`, `regionHeight`→`tile_rows`, `xyResolution`→`xy_pixel_size`, etc.). Flat wavelength lists mapped to (excitation, emission) pairs via `_CODEX_FILTER_SETS` in `project.py`. CODEX lacks RI — default **1.44** for uncleared tissue.

### Channel Names

Create `/meta/CHANNELNAMES.txt` with marker names for each channel. Supports multiple formats:

**Simple list (CODEX format):**
```
DAPI-01
Blank
Blank
Blank
DAPI-02
CD31
CD8
CD45
```

**Cycle-prefixed format:**
```
1: DAPI, Blank, Blank, Blank
2: DAPI, CD31, CD8, CD45
```

SLURM scripts and notebooks automatically load channel names from this file.

### Project Initialization Behavior

`kintsugi init` handles existing data intelligently: raw data stays in place, processed data can be selectively deleted, existing projects get `--slurm` added or `--force` to refresh. Use `kintsugi scan /path/to/directory` to preview.

### SLURM Job Submission (HPC)

```bash
kintsugi init /path/to/project --name "My Experiment" --slurm  # New project with SLURM
kintsugi init /path/to/project --slurm    # Add SLURM to existing project
kintsugi slurm init /path/to/project      # Equivalent SLURM-only command
```

Creates `slurm/` with auto-detected HPC settings (`config.sh`), job symlinks, and workflow README. SLURM scripts auto-load metadata from `experiment.json` and `CHANNELNAMES.txt`.

```bash
kintsugi slurm submit .                             # All steps
kintsugi slurm submit . --steps decon,edf           # Specific steps
kintsugi slurm submit . --cycles 1-5                # Specific cycles
kintsugi slurm submit . --dry-run                   # Preview commands
kintsugi slurm status .                             # Check job status
```

See `workflow/CLAUDE.md` for processing modes, multi-account GPU scheduling, Snakemake workflow, registration, and batch processing.

### Snakemake Workflow (Preferred)

```bash
kintsugi workflow config .                         # Generate config from experiment.json
kintsugi workflow check .                          # Verify SLURM resources
kintsugi workflow run .                            # Full pipeline (auto-generates config if missing)
kintsugi workflow run . --dry-run                  # Preview
kintsugi workflow run . --dashboard                # Run with live progress dashboard
kintsugi workflow run . --dashboard -i 15          # Custom refresh interval
kintsugi workflow run . --cycles 1-5               # Specific cycles
kintsugi workflow run . --local --cores 8          # Local (no SLURM)
kintsugi workflow status . --watch                 # Standalone dashboard
```

`workflow run` auto-generates `workflow/config.yaml` from `meta/experiment.json` if the config is missing. The `--dashboard` flag launches snakemake in the background and displays a live progress dashboard; Ctrl+C detaches without killing SLURM jobs.

**Path resolution**: All workflow subcommands use `_resolve_project_dir()` which auto-detects when you're inside a `workflow/` subdirectory and resolves to the parent project root. So `cd project/workflow && kintsugi workflow config .` works correctly.

**Key function**: `generate_workflow_config(project_dir: Path) -> Path` in `cli.py` — standalone, reusable by both `workflow config` and `workflow run`.

### Batch Processing (Multiple Datasets)

**ALWAYS use the CLI command — NEVER write custom batch scripts:**

```bash
kintsugi workflow batch /path/to/projects_dir           # All eligible datasets
kintsugi workflow batch /path/to/projects_dir --dry-run  # Preview what will run
kintsugi workflow batch . -d CX_19-004                   # Single dataset by name
kintsugi workflow batch . -p 2 --detach                  # Background, 2 concurrent
kintsugi workflow batch . --force                        # Reprocess completed datasets
kintsugi workflow stop /path/to/projects_dir             # Stop background batch
```

Eligibility: project has `workflow/config.yaml` + `data/raw/.staged`, and is missing `data/processed/signal_isolated/.snakemake_complete` (unless `--force`). The command validates GPU resources before starting and divides GPU slots across parallel processes.

**NEVER do this:**
- Do NOT write bash scripts that invoke `snakemake` directly
- Do NOT use `subprocess.run(["snakemake", ...])` — use `kintsugi workflow run`
- Do NOT omit `--profile profiles/slurm` — this causes CPU-only execution
- Do NOT use `-j 1` — GPU pipelining requires `-j` >= total GPU slots (currently 2 on maigan)

### GPU-Only Scheduling

All processing (stitch, decon, EDF, registration) runs on GPU. CPU is 5-25x slower per step. QC rules (stitch/decon/edf/registration QC) run on CPU partition — they use pyvips + matplotlib, not GPU compute. The `-j` flag accounts for both GPU processing slots and CPU QC slots.

`workflow run` will **hard-fail** if no GPU slots are detected, if no SLURM profile exists (unless `--local` is explicitly passed), or if every account has zero live availability. There is no silent CPU fallback.

### Multi-Account State (Apr 8 2026)

`BLOCKED_ACCOUNTS` in `src/kintsugi/hpc.py` contains `{brusko, clive}`. **Maigan is the only active account** (gpu_slots=2, cpu_slots=8). Clive was blocked because its QOS pool was throttled to 1 GPU/312.5 GB and is regularly saturated by other users on the shared investment QOS, causing chronic `QOSGrpMemLimit` blockages.

`workflow run` injects per-account live availability into Snakemake via `--config live_accounts=<json>`. The Snakefile's `_build_cycle_assignment`, `_registration_assignment`, and `_qc_cpu_assignment` use the live data to skip accounts saturated by other users. See `workflow/CLAUDE.md` "Resource Pool Calculation" and the `live-aware-account-routing` skill for details.

## Development Workspace

Use the VS Code multi-root workspace (`kintsugi-dev.code-workspace`) which combines:

| Workspace Folder | Path | Purpose |
|------------------|------|---------|
| KINTSUGI (repo) | `KINTSUGI/` | Source code |
| Test Project | `KINTSUGI/test_data/mini_project/` | 2x2 test dataset |
| Full Project | `KINTSUGI_Projects/.../1904CC1-1L/` | Real data |

Notebooks default to mini_project for testing. Change `PROJECT_DIR` for other projects.

### Automatic Project Sync

**Always edit `KINTSUGI/notebooks/` first** — never edit project folders directly. A git post-commit hook auto-syncs notebook modules and Python files to all project folders via `scripts/sync_to_projects.py` (MD5 checksum comparison). Manual sync: `python scripts/sync_to_projects.py [--dry-run|--force]`.

**Auto-discovery**: `sync_to_projects.py` auto-discovers project folders by globbing `KINTSUGI_Projects/*/notebooks`. No manual updates to `DEFAULT_PROJECT_FOLDERS` are needed when new batch projects are added.

### Adding MCP Support to Existing Projects

If you're not using `kintsugi init`, or need to add MCP support to an existing project, run:

```bash
kintsugi mcp config /path/to/project
```

This creates `.mcp.json` at the project root with the MCP server definition (using an absolute path to the `kintsugi` executable) and `.claude/settings.local.json` to enable it. If `.mcp.json` already exists, it will add the KINTSUGI server to your existing configuration.

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
- `subtract_blank` - Autofluorescence subtraction (`method="global"|"weighted"`)
- `analyze_weighted_subtraction` - Preview per-range weights before applying
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

### TissUUmaps Export

Export processed images to web-optimized DZI tile pyramids for [TissUUmaps](https://github.com/TissUUmaps/TissUUmaps) visualization:

```bash
kintsugi workflow export prepare .              # Convert to DZI + generate .tmap
kintsugi workflow export prepare . --dry-run    # Preview what would be exported
kintsugi workflow export prepare . --force      # Reconvert everything
kintsugi workflow export deploy . user@host:/path         # Deploy via rsync
kintsugi workflow export deploy . user@host:/path -n      # Dry-run rsync
kintsugi workflow export status .               # Show export status
```

**What it does:**
1. Discovers registered images (`data/processed/registered/cyc##/*.tif`) and segmentation masks
2. Converts each to DZI tile pyramid (256x256 PNG tiles) via `pyvips.dzsave()`
3. Copies CSV overlays and h5ad files from `segmented/` and `analysis/`
4. Generates `.tmap` project JSON with per-channel color mapping and `compositeMode: "lighter"`
5. Writes `export_manifest.json` for skip-existing on re-runs (mtime + size comparison)

**Output structure:** `data/exported/` with `images/`, `data/`, `regions/`, and `{project}.tmap`

**Key files:** `src/kintsugi/export.py` (core logic), CLI commands in `cli.py` under `@workflow.group("export")`

Blank/autofluorescence channels are excluded by default (use `--include-blanks` to include). 22 common immune markers have predefined color mappings (DAPI=blue, CD3=green, CD8=red, etc.).

## Build and Development Commands

```bash
# Install in development mode
pip install -e ".[dev]"

# Install optional dependency groups via CLI
kintsugi install --list            # Show all 12 available groups
kintsugi install gpu               # Install GPU acceleration
kintsugi install gpu --conda       # Use conda instead of pip (where available)
kintsugi install all               # Unified pip resolve via pyproject.toml extras

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
- `cli.py` - Click-based CLI with subcommands: check, register, template, info, mcp (config/start/tools/pretrain), init, install, workflow (config/check/run/export). `_resolve_project_dir()` auto-detects `workflow/` subdirectories for all workflow subcommands
- `kreg.py`, `kstitch.py`, `kview2.py` - Bridge modules to notebook implementations
- `deps.py` - Runtime dependency validation (Python packages, libvips, CUDA). `OPTIONAL_GROUPS` dict is the single source of truth for all 12 install groups (gpu, viz, dl, analysis, bio, claude, dev, docs, kronos, denoise, rapids, full)
- `edf.py` - Extended depth of focus processing (CuPy GPU or NumPy CPU)
- `project.py` - Project structure management (`KintsugiProject` class), `_resolve_kintsugi_executable()` for absolute-path MCP config
- `deprecation.py` - Deprecation warnings and migration guidance
- `kcorrect_gpu.py` - GPU-accelerated BaSiC illumination correction (CuPy/SciPy)
- `zarr_io.py` - OME-Zarr format I/O with Dask lazy loading
- `parallel_io.py` - Parallel image loading/saving utilities
- `dl_refinement.py` - Deep learning channel refinement
- `export.py` - TissUUmaps export: DZI tile conversion, .tmap generation, rsync deploy
- `vessel3d.py` - 3D vessel segmentation: Frangi vesselness, presets, marker discovery, multichannel combination, skeletonization, graph morphometry
- `vessel3d_viz.py` - Vessel visualization: ortho views, mask overlays, feature plots

**MCP Server** (`src/kintsugi/mcp/`):
- `server.py` - FastMCP server implementation
- `tools/signal_isolation.py` - Signal processing tools
- `tools/quality_assessment.py` - QC tools
- `tools/visualization.py` - Image preview tools
- `tools/workflow.py` - Workflow management tools
- `tools/learning.py` - Parameter learning with SQLite storage

**Pure Python Modules** (`src/kintsugi/`):
- `signal/` - Autofluorescence subtraction (global + weighted multi-range), bootstrap learning (`bootstrap.py`), multi-project batch orchestration (`batch_multi.py`)
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
- `Kprocess.py` - GPU-accelerated QC statistics and visualization: `compute_*_stats_gpu()`, `collect_*_stats_parallel()`, `plot_summary_heatmaps()`, `run_*_qc()` convenience functions
- `instanseg/` - Instance segmentation models (local customized version)
- `Kutils.py` - Signal isolation utilities (legacy, used by interactive tuners)

**Processing Notebooks** (sequential workflow):
1. `1_Single_Channel_Eval.ipynb` - Parameter tuning and setup
2. `2_Cycle_Processing.ipynb` - Batch cycle processing
3. `3_Signal_Isolation_QC.ipynb` - Signal isolation with integrated QC (Claude-guided or interactive)
4. `4_Segmentation_Analysis.ipynb` - Segmentation, feature extraction, SOM clustering & spatial analysis

## Jupyter Notebook Workflow

**CRITICAL**: Autoreload is enabled in all KINTSUGI notebooks — **NEVER tell the user to restart the kernel** after editing Python modules. Just re-run the relevant cell.

See `notebooks/CLAUDE.md` for full details on autoreload behavior, troubleshooting "Function Not Found" errors, and cell execution order.

## Key Patterns

- **Lazy loading**: `__init__.py` uses lazy imports to avoid loading heavy dependencies at startup
- **Bridge pattern**: `src/kintsugi/*.py` modules wrap notebook implementations for package use
- **Configuration-driven**: Registration and processing pipelines use JSON config files (see `notebooks/config_example.json`)
- **Platform-specific**: Windows/Linux/macOS have different environment files (`envs/`) and dependency handling
- **Parameter learning**: Successful parameters are stored in SQLite databases, indexed by tissue type and marker name
- **MCP integration**: Claude Code can access image processing tools via the MCP server

## Critical Fixes & Enhancements (Jan–Feb 2026)

These bugs have been fixed in the current codebase. Listed here for context when investigating image quality issues.

| Component | Issue | Fix | Key Detail |
|-----------|-------|-----|------------|
| Deconvolution | Horizontal banding (Gibbs ringing) | Enabled Tukey window edge apodization | `_create_tukey_window_3d` existed but was never called (`Kdecon/deconvolution.py`) |
| Deconvolution | Intensity scaling saturated output | Scale to `raw_max` instead of 65535 | Preserves relative channel intensities (`Kdecon/main.py`) |
| BaSiC Correction | Negative flatfield on sparse markers | Falls back to uniform flatfield | `flatfield_min=0.1` clamp in `transform()` (`kcorrect_gpu.py`) |
| Stitching | Vertical stripe artifacts (~30px) | Reprocess with current code | Stripes from old pipeline code, NOT raw data. Script: `scripts/reprocess_striped_zplanes.py` |
| EDF | Detail loss from sigma smoothing | Smooth variance map, not input image | CLIJ2-matching defaults: radius=2, sigma=10 (`edf.py`) |
| Vessel3D | Frangi Ra formula wrong (Ra=\|λ₁\|/\|λ₂\|) | Ra=\|λ₂\|/\|λ₃\| per Frangi 1998 eq. 11 | Tube response was suppressed — small/faint vessels missed (`vessel3d.py`) |
| Vessel3D | L4 GPU OOM on isotropic volumes | VRAM guard + CPU fallback when <40 GB | Queries `cp.cuda.Device().mem_info` before GPU path (`vessel3d.py`) |

**GPU Processing:**
- Multi-GPU with `device_id` parameter, queue-based allocation, automatic OOM retry
- EDF: `blend_depth` and `z_smooth_sigma` for smooth z-transitions

**Quality Control:**
- `ArtifactScanner` — unified artifact detection for TIFF/Zarr with JSON reports
- `QualityGate` — pre-processing validation (stripe, saturation, low signal detection)
- `validate_stats_cache()` / `invalidate_phase_cache()` in `Kio.py` — stale cache prevention
- `quantify_processing_pipeline()` in `Kutils.py` — SNR/uniformity/sharpness metrics

**When investigating artifacts**: Always verify raw data is clean before attributing to data quality. Most historical "artifacts" were caused by the deconvolution apodization bug (now fixed). Compare fresh reprocessing vs existing files first.

**Skills Registry:** `/advise`, `/retrospective`, `/skills` — search and save learnings across sessions

## Registration

See `workflow/CLAUDE.md` for comprehensive registration documentation: tuned non-rigid parameters, three critical bug fixes (params unpacking, coarse NR, Valis init), VALIS dimension parameters, sigma math, smoothing methods, validated results, batch re-registration, and tuning guide.

## Signal Isolation

See `src/kintsugi/signal/CLAUDE.md` for weighted autofluorescence subtraction (per-intensity-range weights), single-project batch signal isolation (recipe-driven multi-step processing, auto method selection, background cleaning, parameter learning, QC reporting), and multi-project batch orchestration (cascading recipe resolution, cross-project reports).

## Batch Processing (Multi-Dataset)

See `workflow/CLAUDE.md` for batch processing documentation: data staging, cleanup with QC guard, full pipeline lifecycle, and sentinel validation.

**Sentinel validation** (`scripts/create_si_sentinels.py`): Validates signal isolation output (manifest JSON + TIF existence/size) and creates missing `.snakemake_complete` sentinels for projects processed outside Snakemake (e.g., by batch scripts). Supports `--dry-run`. Used to promote 25 batch-processed projects to "signal_isolated" status (Feb 2026).

## Dependencies

Core: numpy<2.0, scipy, pandas, scikit-image, opencv-contrib-python-headless, pyvips, seaborn, pqdm

13 optional install groups defined in `deps.py` `OPTIONAL_GROUPS` (single source of truth, also used by `pyproject.toml` extras). Install via `kintsugi install <group>` (supports `--conda` for groups with conda recipes):
- `gpu` - PyTorch + CuPy for GPU acceleration
- `viz` - Napari visualization
- `dl` - Deep learning segmentation (InstanSeg)
- `analysis` - scanpy, scimap for spatial analysis
- `bio` - Bio formats I/O (OME-TIFF, LIF, etc.)
- `claude` - MCP server for Claude Code integration
- `dev` - Development tools (pytest, ruff, black, mypy)
- `docs` - Documentation (Sphinx)
- `kronos` - KRONOS foundation model (prints post-install note about cloning KRONOS repo)
- `denoise` - Advanced denoising (N2V, CARE)
- `optimize` - Parameter optimization (Optuna + SMO)
- `workflow` - Snakemake workflow orchestration (HPC/SLURM)
- `rapids` - RAPIDS GPU-accelerated data science
- `full` - All optional features (composite)

`kintsugi install all` uses a single `pip install -e ".[gpu,viz,dl,...]"` pass via pyproject.toml extras to avoid cascading dependency conflicts. Skips `full` (composite) and `rapids` (needs NVIDIA channel). Falls back to sequential install + constraint repair if pyproject.toml not found. Conda mode (`--conda`) installs conda groups first, then remaining via pip extras.

**External requirements**: libvips (native library). Java/Maven no longer required.

### Dependency Safety

**Constraint guards** prevent known-bad dependency combinations:
- `constraints.txt` at repo root enforces `numpy<2.0` — auto-injected into all `kintsugi install` pip commands
- All torch-using groups (`gpu`, `dl`, `denoise`, `kronos`) use `--index-url https://download.pytorch.org/whl/cu128` to avoid CPU-only torch (cu128 supports Blackwell/B200 sm_100)
- `analysis`, `bio`, `kronos` groups have explicit `numpy>=1.24.0,<2.0.0` in pyproject.toml extras
- Pre-install guards warn before installing groups that would break existing packages

**Validation**: `kintsugi check --strict` catches:
- numpy >= 2.0 installed (ERROR)
- CPU-only PyTorch build (ERROR)
- CuPy missing CUDA runtime libraries (ERROR)
- SLURM TRES patch not applied (ERROR on HPC)

**SLURM TRES patch**: `kintsugi patch slurm` auto-applies the `SLURM_TRES_PER_TASK` fix. Also applied automatically by `kintsugi install all` and `kintsugi install workflow`.

See `docs/DEPENDENCY_GUIDE.md` for full troubleshooting guide.

### CuPy & HPC Cache (IMPORTANT)

**CuPy IS INSTALLED** (`cupy-cuda12x`). Do NOT reinstall. `gpu.cupy_available` returns `False` on login nodes (no GPU hardware) — this is normal. Use `gpu.cupy_installed` for import-only check. Only reinstall via `kintsugi install gpu` after a full env rebuild.

**HPC cache**: SLURM jobs use account-specific caches (`~/.use_conda_{account}.sh` → `~/.cache_redirect.sh`). Never install packages or modify caches directly.

## Git Submodules

`Skills_Registry/` is a Git submodule. Clone with `--recurse-submodules` or run `git submodule update --init --recursive`.

## Testing

Tests are in `tests/` with fixtures in `conftest.py`. Key fixtures: `sample_image`, `sample_multichannel_image`, `sample_stack`, `sample_tiff`, `sample_config`, `temp_dir`.

CI runs on Linux (ubuntu-latest) with Python 3.10-3.12.

**Suite status** (Mar 2026): 814 collected, 10 skipped (GPU hardware + optional deps).

**GPU skip pattern**: On HPC login nodes, CuPy is installed (`CUPY_AVAILABLE=True`) but no GPU hardware exists. Tests that need actual GPU hardware must use `check_gpu()` from `kcorrect_gpu.py`, not `CUPY_AVAILABLE`:
```python
from kintsugi.kcorrect_gpu import check_gpu
GPU_HARDWARE_AVAILABLE, _reason = check_gpu()

@pytest.mark.skipif(not GPU_HARDWARE_AVAILABLE, reason="No GPU hardware available")
```

**Dev dependencies**: `pytest-asyncio` is required for MCP tool tests (27 async tests). Install with `pip install pytest-asyncio` or `pip install -e ".[dev]"`.

## Code Style

- Line length: 100 characters
- Formatter: Black
- Linter: Ruff (with isort, pyupgrade, flake8-bugbear)
- Python 3.10+ required

## Release and Versioning

Uses [Conventional Commits](https://www.conventionalcommits.org/) with automatic semantic versioning. Format: `<type>(<scope>): <description>` — types: `feat` (MINOR), `fix` (PATCH), `docs`, `refactor`, `perf`, `test`, `build`, `ci`, `chore`. Add `!` for breaking changes (MAJOR).

```bash
pre-commit install && pre-commit install --hook-type commit-msg  # Enable hooks
python scripts/release.py --auto                                  # Create release
python scripts/release.py --dry-run --auto                       # Preview
```

## Notebook 4: Segmentation & Spatial Analysis

See `notebooks/CLAUDE.md` for Notebook 4 capabilities (InstanSeg segmentation, feature extraction, SOM clustering, spatial analysis with scanpy/scimap), key dependencies, and migration notes.

## Subdirectory Documentation

Feature-specific documentation lives in CLAUDE.md files near the code:
- `src/kintsugi/CLAUDE.md` — KRONOS foundation model, 3D vessel segmentation, pipeline-aware cleanup
- `src/kintsugi/signal/CLAUDE.md` — Weighted AF subtraction, single-project batch isolation, multi-project batch orchestration
- `workflow/CLAUDE.md` — Snakemake workflow, registration, batch processing, SLURM scheduling
- `notebooks/CLAUDE.md` — Jupyter autoreload, Notebook 4 segmentation/analysis
