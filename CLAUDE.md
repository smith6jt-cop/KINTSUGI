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
├── .claude/           ← Claude Code MCP config (auto-generated)
└── .vscode/           ← VS Code settings (auto-generated)
```

The `.claude/settings.local.json` file is **automatically created** with the MCP server configuration.

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

For HPC clusters with SLURM, add job submission support during project creation:

```bash
kintsugi init /path/to/project --name "My Experiment" --slurm
```

This creates an additional `slurm/` directory:
```
my_project/
└── slurm/
    ├── config.sh      ← Pre-populated with auto-detected HPC settings
    ├── jobs/          ← Symlink to KINTSUGI job scripts
    └── README.md      ← Complete 10-step workflow guide
```

**Auto-detection**: The system automatically detects:
- SLURM account from `SLURM_ACCOUNT` or `SBATCH_ACCOUNT` environment variables
- GPU type from `nvidia-smi`
- Conda environment from `CONDA_DEFAULT_ENV`

**Metadata loading**: SLURM job scripts automatically load parameters from:
1. `/meta/experiment.json` - Microscope parameters (tile grid, pixel sizes, wavelengths)
2. `/meta/CHANNELNAMES.txt` - Channel/marker names (used for EDF output file naming)
3. `slurm/config.sh` environment variables - Fallback if metadata files don't exist

This eliminates the need to manually configure `config.sh` for most parameters.

**SLURM output naming**: EDF uses marker names from `CHANNELNAMES.txt` (e.g., `CD3.tif`), falling back to `CH#`. Registration preserves those names with warped images. SLURM jobs need `KINTSUGI_DIR/notebooks` on `sys.path` for `Kio` imports.

**Memory**: GPU jobs = 48 GB RAM (CuPy does FFT in GPU memory). CPU jobs = 128 GB (SciPy float64). Snakefile lambdas auto-route.

**Add SLURM to existing project** (two equivalent methods):
```bash
kintsugi init /path/to/project --slurm    # Auto-detects existing project, adds SLURM
kintsugi slurm init /path/to/project      # Explicit SLURM-only command
```

**Submit jobs**:
```bash
kintsugi slurm submit /path/to/project              # All steps
kintsugi slurm submit . --steps decon,edf           # Specific steps
kintsugi slurm submit . --cycles 1-5                # Specific cycles
kintsugi slurm submit . --dry-run                   # Preview commands
```

**Check status**:
```bash
kintsugi slurm status .
```

See `workflow/CLAUDE.md` for detailed documentation on processing modes (notebook vs SLURM), multi-account resource pool calculation, concurrent GPU/CPU architecture, Snakemake workflow design decisions, and registration rule configuration.

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
- `cli.py` - Click-based CLI with subcommands: check, register, template, info, mcp (config/start/tools/pretrain), init, workflow (config/check/run/export)
- `kreg.py`, `kstitch.py`, `kview2.py` - Bridge modules to notebook implementations
- `deps.py` - Runtime dependency validation (Python packages, libvips, CUDA)
- `edf.py` - Extended depth of focus processing (CuPy GPU or NumPy CPU)
- `project.py` - Project structure management (`KintsugiProject` class)
- `deprecation.py` - Deprecation warnings and migration guidance
- `kcorrect_gpu.py` - GPU-accelerated BaSiC illumination correction (CuPy/SciPy)
- `zarr_io.py` - OME-Zarr format I/O with Dask lazy loading
- `parallel_io.py` - Parallel image loading/saving utilities
- `dl_refinement.py` - Deep learning channel refinement
- `export.py` - TissUUmaps export: DZI tile conversion, .tmap generation, rsync deploy
- `vessel3d.py` - 3D vessel segmentation: Frangi vesselness, skeletonization, graph morphometry
- `vessel3d_viz.py` - Vessel visualization: ortho views, mask overlays, feature plots

**MCP Server** (`src/kintsugi/mcp/`):
- `server.py` - FastMCP server implementation
- `tools/signal_isolation.py` - Signal processing tools
- `tools/quality_assessment.py` - QC tools
- `tools/visualization.py` - Image preview tools
- `tools/workflow.py` - Workflow management tools
- `tools/learning.py` - Parameter learning with SQLite storage

**Pure Python Modules** (`src/kintsugi/`):
- `signal/` - Autofluorescence subtraction (global + weighted multi-range), bootstrap learning (`bootstrap.py`)
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

## Critical Fixes & Enhancements (Jan 2026)

These bugs have been fixed in the current codebase. Listed here for context when investigating image quality issues.

| Component | Issue | Fix | Key Detail |
|-----------|-------|-----|------------|
| Deconvolution | Horizontal banding (Gibbs ringing) | Enabled Tukey window edge apodization | `_create_tukey_window_3d` existed but was never called (`Kdecon/deconvolution.py`) |
| Deconvolution | Intensity scaling saturated output | Scale to `raw_max` instead of 65535 | Preserves relative channel intensities (`Kdecon/main.py`) |
| BaSiC Correction | Negative flatfield on sparse markers | Falls back to uniform flatfield | `flatfield_min=0.1` clamp in `transform()` (`kcorrect_gpu.py`) |
| Stitching | Vertical stripe artifacts (~30px) | Reprocess with current code | Stripes from old pipeline code, NOT raw data. Script: `scripts/reprocess_striped_zplanes.py` |
| EDF | Detail loss from sigma smoothing | Smooth variance map, not input image | CLIJ2-matching defaults: radius=2, sigma=10 (`edf.py`) |

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

## Registration: Tuned Non-Rigid for CODEX (Feb 2026)

Registration uses rigid + non-rigid alignment with tuned smoothing parameters. Three critical bugs were fixed:

1. **Params not unpacked** (ROOT CAUSE): `serial_non_rigid.py:438` passed `non_rigid_reg_class(params=init_kwargs)` instead of `**init_kwargs` — all smoothing params silently ignored. The earlier "rigid-only is best" conclusion (11/12 datasets worse with non-rigid) was based on this broken code — all tests ran with unsmoothed OpticalFlowWarper.
2. **Coarse non-rigid at 1024px**: VALIS `max_non_rigid_registration_dim_px` defaults to `DEFAULT_MAX_PROCESSED_IMG_SIZE = 1024` (NOT `DEFAULT_MAX_NON_RIGID_REG_SIZE = 3000`). This coarse pass poisoned alignment before the 4096px micro pass. Fix: disable coarse NR by not passing `non_rigid_registrar_cls` to Valis init.
3. **Missing non-rigid config in Valis init**: `registration.py` didn't pass `non_rigid_registrar_cls` or `non_rigid_reg_params` to Valis, so `register()` used defaults (unsmoothed OpticalFlowWarper at 1024px).

**Three VALIS dimension parameters** (often confused):
- `max_image_dim_px` (default 1024) — max saved output image size
- `max_processed_image_dim_px` (default 1024) — rigid feature detection resolution
- `max_non_rigid_registration_dim_px` (default 1024) — coarse non-rigid field resolution

**Current architecture**: Coarse NR is **disabled** (`non_rigid_registrar_cls=None` in Valis init). Non-rigid runs exclusively via `register_micro()` at 4096px with tuned smoothing. Rigid feature detection at 4096px for better matching.

**Optical flow algorithm**: CUDA TVL1 (default) or CPU DeepFlow (fallback). Both use default OpenCV parameters — no tuning knobs exposed.

**Sigma math**: `sigma_pixels = sigma_ratio * max(image_dim_px)`. At 4096px: ratio 0.01 = 41px sigma, ratio 0.005 = 20px sigma, ratio 0.05 = 205px sigma (extreme over-smoothing, washes out all local corrections).

**Smoothing methods**: Only `"gauss"` and `None` work. `"inpaint"` and `"regularize"` are broken in the VALIS code (untested paths). With `smoothing_method="gauss"`, `n_grid_pts` and `fold_penalty` are **no-ops** — they only apply to the broken "regularize"/"inpaint" methods.

**Current parameters** (`workflow/config.yaml`):
```yaml
registration:
  rigid_only: false           # Non-rigid enabled with tuned smoothing
  imgs_ordered: true          # Keep sequential cycle order (VALIS default reorders by similarity)
  align_to_reference: true    # Direct alignment to reference (prevents serial error accumulation)
  reference_cycle: 1
  feature_detector: "VggFD"
  max_image_dim: 4096         # Rigid feature detection resolution (was 2048, improved feature matching)
  non_rigid:
    max_dim: 4096             # Displacement field resolution for register_micro()
    smoothing_method: "gauss" # Prevents noisy displacement fields
    sigma_ratio: 0.01         # sigma = 0.01 * 4096 = ~41px. Default in OpticalFlowWarper is 0.005 (20px)
    n_grid_pts: 50            # No-op with "gauss" (only used by broken "regularize"/"inpaint")
    fold_penalty: 1.0e-6      # No-op with "gauss" (only used by broken "regularize")
    normalize_dimensions: true # Normalize anisotropic images to square aspect ratio for uniform flow
```

**Validated Results** (batch re-registration in progress, 4/16 completed):

| Dataset | Cycles | TIFFs Warped | Time (min) |
|---------|--------|-------------|------------|
| CX_19-001_SP_CC2-A28 | 13 | 52 | 67.9 |
| CX_19-002_lymph-node_R1 | 9 | 36 | 29.5 |
| CX_19-002_lymph-node_R3 | 9 | 36 | 27.6 |
| CX_19-003_lymph-node_R1 | 9 | 36 | 40.8 |

Non-rigid improves on rigid for the vast majority of cycles. Confirms the params unpacking fix (bug #1 above) was the ROOT CAUSE of the earlier "rigid-only is best" conclusion across all 47 datasets.

**Key files:** `workflow/scripts/registration.py` (wrapper), `notebooks/Kreg/registration.py` (VALIS Valis class, line 1771 for dimension defaults), `notebooks/Kreg/serial_non_rigid.py` (line 438, params bug), `notebooks/Kreg/non_rigid_registrars.py` (OpticalFlowWarper)

**Registration QC**: Green/magenta DAPI overlay at **full resolution** (1000x1000 crop) — whole-image thumbnails are too zoomed out to evaluate nuclei overlap. White/gray = aligned, color fringing = misaligned. **Average metrics are meaningless** — poor registration is localized; use spatial NCC heatmaps or targeted overlays.

**Error handling**: Registration failures now raise loudly (no silent fallback copy of EDF images). The old fallback_copy behavior was removed because it masked registration problems.

**VALIS metrics**: `D` values are at processing resolution, not original. Use `rTRE = D / max(processed_shape)` for cross-experiment comparison at different `max_image_dim` settings.

**Batch re-registration**: Remove sentinel + `registration_data/` + `cyc*/` under `registered/`, update config, re-run Snakemake. Script: `/blue/maigan/smith6jt/reregister_all.sh`. Key features:
- Wave-based parallel execution: 5 projects per wave (3 clive + 2 maigan GPUs)
- Targets only registration + qc_registration (CPU-only) to avoid GPU QC contention
- Account distribution via config.yaml reordering (controls `_registration_assignment()`)
- Idempotent: skips projects with `non_rigid_smoothing=gauss` in sentinel

**Tuning guide** (`sigma_pixels = sigma_ratio * max_dim`): Over-smoothing → decrease `sigma_ratio` to 0.005 (20px at 4096). Under-correcting → increase `max_dim` to 8192. Over-warping → increase `sigma_ratio` to 0.02 (82px at 4096). Do NOT use `smoothing_method: "regularize"` — it's broken.

## Weighted Autofluorescence Subtraction (Feb 2026)

Replaces the single global `blank_scale_factor` with per-intensity-range weights. Protects dim signal (FOXP3, CD163) while aggressively removing bright autofluorescence (collagen, lipofuscin).

**Algorithm**: Segment signal histogram into N ranges (default 5: background, very_dim, dim, medium, bright). Per-range weight is computed from signal-to-AF ratio — dim regions where AF > signal get weight 0.3-0.5 (gentle subtraction), bright signal regions get weight 1.0+ (aggressive removal). Cosine transitions between ranges prevent discontinuities.

**Key files:**
- `signal/autofluorescence.py` — `subtract_autofluorescence_weighted()`, `compute_intensity_ranges()`, `build_weight_map()`, `analyze_for_weighted_subtraction()`, `compute_weighted_subtraction_quality()`, Dask variant
- `signal/utils.py` — `adaptive_range_boundaries()`, `smooth_membership()`
- `signal/subtractor.py` — `IntensityRange`, `WeightedSubtractionParameters` dataclasses, `method="weighted"` in `AutofluorescenceSubtractor`
- `signal/bootstrap.py` — `bootstrap_from_project()` for batch pretraining
- `claude/parameter_learning.py` — `algorithm_version` column, range aggregation
- `mcp/tools/signal_isolation.py` — `analyze_weighted_subtraction()` tool, `method` param on `subtract_blank()`

**CLI**: `kintsugi mcp pretrain <project_dir> [--tissue-type TYPE] [--dry-run]` — scans EDF outputs for signal/blank pairs, computes weighted parameters, records to learning database with `algorithm_version="weighted_v1"`.

**Backwards compatible**: `subtract_blank()` defaults to `method="global"` (original behavior). Uniform weights=1.0 produces identical output to global method.

## Batch Processing (Multi-Dataset)

KINTSUGI supports automated batch processing of multiple CODEX datasets (47 in manifest: spleen, lymph node, thymus). See `workflow/CLAUDE.md` for batch processing details including data staging, storage architecture, and sentinel file reference.

**Data staging**: Two methods depending on source:
- **Orange storage** (spleen/LN datasets): `stage_datasets.sh` via rsync SLURM jobs (CPU partition, no GPU conflict)
- **Globus** (thymus datasets from PATH lab): `stage_datasets_globus.py` — transfers from `path.ahc.ufl.edu` SMB share via cifs mount at `/mnt/ahc_share/SHARE/HuBMAP/`. **CRITICAL: Use cifs mount, NOT GVFS** — GVFS causes EOF errors under sustained I/O. Always `globus endpoint activate <UUID>` before transfers.

**Quick reference — processing all staged datasets:**
```bash
tmux new -s batch
bash /blue/maigan/smith6jt/run_all_workflows.sh           # all staged datasets
bash /blue/maigan/smith6jt/run_all_workflows.sh --dry-run  # preview
bash /blue/maigan/smith6jt/run_all_workflows.sh --dataset CX_19-002_lymph-node_R1  # single
```

Processing is sequential (one dataset at a time) because all datasets share the same 5 GPU slots (across clive and maigan). GPU-only scheduling is used — CPU is 5-25x slower per step (~13x full cycle), queuing for GPU always wins. Each Snakemake instance uses the full GPU resource pool. Completed datasets and channels are automatically skipped on re-run.

**Cleanup with QC guard:**
```bash
bash /blue/maigan/smith6jt/cleanup_datasets.sh              # checks QC sentinels, prompts review
bash /blue/maigan/smith6jt/cleanup_datasets.sh --force      # skip prompts (after initial review)
bash /blue/maigan/smith6jt/cleanup_datasets.sh --dry-run    # preview without deleting
```

Cleanup verifies stitch/decon/edf QC sentinel files exist before allowing deletion of intermediates and raw data. Registration QC is NOT required — registration reads EDF outputs, not intermediates. `--force` only skips the interactive review prompt; sentinel checks are always enforced. Never use `--force` without confirming QC has been reviewed.

**EDF blank channel workaround**: Some datasets have all-zeros channels (e.g., CC3-C CH2). EDF script exits non-zero when a channel has 0 valid z-slices, even if other channels succeed. Workaround: manually create `.snakemake_complete` sentinels for cycles where all non-blank channels processed correctly, then run QC.

See `workflow/CLAUDE.md` for full batch pipeline lifecycle.

## Dependencies

Core: numpy<2.0, scipy, pandas, scikit-image, opencv-contrib-python-headless, pyvips, valis-wsi

Optional groups in pyproject.toml:
- `[gpu]` - PyTorch + CuPy for GPU acceleration
- `[claude]` - MCP server for Claude Code integration
- `[denoise]` - PyTorch for N2V/CARE denoising
- `[viz]` - Napari visualization
- `[analysis]` - scanpy, scimap for spatial analysis
- `[kronos]` - KRONOS foundation model integration (PyTorch, h5py, umap-learn, scanpy)
- `[bio]` - Bio formats I/O (OME-TIFF, LIF, etc.)
- `[dl]` - Deep learning segmentation (InstanSeg)
- `[full]` - All optional dependencies
- `[java]` - (DEPRECATED) JPype + PyImageJ for BioFormats

**External requirements**: libvips (native library). Java/Maven no longer required.

### CuPy & HPC Cache (IMPORTANT)

**CuPy IS INSTALLED** (`cupy-cuda12x`). Do NOT reinstall. `gpu.cupy_available` returns `False` on login nodes (no GPU hardware) — this is normal. Use `gpu.cupy_installed` for import-only check. Only reinstall via `kintsugi install gpu` after a full env rebuild.

**HPC cache**: SLURM jobs use account-specific caches (`~/.use_conda_{account}.sh` → `~/.cache_redirect.sh`). Never install packages or modify caches directly.

## Git Submodules

`Skills_Registry/` is a Git submodule. Clone with `--recurse-submodules` or run `git submodule update --init --recursive`.

## Testing

Tests are in `tests/` with fixtures in `conftest.py`. Key fixtures: `sample_image`, `sample_multichannel_image`, `sample_stack`, `sample_tiff`, `sample_config`, `temp_dir`.

CI runs on Windows/Linux/macOS with Python 3.10-3.12.

**Suite status** (Feb 2026): 390 passed, 16 skipped (GPU hardware + dask_image optional dep).

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

## KRONOS Foundation Model Integration (Feb 2026)

[KRONOS](https://github.com/mahmoodlab/KRONOS) is a panel-agnostic foundation model for spatial proteomics (Harvard/BWH, CC-BY-NC-ND-4.0). Trained via self-supervised learning on 47M single-marker patches across 175 protein markers.

**Pipeline position**: Post-registration analysis — runs after EDF/registration produces aligned, marker-named 2D TIFFs.

**Key files:**
- `src/kintsugi/kronos/__init__.py` — Lazy loading, dependency checks
- `src/kintsugi/kronos/model.py` — `KronosModel`, `KronosConfig`, `EmbeddingResult` (HDF5 save/load)
- `src/kintsugi/kronos/markers.py` — `MarkerMapper` maps CHANNELNAMES.txt to KRONOS's 175 markers
- `src/kintsugi/kronos/embeddings.py` — `KronosEmbedder` loads registered images, tiles, normalizes, runs inference
- `src/kintsugi/kronos/analysis.py` — Clustering (Leiden/KMeans), UMAP/PCA, spatial search, cross-dataset comparison
- `notebooks/5_KRONOS_Analysis.ipynb` — Standalone analysis notebook

**Usage:**
```python
from kintsugi.kronos import KronosEmbedder, MarkerMapper
mapper = MarkerMapper.from_project("./my_project", cache_dir="./model_assets")
embedder = KronosEmbedder()
result = embedder.embed_project("./my_project", mapper)
result.save("./embeddings.h5")
```

**Prerequisites**: `pip install -e KRONOS` (from GitHub clone), HuggingFace access for model weights.

**Model specs**: ViT-S/16, 175 MB weights, input `(batch, markers, 256, 256)`, output 3 embedding types (patch/marker/token). Per-marker z-score normalization via `marker_metadata.csv`.

**KRONOS does NOT replace**: stitching, deconvolution, EDF, registration, segmentation masks, or signal isolation — it consumes their outputs.

## 3D Vessel Segmentation (Feb 2026)

Segments blood and lymphatic vessels from deconvolved 3D z-stacks using Frangi vesselness filtering, morphological cleanup, skeletonization, and graph-based morphometry.

**Pipeline position**: Step 3.5 — after deconvolution, before EDF. Uses the full 3D z-stack (EDF collapses z). Runs in parallel with EDF; does not block downstream steps.

**Key files:**
- `src/kintsugi/vessel3d.py` — Core processing: `compute_vesselness_frangi()`, `binarize_vessel_mask()`, `skeletonize_vessels()`, `prune_skeleton()`, `analyze_vessel_graph()`, `export_vessel_results()`, `segment_vessels_3d()` (high-level pipeline)
- `src/kintsugi/vessel3d_viz.py` — Visualization: `ortho_view()`, `overlay_mask_on_raw()`, `render_skeleton_mip()`, `plot_vessel_features()`
- `notebooks/2.5_Vessel_3D_Segmentation.ipynb` — Interactive notebook (7 sections with verification checkpoints)
- `slurm/jobs/03b_vessel3d.sh` — HPC batch job script (Frangi + skeleton, automated)

**Data classes:**
- `VesselSpacing(xy, z)` — Physical voxel spacing; `from_experiment()` loads from experiment.json
- `VesselSegmentationResult` — Container for mask, skeleton, features, vesselness map

**Usage:**
```python
from kintsugi.vessel3d import segment_vessels_3d, VesselSpacing, export_vessel_results

spacing = VesselSpacing(xy=0.377, z=1.5)  # or VesselSpacing.from_experiment(config)
result = segment_vessels_3d(
    volume, spacing=spacing, sigmas=[1, 2, 4, 8],
    device='auto', marker_name='CD31',
)
export_vessel_results(result, "data/processed/vessel_3d/")
```

**Pipeline steps:**
1. Preprocess: isotropic z-resampling (scipy/CuPy zoom) + light Gaussian denoising
2. Frangi vesselness: multi-scale Hessian eigenvalue analysis (GPU-accelerated)
3. Binarize: Otsu threshold + `remove_small_objects` + morphological closing + hole filling
4. Skeletonize: Lee (1994) 3D thinning via `skimage.morphology.skeletonize`
5. Graph analysis: `skan.Skeleton` + `skan.summarize` + distance transform for radius

**Output files** (in `data/processed/vessel_3d/`):
- `binary_mask_{marker}.tif` — 3D vessel mask (isotropic)
- `binary_mask_{marker}_native.tif` — 3D mask at native resolution
- `skeleton_{marker}.tif` — 3D skeleton
- `vessel_features_{marker}.csv` — Per-segment morphometry (length, diameter, tortuosity)
- `vessel_graph_{marker}.graphml` — NetworkX graph
- `vessel_2d_projection_{marker}.tif` — MIP for Notebook 4 spatial analysis

**Dependencies**: `skan>=0.11.0` and `networkx>=3.0` added to the `analysis` optional group. GPU acceleration via CuPy (existing `gpu` group). Skeletonization is CPU-only.

**SLURM**: `03b_vessel3d.sh` runs between `03_deconvolution.sh` and `04_edf.sh`. Environment variables: `VESSEL_CYCLE`, `VESSEL_CHANNEL`, `VESSEL_MARKER`, `VESSEL_SIGMAS`, `VESSEL_MIN_SIZE`. Recommended: 256 GB RAM, 4 h wall time.
