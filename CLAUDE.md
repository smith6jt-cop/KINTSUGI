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

The `/meta/experiment.json` file stores microscope and acquisition parameters. It is **auto-generated** during project creation with sensible defaults and auto-detected values (cycles, z-planes).

**Example experiment.json:**
```json
{
  "tile_rows": 5,
  "tile_cols": 5,
  "tile_overlap": 0.3,
  "xy_pixel_size": 377.0,
  "z_step_size": 1500.0,
  "numerical_aperture": 0.75,
  "tissue_refractive_index": 1.44,
  "wavelengths": {
    "1": [358.0, 461.0],
    "2": [753.0, 775.0],
    "3": [560.0, 575.0],
    "4": [648.0, 668.0]
  },
  "channels_per_cycle": 4,
  "n_cycles": 7,
  "n_zplanes": 15
}
```

**CLI options** to customize during init:
```bash
kintsugi init /path/to/project --name "My Experiment" \
    --tile-rows 5 --tile-cols 5 \
    --xy-pixel-size 377 --z-step-size 1500 \
    --numerical-aperture 0.75 --tissue-ri 1.44
```

**CODEX/Akoya Format Compatibility:**
`ExperimentConfig.from_dict()` automatically translates CODEX field names to KINTSUGI equivalents:

| CODEX Field | KINTSUGI Field |
|-------------|---------------|
| `numCycles` | `n_cycles` |
| `numZPlanes` | `n_zplanes` |
| `regionHeight` | `tile_rows` |
| `regionWidth` | `tile_cols` |
| `numChannels` | `channels_per_cycle` |
| `xyResolution` | `xy_pixel_size` |
| `zPitch` | `z_step_size` |
| `aperture` | `numerical_aperture` |

CODEX flat wavelength lists (e.g., `[358, 488, 550, 650]`) are mapped to proper (excitation, emission) pairs using `_CODEX_FILTER_SETS` in `project.py`. The mapping uses known filter sets:

| Excitation | Filter | Emission |
|------------|--------|----------|
| 358 nm | DAPI | 461 nm |
| 488 nm | FITC/Alexa 488 | 525 nm |
| 550 nm | TRITC/Cy3 | 575 nm |
| 650 nm | Cy5 | 668 nm |
| 750 nm | Cy7 | 775 nm |

CODEX experiment.json does NOT include refractive index. For uncleared tissue (all HuBMAP CODEX datasets), RI is always **1.44**.

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

When running `kintsugi init` on a directory with existing data, the command intelligently handles different scenarios:

**Raw data only** (in `data/raw/`):
- Shows data summary with cycle count
- Options: Continue or Cancel
- Raw data stays in place; project structure is created around it

**Processed data exists** (in `data/processed/`):
- Shows breakdown by processing stage (stitched, deconvolved, edf, etc.)
- Options: Delete processed data, Keep processed data, or Cancel
- If deleting, can select which stages to remove

**Existing project** (has `kintsugi_project.json`):
- If `--slurm` requested but not configured: Automatically adds SLURM
- Otherwise: Shows status and suggests `--force` to refresh templates

**Scanning data**:
```bash
kintsugi scan /path/to/directory   # Preview what init will find
```

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

**SLURM output naming**: SLURM job scripts produce output that matches the notebook conventions:
- **Deconvolution**: `data/processed/deconvolved/cyc##/CH#/*.tif` (z-plane TIFFs per channel)
- **EDF**: `data/processed/edf/cyc##/{channel_name}.tif` (marker-named, e.g., `CD3.tif`, `DAPI-01.tif`)

EDF output files use marker names from `CHANNELNAMES.txt` (loaded via `Kio.load_channel_names()`), falling back to `CH#` if the file is missing. This is critical for downstream compatibility with `_find_edf_file()` in `Kview_qc.py`.

**PYTHONPATH in SLURM jobs**: Job scripts must include `KINTSUGI_DIR/notebooks` on `sys.path` to import notebook modules like `Kio`. Project notebooks directories may only have a subset of synced files. The required path setup is:
```python
sys.path.insert(0, str(PROJECT_DIR / 'notebooks'))
sys.path.insert(0, str(KINTSUGI_DIR / 'notebooks'))  # Required for Kio, Kprocess, etc.
sys.path.insert(0, str(KINTSUGI_DIR))
```

**Memory allocation**: GPU jobs use CuPy (float32 in GPU memory), so 48 GB CPU RAM is sufficient. CPU jobs use SciPy (float64 in system memory) and need more. Snakefile lambdas route `mem_decon`/`cpu_mem_decon` automatically based on GPU vs CPU mode.
```bash
# GPU jobs (CuPy does heavy lifting in GPU memory)
export MEM_DECON=48
export MEM_EDF=48
# CPU jobs (float64 in system memory, needs ~2.5x more)
export CPU_MEM_DECON=128
export CPU_MEM_EDF=96
```

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

### Processing Modes: Notebooks vs SLURM

KINTSUGI has two distinct processing modes with different resource utilization strategies:

| Mode | Context | GPU Policy | CPU Policy | Rationale |
|------|---------|------------|------------|-----------|
| **Notebook** | Interactive, user watching | GPU required, no fallback | Not used | Quality-first; user can intervene if GPU unavailable |
| **SLURM** | Headless batch processing | GPU + CPU concurrent | CPU concurrent | Maximize throughput; use ALL available resources |

**Notebook Processing (Interactive)**:
- GPU is **required** - fails loudly if unavailable
- User is present to diagnose and fix GPU issues
- Quality parameters are non-negotiable (e.g., BaSiC iterations=500)
- See `gpu-quality-priority` skill for details

**SLURM Processing (Headless)**:
- Maximizes throughput by running GPU and CPU jobs **concurrently**
- Uses a **dual-pool architecture** to calculate total concurrent slots from both GPU and CPU resources
- Jobs automatically adapt via `KINTSUGI_DEVICE_MODE` environment variable
- Quality parameters remain unchanged - only the compute device differs

### Resource Pool Calculation (Multi-Account Architecture)

KINTSUGI calculates total concurrent jobs from **all** non-blocked SLURM accounts. Each account contributes **both** GPU slots and CPU slots independently:

- **GPU slots** per account = QOS `gres/gpu` limit (via `sacctmgr`)
- **CPU slots** per account = `floor(0.85 * account_cpus / cpus_per_job)` (85% cap avoids starving other users)
- GPU and CPU partitions are separate resource pools — GPU jobs on `hpg-b200` do NOT consume CPU allocation on `hpg-default`
- **`brusko` account is permanently blocked** — hard-coded in `BLOCKED_ACCOUNTS` frozenset in `hpc.py`

**Total Concurrent** = sum(GPU slots) + sum(CPU slots) across all accounts

**Example calculation** (two accounts, each with GPUs AND CPUs):
| Account | CPUs | Memory | GPUs | GPU Slots | CPU Slots | Calculation |
|---------|------|--------|------|-----------|-----------|-------------|
| `clive` | 104 | 812 GB | 3 | 3 | 11 | GPUs: 3/1; CPUs: floor(0.85*104/8) |
| `maigan` | 80 | 625 GB | 2 | 2 | 8 | GPUs: 2/1; CPUs: floor(0.85*80/8) |
| **Total** | | | **5** | **5** | **19** | **24 concurrent jobs** |

### How Concurrent GPU/CPU Works

1. `detect_multi_account_resources()` / `detect_live_multi_account()` query `sacctmgr show associations` for all user accounts (filtering burst `-b` accounts and blocked accounts)
2. Each account's GPU and CPU slots are calculated independently from its QOS limits
3. GPU jobs are submitted on each account's GPU partition (e.g., `hpg-b200`) with `KINTSUGI_DEVICE_MODE=gpu`
4. CPU jobs are submitted on each account's CPU partition (e.g., `hpg-default`) with `KINTSUGI_DEVICE_MODE=cpu` and guaranteed (non-preemptible) resources
5. Job scripts read `KINTSUGI_DEVICE_MODE` and use appropriate backend (CuPy for GPU, NumPy/SciPy for CPU)
6. CPU jobs get a 5x time multiplier (automatic) but have guaranteed memory allocation

**Important**: `sacctmgr show associations user=USERNAME format=account -n -P` is the correct query format. The older `sacctmgr show user USERNAME format=account` returns empty results.

### Snakemake Workflow (Alternative to submit.sh)

KINTSUGI includes a Snakemake-based workflow as an alternative to `submit.sh` for SLURM batch processing. Both produce files in the same `data/processed/` tree but should not run simultaneously on the same project.

**Setup and usage:**
```bash
kintsugi workflow config .    # Generate workflow/config.yaml + Snakefile
kintsugi workflow check .     # Show live per-account availability
kintsugi workflow run .       # Submit via Snakemake (auto-sets -j)
```

**Key design decisions:**

| Design Choice | Rationale |
|---------------|-----------|
| 3 rules with lambda resources | `ruleorder` does NOT fall back — it always picks the preferred rule, so 6 rules + ruleorder was fundamentally broken |
| Cycle pre-assignment at DAG creation | `_build_cycle_assignment()` distributes cycles proportionally across accounts/modes before any job runs |
| Sentinel files (`.snakemake_complete`) | Stitching/deconvolution produce hundreds of files; declaring every output would create an enormous DAG |
| Per-cycle dependencies | `stitch cyc01 → decon cyc01 → edf cyc01` allows pipelining across cycles |
| No `--resources gpus=N` needed | GPU budget is baked into cycle pre-assignment |

**Config format** (`workflow/config.yaml`):
```yaml
resources:
  accounts:
    - name: clive
      partition_gpu: "hpg-b200,hpg-turin"
      partition_cpu: hpg-default
      gpu_slots: 3
      cpu_slots: 11
    - name: maigan
      partition_gpu: "hpg-b200,hpg-turin"
      partition_cpu: hpg-default
      gpu_slots: 2
      cpu_slots: 8
  total_slots: 24
  cpu_time_multiplier: 5
  cpu_cpus_per_task: 8
```

**Cycle assignment example** (13 cycles, clive 3G+11C, maigan 2G+8C):
- Cycles 1-3 → clive GPU, Cycles 4-5 → maigan GPU
- Cycles 6-16 → clive CPU, Cycles 17+ → maigan CPU (round-robin overflow)

**Lambda resource routing** (per-rule in Snakefile):
```python
resources:
    slurm_partition=lambda wc: _assign(wc)["partition"],
    slurm_account=lambda wc: _assign(wc)["account"],
    gres=lambda wc: "gpu:1" if _is_gpu(wc) else "",
    runtime=lambda wc: base_time * (1 if _is_gpu(wc) else CPU_TIME_MULT),
```

**GPU resource naming**: Use `gres="gpu:1"` (maps to `--gres=gpu:1`). Do NOT use `gpus=1` (maps to `--gpus=1`) or `gpu=1` — both trigger `SLURM_TRES_PER_TASK` conflicts on SLURM >= 24.11.

**SLURM profile** (`workflow/profiles/slurm/config.yaml`):
- `precommand: "module load conda && conda activate KINTSUGI"` — activates the existing conda environment on compute nodes so that conda-installed packages (cupy, torch, etc.) are importable inside srun jobs

**SLURM_TRES_PER_TASK fix** (critical for SLURM >= 24.11):
- SLURM 24.11+ sets `SLURM_TRES_PER_TASK` in GPU job environments, which conflicts with the jobstep plugin's `srun` call
- Error: `srun: fatal: SLURM_TRES_PER_TASK is mutually exclusive with --ntasks-per-gpu`
- Fix: patched `snakemake_executor_plugin_slurm_jobstep/__init__.py` to unset this variable in `__post_init__()`: `os.environ.pop("SLURM_TRES_PER_TASK", None)`
- **This patch must be re-applied after any pip upgrade of `snakemake-executor-plugin-slurm-jobstep`**

**Cycle directory resolution**: `_resolve_raw_cycle_dir()` handles `cyc001_reg001_*`, `cyc001`, `cyc01`, and `Cyc01` naming conventions at DAG creation time. Accepts int or str (Snakemake CLI `--config` passes strings).

**`workflow config` behavior**: Always overwrites the Snakefile and profiles (so pipeline logic and SLURM precommand updates propagate); only copies scripts if they don't already exist.

**What Snakemake replaces vs keeps:**
- Replaces: `submit.sh` orchestration, dependency wiring, `.complete` polling, skip-existing logic, `--array` limit calculation
- Keeps: Python processing code in wrapper scripts (`workflow/scripts/*.py`), `KINTSUGI_DEVICE_MODE` environment variable pattern
- Coexists with: Existing `slurm/` job scripts (run one system or the other, not both)

## Development Workspace

Use the VS Code multi-root workspace (`kintsugi-dev.code-workspace`) which combines:

| Workspace Folder | Path | Purpose |
|------------------|------|---------|
| KINTSUGI (repo) | `KINTSUGI/` | Source code |
| Test Project | `KINTSUGI/test_data/mini_project/` | 2x2 test dataset |
| Full Project | `KINTSUGI_Projects/.../1904CC1-1L/` | Real data |

Notebooks default to mini_project for testing. Change `PROJECT_DIR` for other projects.

### Automatic Project Sync

**IMPORTANT**: Changes to notebook modules in the main repo are **automatically synced** to all project folders after every commit via a git post-commit hook.

**Always edit the main repo first** (`KINTSUGI/notebooks/`), never edit project folders directly. The sync happens automatically via:
- `scripts/sync_to_projects.py` - Sync script (uses MD5 checksum comparison)
- `.git/hooks/post-commit` - Git hook that runs sync after each commit

**Sync uses checksum comparison** (not timestamps) to detect changes. This ensures:
- Notebooks saved with output in project folders don't block updates
- Content changes are always detected regardless of file timestamps
- Network file system timestamp issues don't cause sync failures

**What gets synced:**
- Notebook modules: `Kdecon/`, `Kstitch/`, `Kreg/`, `Kview/`, `Kview2/`, `Kseg/`
- Python files: `Kio.py`, `Kprocess.py`, `Kutils.py`, `Kview_qc.py`, `Kpipeline.py`, `Kvis.py`
- Notebooks: `1_Single_Channel_Eval.ipynb`, `2_Cycle_Processing.ipynb`, etc.

**Project folders synced to:**
- Test project: `KINTSUGI/test_data/mini_project/notebooks/`
- Full project: `KINTSUGI_Projects/CODEX_SP_LN/1904CC1-1L/notebooks/`

**Manual sync** (if needed):
```bash
python scripts/sync_to_projects.py           # Sync all projects
python scripts/sync_to_projects.py --dry-run # Preview changes
python scripts/sync_to_projects.py --verbose # Show detailed output
python scripts/sync_to_projects.py --force   # Force sync all files
```

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
- `Kprocess.py` - GPU-accelerated QC statistics and visualization: `compute_*_stats_gpu()`, `collect_*_stats_parallel()`, `plot_summary_heatmaps()`, `run_*_qc()` convenience functions
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

### Troubleshooting "Function Not Found" Errors

When users report `NameError: name 'function_name' is not defined`:

1. **Check if function is defined in a notebook cell** (not a module):
   - Some functions like `run_deconvolution` are defined inside the notebook, not imported from modules
   - These require running the definition cell before the cell that uses them

2. **Verify cell execution order**:
   ```python
   # Use this to analyze notebook cell order
   import json
   with open('notebook.ipynb') as f:
       nb = json.load(f)
   for i, cell in enumerate(nb['cells']):
       source = ''.join(cell.get('source', []))
       if 'function_name' in source:
           is_def = 'def function_name' in source
           print(f"Cell {i}: {'DEFINES' if is_def else 'CALLS'} function_name")
   ```

3. **Common patterns in 2_Cycle_Processing.ipynb**:
   | Function | Defined In | Called In |
   |----------|------------|-----------|
   | `run_deconvolution` | Cell 24 | Cell 25 |
   | `process_edf_tiff` | Cell 31 | Cell 32 |
   | `visualize_deconvolution` | Cell 27 | (manual) |

4. **Solution**: Run cells sequentially from the top, or at minimum run the definition cell before the calling cell.

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

**Deconvolution Edge Apodization (Critical Fix - Jan 2026):**
- Fixed severe horizontal banding (Gibbs ringing) artifacts in deconvolution output
- Root cause: `_create_tukey_window_3d` function existed but was NEVER CALLED
- The Tukey window provides cosine-tapered edge apodization to prevent FFT artifacts
- Fix: Added `tukey_window = _create_tukey_window_3d(...)` call before deconvolution
- Investigation showed raw and stitched images were NORMAL; artifacts introduced during deconvolution
- Affected code: `notebooks/Kdecon/deconvolution.py` lines 509-512

**Deconvolution Intensity Scaling (Critical Fix):**
- Output scaling now preserves original intensity relationships between channels
- Previously, histogram clipping always stretched output to full 16-bit range (0-65535)
- This caused artificial saturation and inverted relative intensities between channels
- Fix: Scale to `raw_max` (input maximum) instead of 65535
- Critical for quantitative marker expression analysis where relative intensities matter
- Affected code: `notebooks/Kdecon/main.py` lines 283-288

**BaSiC Illumination Correction (Critical Fixes - Jan 2026):**
- **Negative flatfield fix**: Algorithm can produce negative flatfield values for low-signal images
  (sparse markers, blank channels). Now falls back to uniform flatfield (no correction) when detected.
- **Flatfield minimum threshold**: `flatfield_min` parameter in `transform()` (default 0.1) prevents
  extreme amplification by clamping flatfield values
- **Root cause**: Notebook logs showed "WARNING: 94.0% of flatfield clamped (min=-0.4167)" -
  negative flatfield values were causing severe artifacts
- Affected code: `src/kintsugi/kcorrect_gpu.py` lines 567-580 (fit), 618-630 (transform)

**Stitching Stripe Artifacts (Critical Finding - Jan 2026):**
- Vertical stripe patterns (~30px spacing) can appear in stitched images
- **IMPORTANT**: Stripes are NOT caused by raw data quality - investigation confirmed raw tiles are clean
- Root cause: Stripes introduced during stitching by buggy pipeline code
- Fix: Reprocess affected z-planes with current code using `scripts/reprocess_striped_zplanes.py`
- Typical improvement: 2.8-5.5x reduction in stripe metrics after reprocessing
- Diagnostic scripts: `scripts/check_raw_stripe_pattern.py`, `scripts/check_correction_vs_stitching.py`
- When investigating stripe artifacts, always compare fresh reprocessing vs existing files first

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

**EDF Sigma Fix (Critical Fix - Jan 2026):**
- Fixed significant detail loss caused by incorrect sigma parameter application
- Root cause: sigma was smoothing the INPUT IMAGE before variance calculation (destroying detail)
- Fix: sigma now smooths the VARIANCE MAP after variance calculation (matching CLIJ2 behavior)
- This preserves high-frequency detail while regularizing focus selection
- Default parameters updated to match CLIJ2: radius=2, sigma=10 (was radius=5, sigma=20)
- The "altitude map" (variance surface) is smoothed, NOT the image content
- Affected code: `src/kintsugi/edf.py` lines 284-296

**EDF Smooth Transitions:**
- `blend_depth` parameter - Number of adjacent z-slices to blend for smooth transitions
- `z_smooth_sigma` parameter - Gaussian smoothing for z-index map
- Fixes abrupt transitions in areas of changing contrast
- Usage:
```python
from kintsugi.edf import extended_depth_of_focus_variance

edf = extended_depth_of_focus_variance(
    stack,
    radius_x=2,         # CLIJ2 default (was 5)
    radius_y=2,         # CLIJ2 default (was 5)
    sigma=10.0,         # CLIJ2 default (was 20)
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

**Quality Gate (Pre-Processing Validation):**
- `QualityGate` class detects problematic images before deconvolution/EDF
- Useful for detecting: stripe artifacts, corruption patterns, saturation, low signal, tile grids
- **Note**: Most "artifacts" previously attributed to raw data were actually caused by the
  deconvolution edge apodization bug (now fixed, see above). When investigating image
  quality issues, always verify raw data is actually bad before attributing to data quality.
- Usage:
```python
from kintsugi.qc import QualityGate, check_before_processing

# Create gate and check z-stack
gate = QualityGate()
result = gate.check_zstack(z_stack, channel="CD3", cycle="cyc01")

if result.passed:
    deconvolved = decon.process_array(z_stack)
else:
    print(f"Quality gate failed: {result.issues}")
    if result.can_mitigate:
        mitigated = gate.mitigate(z_stack, result)
        deconvolved = decon.process_array(mitigated)

# Or use convenience function with automatic handling
stack, result = check_before_processing(
    z_stack, channel="CD3", cycle="cyc01",
    fail_action="mitigate"  # Options: 'raise', 'warn', 'skip', 'mitigate'
)
```

## Batch Processing (Multi-Dataset)

KINTSUGI supports automated batch processing of multiple CODEX datasets. See `KINTSUGI_Batch_Processing_Guide.docx` for the complete 8-phase workflow.

**Key scripts** (in `/blue/maigan/smith6jt/`):
- `data_discovery.py` - Scans orange storage, builds `dataset_manifest.csv` from experiment.json files
- `setup_all_projects.sh` - Creates KINTSUGI projects for all datasets in manifest
- `stage_datasets.sh` - Submits SLURM rsync jobs to copy raw data from orange to blue
- `configure_workflows.sh` - Generates Snakemake configs for each project

**Storage architecture** (HiPerGator):
- `/orange/maigan/` - Long-term storage for raw CODEX data (~9.7 TB, 34 datasets)
- `/blue/maigan/` - Processing workspace (~17 TB available). Raw data must be copied here.
- Peak blue usage per dataset: ~3x raw size (raw + stitched + deconvolved intermediates)

**Rsync staging pattern** (preserves cycle directory structure):
```bash
rsync -av '${orange_path}/' '${RAW_DIR}/' \
  --include='*/' --include='*.tif' --exclude='*'
```
**CRITICAL**: Never use bash glob `*/` on rsync source — it flattens cycle directories and overwrites files across cycles. Use the source path with trailing `/` and let rsync handle the tree.

**Batch processing phases**:
1. Discovery → 2. Setup + Staging → 3. Single-cycle validation → 4. SLURM batch (stitch→decon→EDF) → 5. Registration → 6. Cleanup intermediates → 7. Signal isolation (local GPU PCs) → 8. Segmentation (local GPU PCs)

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

### CuPy Installation Status (IMPORTANT — READ BEFORE MODIFYING)

**CuPy IS ALREADY INSTALLED** in the KINTSUGI conda environment. **DO NOT attempt to install, reinstall, or upgrade CuPy.** It is installed as `cupy-cuda12x` and is the correct version for the HiPerGator CUDA 12 stack.

**Why CuPy may appear "unavailable":**
- On **login nodes** (no GPU hardware), `import cupy` succeeds but GPU operations fail with `cudaErrorInsufficientDriver`
- Functions like `_check_cupy_available()`, `gpu.cupy_available`, and `check_gpu()` test **GPU hardware availability**, not package installation
- These returning `False` means "no GPU on this node," NOT "CuPy needs to be installed"
- Use `gpu.cupy_installed` (import-only check) to verify the package is actually installed

**If CuPy genuinely needs reinstalling** (extremely rare — only after environment rebuild):
```bash
kintsugi install gpu    # Installs cupy-cuda12x into the active conda env
```

### HPC Cache Redirection (IMPORTANT)

All SLURM/Snakemake jobs MUST use account-specific cache directories, NOT the home directory. The system uses:
- `~/.use_conda_maigan.sh` / `~/.use_conda_clive.sh` — Account switchers that set `BLUE_BASE` and source cache redirection
- `~/.cache_redirect.sh` — Redirects pip, torch, numba, jupyter, XDG caches to `/blue/{account}/smith6jt/scratch/cache/`
- The Snakemake profile precommand automatically sources the correct switcher based on `SLURM_JOB_ACCOUNT`

**Never install packages or modify caches directly.** Always use the conda environment and account-specific paths.

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
