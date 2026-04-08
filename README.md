# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.

Citation Information:

Smith, J. A. et al. Protocol for processing and analyzing multiplexed images improves lymphatic cell identification and spatial architecture in human tissue. STAR Protocols 6, 103976 (2025).

[![CI](https://github.com/smith6jt-cop/KINTSUGI/actions/workflows/ci.yml/badge.svg)](https://github.com/smith6jt-cop/KINTSUGI/actions/workflows/ci.yml)
[![Documentation Status](https://readthedocs.org/projects/kintsugi/badge/?version=latest)](https://kintsugi.readthedocs.io/en/latest/?badge=latest)

## Table of Contents

- [Documentation](#documentation)
- [Installation](#installation)
  - [Linux](#linux)
  - [Windows](#windows)
  - [macOS](#macos)
  - [Optional Features](#optional-features)
- [Updating Existing Projects](#updating-existing-projects)
- [Verify Installation](#verify-installation)
- [Usage](#usage)
  - [Command Line Interface](#command-line-interface)
  - [Python API](#python-api)
- [Notebooks](#notebooks)
- [HPC/SLURM Job Submission](#hpcslurm-job-submission)
  - [Prerequisites](#prerequisites-hpc)
  - [Snakemake Workflow (Recommended)](#snakemake-workflow-recommended)
  - [Progress Dashboard](#progress-dashboard)
  - [Monitoring and Logs](#monitoring-and-logs)
  - [Reprocessing and Recovery](#reprocessing-and-recovery)
  - [Batch Processing Multiple Projects](#batch-processing-multiple-projects)
  - [Legacy SLURM Submission](#legacy-slurm-submission-submitsh)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

## Documentation

**Full documentation with tutorials is available at [kintsugi.readthedocs.io](https://kintsugi.readthedocs.io/)**

The documentation includes:
- [Installation Guide](https://kintsugi.readthedocs.io/en/latest/installation.html) - Detailed setup instructions for all platforms
- [Quick Start](https://kintsugi.readthedocs.io/en/latest/quickstart.html) - Get up and running quickly
- [Processing Workflows](https://kintsugi.readthedocs.io/en/latest/workflows.html) - Step-by-step workflow guides
- [CLI Reference](https://kintsugi.readthedocs.io/en/latest/cli.html) - Command line interface documentation
- [API Reference](https://kintsugi.readthedocs.io/en/latest/api.html) - Python API documentation
- [Troubleshooting](https://kintsugi.readthedocs.io/en/latest/TROUBLESHOOTING.html) - Common issues and solutions

## Installation

KINTSUGI uses a streamlined base installation with optional feature groups that can be added as needed. This ensures fast environment creation and avoids dependency conflicts.

### Linux

```bash
# 1. Clone the repository
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# 2. Create the base conda environment
conda env create -f envs/env-linux.yml

# 3. Activate and verify
conda activate KINTSUGI
kintsugi check
```

### Windows

```powershell
# 1. Clone the repository
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# 2. Create the base conda environment
conda env create -f envs/env-windows.yml

# 3. Activate
conda activate KINTSUGI

# 4. Download and install libvips (REQUIRED for Windows)
#    Download PyVips-dev from: https://zenodo.org/records/14969214
#    Extract to the KINTSUGI folder

# 5. Verify installation
kintsugi check
```

### macOS

```bash
# 1. Install libvips (required)
brew install vips

# 2. Clone the repository
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# 3. Create the base conda environment
conda env create -f envs/env-macos.yml

# 4. Activate and verify
conda activate KINTSUGI
kintsugi check
```

### Optional Features

After installing the base environment, add optional features as needed. **Install `gpu` first** — other groups that need PyTorch depend on it for CUDA-enabled builds.

> **Important:** Always use `kintsugi install` commands — never run `pip install torch` directly, as it installs a CPU-only build that breaks GPU processing.

```bash
# Step 1: GPU acceleration (install this FIRST — provides CUDA torch + CuPy)
kintsugi install gpu

# Step 2 (HPC only): Snakemake workflow orchestration
kintsugi install workflow
kintsugi patch slurm              # Required for SLURM >= 24.11

# Step 3: Additional features as needed
kintsugi install dl               # Deep learning segmentation (InstanSeg)
kintsugi install analysis         # Spatial analysis (scanpy, scimap)
kintsugi install viz              # Napari interactive visualization
kintsugi install bio              # Bio formats I/O (OME-TIFF, LIF)
kintsugi install claude           # Claude Code MCP integration
kintsugi install dev              # Development tools (pytest, ruff, black)

# Or install everything at once (resolves all constraints together)
kintsugi install all

# Step 4: Verify everything is correct
kintsugi check --strict
```

All `kintsugi install` commands automatically enforce a `constraints.txt` file that prevents known-bad version combinations (e.g., numpy 2.x). See [docs/DEPENDENCY_GUIDE.md](docs/DEPENDENCY_GUIDE.md) for troubleshooting.

### Multi-GPU Acceleration

KINTSUGI automatically detects and uses all available NVIDIA GPUs for image processing:

```python
from kintsugi.gpu import get_gpu_manager
gpu = get_gpu_manager()
print(gpu.summary())  # Shows all detected GPUs
```

**GPU-Accelerated Operations:**
- **Illumination Correction** - BaSiC algorithm with CuPy FFT
- **Stitching** - Phase correlation on GPU
- **Deconvolution** - Lucy-Richardson with GPU FFT
- **Extended Depth of Focus** - Variance projection on GPU

```python
# Example: GPU-accelerated illumination correction
from kintsugi.kcorrect_gpu import KCorrectGPU
corrector = KCorrectGPU(device_id=0)  # Specify GPU
flatfield, darkfield = corrector.estimate(images)

# Example: Multi-GPU stitching
from notebooks.Kstitch._translation_computation import get_multi_gpu_accelerator
accelerator = get_multi_gpu_accelerator(tile_shape)
```

> **Note:** For GPU-accelerated single-cell analysis (clustering, UMAP, spatial analysis),
> see the [rapids_singlecell](https://github.com/smith6jt-cop/rapids_singlecell) repository.

**Feature Requirements by Notebook:**

| Notebook | Required Features |
|----------|-------------------|
| 1_Single_Channel_Eval | `gpu` |
| 2_Cycle_Processing | `gpu` |
| 3_Signal_Isolation_QC | `claude` (optional) |
| 4_Segmentation_Analysis | `dl`, `viz`, `analysis` |

Each notebook will check for required dependencies at startup and provide installation instructions if anything is missing.

## Updating Existing Projects

When KINTSUGI is updated (e.g., `git pull`), existing project directories may contain outdated copies of notebooks, modules, workflow scripts, or Snakefiles. Use these commands to bring them up to date.

### Notebooks and Python Modules

A post-commit hook runs `scripts/sync_to_projects.py` automatically after every `git commit`, syncing notebooks and Python modules (Kreg/, Kstitch/, Kview2/, Kio.py, etc.) to all discovered project directories. The sync uses MD5 checksums — only changed files are copied.

To sync manually (e.g., after `git pull`):

```bash
python scripts/sync_to_projects.py           # Sync all projects
python scripts/sync_to_projects.py --dry-run  # Preview changes
python scripts/sync_to_projects.py --force    # Force overwrite all
```

> **Note:** Autoreload (`%autoreload 2`) is enabled in all KINTSUGI notebooks — you do **not** need to restart kernels after module updates. Just re-run the relevant cell.

### Workflow Scripts and Snakefile

Re-run `workflow config` to update the Snakefile, SLURM profiles, and add any new workflow scripts:

```bash
kintsugi workflow config /path/to/project
```

| File | Update behavior |
|------|-----------------|
| **Snakefile** | Always overwritten (pipeline logic must stay current) |
| **profiles/** | Always overwritten (SLURM settings, precommand) |
| **scripts/*.py** | Only added if missing (existing scripts are preserved) |

**Updating existing workflow scripts:** If a workflow script has been updated in the repo (e.g., `registration.py`, `qc_report.py`), `workflow config` will not overwrite the project's copy. To update, delete the stale script and re-run config:

```bash
# Update a single script
rm /path/to/project/workflow/scripts/registration.py
kintsugi workflow config /path/to/project

# Or update all scripts at once
rm /path/to/project/workflow/scripts/*.py
kintsugi workflow config /path/to/project
```

For bulk updates across many projects:

```bash
# Copy all current scripts to every configured project
for proj_scripts in /path/to/KINTSUGI_Projects/*/workflow/scripts/; do
    [ -d "$proj_scripts" ] || continue
    cp workflow/scripts/*.py "$proj_scripts/"
done
```

## Verify Installation

After installation, verify everything is working:

```bash
# Check all dependencies
kintsugi check

# Show version info
kintsugi info

# Generate config template
kintsugi template -o my_config.json
```

## External Dependencies

| Dependency | Purpose | Installation |
|------------|---------|--------------|
| **libvips** | High-performance image I/O | `conda install libvips` (Linux), `brew install vips` (macOS), or Zenodo (Windows) |
| **VALIS** | Image registration | Included in base install |
| **CuPy** | GPU acceleration (optional) | `kintsugi install gpu` |
| **PyTorch** | Deep learning (optional) | `kintsugi install torch` |
| **RAPIDS** | GPU data science (optional) | See [installation docs](docs/installation.md#rapids-gpu-accelerated-data-science-optional) |

> **Note:** Java, Maven, and FIJI/CLIJ2 are no longer required. KINTSUGI now uses pure Python
> implementations (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF).

## Usage

### Command Line Interface

KINTSUGI provides a CLI for common operations:

```bash
# Check dependencies and show system info
kintsugi check
kintsugi info

# Install optional features
kintsugi install gpu
kintsugi install viz --conda  # Use conda where available

# Generate configuration template
kintsugi template -o config.json

# Run registration workflow
kintsugi register config.json --dry-run
kintsugi register config.json

# Initialize a new project
kintsugi init /path/to/project --name "My Project"
kintsugi scan /path/to/directory   # Preview what init will find

# Snakemake workflow management (recommended for HPC)
kintsugi workflow config .         # Auto-detect accounts, generate config + Snakefile
kintsugi workflow check .          # Show live per-account slot availability
kintsugi workflow run .            # Submit via Snakemake with auto-calculated -j

# Legacy SLURM submission
kintsugi slurm init /path/to/project
kintsugi slurm submit /path/to/project
kintsugi slurm status /path/to/project
kintsugi slurm cancel /path/to/project
```

### Python API

```python
import kintsugi

# Check dependencies
kintsugi.check_dependencies()

# Get configuration template
config = kintsugi.get_config_template()

# Access modules
from kintsugi import Kreg, Kview2, Kstitch

# Registration
from kintsugi.kreg import Valis
registrar = Valis(
    src_dir="/path/to/images",
    dst_dir="/path/to/output",
    reference_img_f="cycle1.tif",
)
registrar.register()

# Visualization
from kintsugi.kview2 import imshow, curtain, crop

# Quality Control
from kintsugi.qc import ImageQC, CellQC, MarkerQC
qc = ImageQC()
result = qc.assess(image)

# Advanced Denoising
from kintsugi.denoise import adaptive_denoise, denoise_n2v
denoised = adaptive_denoise(image, strength="auto")
```

### Claude Code Integration

KINTSUGI includes an MCP (Model Context Protocol) server that enables Claude Code to act as an AI-powered image processing assistant.

**Setup:**

```bash
# Install Claude Code dependencies
pip install kintsugi[claude]
```

**If creating a new project:** Use `kintsugi init` - Claude Code configuration is created automatically.

**If adding to an existing project:**
```bash
kintsugi mcp config /path/to/your/project
```

**Available Tools:**

| Category | Tools |
|----------|-------|
| Signal Isolation | `load_channel`, `subtract_blank`, `denoise`, `denoise_advanced`, `apply_clahe`, `clean_background` |
| Quality Assessment | `assess_quality`, `compute_snr` |
| Workflow | `list_channels`, `save_processed`, `suggest_parameters` |
| Parameter Learning | `get_learned_parameters`, `approve_and_learn`, `suggest_with_learning` |

**Parameter Learning:**

The system learns from successful parameter choices:
- Parameters are stored in SQLite databases indexed by tissue type and marker name
- Future recommendations are weighted by past success
- Use `approve_and_learn` to record successful parameters

**Example Claude Interaction:**

```
User: "Load the CD3 channel and suggest denoising parameters"
Claude: [Uses load_channel and suggest_with_learning tools]
Claude: "Based on the image analysis and learned history for tonsil tissue,
         I recommend NLM denoising with patch_size=7, strength=0.15..."
```

See `notebooks/MIGRATION_GUIDE.md` for detailed migration instructions from legacy notebooks.

### Notebook Dependency Checking

Each notebook that requires optional features should include a dependency check at the top:

```python
# At the top of your notebook
from kintsugi.deps import require

# Check specific feature groups
require('gpu', 'viz')

# Or auto-detect from notebook name
require(notebook='4_Segmentation_Analysis')
```

If dependencies are missing, you'll see a clear error with installation instructions.

## Notebooks

The following Jupyter notebooks provide step-by-step workflows:

### 1. Parameter Tuning and Testing
Test illumination correction, stitching, deconvolution, and EDoF.
- [notebooks/1_Single_Channel_Eval.ipynb](notebooks/1_Single_Channel_Eval.ipynb)
- **Requires:** `gpu`

### 2. Batch Processing
Batch processing for illumination correction, stitching, deconvolution, EDoF, and registration.
- [notebooks/2_Cycle_Processing.ipynb](notebooks/2_Cycle_Processing.ipynb)
- **Requires:** `gpu`
- **QC Output:** Automatically generates PDF plots to `PROJECT_DIR/qc_plots/` including:
  - Summary heatmaps (SNR, CV, intensity by cycle/channel)
  - Z-plane profiles for quality assessment

### 3. Signal Isolation & Quality Control (NEW)
Combined signal isolation and quality assessment with Claude Code integration.
- [notebooks/3_Signal_Isolation_QC.ipynb](notebooks/3_Signal_Isolation_QC.ipynb)
- **Requires:** Base only (optional: `claude` for AI-assisted workflow)
- **Features:**
  - Claude-guided parameter selection
  - Interactive widget-based tuning
  - Integrated quality assessment
  - Parameter learning for future recommendations
  - Advanced denoising (N2V, NLM, BM3D-lite)

### 4. Segmentation Analysis
InstanSeg segmentation, feature extraction, and spatial analysis.
- [notebooks/4_Segmentation_Analysis.ipynb](notebooks/4_Segmentation_Analysis.ipynb)
- **Requires:** `dl`, `viz`, `analysis`

> **Note:** See [notebooks/MIGRATION_GUIDE.md](notebooks/MIGRATION_GUIDE.md) for migration guidance from older workflows.

### Running Notebooks

Launch VS Code from the activated environment:
```bash
conda activate KINTSUGI
code .
```

**Important**: Always launch VS Code from the activated conda environment to ensure all packages are available.

## HPC/SLURM Job Submission

KINTSUGI provides two SLURM submission systems. The **Snakemake workflow** (recommended) handles dependency management, skip-existing, and multi-account scheduling declaratively. The legacy `submit.sh` script is still available for simpler setups.

### Prerequisites (HPC)

Snakemake and its SLURM executor plugin are included in the base conda environment. Verify they are available:

```bash
conda activate KINTSUGI
snakemake --version        # Should print 8.x.x or higher
```

GPU support must be installed for HPC processing:

```bash
kintsugi install gpu
kintsugi check
```

### Snakemake Workflow (Recommended)

The Snakemake pipeline runs three processing stages per cycle, with automatic per-cycle dependencies:

| Stage | Description | GPU or CPU |
|-------|-------------|------------|
| `stitch` | BaSiC illumination correction + tile stitching | Both |
| `deconvolve` | Richardson-Lucy deconvolution | Both |
| `edf` | Extended depth of focus (variance projection) | Both |

Dependencies flow per-cycle, enabling pipelining: `stitch cyc01 -> decon cyc01 -> edf cyc01` runs in parallel with `stitch cyc02 -> decon cyc02 -> edf cyc02`.

#### Quick Start

```bash
# 1. Create project with metadata
kintsugi init /path/to/project --name "My Experiment" \
    --tile-rows 9 --tile-cols 7 \
    --xy-pixel-size 377 --z-step-size 1500 \
    --numerical-aperture 0.75 --tissue-ri 1.44

# 2. Copy raw data to data/raw/ (preserving cycle directories)
#    Supported naming: cyc001/, cyc01/, cyc001_reg001_*/, Cyc01/

# 3. Create channel names file (meta/CHANNELNAMES.txt)
#    Simple list format (one marker per line, cycles sequential):
#      DAPI-01, Blank, Blank, Blank, DAPI-02, CD31, CD8, CD45, ...
#    Or cycle-prefixed format:
#      1: DAPI, Blank, Blank, Blank
#      2: DAPI, CD31, CD8, CD45

# 4. Generate workflow config (auto-detects accounts, resources, cycles)
kintsugi workflow config /path/to/project

# 5. Check live resource availability across all accounts
kintsugi workflow check /path/to/project

# 6. Preview (always dry-run first)
kintsugi workflow run /path/to/project --dry-run

# 7. Submit (auto-calculates -j from live availability)
kintsugi workflow run /path/to/project
```

`workflow config` auto-detects SLURM accounts via `sacctmgr show associations`, calculates GPU and CPU slots per account, reads microscope parameters from `meta/experiment.json`, and generates `workflow/config.yaml` + copies the Snakefile and wrapper scripts.

#### Multi-Account Architecture

KINTSUGI distributes jobs across non-blocked SLURM accounts, each with its own GPU and CPU pools. As of Apr 8 2026, `clive` is in `BLOCKED_ACCOUNTS` (its QOS pool was throttled and is regularly saturated by other group members), so the active configuration is **maigan only**:

| Account | GPUs | GPU Slots | CPUs | CPU Slots | Calculation |
|---------|------|-----------|------|-----------|-------------|
| `maigan` | 2 | 2 | 80 | 8 | floor(0.85 * 80 / 8) |
| **Total** | | **2** | | **8** | **10 concurrent jobs** |

Cycles are **pre-assigned** to accounts at DAG creation time via `_build_cycle_assignment()` for deterministic scheduling. Each rule uses lambda resource functions to route jobs to the correct account, partition, and resource allocation. GPU jobs use `gres="gpu:1"` (not `gpus=1`, which triggers SLURM_TRES_PER_TASK conflicts on SLURM >= 24.11).

**Live-aware routing** (Apr 8 2026): `kintsugi workflow run` queries `detect_live_multi_account()` and forwards per-account `gpu_avail`/`cpu_avail`/`mem_avail_gb` to Snakemake via `--config live_accounts=<json>`. The Snakefile assignment helpers use the live data to skip accounts saturated by **other users on the same QOS investment pool**, which is the failure mode that previously caused jobs to get stuck on `QOSGrpMemLimit`. Hard-fails if every account has zero live availability instead of queueing forever.

#### Workflow Configuration (`workflow/config.yaml`)

Auto-generated by `kintsugi workflow config`. Includes cycles, channels, microscope parameters, and per-account resource allocations:

```yaml
project_dir: /path/to/my_project
kintsugi_dir: /path/to/KINTSUGI

cycles: [1, 2, 3, 4, 5, 6, 7]
channels: [1, 2, 3, 4]

tile_rows: 9
tile_cols: 7
tile_overlap: 0.3

resources:
  accounts:
    - name: maigan
      partition_gpu: "hpg-b200,hpg-turin"
      partition_cpu: hpg-default
      gpu_slots: 2
      cpu_slots: 8
  total_gpu_slots: 2
  total_cpu_slots: 8
  total_slots: 10
  cpu_time_multiplier: 5
  cpu_cpus_per_task: 8
```

#### Project Workflow Structure

```
my_project/workflow/
├── config.yaml              # Auto-generated by kintsugi workflow config
├── Snakefile                # Pipeline definition (always overwritten on config)
├── scripts/
│   ├── stitch.py            # BaSiC correction + stitching
│   ├── deconvolve.py        # Richardson-Lucy deconvolution
│   └── edf.py               # Extended depth of focus
└── profiles/slurm/
    └── config.yaml          # SLURM executor profile (precommand, retries, etc.)
```

`workflow config` always overwrites the Snakefile and SLURM profiles (so pipeline logic updates propagate). Scripts are only copied if they don't already exist — see [Updating Existing Projects](#updating-existing-projects) for how to refresh stale scripts.

#### Submitting Specific Cycles or Steps

```bash
# Specific cycles
kintsugi workflow run /path/to/project --cycles 1-3
kintsugi workflow run /path/to/project --cycles 1,3,5

# Force re-run a specific step
kintsugi workflow run /path/to/project --forcerun stitch

# Override concurrent job count
kintsugi workflow run /path/to/project -j 16

# Run locally (no SLURM, for testing)
kintsugi workflow run /path/to/project --local --cores 4

# Direct Snakemake (more control)
cd /path/to/project/workflow
snakemake --profile profiles/slurm -n              # Dry run
snakemake --profile profiles/slurm --allowed-rules deconvolve edf  # Skip stitching
snakemake --dag | dot -Tpng > dag.png              # Visualize DAG
snakemake --report report.html                     # HTML report
```

### Progress Dashboard

KINTSUGI ships a Rich-based progress dashboard that scans sentinel files, queries SLURM via `squeue`/`sacct`, parses log timings, and reports per-cycle, per-stage completion alongside live GPU/memory utilization and ETA estimates. It is exposed two ways:

```bash
# Standalone snapshot of a single project
kintsugi workflow status /path/to/project

# Live auto-refreshing dashboard (Ctrl+C to exit)
kintsugi workflow status /path/to/project --watch
kintsugi workflow status /path/to/project --watch -i 15   # Refresh every 15s

# Scan every project under a parent directory
kintsugi workflow status /path/to/KINTSUGI_Projects --all-projects --watch

# Machine-readable output for scripting
kintsugi workflow status /path/to/project --json

# Hide individual sections
kintsugi workflow status /path/to/project --no-hardware --no-estimates --no-jobs
```

Or attach the dashboard while submitting a run — Snakemake runs in the background and the dashboard refreshes in the foreground. Pressing Ctrl+C **detaches** without killing SLURM jobs:

```bash
kintsugi workflow run /path/to/project --dashboard
kintsugi workflow run /path/to/project --dashboard --dashboard-interval 15
```

**What it shows per project:**

| Section | Content |
|---------|---------|
| Cycle table | Stitch / Decon / EDF status per cycle, per-stage timing, active job & node |
| Aggregate stages | Registration and signal-isolation status with timings |
| QC sentinels | Whether `qc_stitch`, `qc_decon`, `qc_edf`, `qc_registration`, `qc_signal_isolation` reports have been generated |
| SLURM jobs | Job ID, state, elapsed time, MaxRSS, GPU id, partition, account |
| Hardware | Per-account GPU allocated/used/available, memory pools (from `squeue` + `sacctmgr`) |
| Estimate | Wall-clock ETA derived from historical log timings + current parallelism |

The standalone command lives at `kintsugi workflow status` and the embedded variant at `kintsugi workflow run --dashboard`. Both call into `src/kintsugi/dashboard.py` (`scan_project_progress`, `render_dashboard`, `watch_dashboard`).

### Monitoring and Logs

```bash
# View queued/running SLURM jobs
squeue -u $USER

# Filter to KINTSUGI jobs for a project
squeue -u $USER | grep kintsugi

# Cancel all KINTSUGI jobs for a project
kintsugi slurm cancel /path/to/project

# Per-job logs (one per cycle per step)
tail -f /path/to/project/slurm/logs/snakemake/stitch_cyc01.log
tail -f /path/to/project/slurm/logs/snakemake/decon_cyc03.log
tail -f /path/to/project/slurm/logs/snakemake/edf_cyc05.log

# Snakemake coordinator log
ls /path/to/project/workflow/.snakemake/log/
tail -f /path/to/project/workflow/.snakemake/log/*.snakemake.log

# List completed sentinel files
find /path/to/project/data/processed/ -name ".snakemake_complete"

# Count output files per stage
echo "Stitched:"; find data/processed/stitched -name "*.tif" 2>/dev/null | wc -l
echo "Deconvolved:"; find data/processed/deconvolved -name "*.tif" 2>/dev/null | wc -l
echo "EDF:"; find data/processed/edf -name "*.tif" 2>/dev/null | wc -l

# Show remaining work (dry run only shows incomplete jobs)
cd /path/to/project/workflow && snakemake --profile profiles/slurm -n
```

### Reprocessing and Recovery

**Automatic recovery:** Snakemake retries failed jobs automatically (default: 2 retries). No manual intervention needed for transient failures (OOM, node issues).

**Resume after interruption:** If the terminal disconnects or the coordinator dies, re-run the same command. Snakemake checks sentinel files and only submits jobs for incomplete cycles:

```bash
kintsugi workflow run /path/to/project
```

Per-channel skip-existing logic inside the wrapper scripts means even partially-completed cycles resume efficiently — only unfinished channels are reprocessed.

**Force reprocessing:**

```bash
# Force re-stitch a specific cycle (and everything downstream)
cd /path/to/project/workflow
snakemake --profile profiles/slurm \
    --forcerun data/processed/stitched/cyc03/.snakemake_complete

# Force re-run all deconvolution
kintsugi workflow run /path/to/project --forcerun deconvolve

# Delete outputs and reprocess from scratch
rm -rf data/processed/stitched/cyc03/ data/processed/deconvolved/cyc03/ data/processed/edf/cyc03/
kintsugi workflow run /path/to/project  # Snakemake detects missing outputs
```

### Batch Processing Multiple Projects

Process multiple staged datasets with the `kintsugi workflow batch` command:

```bash
kintsugi workflow batch /path/to/KINTSUGI_Projects               # All eligible datasets
kintsugi workflow batch /path/to/KINTSUGI_Projects --dry-run      # Preview eligible
kintsugi workflow batch /path/to/KINTSUGI_Projects -d CX_19-004   # Single dataset
kintsugi workflow batch /path/to/KINTSUGI_Projects -p 2 --detach  # Background, 2 concurrent
kintsugi workflow batch /path/to/KINTSUGI_Projects --force         # Reprocess completed
kintsugi workflow stop /path/to/KINTSUGI_Projects                  # Stop background batch
```

A dataset is eligible if it has `workflow/config.yaml` + `data/raw/.staged` and is missing `data/processed/signal_isolated/.snakemake_complete` (unless `--force`). Sequential processing is the default since all datasets share the same GPU slots.

> **Important:** Always use `kintsugi workflow batch` — do not write custom batch scripts that invoke `snakemake` directly, as they bypass GPU validation and SLURM profile detection.

### Legacy SLURM Submission (`submit.sh`)

The original `submit.sh` script is still available for single-account setups or when Snakemake is not installed:

```bash
./slurm/submit.sh --project /path/to/project
./slurm/submit.sh --project /path/to/project --dry-run
./slurm/submit.sh --project /path/to/project --steps decon,edf
./slurm/submit.sh --project /path/to/project --cycles 1-3
```

Configuration is in `slurm/config.sh` (auto-generated on first run). The legacy system uses a single-account dual-pool architecture.

### Quick Reference

```bash
# === PER-PROJECT SETUP ===
kintsugi init /path/to/project --name "My Experiment"
# Copy raw data to data/raw/, create meta/CHANNELNAMES.txt
kintsugi workflow config /path/to/project

# === RUNNING ===
kintsugi workflow run /path/to/project --dry-run    # Preview
kintsugi workflow run /path/to/project              # Full pipeline via SLURM
kintsugi workflow run /path/to/project --cycles 1-3 # Specific cycles
kintsugi workflow run /path/to/project --local      # Local (no SLURM)

# === MONITORING ===
squeue -u $USER                                     # Job queue
kintsugi workflow check /path/to/project            # Live resource availability
tail -f /path/to/project/slurm/logs/snakemake/*.log # Live logs
```

## Troubleshooting

See [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) for detailed troubleshooting guides.

### Common Issues

**libvips not found (Windows)**
```
Download PyVips-dev from Zenodo and extract to KINTSUGI folder
```

**Import errors**
```bash
# Verify environment is activated
conda activate KINTSUGI

# Check dependencies
kintsugi check
```

**GPU not detected**
```bash
# Check CUDA availability
python -c "import torch; print(torch.cuda.is_available())"

# Check all GPUs with KINTSUGI
python -c "from kintsugi.gpu import get_gpu_manager; print(get_gpu_manager().summary())"

# Install GPU support if missing
kintsugi install gpu
```

**Missing optional dependencies in notebook**
```python
# Run at the top of the notebook to see what's needed
from kintsugi.deps import require
require('gpu', 'viz', strict=False)  # Shows warning instead of error
```

**Conda environment creation hangs**

The base environment is designed to install quickly. If you experience hangs:
```bash
# Use libmamba solver (faster)
conda config --set solver libmamba

# Then retry
conda env create -f envs/env-linux.yml
```

### HPC/SLURM Issues

**"snakemake: command not found"**
```bash
conda activate KINTSUGI
pip install "snakemake>=8.0" snakemake-executor-plugin-slurm
```

**"No workflow/config.yaml found"**

Run `kintsugi workflow config .` from your project directory first.

**Jobs pending indefinitely (QOS limits)**

Check your account allocations:
```bash
sacctmgr show associations user=$(whoami) format=account,partition,qos,grptres -n -P
```
Reduce `-j` when running `kintsugi workflow run` to fit within your allocation.

**Jobs fail with OOM (Out of Memory)**

GPU jobs need ~48 GB RAM (CuPy does FFT in GPU memory). CPU jobs need ~128 GB (SciPy float64 in system memory). Adjust memory in `workflow/config.yaml` if needed:
```yaml
resources:
  mem_stitch: 48000    # MB, GPU jobs
  mem_decon: 48000
  cpu_mem_decon: 128000  # MB, CPU jobs
```

**"Missing output files after job completion" (NFS latency)**

Increase latency wait in `workflow/profiles/slurm/config.yaml`:
```yaml
latency-wait: 300  # Increased from 120 to 300 seconds
```

**"CUDA initialization failed" in job logs**

This is expected on CPU nodes — scripts automatically fall back to CPU mode. GPU nodes are selected by the SLURM partition setting in the workflow config. On login nodes, `gpu.cupy_available` returns False because there is no GPU hardware; this does **not** mean CuPy is missing.

**Stitch model not found for CH2+**

Channel 1 computes the stitching model used by all other channels. If CH1 fails, subsequent channels fail with "No stitch model." Check the CH1 stitching log first.

## Development

### Setup Development Environment

```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
conda env create -f envs/env-linux.yml
conda activate KINTSUGI
kintsugi install dev
```

### Running Tests

```bash
# Run all tests
pytest tests/ -v

# Run with coverage
pytest tests/ --cov=src/kintsugi --cov-report=html
```

### Code Quality

```bash
# Lint code
ruff check src/ tests/

# Format code
black src/ tests/
```

## Data

- **Test Data**: [KINTSUGI Zenodo Community](https://zenodo.org/communities/kintsugi)
- **Processed Results**: [Globus](https://app.globus.org/file-manager?origin_id=10f408d9-f5ee-11ef-bf21-0affeb6b961d&origin_path=%2F)

Create a `data` folder in the KINTSUGI directory and move your image data there:
```
KINTSUGI/
└── data/
    └── [your image files]
```

## License

See [License.txt](License.txt) for license information.

## Citation

```bibtex
@article{smith2025protocol,
  title={Protocol for processing and analyzing multiplexed images improves lymphatic cell identification and spatial architecture in human tissue},
  author={Smith, J. A. and others},
  journal={STAR Protocols},
  volume={6},
  pages={103976},
  year={2025}
}
```
<p align="center">
  <img src="/docs/CD3e.gif" alt="CD31 Autofluorescence Removal" style="float: right; margin-left: 20px;">
