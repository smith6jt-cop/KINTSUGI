# Quick Start Guide

## Verify Installation

After installation, verify everything is working:

```bash
# Check all dependencies
kintsugi check

# Show version info
kintsugi info
```

## Create a Project

The recommended way to start is with `kintsugi init`:

```bash
# Create a new project with standard structure
kintsugi init /path/to/my_project --name "My Experiment" \
    --tile-rows 9 --tile-cols 7 \
    --xy-pixel-size 377 --z-step-size 1500

# Preview what init will find in an existing directory
kintsugi scan /path/to/my_data
```

This creates a complete project directory:
```
my_project/
├── data/
│   ├── raw/           ← Put your raw images here (cyc001/, cyc002/, etc.)
│   └── processed/     ← Outputs go here automatically
├── meta/
│   ├── experiment.json    ← Microscope parameters (auto-generated)
│   └── CHANNELNAMES.txt   ← Channel/marker names (user-provided)
├── notebooks/         ← Working copies of processing notebooks
├── configs/           ← Processing configuration files
└── .claude/           ← Claude Code MCP config (auto-generated)
```

### Channel Names

Create `meta/CHANNELNAMES.txt` with marker names:
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

## Command Line Interface

KINTSUGI provides a CLI for common operations:

```bash
# Check dependencies
kintsugi check

# Show system info
kintsugi info

# Generate configuration template
kintsugi template -o config.json

# Run registration workflow
kintsugi register config.json --dry-run
kintsugi register config.json
```

## Python API

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
from kintsugi.qc import ImageQC
qc = ImageQC()
result = qc.assess(image)

# Denoising
from kintsugi.denoise import adaptive_denoise
denoised = adaptive_denoise(image, strength="auto")

# Segmentation
from kintsugi.segment import segment_nuclei_watershed
nuclei = segment_nuclei_watershed(dapi_image)
```

## Claude Code Integration

KINTSUGI includes an MCP server for Claude Code integration, enabling AI-assisted image processing.

### Setup

```bash
# Install Claude Code dependencies
pip install kintsugi[claude]
```

**If creating a new project:** Use `kintsugi init` - Claude Code configuration is created automatically.

**If adding to an existing project:**
```bash
kintsugi mcp config /path/to/your/project
```

### Usage

Once configured, Claude Code can:
- Load and analyze channels
- Suggest optimal processing parameters
- Apply denoising, CLAHE, and background subtraction
- Learn from successful parameters for future recommendations

Example interaction:
```
User: "Load the CD3 channel and suggest denoising parameters"
Claude: [Analyzes image and provides recommendations based on learned history]
```

See `notebooks/MIGRATION_GUIDE.md` for transitioning from legacy notebooks.

## Jupyter Notebooks

The following Jupyter notebooks provide step-by-step workflows:

### 1. Parameter Tuning and Testing
Test illumination correction, stitching, deconvolution, and EDoF.
- `notebooks/1_Single_Channel_Eval.ipynb`

### 2. Batch Processing
Batch processing for illumination correction, stitching, deconvolution, EDoF, and registration.
- `notebooks/2_Cycle_Processing.ipynb`

### 3. Signal Isolation & Quality Control (NEW)
Combined signal isolation and QC with Claude Code integration.
- `notebooks/3_Signal_Isolation_QC.ipynb`
- Features: Claude-guided workflow, parameter learning, integrated QC

### 4. Segmentation Analysis
InstanSeg segmentation, feature extraction, and spatial analysis.
- `notebooks/4_Segmentation_Analysis.ipynb`

### 5. Vessel Analysis
Specialized analysis for vessel structures.
- `notebooks/Vessel_Analysis.ipynb`

> **Note:** Old notebooks 3 and 5 are deprecated. Use `3_Signal_Isolation_QC.ipynb`.

### Running Notebooks

Launch VS Code from the activated environment:

```bash
conda activate KINTSUGI
code .
```

**Important**: Always launch VS Code from the activated conda environment to ensure all packages are available.

## Data Organization

Place raw image data in the `data/raw/` directory of your project, organized by cycle:

```
my_project/
└── data/
    └── raw/
        ├── cyc001/           ← Cycle 1 tiles
        │   ├── 1_00001_Z001_CH1.tif
        │   ├── 1_00001_Z001_CH2.tif
        │   └── ...
        ├── cyc002/           ← Cycle 2 tiles
        └── ...
```

> **Tip**: If you run `kintsugi init` on a directory that already has raw data, it will detect and organize the data automatically.

> **CODEX long-form names**: Raw folders from CODEX acquisitions (e.g. `cyc004_reg001_211206_201615`) are renamed to short form (`cyc004`) the first time you open Notebook 1 or Notebook 2 via `Kio.normalize_cycle_dirs()`. The Snakemake workflow also resolves long-form names automatically (`workflow/scripts/stitch.py:find_cycle_dir()`), so no manual rename is required either way.

## HPC Quick Start (SLURM)

For large datasets on HPC clusters, use the Snakemake workflow:

```bash
# Add SLURM/workflow support to an existing project
kintsugi init /path/to/project --slurm

# Generate Snakemake config (auto-detects accounts and resources)
kintsugi workflow config /path/to/project

# Check resource availability
kintsugi workflow check /path/to/project

# Preview the pipeline
kintsugi workflow run /path/to/project --dry-run

# Submit (run inside tmux!)
tmux new -s kintsugi
kintsugi workflow run /path/to/project
```

The pipeline runs stitch, deconvolution, and EDF per cycle, distributing jobs across all available GPU and CPU slots.

## Next Steps

- See the [Workflows](workflows.md) guide for detailed processing pipelines
- See the [CLI Reference](cli.md) for all available commands
- Check [Troubleshooting](TROUBLESHOOTING.md) if you encounter issues
- Review the [API Reference](api.md) for programmatic usage
