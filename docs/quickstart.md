# Quick Start Guide

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

# Generate configuration
kintsugi mcp config /path/to/your/project
```

Add the output to `.claude/settings.local.json`:

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

Create a `data` folder in the KINTSUGI directory and move your image data there:

```
KINTSUGI/
└── data/
    └── [your image files]
```

## Next Steps

- See the [Workflows](workflows.md) guide for detailed processing pipelines
- Check [Troubleshooting](TROUBLESHOOTING.md) if you encounter issues
- Review the [API Reference](api.md) for programmatic usage
