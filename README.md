# KINTSUGI: Knowledge Integration with New Technologies for Simplified User-Guided Image processing

<p align="center">
  <img src="/docs/CD3e.gif" alt="CD31 Autofluorescence Removal" style="float: right; margin-left: 20px;">

## Multiplex image processing for challenging datasets with a focus on user integration rather than automation.  This pipeline includes 2D/3D GPU/CPU illumination correction, stitching, deconvolution, extended depth of focus, registration, autofluorescence removal, segmentation, clustering, and spatial analysis.

Citation Information:

Smith, J. A. et al. Protocol for processing and analyzing multiplexed images improves lymphatic cell identification and spatial architecture in human tissue. STAR Protocols 6, 103976 (2025).

</p>

[![DOI](https://zenodo.org/badge/794118146.svg)](https://doi.org/10.5281/zenodo.14984518)
[![CI](https://github.com/smith6jt-cop/KINTSUGI/actions/workflows/ci.yml/badge.svg)](https://github.com/smith6jt-cop/KINTSUGI/actions/workflows/ci.yml)

<div>

## Table of Contents

- [Installation](#installation)
  - [Linux](#linux)
  - [Windows](#windows)
  - [macOS](#macos)
  - [Optional Features](#optional-features)
- [Verify Installation](#verify-installation)
- [Usage](#usage)
  - [Command Line Interface](#command-line-interface)
  - [Python API](#python-api)
- [Notebooks](#notebooks)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

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

After installing the base environment, add optional features as needed using the `kintsugi install` command:

```bash
# GPU acceleration (PyTorch + CuPy for CUDA)
kintsugi install gpu

# Napari interactive visualization
kintsugi install viz

# Deep learning segmentation (InstanSeg)
kintsugi install dl

# Spatial analysis (scanpy, scimap)
kintsugi install analysis

# Bio formats I/O (OME-TIFF, LIF, etc.)
kintsugi install bio

# All optional features
kintsugi install full
```

**Feature Requirements by Notebook:**

| Notebook | Required Features |
|----------|-------------------|
| 1_Single_Channel_Eval | `gpu` |
| 2_Cycle_Processing | `gpu` |
| 3_Signal_Isolation | None (base) |
| 4_Segmentation_Analysis | `dl`, `viz`, `analysis` |
| 5_Cluster_Analysis | `analysis` |
| Image_Registration_Workflow | None (base) |
| Vessel_Analysis | `viz`, `analysis` |

Each notebook will check for required dependencies at startup and provide installation instructions if anything is missing.

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
| **PyTorch** | Deep learning (optional) | `kintsugi install gpu` or `kintsugi install dl` |

> **Note:** Java, Maven, and FIJI/CLIJ2 are no longer required. KINTSUGI now uses pure Python
> implementations (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF).

## Usage

### Command Line Interface

KINTSUGI provides a CLI for common operations:

```bash
# Check dependencies
kintsugi check

# Show system info
kintsugi info

# Install optional features
kintsugi install gpu
kintsugi install viz --conda  # Use conda where available

# Generate configuration template
kintsugi template -o config.json

# Run registration workflow
kintsugi register config.json --dry-run
kintsugi register config.json
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
```

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

### 3. Signal Isolation
Autofluorescence subtraction, filtering, and final processing to isolate signal.
- [notebooks/3_Signal_Isolation.ipynb](notebooks/3_Signal_Isolation.ipynb)
- **Requires:** Base only

### 4. Segmentation Analysis
InstanSeg segmentation, feature extraction, and spatial analysis.
- [notebooks/4_Segmentation_Analysis.ipynb](notebooks/4_Segmentation_Analysis.ipynb)
- **Requires:** `dl`, `viz`, `analysis`

### 5. Vessel Analysis
Specialized analysis for vessel structures.
- [notebooks/Vessel_Analysis.ipynb](notebooks/Vessel_Analysis.ipynb)
- **Requires:** `viz`, `analysis`

### Running Notebooks

Launch VS Code from the activated environment:
```bash
conda activate KINTSUGI
code .
```

**Important**: Always launch VS Code from the activated conda environment to ensure all packages are available.

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

## Development

### Setup Development Environment

```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
conda env create -f envs/env-linux.yml
conda activate KINTSUGI
pip install -e ".[dev]"
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
