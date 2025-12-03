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
  - [Quick Start (All Platforms)](#quick-start-all-platforms)
  - [Windows Installation](#windows-installation)
  - [Linux Installation](#linux-installation)
  - [macOS Installation](#macos-installation)
  - [Verify Installation](#verify-installation)
- [External Dependencies](#external-dependencies)
- [Usage](#usage)
  - [Command Line Interface](#command-line-interface)
  - [Python API](#python-api)
- [Notebooks](#notebooks)
- [Troubleshooting](#troubleshooting)
- [Development](#development)

## Installation

### Quick Start (All Platforms)

1. **Install Miniconda** (if not already installed):
   - Download from [Anaconda](https://www.anaconda.com/download/success#miniconda)
   - Follow platform-specific installation instructions

2. **Clone and Install KINTSUGI**:
```bash
# Clone the repository
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# Create conda environment (choose your platform)
# Linux:
conda env create -f envs/env-linux.yml

# Windows:
conda env create -f envs/env-windows.yml

# macOS:
conda env create -f envs/env-macos.yml

# Activate and install package
conda activate KINTSUGI
pip install -e .
```

3. **Verify Installation**:
```bash
kintsugi check
```

### Windows Installation

#### Option A: Using Installation Script (Recommended)
```powershell
# Open PowerShell as Administrator
cd C:\Users\[your username]
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# Run installation script
.\scripts\install.ps1
```

#### Option B: Manual Installation
```powershell
# Open Anaconda Prompt
conda update -n base conda
conda install -n base conda-libmamba-solver
conda config --set solver libmamba

# Clone repository
cd C:\Users\[your username]
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# Create environment
conda env create -f envs/env-windows.yml
conda activate KINTSUGI

# Install KINTSUGI package
pip install -e .
```

#### Download Windows Dependencies
Windows requires additional binary dependencies from Zenodo:
- Download from: [https://zenodo.org/records/14969214](https://zenodo.org/records/14969214)
- Extract to the KINTSUGI folder:
  - `maven-3.9.9`
  - `java-jdk21`
  - `PyVips-dev-8.16`
  - `FIJI` with Clij2 plugin

#### GPU Acceleration (Optional)
For GPU-accelerated deconvolution, install CuPy for your CUDA version:
```bash
# For CUDA 11.x
pip install cupy-cuda11x

# For CUDA 12.x
pip install cupy-cuda12x
```

### Linux Installation

#### Option A: Using Installation Script (Recommended)
```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
chmod +x scripts/install.sh
./scripts/install.sh
```

#### Option B: Manual Installation
```bash
# Install system dependencies
sudo apt-get update
sudo apt-get install -y libvips-dev openjdk-11-jdk maven

# Clone and setup
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

conda env create -f envs/env-linux.yml
conda activate KINTSUGI
pip install -e .
```

### macOS Installation

```bash
# Install system dependencies
brew install vips openjdk@11 maven

# Clone and setup
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

conda env create -f envs/env-macos.yml
conda activate KINTSUGI
pip install -e .
```

### Verify Installation

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

KINTSUGI relies on several external dependencies:

| Dependency | Purpose | Installation |
|------------|---------|--------------|
| **libvips** | High-performance image I/O | `conda install libvips` (Linux/macOS) or Zenodo (Windows) |
| **Java 11+** | BioFormats support | `conda install openjdk=11` or Zenodo |
| **Maven** | Java dependency management | `conda install maven` or Zenodo |
| **VALIS** | Image registration | `pip install valis-wsi` (included) |
| **FIJI + Clij2** | ImageJ integration | Download from Zenodo |
| **CuPy** | GPU acceleration for deconvolution (optional) | `pip install cupy-cuda11x` |

## Usage

### Command Line Interface

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

## Notebooks

The following Jupyter notebooks provide step-by-step workflows:

### 1. Parameter Tuning and Testing
Test illumination correction, stitching, deconvolution, and EDoF.
- [notebooks/1_Single_Channel_Eval.ipynb](notebooks/1_Single_Channel_Eval.ipynb)

### 2. Batch Processing
Batch processing for illumination correction, stitching, deconvolution, EDoF, and registration.
- [notebooks/2_Cycle_Processing.ipynb](notebooks/2_Cycle_Processing.ipynb)

### 3. Signal Isolation
Autofluorescence subtraction, filtering, and final processing to isolate signal.
- [notebooks/3_Signal_Isolation.ipynb](notebooks/3_Signal_Isolation.ipynb)

### 4. Segmentation Analysis
InstanSeg segmentation, feature extraction, and spatial analysis.
- [notebooks/4_Segmentation_Analysis.ipynb](notebooks/4_Segmentation_Analysis.ipynb)

### 5. Vessel Analysis
Specialized analysis for vessel structures.
- [notebooks/Vessel_Analysis.ipynb](notebooks/Vessel_Analysis.ipynb)

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

**Java not found**
```bash
# Check Java installation
java -version

# If missing, install via conda
conda install openjdk=11
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
