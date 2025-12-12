# Installation

## Quick Start (All Platforms)

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

## Windows Installation

### Option A: Using Installation Script (Recommended)

```powershell
# Open PowerShell as Administrator
cd C:\Users\[your username]
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

# Run installation script
.\scripts\install.ps1
```

### Option B: Manual Installation

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

### Download Windows Dependencies

Windows requires additional binary dependencies from Zenodo:
- Download from: [https://zenodo.org/records/14969214](https://zenodo.org/records/14969214)
- Extract to the KINTSUGI folder:
  - `PyVips-dev-8.16` (required for image I/O)

> **Note:** Java, Maven, and FIJI are no longer required. KINTSUGI now uses
> pure Python implementations (CuPy/NumPy) for all processing including EDF.

## Linux Installation

### Option A: Using Installation Script (Recommended)

```bash
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
chmod +x scripts/install.sh
./scripts/install.sh
```

### Option B: Manual Installation

```bash
# Install system dependencies (Java/Maven no longer required)
sudo apt-get update
sudo apt-get install -y libvips-dev

# Clone and setup
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

conda env create -f envs/env-linux.yml
conda activate KINTSUGI
pip install -e .
```

## macOS Installation

```bash
# Install system dependencies (Java/Maven no longer required)
brew install vips

# Clone and setup
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI

conda env create -f envs/env-macos.yml
conda activate KINTSUGI
pip install -e .
```

## GPU Acceleration (Optional)

KINTSUGI supports multi-GPU acceleration for significantly faster processing.

### Basic GPU Support (CuPy + PyTorch)

For GPU-accelerated deconvolution, stitching, and denoising:

```bash
# Install via kintsugi CLI (recommended)
kintsugi install gpu

# Or manually for CUDA 12.x
pip install cupy-cuda12x torch torchvision

# Or for CUDA 11.x
pip install cupy-cuda11x torch torchvision
```

### Multi-GPU Support

KINTSUGI automatically detects and uses all available non-integrated NVIDIA GPUs:

```python
from kintsugi.gpu import get_gpu_manager

# Check available GPUs
gpu = get_gpu_manager()
print(gpu.summary())
# Output: Found 2 GPU(s):
#   [0] NVIDIA A100 (40.0 GB, CC 8.0)
#   [1] NVIDIA A100 (40.0 GB, CC 8.0)

# PyTorch models automatically use DataParallel
from kintsugi.denoise import CAREDenoiser
denoiser = CAREDenoiser()  # Uses all GPUs for inference

# Multi-GPU stitching
from notebooks.Kstitch._translation_computation import get_multi_gpu_accelerator
accelerator = get_multi_gpu_accelerator((1024, 1024))
```

### RAPIDS GPU-Accelerated Data Science (Optional)

For GPU-accelerated pandas and scikit-learn operations (10-100x speedup on large datasets):

```bash
# Option 1: Via conda (recommended)
# Use cuda-version matching your driver (12.4+ for newer GPUs like B200/H100)
conda install -c rapidsai -c conda-forge -c nvidia \
    cudf=24.10 cuml=24.10 cugraph=24.10 \
    cuda-version=12.6

# Option 2: Via pip
pip install cudf-cu12 cuml-cu12 cugraph-cu12 \
    --extra-index-url=https://pypi.nvidia.com

# Option 3: Via kintsugi optional dependencies
pip install -e ".[rapids]" --extra-index-url=https://pypi.nvidia.com
```

> **Note:** Check your CUDA version with `nvidia-smi`. Newer GPUs (B200, H100, etc.)
> require CUDA 12.4+. Match the `cuda-version` to your driver capabilities.

RAPIDS provides:
- **cuDF**: GPU DataFrame (drop-in pandas replacement)
- **cuML**: GPU machine learning (drop-in scikit-learn replacement)
- **cuDF-pandas**: Transparent pandas acceleration

```python
from kintsugi.rapids import get_rapids_manager, enable_cudf_pandas

# Check RAPIDS status
rapids = get_rapids_manager()
print(rapids.summary())

# GPU-accelerated operations
labels, centers = rapids.kmeans(data, n_clusters=10)
embedding = rapids.umap(data, n_components=2)

# Enable transparent pandas acceleration (call before importing pandas)
enable_cudf_pandas()
import pandas as pd  # Now GPU-accelerated!
```

## External Dependencies

| Dependency | Purpose | Installation |
|------------|---------|--------------|
| **libvips** | High-performance image I/O | `conda install libvips` (Linux/macOS) or Zenodo (Windows) |
| **VALIS** | Image registration | `pip install valis-wsi` (included) |
| **CuPy** | GPU acceleration (optional) | `conda install cupy` or `pip install cupy-cuda12x` |
| **PyTorch** | Deep learning models (optional) | `pip install torch torchvision` |
| **RAPIDS** | GPU data science (optional) | See RAPIDS section above |

> **Note:** Java, Maven, and FIJI/CLIJ2 are no longer required. KINTSUGI now uses pure Python
> implementations (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF).

## Verifying GPU Setup

```bash
# Check all dependencies including GPU
kintsugi check

# Or via Python
python -c "
from kintsugi.gpu import get_gpu_manager
from kintsugi.rapids import get_rapids_manager

print('=== GPU Status ===')
print(get_gpu_manager().summary())

print('\\n=== RAPIDS Status ===')
print(get_rapids_manager().summary())
"
```
