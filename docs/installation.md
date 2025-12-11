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

For GPU-accelerated deconvolution, install CuPy for your CUDA version:

```bash
# For CUDA 11.x
pip install cupy-cuda11x

# For CUDA 12.x
pip install cupy-cuda12x
```

## External Dependencies

| Dependency | Purpose | Installation |
|------------|---------|--------------|
| **libvips** | High-performance image I/O | `conda install libvips` (Linux/macOS) or Zenodo (Windows) |
| **VALIS** | Image registration | `pip install valis-wsi` (included) |
| **CuPy** | GPU acceleration (optional) | `conda install cupy` or `pip install cupy-cuda12x` |

> **Note:** Java, Maven, and FIJI/CLIJ2 are no longer required. KINTSUGI now uses pure Python
> implementations (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF).
