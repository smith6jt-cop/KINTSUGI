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
  - `maven-3.9.9`
  - `java-jdk21`
  - `PyVips-dev-8.16`
  - `FIJI` with Clij2 plugin

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

## macOS Installation

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
| **Java 11+** | BioFormats support | `conda install openjdk=11` or Zenodo |
| **Maven** | Java dependency management | `conda install maven` or Zenodo |
| **VALIS** | Image registration | `pip install valis-wsi` (included) |
| **FIJI + Clij2** | ImageJ integration | Download from Zenodo |
| **CuPy** | GPU acceleration for deconvolution (optional) | `pip install cupy-cuda11x` |
