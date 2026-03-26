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

# Activate and verify
conda activate KINTSUGI
kintsugi check
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

### GPU-Accelerated Image Processing

KINTSUGI uses CuPy for GPU-accelerated image processing operations:

| Operation | GPU Acceleration |
|-----------|------------------|
| Illumination Correction | CuPy FFT-based BaSiC algorithm |
| Stitching | CuPy phase correlation |
| Deconvolution | CuPy Lucy-Richardson with FFT |
| Extended Depth of Focus | CuPy variance projection |

```python
# Check GPU availability
from kintsugi.gpu import get_gpu_manager
gpu = get_gpu_manager()
print(gpu.summary())

# GPU-accelerated illumination correction
from kintsugi.kcorrect_gpu import KCorrectGPU
corrector = KCorrectGPU(device_id=0)
flatfield, darkfield = corrector.estimate(images)

# Multi-GPU stitching
from notebooks.Kstitch._translation_computation import get_multi_gpu_accelerator
accelerator = get_multi_gpu_accelerator(tile_shape)
```

### GPU-Accelerated Single-Cell Analysis

For GPU-accelerated single-cell analysis (clustering, UMAP, spatial analysis with RAPIDS),
see the dedicated repository: [rapids_singlecell](https://github.com/smith6jt-cop/rapids_singlecell)

## External Dependencies

| Dependency | Purpose | Installation |
|------------|---------|--------------|
| **libvips** | High-performance image I/O | `conda install libvips` (Linux/macOS) or Zenodo (Windows) |
| **VALIS** | Image registration | Vendored as `Kreg` (no separate install needed) |
| **CuPy** | GPU image processing | `conda install cupy` or `pip install cupy-cuda12x` |
| **PyTorch** | Deep learning models (optional) | `pip install torch torchvision` |

> **Note:** Java, Maven, and FIJI/CLIJ2 are no longer required. KINTSUGI now uses pure Python
> implementations (CuPy/NumPy) for all processing including Extended Depth of Focus (EDF).

## HPC Installation (SLURM Clusters)

For HPC environments like UF HiPerGator, KINTSUGI provides a dedicated environment file (`envs/env-hpc.yml`) that installs all dependencies upfront — including GPU/CUDA packages, spatial analysis tools, and workflow orchestration. This avoids the common "install base then install extras" pattern that causes silent failures on HPC.

### Quick Start (Recommended)

```bash
module load conda
git clone https://github.com/smith6jt-cop/KINTSUGI.git
cd KINTSUGI
./scripts/install.sh --hpc
```

This single command:
1. Creates the KINTSUGI conda env from `envs/env-hpc.yml`
2. Installs PyTorch + CUDA 12.4 via conda (not pip — avoids version conflicts)
3. Installs CUDA runtime libraries (libcufft, libcublas) needed by CuPy
4. Deploys LD_LIBRARY_PATH activation scripts (fixes libstdc++ mismatch on old OS kernels)
5. Copies CUDA headers for CuPy JIT compilation
6. Applies the SLURM TRES patch for Snakemake
7. Runs `kintsugi check` and fails if required dependencies are broken

### Manual HPC Installation

If you prefer step-by-step control:

```bash
module load conda

# Create environment with full GPU/CUDA/analysis stack
conda env create -f envs/env-hpc.yml
conda activate KINTSUGI

# Deploy LD_LIBRARY_PATH fix + CUDA headers
kintsugi fix-hpc

# CRITICAL: Re-activate to pick up the new LD_LIBRARY_PATH
conda deactivate && conda activate KINTSUGI

# Apply SLURM compatibility patch
kintsugi patch slurm

# Verify everything works
kintsugi check
```

### Fixing a Broken HPC Environment

If `kintsugi check` shows `[MISSING]` for packages that ARE installed (scikit-learn, matplotlib, pyvips, stackview), this is the **libstdc++ version mismatch** — the system's C++ library is too old and loading before conda's. Fix without recreating the env:

```bash
conda activate KINTSUGI
kintsugi fix-hpc
conda deactivate && conda activate KINTSUGI
kintsugi check
```

For comprehensive troubleshooting, see [HPC Troubleshooting](HPC_TROUBLESHOOTING.md).

### env-hpc.yml vs env-linux.yml

| Feature | `env-linux.yml` (Desktop) | `env-hpc.yml` (HPC) |
|---------|---------------------------|----------------------|
| PyTorch + CUDA | Not included (install via `kintsugi install gpu`) | Included via conda pytorch channel |
| CUDA runtime libs | Not included | `cuda-libraries`, `cuda-cudart-dev` |
| libstdcxx-ng | Not pinned | Pinned `>=12.0` (fixes GLIBCXX mismatch) |
| scanpy/analysis | Not included (install via `kintsugi install analysis`) | Included via conda-forge |
| Napari | Not included | Not included (no display on HPC) |
| Snakemake/SLURM | Included | Included |
| Target use | Workstations, laptops | SLURM clusters (HiPerGator, etc.) |

### HPC-Specific Notes

- **Login nodes vs compute nodes**: Install on login nodes (no GPU). The env includes torch and CuPy, which are importable without GPU hardware. `kintsugi check` reporting "CUDA not available" on login nodes is expected.
- **Never install packages inside SLURM jobs.** The conda env must be fully built before submitting.
- **Cache redirection**: SLURM jobs automatically redirect pip/torch/numba caches to `/blue/` storage via account-specific scripts.
- **SLURM TRES patch**: SLURM >= 24.11 requires `kintsugi patch slurm`. This is auto-applied by `install.sh --hpc` and `kintsugi install all`.
- **Verification**: Run `./scripts/verify_hpc_env.sh` for a comprehensive environment check beyond `kintsugi check`.

## Updating Existing Projects

When you update the KINTSUGI repository (`git pull`), project copies of notebooks, modules, and workflow scripts may become stale. See the **Updating Existing Projects** section in the [README](../README.md) for full details.

**Key commands:**

```bash
# Sync notebooks and Python modules to all project folders
python scripts/sync_to_projects.py          # Auto-discovers KINTSUGI_Projects/*/notebooks
python scripts/sync_to_projects.py --dry-run  # Preview changes

# Refresh workflow scripts for a single project
kintsugi workflow config /path/to/project --force

# Bulk-refresh workflow scripts across all projects
for d in /path/to/KINTSUGI_Projects/*/workflow/scripts; do
  cp workflow/scripts/*.py "$d/"
done
```

> **Note:** A git post-commit hook automatically runs `sync_to_projects.py` after each commit.
> Notebooks use `%autoreload 2` — no kernel restart needed after module updates.

## Verifying Installation

```bash
# Check all dependencies including GPU
kintsugi check

# Or via Python
python -c "
from kintsugi.gpu import get_gpu_manager
print(get_gpu_manager().summary())
"
```
