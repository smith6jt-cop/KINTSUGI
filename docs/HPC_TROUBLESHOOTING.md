# HPC Troubleshooting Guide

This guide covers the most common environment and installation issues on HPC clusters (HiPerGator). For general installation instructions, see [installation.md](installation.md).

## Issue 1: Packages Show [MISSING] Despite Being Installed

**Symptoms:**
```
[MISSING]  scikit-learn
           /lib64/libstdc++.so.6: version `CXXABI_1.3.15' not found
[MISSING]  matplotlib
           /lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found
[ERROR]    pyvips
           cannot load library: /lib64/libstdc++.so.6: version `GLIBCXX_3.4.30' not found
[MISSING]  stackview
```

**Root cause:** The system's `/lib64/libstdc++.so.6` is too old and gets loaded before conda's newer version. The packages are installed but fail to import because their compiled extensions need newer C++ ABI symbols.

**Fix:**
```bash
conda activate KINTSUGI
kintsugi fix-hpc
conda deactivate && conda activate KINTSUGI
kintsugi check
```

`kintsugi fix-hpc` deploys conda activation scripts that prepend `$CONDA_PREFIX/lib` to `LD_LIBRARY_PATH`, ensuring conda's `libstdc++.so.6` loads first.

**Prevention:** Use `./scripts/install.sh --hpc` or `envs/env-hpc.yml` which pins `libstdcxx-ng>=12.0` and deploys the activation scripts automatically.

## Issue 2: `kintsugi install all` Completes But Packages Are Missing

**Symptoms:** `kintsugi install all` prints "Installation complete!" but `kintsugi check` shows `[SKIP]` for scanpy, phenograph, instanseg, etc.

**Root cause:** `install all` in pip mode does a single `pip install -e ".[all_extras]"` which can silently drop packages that fail to resolve. Additionally, the pip-only install path skips conda commands needed for CUDA runtime libraries.

**Fix (for future installs):** On HPC, use `install all --conda` or (better) use `envs/env-hpc.yml`:
```bash
# Option A: Conda mode (runs conda commands for GPU groups)
kintsugi install all --conda

# Option B: Start fresh with HPC env file (recommended)
conda env remove -n KINTSUGI -y
conda env create -f envs/env-hpc.yml
conda activate KINTSUGI
kintsugi fix-hpc
conda deactivate && conda activate KINTSUGI
```

**Fix (for missing packages in existing env):** Install individually to see specific errors:
```bash
pip install scanpy
pip install phenograph
pip install instanseg instanseg-torch
```

**Prevention:** The `install all` command now auto-detects HPC environments and switches to conda mode, which installs conda groups first and then pip-only packages individually (avoiding the mega-extras resolution that silently drops packages).

## Issue 3: CuPy Imports But CUDA Operations Fail

**Symptoms:** `import cupy` succeeds but `cupy.array([1.0])` fails with:
```
ImportError: libcufft.so.10: cannot open shared object file
```

**Root cause:** `cupy-cuda12x` (pip) only provides Python bindings. The actual CUDA runtime libraries (libcufft, libcublas, libcusolver) must be installed separately via conda.

**Fix:**
```bash
conda install cuda-libraries cuda-cudart-dev -c nvidia -y
# Copy headers for CuPy JIT
cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/ 2>/dev/null
```

**Prevention:** Use `envs/env-hpc.yml` which includes `cuda-libraries` and `cuda-cudart-dev` in the conda section.

## Issue 4: PyTorch is CPU-Only

**Symptoms:**
```
[ERROR]    torch-build
           CPU-only PyTorch detected! GPU processing will fail.
```

**Root cause:** PyTorch was installed from pip without the CUDA index URL, or a later pip operation overwrite the conda-installed GPU version.

**Fix:**
```bash
# Remove existing torch
pip uninstall torch torchvision -y

# Re-install from PyTorch CUDA channel
pip install torch torchvision --index-url https://download.pytorch.org/whl/cu124

# Or via conda (preferred on HPC):
conda install pytorch torchvision pytorch-cuda=12.4 -c pytorch -c nvidia -c conda-forge
```

**Prevention:** Use `envs/env-hpc.yml` which installs PyTorch from the `pytorch` conda channel with `pytorch-cuda=12.4`. Never run bare `pip install torch` on HPC — always use the CUDA index URL.

## Issue 5: SLURM Jobs Fail With TRES Error

**Symptoms:** Snakemake SLURM jobs fail immediately with:
```
error: Invalid --tres-per-task specification
```

**Root cause:** SLURM >= 24.11 changed the `SLURM_TRES_PER_TASK` environment variable format, which conflicts with the Snakemake jobstep executor plugin.

**Fix:**
```bash
kintsugi patch slurm
```

This patches the `snakemake-executor-plugin-slurm-jobstep` plugin to unset `SLURM_TRES_PER_TASK`. Must be re-applied after upgrading the plugin.

**Prevention:** `./scripts/install.sh --hpc` and `kintsugi install all` auto-apply this patch.

## Issue 6: torch-cuda Shows [SKIP] on Login Node

**Symptoms:**
```
[SKIP]     torch-cuda (optional)
```

**This is expected.** Login nodes have no GPU hardware, so `torch.cuda.is_available()` returns False. PyTorch and CuPy will work correctly on compute nodes.

To verify the torch build is CUDA-enabled (not CPU-only), check:
```bash
python -c "import torch; print(f'CUDA build: {torch.version.cuda}')"
```

If this prints `CUDA build: 12.4` (or similar), you're fine.

## Verification

Run the comprehensive verification script:
```bash
./scripts/verify_hpc_env.sh
```

Or the quick check:
```bash
kintsugi check
```

## Getting Help

If none of the above resolves your issue:
1. Run `kintsugi check` and save the full output
2. Run `conda list > env_packages.txt` to capture installed packages
3. Report the issue at https://github.com/smith6jt-cop/KINTSUGI/issues
