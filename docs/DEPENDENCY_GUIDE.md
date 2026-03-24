# Dependency Safety Guide

KINTSUGI uses a single conda environment with optional pip install groups. This guide explains how to avoid dependency conflicts and recover from broken environments.

## Safe Install Order

Always install groups in this order to prevent conflicts:

```bash
# 1. Create base environment (conda)
conda env create -f envs/env-linux.yml
conda activate KINTSUGI

# 2. GPU acceleration (installs CUDA torch + CuPy + CUDA runtime libs)
kintsugi install gpu

# 3. Workflow orchestration (HPC/SLURM only)
kintsugi install workflow
kintsugi patch slurm              # Required for SLURM >= 24.11

# 4. Additional groups as needed
kintsugi install dl               # Deep learning segmentation
kintsugi install analysis         # Spatial analysis (scanpy, scimap)
kintsugi install viz              # Napari visualization
kintsugi install bio              # Bio formats I/O

# 5. Verify
kintsugi check --strict
```

## Critical Rules

### Never install torch directly with pip
```bash
# WRONG - installs CPU-only torch
pip install torch

# RIGHT - uses CUDA index for GPU-enabled torch
kintsugi install gpu
```

All `kintsugi install` commands for torch-using groups (`gpu`, `dl`, `denoise`, `kronos`) automatically use the correct CUDA index URL.

### The numpy < 2.0 constraint
KINTSUGI requires `numpy >= 1.24, < 2.0`. Many packages in the `analysis` and `bio` groups may pull numpy 2.x as a transitive dependency. The `constraints.txt` file at the repo root prevents this:

```bash
# All kintsugi install commands automatically use constraints.txt
kintsugi install analysis  # numpy stays < 2.0

# Manual pip install should also use constraints
pip install -c constraints.txt scanpy scimap
```

### CuPy needs CUDA runtime libraries
`cupy-cuda12x` from pip only provides Python bindings. The actual CUDA runtime (`libcufft`, `libcublas`) must come from conda:

```bash
# kintsugi install gpu handles this automatically
kintsugi install gpu

# If you installed cupy manually and get libcufft errors:
conda install cuda-libraries cuda-cudart-dev -c nvidia -y
cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/
```

### SLURM TRES patch (HPC only)
SLURM >= 24.11 sets `SLURM_TRES_PER_TASK` in GPU job environments, which conflicts with the Snakemake jobstep plugin. This patch must be applied after every pip upgrade of the plugin:

```bash
kintsugi patch slurm
kintsugi check --strict   # Verifies patch is applied
```

## Validation

```bash
kintsugi check                     # Full dependency report
kintsugi check --strict            # Exit code 1 if any errors
kintsugi check --for pipeline      # Validate GPU processing deps
kintsugi check --for analysis      # Validate analysis deps
kintsugi check --for segmentation  # Validate segmentation deps
```

## Troubleshooting

### "numpy >= 2.0 detected"
```bash
pip install 'numpy>=1.24,<2.0'
```

### "CPU-only PyTorch detected"
```bash
pip uninstall torch torchvision
kintsugi install gpu
```

### "CuPy installed but CUDA libraries missing"
```bash
conda install cuda-libraries cuda-cudart-dev -c nvidia -y
cp -r $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/
```

### "SLURM_TRES_PER_TASK patch NOT applied"
```bash
kintsugi patch slurm
```

### Full environment rebuild
If the environment is too broken to repair:
```bash
conda deactivate
conda env remove -n KINTSUGI
conda env create -f envs/env-linux.yml
conda activate KINTSUGI
kintsugi install gpu
kintsugi install workflow && kintsugi patch slurm  # HPC only
kintsugi check --strict
```

## Feature Requirements by Notebook

| Notebook | Required Groups | Install Command |
|----------|----------------|-----------------|
| 1_Single_Channel_Eval | `gpu` | `kintsugi install gpu` |
| 2_Cycle_Processing | `gpu` | `kintsugi install gpu` |
| 3_Signal_Isolation_QC | (none) | Base install only |
| 4_Segmentation_Analysis | `dl`, `viz`, `analysis` | `kintsugi install dl viz analysis` |
| 5_Cluster_Analysis | `analysis` | `kintsugi install analysis` |
| 2.5_Vessel_3D_Segmentation | `gpu`, `analysis` | `kintsugi install gpu analysis` |
