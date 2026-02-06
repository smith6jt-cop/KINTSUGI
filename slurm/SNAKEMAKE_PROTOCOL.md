# KINTSUGI Snakemake Integration Protocol

This document outlines a protocol for wrapping the existing KINTSUGI SLURM job
scripts with Snakemake to gain declarative dependency management, automatic
re-execution on failure, and fine-grained per-file tracking.

---

## 1. Current SLURM Architecture Summary

The existing pipeline has four processing stages executed as SLURM array jobs:

| Step | Script | Input | Output | Array Index |
|------|--------|-------|--------|-------------|
| 1. Correction | `01_correction.sh` | (passthrough) | — | cycle |
| 2. Stitching | `02_stitching.sh` | `data/raw/cyc{NN}/` | `data/processed/stitched/cyc{NN}/CH{M}/z{PP}.tif` | cycle |
| 3. Deconvolution | `03_deconvolution.sh` | `data/processed/stitched/cyc{NN}/` | `data/processed/deconvolved/cyc{NN}/CH{M}/z{PP}.tif` | cycle |
| 4. EDF | `04_edf.sh` | `data/processed/deconvolved/cyc{NN}/` | `data/processed/edf/cyc{NN}/CH{M}_edf.tif` | cycle |

**Key characteristics of the current system:**
- `submit.sh` calculates dual-pool resource allocation (GPU + CPU slots)
- Jobs use `--dependency=aftercorr` for per-array-task dependency chains
- `.complete` marker files provide inter-job synchronization across NFS
- `KINTSUGI_DEVICE_MODE` env var selects GPU (CuPy) or CPU (NumPy) backend
- Skip-existing logic inside each script checks for completed output
- Each array task processes **one full cycle** (all channels, all z-planes)

**What Snakemake replaces:**
- `submit.sh` orchestration and dependency wiring
- `.complete` marker file polling (Snakemake tracks outputs natively)
- Skip-existing logic (Snakemake's file-based DAG handles this)
- Manual `--array` limit calculation (Snakemake's `--jobs` flag)

**What Snakemake does NOT replace:**
- The Python processing code inside each job script
- `config.sh` resource definitions
- `utils.sh` logging and QC utilities
- The burst monitor (handled differently; see Section 8)

---

## 2. Project Layout

```
my_project/
├── data/
│   ├── raw/
│   │   ├── cyc01/
│   │   ├── cyc02/
│   │   └── ...
│   └── processed/
│       ├── stitched/
│       │   └── cyc{NN}/CH{M}/z{PP}.tif
│       ├── deconvolved/
│       │   └── cyc{NN}/CH{M}/z{PP}.tif
│       └── edf/
│           └── cyc{NN}/CH{M}_edf.tif
├── meta/
│   ├── experiment.json
│   └── CHANNELNAMES.txt
├── slurm/
│   ├── config.sh
│   ├── jobs/            ← Existing job scripts (reused)
│   └── utils.sh
├── workflow/             ← NEW: Snakemake files
│   ├── Snakefile
│   ├── config.yaml
│   ├── profiles/
│   │   └── slurm/
│   │       └── config.yaml
│   ├── envs/
│   │   └── kintsugi.yaml
│   └── scripts/
│       ├── stitch.py
│       ├── deconvolve.py
│       └── edf.py
└── kintsugi_project.json
```

---

## 3. Configuration Translation

Snakemake uses a YAML config file. Translate `config.sh` into `workflow/config.yaml`:

```yaml
# workflow/config.yaml
# Auto-generated from slurm/config.sh — or maintained independently

project_dir: /path/to/my_project
kintsugi_dir: /path/to/KINTSUGI

# Processing scope
cycles: [1, 2, 3, 4, 5, 6, 7]
channels: [1, 2, 3, 4]
# z-planes are auto-detected from raw data

# Output format
output_format: tiff   # or zarr

# Microscope parameters (override experiment.json if needed)
# If omitted, scripts read from meta/experiment.json automatically
microscope:
  xy_pixel_size: 377       # nm
  z_step_size: 1500        # nm
  numerical_aperture: 0.75
  tissue_refractive_index: 1.44
  wavelengths:
    1: [358.0, 461.0]
    2: [753.0, 775.0]
    3: [560.0, 575.0]
    4: [648.0, 668.0]

# Tile grid
tile_rows: 5
tile_cols: 5
tile_overlap: 0.1

# SLURM resource defaults (used in profile)
resources:
  conda_env: kintsugi
  partition_gpu: hpg-b200
  partition_cpu: hpg-default
  qos: maigan-b
  account: maigan
  gpu_type: b200
```

---

## 4. Snakefile Design

The core principle: **each rule wraps one existing job script's Python logic,
operating on a single cycle**. Snakemake's DAG replaces `submit.sh`'s
`aftercorr` dependency chain.

```python
# workflow/Snakefile

import json
from pathlib import Path

configfile: "config.yaml"

PROJECT = config["project_dir"]
KINTSUGI = config["kintsugi_dir"]
CYCLES = config["cycles"]
CHANNELS = config["channels"]

# Load experiment metadata once
EXPERIMENT_JSON = Path(PROJECT) / "meta" / "experiment.json"
if EXPERIMENT_JSON.exists():
    with open(EXPERIMENT_JSON) as f:
        EXPERIMENT = json.load(f)
    N_ZPLANES = EXPERIMENT.get("n_zplanes", 15)
else:
    EXPERIMENT = {}
    N_ZPLANES = 15

ZPLANES = list(range(1, N_ZPLANES + 1))


# =========================================================================
# Target rule: request all final EDF outputs
# =========================================================================
rule all:
    input:
        expand(
            f"{PROJECT}/data/processed/edf/cyc{{cycle:02d}}/CH{{ch}}_edf.tif",
            cycle=CYCLES,
            ch=CHANNELS,
        )


# =========================================================================
# Rule 1: Stitching (includes illumination correction)
# One job per cycle — produces all channels x all z-planes
# =========================================================================
rule stitch:
    input:
        raw_dir=f"{PROJECT}/data/raw/cyc{{cycle}}"
    output:
        # Declare one sentinel output per cycle to represent completion.
        # The script produces all CH*/z*.tif files internally.
        sentinel=f"{PROJECT}/data/processed/stitched/cyc{{cycle}}/.snakemake_complete"
    params:
        project_dir=PROJECT,
        kintsugi_dir=KINTSUGI,
        channels=CHANNELS,
    resources:
        mem_mb=48000,
        runtime=240,       # minutes
        slurm_partition=config["resources"]["partition_gpu"],
        gpus=1,
        cpus_per_task=4,
    conda:
        "envs/kintsugi.yaml"
    log:
        f"{PROJECT}/slurm/logs/snakemake/stitch_cyc{{cycle}}.log"
    script:
        "scripts/stitch.py"


# =========================================================================
# Rule 2: Deconvolution
# One job per cycle — depends on stitching sentinel
# =========================================================================
rule deconvolve:
    input:
        stitch_done=f"{PROJECT}/data/processed/stitched/cyc{{cycle}}/.snakemake_complete"
    output:
        sentinel=f"{PROJECT}/data/processed/deconvolved/cyc{{cycle}}/.snakemake_complete"
    params:
        project_dir=PROJECT,
        kintsugi_dir=KINTSUGI,
        channels=CHANNELS,
    resources:
        mem_mb=48000,
        runtime=240,
        slurm_partition=config["resources"]["partition_gpu"],
        gpus=1,
        cpus_per_task=4,
    conda:
        "envs/kintsugi.yaml"
    log:
        f"{PROJECT}/slurm/logs/snakemake/decon_cyc{{cycle}}.log"
    script:
        "scripts/deconvolve.py"


# =========================================================================
# Rule 3: Extended Depth of Focus
# One job per cycle — depends on deconvolution sentinel
# =========================================================================
rule edf:
    input:
        decon_done=f"{PROJECT}/data/processed/deconvolved/cyc{{cycle}}/.snakemake_complete"
    output:
        # EDF produces one file per channel (2D projection)
        edf_files=expand(
            f"{PROJECT}/data/processed/edf/cyc{{{{cycle}}}}/CH{{ch}}_edf.tif",
            ch=CHANNELS,
        ),
        sentinel=f"{PROJECT}/data/processed/edf/cyc{{cycle}}/.snakemake_complete"
    params:
        project_dir=PROJECT,
        kintsugi_dir=KINTSUGI,
        channels=CHANNELS,
    resources:
        mem_mb=16000,
        runtime=60,
        slurm_partition=config["resources"]["partition_gpu"],
        gpus=1,
        cpus_per_task=4,
    conda:
        "envs/kintsugi.yaml"
    log:
        f"{PROJECT}/slurm/logs/snakemake/edf_cyc{{cycle}}.log"
    script:
        "scripts/edf.py"
```

### Design decisions

**Sentinel files (`.snakemake_complete`)**: Stitching and deconvolution each
produce dozens of output files (channels x z-planes). Declaring every file as a
Snakemake output is possible but creates a very large DAG. Instead, each rule
declares a single sentinel file that is touched after the script verifies all
expected outputs exist. This mirrors the existing `.complete` marker pattern
but is managed by Snakemake.

**One rule per cycle**: The existing scripts process an entire cycle per
invocation (all channels, all z-planes). This matches the SLURM array-task
granularity and keeps the wrapper scripts thin.

**Why not one rule per channel?** The stitching step computes a stitching model
from CH1 and applies it to all channels — the channels are not independent. The
deconvolution step similarly processes channels sequentially on one GPU.
Splitting by channel would require refactoring the core processing code.

---

## 5. Wrapper Scripts

Each `workflow/scripts/*.py` file is a thin adapter between Snakemake and the
existing KINTSUGI processing logic. These extract the same Python code that
lives inside the current bash heredoc scripts.

### `workflow/scripts/stitch.py`

```python
"""Snakemake wrapper for KINTSUGI stitching (replaces 02_stitching.sh)."""
import sys
import os
import json
from pathlib import Path

# Snakemake injects: snakemake.input, snakemake.output, snakemake.params, etc.
PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = snakemake.params.channels

# Setup Python path (same as job scripts)
sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

# Load experiment metadata
experiment_config = {}
config_path = PROJECT_DIR / "meta" / "experiment.json"
if config_path.exists():
    with open(config_path) as f:
        experiment_config = json.load(f)

# Import processing functions
from Kstitch.stitching import stitch_images
from kintsugi.kcorrect_gpu import KCorrectGPU
from kintsugi.stitch_blend import stitch_with_blending

# --- Core processing logic extracted from 02_stitching.sh ---
# (tile loading, BaSiC correction, stitching model, stitch+blend)
# ... (reuse the Python code from the existing heredoc) ...

# After processing, write sentinel
sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(f"cycle={CYCLE}\ncompleted={__import__('datetime').datetime.now().isoformat()}\n")
```

### `workflow/scripts/deconvolve.py`

```python
"""Snakemake wrapper for KINTSUGI deconvolution (replaces 03_deconvolution.sh)."""
import sys
import os
import json
from pathlib import Path

PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = snakemake.params.channels

sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

experiment_config = {}
config_path = PROJECT_DIR / "meta" / "experiment.json"
if config_path.exists():
    with open(config_path) as f:
        experiment_config = json.load(f)

from Kdecon import decon

# --- Core deconvolution logic extracted from 03_deconvolution.sh ---
# (wavelength lookup, per-channel decon loop, QC generation)
# No need for wait_for_input() — Snakemake guarantees input exists.
# ... (reuse the Python code from the existing heredoc) ...

sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(f"cycle={CYCLE}\ncompleted={__import__('datetime').datetime.now().isoformat()}\n")
```

### `workflow/scripts/edf.py`

```python
"""Snakemake wrapper for KINTSUGI EDF (replaces 04_edf.sh)."""
import sys
import os
import json
from pathlib import Path

PROJECT_DIR = Path(snakemake.params.project_dir)
KINTSUGI_DIR = Path(snakemake.params.kintsugi_dir)
CYCLE = int(snakemake.wildcards.cycle)
CHANNELS = snakemake.params.channels

sys.path.insert(0, str(PROJECT_DIR / "notebooks"))
sys.path.insert(0, str(KINTSUGI_DIR))

from kintsugi.edf import EDFProcessor

# --- Core EDF logic extracted from 04_edf.sh ---
# (quality gate, per-channel EDF, QC generation)
# ... (reuse the Python code from the existing heredoc) ...

sentinel = Path(snakemake.output.sentinel)
sentinel.parent.mkdir(parents=True, exist_ok=True)
sentinel.write_text(f"cycle={CYCLE}\ncompleted={__import__('datetime').datetime.now().isoformat()}\n")
```

---

## 6. SLURM Profile

Use Snakemake's SLURM executor plugin to submit each rule invocation as a
SLURM job. Create a profile at `workflow/profiles/slurm/config.yaml`:

```yaml
# workflow/profiles/slurm/config.yaml
# Snakemake >= 8.x with snakemake-executor-plugin-slurm

executor: slurm

# Default SLURM settings (overridden per-rule via resources:)
default-resources:
  slurm_account: maigan
  slurm_partition: hpg-b200
  mem_mb: 48000
  runtime: 240           # minutes
  cpus_per_task: 4
  gpus: 1

# Maximum concurrent jobs (replaces submit.sh's dual-pool calculation)
# Set this to your GPU_SLOTS + CPU_SLOTS from config.sh
jobs: 8

# Latency for NFS propagation (seconds)
latency-wait: 120

# Retry failed jobs (replaces manual resubmission)
retries: 2

# Keep going if independent jobs can still run
keep-going: true

# Use conda for environment management
use-conda: true
```

### Mapping resource directives

The `resources:` block in each rule maps to SLURM sbatch flags:

| Snakemake resource | SLURM flag | Source in config.sh |
|---|---|---|
| `mem_mb` | `--mem` | `MEM_STITCH`, `MEM_DECON`, `MEM_EDF` |
| `runtime` | `--time` | `TIME_STITCH`, `TIME_DECON`, `TIME_EDF` |
| `cpus_per_task` | `--cpus-per-task` | `CPUS_PER_TASK` |
| `gpus` | `--gpus-per-node` | `GPUS_PER_NODE` |
| `slurm_partition` | `--partition` | `PARTITION` |

---

## 7. Execution Commands

### Full pipeline
```bash
cd my_project/workflow
snakemake --profile profiles/slurm
```

### Specific cycles
```bash
snakemake --profile profiles/slurm \
    data/processed/edf/cyc01/CH1_edf.tif \
    data/processed/edf/cyc02/CH1_edf.tif
```

### Dry run (replaces `kintsugi slurm submit --dry-run`)
```bash
snakemake --profile profiles/slurm -n
```

### Just deconvolution + EDF (skip stitching)
```bash
snakemake --profile profiles/slurm \
    --allowed-rules deconvolve edf
```

### Force re-run of a specific cycle
```bash
snakemake --profile profiles/slurm \
    --forcerun data/processed/deconvolved/cyc03/.snakemake_complete
```

### Visualize the DAG
```bash
snakemake --dag | dot -Tpng > dag.png
```

---

## 8. GPU/CPU Dual-Pool Strategy

The existing `submit.sh` submits GPU and CPU jobs concurrently using
`KINTSUGI_DEVICE_MODE`. In Snakemake, there are two approaches:

### Option A: GPU-only (simpler)

Set `gpus: 1` on all rules. Snakemake's `--jobs` flag limits concurrency.
CPU-fallback happens automatically inside the scripts when CUDA init fails.

### Option B: Explicit dual-pool with rule variants (advanced)

Create GPU and CPU variants of each rule using Snakemake's `ruleorder` and
resource constraints:

```python
# GPU variant (preferred)
rule deconvolve_gpu:
    ...
    resources:
        gpus=1,
        mem_mb=48000,
        slurm_partition=config["resources"]["partition_gpu"],
    params:
        device_mode="gpu",

# CPU variant (fallback when GPUs saturated)
rule deconvolve_cpu:
    ...
    resources:
        gpus=0,
        mem_mb=48000,
        runtime=1200,  # 5x GPU time
        cpus_per_task=8,
        slurm_partition=config["resources"]["partition_cpu"],
    params:
        device_mode="cpu",

ruleorder: deconvolve_gpu > deconvolve_cpu
```

Then use Snakemake's `--resources gpus=3` flag to cap GPU usage:

```bash
snakemake --profile profiles/slurm -j 13 --resources gpus=3
```

This allows Snakemake to schedule 3 GPU jobs + 10 CPU jobs concurrently,
matching the existing dual-pool architecture. The `--resources gpus=3` flag
acts as a global GPU budget; once 3 GPU jobs are running, Snakemake
automatically falls back to the CPU variants for remaining cycles.

### Option C: Snakemake group jobs (batch multiple cycles)

Group multiple cycles into a single SLURM job to reduce scheduler overhead:

```python
rule deconvolve:
    ...
    group: "deconvolve"
```

Note: Use a fixed group name (not `"cycle_{cycle}"`) so Snakemake can batch
multiple cycles into one SLURM submission. A per-cycle group name would
create separate groups per cycle, defeating the batching intent.

---

## 9. Handling Burst Resources

The burst monitor (`burst_monitor.sh`) promotes preemptible jobs to guaranteed
QOS. With Snakemake, this is handled differently:

1. **Retries replace burst**: Set `retries: 2` in the profile. If a burst job
   is preempted, Snakemake automatically resubmits it.

2. **Burst profile**: Create a separate Snakemake profile for burst jobs:

```yaml
# workflow/profiles/slurm-burst/config.yaml
executor: slurm
default-resources:
  slurm_account: maigan-b
  slurm_partition: hpg-default
  slurm_extra: "--requeue"
retries: 3
```

3. **Run both profiles in parallel** (from separate terminals):
```bash
# Terminal 1: Primary GPU jobs
snakemake --profile profiles/slurm -j 3

# Terminal 2: Burst overflow (shares lock file)
snakemake --profile profiles/slurm-burst -j 10
```

Snakemake's file locking prevents duplicate work across the two instances.

---

## 10. Migration Checklist

### Phase 1: Extract Python logic from bash heredocs

- [ ] Extract the Python code from `02_stitching.sh` into `workflow/scripts/stitch.py`
- [ ] Extract the Python code from `03_deconvolution.sh` into `workflow/scripts/deconvolve.py`
- [ ] Extract the Python code from `04_edf.sh` into `workflow/scripts/edf.py`
- [ ] Remove `wait_for_input()` calls (Snakemake handles dependencies)
- [ ] Remove skip-existing checks (Snakemake handles this)
- [ ] Replace `.complete` marker writes with sentinel file writes
- [ ] Keep all QC/logging from `utils.sh` (call via subprocess or reimplement in Python)

### Phase 2: Create Snakemake infrastructure

- [ ] Write `workflow/Snakefile` with rules for stitch, deconvolve, edf
- [ ] Write `workflow/config.yaml` from project's `slurm/config.sh`
- [ ] Create SLURM profile at `workflow/profiles/slurm/config.yaml`
- [ ] Create conda environment spec at `workflow/envs/kintsugi.yaml`
- [ ] Test with `snakemake -n` (dry run)

### Phase 3: Validate equivalence

- [ ] Run Snakemake pipeline on test dataset (`test_data/mini_project/`)
- [ ] Compare outputs with existing SLURM pipeline outputs
- [ ] Verify QC images are generated
- [ ] Test failure recovery (kill a job, verify Snakemake reruns it)
- [ ] Test `--forcerun` for selective reprocessing

### Phase 4: Advanced features

- [ ] Implement GPU/CPU dual-pool (Option B above)
- [ ] Add burst profile
- [ ] Add registration step (post-EDF, multi-cycle)
- [ ] Add signal isolation step
- [ ] Integrate with `kintsugi slurm submit` CLI

---

## 11. Benefits Over Current submit.sh

| Feature | Current (submit.sh) | Snakemake |
|---------|---------------------|-----------|
| Dependency tracking | SLURM `aftercorr` + `.complete` polling | File-based DAG (automatic) |
| Skip completed work | Custom skip-existing logic per script | Built-in (timestamp-based) |
| Failure recovery | Manual resubmission | `retries:` + automatic |
| Partial re-runs | Re-submit entire step | `--forcerun` specific targets |
| Visualization | `squeue` | `--dag`, `--summary`, `--report` |
| Concurrency control | Manual `--array=%N` calculation | `--jobs` + `--resources` |
| Portability | SLURM-specific | SLURM, LSF, PBS, local, cloud |
| Reproducibility | Script + config.sh | Snakefile is a complete specification |

---

## 12. Compatibility Notes

- **Snakemake version**: Requires Snakemake >= 8.0 with `snakemake-executor-plugin-slurm`
- **Conda**: Snakemake can manage conda environments via `--use-conda`
- **Module loads**: The `module load conda` call in existing scripts should move
  into the Snakemake profile's `slurm_extra` or be handled by the conda
  environment specification
- **NFS latency**: Set `latency-wait: 120` (or higher) to account for NFS
  metadata propagation delays; this replaces the `wait_for_input()` polling
- **Coexistence**: The Snakemake workflow can coexist with the existing
  `submit.sh` pipeline. Both produce files in the same `data/processed/`
  tree. Run one or the other, not both simultaneously.
