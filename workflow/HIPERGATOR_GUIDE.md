# Running KINTSUGI on HiPerGator with Snakemake

Complete step-by-step protocol for running the KINTSUGI image processing
pipeline on UF HiPerGator using the Snakemake workflow engine.

This replaces the legacy `submit.sh` + `burst_monitor.sh` orchestration with
Snakemake's native dependency management, automatic failure recovery, and
per-file tracking.

---

## Table of Contents

1. [One-Time Setup: Install Snakemake](#step-1-one-time-setup-install-snakemake-in-kintsugi-environment)
2. [Initialize Your Project](#step-2-initialize-your-project)
3. [Prepare Raw Data](#step-3-prepare-raw-data)
4. [Create Channel Names File](#step-4-create-channel-names-file)
5. [Configure Experiment Metadata](#step-5-configure-experiment-metadata)
6. [Generate Workflow Configuration](#step-6-generate-workflow-configuration)
7. [Review and Edit Configuration](#step-7-review-and-edit-configuration)
8. [Dry Run (Preview)](#step-8-dry-run-preview)
9. [Submit to SLURM](#step-9-submit-to-slurm)
10. [Monitor Progress](#step-10-monitor-progress)
11. [Review QC Output](#step-11-review-qc-output)
12. [Reprocessing and Recovery](#step-12-reprocessing-and-recovery)
13. [Advanced: GPU/CPU Dual-Pool](#step-13-advanced-gpucpu-dual-pool)
14. [Advanced: Burst Resources](#step-14-advanced-burst-resources)
15. [Troubleshooting](#troubleshooting)

---

## Step 1: One-Time Setup — Install Snakemake in KINTSUGI Environment

Snakemake and its SLURM executor plugin must be installed **once** in your
KINTSUGI conda environment. This only needs to be done on a fresh environment
or after recreating it.

### 1a. Load conda and activate your environment

```bash
module load conda
conda activate kintsugi
```

### 1b. Install Snakemake and the SLURM executor

Snakemake 8+ uses a plugin architecture. You need both the core package and
the SLURM executor plugin:

```bash
pip install "snakemake>=8.0" snakemake-executor-plugin-slurm
```

> **Conda alternative** (if pip causes conflicts):
> ```bash
> conda install -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm
> ```

### 1c. Verify the installation

```bash
snakemake --version
# Should print 8.x.x or higher

snakemake --list-executor-plugins
# Should list "slurm" among the available executors
```

### 1d. Install KINTSUGI (if not already installed)

```bash
cd /path/to/KINTSUGI
pip install -e ".[gpu]"
```

This installs KINTSUGI in editable mode with GPU acceleration (CuPy + PyTorch).
Verify:

```bash
kintsugi check
# Should show core dependencies as OK
```

### 1e. Verify GPU access (optional but recommended)

```bash
# On a GPU node (via srun or an interactive job):
srun --partition=hpg-b200 --gpus=1 --mem=8gb --time=00:10:00 --pty bash
module load conda && conda activate kintsugi
python -c "import cupy; print(f'CuPy OK: {cupy.cuda.runtime.getDeviceCount()} GPUs')"
exit
```

---

## Step 2: Initialize Your Project

Create the project directory structure with `kintsugi init`:

```bash
module load conda
conda activate kintsugi

kintsugi init /blue/maigan/your_user/projects/my_experiment \
    --name "My Experiment" \
    --tile-rows 5 --tile-cols 5 \
    --xy-pixel-size 377 --z-step-size 1500 \
    --numerical-aperture 0.75 --tissue-ri 1.44
```

This creates:

```
my_experiment/
├── data/
│   ├── raw/           ← Put raw images here
│   └── processed/     ← Outputs go here
├── meta/
│   ├── experiment.json    ← Microscope parameters (auto-generated)
│   └── CHANNELNAMES.txt   ← You create this (Step 4)
├── notebooks/
├── configs/
└── kintsugi_project.json
```

> **Tip**: If you already have a project (from the legacy SLURM pipeline), skip
> this step. The Snakemake workflow can be added to any existing project.

---

## Step 3: Prepare Raw Data

Copy cycle folders into `data/raw/`:

```bash
cd /blue/maigan/your_user/projects/my_experiment

# Copy from your data source
cp -r /path/to/source/cyc001 data/raw/
cp -r /path/to/source/cyc002 data/raw/
# ... repeat for all cycles
```

**Supported naming conventions** (both work automatically):

| Format | Example |
|--------|---------|
| Short | `cyc01/`, `cyc02/` |
| Zero-padded | `cyc001/`, `cyc002/` |
| Long-form with markers | `cyc001_DAPI_Blank_Blank_Blank/` |

**Expected tile naming**: `*_Z{nnn}_CH{n}.tif` (e.g., `tile_001_Z001_CH1.tif`)

**Verify** your data looks correct:

```bash
kintsugi scan data/raw/
```

---

## Step 4: Create Channel Names File

Create `meta/CHANNELNAMES.txt` with your marker names:

```bash
nano meta/CHANNELNAMES.txt
```

**Simple list format** (one marker per line, cycles sequentially):

```
DAPI-01
Blank
Blank
Blank
DAPI-02
CD31
CD8
CD45
DAPI-03
CD3e
Ki67
PanCK
```

**Alternative cycle-prefixed format**:

```
1: DAPI, Blank, Blank, Blank
2: DAPI, CD31, CD8, CD45
3: DAPI, CD3e, Ki67, PanCK
```

---

## Step 5: Configure Experiment Metadata

Edit `meta/experiment.json` to match your microscope acquisition parameters:

```bash
nano meta/experiment.json
```

**Key parameters to verify**:

| Parameter | Description | Default | Check Against |
|-----------|-------------|---------|---------------|
| `tile_rows` | Tile grid rows | 5 | Your acquisition settings |
| `tile_cols` | Tile grid columns | 5 | Your acquisition settings |
| `tile_overlap` | Overlap fraction | 0.1 | Typically 10% |
| `xy_pixel_size` | XY voxel size (nm) | 377.0 | Objective/camera specs |
| `z_step_size` | Z step size (nm) | 1500.0 | Acquisition Z interval |
| `numerical_aperture` | Objective NA | 0.75 | Engraved on objective |
| `tissue_refractive_index` | Tissue RI | 1.44 | Typical for FFPE tissue |
| `n_cycles` | Number of cycles | Auto-detected | Count of raw cycle dirs |
| `n_zplanes` | Z-planes per stack | Auto-detected | Count from raw tiles |
| `wavelengths` | Ex/Em per channel | See defaults | Filter cube specs |

**Example** `experiment.json`:

```json
{
  "tile_rows": 5,
  "tile_cols": 5,
  "tile_overlap": 0.1,
  "xy_pixel_size": 377.0,
  "z_step_size": 1500.0,
  "numerical_aperture": 0.75,
  "tissue_refractive_index": 1.44,
  "wavelengths": {
    "1": [358.0, 461.0],
    "2": [753.0, 775.0],
    "3": [560.0, 575.0],
    "4": [648.0, 668.0]
  },
  "channels_per_cycle": 4,
  "n_cycles": 7,
  "n_zplanes": 15
}
```

---

## Step 6: Generate Workflow Configuration

This is the key step that replaces the old `slurm/config.sh`. It reads your
`experiment.json` and generates a unified Snakemake configuration.

```bash
cd /blue/maigan/your_user/projects/my_experiment
module load conda
conda activate kintsugi

kintsugi workflow config .
```

This creates the `workflow/` directory:

```
my_experiment/
└── workflow/
    ├── Snakefile              ← Pipeline rules (stitch → decon → edf)
    ├── config.yaml            ← All parameters in one file
    ├── scripts/
    │   ├── stitch.py          ← BaSiC correction + stitching
    │   ├── deconvolve.py      ← Richardson-Lucy deconvolution
    │   └── edf.py             ← Extended depth of focus
    ├── profiles/
    │   └── slurm/
    │       └── config.yaml    ← SLURM executor settings
    └── envs/
        └── kintsugi.yaml      ← Conda environment spec
```

> **Existing project with slurm/config.sh**: The command automatically extracts
> your SLURM account and partition from the existing `slurm/config.sh`.

> **Preview without writing**: Use `kintsugi workflow config . --print-only`
> to see what would be generated.

---

## Step 7: Review and Edit Configuration

### 7a. Workflow config (`workflow/config.yaml`)

Open and verify the auto-generated configuration:

```bash
nano workflow/config.yaml
```

**Sections to check**:

```yaml
# Paths — should match your project
project_dir: "/blue/maigan/your_user/projects/my_experiment"
kintsugi_dir: "/blue/maigan/your_user/KINTSUGI"

# Processing scope — which cycles/channels to process
cycles: [1, 2, 3, 4, 5, 6, 7]
channels: [1, 2, 3, 4]

# Tile grid — must match your acquisition
tile_rows: 5
tile_cols: 5
tile_overlap: 0.1

# Stitching parameters (generally fine as defaults)
ncc_threshold: 0.078
pou: 0.5
blend_sigma: 15.0

# BaSiC illumination correction (generally fine as defaults)
basic_flatfield_min: 0.1
basic_max_iterations: 500
basic_optimization_tolerance: 1.0e-6

# EDF parameters (CLIJ2 defaults — generally fine)
edf:
  radius_x: 2
  radius_y: 2
  sigma: 10.0

# SLURM resources — adjust for your allocation
resources:
  account: "maigan"           # ← Your SLURM account
  partition_gpu: "hpg-b200"   # ← GPU partition
  partition_cpu: "hpg-default"
  mem_stitch: 48000           # MB per stitching job
  mem_decon: 48000            # MB per deconvolution job
  mem_edf: 16000              # MB per EDF job
  time_stitch: 240            # minutes
  time_decon: 240             # minutes
  time_edf: 60                # minutes
```

### 7b. SLURM profile (`workflow/profiles/slurm/config.yaml`)

```bash
nano workflow/profiles/slurm/config.yaml
```

```yaml
executor: slurm

default-resources:
  slurm_account: maigan         # ← Must match your account
  slurm_partition: hpg-b200     # ← Must match available partition
  mem_mb: 48000
  runtime: 240
  cpus_per_task: 4
  gpus: 1

jobs: 8            # Max concurrent SLURM jobs
latency-wait: 120  # Seconds to wait for NFS file propagation
retries: 2         # Auto-retry failed jobs
keep-going: true   # Continue independent jobs on failure
use-conda: true    # Use conda environments
```

**Important**: Set `slurm_account` and `slurm_partition` to values valid for
your HiPerGator allocation. Check available accounts with:

```bash
sacctmgr show user $(whoami) format=Account%30 --noheader
```

Check available partitions with:

```bash
sinfo -s --format="%P %a %l %D %G"
```

---

## Step 8: Dry Run (Preview)

**Always preview before submitting**. This validates the DAG, checks for
missing inputs, and shows exactly what will run without executing anything.

### Option A: Via the KINTSUGI CLI

```bash
cd /blue/maigan/your_user/projects/my_experiment
module load conda && conda activate kintsugi

kintsugi workflow run . --dry-run
```

### Option B: Direct Snakemake command

```bash
cd workflow/
snakemake --profile profiles/slurm -n
```

**Expected output** (example with 7 cycles):

```
Job stats:
job          count
---------  -------
all              1
deconvolve       7
edf              7
stitch           7
total           22

Would execute 21 jobs (the "all" rule is local).
```

**Verify**:
- [ ] Correct number of cycles detected
- [ ] All three steps listed (stitch, deconvolve, edf)
- [ ] No "MissingInputException" errors

### Visualize the dependency graph (optional)

```bash
cd workflow/
snakemake --dag | dot -Tpng > dag.png
```

Open `dag.png` to see the visual DAG of all jobs and their dependencies.

---

## Step 9: Submit to SLURM

### Full pipeline

```bash
cd /blue/maigan/your_user/projects/my_experiment
module load conda && conda activate kintsugi

kintsugi workflow run .
```

Or equivalently, using Snakemake directly:

```bash
cd workflow/
snakemake --profile profiles/slurm
```

### Specific cycles only

```bash
# Cycles 1-3
kintsugi workflow run . --cycles 1-3

# Specific cycles
kintsugi workflow run . --cycles 1,3,5
```

### Specific steps only

```bash
cd workflow/

# Only stitching
snakemake --profile profiles/slurm --allowed-rules stitch

# Deconvolution + EDF (skip stitching, assumes already done)
snakemake --profile profiles/slurm --allowed-rules deconvolve edf
```

### Override job count

```bash
# Allow up to 16 concurrent SLURM jobs
kintsugi workflow run . -j 16
```

### What happens behind the scenes

1. Snakemake reads `config.yaml` and builds a DAG of all jobs
2. For each cycle: stitch depends on raw data, deconvolve depends on stitch sentinel, edf depends on decon sentinel
3. Snakemake submits SLURM jobs respecting the dependency order
4. Each job runs in its own SLURM allocation with the requested resources
5. Sentinel files (`.snakemake_complete`) are written after each step
6. If a job fails, Snakemake automatically retries up to `retries` times
7. Completed cycles are **never re-processed** (Snakemake checks for existing outputs)

---

## Step 10: Monitor Progress

### SLURM job queue

```bash
squeue -u $USER
```

### Snakemake progress

If you ran Snakemake from a terminal (not background), progress is printed
live. For background runs, check the Snakemake log:

```bash
# Snakemake creates a hidden log directory
ls -la workflow/.snakemake/log/
tail -f workflow/.snakemake/log/*.snakemake.log
```

### Per-job logs

Each job writes logs to the project's log directory:

```bash
ls workflow/../slurm/logs/snakemake/
# Files: stitch_cyc01.log, decon_cyc02.log, edf_cyc03.log, etc.
```

### Check completion status

```bash
# List all sentinel files (one per completed step per cycle)
find data/processed/ -name ".snakemake_complete" -exec echo {} \; -exec cat {} \;
```

**Example sentinel content**:
```
stage=stitch
cycle=1
completed=2026-02-06T14:30:00
channels=1-4
zplanes=15
duration_minutes=23.4
```

### Summary of outputs

```bash
# Count output files per stage
echo "Stitched:"; ls data/processed/stitched/cyc*/CH*/*.tif 2>/dev/null | wc -l
echo "Deconvolved:"; ls data/processed/deconvolved/cyc*/CH*/*.tif 2>/dev/null | wc -l
echo "EDF:"; ls data/processed/edf/cyc*/*_edf.tif 2>/dev/null | wc -l
```

---

## Step 11: Review QC Output

QC images are saved alongside the Snakemake logs:

```bash
ls slurm/logs/snakemake/*.png
# Files: cyc01_CH1_decon.png, cyc01_CH1_edf.png, etc.
```

**What to check**:

| Stage | QC Shows | Look For |
|-------|----------|----------|
| Stitching | (sentinel metadata) | Correct tile count, no errors in log |
| Deconvolution | Middle z-slice | Clear signal, no ringing artifacts |
| EDF | 2D projection | Sharp features, <5% zero pixels |

---

## Step 12: Reprocessing and Recovery

### Automatic recovery

If a job fails (OOM, timeout, node failure), Snakemake automatically retries
it up to `retries` times (default: 2). No manual intervention needed.

### Force reprocessing of a specific cycle

```bash
cd workflow/

# Force re-stitch cycle 3 (and everything downstream)
snakemake --profile profiles/slurm \
    --forcerun data/processed/stitched/cyc03/.snakemake_complete
```

### Force reprocessing of a specific step

```bash
# Force re-run all deconvolution (all cycles)
snakemake --profile profiles/slurm --forcerun deconvolve
```

### Delete outputs and reprocess from scratch

```bash
# Remove all processed data for a cycle
rm -rf data/processed/stitched/cyc03/
rm -rf data/processed/deconvolved/cyc03/
rm -rf data/processed/edf/cyc03/

# Snakemake will automatically detect missing outputs and reprocess
cd workflow/ && snakemake --profile profiles/slurm
```

### Resume after interruption

If the submission terminal was closed or the login node timed out, simply
re-run the same command. Snakemake checks for existing outputs and only
submits jobs for missing files:

```bash
cd workflow/ && snakemake --profile profiles/slurm
```

---

## Step 13: Advanced — GPU/CPU Dual-Pool

The legacy `submit.sh` ran GPU and CPU jobs concurrently using a dual-pool
architecture. With Snakemake, the simplest approach is:

### Option A: GPU-only (recommended for most users)

This is the default. All jobs request `gpus: 1` and fall back to CPU
automatically if CUDA initialization fails. Snakemake's `jobs: 8` limits
concurrency.

### Option B: Explicit GPU budget with Snakemake resource constraints

Limit GPU usage globally and let Snakemake schedule optimally:

```bash
# Allow 3 GPU jobs concurrently (out of 8 total)
cd workflow/
snakemake --profile profiles/slurm -j 8 --resources gpus=3
```

When 3 GPU jobs are running, Snakemake holds remaining jobs until a GPU slot
frees up. This is simpler and more reliable than the old dual-pool + burst
monitor approach.

### Option C: GPU + CPU rule variants (advanced)

For maximum throughput, create CPU-fallback variants in the Snakefile.
See `slurm/SNAKEMAKE_PROTOCOL.md` section 8 (Option B) for details.

---

## Step 14: Advanced — Burst Resources

HiPerGator burst QOS provides access to idle cluster resources, but jobs are
preemptible (can be killed when the resource owner needs them back).

### Create a burst profile

```bash
mkdir -p workflow/profiles/slurm-burst
cat > workflow/profiles/slurm-burst/config.yaml << 'EOF'
executor: slurm

default-resources:
  slurm_account: maigan-b
  slurm_partition: hpg-default
  slurm_extra: "--requeue"
  mem_mb: 48000
  runtime: 240
  cpus_per_task: 8
  gpus: 0

jobs: 10
latency-wait: 120
retries: 3
keep-going: true
use-conda: true
EOF
```

### Run primary + burst in parallel

Open two terminals:

```bash
# Terminal 1: Primary GPU jobs (guaranteed resources)
cd workflow/
snakemake --profile profiles/slurm -j 3

# Terminal 2: Burst overflow (preemptible, CPU fallback)
cd workflow/
snakemake --profile profiles/slurm-burst -j 10
```

Snakemake's file locking prevents duplicate work. The burst instance will
automatically pick up cycles not already being processed by the primary
instance.

> **Note**: If a burst job is preempted, `--requeue` tells SLURM to requeue
> it, and `retries: 3` tells Snakemake to retry up to 3 times.

---

## Troubleshooting

### "snakemake: command not found"

```bash
module load conda
conda activate kintsugi
pip install "snakemake>=8.0" snakemake-executor-plugin-slurm
```

### "No workflow/config.yaml found"

Run `kintsugi workflow config .` from your project directory first.

### Jobs pending indefinitely (QOS limits)

Check your QOS limits:

```bash
sacctmgr show qos format=Name%20,GrpTRES%30,MaxTRES%30
```

Reduce `jobs:` in `profiles/slurm/config.yaml` to fit within your allocation.

### "CUDA initialization failed" in job logs

The scripts automatically fall back to CPU mode. This is expected on CPU-only
nodes. To force GPU nodes, verify your profile's `slurm_partition` is a GPU
partition (e.g., `hpg-b200`, `hpg-turin`).

### Jobs fail with OOM (Out of Memory)

Increase memory in `workflow/config.yaml`:

```yaml
resources:
  mem_stitch: 96000   # Doubled from 48GB
  mem_decon: 96000
```

### NFS latency — "Missing output files after job completion"

Increase the wait time in `profiles/slurm/config.yaml`:

```yaml
latency-wait: 300  # Increased from 120 to 300 seconds
```

### Stitch model not found for CH2+

Channel 1 computes the stitching model that all other channels use. If CH1
fails, subsequent channels will fail with "No stitch model at ...". Check the
CH1 stitching log first.

### "MissingInputException: Missing input files for rule stitch"

Your raw data directory is missing or has unexpected naming. Verify:

```bash
ls data/raw/
# Should show cyc01/, cyc02/, etc.
```

### How to check which jobs are complete vs pending

```bash
# Dry run shows only remaining work
cd workflow/ && snakemake --profile profiles/slurm -n
```

### Comparing to the legacy SLURM pipeline

| Feature | Legacy (`submit.sh`) | Snakemake |
|---------|---------------------|-----------|
| Config | `slurm/config.sh` (~60 vars) | `workflow/config.yaml` (one file) |
| Dependencies | `aftercorr` + `.complete` polling | File-based DAG (automatic) |
| Skip completed | Custom logic per script | Built-in (checks output timestamps) |
| Failure recovery | Manual resubmission | `retries: 2` (automatic) |
| Partial re-run | Not supported | `--forcerun` specific targets |
| Concurrency | Manual array calculations | `jobs:` + `--resources` |
| Visualization | `squeue` only | `--dag`, `--summary`, `--report` |
| Burst monitor | Separate daemon (`burst_monitor.sh`) | Second Snakemake profile |
| GPU/CPU pools | 12+ SLURM array jobs | 1 job per step per cycle |

---

## Quick Reference Card

```bash
# === SETUP (one-time) ===
module load conda && conda activate kintsugi
pip install "snakemake>=8.0" snakemake-executor-plugin-slurm

# === PER-PROJECT SETUP ===
cd /path/to/my_project
kintsugi workflow config .                # Generate workflow/

# === RUNNING ===
kintsugi workflow run . --dry-run         # Preview
kintsugi workflow run .                   # Full pipeline via SLURM
kintsugi workflow run . --cycles 1-3      # Specific cycles
kintsugi workflow run . --local --cores 4 # Local (no SLURM)

# === DIRECT SNAKEMAKE (more control) ===
cd workflow/
snakemake --profile profiles/slurm -n     # Dry run
snakemake --profile profiles/slurm        # Submit
snakemake --profile profiles/slurm \
    --forcerun stitch                     # Force re-stitch all
snakemake --report report.html            # Generate HTML report

# === MONITORING ===
squeue -u $USER                           # Job queue
tail -f slurm/logs/snakemake/*.log        # Live logs
find data/processed/ -name ".snakemake_complete" | wc -l  # Completed steps
```
