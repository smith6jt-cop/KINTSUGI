# Workflow CLAUDE.md

Context for Claude Code when working in the `workflow/` directory. See root `CLAUDE.md` for project overview, CLI usage, and architecture.

## Processing Modes: Notebooks vs SLURM

KINTSUGI has two distinct processing modes with different resource utilization strategies:

| Mode | Context | GPU Policy | CPU Policy | Rationale |
|------|---------|------------|------------|-----------|
| **Notebook** | Interactive, user watching | GPU required, no fallback | Not used | Quality-first; user can intervene if GPU unavailable |
| **SLURM** | Headless batch processing | GPU + CPU concurrent | CPU concurrent | Maximize throughput; use ALL available resources |

**Notebook Processing (Interactive)**:
- GPU is **required** - fails loudly if unavailable
- User is present to diagnose and fix GPU issues
- Quality parameters are non-negotiable (e.g., BaSiC iterations=500)
- See `gpu-quality-priority` skill for details

**SLURM Processing (Headless)**:
- Maximizes throughput by running GPU and CPU jobs **concurrently**
- Uses a **dual-pool architecture** to calculate total concurrent slots from both GPU and CPU resources
- Jobs automatically adapt via `KINTSUGI_DEVICE_MODE` environment variable
- Quality parameters remain unchanged - only the compute device differs

## Resource Pool Calculation (Multi-Account Architecture)

KINTSUGI calculates total concurrent jobs from **all** non-blocked SLURM accounts. Each account contributes **both** GPU slots and CPU slots independently:

- **GPU slots** per account = QOS `gres/gpu` limit (via `sacctmgr`)
- **CPU slots** per account = `floor(0.85 * account_cpus / cpus_per_job)` (85% cap avoids starving other users)
- GPU and CPU partitions are separate resource pools — GPU jobs on `hpg-b200` do NOT consume CPU allocation on `hpg-default`
- **`brusko` account is permanently blocked** — hard-coded in `BLOCKED_ACCOUNTS` frozenset in `hpc.py`

**Total Concurrent** = sum(GPU slots) + sum(CPU slots) across all accounts

**Example calculation** (two accounts, each with GPUs AND CPUs):
| Account | CPUs | Memory | GPUs | GPU Slots | CPU Slots | Calculation |
|---------|------|--------|------|-----------|-----------|-------------|
| `clive` | 104 | 812 GB | 3 | 3 | 11 | GPUs: 3/1; CPUs: floor(0.85*104/8) |
| `maigan` | 80 | 625 GB | 2 | 2 | 8 | GPUs: 2/1; CPUs: floor(0.85*80/8) |
| **Total** | | | **5** | **5** | **19** | **24 concurrent jobs** |

## How Concurrent GPU/CPU Works

1. `detect_multi_account_resources()` / `detect_live_multi_account()` query `sacctmgr show associations` for all user accounts (filtering burst `-b` accounts and blocked accounts)
2. Each account's GPU and CPU slots are calculated independently from its QOS limits
3. GPU jobs are submitted on each account's GPU partition (e.g., `hpg-b200`) with `KINTSUGI_DEVICE_MODE=gpu`
4. CPU jobs are submitted on each account's CPU partition (e.g., `hpg-default`) with `KINTSUGI_DEVICE_MODE=cpu` and guaranteed (non-preemptible) resources
5. Job scripts read `KINTSUGI_DEVICE_MODE` and use appropriate backend (CuPy for GPU, NumPy/SciPy for CPU)
6. CPU jobs get a 5x time multiplier (automatic) but have guaranteed memory allocation

**Important**: `sacctmgr show associations user=USERNAME format=account -n -P` is the correct query format. The older `sacctmgr show user USERNAME format=account` returns empty results.

## Snakemake Workflow (Alternative to submit.sh)

KINTSUGI includes a Snakemake-based workflow as an alternative to `submit.sh` for SLURM batch processing. Both produce files in the same `data/processed/` tree but should not run simultaneously on the same project.

**Setup and usage:**
```bash
kintsugi workflow config .    # Generate workflow/config.yaml + Snakefile
kintsugi workflow check .     # Show live per-account availability
kintsugi workflow run .       # Submit via Snakemake (auto-sets -j)
```

**Key design decisions:**

| Design Choice | Rationale |
|---------------|-----------|
| 3 rules with lambda resources | `ruleorder` does NOT fall back — it always picks the preferred rule, so 6 rules + ruleorder was fundamentally broken |
| Cycle pre-assignment at DAG creation | `_build_cycle_assignment()` distributes cycles proportionally across accounts/modes before any job runs |
| Sentinel files (`.snakemake_complete`) | Stitching/deconvolution produce hundreds of files; declaring every output would create an enormous DAG |
| Per-cycle dependencies | `stitch cyc01 → decon cyc01 → edf cyc01` allows pipelining across cycles |
| No `--resources gpus=N` needed | GPU budget is baked into cycle pre-assignment |

**Config format** (`workflow/config.yaml`):
```yaml
resources:
  accounts:
    - name: clive
      partition_gpu: "hpg-b200,hpg-turin"
      partition_cpu: hpg-default
      gpu_slots: 3
      cpu_slots: 11
    - name: maigan
      partition_gpu: "hpg-b200,hpg-turin"
      partition_cpu: hpg-default
      gpu_slots: 2
      cpu_slots: 8
  total_slots: 24
  cpu_time_multiplier: 5
  cpu_cpus_per_task: 8
```

**Cycle assignment example** (13 cycles, clive 3G+11C, maigan 2G+8C):
- Cycles 1-3 → clive GPU, Cycles 4-5 → maigan GPU
- Cycles 6-16 → clive CPU, Cycles 17+ → maigan CPU (round-robin overflow)

**Lambda resource routing** (per-rule in Snakefile):
```python
resources:
    slurm_partition=lambda wc: _assign(wc)["partition"],
    slurm_account=lambda wc: _assign(wc)["account"],
    gres=lambda wc: "gpu:1" if _is_gpu(wc) else "",
    runtime=lambda wc: base_time * (1 if _is_gpu(wc) else CPU_TIME_MULT),
```

**GPU resource naming**: Use `gres="gpu:1"` (maps to `--gres=gpu:1`). Do NOT use `gpus=1` (maps to `--gpus=1`) or `gpu=1` — both trigger `SLURM_TRES_PER_TASK` conflicts on SLURM >= 24.11.

**SLURM profile** (`workflow/profiles/slurm/config.yaml`):
- `precommand: "module load conda && conda activate KINTSUGI"` — activates the existing conda environment on compute nodes so that conda-installed packages (cupy, torch, etc.) are importable inside srun jobs

**SLURM_TRES_PER_TASK fix** (critical for SLURM >= 24.11):
- SLURM 24.11+ sets `SLURM_TRES_PER_TASK` in GPU job environments, which conflicts with the jobstep plugin's `srun` call
- Error: `srun: fatal: SLURM_TRES_PER_TASK is mutually exclusive with --ntasks-per-gpu`
- Fix: patched `snakemake_executor_plugin_slurm_jobstep/__init__.py` to unset this variable in `__post_init__()`: `os.environ.pop("SLURM_TRES_PER_TASK", None)`
- **This patch must be re-applied after any pip upgrade of `snakemake-executor-plugin-slurm-jobstep`**

**Cycle directory resolution**: `_resolve_raw_cycle_dir()` handles `cyc001_reg001_*`, `cyc001`, `cyc01`, and `Cyc01` naming conventions at DAG creation time. Accepts int or str (Snakemake CLI `--config` passes strings).

**`workflow config` behavior**: Always overwrites the Snakefile and profiles (so pipeline logic and SLURM precommand updates propagate); only copies scripts if they don't already exist.

**Per-channel skip-existing checks** (inside wrapper scripts):
- Snakemake controls the DAG at cycle level via sentinel files. If a sentinel is missing, the entire cycle reruns.
- Wrapper scripts (`stitch.py`, `deconvolve.py`, `edf.py`) add per-channel skip-existing checks inside the job to avoid re-processing completed channels when a job was interrupted mid-cycle.
- A channel is only skipped when ALL expected output files exist (partially-complete channels are fully reprocessed).
- `stitch.py`: Checks all z-plane TIFs + `result_df.pkl` for CH1
- `deconvolve.py`: Compares decon TIF count against stitched input count
- `edf.py`: Checks if marker-named output file exists
- Sentinel files include `skipped=N` for visibility in logs

**What Snakemake replaces vs keeps:**
- Replaces: `submit.sh` orchestration, dependency wiring, `.complete` polling, `--array` limit calculation
- Keeps: Python processing code in wrapper scripts (`workflow/scripts/*.py`), `KINTSUGI_DEVICE_MODE` environment variable pattern, per-channel skip-existing inside wrapper scripts
- Coexists with: Existing `slurm/` job scripts (run one system or the other, not both)

## Batch Processing (Multi-Dataset)

KINTSUGI supports automated batch processing of multiple CODEX datasets. See `KINTSUGI_Batch_Processing_Guide.docx` for the complete 8-phase workflow.

**Key scripts** (in `/blue/maigan/smith6jt/`):
- `data_discovery.py` - Scans orange storage, builds `dataset_manifest.csv` from experiment.json files
- `setup_all_projects.sh` - Creates KINTSUGI projects for all datasets in manifest
- `stage_datasets.sh` - Submits SLURM rsync jobs to copy raw data from orange to blue
- `configure_workflows.sh` - Generates Snakemake configs for each project

**Storage architecture** (HiPerGator):
- `/orange/maigan/` - Long-term storage for raw CODEX data (~9.7 TB, 34 datasets)
- `/blue/maigan/` - Processing workspace (~17 TB available). Raw data must be copied here.
- Peak blue usage per dataset: ~3x raw size (raw + stitched + deconvolved intermediates)

**Rsync staging pattern** (preserves cycle directory structure):
```bash
rsync -av '${orange_path}/' '${RAW_DIR}/' \
  --include='*/' --include='*.tif' --exclude='*'
```
**CRITICAL**: Never use bash glob `*/` on rsync source — it flattens cycle directories and overwrites files across cycles. Use the source path with trailing `/` and let rsync handle the tree.

**Batch processing phases**:
1. Discovery → 2. Setup + Staging → 3. Single-cycle validation → 4. SLURM batch (stitch→decon→EDF) → 5. Registration → 6. Cleanup intermediates → 7. Signal isolation (local GPU PCs) → 8. Segmentation (local GPU PCs)
