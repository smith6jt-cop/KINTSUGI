# Workflow CLAUDE.md

Context for Claude Code when working in the `workflow/` directory. See root `CLAUDE.md` for project overview, CLI usage, and architecture.

## Processing Modes: Notebooks vs SLURM

KINTSUGI has two distinct processing modes with different resource utilization strategies:

| Mode | Context | GPU Policy | CPU Policy | Rationale |
|------|---------|------------|------------|-----------|
| **Notebook** | Interactive, user watching | GPU required, no fallback | Not used | Quality-first; user can intervene if GPU unavailable |
| **SLURM** | Headless batch processing | GPU-only, overflow queues | Not used (too slow) | GPU ~22 min vs CPU ~282 min per cycle (stitch 25x, decon 5x, EDF 15x) |

**Notebook Processing (Interactive)**:
- GPU is **required** - fails loudly if unavailable
- User is present to diagnose and fix GPU issues
- Quality parameters are non-negotiable (e.g., BaSiC iterations=500)
- See `gpu-quality-priority` skill for details

**SLURM Processing (Headless — GPU-Only Scheduling)**:
- **All cycles route through GPU slots** — no CPU fallback
- Overflow cycles queue in SLURM until a GPU slot frees up
- Measured speedups (CX_19-003, 9x7, 4ch, 20z): stitch 25x, decon 5x, EDF 15x (~13x full cycle)
- Snakemake `-j` is set to total GPU slots (not total slots) to limit concurrent submissions
- Per-cycle pipeline (`stitch→decon→edf`) enables cycles to flow through freed GPU slots
- See `gpu-only-scheduling` skill for the performance data behind this decision

## Resource Pool Calculation (Multi-Account Architecture)

KINTSUGI calculates available GPU slots from **all** non-blocked SLURM accounts:

- **GPU slots** per account = QOS `gres/gpu` limit (via `sacctmgr`)
- **CPU slots** are still detected but **not used for scheduling** (CPU is 5-25x slower per step — see `gpu-only-scheduling` skill)
- **`brusko` account is permanently blocked** — hard-coded in `BLOCKED_ACCOUNTS` frozenset in `hpc.py`

**Total Concurrent** = sum(GPU slots) across all accounts (CPU pool is informational only)

**Example calculation** (two accounts):
| Account | GPUs | GPU Slots | CPU Slots (unused) |
|---------|------|-----------|-------------------|
| `clive` | 3 | 3 | 11 |
| `maigan` | 2 | 2 | 8 |
| **Total** | **5** | **5** | 19 (not used) |

## How GPU-Only Scheduling Works

1. `detect_multi_account_resources()` / `detect_live_multi_account()` query `sacctmgr show associations` for all user accounts (filtering burst `-b` accounts and blocked accounts)
2. Each account's GPU slots are calculated from its QOS limits
3. `_build_cycle_assignment()` pre-assigns ALL cycles to GPU, round-robin across accounts proportional to GPU slot counts
4. Overflow cycles (beyond total GPU slots) still get GPU assignments — they queue in SLURM until a GPU slot frees up
5. Snakemake `-j` = total GPU slots, so concurrent SLURM submissions match available GPUs
6. Per-cycle pipeline (`stitch→decon→edf`) means freed GPU slots are immediately used by the next stage

**Why no CPU fallback**: Measured per-step (CX_19-003, 9x7, 4ch, 20z): stitch GPU ~8 min vs CPU ~200 min (25x), decon ~12 vs ~60 min (5x), EDF ~1.5 vs ~22 min (15x). Full cycle ~22 min GPU vs ~282 min CPU. Queuing for GPU is always faster.

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
| 3 per-cycle rules with lambda resources | `ruleorder` does NOT fall back — it always picks the preferred rule, so 6 rules + ruleorder was fundamentally broken |
| GPU-only cycle assignment | `_build_cycle_assignment()` routes ALL cycles to GPU — CPU is 5-25x slower per step (~13x full cycle), queuing for GPU is faster |
| Sentinel files (`.snakemake_complete`) | Stitching/deconvolution produce hundreds of files; declaring every output would create an enormous DAG |
| Per-cycle dependencies | `stitch cyc01 → decon cyc01 → edf cyc01` allows pipelining across cycles |
| No `--resources gpus=N` needed | GPU budget is baked into cycle pre-assignment |
| Registration as aggregate rule (no `{cycle}` wildcard) | Registration aligns ALL cycles at once — needs all EDF outputs before starting |
| Static account assignment for registration & QC | `_registration_assignment()` picks first GPU account; reused by QC rules |
| Aggregate QC rules (3 stages) | QC runs once across ALL cycles for cross-cycle comparison heatmaps |
| Dynamic worker counts | Wrapper scripts read `cpus_per_task` from Snakemake resources instead of hardcoded 4 |

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

**Cycle assignment example** (9 cycles, clive 3G, maigan 2G = 5 GPU slots):
- Cycles 1-3 → clive GPU, Cycles 4-5 → maigan GPU (fills all 5 GPU slots)
- Cycles 6-9 → overflow: round-robin clive/maigan GPU (queue in SLURM until a GPU frees up)
- No CPU assignments — all cycles run on GPU

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

**QC Report Rules** (aggregate — added 2026-02-13):

The workflow includes 3 aggregate QC rules that produce rich QC reports (summary heatmaps, z-plane profiles, cross-stage comparison PDFs) after each processing stage completes. These replicate the Notebook 2 QC outputs via `Kprocess.py` functions.

| Rule | Depends On | Kprocess Function | Output |
|------|-----------|-------------------|--------|
| `qc_stitch` | All stitch sentinels (all cycles) | `run_stitched_qc()` | `qc_plots/stitched_summary_heatmaps.pdf` + z-profiles for ALL cycles |
| `qc_decon` | All decon sentinels (all cycles) | `run_decon_qc()` | `qc_plots/deconvolved_summary_heatmaps.pdf` + z-profiles for ALL cycles |
| `qc_edf` | All EDF sentinels (all cycles) | `run_edf_qc()` | `qc_plots/edf_summary_heatmaps.pdf` |

Key design decisions:
- **Aggregate (not per-cycle)**: QC runs once across ALL cycles for cross-cycle comparison heatmaps
- **Single dispatch script**: `scripts/qc_report.py` dispatches to the correct Kprocess function based on `snakemake.params.stage`
- **Cross-stage comparison**: Each stage loads the previous stage's cached stats pickle (`cache/stitch_stats.pkl` → `decon_stats.pkl` → `edf_stats.pkl`)
- **Headless matplotlib**: `matplotlib.use("Agg")` is set before any Kprocess import to prevent display issues on compute nodes
- **GPU with CPU fallback**: Same pattern as processing scripts — tries CuPy, falls back to CPU
- **Resources**: 64 GB RAM, 120 min, 1 GPU on first account (reuses `_REG_ASSIGN` from registration)
- **`rule all` includes QC**: `all_qc_sentinels()` is included alongside registration sentinel in the default target
- **`rule qc` standalone target**: Run `snakemake qc` to generate only QC reports without reprocessing

```bash
snakemake -n qc           # Dry run — show QC DAG
snakemake qc --profile profiles/slurm -j 1  # Run QC only
```

**Config resources** (`workflow/config.yaml`):
```yaml
resources:
  mem_qc: 64000    # 64 GB for QC report generation
  time_qc: 120     # 2 hours
```

**Registration Rule** (aggregate — added 2026-02-13):

Registration (multi-cycle alignment via VALIS) is Rule 4 in the pipeline. Unlike stitch/decon/edf which process one cycle each, registration processes ALL cycles in a single SLURM job.

**Full pipeline**: `stitch → deconvolve → edf` (per-cycle, pipelined) → `registration` (all cycles) + `qc_stitch/qc_decon/qc_edf` (aggregate after each stage)

| Aspect | Details |
|--------|---------|
| **Input** | All EDF sentinels (waits for every cycle to complete) |
| **Output** | `data/processed/registered/.snakemake_complete` sentinel |
| **Account** | Static via `_registration_assignment()` — first GPU account preferred |
| **GPU usage** | VggFD feature detector (GPU); falls back to OrbFD on CPU |
| **Script** | `scripts/registration.py` |
| **Resources** | 64 GB RAM, 120 min (GPU) / 600 min (CPU), 4 CPUs |

Key design decisions:
- **No `{cycle}` wildcard**: Registration is a single job, not per-cycle. Uses static `_REG_ASSIGN` dict instead of `_assign(wc)` lambdas
- **`_registration_assignment()`**: Prefers first account with `gpu_slots > 0` for VggFD; falls back to CPU with OrbFD
- **Single-cycle handling**: If `len(CYCLES) < 2`, copies EDF → registered without alignment (no VALIS needed)
- **Error fallback**: On VALIS failure, copies EDF images unchanged (matches notebook behavior)
- **Skip-existing**: All-or-nothing — checks if ALL registered files exist for ALL cycles before skipping
- **Slide matching**: Builds `cycle_to_slide_name` dict from DAPI filenames for robust `slide_dict` lookup

**Registration config** (`workflow/config.yaml`):
```yaml
registration:
  reference_cycle: 1          # Cycle used as fixed reference
  reference_channel: 1        # Channel for feature detection (usually DAPI)
  max_image_dim: 2048         # Max dimension for registration computation
  rigid_only: false           # Set true to skip non-rigid registration
  feature_detector: "VggFD"   # "VggFD" (GPU, better) or "OrbFD" (CPU-only)

resources:
  mem_registration: 64000     # 64 GB
  time_registration: 120      # 2 hours
  cpu_mem_registration: 64000 # Same for CPU (loads one image per cycle, not z-stacks)
```

**Per-channel skip-existing**: Wrapper scripts (`stitch.py`, `deconvolve.py`, `edf.py`) check per-channel completion inside jobs to resume interrupted cycles. Sentinel files include `skipped=N` in logs. Dynamic worker counts read from `snakemake.resources.cpus_per_task`.

**Snakemake replaces** `submit.sh` orchestration/dependency wiring; **keeps** Python wrapper scripts (`workflow/scripts/*.py`). Run one system or the other, not both.

## Batch Processing (Multi-Dataset)

KINTSUGI supports automated batch processing of multiple CODEX datasets. See `/blue/maigan/smith6jt/README.md` for the complete pipeline reference and `KINTSUGI_Batch_Processing_Guide.docx` for the 8-phase workflow.

**Key scripts** (in `/blue/maigan/smith6jt/`):
- `dataset_manifest.csv` - Central registry of all 47 datasets (34 spleen/LN + 13 thymus; source paths, parameters, sizes)
- `setup_all_projects.sh` - Creates KINTSUGI projects for all datasets in manifest
- `configure_all_workflows.sh` - Generates Snakemake configs for each project
- `stage_datasets.sh` - Submits SLURM rsync jobs to copy raw data from orange to blue (wave-based)
- `stage_datasets_globus.py` - Globus-based staging for thymus datasets from PATH lab SMB share
- `thymus_manifest.csv` - Standalone manifest for 13 thymus datasets (subset of main manifest)
- `run_all_workflows.sh` - Runs Snakemake for all staged datasets sequentially
- `cleanup_datasets.sh` - Verifies EDF outputs + QC sentinels, prompts for QC review, deletes intermediates and raw data (`--force` skips prompt)
- `pipeline_status.sh` - Shows current state of every dataset

**Pipeline lifecycle:**
```
1. setup_all_projects.sh      Create project dirs + copy metadata from orange
2. configure_all_workflows.sh Generate Snakemake configs for all projects
3. stage_datasets.sh 5        Stage next wave of raw data (orange → blue)
4. run_all_workflows.sh       Process all staged datasets (stitch → decon → EDF → registration + QC)
5. cleanup_datasets.sh        Verify outputs + QC, prompt review, delete intermediates + raw
6. Repeat 3-5 for next wave
```

**Running batch workflows:**
```bash
# Process all staged datasets (run inside tmux/screen!)
tmux new -s batch
bash run_all_workflows.sh               # all staged datasets, sequential
bash run_all_workflows.sh --dry-run     # preview without executing
bash run_all_workflows.sh --dataset CX_19-002_lymph-node_R1  # single dataset
```

**Why sequential**: All datasets share 5 GPU slots — parallel Snakemake instances cause contention. Run inside tmux. Re-runs auto-skip completed datasets via sentinel files + per-channel skip-existing.

**Storage**: `/orange/maigan/` (long-term, ~9.7 TB spleen/LN), PATH lab SMB (thymus, ~2.4 TB via Globus), `/blue/maigan/` (processing workspace). Peak blue usage per dataset: ~3x raw size.

**CRITICAL rsync rule**: Never use bash glob `*/` on rsync source — it flattens cycle directories. Use trailing `/` and let rsync handle the tree: `rsync -av '${orange_path}/' '${RAW_DIR}/' --include='*/' --include='*.tif' --exclude='*'`

**Globus staging** (thymus datasets from PATH lab): Uses `stage_datasets_globus.py` with PATH lab GCP endpoint → HiPerGator. Transfer log: `/blue/maigan/smith6jt/globus_transfers.json`. Requires `module load globus` and session consent. See `stage_datasets_globus.py` docstring for endpoint IDs and path quirks.

**Sentinel files:**
| File | Created by | Meaning |
|------|-----------|---------|
| `data/raw/.staged` | `stage_datasets.sh` | Raw data rsync complete |
| `data/processed/stitched/cyc{NN}/.snakemake_complete` | Snakemake stitch rule | Stitching done for cycle |
| `data/processed/deconvolved/cyc{NN}/.snakemake_complete` | Snakemake decon rule | Deconvolution done for cycle |
| `data/processed/edf/cyc{NN}/.snakemake_complete` | Snakemake edf rule | EDF done for cycle |
| `data/processed/registered/.snakemake_complete` | Snakemake registration rule | All cycles registered |
| `qc_plots/.snakemake_complete_stitch` | Snakemake qc_stitch rule | Stitch QC report done |
| `qc_plots/.snakemake_complete_decon` | Snakemake qc_decon rule | Decon QC report done |
| `qc_plots/.snakemake_complete_edf` | Snakemake qc_edf rule | EDF QC report done |
| `data/processed/edf/.complete` | `cleanup_datasets.sh` | Dataset fully processed and cleaned |

**Cleanup QC Guard** (added 2026-02-16): `cleanup_datasets.sh` checks all 3 QC sentinels (`qc_plots/.snakemake_complete_{stitch,decon,edf}`) before proceeding. If any are missing, the dataset is skipped with "QC not complete". If all are present, the user is prompted to confirm QC review before deletion. Use `--force` to skip the interactive prompt (for re-runs after initial review). In `--dry-run` mode, shows what would be prompted without blocking.

```bash
./cleanup_datasets.sh              # Interactive: checks QC, prompts before each dataset
./cleanup_datasets.sh --force      # Skip prompts (after initial QC review)
./cleanup_datasets.sh --dry-run    # Preview without deleting
```

**Phases**: Discovery → Setup/Staging → Validation → SLURM batch → **QC Review** → Cleanup → Signal isolation → Segmentation
