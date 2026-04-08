# Command Line Interface

KINTSUGI provides a command-line interface for common operations.

## Commands

### kintsugi check

Check all dependencies and report their status.

```bash
kintsugi check
kintsugi check --verbose  # More detailed output
kintsugi check --json     # Output as JSON
```

### kintsugi info

Display system information and KINTSUGI version.

```bash
kintsugi info
```

### kintsugi template

Generate a configuration template file.

```bash
kintsugi template -o config.json
kintsugi template --output my_config.json
```

### kintsugi register

Run the registration workflow.

```bash
# Dry run (show what would be done)
kintsugi register config.json --dry-run

# Run registration
kintsugi register config.json

# Override config options
kintsugi register config.json --src /path/to/images --dst /path/to/output
```

### kintsugi scan

Preview what `kintsugi init` will find in a directory.

```bash
kintsugi scan /path/to/directory
```

### kintsugi init

Initialize a new KINTSUGI project with standard directory structure.

```bash
kintsugi init /path/to/project --name "My Project"
kintsugi init /path/to/project --name "My Project" --slurm  # Add SLURM support

# With microscope parameters
kintsugi init /path/to/project --name "My Project" \
    --tile-rows 9 --tile-cols 7 \
    --xy-pixel-size 377 --z-step-size 1500 \
    --numerical-aperture 0.75 --tissue-ri 1.44
```

**Options:**
- `--name`: Project name (required)
- `--slurm`: Add SLURM/Snakemake workflow support
- `--force`: Refresh templates on existing project
- `--tile-rows`, `--tile-cols`: Tile grid dimensions
- `--xy-pixel-size`, `--z-step-size`: Pixel sizes in nm
- `--numerical-aperture`, `--tissue-ri`: Optical parameters

### kintsugi install

Install optional dependency groups.

```bash
kintsugi install gpu        # GPU acceleration (CuPy for CUDA)
kintsugi install torch      # PyTorch for deep learning models
kintsugi install bio        # Spatial biology analysis (scanpy, scimap, squidpy)
kintsugi install viz        # Napari visualization
kintsugi install claude     # Claude Code MCP integration
kintsugi install dev        # Development tools (pytest, ruff, black, mypy)
kintsugi install all        # All optional features
```

## MCP Server Commands

KINTSUGI includes an MCP (Model Context Protocol) server for Claude Code integration.

### kintsugi mcp start

Start the MCP server for Claude Code integration.

```bash
kintsugi mcp start
```

The server exposes image processing tools that Claude Code can use for signal isolation and quality assessment.

### kintsugi mcp tools

List all available MCP tools.

```bash
kintsugi mcp tools
```

**Available Tools:**

| Category | Tools |
|----------|-------|
| **Signal Isolation** | `load_channel`, `subtract_blank`, `denoise`, `denoise_advanced`, `apply_clahe`, `clean_background`, `gaussian_subtract` |
| **Quality Assessment** | `assess_quality`, `compute_snr` |
| **Visualization** | `get_image_stats`, `get_thumbnail` |
| **Workflow** | `list_channels`, `save_processed`, `suggest_parameters`, `generate_jupyter_cell` |
| **Parameter Learning** | `get_learned_parameters`, `record_successful_parameters`, `suggest_with_learning`, `approve_and_learn`, `get_learning_statistics` |

### kintsugi mcp config

Generate and create Claude Code MCP configuration for a project.

```bash
kintsugi mcp config /path/to/project
```

This automatically creates `.claude/settings.local.json` with the MCP server configuration. If the file already exists, it adds the KINTSUGI server to your existing configuration.

Use `--print-only` to display the configuration without creating files:

```bash
kintsugi mcp config /path/to/project --print-only
```

## Configuration File Format

The configuration file is a JSON file with the following structure:

```json
{
    "src_dir": "/path/to/source/images",
    "dst_dir": "/path/to/output",
    "reference_image": "cycle1.tif",
    "image_type": "tif",
    "series": 0,
    "max_image_dim_px": 2048,
    "max_processed_image_dim_px": 2048,
    "micro_rigid_registrar_cls": "RigidRegistrar",
    "align_to_reference": true,
    "create_masks": true,
    "resolution_xyu": [0.325, 0.325, "um"],
    "channel_names": [],
    "compose_non_rigid": true,
    "crop_to_overlap": true
}
```

## Workflow Commands (HPC/Snakemake)

Manage the Snakemake-based processing pipeline for SLURM clusters.

### kintsugi workflow config

Generate Snakemake workflow configuration. Auto-detects SLURM accounts, calculates GPU and CPU slots, reads microscope parameters from `meta/experiment.json`.

```bash
kintsugi workflow config /path/to/project
kintsugi workflow config /path/to/project --print-only  # Preview without writing
```

### kintsugi workflow check

Show live per-account resource availability (GPU and CPU slots).

```bash
kintsugi workflow check /path/to/project
```

### kintsugi workflow run

Run the Snakemake processing pipeline. Auto-calculates concurrent job count from live resource availability.

```bash
kintsugi workflow run /path/to/project              # Full pipeline
kintsugi workflow run /path/to/project --dry-run     # Preview
kintsugi workflow run /path/to/project --cycles 1-3  # Specific cycles
kintsugi workflow run /path/to/project --forcerun stitch  # Force re-run
kintsugi workflow run /path/to/project --local --cores 4  # Local execution
kintsugi workflow run /path/to/project -j 16         # Override job count
kintsugi workflow run /path/to/project --dashboard   # Run with live progress dashboard
kintsugi workflow run /path/to/project --dashboard --dashboard-interval 15
```

When `--dashboard` is supplied, Snakemake runs in the background and the dashboard refreshes in the foreground. Pressing Ctrl+C **detaches** without killing the SLURM jobs.

### kintsugi workflow status

Display a live pipeline progress dashboard for one or more projects. Scans sentinel files, queries SLURM via `squeue`/`sacct`, parses log timings, and reports per-cycle, per-stage completion alongside live GPU/memory utilization and ETA estimates.

```bash
kintsugi workflow status /path/to/project              # One-shot snapshot
kintsugi workflow status /path/to/project --watch      # Auto-refresh
kintsugi workflow status /path/to/project -w -i 15     # Custom refresh interval (seconds)
kintsugi workflow status /path/to/project --json       # Machine-readable output
kintsugi workflow status /path/to/parent_dir --all-projects --watch  # Scan every project
kintsugi workflow status /path/to/project --no-hardware --no-estimates --no-jobs
```

**Flags:**

| Flag | Description |
|------|-------------|
| `--watch`, `-w` | Auto-refresh the dashboard until Ctrl+C |
| `--interval`, `-i` | Refresh interval in seconds (default: 30) |
| `--all-projects` | Recursively discover every directory under `PROJECT_DIR` containing `workflow/config.yaml` |
| `--json` | Emit a JSON snapshot instead of the Rich table (suitable for scripting) |
| `--no-hardware` | Hide the per-account GPU/CPU/memory utilization section |
| `--no-estimates` | Hide the wall-clock ETA estimate |
| `--no-jobs` | Hide the active SLURM job details column |

**What it reports per project:**

- **Cycle table** — stitch / decon / EDF status per cycle, per-stage timing, active job, and node assignment
- **Aggregate stages** — registration and signal-isolation status with timings
- **QC sentinels** — completion of `qc_stitch`, `qc_decon`, `qc_edf`, `qc_registration`, `qc_signal_isolation`
- **SLURM jobs** — job id, state, elapsed time, MaxRSS, GPU id, partition, account
- **Hardware utilization** — per-account GPU allocated/used/available and memory pools
- **Completion estimate** — wall-clock ETA derived from historical log timings and current parallelism

The dashboard is also available embedded in `kintsugi workflow run --dashboard` (see above). Both entry points share `src/kintsugi/dashboard.py` (`scan_project_progress`, `render_dashboard`, `watch_dashboard`).

## SLURM Commands

Manage legacy SLURM job submission.

### kintsugi slurm init

Add SLURM support to an existing project.

```bash
kintsugi slurm init /path/to/project
```

### kintsugi slurm submit

Submit processing jobs to SLURM.

```bash
kintsugi slurm submit /path/to/project
kintsugi slurm submit /path/to/project --steps decon,edf
kintsugi slurm submit /path/to/project --cycles 1-5
kintsugi slurm submit /path/to/project --dry-run
```

### kintsugi slurm status

Check SLURM job status for a project.

```bash
kintsugi slurm status /path/to/project
```

### kintsugi slurm cancel

Cancel all KINTSUGI SLURM jobs for a project.

```bash
kintsugi slurm cancel /path/to/project
```

## Environment Variables

| Variable | Description |
|----------|-------------|
| `KINTSUGI_DATA_DIR` | Default data directory |
| `KINTSUGI_DEVICE_MODE` | Processing backend: `gpu` or `cpu` (set by Snakemake) |
| `VIPS_PATH` | libvips binary directory (Windows) |
| `SLURM_ACCOUNT` | Default SLURM account (auto-detected) |
| `SLURM_JOB_ACCOUNT` | Current job's SLURM account (set by SLURM) |
