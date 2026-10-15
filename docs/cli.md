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
kintsugi install --list     # Show every optional group (14 in deps.py OPTIONAL_GROUPS)
kintsugi install gpu        # GPU acceleration (PyTorch + CuPy for CUDA)
kintsugi install analysis   # Spatial analysis (scanpy, scimap)
kintsugi install bio        # Bio formats I/O (OME-TIFF, LIF, etc.)
kintsugi install viz        # Napari visualization
kintsugi install mcp        # MCP server for AI assistants (mcp>=2.0,<3)
kintsugi install dev        # Development tools (pytest, ruff, black, mypy)
kintsugi install all        # All optional features
```

The groups are defined in `src/kintsugi/deps.py` (`OPTIONAL_GROUPS`, 14 groups:
`gpu`, `viz`, `dl`, `analysis`, `bio`, `claude`, `dev`, `docs`, `kronos`, `denoise`,
`optimize`, `rapids`, `workflow`, `full`). There is no `torch` group — PyTorch is
installed by `gpu` (and by `dl`, `denoise`, `kronos`). For development installs,
`pip install -e ".[mcp]"` is equivalent to `kintsugi install mcp`; do **not** run
`pip install kintsugi[mcp]`, which installs an unrelated PyPI package.

## MCP Server Commands

KINTSUGI includes an MCP (Model Context Protocol) server for Claude Code integration.
It requires the 2.x Python SDK (`mcp>=2.0,<3`): `pip install -e ".[mcp]"` or
`kintsugi install mcp`. The full guide — tools, resources, the prompt, the learning
loop, security model and troubleshooting — is [Claude Code Integration (MCP Server)](mcp.md).

### kintsugi mcp start

Start the MCP server for Claude Code integration.

```bash
kintsugi mcp start                                   # project from KINTSUGI_PROJECT_DIR or a project cwd
kintsugi mcp start --project-dir /path/to/project    # what generated client configs pass
```

The server speaks JSON-RPC over stdio, so this command writes nothing to stdout itself;
status and errors go to stderr. MCP clients start it from their configuration (see
`kintsugi mcp config`). `--project-dir` (which must contain `kintsugi_project.json`) sets the
project that `project_path` arguments default to; `KINTSUGI_PROJECT_DIR` does the same, and
without either the working directory is used only if it is a KINTSUGI project.

### kintsugi mcp tools

List the registered MCP tools, or sync the generated tool reference in the docs.

```bash
kintsugi mcp tools                                   # Rich table (name, title, read-only, description)
kintsugi mcp tools --markdown                        # Same as a Markdown table
kintsugi mcp tools --write docs/cli.md CLAUDE.md     # Regenerate the <!-- mcp-tools:* --> blocks
kintsugi mcp tools --check docs/cli.md CLAUDE.md     # Exit non-zero if a block is stale (CI runs this)
```

The registry (`kintsugi.mcp.tool_specs.TOOL_SPECS`) is the single source of truth. Docs
embed the generated table between `<!-- mcp-tools:table:begin -->` /
`<!-- mcp-tools:table:end -->` markers (or a names-only list between
`<!-- mcp-tools:list:begin -->` / `<!-- mcp-tools:list:end -->`); `--write` regenerates
them and `--check` fails on drift. To add a block to a page, write the markers with a
placeholder line between them and run `--write` — never hand-edit the table.

**Available tools** (generated):

<!-- mcp-tools:table:begin -->
| Tool | Title | Category | Read-only | Destructive | Description |
|------|-------|----------|:---------:|:------------:|-------------|
| `load_channel` | Load channel | Signal isolation | no | no | Load a channel image from a KINTSUGI project. Returns image metadata and prepares it for processing. |
| `subtract_blank` | Subtract blank (autofluorescence) | Signal isolation | no | no | Subtract autofluorescence/blank channel from signal. Uses optimized Dask operations for large images. |
| `analyze_weighted_subtraction` | Preview weighted subtraction | Signal isolation | yes | no | Preview the per-intensity-range weights for weighted autofluorescence subtraction WITHOUT applying it. Use this to verify the plan before calling subtract_blank with method='weighted'. |
| `denoise` | Denoise (basic filters) | Signal isolation | no | no | Apply denoising filters (percentile, uniform, median) to reduce noise. |
| `apply_clahe` | Apply CLAHE contrast | Signal isolation | no | no | Apply Contrast Limited Adaptive Histogram Equalization to enhance local contrast. |
| `clean_background` | Clean background | Signal isolation | no | no | Remove background and small objects from the image. |
| `gaussian_subtract` | Gaussian background subtraction | Signal isolation | no | no | Subtract Gaussian-blurred version to remove structured background/autofluorescence. |
| `denoise_advanced` | Denoise (advanced) | Signal isolation | no | no | Apply advanced denoising: adaptive (auto-select best filter), bilateral (edge- preserving), nlm (non-local means), bm3d (block matching 3D), n2v (self-supervised Noise2Void), or patch_svd. |
| `assess_quality` | Assess quality | Quality assessment | yes | no | Assess channel quality with heuristic metrics (SNR, contrast, dynamic range, background) and flag channels that need manual review. Returns a 0-1 quality score, the metrics and recommendations. |
| `compute_snr` | Compute SNR | Quality assessment | yes | no | Compute Signal-to-Noise Ratio for a channel. |
| `batch_assess_quality` | Batch quality assessment | Quality assessment | yes | no | Assess several loaded channels at once: scores sorted worst-first and a review queue of the channels below the threshold. |
| `get_image_stats` | Image statistics | Visualization | yes | no | Get statistics for a loaded image (min, max, mean, std, histogram). |
| `get_thumbnail` | Thumbnail preview | Visualization | yes | no | Get a downsampled thumbnail of the image for preview. |
| `compare_channels` | Compare two channels | Visualization | yes | no | Compare two loaded channels statistically (means, difference, correlation, SSIM) with a one-line interpretation; useful for before/after checks. |
| `get_intensity_profile` | Intensity profile | Visualization | yes | no | Intensity values along a horizontal or vertical line through a loaded channel, for spotting gradients or uneven illumination. |
| `list_channels` | List project channels | Workflow | yes | no | List all available channels in the current project/cycle. |
| `save_processed` | Save processed channel | Workflow | no | yes | Save the processed channel to the project output directory. |
| `get_processing_history` | Processing history | Workflow | yes | no | Get the processing history for a channel, showing all operations applied. |
| `suggest_parameters` | Suggest parameters (heuristic) | Workflow | yes | no | Rule-based suggestions for every operation from the channel's intensity statistics, without the learning database. Prefer suggest_with_learning when the tissue type and marker are known. |
| `generate_jupyter_cell` | Generate Jupyter tuning cell | Workflow | yes | no | Generate a Jupyter notebook cell for interactive parameter tuning with Kview2 widgets. |
| `cluster_channels` | Cluster channels | Signal isolation | no | no | Cluster channels by image feature similarity. Groups channels that share similar intensity profiles, SNR, and texture so parameters can be tuned once per cluster. |
| `propagate_parameters` | Propagate cluster parameters | Signal isolation | no | yes | Apply tuned signal isolation parameters from a representative channel to all members of a cluster. |
| `get_session_state` | Session state | Workflow | yes | no | Show what this server session holds: the project root it serves, every loaded channel with shape/dtype/metadata, per-channel processing step counts, and memory usage. |
| `unload_channel` | Unload channel | Workflow | no | no | Free server memory by dropping a loaded channel (or everything with clear_all). Files on disk are untouched. |
| `get_learned_parameters` | Learned parameters | Parameter learning | yes | no | Get recommended parameters based on learned history from similar tissue/marker combinations. Returns parameters with confidence scores. |
| `record_successful_parameters` | Record parameters | Parameter learning | no | no | Record a successful parameter set for future learning. Call after user approves processing results. |
| `suggest_with_learning` | Suggest parameters (learned + heuristic) | Parameter learning | yes | no | Get comprehensive parameter suggestions combining image analysis AND learned history. This is the primary recommendation tool. |
| `approve_and_learn` | Approve and learn | Parameter learning | no | no | Approve processing results and record ALL parameters for learning. Call when user approves final processed result. |
| `get_learning_statistics` | Learning database statistics | Parameter learning | yes | no | Get statistics about the parameter learning database. |
| `optimize_parameters` | Optimize parameters (Optuna) | Signal isolation | no | no | Find optimal signal isolation parameters via Bayesian optimization (Optuna). Automatically searches the parameter space — no manual slider adjustment needed. |
| `predict_parameters` | Predict parameters (Random Forest) | Signal isolation | yes | no | Predict signal isolation parameters from image features using a trained Random Forest model. Includes confidence scores and uncertainty estimates. |
| `estimate_background` | Estimate background (SMO) | Signal isolation | no | no | Estimate and subtract background using SMO (Silver Mountain Operator). Parameter- free — no blank channel needed. |

32 tools. Generated by `kintsugi mcp tools --write`; do not edit by hand.
<!-- mcp-tools:table:end -->

### kintsugi mcp config

Write MCP client configuration (and `AGENTS.md`) for a project.

```bash
kintsugi mcp config /path/to/project                   # every project-scoped client + snippets
kintsugi mcp config /path/to/project --client cursor   # one client (repeatable)
kintsugi mcp config /path/to/project --print-only      # show everything, write nothing
kintsugi mcp config /path/to/project --force           # replace existing kintsugi entries
```

By default it writes or merges the configuration of every client that keeps it in the
project (Claude Code, VS Code, Cursor, Gemini CLI, Antigravity, Codex, Zed), writes
`AGENTS.md` if the project has none, and prints the entry to paste for Claude Desktop,
Windsurf, Cline and JetBrains; nothing is written outside the project. Other servers and
settings in existing files are kept, and the generated files are added to the project's
`.gitignore` (they hold machine-specific absolute paths). A directory without
`kintsugi_project.json` is refused unless `--force` is given. `kintsugi init` writes the
same files (`--mcp-client NAME`, repeatable, or `--mcp-client none` to choose). For
Claude Code, `.mcp.json` at the project root holds the server definition:

```json
{
  "mcpServers": {
    "kintsugi": {
      "type": "stdio",
      "command": "/absolute/path/to/kintsugi",
      "args": ["mcp", "start", "--project-dir", "/path/to/project"]
    }
  }
}
```

The project travels in `args` rather than a `cwd` entry, because several clients ignore
`cwd`. `.claude/settings.local.json` approves only that server:
`{"enabledMcpjsonServers": ["kintsugi"]}`. The equivalent Claude Code CLI command is:

```bash
claude mcp add --scope project --transport stdio kintsugi -- \
  kintsugi mcp start --project-dir /path/to/project
```

### kintsugi mcp clients

List the MCP clients `kintsugi mcp config` knows: name, scope (`project` = written into
the project, `user` / `ide` = printed as a snippet), configuration file, instructions file
and documentation link.

```bash
kintsugi mcp clients
```

### kintsugi mcp pretrain

Pre-populate the parameter-learning database from existing processed data.

```bash
kintsugi mcp pretrain /path/to/project --tissue-type "lymph node"
kintsugi mcp pretrain /path/to/project --dry-run --verbose   # Preview without recording
```

Scans `data/processed/edf/` for signal/blank pairs, computes weighted subtraction
parameters for each, and records them with `algorithm_version='weighted_v1'` and
`user_approved=False` (bootstrap records count at reduced weight in recommendations).
Options: `-t/--tissue-type`, `--dry-run`, `-v/--verbose`. The databases live at
`<project>/.kintsugi/parameter_learning.db` and `~/.kintsugi/parameter_learning.db`.
An existing `<project>/.claude_workflow/` from an earlier release is read as is and moved
to `.kintsugi/` on the first write. Both databases are created when parameters are first
recorded; looking up recommendations or statistics never creates them.

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
