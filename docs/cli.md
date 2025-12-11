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

### kintsugi init

Initialize a new KINTSUGI project.

```bash
kintsugi init /path/to/project
kintsugi init /path/to/project --name "My Project" --description "Project description"
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

Generate Claude Code MCP configuration for a project.

```bash
kintsugi mcp config /path/to/project
```

This outputs a JSON configuration that you can add to `.claude/settings.local.json`:

```json
{
    "mcpServers": {
        "kintsugi": {
            "command": "kintsugi",
            "args": ["mcp", "start"],
            "cwd": "/path/to/project"
        }
    }
}
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

## Environment Variables

| Variable | Description |
|----------|-------------|
| `KINTSUGI_DATA_DIR` | Default data directory |
| `VIPS_PATH` | libvips binary directory (Windows) |

## Installing Optional Features

Install optional dependency groups:

```bash
# GPU acceleration (PyTorch + CuPy)
pip install kintsugi[gpu]

# Napari visualization
pip install kintsugi[viz]

# Spatial analysis (scanpy, scimap)
pip install kintsugi[analysis]

# Deep learning segmentation
pip install kintsugi[dl]

# Claude Code MCP integration
pip install kintsugi[claude]

# Advanced denoising (N2V, CARE)
pip install kintsugi[denoise]

# All optional features
pip install kintsugi[full]
```
