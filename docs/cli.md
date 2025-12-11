# Command Line Interface

KINTSUGI provides a command-line interface for common operations.

## Commands

### kintsugi check

Check all dependencies and report their status.

```bash
kintsugi check
kintsugi check --verbose  # More detailed output
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
