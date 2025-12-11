"""
Command-line interface for KINTSUGI.

Provides CLI entry points for:
- kintsugi: Main entry point
- kintsugi-register: Registration workflow
- kintsugi-check: Dependency checking
"""

import json
import sys
from pathlib import Path

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()


@click.group()
@click.version_option()
def main():
    """
    KINTSUGI: Multi-cycle immunofluorescence image processing.

    Knowledge Integration with New Technologies for Simplified
    User-Guided Image processing.
    """
    pass


@main.command()
@click.option("--verbose/--quiet", "-v/-q", default=True, help="Verbose output")
@click.option("--json", "output_json", is_flag=True, help="Output as JSON")
def check(verbose: bool, output_json: bool):
    """Check all KINTSUGI dependencies."""
    from kintsugi.deps import DependencyChecker

    checker = DependencyChecker()

    if output_json:
        result = checker.check_all(verbose=False)
        click.echo(json.dumps(result, indent=2))
    else:
        checker.check_all(verbose=verbose)


@main.command()
@click.argument("config", type=click.Path(exists=True))
@click.option("--src", "-s", type=click.Path(exists=True), help="Source directory")
@click.option("--dst", "-d", type=click.Path(), help="Destination directory")
@click.option("--reference", "-r", help="Reference image filename")
@click.option("--dry-run", is_flag=True, help="Show what would be done")
def register(config: str, src: str | None, dst: str | None, reference: str | None, dry_run: bool):
    """
    Run image registration workflow.

    CONFIG is a JSON configuration file with registration parameters.
    Command-line options override config file values.
    """
    config_path = Path(config)

    with open(config_path) as f:
        cfg = json.load(f)

    # Override with CLI options
    if src:
        cfg["src_dir"] = src
    if dst:
        cfg["dst_dir"] = dst
    if reference:
        cfg["reference_image"] = reference

    # Validate required fields
    required_fields = ["src_dir", "dst_dir"]
    missing = [f for f in required_fields if f not in cfg or not cfg[f]]

    if missing:
        console.print(f"[red]Missing required fields: {', '.join(missing)}[/red]")
        raise SystemExit(1)

    src_path = Path(cfg["src_dir"])
    dst_path = Path(cfg["dst_dir"])

    if not src_path.exists():
        console.print(f"[red]Source directory does not exist: {src_path}[/red]")
        raise SystemExit(1)

    if dry_run:
        console.print(
            Panel.fit(
                f"[bold]Registration Workflow (Dry Run)[/bold]\n\n"
                f"Source: {src_path}\n"
                f"Destination: {dst_path}\n"
                f"Reference: {cfg.get('reference_image', 'auto')}\n"
                f"Image type: {cfg.get('image_type', 'tif')}\n"
                f"Max dimension: {cfg.get('max_image_dim_px', 2048)}px",
                title="KINTSUGI",
            )
        )
        return

    # Create destination directory
    dst_path.mkdir(parents=True, exist_ok=True)

    console.print("[bold]Starting registration workflow...[/bold]")

    try:
        # Import Kreg and run registration
        # Note: Uses the notebooks/Kreg module via path manipulation
        # This will be improved when notebooks are properly packaged
        notebooks_path = Path(__file__).parent.parent.parent / "notebooks"
        if str(notebooks_path) not in sys.path:
            sys.path.insert(0, str(notebooks_path))

        from Kreg.registration import Valis

        registrar = Valis(
            src_dir=str(src_path),
            dst_dir=str(dst_path),
            reference_img_f=cfg.get("reference_image"),
            img_type=cfg.get("image_type", "tif"),
            series=cfg.get("series", 0),
            max_image_dim_px=cfg.get("max_image_dim_px", 2048),
            max_processed_image_dim_px=cfg.get("max_processed_image_dim_px", 2048),
            align_to_reference=cfg.get("align_to_reference", True),
            create_masks=cfg.get("create_masks", True),
            resolution_xyu=cfg.get("resolution_xyu"),
            channel_names=cfg.get("channel_names"),
            compose_non_rigid=cfg.get("compose_non_rigid", True),
            crop_to_overlap=cfg.get("crop_to_overlap", True),
        )

        # Run registration
        rigid_registrar, non_rigid_registrar, error_df = registrar.register()

        console.print("[green]Registration complete![/green]")

        if error_df is not None and len(error_df) > 0:
            console.print("\n[yellow]Registration errors:[/yellow]")
            console.print(error_df.to_string())

    except ImportError as e:
        console.print(f"[red]Import error: {e}[/red]")
        console.print("Make sure KINTSUGI is properly installed.")
        raise SystemExit(1)
    except Exception as e:
        console.print(f"[red]Registration failed: {e}[/red]")
        raise SystemExit(1)


@main.command()
@click.option("--output", "-o", type=click.Path(), help="Output file path")
def template(output: str | None):
    """Generate a template configuration file."""
    from kintsugi import get_config_template

    cfg = get_config_template()

    if output:
        with open(output, "w") as f:
            json.dump(cfg, f, indent=2)
        console.print(f"[green]Template saved to: {output}[/green]")
    else:
        console.print(json.dumps(cfg, indent=2))


@main.command()
def info():
    """Show KINTSUGI version and environment information."""
    import platform

    from kintsugi import __version__

    table = Table(title="KINTSUGI Information")
    table.add_column("Property", style="cyan")
    table.add_column("Value", style="green")

    table.add_row("Version", __version__)
    table.add_row("Python", platform.python_version())
    table.add_row("Platform", platform.platform())
    table.add_row("Architecture", platform.machine())

    # Check for GPU
    try:
        import torch

        if torch.cuda.is_available():
            gpu_info = f"{torch.cuda.device_count()} GPU(s): {torch.cuda.get_device_name(0)}"
        else:
            gpu_info = "Not available (CPU mode)"
    except ImportError:
        gpu_info = "PyTorch not installed"

    table.add_row("GPU", gpu_info)

    console.print(table)


# ============================================================================
# MCP Server Commands for Claude Code Integration
# ============================================================================


@main.group()
def mcp():
    """Claude Code MCP server commands."""
    pass


@mcp.command()
def start():
    """
    Start the KINTSUGI MCP server for Claude Code integration.

    The MCP server exposes image processing tools that Claude can use
    for signal isolation and quality assessment.

    Configure in Claude Code settings:

    \b
    {
        "mcpServers": {
            "kintsugi": {
                "command": "kintsugi",
                "args": ["mcp", "start"]
            }
        }
    }
    """
    try:
        from kintsugi.mcp import run_server
        import asyncio

        console.print("[bold green]Starting KINTSUGI MCP server...[/bold green]")
        asyncio.run(run_server())

    except ImportError as e:
        console.print(f"[red]MCP dependencies not installed: {e}[/red]")
        console.print("Install with: [bold]pip install kintsugi[claude][/bold]")
        raise SystemExit(1)


@mcp.command()
def tools():
    """List available MCP tools."""
    table = Table(title="KINTSUGI MCP Tools")
    table.add_column("Tool", style="cyan")
    table.add_column("Description", style="white")

    tool_list = [
        # Signal Isolation
        ("load_channel", "Load a channel image from a KINTSUGI project"),
        ("subtract_blank", "Subtract autofluorescence/blank channel from signal"),
        ("denoise", "Apply denoising filters (percentile, uniform, median)"),
        ("denoise_advanced", "Advanced denoising (N2V, NLM, BM3D, bilateral, adaptive)"),
        ("apply_clahe", "Apply Contrast Limited Adaptive Histogram Equalization"),
        ("clean_background", "Remove background and small objects"),
        ("gaussian_subtract", "Subtract Gaussian-blurred version for background removal"),
        # Quality Assessment
        ("assess_quality", "Assess channel quality using DL/heuristic metrics"),
        ("compute_snr", "Compute Signal-to-Noise Ratio"),
        # Visualization
        ("get_image_stats", "Get statistics for a loaded image"),
        ("get_thumbnail", "Get a downsampled thumbnail for preview"),
        # Workflow
        ("list_channels", "List available channels in the project"),
        ("save_processed", "Save processed channel to output directory"),
        ("suggest_parameters", "Analyze channel and suggest processing parameters"),
        ("generate_jupyter_cell", "Generate Jupyter cell for interactive tuning"),
        # Parameter Learning
        ("get_learned_parameters", "Get recommended params from tissue/marker history"),
        ("record_successful_parameters", "Record approved params for future learning"),
        ("suggest_with_learning", "Get suggestions combining analysis + learned history"),
        ("approve_and_learn", "Approve results and record all params for learning"),
        ("get_learning_statistics", "Get statistics about the learning database"),
    ]

    for name, desc in tool_list:
        table.add_row(name, desc)

    console.print(table)
    console.print("\n[dim]Use 'kintsugi mcp start' to start the server[/dim]")
    console.print("[dim]Learning tools remember successful parameters by tissue/marker[/dim]")


@mcp.command()
@click.argument("project_path", type=click.Path(exists=True))
def config(project_path: str):
    """
    Generate Claude Code MCP configuration for a project.

    PROJECT_PATH is the path to your KINTSUGI project directory.
    """
    import json as json_mod
    from pathlib import Path

    project_path = Path(project_path).resolve()

    config_json = {
        "mcpServers": {
            "kintsugi": {
                "command": "kintsugi",
                "args": ["mcp", "start"],
                "cwd": str(project_path),
            }
        }
    }

    console.print(Panel.fit(
        json_mod.dumps(config_json, indent=2),
        title="Claude Code MCP Configuration",
        subtitle="Add to .claude/settings.local.json",
    ))

    console.print("\n[bold]To configure Claude Code:[/bold]")
    console.print("1. Copy the JSON above")
    console.print("2. Add to your project's .claude/settings.local.json")
    console.print("3. Restart Claude Code")


# ============================================================================
# Project Commands
# ============================================================================


@main.command()
@click.argument("project_path", type=click.Path())
@click.option("--name", "-n", help="Project name")
@click.option("--description", "-d", default="", help="Project description")
def init(project_path: str, name: str | None, description: str):
    """
    Initialize a new KINTSUGI project.

    Creates the project directory structure and configuration.
    """
    from kintsugi.project import KintsugiProject

    try:
        project = KintsugiProject.create(
            project_path,
            name=name,
            description=description,
        )
        console.print(f"\n[green]Project created successfully![/green]")

    except Exception as e:
        console.print(f"[red]Failed to create project: {e}[/red]")
        raise SystemExit(1)


# Standalone entry points
def check_dependencies():
    """Entry point for kintsugi-check command."""
    from kintsugi.deps import check_dependencies as _check

    _check(verbose=True)


def register_cli():
    """Entry point for kintsugi-register command."""
    # This wraps the main CLI for direct invocation
    main(["register"])


if __name__ == "__main__":
    main()
