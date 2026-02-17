"""
Command-line interface for KINTSUGI.

Provides CLI entry points for:
- kintsugi: Main entry point
- kintsugi-register: Registration workflow
- kintsugi-check: Dependency checking
- kintsugi install: Install optional dependencies
"""

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

console = Console()

# Optional dependency groups - use direct package names to avoid PyPI collision
# (the 'kintsugi' package exists on PyPI as an unrelated project)
OPTIONAL_GROUPS: dict[str, dict[str, Any]] = {
    "gpu": {
        "description": "GPU acceleration (CuPy for CUDA)",
        "install_cmd": "pip install cupy-cuda12x",
        "packages": ["cupy"],
    },
    "torch": {
        "description": "PyTorch for deep learning models",
        "install_cmd": "pip install torch torchvision",
        "packages": ["torch", "torchvision"],
    },
    "bio": {
        "description": "Spatial biology analysis (scanpy, scimap, squidpy)",
        "install_cmd": "pip install scanpy scimap squidpy anndata",
        "packages": ["scanpy", "scimap", "squidpy", "anndata"],
    },
    "viz": {
        "description": "Napari visualization",
        "install_cmd": "pip install napari[all]",
        "packages": ["napari"],
    },
    "claude": {
        "description": "Claude Code MCP integration",
        "install_cmd": "pip install mcp anthropic",
        "packages": ["mcp", "anthropic"],
    },
    "dev": {
        "description": "Development tools (pytest, ruff, black, mypy)",
        "install_cmd": "pip install pytest pytest-cov ruff black mypy",
        "packages": ["pytest", "ruff", "black", "mypy"],
    },
}


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
    except (ImportError, OSError):
        gpu_info = "PyTorch not installed"

    table.add_row("GPU", gpu_info)

    console.print(table)


# ============================================================================
# Install Command for Optional Dependencies
# ============================================================================


@main.command()
@click.argument("group", required=False)
@click.option("--list", "-l", "list_groups", is_flag=True, help="List available groups")
def install(group: str | None, list_groups: bool):
    """
    Install optional dependency groups.

    GROUP is the name of the dependency group to install (gpu, torch, bio, viz, claude, dev).
    Use 'all' to install all optional dependencies.

    \b
    Examples:
        kintsugi install --list     # Show available groups
        kintsugi install gpu        # Install GPU acceleration
        kintsugi install bio        # Install spatial biology tools
        kintsugi install all        # Install everything
    """
    if list_groups or group is None:
        table = Table(title="Optional Dependency Groups")
        table.add_column("Group", style="cyan")
        table.add_column("Description", style="white")
        table.add_column("Packages", style="dim")

        for name, info in OPTIONAL_GROUPS.items():
            table.add_row(name, info["description"], ", ".join(info["packages"]))

        console.print(table)
        console.print("\n[dim]Usage: kintsugi install <group>[/dim]")
        console.print("[dim]       kintsugi install all[/dim]")
        return

    if group == "all":
        groups_to_install = list(OPTIONAL_GROUPS.keys())
    elif group in OPTIONAL_GROUPS:
        groups_to_install = [group]
    else:
        console.print(f"[red]Unknown group: {group}[/red]")
        console.print(f"Available groups: {', '.join(OPTIONAL_GROUPS.keys())}, all")
        raise SystemExit(1)

    for grp in groups_to_install:
        info = OPTIONAL_GROUPS[grp]
        console.print(f"\n[bold]Installing {grp}:[/bold] {info['description']}")
        console.print(f"[dim]Running: {info['install_cmd']}[/dim]")

        try:
            subprocess.run(
                info["install_cmd"],
                shell=True,
                check=True,
                capture_output=False,
            )
            console.print(f"[green]✓ {grp} installed successfully[/green]")
        except subprocess.CalledProcessError as e:
            console.print(f"[red]✗ Failed to install {grp}: {e}[/red]")
            if group != "all":
                raise SystemExit(1)

    console.print("\n[bold green]Installation complete![/bold green]")


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
        import asyncio

        from kintsugi.mcp import run_server

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
@click.option("--print-only", is_flag=True, help="Only print config, don't create file")
def config(project_path: str, print_only: bool):
    """
    Generate Claude Code MCP configuration for a project.

    PROJECT_PATH is the path to your KINTSUGI project directory.

    By default, creates .claude/settings.local.json in the project directory.
    Use --print-only to just display the configuration without creating files.
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

    if print_only:
        console.print(
            Panel.fit(
                json_mod.dumps(config_json, indent=2),
                title="Claude Code MCP Configuration",
            )
        )
        return

    # Create .claude directory and settings file
    claude_dir = project_path / ".claude"
    claude_dir.mkdir(exist_ok=True)

    settings_file = claude_dir / "settings.local.json"

    if settings_file.exists():
        # Check if we need to update existing config
        with open(settings_file) as f:
            existing = json_mod.load(f)

        if "mcpServers" not in existing:
            existing["mcpServers"] = {}

        if "kintsugi" in existing.get("mcpServers", {}):
            console.print(f"[yellow]KINTSUGI MCP config already exists in {settings_file}[/yellow]")
            console.print("Use --print-only to see the configuration.")
            return

        # Add kintsugi to existing config
        existing["mcpServers"]["kintsugi"] = config_json["mcpServers"]["kintsugi"]
        with open(settings_file, "w") as f:
            json_mod.dump(existing, f, indent=2)
        console.print(f"[green]Added KINTSUGI MCP config to existing {settings_file}[/green]")
    else:
        # Create new config file
        with open(settings_file, "w") as f:
            json_mod.dump(config_json, f, indent=2)
        console.print(f"[green]Created {settings_file}[/green]")

    console.print(
        Panel.fit(
            json_mod.dumps(config_json, indent=2),
            title="Claude Code MCP Configuration",
        )
    )
    console.print("\n[bold]Next steps:[/bold]")
    console.print("1. Open this project folder in VS Code")
    console.print("2. Start Claude Code - the MCP server will be available automatically")


@mcp.command()
@click.argument("project_dir", type=click.Path(exists=True))
@click.option("--tissue-type", "-t", default="unknown", help="Tissue type to record")
@click.option("--dry-run", is_flag=True, help="Preview what would be done without recording")
@click.option("--verbose", "-v", is_flag=True, help="Show detailed progress")
def pretrain(project_dir: str, tissue_type: str, dry_run: bool, verbose: bool):
    """
    Pre-populate the learning database from existing processed data.

    Scans a KINTSUGI project for signal/blank pairs in data/processed/edf/
    and computes weighted subtraction parameters for each. Records them
    in the learning database with algorithm_version='weighted_v1' and
    user_approved=False (bootstrap, lower weight in recommendations).

    \b
    Example:
        kintsugi mcp pretrain /path/to/project --tissue-type "lymph node"
        kintsugi mcp pretrain . --dry-run --verbose
    """
    from kintsugi.signal.bootstrap import bootstrap_from_project

    console.print(f"[bold]Bootstrapping weighted subtraction parameters from:[/bold] {project_dir}")
    if dry_run:
        console.print("[yellow]DRY RUN — no parameters will be recorded[/yellow]")

    result = bootstrap_from_project(
        project_dir=project_dir,
        tissue_type=tissue_type,
        dry_run=dry_run,
        verbose=verbose,
    )

    if result["status"] == "no_data":
        console.print(f"[yellow]{result['message']}[/yellow]")
        return

    if dry_run:
        console.print(f"\n[bold]Found {result['pairs_found']} signal/blank pairs[/bold]")
        if result.get("pairs"):
            table = Table(title="Pairs Found")
            table.add_column("Cycle", style="cyan")
            table.add_column("Marker", style="white")
            for p in result["pairs"]:
                table.add_row(p["cycle"], p["marker"])
            console.print(table)
        return

    # Show results
    console.print("\n[bold]Bootstrap Results:[/bold]")
    console.print(f"  Pairs found: {result.get('pairs_found', 0)}")
    console.print(f"  Pairs processed: {result.get('pairs_processed', 0)}")
    console.print(f"  Pairs recorded: {result.get('pairs_recorded', 0)}")

    if result.get("results"):
        table = Table(title="Processing Results")
        table.add_column("Cycle", style="cyan")
        table.add_column("Marker", style="white")
        table.add_column("Quality", style="green")
        table.add_column("Recorded", style="bold")

        for r in result["results"]:
            quality_str = f"{r.get('quality', 0):.3f}" if r.get("recorded") else "—"
            recorded_str = "✓" if r.get("recorded") else f"✗ {r.get('error', '')}"
            table.add_row(
                r.get("cycle", ""),
                r.get("marker", ""),
                quality_str,
                recorded_str,
            )
        console.print(table)


# ============================================================================
# SLURM Commands for HPC Job Submission
# ============================================================================


@main.group()
def slurm():
    """SLURM job submission commands for HPC clusters."""
    pass


@slurm.command("init")
@click.argument("project_dir", type=click.Path(exists=True))
@click.option("--account", "-a", help="HPC account name (auto-detected if not provided)")
@click.option("--partition", "-p", default=None, help="SLURM partition (default: gpu)")
@click.option("--qos", "-q", default=None, help="SLURM Quality of Service")
@click.option("--gpu-type", "-g", default=None, help="GPU type for constraints (auto-detected)")
def slurm_init(
    project_dir: str,
    account: str | None,
    partition: str | None,
    qos: str | None,
    gpu_type: str | None,
):
    """
    Add SLURM support to an existing KINTSUGI project.

    Creates the slurm/ directory with configuration files and job script
    symlinks. Auto-detects HPC settings from the environment when possible.

    PROJECT_DIR is the path to your existing KINTSUGI project directory.

    \b
    Examples:
        kintsugi slurm init .                     # Current directory
        kintsugi slurm init /path/to/project      # Specific project
        kintsugi slurm init . --account myacct    # With explicit account
    """
    from pathlib import Path

    from kintsugi.project import KintsugiProject

    project_dir = Path(project_dir).resolve()
    config_file = project_dir / "kintsugi_project.json"

    if not config_file.exists():
        console.print(f"[red]Not a KINTSUGI project: {project_dir}[/red]")
        console.print("Run 'kintsugi init' first to create a project.")
        raise SystemExit(1)

    try:
        project = KintsugiProject.load(project_dir)
        slurm_dir = project.setup_slurm(
            account=account,
            partition=partition,
            qos=qos,
            gpu_type=gpu_type,
        )

        console.print(f"\n[green]SLURM support added to {project_dir}[/green]")
        console.print("\n[bold]Next steps:[/bold]")
        console.print(f"  1. Review and edit: {slurm_dir / 'config.sh'}")
        console.print(f"  2. Submit jobs: kintsugi slurm submit {project_dir}")
        console.print(
            f"  3. Or directly: {project._kintsugi_path}/slurm/submit.sh --project {project_dir}"
        )

    except Exception as e:
        console.print(f"[red]Failed to add SLURM support: {e}[/red]")
        raise SystemExit(1)


@slurm.command("submit")
@click.argument("project_dir", type=click.Path(exists=True))
@click.option(
    "--steps",
    "-s",
    default="all",
    help="Comma-separated steps: correction,stitch,decon,edf,all (default: all)",
)
@click.option("--cycles", "-c", default=None, help="Override cycles: '1-7' or '1,2,5'")
@click.option("--dry-run", "-n", is_flag=True, help="Show commands without submitting")
@click.option(
    "--use-burst",
    "-b",
    is_flag=True,
    help="Also submit to burst QOS for faster processing (preemptible)",
)
def slurm_submit(project_dir: str, steps: str, cycles: str | None, dry_run: bool, use_burst: bool):
    """
    Submit SLURM jobs for a KINTSUGI project.

    This is a convenience wrapper around the submit.sh script. It will use
    the configuration from PROJECT_DIR/slurm/config.sh.

    PROJECT_DIR is the path to your KINTSUGI project directory.

    \b
    Examples:
        kintsugi slurm submit .                   # Submit all steps
        kintsugi slurm submit . --steps decon     # Just deconvolution
        kintsugi slurm submit . --cycles 1-3      # Specific cycles
        kintsugi slurm submit . --dry-run         # Preview commands
        kintsugi slurm submit . --use-burst       # Also submit burst jobs
    """
    from pathlib import Path

    project_dir = Path(project_dir).resolve()

    # Validate project
    if not (project_dir / "kintsugi_project.json").exists():
        console.print(f"[red]Not a KINTSUGI project: {project_dir}[/red]")
        console.print("Run 'kintsugi init' first to create a project.")
        raise SystemExit(1)

    # Check for SLURM config
    config_file = project_dir / "slurm" / "config.sh"
    if not config_file.exists():
        console.print("[red]SLURM not configured for this project.[/red]")
        console.print("Run 'kintsugi slurm init' first to add SLURM support.")
        raise SystemExit(1)

    # Find submit.sh script
    try:
        from kintsugi.project import KintsugiProject

        project = KintsugiProject.load(project_dir)
        submit_script = project._kintsugi_path / "slurm" / "submit.sh"
    except Exception:
        # Fallback: try to find it relative to this file
        submit_script = Path(__file__).parent.parent.parent / "slurm" / "submit.sh"

    if not submit_script.exists():
        console.print(f"[red]Submit script not found: {submit_script}[/red]")
        raise SystemExit(1)

    # Build command
    cmd = [str(submit_script), "--project", str(project_dir), "--steps", steps]

    if cycles:
        cmd.extend(["--cycles", cycles])

    if dry_run:
        cmd.append("--dry-run")

    if use_burst:
        cmd.append("--use-burst")

    console.print(f"[bold]Running:[/bold] {' '.join(cmd)}")
    console.print()

    try:
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            raise SystemExit(result.returncode)
    except FileNotFoundError:
        console.print("[red]Could not execute submit script[/red]")
        raise SystemExit(1)


@slurm.command("status")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def slurm_status(project_dir: str):
    """
    Show status of SLURM jobs for a project.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from pathlib import Path

    project_dir = Path(project_dir).resolve()
    project_name = project_dir.name

    # Run squeue to get job status
    cmd = [
        "squeue",
        "-u",
        os.environ.get("USER", ""),
        "-n",
        f"kintsugi_*_{project_name}",
        "--format=%.18i %.9P %.30j %.8u %.2t %.10M %.6D %R",
    ]

    console.print(f"[bold]SLURM jobs for:[/bold] {project_name}")
    console.print()

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.stdout.strip():
            console.print(result.stdout)
        else:
            console.print("[dim]No running jobs found.[/dim]")

        # Also show recent runs
        runs_dir = project_dir / "slurm" / "runs"
        if runs_dir.exists():
            runs = sorted(runs_dir.iterdir(), reverse=True)[:5]
            if runs:
                console.print("\n[bold]Recent runs:[/bold]")
                for run in runs:
                    console.print(f"  {run.name}")

    except FileNotFoundError:
        console.print("[yellow]squeue not available (not on a SLURM cluster?)[/yellow]")
    except Exception as e:
        console.print(f"[red]Error checking status: {e}[/red]")


@slurm.command("cancel")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--run-id", default=None, help="Cancel jobs from a specific run ID only")
def slurm_cancel(project_dir: str, run_id: str | None):
    """
    Cancel all SLURM jobs for a project.

    If --run-id is provided, cancels only the jobs from that specific run
    (reads from slurm/runs/<run-id>/job_ids.txt). Otherwise, cancels all
    running kintsugi jobs matching this project name.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi slurm cancel .                          # Cancel all project jobs
        kintsugi slurm cancel . --run-id 20260206_143000 # Cancel specific run
    """
    from pathlib import Path

    project_dir = Path(project_dir).resolve()
    project_name = project_dir.name

    if run_id:
        # Cancel specific run's jobs from job_ids.txt
        job_ids_file = project_dir / "slurm" / "runs" / run_id / "job_ids.txt"
        if not job_ids_file.exists():
            console.print(f"[red]No job_ids.txt found for run {run_id}[/red]")
            console.print(f"  Expected: {job_ids_file}")
            raise SystemExit(1)

        job_ids = [
            line.strip() for line in job_ids_file.read_text().strip().splitlines() if line.strip()
        ]

        if not job_ids:
            console.print(f"[yellow]No job IDs recorded for run {run_id}[/yellow]")
            return

        console.print(f"[bold]Cancelling {len(job_ids)} job(s) from run {run_id}:[/bold]")
        for jid in job_ids:
            console.print(f"  {jid}")

        try:
            result = subprocess.run(
                ["scancel"] + job_ids,
                capture_output=True,
                text=True,
            )
            if result.returncode == 0:
                console.print(f"[green]Cancelled {len(job_ids)} job(s)[/green]")
            else:
                # scancel may warn about already-completed jobs, which is fine
                stderr = result.stderr.strip()
                if stderr:
                    console.print(f"[yellow]{stderr}[/yellow]")
                console.print("[green]Cancel request sent[/green]")
        except FileNotFoundError:
            console.print("[yellow]scancel not available (not on a SLURM cluster?)[/yellow]")
    else:
        # Cancel all kintsugi jobs for this project via squeue pattern
        try:
            result = subprocess.run(
                [
                    "squeue",
                    "-u",
                    os.environ.get("USER", ""),
                    "--name",
                    f"kintsugi_%_{project_name}",
                    "-h",
                    "-o",
                    "%i",
                ],
                capture_output=True,
                text=True,
            )
            job_ids = [line.strip() for line in result.stdout.strip().splitlines() if line.strip()]

            if not job_ids:
                console.print(f"[dim]No running kintsugi jobs found for {project_name}[/dim]")
                return

            console.print(f"[bold]Cancelling {len(job_ids)} job(s) for {project_name}:[/bold]")
            for jid in job_ids:
                console.print(f"  {jid}")

            cancel_result = subprocess.run(
                ["scancel"] + job_ids,
                capture_output=True,
                text=True,
            )
            if cancel_result.returncode == 0:
                console.print(f"[green]Cancelled {len(job_ids)} job(s)[/green]")
            else:
                stderr = cancel_result.stderr.strip()
                if stderr:
                    console.print(f"[yellow]{stderr}[/yellow]")
                console.print("[green]Cancel request sent[/green]")

        except FileNotFoundError:
            console.print("[yellow]squeue/scancel not available (not on a SLURM cluster?)[/yellow]")
        except Exception as e:
            console.print(f"[red]Error cancelling jobs: {e}[/red]")


# ============================================================================
# Workflow Commands (Snakemake)
# ============================================================================


@main.group()
def workflow():
    """Manage Snakemake-based processing workflows."""
    pass


@workflow.command("config")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--print-only", is_flag=True, help="Print config to stdout without writing")
def workflow_config(project_dir: str, print_only: bool):
    """
    Generate workflow/config.yaml for a project.

    Reads meta/experiment.json (auto-translating CODEX field names) and
    slurm/config.sh to create a unified Snakemake configuration file.

    All microscope and acquisition parameters are written to config.yaml so
    that workflow scripts never need to read experiment.json directly.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    import json
    import shutil
    from pathlib import Path

    import yaml

    from kintsugi.project import ExperimentConfig

    project_dir = Path(project_dir).resolve()

    # Load experiment config via ExperimentConfig.from_dict() for CODEX translation
    exp_path = project_dir / "meta" / "experiment.json"
    if not exp_path.exists():
        console.print(f"[red]Missing: {exp_path}[/red]")
        console.print("Cannot generate workflow config without experiment metadata.")
        raise SystemExit(1)

    with open(exp_path) as f:
        raw_experiment = json.load(f)

    exp = ExperimentConfig.from_dict(raw_experiment.copy())

    # Load channel names from CHANNELNAMES.txt
    channel_names_dict: dict[int, list[str]] = {}
    try:
        sys.path.insert(0, str(Path(__file__).parent.parent.parent.resolve() / "notebooks"))
        from Kio import load_channel_names

        channel_names_dict = load_channel_names(
            project_dir / "meta",
            channels_per_cycle=raw_experiment.get(
                "channels_per_cycle", raw_experiment.get("numChannels", 4)
            ),
            make_unique=True,
        )
        console.print(
            f"  Loaded channel names: {sum(len(v) for v in channel_names_dict.values())} markers"
        )
    except FileNotFoundError:
        console.print(
            "[yellow]  Warning: meta/CHANNELNAMES.txt not found — EDF will use CH# names[/yellow]"
        )
    except Exception as e:
        console.print(f"[yellow]  Warning: Could not load channel names: {e}[/yellow]")

    # Validate critical fields — fail loudly, no silent defaults
    if exp.n_cycles is None:
        console.print("[red]ERROR: n_cycles not found in experiment.json[/red]")
        console.print("Expected field: 'n_cycles' (KINTSUGI) or 'numCycles' (CODEX)")
        raise SystemExit(1)

    if exp.n_zplanes is None:
        console.print("[red]ERROR: n_zplanes not found in experiment.json[/red]")
        console.print("Expected field: 'n_zplanes' (KINTSUGI) or 'numZPlanes' (CODEX)")
        raise SystemExit(1)

    # Detect KINTSUGI installation path
    kintsugi_dir = Path(__file__).parent.parent.parent.resolve()

    # Build cycles and channels lists
    cycles_list = list(range(1, exp.n_cycles + 1))
    channels_list = list(range(1, exp.channels_per_cycle + 1))

    # Build wavelengths dict for YAML (string keys for readability)
    wavelengths_dict = {}
    for ch in sorted(exp.wavelengths.keys()):
        ex, em = exp.wavelengths[ch]
        wavelengths_dict[str(ch)] = [float(ex), float(em)]

    # Discover ALL user SLURM accounts (excluding burst -b accounts)
    from kintsugi.hpc import detect_multi_account_resources, get_user_slurm_accounts

    all_accounts = get_user_slurm_accounts()
    # Filter: exclude burst accounts
    accounts = [a for a in all_accounts if not a.endswith("-b")]
    if not accounts:
        console.print("[yellow]  Warning: No SLURM accounts detected[/yellow]")

    # GPU partition: comma-separated list for SLURM native fallback
    partition_gpu = "hpg-b200,hpg-turin"
    partition_cpu = "hpg-default"

    # Detect multi-account resources via sacctmgr
    pool = detect_multi_account_resources(
        accounts,
        resources_config={
            "cpus_per_cpu_job": 8,
            "cpu_utilization_cap": 0.85,
            "partition_gpu": partition_gpu,
            "partition_cpu": partition_cpu,
        },
    )

    # Build config dict and dump via PyYAML for correct types
    # Convert channel_names_dict keys to int for YAML (Snakemake reads as int)
    channel_names_yaml = (
        {int(k): v for k, v in channel_names_dict.items()} if channel_names_dict else {}
    )

    config_dict = {
        # Paths
        "project_dir": project_dir.as_posix(),
        "kintsugi_dir": kintsugi_dir.as_posix(),
        # Processing scope
        "cycles": cycles_list,
        "channels": channels_list,
        "output_format": "tiff",
        # Channel names (from CHANNELNAMES.txt)
        "channel_names": channel_names_yaml,
        # Tile grid
        "tile_rows": exp.tile_rows,
        "tile_cols": exp.tile_cols,
        "tile_overlap": exp.tile_overlap,
        # Microscope parameters
        "xy_pixel_size": exp.xy_pixel_size,
        "z_step_size": exp.z_step_size,
        "numerical_aperture": exp.numerical_aperture,
        "tissue_refractive_index": exp.tissue_refractive_index,
        "wavelengths": wavelengths_dict,
        "n_zplanes": exp.n_zplanes,
        # Stitching parameters
        "ncc_threshold": 0.078,
        "pou": 0.5,
        "blend_sigma": 15.0,
        "basic_flatfield_min": 0.1,
        "basic_max_iterations": 500,
        "basic_optimization_tolerance": 1.0e-6,
        # EDF parameters
        "edf": {
            "radius_x": 2,
            "radius_y": 2,
            "sigma": 10.0,
            "tiles": [3, 3],
            "blend_depth": 0,
            "z_smooth_sigma": 1.0,
        },
        # Registration parameters
        "registration": {
            "reference_cycle": 1,
            "reference_channel": 1,
            "max_image_dim": 2048,
            "rigid_only": True,
            "feature_detector": "VggFD",
            "imgs_ordered": True,
            "align_to_reference": True,
        },
        "quality_gate": {
            "max_zero_pct": 0.10,
            "max_sat_pct": 0.05,
            "min_valid_slices": 3,
        },
        # SLURM resources
        "resources": {
            "accounts": [
                {
                    "name": acct["name"],
                    "partition_gpu": acct["partition_gpu"],
                    "partition_cpu": acct["partition_cpu"],
                    "gpu_slots": acct["gpu_slots"],
                    "cpu_slots": acct["cpu_slots"],
                }
                for acct in pool["accounts"]
            ],
            "total_gpu_slots": pool["total_gpu_slots"],
            "total_cpu_slots": pool["total_cpu_slots"],
            "total_slots": pool["total_slots"],
            "cpu_utilization_cap": 0.85,
            "qos": "",
            "conda_env": "kintsugi",
            # GPU job resources (CuPy processes in GPU memory)
            "mem_stitch": 48000,
            "mem_decon": 48000,
            "mem_edf": 48000,
            "time_stitch": 240,
            "time_decon": 240,
            "time_edf": 60,
            "mem_registration": 64000,
            "time_registration": 120,
            "mem_qc_registration": 16000,
            "time_qc_registration": 30,
            # CPU job resources (float64 in system memory, needs more)
            "cpu_cpus_per_task": 8,
            "cpu_time_multiplier": 5,
            "cpu_mem_stitch": 48000,
            "cpu_mem_decon": 128000,
            "cpu_mem_edf": 96000,
            "cpu_mem_registration": 64000,
        },
    }

    # Generate YAML with header comment
    header = (
        "# =============================================================================\n"
        "# KINTSUGI Snakemake Pipeline Configuration\n"
        "# =============================================================================\n"
        "# Auto-generated from meta/experiment.json via ExperimentConfig.from_dict()\n"
        "# Edit paths and resource settings as needed.\n"
        "# =============================================================================\n\n"
    )
    config_content = header + yaml.dump(
        config_dict, default_flow_style=False, sort_keys=False, allow_unicode=True
    )

    if print_only:
        console.print(config_content)
        return

    # Write config
    wf_dir = project_dir / "workflow"
    wf_dir.mkdir(parents=True, exist_ok=True)
    config_file = wf_dir / "config.yaml"
    config_file.write_text(config_content)
    console.print(f"[green]Created workflow config:[/green] {config_file}")

    # Show key parameters for verification
    console.print(f"  n_cycles: {exp.n_cycles}")
    console.print(f"  n_zplanes: {exp.n_zplanes}")
    console.print(f"  tile_grid: {exp.tile_rows}x{exp.tile_cols}")
    console.print(f"  tile_overlap: {exp.tile_overlap}")
    console.print(f"  xy_pixel_size: {exp.xy_pixel_size}")
    console.print(f"  wavelengths: {len(exp.wavelengths)} channels")
    if channel_names_dict:
        console.print(
            f"  channel_names: {sum(len(v) for v in channel_names_dict.values())} markers"
        )
    if pool["total_gpu_slots"] > 0:
        acct_summary = ", ".join(f"{a['name']}({a['gpu_slots']}G)" for a in pool["accounts"])
        console.print(
            f"  multi-account: {acct_summary} = "
            f"{pool['total_gpu_slots']} GPU slots (GPU-only scheduling)"
        )

    # Copy Snakefile, profiles, and scripts from KINTSUGI installation.
    # Always overwrite Snakefile (pipeline logic) and profiles (precommand,
    # cache redirection, SLURM settings). Scripts are copied per-file: new
    # scripts are always added, but existing scripts are not overwritten
    # (users may customize wrapper scripts per-project).
    src_wf = kintsugi_dir / "workflow"
    if src_wf.exists():
        for item in ["Snakefile", "profiles"]:
            src = src_wf / item
            dst = wf_dir / item
            if src.exists():
                if src.is_dir():
                    if dst.exists():
                        shutil.rmtree(str(dst))
                    shutil.copytree(str(src), str(dst))
                else:
                    shutil.copy2(str(src), str(dst))
                console.print(f"  Updated {item}")

        # Scripts: per-file copy (add new scripts, don't overwrite existing)
        src_scripts = src_wf / "scripts"
        dst_scripts = wf_dir / "scripts"
        if src_scripts.exists():
            dst_scripts.mkdir(parents=True, exist_ok=True)
            added = []
            for src_file in src_scripts.iterdir():
                if src_file.is_file():
                    dst_file = dst_scripts / src_file.name
                    if not dst_file.exists():
                        shutil.copy2(str(src_file), str(dst_file))
                        added.append(src_file.name)
            if added:
                console.print(f"  Added scripts: {', '.join(added)}")
            else:
                console.print("  Scripts up to date")

    console.print()
    console.print("[bold]Next steps:[/bold]")
    console.print(f"  cd {wf_dir}")
    console.print("  kintsugi workflow check .                  # Verify resources")
    console.print("  snakemake -n                               # Dry run")
    console.print("  kintsugi workflow run .                    # Submit via SLURM")


@workflow.command("check")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def workflow_check(project_dir: str):
    """
    Check SLURM resource availability for multi-account execution.

    Reads workflow/config.yaml, queries all configured accounts via
    sacctmgr, and prints recommended snakemake command with correct -j flag.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from pathlib import Path

    import yaml

    from kintsugi.hpc import detect_live_multi_account

    project_dir = Path(project_dir).resolve()
    config_path = project_dir / "workflow" / "config.yaml"

    if not config_path.exists():
        console.print("[red]No workflow/config.yaml found.[/red]")
        console.print("Run 'kintsugi workflow config .' first.")
        raise SystemExit(1)

    with open(config_path) as f:
        wf_config = yaml.safe_load(f)

    res = wf_config.get("resources", {})
    account_list = [a["name"] for a in res.get("accounts", [])]

    if not account_list:
        console.print("[yellow]No accounts configured in workflow/config.yaml[/yellow]")
        console.print("Run 'kintsugi workflow config .' to regenerate.")
        raise SystemExit(1)

    util_cap = res.get("cpu_utilization_cap", 0.85)

    console.print(f"\n[bold]Resource Check[/bold]  ({project_dir.name})\n")

    pool = detect_live_multi_account(
        account_list,
        resources_config={
            "cpus_per_cpu_job": res.get("cpu_cpus_per_task", 8),
            "mem_per_cpu_job": res.get("cpu_mem_decon", 48000) // 1000,
            "cpu_utilization_cap": util_cap,
        },
    )

    # Per-account breakdown
    for acct in pool["accounts"]:
        alloc = acct["alloc"]
        usage = acct["usage"]
        console.print(f"  [bold]Account: {acct['name']}[/bold]")
        console.print(
            f"    Allocation: {alloc['cpus']} CPUs, {alloc['mem_gb']}GB, {alloc['gpus']} GPUs"
        )
        console.print(
            f"    In use:     {usage['cpus']} CPUs, {usage['mem_gb']}GB, {usage['gpus']} GPUs"
        )
        console.print(f"    GPU slots: {acct['gpu_slots']} max, {acct['gpu_avail']} available")
        console.print(
            f"    CPU slots: {acct['cpu_slots']} max ({int(util_cap * 100)}% cap), "
            f"{acct['cpu_avail']} available"
        )
        console.print()

    gpu_avail = pool["total_gpu_avail"]
    gpu_total = pool["total_gpu_slots"]

    console.print(f"  Total GPU: {gpu_total} max, {gpu_avail} available")
    console.print(f"  (CPU pool: {pool['total_cpu_slots']} slots — not used; GPU-only scheduling)")

    if gpu_avail > 0:
        console.print("\n  [bold]Recommended:[/bold]")
        console.print(f"    snakemake --profile profiles/slurm -j {gpu_avail}")
    else:
        console.print(
            "\n  [yellow]No GPU slots available right now. Other users are using the full allocation.[/yellow]"
        )
        console.print(f"  Max allocation would give {gpu_total} GPU slots when free.")


@workflow.command("run")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--jobs", "-j", default=8, type=int, help="Max concurrent SLURM jobs (default: 8)")
@click.option("--dry-run", "-n", is_flag=True, help="Preview without executing")
@click.option("--local", is_flag=True, help="Run locally instead of via SLURM")
@click.option("--cores", default=4, type=int, help="CPU cores for local execution (default: 4)")
@click.option(
    "--forcerun",
    default=None,
    help="Force re-run a specific rule (stitch/deconvolve/edf/registration)",
)
@click.option("--cycles", "-c", default=None, help="Override cycles: '1-3' or '1,2,5'")
def workflow_run(
    project_dir: str,
    jobs: int,
    dry_run: bool,
    local: bool,
    cores: int,
    forcerun: str | None,
    cycles: str | None,
):
    """
    Run the Snakemake processing pipeline.

    Requires workflow/config.yaml to exist. Create it with:
      kintsugi workflow config .

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow run .                    # Full pipeline via SLURM
        kintsugi workflow run . --dry-run          # Preview
        kintsugi workflow run . --local --cores 4  # Local execution
        kintsugi workflow run . --forcerun stitch  # Force re-stitch
    """
    from pathlib import Path

    project_dir = Path(project_dir).resolve()
    wf_dir = project_dir / "workflow"

    if not (wf_dir / "config.yaml").exists():
        console.print("[red]No workflow/config.yaml found.[/red]")
        console.print("Run 'kintsugi workflow config .' first.")
        raise SystemExit(1)

    if not (wf_dir / "Snakefile").exists():
        console.print("[red]No workflow/Snakefile found.[/red]")
        console.print("Run 'kintsugi workflow config .' to set up the workflow.")
        raise SystemExit(1)

    # Build snakemake command
    cmd = ["snakemake", "--directory", str(wf_dir), "--snakefile", str(wf_dir / "Snakefile")]

    if local:
        cmd.extend(["--cores", str(cores)])
    else:
        profile_dir = wf_dir / "profiles" / "slurm"
        if profile_dir.exists():
            # Query live resource availability (allocation minus current usage)
            import yaml

            from kintsugi.hpc import detect_live_multi_account

            with open(wf_dir / "config.yaml") as f:
                wf_config = yaml.safe_load(f)
            res = wf_config.get("resources", {})
            account_list = [a["name"] for a in res.get("accounts", [])]

            pool = detect_live_multi_account(
                account_list,
                resources_config={
                    "cpus_per_cpu_job": res.get("cpu_cpus_per_task", 8),
                    "mem_per_cpu_job": res.get("cpu_mem_decon", 48000) // 1000,
                    "cpu_utilization_cap": res.get("cpu_utilization_cap", 0.85),
                },
            )
            # GPU-only scheduling: -j is limited to available GPU slots.
            # Measured speedups (CX_19-003): stitch 25x, decon 5x, EDF 15x.
            # Full cycle ~22 min GPU vs ~282 min CPU — queuing always wins.
            total = pool["total_gpu_avail"] or pool["total_gpu_slots"] or jobs

            console.print(
                f"  Resources: {pool['total_gpu_avail']} GPU available "
                f"(of {pool['total_gpu_slots']} total GPU slots)"
            )

            cmd.extend(["--profile", str(profile_dir), "-j", str(total)])
        else:
            console.print("[yellow]No SLURM profile found, running locally.[/yellow]")
            cmd.extend(["--cores", str(cores)])

    if dry_run:
        cmd.append("-n")

    if forcerun:
        cmd.extend(["--forcerun", forcerun])

    if cycles:
        # Parse cycle range and build specific targets
        import re

        cycle_nums = []
        if "-" in cycles:
            m = re.match(r"(\d+)-(\d+)", cycles)
            if m:
                cycle_nums = list(range(int(m.group(1)), int(m.group(2)) + 1))
        else:
            cycle_nums = [int(c.strip()) for c in cycles.split(",")]

        if cycle_nums:
            # Read project_dir from config to build target paths
            import yaml

            with open(wf_dir / "config.yaml") as f:
                wf_config = yaml.safe_load(f)
            proj = wf_config["project_dir"]
            targets = [
                f"{proj}/data/processed/edf/cyc{c:02d}/.snakemake_complete" for c in cycle_nums
            ]
            cmd.extend(targets)

    console.print(f"[bold]Running:[/bold] {' '.join(cmd)}")
    console.print()

    try:
        result = subprocess.run(cmd, check=False)
        if result.returncode != 0:
            raise SystemExit(result.returncode)
    except FileNotFoundError:
        console.print("[red]snakemake not found. Install with:[/red]")
        console.print("  pip install snakemake snakemake-executor-plugin-slurm")
        raise SystemExit(1)


# ============================================================================
# Export Commands (TissUUmaps)
# ============================================================================


@workflow.group("export")
def export_group():
    """Export processed data for TissUUmaps web visualization.

    \b
    Commands:
        prepare  - Convert registered images to DZI tiles + generate .tmap
        deploy   - Rsync exported data to a remote server
        status   - Show export status for a project
    """
    pass


@export_group.command("prepare")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--force", "-f", is_flag=True, help="Force reconversion of all files")
@click.option("--dry-run", "-n", is_flag=True, help="Preview what would be exported")
@click.option("--include-blanks", is_flag=True, help="Include blank/autofluorescence channels")
def export_prepare(project_dir: str, force: bool, dry_run: bool, include_blanks: bool):
    """
    Convert processed images to DZI tiles and generate a .tmap project file.

    Reads registered images from data/processed/registered/, converts each to
    a DZI tile pyramid (256x256 PNG tiles), and creates a TissUUmaps project
    file (.tmap) with per-channel color mapping.

    Also copies segmentation masks, CSV overlays, and GeoJSON regions if present.

    Re-runs skip unchanged files unless --force is used.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow export prepare .              # Export current project
        kintsugi workflow export prepare . --dry-run    # Preview what would be exported
        kintsugi workflow export prepare . --force      # Reconvert everything
    """
    from pathlib import Path

    from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn

    from kintsugi.export import prepare_export

    project_dir = Path(project_dir).resolve()

    if not (project_dir / "kintsugi_project.json").exists():
        # Also allow projects that have registered/ without kintsugi_project.json
        if not (project_dir / "data" / "processed" / "registered").exists():
            console.print(f"[red]Not a KINTSUGI project: {project_dir}[/red]")
            console.print("No kintsugi_project.json or data/processed/registered/ found.")
            raise SystemExit(1)

    console.print(f"[bold]TissUUmaps Export[/bold]  ({project_dir.name})\n")

    if dry_run:
        result = prepare_export(
            project_dir,
            dry_run=True,
            exclude_blanks=not include_blanks,
        )

        if result["status"] == "error":
            console.print(f"[red]{result['message']}[/red]")
            raise SystemExit(1)

        console.print(f"  Export directory: {result['export_dir']}")
        console.print(
            f"  Registered: {result['registered_images']} images across {result['registered_cycles']} cycles"
        )
        console.print(f"  Segmentation masks: {result['segmentation_masks']}")
        console.print(f"  Overlays: {result['overlays']}")
        console.print(f"  Total DZI conversions: {result['total_dzi_conversions']}")

        if result.get("cycles"):
            console.print("\n  [bold]Channels per cycle:[/bold]")
            for cycle, markers in result["cycles"].items():
                console.print(f"    {cycle}: {', '.join(markers)}")

        if result.get("overlay_files"):
            console.print(f"\n  [bold]Overlay files:[/bold] {', '.join(result['overlay_files'])}")

        if result.get("masks"):
            console.print(f"  [bold]Masks:[/bold] {', '.join(result['masks'])}")
        return

    # Real export with progress bar
    progress = Progress(
        TextColumn("[bold]{task.description}"),
        BarColumn(),
        TextColumn("{task.completed}/{task.total}"),
        TimeRemainingColumn(),
        console=console,
    )

    task_id = None

    def on_progress(current: int, total: int, message: str) -> None:
        nonlocal task_id
        if task_id is None:
            task_id = progress.add_task("Converting", total=total)
        progress.update(task_id, completed=current, description=message)

    with progress:
        result = prepare_export(
            project_dir,
            force=force,
            exclude_blanks=not include_blanks,
            progress_callback=on_progress,
        )

    if result["status"] == "error":
        console.print(f"\n[red]{result['message']}[/red]")
        raise SystemExit(1)

    console.print("\n[green]Export complete![/green]")
    console.print(f"  Images converted: {result['images_converted']}")
    if result["images_skipped"]:
        console.print(f"  Images skipped (unchanged): {result['images_skipped']}")
    if result["segmentation_masks"]:
        console.print(f"  Segmentation masks: {result['segmentation_masks']}")
    if result["overlays_copied"]:
        console.print(f"  Overlays copied: {result['overlays_copied']}")
    if result["regions_copied"]:
        console.print(f"  Regions copied: {result['regions_copied']}")
    console.print(f"  .tmap file: {result['tmap_path']}")
    console.print(
        f"\n[dim]Deploy with: kintsugi workflow export deploy {project_dir} user@host:/path[/dim]"
    )


@export_group.command("deploy")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.argument("remote", type=str)
@click.option("--dry-run", "-n", is_flag=True, help="Preview rsync transfer without sending")
@click.option("--bwlimit", type=int, default=None, help="Bandwidth limit in KB/s")
@click.option("--ssh-key", type=click.Path(exists=True), default=None, help="SSH private key path")
def export_deploy(
    project_dir: str, remote: str, dry_run: bool, bwlimit: int | None, ssh_key: str | None
):
    """
    Deploy exported data to a remote server via rsync.

    REMOTE is the rsync destination (e.g., user@host:/var/www/tissuumaps/projects/).

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow export deploy . user@server:/path
        kintsugi workflow export deploy . user@server:/path --dry-run
        kintsugi workflow export deploy . user@server:/path --bwlimit 10000
    """
    from pathlib import Path

    from kintsugi.export import deploy_export

    project_dir = Path(project_dir).resolve()
    export_dir = project_dir / "data" / "exported"

    if not export_dir.exists():
        console.print(f"[red]No export directory found: {export_dir}[/red]")
        console.print("Run 'kintsugi workflow export prepare .' first.")
        raise SystemExit(1)

    console.print(f"[bold]Deploying export[/bold]  ({project_dir.name})")
    console.print(f"  Source: {export_dir}")
    console.print(f"  Destination: {remote}")
    if dry_run:
        console.print("  [yellow]DRY RUN — no files will be transferred[/yellow]")
    console.print()

    rc = deploy_export(
        export_dir,
        remote,
        dry_run=dry_run,
        bwlimit=bwlimit,
        ssh_key=ssh_key,
    )

    if rc == 0:
        if not dry_run:
            console.print("\n[green]Deploy complete![/green]")
    else:
        console.print(f"\n[red]rsync exited with code {rc}[/red]")
        raise SystemExit(rc)


@export_group.command("status")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def export_status(project_dir: str):
    """
    Show export status for a project.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from pathlib import Path

    from kintsugi.export import get_export_status

    project_dir = Path(project_dir).resolve()

    console.print(f"[bold]Export Status[/bold]  ({project_dir.name})\n")

    status = get_export_status(project_dir)

    if status["status"] == "not_exported":
        console.print("[dim]No export found. Run 'kintsugi workflow export prepare .' first.[/dim]")
        return

    console.print(f"  Export directory: {status['export_dir']}")
    console.print(f"  DZI images: {status['dzi_count']}")
    console.print(f"  Total size: {status['total_size_mb']} MB")

    if status.get("tmap_files"):
        console.print(f"  .tmap files: {', '.join(status['tmap_files'])}")

    if status.get("manifest_entries"):
        console.print(f"  Manifest entries: {status['manifest_entries']}")

    if status.get("last_updated"):
        console.print(f"  Last updated: {status['last_updated']}")

    if status.get("overlay_files"):
        console.print(f"  Overlay files: {', '.join(status['overlay_files'])}")

    if status.get("geojson_regions"):
        console.print(f"  GeoJSON regions: {status['geojson_regions']}")


# ============================================================================
# Project Commands
# ============================================================================


@main.command()
@click.argument("project_path", type=click.Path())
@click.option("--name", "-n", help="Project name")
@click.option("--description", "-d", default="", help="Project description")
@click.option(
    "--force", "-f", is_flag=True, help="Skip confirmation prompts and directory scanning"
)
@click.option("--slurm", is_flag=True, help="Initialize SLURM job submission support")
@click.option("--slurm-account", help="HPC account for SLURM (auto-detected if not provided)")
@click.option("--slurm-partition", default=None, help="SLURM partition (default: gpu)")
@click.option("--slurm-qos", default=None, help="SLURM Quality of Service")
@click.option("--slurm-gpu-type", default=None, help="GPU type for SLURM (auto-detected)")
# Microscope parameters (stored in /meta/experiment.json)
@click.option("--tile-rows", type=int, default=None, help="Number of tile rows (auto-detected)")
@click.option("--tile-cols", type=int, default=None, help="Number of tile columns (auto-detected)")
@click.option(
    "--xy-pixel-size", type=float, default=377.0, help="XY pixel size in nm (default: 377)"
)
@click.option("--z-step-size", type=float, default=1500.0, help="Z step size in nm (default: 1500)")
@click.option("--numerical-aperture", type=float, default=0.75, help="Objective NA (default: 0.75)")
@click.option(
    "--tissue-ri", type=float, default=1.44, help="Tissue refractive index (default: 1.44)"
)
def init(
    project_path: str,
    name: str | None,
    description: str,
    force: bool,
    slurm: bool,
    slurm_account: str | None,
    slurm_partition: str | None,
    slurm_qos: str | None,
    slurm_gpu_type: str | None,
    tile_rows: int | None,
    tile_cols: int | None,
    xy_pixel_size: float,
    z_step_size: float,
    numerical_aperture: float,
    tissue_ri: float,
):
    """
    Initialize a new KINTSUGI project.

    Creates the project directory structure and configuration.

    SAFETY: This command will NEVER delete existing files. If data exists in the
    directory, you will be prompted to review it before proceeding.

    Use --force to skip directory scanning (faster for large existing datasets).

    Use --slurm to add SLURM job submission support (creates slurm/ directory
    with config.sh, job script symlinks, and README).
    """
    from pathlib import Path

    from kintsugi.project import ExistingDataReport, KintsugiProject, scan_existing_data

    project_path = Path(project_path).resolve()
    config_file = project_path / "kintsugi_project.json"

    # Handle existing projects
    if config_file.exists():
        project = KintsugiProject.load(project_path)
        slurm_dir = project_path / "slurm"

        if force:
            # Fast path: refresh all assets from the current repo
            console.print(f"\n[bold]Updating existing project:[/bold] {project_path}")
            project._kintsugi_path = project._detect_kintsugi_path()
            project.config.kintsugi_path = str(project._kintsugi_path)
            project.paths.create_all()
            project.setup_notebooks(overwrite=True)
            project._create_claude_config()
            project._create_vscode_config()
            if slurm:
                project.setup_slurm(
                    account=slurm_account,
                    partition=slurm_partition,
                    qos=slurm_qos,
                    gpu_type=slurm_gpu_type,
                )
            project.save()
            console.print("\n[green]Project updated from KINTSUGI templates.[/green]")
            return

        # Check if SLURM was requested but not set up
        if slurm and not slurm_dir.exists():
            console.print(f"\n[bold]Project exists:[/bold] {project_path}")
            console.print("[yellow]SLURM support not configured. Adding now...[/yellow]")
            project.setup_slurm(
                account=slurm_account,
                partition=slurm_partition,
                qos=slurm_qos,
                gpu_type=slurm_gpu_type,
            )
            project.save()
            console.print("\n[green]SLURM support added to existing project.[/green]")
            return

        # Project exists, no special action needed
        console.print(f"\n[bold]Project already exists:[/bold] {project_path}")
        if slurm_dir.exists():
            console.print("[dim]SLURM support is already configured.[/dim]")
        else:
            console.print(
                "[dim]Use --slurm to add SLURM support, or run: kintsugi slurm init .[/dim]"
            )
        console.print("[dim]Use --force to refresh notebooks and configs from KINTSUGI.[/dim]")
        return

    try:
        # Skip scanning if --force is used (much faster for large datasets)
        if force:
            console.print(f"\n[bold]Initializing project:[/bold] {project_path}")
            console.print("[dim]Skipping directory scan (--force)[/dim]")
            report = ExistingDataReport()  # Empty report
        else:
            # Scan for existing data first
            console.print(f"\n[bold]Scanning directory:[/bold] {project_path}")
            report = scan_existing_data(project_path)

        # Track which processed stages to delete (if any)
        stages_to_delete: list[str] = []

        if report.has_data:
            # Build data summary based on raw vs processed
            summary_lines = []

            if report.has_raw_data:
                summary_lines.append(
                    f"[green]Raw data:[/green] {report.raw_image_count} images "
                    f"({report.raw_size_mb:.1f} MB)"
                )
                if report.raw_cycle_folders:
                    summary_lines.append(
                        f"  Cycles: {', '.join(report.raw_cycle_folders[:10])}"
                        + ("..." if len(report.raw_cycle_folders) > 10 else "")
                    )

            if report.has_processed_data:
                summary_lines.append(
                    f"[yellow]Processed data:[/yellow] {report.processed_size_mb:.1f} MB"
                )
                for stage, count in report.processed_stages.items():
                    summary_lines.append(f"  {stage}: {count} files")

            # Determine panel style based on data type
            if report.has_processed_data:
                panel_title = "Existing Processed Data Found"
                panel_style = "yellow"
            else:
                panel_title = "Raw Data Found"
                panel_style = "green"

            console.print(
                Panel(
                    "\n".join(summary_lines),
                    title=panel_title,
                    border_style=panel_style,
                )
            )

            # Show metadata from sample images
            if report.metadata_samples:
                console.print("\n[bold]Sample Image Metadata:[/bold]")
                meta_table = Table(show_header=True, header_style="bold")
                meta_table.add_column("File")
                meta_table.add_column("Dimensions")
                meta_table.add_column("Channels")
                meta_table.add_column("Pixel Size")
                meta_table.add_column("OME")

                for sample in report.metadata_samples[:5]:
                    dims = str(sample.get("dimensions", "N/A"))
                    if len(dims) > 25:
                        dims = dims[:22] + "..."
                    ch_names = sample.get("channel_names", [])
                    ch_str = str(sample.get("channels") or "N/A")
                    if ch_names:
                        ch_str = f"{len(ch_names)}: {', '.join(ch_names[:3])}"
                        if len(ch_names) > 3:
                            ch_str += "..."
                    px = sample.get("pixel_size")
                    px_str = f"{px:.4f} {sample.get('pixel_unit', 'um')}" if px else "N/A"

                    meta_table.add_row(
                        sample.get("file", "?")[:30],
                        dims,
                        ch_str,
                        px_str,
                        "Yes" if sample.get("is_ome") else "No",
                    )
                console.print(meta_table)

            if not force:
                # Different options depending on what data exists
                if report.has_processed_data:
                    # Processed data exists - offer delete or keep options
                    console.print("\n[bold]Options:[/bold]")
                    console.print(
                        "  1. Delete processed - Remove all processed data and start fresh"
                    )
                    console.print(
                        "  2. Keep processed - Initialize project keeping existing processed data"
                    )
                    console.print("  3. Cancel - Exit without making changes")

                    choice = click.prompt(
                        "\nSelect option",
                        type=click.Choice(["1", "2", "3", "delete", "keep", "cancel"]),
                        default="2",
                    )

                    if choice in ["3", "cancel"]:
                        console.print("[yellow]Cancelled. No changes made.[/yellow]")
                        return

                    if choice in ["1", "delete"]:
                        # Ask which stages to delete
                        stages = list(report.processed_stages.keys())
                        if len(stages) == 1:
                            stages_to_delete = stages
                            console.print(
                                f"[yellow]Will delete {stages[0]} data after project creation.[/yellow]"
                            )
                        else:
                            console.print(
                                "\n[bold]Select stages to delete (comma-separated, or 'all'):[/bold]"
                            )
                            for i, stage in enumerate(stages, 1):
                                count = report.processed_stages[stage]
                                console.print(f"  {i}. {stage} ({count} files)")

                            stage_choice = click.prompt(
                                "\nStages to delete",
                                default="all",
                            )

                            if stage_choice.lower() == "all":
                                stages_to_delete = stages
                            else:
                                # Parse comma-separated list of numbers or names
                                for item in stage_choice.split(","):
                                    item = item.strip()
                                    if item.isdigit():
                                        idx = int(item) - 1
                                        if 0 <= idx < len(stages):
                                            stages_to_delete.append(stages[idx])
                                    elif item in stages:
                                        stages_to_delete.append(item)

                            if stages_to_delete:
                                console.print(
                                    f"[yellow]Will delete: {', '.join(stages_to_delete)}[/yellow]"
                                )
                else:
                    # Only raw data - simple continue/cancel
                    console.print(
                        "\n[bold green]Raw data found in expected location.[/bold green] "
                        "Project structure will be created around it."
                    )
                    console.print("\n[bold]Options:[/bold]")
                    console.print("  1. Continue - Create project structure")
                    console.print("  2. Cancel - Exit without making changes")

                    choice = click.prompt(
                        "\nSelect option",
                        type=click.Choice(["1", "2", "continue", "cancel"]),
                        default="1",
                    )

                    if choice in ["2", "cancel"]:
                        console.print("[yellow]Cancelled. No changes made.[/yellow]")
                        return

        # Create the project
        KintsugiProject.create(
            project_path,
            name=name,
            description=description,
            existing_data_report=report,
            slurm=slurm,
            slurm_account=slurm_account,
            slurm_partition=slurm_partition,
            slurm_qos=slurm_qos,
            slurm_gpu_type=slurm_gpu_type,
            # Microscope parameters
            tile_rows=tile_rows,
            tile_cols=tile_cols,
            xy_pixel_size=xy_pixel_size,
            z_step_size=z_step_size,
            numerical_aperture=numerical_aperture,
            tissue_refractive_index=tissue_ri,
        )

        # Delete processed stages if requested
        if stages_to_delete:
            import shutil

            processed_dir = project_path / "data" / "processed"
            for stage in stages_to_delete:
                stage_dir = processed_dir / stage
                if stage_dir.exists():
                    console.print(f"  Deleting {stage}...")
                    shutil.rmtree(stage_dir)
            console.print("[green]Processed data deleted.[/green]")

        console.print("\n[green]Project created successfully![/green]")

        # Remind about Claude Code
        console.print(
            "\n[dim]Claude Code is configured. Open this folder in VS Code to start.[/dim]"
        )

    except Exception as e:
        console.print(f"[red]Failed to create project: {e}[/red]")
        raise SystemExit(1)


@main.command()
@click.argument("directory", type=click.Path(exists=True))
@click.option("--depth", "-d", default=3, help="Maximum scan depth")
@click.option("--samples", "-s", default=5, help="Number of images to sample for metadata")
def scan(directory: str, depth: int, samples: int):
    """
    Scan a directory for image data and extract metadata.

    Useful for inspecting data before initializing a project.
    Shows image counts, dimensions, channel names, and pixel sizes.
    """
    from pathlib import Path

    from kintsugi.project import scan_existing_data

    directory = Path(directory).resolve()

    console.print(f"\n[bold]Scanning:[/bold] {directory}")
    console.print(f"[dim]Depth: {depth}, Samples: {samples}[/dim]\n")

    report = scan_existing_data(directory, max_depth=depth, sample_count=samples)

    if not report.has_data:
        console.print("[yellow]No image data found.[/yellow]")
        return

    # Build summary based on raw vs processed data
    summary_lines = []

    if report.has_raw_data:
        summary_lines.append(
            f"[green]Raw data:[/green] {report.raw_image_count} images "
            f"({report.raw_size_mb:.1f} MB)"
        )
        if report.raw_cycle_folders:
            cycles_str = ", ".join(report.raw_cycle_folders[:10])
            if len(report.raw_cycle_folders) > 10:
                cycles_str += f"... (+{len(report.raw_cycle_folders) - 10})"
            summary_lines.append(f"  Cycles: {cycles_str}")

    if report.has_processed_data:
        summary_lines.append(f"[yellow]Processed data:[/yellow] {report.processed_size_mb:.1f} MB")
        for stage, count in report.processed_stages.items():
            summary_lines.append(f"  {stage}: {count} files")

    # Add overall totals
    summary_lines.append("")
    summary_lines.append(f"Total: {report.image_count} image files ({report.total_size_mb:.1f} MB)")
    if report.filename_patterns:
        summary_lines.append(
            f"Patterns: {', '.join(report.filename_patterns[:5])}"
            + ("..." if len(report.filename_patterns) > 5 else "")
        )

    # Determine panel style
    if report.has_processed_data:
        panel_title = "Data Summary (Raw + Processed)"
        panel_style = "yellow"
    elif report.has_raw_data:
        panel_title = "Data Summary (Raw)"
        panel_style = "green"
    else:
        panel_title = "Data Summary"
        panel_style = "blue"

    console.print(
        Panel(
            "\n".join(summary_lines),
            title=panel_title,
            border_style=panel_style,
        )
    )

    # Detailed metadata table
    if report.metadata_samples:
        console.print("\n[bold]Image Metadata (sampled):[/bold]")
        meta_table = Table(show_header=True, header_style="bold")
        meta_table.add_column("File", max_width=35)
        meta_table.add_column("Dimensions")
        meta_table.add_column("Type")
        meta_table.add_column("Channels")
        meta_table.add_column("Pixel Size")
        meta_table.add_column("OME")

        for sample in report.metadata_samples:
            dims = str(sample.get("dimensions", "N/A"))
            ch_names = sample.get("channel_names", [])
            ch_str = str(sample.get("channels") or "N/A")
            if ch_names:
                ch_str = f"{len(ch_names)}: {', '.join(ch_names[:3])}"
                if len(ch_names) > 3:
                    ch_str += f"... (+{len(ch_names) - 3})"
            px = sample.get("pixel_size")
            px_str = f"{px:.4f} {sample.get('pixel_unit', 'um')}" if px else "N/A"

            meta_table.add_row(
                sample.get("file", "?"),
                dims,
                sample.get("dtype", "N/A") if "dtype" in sample else "N/A",
                ch_str,
                px_str,
                "Yes" if sample.get("is_ome") else "No",
            )
        console.print(meta_table)

        # Show channel names if OME
        all_channels = set()
        for sample in report.metadata_samples:
            all_channels.update(sample.get("channel_names", []))
        if all_channels:
            console.print(f"\n[bold]Channel names found:[/bold] {', '.join(sorted(all_channels))}")

    console.print(f"\n[dim]To create a project here: kintsugi init {directory}[/dim]")


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
