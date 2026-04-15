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

import click
from rich.console import Console
from rich.panel import Panel
from rich.table import Table

from kintsugi.deps import OPTIONAL_GROUPS, _find_project_root

console = Console()


def _resolve_project_dir(path: str | Path) -> Path:
    """Resolve project directory, auto-detecting if inside a workflow/ subdirectory.

    If the resolved path is named "workflow", has a config.yaml, and its parent
    has a meta/ directory, return the parent (the actual project root).
    """
    resolved = Path(path).resolve()
    if (
        resolved.name == "workflow"
        and (resolved / "config.yaml").exists()
        and (resolved.parent / "meta").is_dir()
    ):
        return resolved.parent
    return resolved


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
@click.option("--strict", is_flag=True, help="Exit code 1 if any errors detected")
@click.option(
    "--for",
    "workflow_name",
    type=click.Choice(
        ["pipeline", "analysis", "segmentation", "vessel3d", "kronos", "all"],
        case_sensitive=False,
    ),
    default=None,
    help="Validate deps for a specific workflow",
)
def check(verbose: bool, output_json: bool, strict: bool, workflow_name: str | None):
    """Check all KINTSUGI dependencies."""
    from kintsugi.deps import DependencyChecker, DependencyStatus, require

    # If --for is specified, check workflow-specific deps first
    if workflow_name:
        workflow_groups = {
            "pipeline": ["gpu"],
            "analysis": ["analysis", "dl", "viz"],
            "segmentation": ["dl", "viz"],
            "vessel3d": ["gpu", "analysis"],
            "kronos": ["kronos"],
            "all": ["gpu", "viz", "dl", "analysis"],
        }
        groups = workflow_groups.get(workflow_name, [])
        if groups:
            try:
                require(*groups, strict=True)
                if not output_json:
                    console.print(
                        f"[green]All dependencies for '{workflow_name}' workflow are installed.[/green]"
                    )
            except Exception as e:
                if not output_json:
                    console.print(f"[red]{e}[/red]")
                if strict:
                    raise SystemExit(1)

    checker = DependencyChecker()

    if output_json:
        result = checker.check_all(verbose=False)
        click.echo(json.dumps(result, indent=2))
    else:
        checker.check_all(verbose=verbose)

    if strict:
        errors = [r for r in checker.results if r.status == DependencyStatus.ERROR]
        if errors:
            if not output_json:
                console.print(f"\n[red]Strict check failed: {len(errors)} error(s) found[/red]")
                for err in errors:
                    console.print(f"  [red]{err.name}: {err.message}[/red]")
            raise SystemExit(1)


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
@click.option("--conda", "use_conda", is_flag=True, help="Use conda instead of pip where available")
def install(group: str | None, list_groups: bool, use_conda: bool):
    """
    Install optional dependency groups.

    GROUP is the name of the dependency group to install.
    Use 'all' to install all optional dependencies.

    \b
    Examples:
        kintsugi install --list     # Show available groups
        kintsugi install gpu        # Install GPU acceleration
        kintsugi install analysis   # Install spatial analysis tools
        kintsugi install all        # Install everything
        kintsugi install gpu --conda  # Use conda for GPU packages
    """
    if list_groups or group is None:
        table = Table(title="Optional Dependency Groups")
        table.add_column("Group", style="cyan")
        table.add_column("Description", style="white")
        table.add_column("Packages", style="dim")

        for name, info in OPTIONAL_GROUPS.items():
            if name == "full":
                continue  # Skip composite group in listing
            table.add_row(name, info["description"], ", ".join(info["packages"]))

        console.print(table)
        console.print("\n[dim]Usage: kintsugi install <group>[/dim]")
        console.print("[dim]       kintsugi install all[/dim]")
        return

    if group == "all":
        _install_all(use_conda)
        return

    if group not in OPTIONAL_GROUPS:
        console.print(f"[red]Unknown group: {group}[/red]")
        available = [k for k in OPTIONAL_GROUPS if k != "full"]
        console.print(f"Available groups: {', '.join(available)}, all")
        raise SystemExit(1)

    # Single-group install
    from kintsugi.deps import _inject_constraints, _pre_install_guard

    # Pre-install guard checks
    guard_warnings = _pre_install_guard(group)
    for warning in guard_warnings:
        console.print(f"\n[yellow]Warning: {warning}[/yellow]")

    info = OPTIONAL_GROUPS[group]
    console.print(f"\n[bold]Installing {group}:[/bold] {info['description']}")

    if use_conda and "conda_cmd" in info:
        cmd = info["conda_cmd"]
    else:
        cmd = info["install_cmd"]
        # Inject constraints file for pip commands
        cmd = _inject_constraints(cmd)

    console.print(f"[dim]Running: {cmd}[/dim]")

    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=False)
        console.print(f"[green]✓ {group} installed successfully[/green]")
    except subprocess.CalledProcessError as e:
        console.print(f"[red]✗ Failed to install {group}: {e}[/red]")
        raise SystemExit(1)

    if "note" in info:
        console.print(f"[yellow]Note: {info['note']}[/yellow]")

    # Post-install validation
    if not _post_install_validate():
        console.print(
            "\n[red bold]Installation completed with ERRORS. Environment may be broken.[/red bold]"
        )
        raise SystemExit(1)

    console.print("\n[bold green]Installation complete![/bold green]")


def _post_install_validate() -> bool:
    """Run post-install validation. Returns False if critical errors found."""
    from kintsugi.deps import DependencyChecker, DependencyStatus

    console.print("\n[dim]Validating environment...[/dim]")
    checker = DependencyChecker()
    checker.check_all(verbose=False)

    errors = [r for r in checker.results if r.status == DependencyStatus.ERROR]
    missing_required = [
        r
        for r in checker.results
        if not r.is_optional and r.status in (DependencyStatus.MISSING, DependencyStatus.ERROR)
    ]

    if missing_required:
        console.print("[red bold]CRITICAL: Required dependencies broken after install:[/red bold]")
        for err in missing_required:
            msg = err.message or "not found"
            console.print(f"  [red]✗ {err.name}: {msg}[/red]")
        console.print("[dim]Run 'kintsugi check' for full details[/dim]")
        return False

    if errors:
        console.print("[yellow]Post-install warnings:[/yellow]")
        for err in errors:
            console.print(f"  [yellow]{err.name}: {err.message}[/yellow]")
        console.print("[dim]Run 'kintsugi check --strict' for full details[/dim]")

    return True


def _install_all(use_conda: bool) -> None:
    """Install all optional dependency groups with unified resolution.

    Uses pyproject.toml extras in a single pip pass to avoid cascading
    dependency conflicts from sequential per-group pip invocations.
    Falls back to sequential install + constraint repair if pyproject.toml
    is not found (e.g., wheel install without source checkout).

    On HPC (detected automatically), switches to conda mode to ensure
    CUDA runtime libraries and conda-compiled packages are installed
    correctly via conda channels rather than pip.
    """
    from kintsugi.deps import _inject_constraints, detect_hpc

    # Auto-detect HPC and switch to conda mode
    if detect_hpc() and not use_conda:
        console.print(
            "[yellow]HPC environment detected. Switching to conda mode "
            "for GPU/CUDA packages.[/yellow]"
        )
        use_conda = True

    # 'docs' is excluded from 'all' because Sphinx/myst-parser are only
    # needed to build documentation, not to run KINTSUGI. Install with
    # `kintsugi install docs` if you need them.
    _SKIP_FROM_ALL = {"full", "rapids", "docs"}
    extras = [k for k in OPTIONAL_GROUPS if k not in _SKIP_FROM_ALL]
    project_root = _find_project_root()
    failed_groups: list[str] = []

    if project_root and not use_conda:
        # Pre-install pims to work around setuptools build failure (needed by dask-image)
        from kintsugi.deps import _pre_install_pims

        _pre_install_pims()

        # Single pip resolution pass using pyproject.toml extras
        extras_str = ",".join(extras)
        cmd = f'pip install -e "{project_root}[{extras_str}]"'
        console.print(f"\n[bold]Installing all groups via extras:[/bold] {extras_str}")
        console.print(f"[dim]Running: {cmd}[/dim]")
        try:
            subprocess.run(cmd, shell=True, check=True, capture_output=False)
            console.print("[green]✓ All groups installed successfully[/green]")
        except subprocess.CalledProcessError as e:
            console.print(f"[red]✗ Failed: {e}[/red]")
            raise SystemExit(1)
    elif use_conda:
        # Conda mode: run conda_cmd groups first, then pip-only groups individually.
        # IMPORTANT: We do NOT use pyproject.toml extras for the pip phase because
        # that would re-install torch from PyPI, overwriting the conda version.
        conda_groups = [g for g in extras if "conda_cmd" in OPTIONAL_GROUPS[g]]
        pip_only_groups = [g for g in extras if g not in conda_groups]

        # Phase 1: Install all conda groups
        for grp in conda_groups:
            info = OPTIONAL_GROUPS[grp]
            console.print(f"\n[bold]Installing {grp} (conda):[/bold] {info['description']}")
            cmd = info["conda_cmd"]
            console.print(f"[dim]Running: {cmd}[/dim]")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=False)
                console.print(f"[green]✓ {grp} installed[/green]")
            except subprocess.CalledProcessError as e:
                console.print(f"[red]✗ FAILED to install {grp}: {e}[/red]")
                failed_groups.append(grp)

        # Phase 2: Install pip-only groups individually (not via pyproject extras
        # to avoid re-pulling torch from PyPI and overwriting conda version)
        for grp in pip_only_groups:
            info = OPTIONAL_GROUPS[grp]
            cmd = _inject_constraints(info["install_cmd"])
            console.print(f"\n[bold]Installing {grp} (pip):[/bold] {info['description']}")
            console.print(f"[dim]Running: {cmd}[/dim]")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=False)
                console.print(f"[green]✓ {grp} installed[/green]")
            except subprocess.CalledProcessError as e:
                console.print(f"[red]✗ FAILED to install {grp}: {e}[/red]")
                failed_groups.append(grp)

        if failed_groups:
            console.print(
                f"\n[red bold]Installation FAILED for groups: {', '.join(failed_groups)}[/red bold]"
            )
    else:
        # Fallback: sequential install + repair (no pyproject.toml found)
        console.print(
            "[yellow]Warning: pyproject.toml not found. "
            "Installing sequentially (may cause conflicts).[/yellow]"
        )
        for grp in extras:
            info = OPTIONAL_GROUPS[grp]
            console.print(f"\n[bold]Installing {grp}:[/bold] {info['description']}")
            cmd = info["install_cmd"]
            console.print(f"[dim]Running: {cmd}[/dim]")
            try:
                subprocess.run(cmd, shell=True, check=True, capture_output=False)
                console.print(f"[green]✓ {grp} installed[/green]")
            except subprocess.CalledProcessError as e:
                console.print(f"[red]✗ FAILED to install {grp}: {e}[/red]")
                failed_groups.append(grp)

        # Repair core constraints
        console.print("\n[bold]Repairing core constraints...[/bold]")
        subprocess.run("pip install -e .", shell=True, check=False)

    # Print notes for all groups
    for grp in extras:
        if "note" in OPTIONAL_GROUPS[grp]:
            console.print(f"\n[yellow]Note ({grp}): {OPTIONAL_GROUPS[grp]['note']}[/yellow]")

    # Auto-apply SLURM TRES patch if jobstep plugin is installed
    _auto_patch_slurm_tres()

    # Post-install validation (fails loudly on critical errors)
    validation_ok = _post_install_validate()

    # Remind about rapids (excluded from 'all')
    console.print(
        "\n[dim]Note: 'rapids' requires separate installation "
        "(NVIDIA channel). Run: kintsugi install rapids[/dim]"
    )

    # Remind about docs (excluded from 'all')
    console.print(
        "[dim]Note: 'docs' is excluded from 'all' (build-time only). "
        "Run: kintsugi install docs[/dim]"
    )

    if failed_groups or not validation_ok:
        console.print(
            "\n[red bold]Installation completed with ERRORS. Environment may be broken.[/red bold]"
        )
        console.print("[dim]Run 'kintsugi check' for details[/dim]")
        raise SystemExit(1)

    console.print("\n[bold green]Installation complete![/bold green]")


def _auto_patch_slurm_tres() -> None:
    """Automatically apply the SLURM_TRES_PER_TASK patch if needed."""
    try:
        import snakemake_executor_plugin_slurm_jobstep

        init_file = Path(snakemake_executor_plugin_slurm_jobstep.__file__)
        source = init_file.read_text()
        if "SLURM_TRES_PER_TASK" not in source:
            console.print(
                "\n[yellow]Applying SLURM_TRES_PER_TASK patch to jobstep plugin...[/yellow]"
            )
            _apply_tres_patch(init_file, source)
    except ImportError:
        pass  # Not installed
    except Exception:
        pass  # Can't read/write


def _apply_tres_patch(init_file: Path, source: str) -> bool:
    """Apply the SLURM_TRES_PER_TASK patch to the jobstep plugin.

    Inserts os.environ.pop("SLURM_TRES_PER_TASK", None) into the
    __post_init__ method of the executor class.
    """
    # Find __post_init__ method
    patch_line = '        os.environ.pop("SLURM_TRES_PER_TASK", None)  # KINTSUGI patch\n'

    # Look for the __post_init__ method and insert after it
    lines = source.splitlines(keepends=True)
    for i, line in enumerate(lines):
        if "def __post_init__" in line:
            # Insert after the next line (which is usually the docstring or first statement)
            insert_at = i + 1
            # Skip docstring if present
            while insert_at < len(lines) and (
                lines[insert_at].strip().startswith('"""')
                or lines[insert_at].strip().startswith("'''")
                or lines[insert_at].strip() == ""
            ):
                insert_at += 1
                # Handle multi-line docstrings
                if insert_at > 0 and (
                    lines[insert_at - 1].strip().endswith('"""')
                    or lines[insert_at - 1].strip().endswith("'''")
                ):
                    break

            # Also need 'import os' at the top
            if "import os" not in source:
                lines.insert(0, "import os\n")
                insert_at += 1

            lines.insert(insert_at, patch_line)
            init_file.write_text("".join(lines))
            console.print("[green]✓ SLURM_TRES_PER_TASK patch applied[/green]")
            return True

    console.print("[yellow]Could not locate __post_init__ method in jobstep plugin[/yellow]")
    return False


@main.command("patch")
@click.argument("target", type=click.Choice(["slurm"]))
def patch_command(target: str):
    """Apply patches to third-party plugins.

    \b
    Targets:
        slurm   Apply SLURM_TRES_PER_TASK fix to snakemake-executor-plugin-slurm-jobstep
    """
    if target == "slurm":
        try:
            import snakemake_executor_plugin_slurm_jobstep

            init_file = Path(snakemake_executor_plugin_slurm_jobstep.__file__)
            source = init_file.read_text()

            if "SLURM_TRES_PER_TASK" in source:
                console.print("[green]SLURM_TRES_PER_TASK patch is already applied.[/green]")
                return

            if _apply_tres_patch(init_file, source):
                console.print(
                    "\n[dim]This patch must be re-applied after upgrading "
                    "snakemake-executor-plugin-slurm-jobstep.[/dim]"
                )
            else:
                raise SystemExit(1)

        except ImportError:
            console.print(
                "[red]snakemake-executor-plugin-slurm-jobstep not installed.[/red]\n"
                "[dim]Install with: kintsugi install workflow[/dim]"
            )
            raise SystemExit(1)


@main.command("fix-hpc")
def fix_hpc_command():
    """Fix common HPC environment issues (libstdc++, CUDA paths).

    \b
    Deploys conda activation scripts that:
      - Prepend $CONDA_PREFIX/lib to LD_LIBRARY_PATH (fixes libstdc++ mismatch)
      - Set CUDA_PATH/CUDA_HOME for CuPy JIT compilation

    Also copies CUDA headers and installs libstdcxx-ng if needed.

    \b
    After running this command, deactivate and reactivate the environment:
      conda deactivate && conda activate KINTSUGI
    """
    import os

    from kintsugi.deps import deploy_activation_scripts, detect_hpc

    if not detect_hpc():
        console.print("[yellow]Warning: HPC environment not detected. Running anyway.[/yellow]")

    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    if not conda_prefix:
        console.print("[red]Error: No conda environment is active.[/red]")
        raise SystemExit(1)

    # Step 1: Deploy activation scripts
    console.print("\n[bold]Step 1: Deploying LD_LIBRARY_PATH activation scripts[/bold]")
    if deploy_activation_scripts():
        console.print("[green]✓ Activation scripts deployed[/green]")
    else:
        console.print("[red]✗ Failed to deploy activation scripts[/red]")
        console.print(
            "[dim]Ensure KINTSUGI is installed in dev mode "
            "(pip install -e .) and envs/activate.d/ exists[/dim]"
        )

    # Step 2: Install libstdcxx-ng if needed
    console.print("\n[bold]Step 2: Ensuring conda libstdcxx-ng is up to date[/bold]")
    try:
        subprocess.run(
            [
                "conda",
                "install",
                "libstdcxx-ng>=12.0",
                "libgcc-ng>=12.0",
                "-c",
                "conda-forge",
                "-y",
            ],
            check=True,
            capture_output=False,
        )
        console.print("[green]✓ libstdcxx-ng up to date[/green]")
    except subprocess.CalledProcessError:
        console.print("[yellow]Warning: Could not update libstdcxx-ng[/yellow]")

    # Step 3: Copy CUDA headers for CuPy JIT compilation
    console.print("\n[bold]Step 3: Copying CUDA headers[/bold]")
    targets_include = Path(conda_prefix) / "targets" / "x86_64-linux" / "include"
    conda_include = Path(conda_prefix) / "include"
    if targets_include.is_dir():
        try:
            import shutil as _shutil

            conda_include.mkdir(parents=True, exist_ok=True)
            for item in targets_include.iterdir():
                dst = conda_include / item.name
                if item.is_dir():
                    _shutil.copytree(item, dst, dirs_exist_ok=True)
                else:
                    _shutil.copy2(item, dst)
            console.print("[green]✓ CUDA headers copied[/green]")
        except OSError as e:
            console.print(f"[yellow]Warning: Could not copy CUDA headers: {e}[/yellow]")
    else:
        console.print("[dim]CUDA headers not found (cuda-cudart-dev may not be installed)[/dim]")

    # Instructions
    console.print("\n[bold yellow]IMPORTANT: Deactivate and reactivate to apply:[/bold yellow]")
    console.print("  conda deactivate && conda activate KINTSUGI")
    console.print("\n[dim]Then verify: kintsugi check[/dim]")


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

    Run 'kintsugi mcp config <project_path>' to generate the .mcp.json
    configuration file, or use 'kintsugi init' which creates it automatically.
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

    Creates .mcp.json (MCP server definition) at the project root and
    .claude/settings.local.json (enables the server). Use --print-only
    to just display the configuration without creating files.
    """
    import json as json_mod
    from pathlib import Path

    project_path = Path(project_path).resolve()

    from kintsugi.project import _resolve_kintsugi_executable

    mcp_config = {
        "mcpServers": {
            "kintsugi": {
                "command": _resolve_kintsugi_executable(),
                "args": ["mcp", "start"],
                "cwd": str(project_path),
            }
        }
    }

    if print_only:
        console.print(
            Panel.fit(
                json_mod.dumps(mcp_config, indent=2),
                title="Claude Code MCP Configuration (.mcp.json)",
            )
        )
        return

    # Create .mcp.json at project root
    mcp_file = project_path / ".mcp.json"

    if mcp_file.exists():
        with open(mcp_file) as f:
            existing = json_mod.load(f)

        if "kintsugi" in existing.get("mcpServers", {}):
            console.print(f"[yellow]KINTSUGI MCP config already exists in {mcp_file}[/yellow]")
            console.print("Use --print-only to see the configuration.")
            return

        # Add kintsugi to existing .mcp.json
        if "mcpServers" not in existing:
            existing["mcpServers"] = {}
        existing["mcpServers"]["kintsugi"] = mcp_config["mcpServers"]["kintsugi"]
        with open(mcp_file, "w") as f:
            json_mod.dump(existing, f, indent=2)
            f.write("\n")
        console.print(f"[green]Added KINTSUGI MCP config to existing {mcp_file}[/green]")
    else:
        with open(mcp_file, "w") as f:
            json_mod.dump(mcp_config, f, indent=2)
            f.write("\n")
        console.print(f"[green]Created {mcp_file}[/green]")

    # Ensure .claude/settings.local.json enables the MCP server
    claude_dir = project_path / ".claude"
    claude_dir.mkdir(exist_ok=True)
    settings_file = claude_dir / "settings.local.json"

    if settings_file.exists():
        with open(settings_file) as f:
            settings = json_mod.load(f)
        enabled = settings.get("enabledMcpjsonServers", [])
        if "kintsugi" not in enabled:
            enabled.append("kintsugi")
            settings["enabledMcpjsonServers"] = enabled
            with open(settings_file, "w") as f:
                json_mod.dump(settings, f, indent=2)
                f.write("\n")
            console.print(
                f"[green]Added kintsugi to enabledMcpjsonServers in {settings_file}[/green]"
            )
    else:
        settings = {
            "enableAllProjectMcpServers": True,
            "enabledMcpjsonServers": ["kintsugi"],
        }
        with open(settings_file, "w") as f:
            json_mod.dump(settings, f, indent=2)
            f.write("\n")
        console.print(f"[green]Created {settings_file}[/green]")

    console.print(
        Panel.fit(
            json_mod.dumps(mcp_config, indent=2),
            title="Claude Code MCP Configuration (.mcp.json)",
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

    DEPRECATED: This command uses the legacy submit.sh script. For new
    projects, use the Snakemake-based workflow instead:

    \b
        kintsugi workflow config .   # Generate Snakemake config
        kintsugi workflow run .      # Submit via Snakemake

    The Snakemake workflow provides automatic dependency management, failure
    recovery, skip-existing, and aggregate QC reports. This command will be
    removed in KINTSUGI v3.0.0.

    PROJECT_DIR is the path to your KINTSUGI project directory.

    \b
    Examples:
        kintsugi slurm submit .                   # Submit all steps
        kintsugi slurm submit . --steps decon     # Just deconvolution
        kintsugi slurm submit . --cycles 1-3      # Specific cycles
        kintsugi slurm submit . --dry-run         # Preview commands
        kintsugi slurm submit . --use-burst       # Also submit burst jobs
    """
    import warnings
    from pathlib import Path

    warnings.warn(
        "kintsugi slurm submit is deprecated. Use 'kintsugi workflow run' instead. "
        "See 'kintsugi workflow --help' for details. "
        "This command will be removed in KINTSUGI v3.0.0.",
        DeprecationWarning,
        stacklevel=2,
    )
    console.print(
        "[yellow]WARNING: 'kintsugi slurm submit' is deprecated.[/yellow]\n"
        "Use the Snakemake-based workflow instead:\n"
        "  kintsugi workflow config .   # Generate config\n"
        "  kintsugi workflow run .      # Submit jobs\n"
    )

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


def _load_channel_names_lightweight(
    meta_dir: Path,
    channels_per_cycle: int = 4,
) -> dict[int, list[str]]:
    """Load channel names from CHANNELNAMES.txt using only stdlib.

    This is a lightweight reimplementation of ``Kio.load_channel_names`` +
    ``Kio.make_channel_names_unique`` that avoids importing numpy/scipy/skimage.
    Used by ``generate_workflow_config`` which only needs parsed channel names.

    Raises ``FileNotFoundError`` if no channel names file is found.
    """
    import re as _re

    alt_names = [
        "CHANNELNAMES.txt",
        "channelnames.txt",
        "channel_names.txt",
        "channel_names.csv",
        "channels.txt",
        "markers.txt",
    ]
    channel_file = None
    for name in alt_names:
        candidate = Path(meta_dir) / name
        if candidate.exists():
            channel_file = candidate
            break
    if channel_file is None:
        raise FileNotFoundError(f"No channel names file found in {meta_dir}")

    lines = [
        ln.strip()
        for ln in channel_file.read_text().splitlines()
        if ln.strip() and not ln.strip().startswith("#")
    ]
    if not lines:
        raise ValueError(f"Channel names file is empty: {channel_file}")

    channel_dict: dict[int, list[str]] = {}
    first_line = lines[0]

    # Detect cycle-prefixed format
    is_cycle_prefixed = (
        ":" in first_line
        or "\t" in first_line
        or (first_line.split(",")[0].strip().isdigit() and len(first_line.split(",")) > 2)
    )

    if is_cycle_prefixed:
        for line in lines:
            try:
                if ":" in line:
                    cycle_str, names_str = line.split(":", 1)
                    cycle = int(cycle_str.strip())
                    names = [n.strip() for n in names_str.split(",")]
                elif "\t" in line:
                    parts = line.split("\t")
                    cycle = int(parts[0].strip())
                    names = [n.strip() for n in parts[1:]]
                else:
                    parts = line.split(",")
                    cycle = int(parts[0].strip())
                    names = [n.strip() for n in parts[1:]]
                channel_dict[cycle] = names
            except (ValueError, IndexError):
                continue
    else:
        current_cycle = 0
        cycle_channels: list[str] = []
        for line in lines:
            dapi_match = _re.match(r"DAPI[-_]?(\d+)", line, _re.IGNORECASE)
            dapi_bare = _re.match(r"DAPI\s*$", line, _re.IGNORECASE)
            if dapi_match:
                if cycle_channels and current_cycle > 0:
                    channel_dict[current_cycle] = cycle_channels
                current_cycle = int(dapi_match.group(1))
                cycle_channels = [line]
            elif dapi_bare:
                if cycle_channels and current_cycle > 0:
                    channel_dict[current_cycle] = cycle_channels
                current_cycle += 1
                cycle_channels = [line]
            elif current_cycle > 0:
                cycle_channels.append(line)
                if len(cycle_channels) == channels_per_cycle:
                    channel_dict[current_cycle] = cycle_channels
                    cycle_channels = []
        if cycle_channels and current_cycle > 0:
            channel_dict[current_cycle] = cycle_channels

    # --- make unique (inline reimplementation of Kio.make_channel_names_unique) ---
    global_name_counts: dict[str, int] = {}
    for cycle_names in channel_dict.values():
        for name in cycle_names:
            if name.upper().startswith("DAPI"):
                continue
            key = name.upper()
            global_name_counts[key] = global_name_counts.get(key, 0) + 1
    needs_suffix = {n for n, c in global_name_counts.items() if c > 1}

    unique_dict: dict[int, list[str]] = {}
    for cycle_num in sorted(channel_dict.keys()):
        new_names = []
        for ch_idx, name in enumerate(channel_dict[cycle_num]):
            if name.upper().startswith("DAPI"):
                if cycle_num > 1:
                    new_names.append(f"DAPI_cyc{cycle_num:02d}")
                else:
                    new_names.append(name)
            elif name.upper() in needs_suffix:
                letter = chr(ord("a") + ch_idx - 1) if ch_idx > 0 else "a"
                new_names.append(f"{name}_{cycle_num}{letter}")
            else:
                new_names.append(name)
        unique_dict[cycle_num] = new_names

    return unique_dict


def generate_workflow_config(project_dir: Path, print_only: bool = False) -> Path | None:
    """Generate workflow/config.yaml for a project.

    Reads meta/experiment.json (auto-translating CODEX field names) and
    detects SLURM resources to create a unified Snakemake configuration.

    Parameters
    ----------
    project_dir : Path
        Resolved path to the KINTSUGI project directory.
    print_only : bool
        If True, print config to stdout without writing files.

    Returns
    -------
    Path or None
        Path to the created config file, or None if ``print_only``.

    Raises
    ------
    SystemExit
        If experiment.json is missing or invalid.
    """
    import json
    import shutil

    import yaml

    from kintsugi.project import ExperimentConfig

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
    # Use lightweight loader to avoid heavy scipy/skimage imports from Kio.py
    channel_names_dict: dict[int, list[str]] = {}
    try:
        channel_names_dict = _load_channel_names_lightweight(
            project_dir / "meta",
            channels_per_cycle=raw_experiment.get(
                "channels_per_cycle", raw_experiment.get("numChannels", 4)
            ),
        )
        console.print(
            f"  Loaded channel names: {sum(len(v) for v in channel_names_dict.values())} markers"
        )
    except FileNotFoundError:
        console.print(
            "[yellow]  Warning: meta/CHANNELNAMES.txt not found — EDF will use CH# names[/yellow]"
        )
    except ValueError as e:
        console.print(f"[yellow]  Warning: {e}[/yellow]")
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

    # Auto-detect tissue type from project name
    # Use importlib to load batch_multi.py directly, bypassing signal/__init__.py
    # which imports scipy-dependent modules that crash with numpy version mismatches
    try:
        import importlib.util

        _mod_name = "_kintsugi_batch_multi"
        _batch_multi_path = Path(__file__).parent / "signal" / "batch_multi.py"
        _spec = importlib.util.spec_from_file_location(_mod_name, _batch_multi_path)
        if _spec is None or _spec.loader is None:
            raise ImportError(f"Cannot load spec for {_batch_multi_path}")
        _mod = importlib.util.module_from_spec(_spec)
        sys.modules[_mod_name] = _mod  # Required for dataclass processing
        _loaded = False
        try:
            _spec.loader.exec_module(_mod)
            _loaded = True
        finally:
            if not _loaded:
                sys.modules.pop(_mod_name, None)
        tissue_type = _mod.parse_tissue_type(project_dir.name, project_dir)
    except Exception:
        tissue_type = "unknown"

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
        # Post-stitch boundary QC (detects rare misalignment and retries
        # with alternate z-planes for the stitch model)
        "boundary_qc_enabled": True,
        "boundary_qc_threshold": 0.65,
        "boundary_qc_max_alternates": 4,
        # EDF parameters
        "edf": {
            "radius_x": 2,
            "radius_y": 2,
            "sigma": 10.0,
            "tiles": [3, 3],
            "blend_depth": 0,
            "z_smooth_sigma": 1.0,
        },
        # Signal isolation parameters
        "signal_isolation": {
            "method": "auto",
            "tissue_type": tissue_type,
            "tile_smooth_sigma": 0.0,
            "recipe_dir": "",
            "learn": True,
            "force": False,
        },
        # Registration parameters
        "registration": {
            "reference_cycle": 1,
            "reference_channel": 1,
            "max_image_dim": 4096,
            "rigid_only": False,
            "feature_detector": "VggFD",
            "imgs_ordered": True,
            "align_to_reference": True,
            "normalize_dimensions": True,
            "non_rigid": {
                "max_dim": 4096,
                "smoothing_method": "gauss",
                "sigma_ratio": 0.01,
                "n_grid_pts": 50,
                "fold_penalty": 1e-6,
            },
        },
        "quality_gate": {
            "max_zero_pct": 0.10,
            "max_sat_pct": 0.05,
            "min_valid_slices": 3,
        },
        # Spillover correction (REDSEA)
        "redsea": {
            "element_shape": "star",
            "element_size": 0,
            "markers_of_interest": [],
        },
        # 3D Vessel Segmentation
        "vessel3d": {
            "channel": 2,
            "marker": "CD31",
            "sigmas": [1, 2, 4, 8],
            "min_size": 500,
            "prune_min_um": 5.0,
        },
        # Optional pipeline stages (cleanup safety system)
        # Declare which optional stages are planned for this project.
        # Cleanup will BLOCK deletion of intermediate data consumed by
        # enabled optional stages until their sentinels are present.
        # If this section is absent, cleanup conservatively blocks deletion
        # of deconvolved/ data with a warning to update the config.
        "optional_stages": {
            "vessel3d": {
                "enabled": False,  # Set true if 3D vessel segmentation is planned
                "cycles": [],  # Which cycles need it (empty = none planned)
            },
            "spillover": {
                "enabled": False,  # Set true if spillover correction is planned
            },
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
            # Signal isolation resources (CPU-only)
            "mem_signal_isolation": 32000,
            "time_signal_isolation": 120,
            "cpus_signal_isolation": 4,
            "time_qc_signal_isolation": 60,
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
        return None

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

    return config_file


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
    config_file = generate_workflow_config(_resolve_project_dir(project_dir), print_only=print_only)

    if config_file is not None:
        project_path = config_file.parent.parent
        display_path = "." if Path.cwd().resolve() == project_path else str(project_path)
        console.print()
        console.print("[bold]Next steps:[/bold]")
        console.print(f"  kintsugi workflow check {display_path}    # Verify resources")
        console.print(f"  kintsugi workflow run {display_path} -n   # Dry run")
        console.print(f"  kintsugi workflow run {display_path}      # Submit via SLURM")


@workflow.command("check")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def workflow_check(project_dir: str):
    """
    Check SLURM resource availability for multi-account execution.

    Reads workflow/config.yaml, queries all configured accounts via
    sacctmgr, and prints recommended snakemake command with correct -j flag.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """

    import yaml

    from kintsugi.hpc import detect_live_multi_account

    project_dir = _resolve_project_dir(project_dir)
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


def _run_with_dashboard(cmd: list[str], project_dir: Path, interval: int = 30) -> None:
    """Launch snakemake in background and display a live progress dashboard.

    Parameters
    ----------
    cmd : list[str]
        Full snakemake command to execute.
    project_dir : Path
        Resolved project directory for dashboard scanning.
    interval : int
        Dashboard refresh interval in seconds.
    """
    import time

    from kintsugi.dashboard import render_dashboard, scan_project_progress

    proc = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    console.print(f"  Snakemake launched (PID {proc.pid})")
    console.print(f"  Dashboard refreshing every {interval}s. Press Ctrl+C to detach.\n")

    try:
        while True:
            console.clear()
            try:
                status = scan_project_progress(project_dir)
                render_dashboard([status])
            except Exception as e:
                console.print(f"[yellow]Dashboard error: {e}[/yellow]")

            rc = proc.poll()
            if rc is not None:
                # Snakemake finished — render one final dashboard
                console.clear()
                try:
                    status = scan_project_progress(project_dir)
                    render_dashboard([status])
                except Exception:
                    pass
                if rc == 0:
                    console.print("[green]Snakemake completed successfully.[/green]")
                else:
                    stderr_out = proc.stderr.read().decode() if proc.stderr else ""
                    console.print(f"[red]Snakemake exited with code {rc}.[/red]")
                    if stderr_out:
                        console.print(stderr_out[-500:])
                raise SystemExit(rc)

            time.sleep(interval)
    except KeyboardInterrupt:
        console.print(f"\n  Detached from dashboard. Snakemake still running (PID {proc.pid}).")
        console.print("  SLURM jobs continue independently.")
        console.print(f"  Check progress: kintsugi workflow status {project_dir} --watch")


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
@click.option("--dashboard", "-d", is_flag=True, help="Show live progress dashboard while running")
@click.option(
    "--dashboard-interval",
    default=30,
    type=int,
    help="Dashboard refresh interval in seconds (default: 30)",
)
@click.option(
    "--detach",
    is_flag=True,
    help="Run in background (survives SSH disconnect and terminal exit)",
)
def workflow_run(
    project_dir: str,
    jobs: int,
    dry_run: bool,
    local: bool,
    cores: int,
    forcerun: str | None,
    cycles: str | None,
    dashboard: bool,
    dashboard_interval: int,
    detach: bool,
):
    """
    Run the Snakemake processing pipeline.

    Auto-generates workflow/config.yaml from meta/experiment.json if missing.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow run .                    # Full pipeline via SLURM
        kintsugi workflow run . --dry-run          # Preview
        kintsugi workflow run . --local --cores 4  # Local execution
        kintsugi workflow run . --forcerun stitch  # Force re-stitch
        kintsugi workflow run . --dashboard        # Run with live dashboard
        kintsugi workflow run . --detach           # Background (survives disconnect)
    """

    project_dir = _resolve_project_dir(project_dir)
    wf_dir = project_dir / "workflow"

    # Auto-generate config if missing but experiment.json exists
    if not (wf_dir / "config.yaml").exists():
        exp_path = project_dir / "meta" / "experiment.json"
        if exp_path.exists():
            console.print("[yellow]No workflow/config.yaml found. Auto-generating...[/yellow]")
            try:
                generate_workflow_config(project_dir)
            except SystemExit:
                console.print(
                    "[red]Auto-generation failed. Run 'kintsugi workflow config .' manually.[/red]"
                )
                raise SystemExit(1)
        else:
            console.print("[red]No workflow/config.yaml and no meta/experiment.json.[/red]")
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

            if pool["total_gpu_slots"] == 0:
                console.print("[red bold]ERROR: No GPU slots detected.[/red bold]")
                console.print("  Config has no accounts with gpu_slots > 0")
                console.print("  Run: kintsugi workflow config .")
                raise SystemExit(1)

            # Build per-account live availability for the Snakefile to consume.
            # When `clive` is saturated by another user (e.g. QOSGrpMemLimit),
            # gpu_avail/cpu_avail/mem_avail_gb will be 0 and the assignment
            # helpers in Snakefile will route everything to a healthier account.
            live_accounts = [
                {
                    "name": a["name"],
                    "gpu_avail": a.get("gpu_avail", 0),
                    "cpu_avail": a.get("cpu_avail", 0),
                    "mem_avail_gb": max(
                        0,
                        a.get("alloc", {}).get("mem_gb", 0) - a.get("usage", {}).get("mem_gb", 0),
                    ),
                }
                for a in pool["accounts"]
            ]

            # Hard-fail if every account has zero live headroom — better than
            # queuing forever waiting for another user's job to finish.
            if pool["total_gpu_avail"] == 0 and pool.get("total_cpu_avail", 0) == 0:
                console.print(
                    "[red bold]ERROR: All configured SLURM accounts are saturated.[/red bold]"
                )
                for a in pool["accounts"]:
                    used_mem = a.get("usage", {}).get("mem_gb", 0)
                    used_gpu = a.get("usage", {}).get("gpus", 0)
                    console.print(
                        f"  {a['name']}: 0 GPU avail, 0 CPU avail "
                        f"(in use: {used_gpu} GPU, {used_mem} GB mem — "
                        f"may include other users on the QOS)"
                    )
                console.print(
                    "[dim]Wait for other jobs to drain or check `squeue -A <account>`[/dim]"
                )
                raise SystemExit(1)

            if total < 2:
                console.print(
                    f"[yellow]WARNING: Only {total} GPU slot — "
                    f"pipeline stages cannot overlap[/yellow]"
                )

            # Allow QC rules (CPU partition, no GPU) to run alongside GPU jobs.
            # Cap CPU overhead at 4 to avoid overwhelming the scheduler.
            total_with_cpu = total + min(4, pool.get("total_cpu_slots", 0))

            console.print(
                f"  Resources: {pool['total_gpu_avail']} GPU available "
                f"(of {pool['total_gpu_slots']} total GPU slots)"
            )
            for a in live_accounts:
                console.print(
                    f"    {a['name']}: {a['gpu_avail']} GPU avail, "
                    f"{a['cpu_avail']} CPU avail, "
                    f"{a['mem_avail_gb']} GB mem avail"
                )

            cmd.extend(["--profile", str(profile_dir), "-j", str(total_with_cpu)])
            # Forward live availability to Snakemake so the Snakefile's
            # assignment helpers can route around saturated accounts.
            cmd.extend(["--config", f"live_accounts={json.dumps(live_accounts)}"])
        else:
            if not local:
                console.print("[red bold]ERROR: No SLURM profile found.[/red bold]")
                console.print("  Expected: workflow/profiles/slurm/config.yaml")
                console.print("  Run: kintsugi workflow config .")
                raise SystemExit(1)
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

    # ── Detach mode: relaunch in background and return immediately ──
    if detach and not dry_run:
        log_dir = project_dir / "slurm" / "logs" / "snakemake"
        log_dir.mkdir(parents=True, exist_ok=True)
        log_path = log_dir / "workflow_run_detached.log"
        pid_path = project_dir / ".kintsugi_run.pid"

        # Re-invoke ourselves without --detach
        relaunch_cmd = [sys.executable, "-m", "kintsugi.cli", "workflow", "run",
                        str(project_dir)]
        if jobs != 8:
            relaunch_cmd.extend(["-j", str(jobs)])
        if local:
            relaunch_cmd.append("--local")
            relaunch_cmd.extend(["--cores", str(cores)])
        if forcerun:
            relaunch_cmd.extend(["--forcerun", forcerun])
        if cycles:
            relaunch_cmd.extend(["--cycles", cycles])

        with open(log_path, "w") as log_file:
            proc = subprocess.Popen(
                relaunch_cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )

        pid_path.write_text(str(proc.pid))
        console.print(f"  Pipeline launched in background (PID {proc.pid})")
        console.print(f"  Log:    {log_path}")
        console.print(f"  PID:    {pid_path}")
        console.print(f"\n  Monitor: tail -f {log_path}")
        console.print(f"  Status:  kintsugi workflow status {project_dir}")
        console.print(f"  Stop:    kill $(cat {pid_path})")
        return

    # ── SLURM context warning ──
    # Skip if we were relaunched by --detach (start_new_session makes us
    # persistent, even though SLURM_JOB_ID is still inherited).
    _is_session_leader = hasattr(os, "getsid") and os.getsid(0) == os.getpid()
    if (
        not local
        and not dry_run
        and not dashboard
        and not _is_session_leader
        and os.environ.get("SLURM_JOB_ID")
    ):
        console.print(
            "[yellow bold]WARNING: Running inside a SLURM job without "
            "--detach or --dashboard.[/yellow bold]"
        )
        console.print(
            "  The Snakemake coordinator must stay alive to submit "
            "dependent jobs."
        )
        console.print(
            "  If this terminal/process is killed, the pipeline will "
            "stall with no new jobs submitted."
        )
        console.print(
            "  Recommended: add --detach (background) or "
            "--dashboard (Ctrl+C to detach)\n"
        )

    if dashboard and not dry_run:
        _run_with_dashboard(cmd, project_dir, interval=dashboard_interval)
    else:
        try:
            result = subprocess.run(cmd, check=False)
            if result.returncode != 0:
                raise SystemExit(result.returncode)
        except FileNotFoundError:
            console.print("[red]snakemake not found. Install with:[/red]")
            console.print("  pip install snakemake snakemake-executor-plugin-slurm")
            raise SystemExit(1)


# ============================================================================
# Status Dashboard
# ============================================================================


@workflow.command("status")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--watch", "-w", is_flag=True, help="Auto-refresh dashboard")
@click.option(
    "--interval", "-i", default=30, type=int, help="Refresh interval in seconds (default: 30)"
)
@click.option("--no-hardware", is_flag=True, help="Hide hardware utilization section")
@click.option("--no-estimates", is_flag=True, help="Hide completion time estimates")
@click.option("--no-jobs", is_flag=True, help="Hide active SLURM job details")
@click.option("--all-projects", is_flag=True, help="Scan all projects under parent directory")
@click.option("--json", "json_output", is_flag=True, help="Output as JSON instead of table")
def workflow_status(
    project_dir: str,
    watch: bool,
    interval: int,
    no_hardware: bool,
    no_estimates: bool,
    no_jobs: bool,
    all_projects: bool,
    json_output: bool,
):
    """
    Show pipeline progress dashboard with real-time status.

    Displays per-cycle, per-stage completion status with SLURM job details,
    memory/GPU utilization, and estimated completion times.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow status .                # One-shot status
        kintsugi workflow status . --watch        # Auto-refresh dashboard
        kintsugi workflow status . -w -i 15       # Refresh every 15s
        kintsugi workflow status . --json         # JSON output for scripting
        kintsugi workflow status /blue/maigan/smith6jt --all-projects
    """

    from kintsugi.dashboard import (
        render_dashboard,
        scan_project_progress,
        watch_dashboard,
    )

    project_dir = _resolve_project_dir(project_dir)

    if all_projects:
        project_dirs = []
        # Discover projects under this directory
        for config_file in sorted(project_dir.rglob("workflow/config.yaml")):
            project_dirs.append(config_file.parent.parent)
        if not project_dirs:
            console.print("[yellow]No projects with workflow/config.yaml found.[/yellow]")
            raise SystemExit(1)
        console.print(f"Found {len(project_dirs)} projects")
    else:
        if not (project_dir / "workflow" / "config.yaml").exists():
            console.print("[red]No workflow/config.yaml found.[/red]")
            console.print("Run 'kintsugi workflow config .' first.")
            raise SystemExit(1)
        project_dirs = [project_dir]

    if json_output:
        import json as json_mod

        statuses = []
        for pdir in project_dirs:
            try:
                statuses.append(scan_project_progress(pdir))
            except Exception as e:
                console.print(f"[red]Error scanning {pdir}: {e}[/red]")

        # Serialize to JSON
        output = []
        for s in statuses:
            from kintsugi.dashboard import estimate_completion

            est = estimate_completion(s)
            proj = {
                "project_dir": str(s.project_dir),
                "project_name": s.project_name,
                "cycles": s.cycles,
                "cycle_statuses": {
                    cyc: {
                        "stages": cs.stages,
                        "timings": cs.timings,
                    }
                    for cyc, cs in s.cycle_statuses.items()
                },
                "registration_status": s.registration_status,
                "registration_timing": s.registration_timing,
                "qc_statuses": s.qc_statuses,
                "estimate": {
                    "wall_clock_seconds": est["wall_clock_seconds"],
                    "eta": est["eta"].isoformat() if est["eta"] else None,
                    "confidence": est["confidence"],
                    "parallelism": est["parallelism"],
                },
                "scan_time": s.scan_time.isoformat(),
            }
            if s.hardware:
                proj["hardware"] = {
                    "accounts": s.hardware.accounts,
                    "total_gpu_alloc": s.hardware.total_gpu_alloc,
                    "total_gpu_used": s.hardware.total_gpu_used,
                    "total_gpu_avail": s.hardware.total_gpu_avail,
                }
            output.append(proj)

        click.echo(json_mod.dumps(output, indent=2, default=str))
        return

    if watch:
        watch_dashboard(
            project_dirs,
            interval=interval,
            show_hardware=not no_hardware,
            show_estimates=not no_estimates,
            show_jobs=not no_jobs,
        )
    else:
        statuses = []
        for pdir in project_dirs:
            try:
                statuses.append(scan_project_progress(pdir))
            except Exception as e:
                console.print(f"[red]Error scanning {pdir}: {e}[/red]")

        if statuses:
            render_dashboard(
                statuses,
                show_hardware=not no_hardware,
                show_estimates=not no_estimates,
                show_jobs=not no_jobs,
            )


# ============================================================================
# Batch Processing (Multi-Dataset)
# ============================================================================


def _discover_batch_projects(
    projects_dir: Path, force: bool = False, dataset: str | None = None
) -> list[Path]:
    """Find projects eligible for batch processing.

    Eligible = has workflow/config.yaml + data/raw/.staged,
    and missing signal_isolated/.snakemake_complete (unless force).
    """
    eligible = []
    for p in sorted(projects_dir.iterdir()):
        if not p.is_dir():
            continue
        if dataset and dataset.lower() not in p.name.lower():
            continue
        if not (p / "workflow" / "config.yaml").exists():
            continue
        if not (p / "data" / "raw" / ".staged").exists():
            continue
        if (
            not force
            and (p / "data" / "processed" / "signal_isolated" / ".snakemake_complete").exists()
        ):
            continue
        eligible.append(p)
    return eligible


def _detect_project_stage(project_dir: Path) -> str:
    """Detect the furthest completed stage for a project."""
    if (project_dir / "data" / "processed" / "signal_isolated" / ".snakemake_complete").exists():
        return "[green]signal_isolated[/green]"

    if (project_dir / "data" / "processed" / "registered" / ".snakemake_complete").exists():
        return "[cyan]registered[/cyan]"

    # Check EDF - are all cycles done?
    edf_dir = project_dir / "data" / "processed" / "edf"
    if edf_dir.exists() and list(edf_dir.glob("cyc*/.snakemake_complete")):
        return "edf (partial)"

    decon_dir = project_dir / "data" / "processed" / "deconvolved"
    if decon_dir.exists() and list(decon_dir.glob("cyc*/.snakemake_complete")):
        return "deconvolved (partial)"

    stitch_dir = project_dir / "data" / "processed" / "stitched"
    if stitch_dir.exists() and list(stitch_dir.glob("cyc*/.snakemake_complete")):
        return "stitched (partial)"

    if (project_dir / "data" / "raw" / ".staged").exists():
        return "[yellow]staged (not started)[/yellow]"

    return "[dim]not staged[/dim]"


def _run_single_project(
    project_dir: Path, gpu_budget: int, force: bool, phases: list[str] | None = None
) -> tuple[str, int, float, str]:
    """Run kintsugi workflow run for a single project."""
    import time

    cmd = ["kintsugi", "workflow", "run", str(project_dir), "-j", str(gpu_budget)]
    if force:
        cmd.extend(["--forcerun", "all"])

    start = time.time()
    try:
        result = subprocess.run(
            cmd,
            check=False,
            capture_output=True,
            text=True,
            timeout=86400,  # 24h max
        )
        duration = time.time() - start
        return (
            project_dir.name,
            result.returncode,
            duration,
            result.stderr[-1000:] if result.stderr else "",
        )
    except subprocess.TimeoutExpired:
        duration = time.time() - start
        return (project_dir.name, -1, duration, "Timeout after 24h")
    except Exception as e:
        duration = time.time() - start
        return (project_dir.name, -1, duration, str(e))


@workflow.command("batch")
@click.argument("projects_dir", type=click.Path(exists=True))
@click.option("--parallel", "-p", default=1, type=int, help="Concurrent datasets (default: 1)")
@click.option("--dataset", "-d", default=None, help="Process single dataset by name substring")
@click.option("--dry-run", "-n", is_flag=True, help="Preview eligible datasets")
@click.option("--detach", is_flag=True, help="Run in background (survives SSH disconnect)")
@click.option("--force", "-f", is_flag=True, help="Reprocess even if completed")
@click.option("--phases", default=None, help="Comma-separated stages (e.g., stitch,decon,edf)")
def workflow_batch(
    projects_dir: str,
    parallel: int,
    dataset: str | None,
    dry_run: bool,
    detach: bool,
    force: bool,
    phases: str | None,
):
    """Run processing pipeline across multiple datasets.

    Discovers eligible projects (with config + staged data), validates GPU
    resources, and runs kintsugi workflow run for each dataset. This is a
    thin orchestrator -- each dataset runs through the full existing pipeline
    with all GPU routing, SLURM profile, and sentinel tracking intact.

    PROJECTS_DIR is the parent directory containing project subdirectories.

    \b
    Examples:
        kintsugi workflow batch /path/to/projects           # All eligible
        kintsugi workflow batch /path/to/projects --dry-run  # Preview
        kintsugi workflow batch . -d CX_19-004               # Single dataset
        kintsugi workflow batch . -p 2 --detach              # Background, 2 concurrent
        kintsugi workflow batch . --force                    # Reprocess completed
    """
    import time

    import yaml

    from kintsugi.hpc import detect_live_multi_account

    projects_dir_path = Path(projects_dir).resolve()

    # Discover eligible projects
    eligible = _discover_batch_projects(projects_dir_path, force=force, dataset=dataset)

    if not eligible:
        console.print("[yellow]No eligible projects found.[/yellow]")
        console.print("[dim]Projects need workflow/config.yaml + data/raw/.staged")
        if not force:
            console.print("and must not have data/processed/registered/.snakemake_complete")
            console.print("Use --force to reprocess completed projects.[/dim]")
        return

    # Validate GPU resources from first project's config
    config_path = eligible[0] / "workflow" / "config.yaml"
    with open(config_path) as f:
        wf_config = yaml.safe_load(f)

    res = wf_config.get("resources", {})
    account_list = [a["name"] for a in res.get("accounts", [])]

    if not account_list:
        console.print("[red bold]ERROR: No SLURM accounts in config.yaml[/red bold]")
        console.print("  Run: kintsugi workflow config <project_dir>")
        raise SystemExit(1)

    pool = detect_live_multi_account(
        account_list,
        resources_config={
            "cpus_per_cpu_job": res.get("cpu_cpus_per_task", 8),
            "mem_per_cpu_job": res.get("cpu_mem_decon", 48000) // 1000,
            "cpu_utilization_cap": res.get("cpu_utilization_cap", 0.85),
        },
    )

    if pool["total_gpu_slots"] == 0:
        console.print("[red bold]ERROR: No GPU slots detected.[/red bold]")
        console.print("  Config has no accounts with gpu_slots > 0")
        console.print("  Run: kintsugi workflow config <project_dir>")
        raise SystemExit(1)

    gpu_per_proc = max(
        1, (pool.get("total_gpu_avail") or pool["total_gpu_slots"]) // max(1, parallel)
    )

    # Parse phases
    phase_list = [p.strip() for p in phases.split(",")] if phases else None

    # Show discovery summary
    console.print(f"\n[bold]Batch Processing[/bold]  ({len(eligible)} datasets)\n")
    console.print(
        f"  GPU slots: {pool.get('total_gpu_avail', '?')} available "
        f"(of {pool['total_gpu_slots']} total)"
    )
    console.print(f"  Parallel: {parallel}  |  GPU budget per process: {gpu_per_proc}")
    if phase_list:
        console.print(f"  Phases: {', '.join(phase_list)}")
    console.print()

    # Project table
    table = Table(title="Eligible Datasets")
    table.add_column("#", justify="right", style="dim")
    table.add_column("Project", style="cyan")
    table.add_column("Cycles", justify="right")
    table.add_column("Current Stage", style="bold")

    for i, proj in enumerate(eligible, 1):
        # Read cycle count from config
        try:
            with open(proj / "workflow" / "config.yaml") as f:
                pcfg = yaml.safe_load(f)
            n_cycles = str(len(pcfg.get("cycles", [])))
        except Exception:
            n_cycles = "?"

        stage = _detect_project_stage(proj)
        table.add_row(str(i), proj.name, n_cycles, stage)

    console.print(table)
    console.print()

    if dry_run:
        console.print("[dim]Dry run -- no processing started.[/dim]")
        return

    # Detach mode: relaunch in background
    if detach:
        log_path = projects_dir_path / ".kintsugi_batch.log"
        pid_path = projects_dir_path / ".kintsugi_batch.pid"

        # Build command to re-invoke ourselves without --detach
        cmd = [
            sys.executable,
            "-m",
            "kintsugi.cli",
            "workflow",
            "batch",
            str(projects_dir_path),
        ]
        cmd.extend(["-p", str(parallel)])
        if dataset:
            cmd.extend(["-d", dataset])
        if force:
            cmd.append("-f")
        if phases:
            cmd.extend(["--phases", phases])

        with open(log_path, "w") as log_file:
            proc = subprocess.Popen(
                cmd,
                stdout=log_file,
                stderr=subprocess.STDOUT,
                start_new_session=True,
            )

        pid_path.write_text(str(proc.pid))
        console.print(f"  Batch launched in background (PID {proc.pid})")
        console.print(f"  Log: {log_path}")
        console.print(f"  PID: {pid_path}")
        console.print(f"\n  Monitor: tail -f {log_path}")
        console.print(f"  Stop:    kintsugi workflow stop {projects_dir_path}")
        return

    # Execute
    results: list[tuple[str, int, float]] = []
    start_total = time.time()

    if parallel <= 1:
        # Sequential execution
        for i, proj in enumerate(eligible, 1):
            console.print(f"[bold][{i}/{len(eligible)}] Processing {proj.name}...[/bold]")
            name, rc, dur, stderr = _run_single_project(proj, gpu_per_proc, force, phase_list)
            status = "[green]OK[/green]" if rc == 0 else f"[red]FAILED (rc={rc})[/red]"
            console.print(f"  {status}  ({dur / 60:.1f} min)")
            if rc != 0 and stderr:
                console.print(f"  [dim]{stderr[-200:]}[/dim]")
            results.append((name, rc, dur))
    else:
        # Parallel execution
        from concurrent.futures import ThreadPoolExecutor, as_completed

        console.print(f"[bold]Running {len(eligible)} datasets with {parallel} workers[/bold]")
        with ThreadPoolExecutor(max_workers=parallel) as executor:
            futures = {
                executor.submit(_run_single_project, p, gpu_per_proc, force, phase_list): p
                for p in eligible
            }
            for future in as_completed(futures):
                name, rc, dur, stderr = future.result()
                status = "[green]OK[/green]" if rc == 0 else f"[red]FAILED (rc={rc})[/red]"
                console.print(f"  {name}: {status}  ({dur / 60:.1f} min)")
                if rc != 0 and stderr:
                    console.print(f"    [dim]{stderr[-200:]}[/dim]")
                results.append((name, rc, dur))

    # Summary
    total_time = time.time() - start_total
    succeeded = sum(1 for _, rc, _ in results if rc == 0)
    failed = sum(1 for _, rc, _ in results if rc != 0)

    console.print()
    summary = Table(title="Batch Summary")
    summary.add_column("Project", style="cyan")
    summary.add_column("Status")
    summary.add_column("Duration", justify="right")
    for name, rc, dur in results:
        status = "[green]OK[/green]" if rc == 0 else "[red]FAILED[/red]"
        summary.add_row(name, status, f"{dur / 60:.1f} min")
    console.print(summary)
    console.print(f"\n  {succeeded} succeeded, {failed} failed, total {total_time / 60:.1f} min")


@workflow.command("stop")
@click.argument("projects_dir", type=click.Path(exists=True))
def workflow_stop(projects_dir: str):
    """Stop a running batch process.

    Sends SIGTERM to the batch orchestrator. SLURM jobs already submitted
    continue running -- use scancel to cancel them.
    """
    import signal

    pid_path = Path(projects_dir).resolve() / ".kintsugi_batch.pid"
    if not pid_path.exists():
        console.print("[yellow]No batch PID file found.[/yellow]")
        return

    try:
        pid = int(pid_path.read_text().strip())
        os.kill(pid, signal.SIGTERM)
        console.print(f"  Sent SIGTERM to batch process (PID {pid})")
        console.print("  SLURM jobs already submitted continue running.")
        console.print("  Use scancel to cancel individual SLURM jobs.")
        pid_path.unlink()
    except ProcessLookupError:
        console.print(f"  Process {pid} not running. Cleaning up PID file.")
        pid_path.unlink()
    except ValueError:
        console.print(f"[red]Invalid PID in {pid_path}[/red]")


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

    from rich.progress import BarColumn, Progress, TextColumn, TimeRemainingColumn

    from kintsugi.export import prepare_export

    project_dir = _resolve_project_dir(project_dir)

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

    from kintsugi.export import deploy_export

    project_dir = _resolve_project_dir(project_dir)
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

    from kintsugi.export import get_export_status

    project_dir = _resolve_project_dir(project_dir)

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
# Cleanup Commands (Pipeline-Aware Safety System)
# ============================================================================


@workflow.group("cleanup")
def cleanup_group():
    """Pipeline-aware cleanup of intermediate data.

    Checks ALL downstream consumers before deleting intermediate files.
    Stages data to a recoverable trash directory instead of permanent deletion.

    \b
    Commands:
        status   - Show per-directory safe/blocked status with reasons
        plan     - Dry-run showing what would be deleted and sizes
        execute  - Run cleanup with full safety checks
        recover  - List and restore trashed data
        purge    - Permanently delete old trash entries
    """
    pass


@cleanup_group.command("status")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def cleanup_status(project_dir: str):
    """
    Show cleanup safety status for each intermediate directory.

    Checks pipeline dependency graph and reports which directories can
    be safely deleted and which are blocked by incomplete consumers.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """

    from kintsugi.cleanup import assess_cleanup_safety

    project_dir = _resolve_project_dir(project_dir)

    console.print(f"[bold]Cleanup Status[/bold]  ({project_dir.name})\n")

    manifest = assess_cleanup_safety(project_dir)

    if manifest.warnings:
        for warning in manifest.warnings:
            console.print(f"  [yellow]Warning: {warning}[/yellow]")
        console.print()

    if not manifest.entries:
        console.print("[dim]No intermediate directories found.[/dim]")
        return

    table = Table(title="Directory Status")
    table.add_column("Directory", style="cyan")
    table.add_column("Status", style="bold")
    table.add_column("Size", style="dim", justify="right")
    table.add_column("Consumers", style="white")
    table.add_column("Reason", style="white", max_width=60)

    for entry in manifest.entries:
        if entry.status == "safe":
            status_str = "[green]SAFE[/green]"
        elif entry.status == "blocked":
            status_str = "[red]BLOCKED[/red]"
        elif entry.status == "already_deleted":
            status_str = "[dim]GONE[/dim]"
        else:
            status_str = entry.status

        size_str = ""
        if entry.size_bytes > 0:
            if entry.size_gb >= 1:
                size_str = f"{entry.size_gb:.1f} GB"
            else:
                size_str = f"{entry.size_mb:.0f} MB"

        consumers_str = ", ".join(entry.consumers)
        table.add_row(
            entry.directory + "/",
            status_str,
            size_str,
            consumers_str,
            entry.reason,
        )

    console.print(table)

    if manifest.safe_entries:
        total_gb = manifest.total_reclaimable_gb
        console.print(f"\n  [green]Reclaimable: {total_gb:.1f} GB[/green]")
        console.print(f"  Run: kintsugi workflow cleanup execute {project_dir}")

    if manifest.blocked_entries:
        console.print(
            f"\n  [yellow]{len(manifest.blocked_entries)} director(ies) blocked. "
            "See reasons above.[/yellow]"
        )


@cleanup_group.command("plan")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def cleanup_plan(project_dir: str):
    """
    Dry-run: show what would be deleted and sizes.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """

    from kintsugi.cleanup import _dir_size, assess_cleanup_safety

    project_dir = _resolve_project_dir(project_dir)

    console.print(f"[bold]Cleanup Plan[/bold]  ({project_dir.name})\n")

    manifest = assess_cleanup_safety(project_dir)

    if manifest.warnings:
        for warning in manifest.warnings:
            console.print(f"  [yellow]Warning: {warning}[/yellow]")
        console.print()

    # Show what would be deleted
    will_delete = []
    will_skip = []

    for entry in manifest.entries:
        if entry.status == "safe":
            will_delete.append(entry)
        elif entry.status == "blocked":
            will_skip.append(entry)

    # Also check raw/
    raw_dir = project_dir / "data" / "raw"
    if raw_dir.exists():
        from kintsugi.cleanup import _check_stage_sentinels, _get_cycles, _get_stage

        config = {}
        try:
            import yaml

            cfg_path = project_dir / "workflow" / "config.yaml"
            if cfg_path.exists():
                with open(cfg_path) as f:
                    config = yaml.safe_load(f) or {}
        except Exception:
            pass
        cycles = _get_cycles(config, project_dir)
        stitch_stage = _get_stage("stitch")
        if stitch_stage and cycles:
            complete, _ = _check_stage_sentinels(project_dir, stitch_stage, cycles)
            if complete:
                raw_size = _dir_size(raw_dir)
                will_delete.append(
                    type(
                        "Entry",
                        (),
                        {
                            "directory": "raw",
                            "size_bytes": raw_size,
                            "size_gb": raw_size / (1024**3),
                            "reason": "All stitch sentinels present",
                        },
                    )()
                )

    if will_delete:
        console.print("  [bold green]WILL DELETE:[/bold green]")
        total = 0
        for entry in will_delete:
            size = entry.size_bytes
            total += size
            if size >= 1024**3:
                size_str = f"{size / 1024**3:.1f} GB"
            elif size >= 1024**2:
                size_str = f"{size / 1024**2:.0f} MB"
            else:
                size_str = f"{size} bytes"
            console.print(f"    {entry.directory}/  ({size_str})")
        console.print(f"\n    Total: {total / 1024**3:.1f} GB")
    else:
        console.print("  [dim]Nothing to delete.[/dim]")

    if will_skip:
        console.print("\n  [bold yellow]BLOCKED (will skip):[/bold yellow]")
        for entry in will_skip:
            console.print(f"    {entry.directory}/  — {entry.reason}")


@cleanup_group.command("execute")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--no-trash", is_flag=True, help="Permanent delete (skip trash staging)")
@click.option(
    "--force", is_flag=True, help="Skip interactive prompt (safety checks still enforced)"
)
@click.option("--skip-vessel3d", is_flag=True, help="Override vessel3d block (escape hatch)")
def cleanup_execute(project_dir: str, no_trash: bool, force: bool, skip_vessel3d: bool):
    """
    Execute cleanup with full pipeline safety checks.

    By default, moves data to data/.trash/ for recovery. Use --no-trash
    for permanent deletion. Safety checks (sentinel verification) are
    ALWAYS enforced — --force only skips the interactive confirmation.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow cleanup execute .                 # Interactive + trash
        kintsugi workflow cleanup execute . --force         # Skip prompt + trash
        kintsugi workflow cleanup execute . --no-trash      # Permanent delete
        kintsugi workflow cleanup execute . --skip-vessel3d # Override vessel3d block
    """

    from kintsugi.cleanup import assess_cleanup_safety, execute_cleanup

    project_dir = _resolve_project_dir(project_dir)

    console.print(f"[bold]Cleanup Execute[/bold]  ({project_dir.name})\n")

    manifest = assess_cleanup_safety(project_dir)

    if manifest.warnings:
        for warning in manifest.warnings:
            console.print(f"  [yellow]Warning: {warning}[/yellow]")
        console.print()

    safe_count = len(manifest.safe_entries)
    blocked_count = len(manifest.blocked_entries)

    if safe_count == 0 and not skip_vessel3d:
        console.print("[dim]Nothing safe to delete.[/dim]")
        if blocked_count > 0:
            console.print(
                f"  {blocked_count} director(ies) blocked. "
                "Run 'kintsugi workflow cleanup status .' for details."
            )
        return

    # Show summary
    console.print(f"  Safe to delete: {safe_count} directories")
    if blocked_count > 0:
        console.print(f"  Blocked: {blocked_count} directories")
    if manifest.total_reclaimable_gb > 0:
        console.print(f"  Reclaimable: {manifest.total_reclaimable_gb:.1f} GB")

    use_trash = not no_trash
    if use_trash:
        console.print("  Method: move to data/.trash/ (recoverable)")
    else:
        console.print("  Method: [red]permanent delete[/red]")

    if skip_vessel3d:
        console.print("  [yellow]Override: skipping vessel3d safety check[/yellow]")

    # Interactive confirmation
    if not force:
        console.print()
        response = click.prompt(
            "  Proceed with cleanup?",
            type=click.Choice(["y", "n"], case_sensitive=False),
            default="n",
        )
        if response.lower() != "y":
            console.print("[yellow]Cancelled.[/yellow]")
            return

    # Execute
    result = execute_cleanup(
        project_dir,
        manifest=manifest,
        use_trash=use_trash,
        skip_vessel3d=skip_vessel3d,
    )

    # Report results
    if result["deleted"]:
        console.print(f"\n  [green]Deleted {len(result['deleted'])} director(ies):[/green]")
        for d in result["deleted"]:
            size = d.get("size_bytes", 0)
            size_str = f"{size / 1024**3:.1f} GB" if size >= 1024**3 else f"{size / 1024**2:.0f} MB"
            if d.get("trash_path"):
                console.print(f"    {d['directory']}/  ({size_str}) -> .trash/")
            else:
                console.print(f"    {d['directory']}/  ({size_str}) [deleted]")

    if result["skipped"]:
        console.print(f"\n  [yellow]Skipped {len(result['skipped'])} director(ies):[/yellow]")
        for s in result["skipped"]:
            console.print(f"    {s['directory']}/  — {s['reason']}")

    if result["errors"]:
        console.print(f"\n  [red]Errors: {len(result['errors'])}[/red]")
        for e in result["errors"]:
            console.print(f"    {e['directory']}/  — {e['error']}")


@cleanup_group.command("recover")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--entry", default=None, help="Name of trash entry to recover")
def cleanup_recover(project_dir: str, entry: str | None):
    """
    List and restore trashed data.

    Without --entry, lists all items in data/.trash/.
    With --entry, restores the specified item to its original location.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow cleanup recover .                          # List trash
        kintsugi workflow cleanup recover . --entry deconvolved_20260219_120000  # Restore
    """
    from pathlib import Path

    from kintsugi.cleanup import list_trash, recover_trash

    project_dir = _resolve_project_dir(project_dir)

    if entry is None:
        # List trash contents
        entries = list_trash(project_dir)
        if not entries:
            console.print("[dim]Trash is empty.[/dim]")
            return

        console.print(f"[bold]Trash Contents[/bold]  ({project_dir.name})\n")
        table = Table()
        table.add_column("Entry", style="cyan")
        table.add_column("Original Path", style="white")
        table.add_column("Timestamp", style="dim")
        table.add_column("Size", style="dim", justify="right")

        for e in entries:
            name = Path(e["current_path"]).name
            size = e.get("size_bytes", 0)
            size_str = f"{size / 1024**3:.1f} GB" if size >= 1024**3 else f"{size / 1024**2:.0f} MB"
            table.add_row(
                name,
                e.get("original_path", "unknown"),
                e.get("timestamp", "unknown"),
                size_str,
            )

        console.print(table)
        console.print(
            "\n[dim]Recover with: kintsugi workflow cleanup recover . --entry <name>[/dim]"
        )
    else:
        result = recover_trash(project_dir, entry)
        if "error" in result:
            console.print(f"[red]{result['error']}[/red]")
            raise SystemExit(1)
        console.print(f"[green]Recovered: {result['recovered']}[/green]")


@cleanup_group.command("purge")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--days", default=7, type=int, help="Delete entries older than N days (default: 7)")
@click.option("--all", "purge_all", is_flag=True, help="Delete ALL trash entries")
@click.option("--force", is_flag=True, help="Skip confirmation prompt")
def cleanup_purge(project_dir: str, days: int, purge_all: bool, force: bool):
    """
    Permanently delete old trash entries.

    By default, deletes entries older than 7 days. Use --all to purge everything.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow cleanup purge .              # Entries > 7 days
        kintsugi workflow cleanup purge . --days 1     # Entries > 1 day
        kintsugi workflow cleanup purge . --all        # Everything
    """

    from kintsugi.cleanup import list_trash, purge_trash

    project_dir = _resolve_project_dir(project_dir)

    entries = list_trash(project_dir)
    if not entries:
        console.print("[dim]Trash is empty.[/dim]")
        return

    console.print(f"[bold]Purge Trash[/bold]  ({project_dir.name})")
    console.print(f"  Entries in trash: {len(entries)}")
    if purge_all:
        console.print("  [red]Will delete ALL trash entries permanently[/red]")
    else:
        console.print(f"  Will delete entries older than {days} days")

    if not force:
        response = click.prompt(
            "\n  Proceed?",
            type=click.Choice(["y", "n"], case_sensitive=False),
            default="n",
        )
        if response.lower() != "y":
            console.print("[yellow]Cancelled.[/yellow]")
            return

    result = purge_trash(project_dir, max_age_days=days, purge_all=purge_all)

    if result["purged"]:
        total_gb = result["total_bytes"] / (1024**3)
        console.print(
            f"\n  [green]Purged {len(result['purged'])} entries ({total_gb:.1f} GB)[/green]"
        )
    else:
        console.print("\n  [dim]No entries matched the age criteria.[/dim]")


# ============================================================================
# Signal Isolation Commands
# ============================================================================


@workflow.group("isolate")
def isolate_group():
    """Batch signal isolation with per-marker optimization.

    Automatically selects global vs weighted subtraction per channel,
    smooths blanks to remove tile-grid artifacts, and generates QC reports.

    \b
    Commands:
        plan     - Preview per-channel parameters (dry-run)
        run      - Process all channels with optimized parameters
        qc       - Generate QC visualization pages
        status   - Show summary table from manifest
    """
    pass


@isolate_group.command("plan")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option(
    "--method",
    type=click.Choice(["auto", "global", "weighted"]),
    default="auto",
    help="Subtraction method (default: auto per-marker selection)",
)
@click.option("--tissue-type", default=None, help="Tissue type for parameter tuning")
@click.option(
    "--tile-smooth-sigma",
    type=float,
    default=0.0,
    help="Gaussian sigma for blank smoothing (0 to disable, default: 0)",
)
@click.option("--channels", default=None, help="Comma-separated marker names to process")
@click.option(
    "--recipe-dir",
    default=None,
    type=click.Path(exists=True),
    help="Directory with *_param.txt legacy recipe files",
)
def isolate_plan(
    project_dir: str,
    method: str,
    tissue_type: str | None,
    tile_smooth_sigma: float,
    channels: str | None,
    recipe_dir: str | None,
):
    """
    Preview per-channel parameters without processing.

    Analyzes each signal/blank pair and shows the method that would be selected,
    suggested parameters, and blank-to-signal ratio. With --recipe-dir, shows
    legacy recipe parameters instead.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from kintsugi.signal.batch import process_batch

    channel_list = [c.strip() for c in channels.split(",")] if channels else None

    console.print("[bold]Signal Isolation Plan[/bold]  (dry-run)\n")
    if recipe_dir:
        console.print(f"[dim]Using recipes from: {recipe_dir}[/dim]\n")
    result = process_batch(
        project_dir,
        method=method,
        tissue_type=tissue_type,
        tile_smooth_sigma=tile_smooth_sigma,
        dry_run=True,
        channels=channel_list,
        recipe_dir=recipe_dir,
    )

    # Print per-channel summary
    table = Table(title="Per-Channel Parameter Preview")
    table.add_column("Marker", style="cyan")
    table.add_column("Cycle", justify="right")
    table.add_column("Blank")
    table.add_column("Method", style="bold")
    table.add_column("Scale", justify="right")
    table.add_column("B/S Ratio", justify="right", style="dim")
    table.add_column("Correlation", justify="right")

    for marker, ch in sorted(result.channels.items()):
        analysis = ch.analysis or {}
        signal_p99 = analysis.get("signal_p99", 0)
        blank_p99 = analysis.get("blank_p99", 0)
        ratio = f"{blank_p99 / signal_p99:.2f}" if signal_p99 > 0 else "-"
        scale = ch.parameters.get("blank_scale_factor", ch.parameters.get("base_scale_factor", "?"))
        corr = analysis.get("correlation", "?")
        method_str = f"[magenta]{ch.method}[/magenta]" if ch.method == "weighted" else ch.method
        table.add_row(
            marker,
            str(ch.cycle),
            ch.blank_used,
            method_str,
            f"{scale:.2f}" if isinstance(scale, (int, float)) else str(scale),
            ratio,
            f"{corr:.3f}" if isinstance(corr, (int, float)) else str(corr),
        )

    console.print(table)
    recipe_count = result.summary.get("recipe", 0)
    summary_parts = [
        f"\n[bold]Summary:[/bold] {result.summary.get('total', 0)} channels —",
        f"global: {result.summary.get('global', 0)},",
        f"weighted: {result.summary.get('weighted', 0)}",
    ]
    if recipe_count:
        summary_parts.append(f", recipe: {recipe_count}")
    console.print(" ".join(summary_parts))
    if tile_smooth_sigma > 0:
        console.print(f"[dim]Tile smooth sigma: {tile_smooth_sigma}[/dim]")


@isolate_group.command("run")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option(
    "--method",
    type=click.Choice(["auto", "global", "weighted"]),
    default="auto",
    help="Subtraction method (default: auto per-marker selection)",
)
@click.option("--tissue-type", default=None, help="Tissue type for parameter tuning")
@click.option(
    "--tile-smooth-sigma",
    type=float,
    default=0.0,
    help="Gaussian sigma for blank smoothing (0 to disable, default: 0)",
)
@click.option("--channels", default=None, help="Comma-separated marker names to process")
@click.option("--output-dir", default=None, type=click.Path(), help="Output directory")
@click.option("--force", "-f", is_flag=True, help="Overwrite existing outputs")
@click.option(
    "--recipe-dir",
    default=None,
    type=click.Path(exists=True),
    help="Directory with *_param.txt legacy recipe files",
)
@click.option(
    "--learn/--no-learn",
    default=True,
    help="Record outcomes to parameter learning DB (default: --learn)",
)
def isolate_run(
    project_dir: str,
    method: str,
    tissue_type: str | None,
    tile_smooth_sigma: float,
    channels: str | None,
    output_dir: str | None,
    force: bool,
    recipe_dir: str | None,
    learn: bool,
):
    """
    Process all channels with optimized parameters.

    For each signal channel, analyzes the signal/blank pair, selects
    global or weighted subtraction, and saves the result. With --recipe-dir,
    uses legacy Notebook 3 parameter files for multi-step processing.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from kintsugi.signal.batch import process_batch

    channel_list = [c.strip() for c in channels.split(",")] if channels else None

    console.print("[bold]Signal Isolation — Batch Processing[/bold]\n")
    if recipe_dir:
        console.print(f"[dim]Using recipes from: {recipe_dir}[/dim]")
    result = process_batch(
        project_dir,
        method=method,
        tissue_type=tissue_type,
        tile_smooth_sigma=tile_smooth_sigma,
        dry_run=False,
        output_dir=output_dir,
        force=force,
        channels=channel_list,
        recipe_dir=recipe_dir,
        learn=learn,
    )

    # Print summary
    s = result.summary
    console.print(f"\n[bold]Done:[/bold] {s.get('total', 0)} channels processed")
    method_parts = []
    if s.get("recipe", 0):
        method_parts.append(f"recipe: {s['recipe']}")
    if s.get("global", 0):
        method_parts.append(f"global: {s['global']}")
    if s.get("weighted", 0):
        method_parts.append(f"weighted: {s['weighted']}")
    method_parts.append(f"skipped: {s.get('skipped', 0)}")
    method_parts.append(f"errors: {s.get('error', 0)}")
    console.print(f"  {' | '.join(method_parts)}")
    if s.get("mean_quality"):
        console.print(f"  mean quality: {s['mean_quality']:.3f}")

    # Show any warnings
    warned = [(m, ch.warning) for m, ch in result.channels.items() if ch.warning]
    if warned:
        console.print("\n[yellow]Warnings:[/yellow]")
        for marker, warning in warned:
            console.print(f"  {marker}: {warning}")


@isolate_group.command("qc")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--page-size", type=int, default=6, help="Channels per page (default: 6)")
@click.option("--dpi", type=int, default=120, help="Output DPI (default: 120)")
@click.option("--downsample", type=int, default=16, help="Downsample factor (default: 16)")
@click.option(
    "--normalize-mode",
    type=click.Choice(["independent", "matched"]),
    default="independent",
    help="Normalization mode: 'independent' (each image uses own range) or "
    "'matched' (both use Before's range to reduce tile artifact visibility).",
)
def isolate_qc(project_dir: str, page_size: int, dpi: int, downsample: int, normalize_mode: str):
    """
    Generate QC visualization pages for signal isolation results.

    Creates multi-page PNG reports with three columns per channel:
    Before (normalized), After (normalized), Difference (inferno).

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from kintsugi.signal.isolation_qc import generate_qc_pages

    console.print("[bold]Signal Isolation QC[/bold]\n")
    pages = generate_qc_pages(
        project_dir,
        page_size=page_size,
        dpi=dpi,
        downsample=downsample,
        normalize_mode=normalize_mode,
    )

    if pages:
        console.print(f"[green]Generated {len(pages)} QC pages:[/green]")
        for p in pages:
            console.print(f"  {p}")
    else:
        console.print("[yellow]No signal-isolated images found.[/yellow]")


@isolate_group.command("status")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
def isolate_status(project_dir: str):
    """
    Show summary table from signal isolation manifest.

    Displays per-channel method, parameters, quality score, and warnings.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).
    """
    from kintsugi.signal.isolation_qc import generate_summary_table

    output = generate_summary_table(project_dir)
    # generate_summary_table already prints via Rich console, but also returns text
    # Print the text for non-Rich terminals
    if "No manifest found" in output:
        console.print(f"[yellow]{output}[/yellow]")
    else:
        # Re-render with our console for consistent formatting

        project_dir = _resolve_project_dir(project_dir)
        manifest_path = (
            project_dir
            / "data"
            / "processed"
            / "signal_isolated"
            / "signal_isolation_manifest.json"
        )
        if manifest_path.exists():
            import json

            from rich.table import Table as RichTable

            with open(manifest_path) as f:
                manifest = json.load(f)

            table = RichTable(title="Signal Isolation Summary")
            table.add_column("Marker", style="cyan")
            table.add_column("Cycle", justify="right")
            table.add_column("Blank", style="dim")
            table.add_column("Method", style="bold")
            table.add_column("Scale", justify="right")
            table.add_column("B/S Ratio", justify="right", style="dim")
            table.add_column("Q Score", justify="right")
            table.add_column("Zero%", justify="right")
            table.add_column("Flags", style="yellow", max_width=30)

            for marker, info in sorted(manifest.get("channels", {}).items()):
                if info.get("status") == "skipped":
                    continue
                params = info.get("parameters", {})
                analysis = info.get("analysis", {})
                quality = info.get("quality_metrics", {})
                scale = params.get("blank_scale_factor", params.get("base_scale_factor", ""))
                signal_p99 = analysis.get("signal_p99", 0)
                blank_p99 = analysis.get("blank_p99", 0)
                ratio = f"{blank_p99 / signal_p99:.2f}" if signal_p99 > 0 else "-"
                q_score = quality.get("quality_score", "")
                if isinstance(q_score, (int, float)):
                    if q_score >= 0.7:
                        q_str = f"[green]{q_score:.3f}[/green]"
                    elif q_score >= 0.5:
                        q_str = f"[yellow]{q_score:.3f}[/yellow]"
                    else:
                        q_str = f"[red]{q_score:.3f}[/red]"
                else:
                    q_str = str(q_score)
                method = info.get("method", "")
                method_str = f"[magenta]{method}[/magenta]" if method == "weighted" else method
                table.add_row(
                    marker,
                    str(info.get("cycle", "")),
                    info.get("blank_used", ""),
                    method_str,
                    f"{scale:.2f}" if isinstance(scale, (int, float)) else str(scale),
                    ratio,
                    q_str,
                    (
                        f"{info.get('zero_percent', ''):.1f}"
                        if isinstance(info.get("zero_percent"), (int, float))
                        else ""
                    ),
                    info.get("warning", ""),
                )

            summary = manifest.get("summary", {})
            if summary:
                table.add_section()
                table.add_row(
                    f"[bold]Total: {summary.get('total', 0)}[/bold]",
                    "",
                    "",
                    f"G:{summary.get('global', 0)} W:{summary.get('weighted', 0)}",
                    "",
                    "",
                    f"avg={summary.get('mean_quality', 0):.3f}",
                    "",
                    f"skip={summary.get('skipped', 0)} err={summary.get('error', 0)}",
                )

            console.print(table)


@isolate_group.command("batch")
@click.argument("projects_dir", type=click.Path(exists=True), default=".")
@click.option("--dry-run", "-n", is_flag=True, help="Preview without processing")
@click.option("--dataset", "-d", default=None, help="Process single project by name")
@click.option("--force", "-f", is_flag=True, help="Reprocess completed projects")
@click.option(
    "--template-recipe-dir",
    default=None,
    type=click.Path(exists=True),
    help="Override template recipe directory (default: CX_19-001 recipes)",
)
@click.option("--skip-qc", is_flag=True, help="Skip QC page generation")
@click.option(
    "--learn/--no-learn",
    default=True,
    help="Record outcomes to parameter learning DB (default: --learn)",
)
@click.option(
    "--method",
    type=click.Choice(["auto", "global", "weighted"]),
    default="auto",
    help="Override method for auto-analysis channels (default: auto)",
)
@click.option(
    "--min-confidence",
    type=float,
    default=0.6,
    help="Min confidence for learned DB params (default: 0.6)",
)
@click.option(
    "--tile-smooth-sigma",
    type=float,
    default=0.0,
    help="Gaussian sigma for blank smoothing (0 to disable, default: 0)",
)
def isolate_batch(
    projects_dir: str,
    dry_run: bool,
    dataset: str | None,
    force: bool,
    template_recipe_dir: str | None,
    skip_qc: bool,
    learn: bool,
    method: str,
    min_confidence: float,
    tile_smooth_sigma: float,
):
    """
    Run signal isolation across multiple projects.

    Discovers all projects with registered images and processes them
    sequentially. Uses cascading recipe resolution per marker:

    \b
      1. Project's own recipes (Processing_parameters/)
      2. Learned DB params (from previous successful runs)
      3. Template recipes (default: CX_19-001)
      4. Auto-analysis (global vs weighted selection)

    PROJECTS_DIR is the parent directory containing project subdirectories
    (default: current directory).
    """
    from kintsugi.signal.batch_multi import (
        discover_projects,
        process_all_projects,
        resolve_recipes_for_project,
    )

    if dry_run:
        console.print("[bold]Multi-Project Signal Isolation — Preview[/bold]\n")
    else:
        console.print("[bold]Multi-Project Signal Isolation — Batch[/bold]\n")

    # Discovery phase
    from pathlib import Path

    projects = discover_projects(projects_dir, force=force, dataset=dataset)

    if not projects:
        console.print("[yellow]No projects found ready for signal isolation.[/yellow]")
        console.print(
            "[dim]Projects need workflow/config.yaml with channel_names\n"
            "and data/processed/registered/ with .tif files.[/dim]"
        )
        return

    # Show discovery table
    disc_table = Table(title=f"Discovered Projects ({len(projects)})")
    disc_table.add_column("#", justify="right", style="dim")
    disc_table.add_column("Project", style="cyan")
    disc_table.add_column("Tissue", style="bold")
    disc_table.add_column("Channels", justify="right")
    disc_table.add_column("Own Recipes", justify="center")
    disc_table.add_column("Recipe Sources", style="dim")

    for i, proj in enumerate(projects, 1):
        # Resolve recipes for preview
        resolved = resolve_recipes_for_project(
            proj,
            template_recipe_dir=template_recipe_dir,
            min_confidence=min_confidence,
            projects_dir=Path(projects_dir).resolve(),
        )
        source_summary = resolved.summary()
        source_parts = []
        for src, count in sorted(source_summary.items()):
            source_parts.append(f"{src}:{count}")
        disc_table.add_row(
            str(i),
            proj.name,
            proj.tissue_type,
            str(proj.n_registered_channels),
            "[green]Yes[/green]" if proj.has_own_recipes else "-",
            " ".join(source_parts),
        )

    console.print(disc_table)

    if dry_run:
        console.print(f"\n[dim]Dry run — {len(projects)} projects would be processed.[/dim]")
        return

    # Process
    result = process_all_projects(
        projects_dir,
        dry_run=False,
        dataset=dataset,
        force=force,
        template_recipe_dir=template_recipe_dir,
        skip_qc=skip_qc,
        learn=learn,
        method=method,
        min_confidence=min_confidence,
        tile_smooth_sigma=tile_smooth_sigma,
    )

    # Print summary
    console.print("\n[bold]Batch Complete[/bold]")
    console.print(
        f"  Total: {result.total_projects} | "
        f"Completed: {result.completed} | "
        f"Errors: {result.errors} | "
        f"Skipped: {result.skipped}"
    )

    # Per-project results
    if result.project_results:
        res_table = Table(title="Per-Project Results")
        res_table.add_column("Project", style="cyan")
        res_table.add_column("Tissue")
        res_table.add_column("Status", style="bold")
        res_table.add_column("Channels", justify="right")
        res_table.add_column("Quality", justify="right")
        res_table.add_column("Warnings", style="yellow", max_width=40)

        for name, info in result.project_results.items():
            status = info.get("status", "")
            status_str = (
                f"[green]{status}[/green]"
                if status == "completed"
                else f"[red]{status}[/red]"
                if status == "error"
                else status
            )
            summary = info.get("summary", {})
            n_channels = summary.get("total", 0)
            mean_q = summary.get("mean_quality", 0)
            q_str = f"{mean_q:.3f}" if mean_q else "-"
            warns = info.get("warnings", [])
            warn_str = f"{len(warns)} warnings" if warns else ""
            if info.get("error"):
                warn_str = info["error"][:40]

            res_table.add_row(
                name,
                info.get("tissue_type", ""),
                status_str,
                str(n_channels),
                q_str,
                warn_str,
            )

        console.print(res_table)


# ============================================================================
# Channel Clustering Commands
# ============================================================================


@isolate_group.command("cluster")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--n-clusters", type=int, default=None, help="Number of clusters (auto if omitted)")
@click.option(
    "--no-wavelength-aware",
    is_flag=True,
    help="Disable wavelength-aware clustering (allow cross-wavelength clusters)",
)
@click.option("--output", "-o", type=click.Path(), default=None, help="Output JSON path")
@click.option("--plot", "-p", is_flag=True, help="Save PCA cluster plot")
def isolate_cluster(
    project_dir: str,
    n_clusters: int | None,
    no_wavelength_aware: bool,
    output: str | None,
    plot: bool,
):
    """
    Cluster channels by image feature similarity.

    Extracts intensity, texture, and SNR features from registered/EDF images,
    then groups similar channels using agglomerative clustering. Channels in
    the same cluster can share signal isolation parameters.

    \b
    Output JSON includes:
      - cluster assignments (channel → cluster_id)
      - cluster representatives (closest to centroid)
      - feature summary per cluster
    """
    import json
    from pathlib import Path

    from kintsugi.signal.clustering import (
        cluster_channels,
        get_cluster_representatives,
        plot_cluster_summary,
    )
    from kintsugi.signal.features import batch_extract_features

    project_path = _resolve_project_dir(project_dir)

    # Discover channel images
    search_dirs = [
        project_path / "data" / "processed" / "registered",
        project_path / "data" / "processed" / "edf",
        project_path / "data" / "processed" / "stitched",
    ]

    import tifffile

    marker_dict = {}
    source_dir = None
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        for tif in sorted(search_dir.rglob("*.tif")):
            name = tif.stem
            if name not in marker_dict:
                try:
                    marker_dict[name] = tifffile.imread(str(tif))
                except Exception as e:
                    console.print(f"[yellow]Warning: could not load {tif.name}: {e}[/yellow]")
        if marker_dict:
            source_dir = search_dir
            break

    if not marker_dict:
        console.print("[red]No channel images found in project.[/red]")
        return

    console.print(f"Found [bold]{len(marker_dict)}[/bold] channels in {source_dir}")

    # Extract features
    console.print("Extracting channel features...")
    features = batch_extract_features(marker_dict, progress=True)

    # Cluster
    console.print("Clustering channels...")
    assignments = cluster_channels(
        features,
        n_clusters=n_clusters,
        wavelength_aware=not no_wavelength_aware,
    )
    representatives = get_cluster_representatives(features, assignments)

    # Display results
    n_clusters_found = len(representatives)
    console.print(f"\n[bold]{n_clusters_found} clusters[/bold] from {len(assignments)} channels:\n")

    table = Table(title="Channel Clusters")
    table.add_column("Cluster", justify="right", style="bold")
    table.add_column("Representative", style="cyan")
    table.add_column("Members", justify="right")
    table.add_column("Wavelength", style="dim")
    table.add_column("All Members")

    for cid, (rep_name, count) in sorted(representatives.items()):
        members = [n for n, c in assignments.items() if c == cid]
        wl = features[rep_name]["wavelength_group"]
        table.add_row(
            str(cid),
            rep_name,
            str(count),
            wl,
            ", ".join(sorted(members)),
        )

    console.print(table)

    # Save results
    out_data = {
        "n_channels": len(assignments),
        "n_clusters": n_clusters_found,
        "assignments": assignments,
        "representatives": {
            str(k): {"name": v[0], "member_count": v[1]} for k, v in representatives.items()
        },
    }

    if output:
        out_path = Path(output)
    else:
        out_path = project_path / "configs" / "channel_clusters.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)

    with open(out_path, "w") as f:
        json.dump(out_data, f, indent=2)
    console.print(f"\nCluster assignments saved to [bold]{out_path}[/bold]")

    # Plot
    if plot:
        plot_path = out_path.with_suffix(".png")
        plot_cluster_summary(features, assignments, save_path=plot_path)
        console.print(f"Cluster plot saved to [bold]{plot_path}[/bold]")


@isolate_group.command("propagate")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--cluster-id", type=int, required=True, help="Cluster ID to process")
@click.option(
    "--params",
    type=click.Path(exists=True),
    required=True,
    help="JSON file with processing parameters",
)
@click.option("--output-dir", type=click.Path(), default=None, help="Output directory")
def isolate_propagate(
    project_dir: str,
    cluster_id: int,
    params: str,
    output_dir: str | None,
):
    """
    Propagate tuned parameters to all members of a cluster.

    Loads cluster assignments from configs/channel_clusters.json (or the
    output of 'isolate cluster'), applies the given parameters to every
    channel in the specified cluster, and saves processed TIFFs with
    quality scores.
    """
    import json
    from pathlib import Path

    import tifffile

    from kintsugi.signal.clustering import propagate_cluster_parameters

    project_path = _resolve_project_dir(project_dir)

    # Load cluster assignments
    clusters_path = project_path / "configs" / "channel_clusters.json"
    if not clusters_path.exists():
        console.print("[red]No channel_clusters.json found. Run 'isolate cluster' first.[/red]")
        return

    with open(clusters_path) as f:
        cluster_data = json.load(f)

    assignments = cluster_data["assignments"]
    members = [name for name, cid in assignments.items() if cid == cluster_id]

    if not members:
        console.print(f"[red]No channels in cluster {cluster_id}.[/red]")
        return

    console.print(
        f"Propagating to [bold]{len(members)}[/bold] channels in cluster {cluster_id}: "
        f"{', '.join(sorted(members))}"
    )

    # Load parameters
    with open(params) as f:
        param_dict = json.load(f)

    # Load channel images
    search_dirs = [
        project_path / "data" / "processed" / "registered",
        project_path / "data" / "processed" / "edf",
        project_path / "data" / "processed" / "stitched",
    ]

    marker_dict = {}
    for search_dir in search_dirs:
        if not search_dir.exists():
            continue
        for tif in sorted(search_dir.rglob("*.tif")):
            name = tif.stem
            if name in members or name in param_dict.get("blank_map", {}).values():
                if name not in marker_dict:
                    marker_dict[name] = tifffile.imread(str(tif))
        if all(m in marker_dict for m in members):
            break

    # Output directory
    if output_dir:
        out = Path(output_dir)
    else:
        out = project_path / "data" / "processed" / "signal_isolated"

    # Blank map from params or empty
    blank_map = param_dict.get("blank_map", {})

    # Propagate
    results = propagate_cluster_parameters(marker_dict, members, param_dict, blank_map, out)

    # Display results
    table = Table(title=f"Cluster {cluster_id} — Propagation Results")
    table.add_column("Channel", style="cyan")
    table.add_column("Status", style="bold")
    table.add_column("Quality", justify="right")

    for ch, res in sorted(results.items()):
        status = res["status"]
        status_str = f"[green]{status}[/green]" if status == "success" else f"[red]{status}[/red]"
        q = f"{res['quality_score']:.3f}" if res.get("quality_score") else "-"
        table.add_row(ch, status_str, q)

    console.print(table)

    successes = sum(1 for r in results.values() if r["status"] == "success")
    console.print(f"\n{successes}/{len(members)} channels processed successfully → {out}")


# ============================================================================
# Vessel 3D Segmentation Command
# ============================================================================


@workflow.command("vessel3d")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.option("--cycle", type=int, required=True, help="Cycle number with vessel marker")
@click.option("--channel", type=int, default=2, help="Channel index (1-based, default: 2)")
@click.option("--marker", default="CD31", help="Marker name (default: CD31)")
@click.option("--dry-run", "-n", is_flag=True, help="Preview without executing")
def workflow_vessel3d(project_dir: str, cycle: int, channel: int, marker: str, dry_run: bool):
    """
    Run 3D vessel segmentation on a deconvolved z-stack.

    Segments blood/lymphatic vessels using Frangi vesselness filtering,
    morphological cleanup, skeletonization, and graph-based morphometry.

    Requires deconvolution output for the target cycle. Runs as a standalone
    Snakemake rule (not part of the default pipeline).

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow vessel3d . --cycle 3                   # CD31 in cycle 3
        kintsugi workflow vessel3d . --cycle 5 --marker CD34     # CD34 in cycle 5
        kintsugi workflow vessel3d . --cycle 3 --dry-run         # Preview
    """

    project_dir = _resolve_project_dir(project_dir)
    wf_dir = project_dir / "workflow"

    if not (wf_dir / "config.yaml").exists():
        console.print("[red]No workflow/config.yaml found.[/red]")
        console.print("Run 'kintsugi workflow config .' first.")
        raise SystemExit(1)

    # Build snakemake command targeting vessel3d for the specified cycle
    cyc_str = f"{cycle:02d}"
    target = f"{project_dir}/data/processed/vessel_3d/cyc{cyc_str}/.snakemake_complete"

    cmd = [
        "snakemake",
        "--directory",
        str(wf_dir),
        "--snakefile",
        str(wf_dir / "Snakefile"),
        "--configfile",
        str(wf_dir / "config.yaml"),
        "--config",
        f"vessel3d={{channel: {channel}, marker: {marker}}}",
        target,
    ]

    profile_dir = wf_dir / "profiles" / "slurm"
    if profile_dir.exists():
        cmd.extend(["--profile", str(profile_dir), "-j", "1"])
    else:
        cmd.extend(["--cores", "4"])

    if dry_run:
        cmd.append("-n")

    console.print(f"[bold]Vessel 3D Segmentation[/bold]  (cycle {cycle}, {marker})\n")
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
# Spillover Correction Command
# ============================================================================


@workflow.command("spillover")
@click.argument("project_dir", type=click.Path(exists=True), default=".")
@click.argument("sample", required=False, default=None)
@click.option(
    "--element-shape",
    default="star",
    type=click.Choice(["star", "square"]),
    help="Structuring element shape (default: star)",
)
@click.option("--dry-run", "-n", is_flag=True, help="Preview without executing")
def workflow_spillover(project_dir: str, sample: str | None, element_shape: str, dry_run: bool):
    """
    Run REDSEA spillover correction on segmented data.

    Corrects lateral signal spillover between neighboring cells using REDSEA
    (Bai et al., 2021). Requires signal-isolated images and cell segmentation masks.

    SAMPLE is the sample name (matches files in signal_isolated/ and segmented/).
    If omitted, runs for all available samples.

    PROJECT_DIR is the path to your KINTSUGI project directory (default: current).

    \b
    Examples:
        kintsugi workflow spillover .                            # All samples
        kintsugi workflow spillover . sample_01                  # Specific sample
        kintsugi workflow spillover . --element-shape square     # Square kernel
        kintsugi workflow spillover . --dry-run                  # Preview
    """

    project_dir = _resolve_project_dir(project_dir)
    wf_dir = project_dir / "workflow"

    if not (wf_dir / "config.yaml").exists():
        console.print("[red]No workflow/config.yaml found.[/red]")
        console.print("Run 'kintsugi workflow config .' first.")
        raise SystemExit(1)

    cmd = [
        "snakemake",
        "--directory",
        str(wf_dir),
        "--snakefile",
        str(wf_dir / "Snakefile"),
        "--configfile",
        str(wf_dir / "config.yaml"),
        "--config",
        f"redsea={{element_shape: {element_shape}}}",
    ]

    # Target specific sample or all spillover targets
    if sample:
        target = f"{project_dir}/data/processed/spillover_corrected/.{sample}_spillover.complete"
        cmd.append(target)
    else:
        cmd.append("spillover_correction")

    profile_dir = wf_dir / "profiles" / "slurm"
    if profile_dir.exists():
        cmd.extend(["--profile", str(profile_dir), "-j", "1"])
    else:
        cmd.extend(["--cores", "4"])

    if dry_run:
        cmd.append("-n")

    console.print(f"[bold]Spillover Correction[/bold]  ({element_shape} kernel)\n")
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
    "--xy-pixel-size", type=float, default=None, help="XY pixel size in nm (default: 377)"
)
@click.option("--z-step-size", type=float, default=None, help="Z step size in nm (default: 1500)")
@click.option("--numerical-aperture", type=float, default=None, help="Objective NA (default: 0.75)")
@click.option(
    "--tissue-ri", type=float, default=None, help="Tissue refractive index (default: 1.44)"
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
    xy_pixel_size: float | None,
    z_step_size: float | None,
    numerical_aperture: float | None,
    tissue_ri: float | None,
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

            # Update project name if provided
            if name is not None and name != project.config.name:
                old_name = project.config.name
                project.config.name = name
                console.print(f"  Updated project name: '{old_name}' -> '{name}'")

            project.paths.create_all()
            project.setup_notebooks(overwrite=True)
            project._create_claude_config()
            project._create_vscode_config()

            # Regenerate default_parameters.json
            project._create_default_params_file()
            console.print("  Regenerated default_parameters.json")

            # Update experiment.json if any microscope params explicitly provided
            microscope_overrides = {}
            if tile_rows is not None:
                microscope_overrides["tile_rows"] = tile_rows
            if tile_cols is not None:
                microscope_overrides["tile_cols"] = tile_cols
            if xy_pixel_size is not None:
                microscope_overrides["xy_pixel_size"] = xy_pixel_size
            if z_step_size is not None:
                microscope_overrides["z_step_size"] = z_step_size
            if numerical_aperture is not None:
                microscope_overrides["numerical_aperture"] = numerical_aperture
            if tissue_ri is not None:
                microscope_overrides["tissue_refractive_index"] = tissue_ri

            if microscope_overrides:
                exp_config = project.load_experiment_config()
                if exp_config is None:
                    from kintsugi.project import ExperimentConfig

                    exp_config = ExperimentConfig()
                for key, value in microscope_overrides.items():
                    setattr(exp_config, key, value)
                project.save_experiment_config(exp_config)
                console.print(
                    f"  Updated experiment.json: {', '.join(microscope_overrides.keys())}"
                )

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

        # Create the project (apply defaults for None microscope params)
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
            # Microscope parameters (apply defaults for new projects)
            tile_rows=tile_rows,
            tile_cols=tile_cols,
            xy_pixel_size=xy_pixel_size if xy_pixel_size is not None else 377.0,
            z_step_size=z_step_size if z_step_size is not None else 1500.0,
            numerical_aperture=numerical_aperture if numerical_aperture is not None else 0.75,
            tissue_refractive_index=tissue_ri if tissue_ri is not None else 1.44,
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
