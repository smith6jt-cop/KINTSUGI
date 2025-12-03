#
# KINTSUGI Installation Script for Windows
#
# Usage:
#   .\scripts\install.ps1 [OPTIONS]
#
# Options:
#   -EnvName NAME       Environment name (default: KINTSUGI)
#   -NoGpu              Skip GPU/CUDA installation
#   -Minimal            Use minimal environment (fewer dependencies)
#   -Dev                Include development dependencies
#   -SkipValidate       Skip dependency validation after install
#   -Help               Show this help message
#

param(
    [string]$EnvName = "KINTSUGI",
    [switch]$NoGpu,
    [switch]$Minimal,
    [switch]$Dev,
    [switch]$SkipValidate,
    [switch]$Help
)

# Enable strict mode
Set-StrictMode -Version Latest
$ErrorActionPreference = "Stop"

# Colors
function Write-Header($message) {
    Write-Host "============================================" -ForegroundColor Blue
    Write-Host $message -ForegroundColor Blue
    Write-Host "============================================" -ForegroundColor Blue
}

function Write-Success($message) {
    Write-Host "[OK] " -ForegroundColor Green -NoNewline
    Write-Host $message
}

function Write-Warning($message) {
    Write-Host "[WARNING] " -ForegroundColor Yellow -NoNewline
    Write-Host $message
}

function Write-Error($message) {
    Write-Host "[ERROR] " -ForegroundColor Red -NoNewline
    Write-Host $message
}

function Write-Info($message) {
    Write-Host "[INFO] " -ForegroundColor Blue -NoNewline
    Write-Host $message
}

# Show help
if ($Help) {
    Get-Content $PSCommandPath | Select-Object -First 17 | Select-Object -Skip 1 | ForEach-Object { $_ -replace '^# ', '' }
    exit 0
}

# Get script directory
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Path
$ProjectDir = Split-Path -Parent $ScriptDir

Write-Header "KINTSUGI Installation"
Write-Host ""
Write-Info "Project directory: $ProjectDir"
Write-Info "Environment name: $EnvName"
Write-Info "GPU support: $(-not $NoGpu)"
Write-Info "Minimal install: $Minimal"
Write-Host ""

# Check for conda
try {
    $condaVersion = conda --version
    Write-Success "Conda found: $condaVersion"
} catch {
    Write-Error "Conda not found. Please install Miniconda or Anaconda first."
    Write-Info "Download from: https://www.anaconda.com/download/success#miniconda"
    exit 1
}

# Initialize conda for PowerShell
$condaPath = (Get-Command conda).Source | Split-Path -Parent | Split-Path -Parent
. "$condaPath\shell\condabin\conda-hook.ps1"

# Check if environment exists
$envList = conda env list
$UpdateMode = $false
if ($envList -match "^$EnvName\s") {
    Write-Warning "Environment '$EnvName' already exists."
    $response = Read-Host "Do you want to remove and recreate it? (y/N)"
    if ($response -eq 'y' -or $response -eq 'Y') {
        Write-Info "Removing existing environment..."
        conda env remove -n $EnvName -y
    } else {
        Write-Info "Updating existing environment..."
        $UpdateMode = $true
    }
}

# Select environment file
$EnvFile = "envs\env-windows.yml"
if ($Minimal) {
    $EnvFile = "env_streamlined.yml"
}

Write-Header "Creating Conda Environment"
Set-Location $ProjectDir

if ($UpdateMode) {
    Write-Info "Updating environment from $EnvFile..."
    conda env update -n $EnvName -f $EnvFile --prune
} else {
    Write-Info "Creating environment from $EnvFile..."
    conda env create -n $EnvName -f $EnvFile
}

Write-Success "Conda environment created/updated"

# Activate environment
Write-Header "Installing KINTSUGI Package"
conda activate $EnvName
Write-Success "Environment activated: $EnvName"

# Install KINTSUGI as editable package
Write-Info "Installing KINTSUGI package..."
if ($Dev) {
    pip install -e ".[dev]"
} else {
    pip install -e .
}
Write-Success "KINTSUGI package installed"

# Validate installation
if (-not $SkipValidate) {
    Write-Header "Validating Installation"

    # Run dependency check
    Write-Info "Checking dependencies..."
    try {
        python -c "from kintsugi import check_dependencies; check_dependencies()"
    } catch {
        Write-Warning "Some optional dependencies may be missing."
        Write-Info "This is normal for minimal installations."
    }

    # Test basic imports
    Write-Info "Testing basic imports..."
    $testScript = @"
import kintsugi
print(f'KINTSUGI version: {kintsugi.__version__}')

# Test pyvips
try:
    import pyvips
    v = pyvips.version(0)
    print(f'libvips version: {v}.{pyvips.version(1)}.{pyvips.version(2)}')
except Exception as e:
    print(f'libvips: Not available ({e})')

# Test numpy
import numpy as np
print(f'NumPy version: {np.__version__}')

print('Basic imports successful!')
"@
    try {
        python -c $testScript
        Write-Success "Import test passed"
    } catch {
        Write-Warning "Some imports failed"
    }
}

Write-Header "Installation Complete!"
Write-Host ""
Write-Info "To activate the environment:"
Write-Host "    conda activate $EnvName"
Write-Host ""
Write-Info "To verify installation:"
Write-Host "    kintsugi check"
Write-Host ""
Write-Info "To see available commands:"
Write-Host "    kintsugi --help"
Write-Host ""
Write-Info "IMPORTANT: Download additional dependencies from Zenodo:"
Write-Host "    https://zenodo.org/records/14969214"
Write-Host ""
Write-Info "Extract the following to the KINTSUGI folder:"
Write-Host "    - maven-3.9.9"
Write-Host "    - java-jdk21"
Write-Host "    - PyVips-dev-8.16"
Write-Host "    - FIJI with Clij2 plugin"
Write-Host "    - MatlabRuntime 2024a (optional, for deconvolution)"
Write-Host ""
