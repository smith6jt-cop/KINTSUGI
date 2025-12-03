#!/bin/bash
#
# KINTSUGI Installation Script for Linux/macOS
#
# Usage:
#   ./scripts/install.sh [OPTIONS]
#
# Options:
#   --env-name NAME     Environment name (default: KINTSUGI)
#   --no-gpu            Skip GPU/CUDA installation
#   --minimal           Use minimal environment (fewer dependencies)
#   --dev               Include development dependencies
#   --skip-validate     Skip dependency validation after install
#   --help              Show this help message
#

set -e  # Exit on error

# Default values
ENV_NAME="KINTSUGI"
USE_GPU=true
MINIMAL=false
DEV=false
SKIP_VALIDATE=false

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Print functions
print_header() {
    echo -e "${BLUE}============================================${NC}"
    echo -e "${BLUE}$1${NC}"
    echo -e "${BLUE}============================================${NC}"
}

print_success() {
    echo -e "${GREEN}[OK]${NC} $1"
}

print_warning() {
    echo -e "${YELLOW}[WARNING]${NC} $1"
}

print_error() {
    echo -e "${RED}[ERROR]${NC} $1"
}

print_info() {
    echo -e "${BLUE}[INFO]${NC} $1"
}

# Show help
show_help() {
    head -20 "$0" | tail -17 | sed 's/^# //'
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --no-gpu)
            USE_GPU=false
            shift
            ;;
        --minimal)
            MINIMAL=true
            shift
            ;;
        --dev)
            DEV=true
            shift
            ;;
        --skip-validate)
            SKIP_VALIDATE=true
            shift
            ;;
        --help|-h)
            show_help
            ;;
        *)
            print_error "Unknown option: $1"
            show_help
            ;;
    esac
done

# Get script directory
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_DIR="$( cd "$SCRIPT_DIR/.." && pwd )"

print_header "KINTSUGI Installation"
echo ""
print_info "Project directory: $PROJECT_DIR"
print_info "Environment name: $ENV_NAME"
print_info "GPU support: $USE_GPU"
print_info "Minimal install: $MINIMAL"
echo ""

# Detect OS
OS=$(uname -s)
case "$OS" in
    Linux*)
        PLATFORM="linux"
        ENV_FILE="envs/env-linux.yml"
        ;;
    Darwin*)
        PLATFORM="macos"
        ENV_FILE="envs/env-macos.yml"
        ;;
    *)
        print_error "Unsupported operating system: $OS"
        exit 1
        ;;
esac
print_info "Detected platform: $PLATFORM"

# Check for conda
if ! command -v conda &> /dev/null; then
    print_error "Conda not found. Please install Miniconda or Anaconda first."
    print_info "Download from: https://www.anaconda.com/download/success#miniconda"
    exit 1
fi
print_success "Conda found: $(conda --version)"

# Initialize conda for script
eval "$(conda shell.bash hook)"

# Check if environment already exists
if conda env list | grep -q "^$ENV_NAME "; then
    print_warning "Environment '$ENV_NAME' already exists."
    read -p "Do you want to remove and recreate it? (y/N) " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        print_info "Removing existing environment..."
        conda env remove -n "$ENV_NAME" -y
    else
        print_info "Updating existing environment..."
        UPDATE_MODE=true
    fi
fi

# Select environment file
if [ "$MINIMAL" = true ]; then
    ENV_FILE="env_streamlined.yml"
fi

print_header "Creating Conda Environment"
cd "$PROJECT_DIR"

if [ "$UPDATE_MODE" = true ]; then
    print_info "Updating environment from $ENV_FILE..."
    conda env update -n "$ENV_NAME" -f "$ENV_FILE" --prune
else
    print_info "Creating environment from $ENV_FILE..."
    conda env create -n "$ENV_NAME" -f "$ENV_FILE"
fi

print_success "Conda environment created/updated"

# Activate environment
print_header "Installing KINTSUGI Package"
conda activate "$ENV_NAME"
print_success "Environment activated: $ENV_NAME"

# Install KINTSUGI as editable package
print_info "Installing KINTSUGI package..."
if [ "$DEV" = true ]; then
    pip install -e ".[dev]"
else
    pip install -e .
fi
print_success "KINTSUGI package installed"

# Validate installation
if [ "$SKIP_VALIDATE" = false ]; then
    print_header "Validating Installation"

    # Run dependency check
    print_info "Checking dependencies..."
    python -c "from kintsugi import check_dependencies; check_dependencies()" || {
        print_warning "Some optional dependencies may be missing."
        print_info "This is normal for minimal installations."
    }

    # Test basic imports
    print_info "Testing basic imports..."
    python -c "
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
" && print_success "Import test passed" || print_warning "Some imports failed"
fi

print_header "Installation Complete!"
echo ""
print_info "To activate the environment:"
echo "    conda activate $ENV_NAME"
echo ""
print_info "To verify installation:"
echo "    kintsugi check"
echo ""
print_info "To see available commands:"
echo "    kintsugi --help"
echo ""
print_info "For Windows users: Download additional dependencies from Zenodo:"
echo "    https://zenodo.org/records/14969214"
echo ""
