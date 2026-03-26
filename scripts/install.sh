#!/bin/bash
#
# KINTSUGI Installation Script for Linux/macOS
#
# Usage:
#   ./scripts/install.sh [OPTIONS]
#
# Options:
#   --env-name NAME     Environment name (default: KINTSUGI)
#   --features LIST     Comma-separated features to install (gpu,viz,dl,analysis,bio,full)
#   --hpc               Use HPC environment file (envs/env-hpc.yml) with full GPU/CUDA stack
#   --skip-validate     Skip dependency validation after install
#   --help              Show this help message
#

set -e  # Exit on error

# Default values
ENV_NAME="KINTSUGI"
FEATURES=""
SKIP_VALIDATE=false
HPC_MODE=false

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
    cat << EOF
KINTSUGI Installation Script

Usage: ./scripts/install.sh [OPTIONS]

Options:
  --env-name NAME     Environment name (default: KINTSUGI)
  --features LIST     Comma-separated features to install after base
                      Available: gpu, viz, dl, analysis, bio, full
  --hpc               Use HPC environment (includes GPU, CUDA, analysis)
  --skip-validate     Skip dependency validation after install
  --help              Show this help message

Examples:
  ./scripts/install.sh                          # Base install only (desktop Linux)
  ./scripts/install.sh --hpc                    # Full HPC install (HiPerGator)
  ./scripts/install.sh --features gpu           # Base + GPU support
  ./scripts/install.sh --features gpu,viz,dl    # Base + multiple features
  ./scripts/install.sh --features full          # Base + all features
EOF
    exit 0
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --env-name)
            ENV_NAME="$2"
            shift 2
            ;;
        --features)
            FEATURES="$2"
            shift 2
            ;;
        --hpc)
            HPC_MODE=true
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

# Auto-detect HPC if --hpc not explicitly passed
if [ "$HPC_MODE" = false ]; then
    if [ -n "${SLURM_CONF:-}" ] || [ -n "${MODULESHOME:-}" ]; then
        print_info "HPC environment detected (SLURM/modules found)"
        print_info "Use --hpc flag for full HPC install, or proceed with desktop install"
        echo ""
        read -p "Use HPC environment file? (Y/n) " -n 1 -r
        echo
        if [[ ! $REPLY =~ ^[Nn]$ ]]; then
            HPC_MODE=true
        fi
    fi
fi

print_header "KINTSUGI Installation"
echo ""
print_info "Project directory: $PROJECT_DIR"
print_info "Environment name: $ENV_NAME"
if [ "$HPC_MODE" = true ]; then
    print_info "Mode: HPC (comprehensive GPU/CUDA/analysis)"
else
    print_info "Mode: Desktop (base + optional features)"
    print_info "Features: ${FEATURES:-none (base only)}"
fi
echo ""

# Detect OS
OS=$(uname -s)
case "$OS" in
    Linux*)
        PLATFORM="linux"
        if [ "$HPC_MODE" = true ]; then
            ENV_FILE="envs/env-hpc.yml"
        else
            ENV_FILE="envs/env-linux.yml"
        fi
        ;;
    Darwin*)
        PLATFORM="macos"
        ENV_FILE="envs/env-macos.yml"
        if [ "$HPC_MODE" = true ]; then
            print_warning "HPC mode is not supported on macOS. Using standard macOS environment."
            HPC_MODE=false
        fi
        # Check for libvips on macOS
        if ! command -v vips &> /dev/null; then
            print_warning "libvips not found. Installing via Homebrew..."
            if command -v brew &> /dev/null; then
                brew install vips
            else
                print_error "Homebrew not found. Please install libvips manually:"
                print_error "  brew install vips"
                exit 1
            fi
        fi
        ;;
    *)
        print_error "Unsupported operating system: $OS"
        print_error "Use install.ps1 for Windows"
        exit 1
        ;;
esac
print_info "Detected platform: $PLATFORM"
print_info "Environment file: $ENV_FILE"

# Verify env file exists
if [ ! -f "$PROJECT_DIR/$ENV_FILE" ]; then
    print_error "Environment file not found: $PROJECT_DIR/$ENV_FILE"
    exit 1
fi

# Check for conda
if ! command -v conda &> /dev/null; then
    print_error "Conda not found. Please install Miniconda or Miniforge first."
    print_info "Download from: https://github.com/conda-forge/miniforge"
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
print_header "Activating Environment"
conda activate "$ENV_NAME"
print_success "Environment activated: $ENV_NAME"

# HPC-specific post-install steps
if [ "$HPC_MODE" = true ]; then
    print_header "HPC Post-Install Configuration"

    # Deploy conda activation scripts for LD_LIBRARY_PATH
    ACTIVATE_DIR="$CONDA_PREFIX/etc/conda/activate.d"
    DEACTIVATE_DIR="$CONDA_PREFIX/etc/conda/deactivate.d"
    mkdir -p "$ACTIVATE_DIR" "$DEACTIVATE_DIR"

    if [ -f "$PROJECT_DIR/envs/activate.d/env_vars.sh" ]; then
        cp "$PROJECT_DIR/envs/activate.d/env_vars.sh" "$ACTIVATE_DIR/"
        cp "$PROJECT_DIR/envs/deactivate.d/env_vars.sh" "$DEACTIVATE_DIR/"
        chmod +x "$ACTIVATE_DIR/env_vars.sh" "$DEACTIVATE_DIR/env_vars.sh"
        print_success "Deployed LD_LIBRARY_PATH activation scripts"
    else
        print_warning "Activation scripts not found in envs/activate.d/"
    fi

    # Copy CUDA headers for CuPy JIT compilation
    TARGETS_INCLUDE="$CONDA_PREFIX/targets/x86_64-linux/include"
    if [ -d "$TARGETS_INCLUDE" ]; then
        cp -r "$TARGETS_INCLUDE"/* "$CONDA_PREFIX/include/" 2>/dev/null || true
        print_success "CUDA headers copied to conda include directory"
    fi

    # Re-activate to pick up new env vars
    print_info "Re-activating environment to apply LD_LIBRARY_PATH..."
    conda deactivate
    conda activate "$ENV_NAME"
    print_success "Environment re-activated with HPC fixes"

    # Apply SLURM TRES patch
    print_info "Applying SLURM TRES patch..."
    kintsugi patch slurm 2>/dev/null || print_warning "SLURM patch skipped (may not be needed)"
fi

# Install optional features if specified (desktop mode)
if [ -n "$FEATURES" ] && [ "$HPC_MODE" = false ]; then
    print_header "Installing Optional Features"

    # Split features by comma
    IFS=',' read -ra FEATURE_ARRAY <<< "$FEATURES"

    for feature in "${FEATURE_ARRAY[@]}"; do
        feature=$(echo "$feature" | xargs)  # Trim whitespace
        print_info "Installing feature: $feature"
        kintsugi install "$feature" || {
            print_warning "Failed to install $feature, continuing..."
        }
    done
fi

# Validate installation
if [ "$SKIP_VALIDATE" = false ]; then
    print_header "Validating Installation"
    print_info "Checking dependencies..."

    if [ "$HPC_MODE" = true ]; then
        # Strict validation on HPC — fail if required deps are broken
        if ! kintsugi check; then
            print_error "Validation FAILED. Required dependencies are missing or broken."
            print_error "Review the errors above. Common fixes:"
            print_error "  - Run: kintsugi fix-hpc"
            print_error "  - Then: conda deactivate && conda activate $ENV_NAME"
            exit 1
        fi
        print_success "All required dependencies validated"
    else
        kintsugi check || {
            print_warning "Some optional dependencies may be missing."
            print_info "Install them with: kintsugi install <feature>"
        }
    fi
fi

print_header "Installation Complete!"
echo ""
print_info "To activate the environment:"
echo "    conda activate $ENV_NAME"
echo ""
print_info "To verify installation:"
echo "    kintsugi check"
echo ""

if [ "$HPC_MODE" = true ]; then
    print_info "HPC-specific commands:"
    echo "    kintsugi fix-hpc           # Repair HPC env issues"
    echo "    kintsugi patch slurm       # Re-apply SLURM patch after updates"
    echo "    kintsugi workflow run .     # Run Snakemake pipeline"
else
    print_info "To install optional features:"
    echo "    kintsugi install gpu       # GPU acceleration (CuPy for CUDA)"
    echo "    kintsugi install dl        # Deep learning segmentation"
    echo "    kintsugi install analysis  # Spatial analysis (scanpy, scimap)"
    echo "    kintsugi install viz       # Napari visualization"
    echo "    kintsugi install all       # All optional features"
fi
echo ""
