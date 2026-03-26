#!/bin/bash
#
# KINTSUGI HPC Environment Verification Script
#
# Checks that the conda environment is correctly configured for HPC use.
# Run after installation or after running 'kintsugi fix-hpc'.
#
# Usage:
#   conda activate KINTSUGI
#   ./scripts/verify_hpc_env.sh
#

set -e

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

PASS=0
FAIL=0
WARN=0

check_pass() {
    echo -e "  ${GREEN}[OK]${NC} $1"
    PASS=$((PASS + 1))
}

check_fail() {
    echo -e "  ${RED}[FAIL]${NC} $1"
    FAIL=$((FAIL + 1))
}

check_warn() {
    echo -e "  ${YELLOW}[WARN]${NC} $1"
    WARN=$((WARN + 1))
}

echo "============================================"
echo "KINTSUGI HPC Environment Verification"
echo "============================================"
echo ""

# 0. Check conda environment is active
echo "[Conda Environment]"
if [ -z "${CONDA_PREFIX:-}" ]; then
    check_fail "No conda environment active"
    echo "  Run: conda activate KINTSUGI"
    exit 1
fi
check_pass "Conda env active: $(basename "$CONDA_PREFIX")"

# 1. Check LD_LIBRARY_PATH includes conda lib
echo ""
echo "[Library Paths]"
if echo "${LD_LIBRARY_PATH:-}" | grep -q "$CONDA_PREFIX/lib"; then
    check_pass "LD_LIBRARY_PATH includes \$CONDA_PREFIX/lib"
else
    check_fail "LD_LIBRARY_PATH missing \$CONDA_PREFIX/lib"
    echo "         Fix: kintsugi fix-hpc && conda deactivate && conda activate KINTSUGI"
fi

# 2. Check conda's libstdc++ exists
if [ -f "$CONDA_PREFIX/lib/libstdc++.so.6" ]; then
    check_pass "Conda libstdc++.so.6 exists"
else
    check_fail "Conda libstdc++.so.6 not found"
    echo "         Fix: conda install 'libstdcxx-ng>=12.0' -c conda-forge"
fi

# 3. Check activation scripts are deployed
if [ -f "$CONDA_PREFIX/etc/conda/activate.d/env_vars.sh" ]; then
    check_pass "Activation scripts deployed"
else
    check_fail "Activation scripts not deployed"
    echo "         Fix: kintsugi fix-hpc"
fi

# 4. Critical imports
echo ""
echo "[Core Packages]"

for pkg in sklearn matplotlib pyvips stackview scipy pandas; do
    if python -c "import $pkg" 2>/dev/null; then
        check_pass "$pkg imports successfully"
    else
        check_fail "$pkg import FAILED (likely libstdc++ issue)"
    fi
done

# 5. GPU packages
echo ""
echo "[GPU / CUDA]"

if python -c "import torch; assert torch.version.cuda is not None; print(f'CUDA {torch.version.cuda}')" 2>/dev/null; then
    CUDA_VER=$(python -c "import torch; print(torch.version.cuda)" 2>/dev/null)
    check_pass "PyTorch with CUDA $CUDA_VER"
else
    if python -c "import torch" 2>/dev/null; then
        check_fail "PyTorch is CPU-only (no CUDA)"
    else
        check_fail "PyTorch not installed"
    fi
fi

if python -c "import cupy" 2>/dev/null; then
    check_pass "CuPy importable"
else
    check_fail "CuPy import failed"
fi

# Note: torch.cuda.is_available() is expected to be False on login nodes
if python -c "import torch; assert torch.cuda.is_available()" 2>/dev/null; then
    check_pass "CUDA hardware available"
else
    check_warn "CUDA hardware not available (expected on login nodes)"
fi

# 6. Analysis packages
echo ""
echo "[Analysis Packages]"

for pkg in scanpy anndata phenograph scimap; do
    if python -c "import $pkg" 2>/dev/null; then
        check_pass "$pkg"
    else
        check_warn "$pkg not available (optional)"
    fi
done

# 7. numpy version
echo ""
echo "[Constraints]"

NUMPY_VER=$(python -c "import numpy; print(numpy.__version__)" 2>/dev/null)
if python -c "from packaging.version import parse; import numpy; assert parse(numpy.__version__) < parse('2.0')" 2>/dev/null; then
    check_pass "numpy $NUMPY_VER (< 2.0)"
else
    check_fail "numpy $NUMPY_VER (>= 2.0 will break processing!)"
fi

# 8. Full kintsugi check
echo ""
echo "[KINTSUGI Check]"
if kintsugi check > /dev/null 2>&1; then
    check_pass "kintsugi check passed"
else
    check_warn "kintsugi check reported issues (run 'kintsugi check' for details)"
fi

# Summary
echo ""
echo "============================================"
echo "SUMMARY"
echo "============================================"
echo -e "  ${GREEN}Passed: $PASS${NC}"
if [ $WARN -gt 0 ]; then
    echo -e "  ${YELLOW}Warnings: $WARN${NC}"
fi
if [ $FAIL -gt 0 ]; then
    echo -e "  ${RED}Failed: $FAIL${NC}"
    echo ""
    echo "To fix most failures, run:"
    echo "  kintsugi fix-hpc"
    echo "  conda deactivate && conda activate KINTSUGI"
    exit 1
else
    echo ""
    echo -e "  ${GREEN}Environment is ready for HPC use!${NC}"
fi
