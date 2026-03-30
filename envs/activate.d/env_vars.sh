#!/bin/bash
# KINTSUGI conda activation script
# Ensures conda's libstdc++ takes priority over system libs on HPC
#
# On HPC systems (e.g., HiPerGator), the system /lib64/libstdc++.so.6
# is often too old (missing GLIBCXX_3.4.30, CXXABI_1.3.15), causing
# import failures for scikit-learn, matplotlib, pyvips, stackview, etc.
# This script prepends conda's lib directory so the newer libstdc++
# from conda-forge is loaded first.
#
# Deployed to: $CONDA_PREFIX/etc/conda/activate.d/env_vars.sh

# Save current state for clean deactivation
export _KINTSUGI_OLD_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

# Prepend conda lib so conda's libstdc++.so.6 loads before system's
_CONDA_LIB="${CONDA_PREFIX}/lib"
if [ -d "${_CONDA_LIB}" ]; then
    export LD_LIBRARY_PATH="${_CONDA_LIB}${LD_LIBRARY_PATH:+:${LD_LIBRARY_PATH}}"
fi

# Also add CUDA target libs if present (for CuPy/CUDA runtime)
_CONDA_TARGETS_LIB="${CONDA_PREFIX}/targets/x86_64-linux/lib"
if [ -d "${_CONDA_TARGETS_LIB}" ]; then
    export LD_LIBRARY_PATH="${_CONDA_TARGETS_LIB}:${LD_LIBRARY_PATH}"
fi

# Set CUDA paths for CuPy JIT compilation and torch
if [ -f "${CONDA_PREFIX}/include/cuda.h" ] || [ -d "${CONDA_PREFIX}/targets/x86_64-linux/include" ]; then
    export CUDA_PATH="${CONDA_PREFIX}"
    export CUDA_HOME="${CONDA_PREFIX}"
fi

unset _CONDA_LIB _CONDA_TARGETS_LIB
