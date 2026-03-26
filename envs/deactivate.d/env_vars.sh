#!/bin/bash
# KINTSUGI conda deactivation script
# Restores environment variables to pre-activation state
#
# Deployed to: $CONDA_PREFIX/etc/conda/deactivate.d/env_vars.sh

# Restore LD_LIBRARY_PATH to pre-activation state
if [ -n "${_KINTSUGI_OLD_LD_LIBRARY_PATH+x}" ]; then
    if [ -z "${_KINTSUGI_OLD_LD_LIBRARY_PATH}" ]; then
        unset LD_LIBRARY_PATH
    else
        export LD_LIBRARY_PATH="${_KINTSUGI_OLD_LD_LIBRARY_PATH}"
    fi
    unset _KINTSUGI_OLD_LD_LIBRARY_PATH
fi

# Unset CUDA paths set during activation
unset CUDA_PATH CUDA_HOME 2>/dev/null
