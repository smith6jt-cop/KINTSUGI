import warnings
try:
    import cupy as cp
    _HAS_CUPY = True
except Exception:
    cp = None
    _HAS_CUPY = False

import numpy as np

try:
    import KCorrect
except Exception:
    KCorrect = None


def to_device(x):
    if _HAS_CUPY:
        return cp.asarray(x)
    return x


def to_host(x):
    if _HAS_CUPY and isinstance(x, cp.ndarray):
        return cp.asnumpy(x)
    return x


def correct_stack(arr, flat, dark, clip_min=0.0, clip_max=1.0, out_dtype=None):
    """
    Elementwise correction: (arr - dark) / flat with clipping.
    Accepts NumPy or CuPy arrays; returns same-device array.
    arr, flat, dark should be broadcastable; arr in [0,1].
    """
    xp = cp if (_HAS_CUPY and isinstance(arr, cp.ndarray)) else np
    eps_mask = (flat != 0)
    arr = arr - dark
    arr = xp.where(eps_mask, arr / flat, 0)
    arr = xp.clip(arr, clip_min, clip_max)
    if out_dtype is not None:
        arr = arr.astype(out_dtype, copy=False)
    return arr


def compute_fields(arr,
                   if_darkfield=True,
                   max_iterations=500,
                   optimization_tolerance=1e-6,
                   max_reweight_iterations=25,
                   reweight_tolerance=1e-3):
    """
    Compute (flat, dark) with KCorrect on CPU (NumPy arrays).
    `arr` may be CuPy or NumPy; converted to NumPy for KCorrect.
    """
    if KCorrect is None:
        raise ImportError("KCorrect is required to compute illumination fields")
    if _HAS_CUPY and isinstance(arr, cp.ndarray):
        arr_np = cp.asnumpy(arr)
    else:
        arr_np = arr
    flat, dark = KCorrect.KCorrect(
        arr_np,
        if_darkfield=if_darkfield,
        max_iterations=max_iterations,
        optimization_tolerance=optimization_tolerance,
        max_reweight_iterations=max_reweight_iterations,
        reweight_tolerance=reweight_tolerance,
    )
    return flat, dark
