"""
Signal Isolation Quality Scoring

Provides a computable composite quality score for signal isolation results,
replacing or supplementing visual inspection. Used by channel clustering
for automated QC flagging and parameter propagation validation.

Components:
    - SNR (sigmoid-normalized)
    - Dynamic range utilization
    - Background uniformity
    - Structure preservation (SSIM)
    - Residual autofluorescence correlation
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger("kintsugi.qc.signal_quality")


def compute_signal_isolation_quality(
    original: np.ndarray,
    processed: np.ndarray,
    blank: np.ndarray | None = None,
    mask: np.ndarray | None = None,
    weights: dict | None = None,
) -> dict:
    """
    Composite quality score for signal isolation results.

    Parameters
    ----------
    original : np.ndarray
        Original (pre-processing) single-channel image.
    processed : np.ndarray
        Processed (post-signal-isolation) image.
    blank : np.ndarray, optional
        Blank/autofluorescence channel used in subtraction.
    mask : np.ndarray, optional
        Binary tissue mask. If None, estimated from Otsu on original.
    weights : dict, optional
        Custom weights for composite score components.
        Keys: 'snr', 'dynamic_range', 'bg_uniformity', 'structure', 'residual_af'
        Default: {'snr': 0.25, 'dynamic_range': 0.15, 'bg_uniformity': 0.20,
                  'structure': 0.25, 'residual_af': 0.15}

    Returns
    -------
    dict with:
        quality_score : float 0-1 (composite, higher is better)
        snr : float (raw SNR value)
        dynamic_range_utilization : float 0-1
        background_uniformity : float 0-1
        structure_preservation : float 0-1
        residual_af_score : float 0-1
        artifact_flag : bool
        components : dict (individual sub-scores before weighting)
        pass_threshold : float (0.4 default)
        passed : bool (quality_score >= pass_threshold)
    """
    default_weights = {
        "snr": 0.25,
        "dynamic_range": 0.15,
        "bg_uniformity": 0.20,
        "structure": 0.25,
        "residual_af": 0.15,
    }
    w = weights or default_weights

    orig_f = original.astype(np.float64)
    proc_f = processed.astype(np.float64)

    if np.issubdtype(original.dtype, np.integer):
        dtype_max = float(np.iinfo(original.dtype).max)
    else:
        dtype_max = float(orig_f.max()) if orig_f.max() > 0 else 1.0

    # Estimate tissue mask if not provided
    if mask is None:
        from skimage.filters import threshold_otsu

        try:
            thresh = threshold_otsu(original)
        except ValueError:
            thresh = np.mean(orig_f)
        mask = orig_f > thresh

    bg_mask = ~mask
    fg_mask = mask

    # 1. SNR (sigmoid-normalized to 0-1, midpoint at SNR=5)
    if bg_mask.any() and fg_mask.any():
        bg_std = float(np.std(proc_f[bg_mask]))
        fg_mean = float(np.mean(proc_f[fg_mask]))
        bg_mean = float(np.mean(proc_f[bg_mask]))
        raw_snr = (fg_mean - bg_mean) / bg_std if bg_std > 0 else 0.0
    else:
        raw_snr = 0.0
    snr_score = float(1.0 / (1.0 + np.exp(-0.5 * (raw_snr - 5.0))))

    # 2. Dynamic range utilization
    p01, p99 = np.percentile(proc_f, [1, 99])
    dr_score = float((p99 - p01) / dtype_max) if dtype_max > 0 else 0.0
    dr_score = min(dr_score, 1.0)

    # 3. Background uniformity
    if bg_mask.any():
        bg_mean_val = float(np.mean(proc_f[bg_mask]))
        bg_cv = float(np.std(proc_f[bg_mask])) / (bg_mean_val + 1e-10)
        bg_uniformity = float(max(0.0, 1.0 - bg_cv))
    else:
        bg_uniformity = 0.5

    # 4. Structure preservation (SSIM on downsampled images)
    try:
        from skimage.metrics import structural_similarity as ssim

        ds = max(1, orig_f.shape[0] // 512)
        orig_ds = orig_f[::ds, ::ds]
        proc_ds = proc_f[::ds, ::ds]
        data_range = max(orig_ds.max(), proc_ds.max()) - min(orig_ds.min(), proc_ds.min())
        if data_range > 0:
            structure_score = float(ssim(orig_ds, proc_ds, data_range=data_range))
        else:
            structure_score = 1.0
    except ImportError:
        structure_score = 0.5

    # 5. Residual autofluorescence (correlation with blank in background)
    if blank is not None and bg_mask.any():
        blank_f = blank.astype(np.float64)
        proc_bg = proc_f[bg_mask].ravel()
        blank_bg = blank_f[bg_mask].ravel()
        # Limit sample size for speed
        n_sample = min(10000, len(proc_bg))
        if proc_bg[:n_sample].std() > 0 and blank_bg[:n_sample].std() > 0:
            corr = float(np.corrcoef(proc_bg[:n_sample], blank_bg[:n_sample])[0, 1])
            residual_af = float(max(0.0, 1.0 - abs(corr)))
        else:
            residual_af = 1.0
    else:
        residual_af = 0.5  # neutral if no blank

    # Artifact detection
    artifact_flag = False
    if np.any(proc_f < 0):
        artifact_flag = True
    if proc_f.size > 0 and np.mean(proc_f == 0) > 0.95:
        artifact_flag = True
    if proc_f.size > 0 and np.mean(proc_f >= dtype_max * 0.99) > 0.01:
        artifact_flag = True

    # Composite score (weighted geometric mean)
    components = {
        "snr": snr_score,
        "dynamic_range": dr_score,
        "bg_uniformity": bg_uniformity,
        "structure": structure_score,
        "residual_af": residual_af,
    }

    log_sum = sum(w[k] * np.log(max(v, 1e-10)) for k, v in components.items())
    composite = float(np.exp(log_sum))

    if artifact_flag:
        composite *= 0.5

    pass_threshold = 0.4

    return {
        "quality_score": composite,
        "snr": raw_snr,
        "dynamic_range_utilization": dr_score,
        "background_uniformity": bg_uniformity,
        "structure_preservation": structure_score,
        "residual_af_score": residual_af,
        "artifact_flag": artifact_flag,
        "components": components,
        "pass_threshold": pass_threshold,
        "passed": composite >= pass_threshold,
    }
