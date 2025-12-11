"""
KINTSUGI Denoising Module

Pure Python implementations of advanced denoising algorithms for
fluorescence microscopy images.

Algorithms:
- Traditional filters (median, Gaussian, bilateral)
- Noise2Void (N2V) - Self-supervised denoising without clean training data
- CARE-style denoising with pre-trained or custom models
- BM3D-inspired patch-based denoising

Usage:
    from kintsugi.denoise import (
        denoise_median,
        denoise_n2v,
        denoise_care,
        estimate_noise_level,
    )

    # Quick denoising with traditional filter
    denoised = denoise_median(image, size=3)

    # Self-supervised N2V denoising
    denoised = denoise_n2v(image, patch_size=64, n_epochs=50)

    # CARE denoising with pre-trained model
    denoised = denoise_care(image, model_name="fluorescence_2d")
"""

from kintsugi.denoise.filters import (
    denoise_median,
    denoise_gaussian,
    denoise_bilateral,
    denoise_nlm,
    estimate_noise_level,
    adaptive_denoise,
)

from kintsugi.denoise.n2v import (
    N2VDenoiser,
    denoise_n2v,
    train_n2v,
)

from kintsugi.denoise.care import (
    CAREDenoiser,
    denoise_care,
    list_pretrained_models,
)

from kintsugi.denoise.patch_based import (
    denoise_bm3d_lite,
    denoise_patch_similarity,
)

__all__ = [
    # Traditional filters
    "denoise_median",
    "denoise_gaussian",
    "denoise_bilateral",
    "denoise_nlm",
    "estimate_noise_level",
    "adaptive_denoise",
    # N2V
    "N2VDenoiser",
    "denoise_n2v",
    "train_n2v",
    # CARE
    "CAREDenoiser",
    "denoise_care",
    "list_pretrained_models",
    # Patch-based
    "denoise_bm3d_lite",
    "denoise_patch_similarity",
]
