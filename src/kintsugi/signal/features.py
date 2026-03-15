"""
Channel Feature Extraction for Clustering and Parameter Prediction

Extracts a fixed-length numerical feature vector from a single-channel image.
Used for channel clustering (Tier 1), parameter prediction (Tier 2), and
similarity-based parameter transfer across datasets.

Features include intensity statistics, SNR estimates, texture entropy,
gradient magnitude, blank correlation, and wavelength group inference.
"""

from __future__ import annotations

import logging

import numpy as np

logger = logging.getLogger("kintsugi.signal.features")


# --- Wavelength group inference ---

WAVELENGTH_KEYWORDS: dict[str, list[str]] = {
    "405": ["dapi", "hoechst", "405"],
    "488": ["488", "af488", "fitc", "gfp", "alexa488"],
    "550": ["550", "cy3", "tritc", "af555", "pe", "atto550"],
    "594": ["594", "af594", "texas", "texasred"],
    "650": ["647", "650", "cy5", "af647", "apc"],
    "750": ["750", "cy7", "af750", "ir800"],
}

WAVELENGTH_GROUPS = ["405", "488", "550", "594", "650", "750", "unknown"]


def infer_wavelength_group(channel_name: str) -> str:
    """
    Infer excitation/emission wavelength group from channel name string.

    Matches against known fluorophore naming patterns.
    Returns wavelength string (e.g., '488') or 'unknown'.

    Used to enforce same-wavelength clustering, since autofluorescence
    varies primarily by excitation wavelength.

    Parameters
    ----------
    channel_name : str
        Channel or marker name (e.g., 'CD3-AF488', 'DAPI-01').

    Returns
    -------
    str
        Wavelength group ('405', '488', '550', '594', '650', '750', or 'unknown').
    """
    name_lower = channel_name.lower()
    for group, keywords in WAVELENGTH_KEYWORDS.items():
        if any(kw in name_lower for kw in keywords):
            return group
    return "unknown"


def extract_channel_features(
    image: np.ndarray,
    blank: np.ndarray | None = None,
    channel_name: str | None = None,
    downsample: int = 4,
) -> dict:
    """
    Extract statistical features for channel clustering and parameter prediction.

    Parameters
    ----------
    image : np.ndarray
        2D single-channel image (uint8, uint16, or float).
    blank : np.ndarray, optional
        Corresponding blank/autofluorescence channel.
    channel_name : str, optional
        Channel name for wavelength group inference.
    downsample : int
        Downsample factor for texture/gradient features (speed optimization).
        Statistics (percentiles, moments) are computed on full resolution.

    Returns
    -------
    dict
        Feature dictionary with keys:
        - intensity_percentiles: list[float] — [p1, p5, p10, p25, p50, p75, p90, p95, p99]
        - mean, std, skewness, kurtosis: float
        - snr_estimate: float — estimated from Otsu-separated foreground/background
        - background_level: float — estimated background from Otsu threshold
        - foreground_fraction: float — fraction of pixels above Otsu threshold
        - texture_entropy: float — Shannon entropy of 256-bin intensity histogram
        - gradient_magnitude_mean: float — mean Sobel gradient magnitude
        - blank_correlation: float or None — Pearson r with blank channel
        - blank_ratio_mean: float or None — mean(blank) / mean(signal) ratio
        - wavelength_group: str — inferred from channel_name, or 'unknown'
        - dynamic_range: float — (p99 - p01) as fraction of dtype max
    """
    from scipy import stats as sp_stats
    from scipy.ndimage import sobel

    img = image.astype(np.float64)

    if np.issubdtype(image.dtype, np.integer):
        dtype_max = float(np.iinfo(image.dtype).max)
    else:
        dtype_max = float(img.max()) if img.max() > 0 else 1.0

    # Basic intensity statistics (full resolution)
    percentiles = [1, 5, 10, 25, 50, 75, 90, 95, 99]
    pvals = np.percentile(img, percentiles).tolist()
    mean_val = float(np.mean(img))
    std_val = float(np.std(img))
    skew_val = float(sp_stats.skew(img.ravel()))
    kurt_val = float(sp_stats.kurtosis(img.ravel()))

    # Background estimation via Otsu
    from skimage.filters import threshold_otsu

    try:
        thresh = threshold_otsu(image)
    except ValueError:
        thresh = mean_val

    bg_mask = img <= thresh
    fg_mask = ~bg_mask

    bg_mean = float(np.mean(img[bg_mask])) if bg_mask.any() else 0.0
    bg_std = float(np.std(img[bg_mask])) if bg_mask.any() else 1.0
    fg_mean = float(np.mean(img[fg_mask])) if fg_mask.any() else mean_val
    fg_fraction = float(fg_mask.sum() / img.size)

    snr = (fg_mean - bg_mean) / bg_std if bg_std > 0 else 0.0

    # Texture entropy (256-bin histogram)
    hist, _ = np.histogram(img, bins=256, density=True)
    hist = hist[hist > 0]
    entropy = float(-np.sum(hist * np.log2(hist)))

    # Gradient magnitude (downsampled for speed)
    ds = max(1, downsample)
    img_ds = img[::ds, ::ds]
    gx = sobel(img_ds, axis=0)
    gy = sobel(img_ds, axis=1)
    grad_mean = float(np.mean(np.sqrt(gx**2 + gy**2)))

    # Blank correlation
    blank_corr = None
    blank_ratio = None
    if blank is not None:
        blank_f = blank.astype(np.float64)
        if blank_f.std() > 0 and img.std() > 0:
            s_flat = img[::ds, ::ds].ravel()
            b_flat = blank_f[::ds, ::ds].ravel()
            blank_corr = float(np.corrcoef(s_flat, b_flat)[0, 1])
        blank_mean = float(np.mean(blank_f))
        blank_ratio = blank_mean / mean_val if mean_val > 0 else None

    # Wavelength group
    wl_group = infer_wavelength_group(channel_name) if channel_name else "unknown"

    # Dynamic range
    dynamic_range = (pvals[-1] - pvals[0]) / dtype_max if dtype_max > 0 else 0.0

    return {
        "intensity_percentiles": pvals,
        "mean": mean_val,
        "std": std_val,
        "skewness": skew_val,
        "kurtosis": kurt_val,
        "snr_estimate": snr,
        "background_level": bg_mean,
        "foreground_fraction": fg_fraction,
        "texture_entropy": entropy,
        "gradient_magnitude_mean": grad_mean,
        "blank_correlation": blank_corr,
        "blank_ratio_mean": blank_ratio,
        "wavelength_group": wl_group,
        "dynamic_range": dynamic_range,
    }


def features_to_vector(features: dict) -> np.ndarray:
    """
    Convert feature dict to fixed-length numpy array for clustering/prediction.

    Returns 1D array of length 22. Missing values (None) replaced with 0.
    Wavelength group is one-hot encoded as a 7-element suffix
    (405, 488, 550, 594, 650, 750, unknown).

    Parameters
    ----------
    features : dict
        Feature dictionary from extract_channel_features.

    Returns
    -------
    np.ndarray
        1D float64 array of length 22.
    """
    numeric = [
        features["mean"],
        features["std"],
        features["skewness"],
        features["kurtosis"],
        features["snr_estimate"],
        features["background_level"],
        features["foreground_fraction"],
        features["texture_entropy"],
        features["gradient_magnitude_mean"],
        features.get("blank_correlation") or 0.0,
        features.get("blank_ratio_mean") or 0.0,
        features["dynamic_range"],
    ]
    # Percentile summary (reduce 9 to 3 for compactness)
    p = features["intensity_percentiles"]
    numeric.extend([p[0], p[4], p[8]])  # p1, p50, p99

    # One-hot wavelength group
    wl_one_hot = [1.0 if features["wavelength_group"] == g else 0.0 for g in WAVELENGTH_GROUPS]

    return np.array(numeric + wl_one_hot, dtype=np.float64)


def batch_extract_features(
    marker_dict: dict,
    blank_map: dict | None = None,
    progress: bool = True,
) -> dict:
    """
    Extract features for all channels in marker_dict.

    Parameters
    ----------
    marker_dict : dict
        Channel name -> image array (numpy or dask).
    blank_map : dict, optional
        Channel name -> blank channel name mapping.
    progress : bool
        Show tqdm progress bar if available.

    Returns
    -------
    dict
        channel_name -> feature dict
    """
    try:
        from tqdm import tqdm

        iterator = (
            tqdm(marker_dict.items(), desc="Extracting features")
            if progress
            else marker_dict.items()
        )
    except ImportError:
        iterator = marker_dict.items()

    features = {}
    for name, data in iterator:
        img = data.compute() if hasattr(data, "compute") else np.array(data)

        blank = None
        if blank_map and name in blank_map:
            blank_name = blank_map[name]
            if blank_name in marker_dict:
                b = marker_dict[blank_name]
                blank = b.compute() if hasattr(b, "compute") else np.array(b)

        features[name] = extract_channel_features(img, blank=blank, channel_name=name)
    return features
