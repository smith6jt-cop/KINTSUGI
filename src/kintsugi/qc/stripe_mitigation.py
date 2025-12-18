"""
Stripe artifact mitigation for microscopy images.

Provides multiple strategies for removing or reducing periodic horizontal
stripe artifacts detected by stripe_artifact.py.
"""

from dataclasses import dataclass
from typing import Literal

import numpy as np
from scipy.ndimage import gaussian_filter


@dataclass
class MitigationResult:
    """Result of artifact mitigation."""

    method: Literal["notch", "directional", "interpolate", "keep_as_is"]
    success: bool
    original_score: float
    mitigated_score: float
    improvement: float  # Reduction in stripe score (positive = better)
    image: np.ndarray  # Mitigated image
    message: str = ""


def apply_notch_filter(
    image: np.ndarray,
    frequencies: list[float],
    bandwidth: float = 0.02,
    preserve_dc: bool = True,
) -> np.ndarray:
    """
    Remove specific periodic frequencies using FFT notch filter.

    This is the recommended method for periodic stripe artifacts as it
    precisely targets the artifact frequencies while preserving other
    image content.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    frequencies : list[float]
        Normalized frequencies to remove (0-0.5, from detect_stripe_artifacts)
    bandwidth : float
        Width of the notch in normalized frequency units
    preserve_dc : bool
        If True, preserve the DC (mean) component

    Returns
    -------
    np.ndarray
        Filtered image with same dtype as input
    """
    original_dtype = image.dtype
    height, width = image.shape

    # Compute FFT
    fft = np.fft.fft2(image.astype(np.float64))
    fft_shifted = np.fft.fftshift(fft)

    center_y = height // 2

    # Create notch filter mask
    mask = np.ones_like(fft_shifted, dtype=np.float64)

    for freq in frequencies:
        # Convert normalized frequency to pixel position
        pixel_pos = int(freq * height)

        if pixel_pos <= 0:
            continue

        # Notch width in pixels
        notch_width = max(1, int(bandwidth * height))

        # Create smooth notch (Gaussian roll-off)
        for dy in range(-notch_width * 2, notch_width * 2 + 1):
            y_pos_upper = center_y + pixel_pos + dy
            y_pos_lower = center_y - pixel_pos + dy

            if 0 <= y_pos_upper < height:
                # Gaussian roll-off for smooth transition
                weight = np.exp(-(dy ** 2) / (2 * (notch_width / 2) ** 2))
                attenuation = 1.0 - weight
                mask[y_pos_upper, :] *= attenuation

            if 0 <= y_pos_lower < height:
                weight = np.exp(-(dy ** 2) / (2 * (notch_width / 2) ** 2))
                attenuation = 1.0 - weight
                mask[y_pos_lower, :] *= attenuation

    # Preserve DC if requested
    if preserve_dc:
        mask[center_y, :] = 1.0

    # Apply filter
    filtered_fft = fft_shifted * mask
    filtered = np.fft.ifft2(np.fft.ifftshift(filtered_fft))
    filtered = np.real(filtered)

    # Clip to valid range and restore dtype
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        filtered = np.clip(filtered, info.min, info.max).astype(original_dtype)
    else:
        filtered = filtered.astype(original_dtype)

    return filtered


def apply_directional_filter(
    image: np.ndarray,
    sigma_horizontal: float = 2.0,
    sigma_vertical: float = 0.0,
) -> np.ndarray:
    """
    Apply directional Gaussian blur to reduce stripe artifacts.

    Horizontal stripes can be reduced by blurring along the horizontal
    direction (reducing vertical frequency content). This is a simpler
    but less precise method than notch filtering.

    Parameters
    ----------
    image : np.ndarray
        2D image array
    sigma_horizontal : float
        Blur strength in horizontal direction (perpendicular to stripes)
    sigma_vertical : float
        Blur strength in vertical direction (along stripes, usually 0)

    Returns
    -------
    np.ndarray
        Filtered image with same dtype as input
    """
    original_dtype = image.dtype

    # Convert to float for processing
    img_float = image.astype(np.float64)

    # Apply anisotropic Gaussian filter
    # sigma order is (row, col) = (vertical, horizontal)
    filtered = gaussian_filter(img_float, sigma=(sigma_vertical, sigma_horizontal))

    # Restore dtype
    if np.issubdtype(original_dtype, np.integer):
        info = np.iinfo(original_dtype)
        filtered = np.clip(filtered, info.min, info.max).astype(original_dtype)
    else:
        filtered = filtered.astype(original_dtype)

    return filtered


def interpolate_zplane(
    z_stack: np.ndarray,
    bad_plane: int,
    method: Literal["linear", "nearest"] = "linear",
) -> np.ndarray:
    """
    Replace a bad z-plane with interpolation from adjacent planes.

    This is useful when the artifact is severe and cannot be filtered,
    and adjacent planes are clean.

    Parameters
    ----------
    z_stack : np.ndarray
        3D array of shape (nz, height, width)
    bad_plane : int
        Index of the z-plane to replace
    method : str
        Interpolation method: "linear" or "nearest"

    Returns
    -------
    np.ndarray
        The interpolated z-plane (2D array)

    Raises
    ------
    ValueError
        If bad_plane is at the edge (no neighbors on one side)
    """
    nz = z_stack.shape[0]

    if bad_plane <= 0 or bad_plane >= nz - 1:
        raise ValueError(
            f"Cannot interpolate z-plane {bad_plane}: must have neighbors on both sides. "
            f"Stack has {nz} planes (indices 0 to {nz - 1})."
        )

    plane_below = z_stack[bad_plane - 1].astype(np.float64)
    plane_above = z_stack[bad_plane + 1].astype(np.float64)

    if method == "linear":
        # Simple average of adjacent planes
        interpolated = (plane_below + plane_above) / 2.0
    elif method == "nearest":
        # Use the plane with better quality (assume they're similar, use average)
        interpolated = (plane_below + plane_above) / 2.0
    else:
        raise ValueError(f"Unknown method: {method}")

    # Restore original dtype
    if np.issubdtype(z_stack.dtype, np.integer):
        info = np.iinfo(z_stack.dtype)
        interpolated = np.clip(interpolated, info.min, info.max).astype(z_stack.dtype)
    else:
        interpolated = interpolated.astype(z_stack.dtype)

    return interpolated


def mitigate_artifact(
    image: np.ndarray,
    method: Literal["notch", "directional", "interpolate", "keep_as_is"],
    frequencies: list[float] | None = None,
    z_stack: np.ndarray | None = None,
    z_index: int | None = None,
    notch_bandwidth: float = 0.02,
    directional_sigma: float = 2.0,
) -> MitigationResult:
    """
    Apply the selected mitigation method to an artifact-containing image.

    Parameters
    ----------
    image : np.ndarray
        2D image with artifact
    method : str
        Mitigation method: "notch", "directional", "interpolate", or "keep_as_is"
    frequencies : list[float], optional
        Frequencies to notch (required for "notch" method)
    z_stack : np.ndarray, optional
        Full z-stack (required for "interpolate" method)
    z_index : int, optional
        Index of this plane in z_stack (required for "interpolate" method)
    notch_bandwidth : float
        Bandwidth for notch filter
    directional_sigma : float
        Sigma for directional blur

    Returns
    -------
    MitigationResult
        Result containing mitigated image and metrics
    """
    from .stripe_artifact import detect_stripe_artifact

    # Get original score
    original_result = detect_stripe_artifact(image)
    original_score = original_result.stripe_score

    if method == "keep_as_is":
        return MitigationResult(
            method="keep_as_is",
            success=True,
            original_score=original_score,
            mitigated_score=original_score,
            improvement=0.0,
            image=image,
            message="Image kept as-is (no mitigation applied)",
        )

    elif method == "notch":
        if frequencies is None or len(frequencies) == 0:
            # Use frequencies from detection
            frequencies = original_result.frequency_peaks

        if len(frequencies) == 0:
            return MitigationResult(
                method="notch",
                success=False,
                original_score=original_score,
                mitigated_score=original_score,
                improvement=0.0,
                image=image,
                message="No frequencies to notch (detection may have failed)",
            )

        mitigated = apply_notch_filter(image, frequencies, bandwidth=notch_bandwidth)

    elif method == "directional":
        mitigated = apply_directional_filter(image, sigma_horizontal=directional_sigma)

    elif method == "interpolate":
        if z_stack is None or z_index is None:
            return MitigationResult(
                method="interpolate",
                success=False,
                original_score=original_score,
                mitigated_score=original_score,
                improvement=0.0,
                image=image,
                message="Z-stack and z_index required for interpolation",
            )

        try:
            mitigated = interpolate_zplane(z_stack, z_index)
        except ValueError as e:
            return MitigationResult(
                method="interpolate",
                success=False,
                original_score=original_score,
                mitigated_score=original_score,
                improvement=0.0,
                image=image,
                message=str(e),
            )

    else:
        raise ValueError(f"Unknown method: {method}")

    # Check effectiveness
    mitigated_result = detect_stripe_artifact(mitigated)
    mitigated_score = mitigated_result.stripe_score
    improvement = original_score - mitigated_score

    success = improvement > 0 or mitigated_score < 0.08

    if success:
        message = f"Mitigation successful: score reduced from {original_score:.3f} to {mitigated_score:.3f}"
    else:
        message = f"Mitigation may not be effective: score {original_score:.3f} → {mitigated_score:.3f}"

    return MitigationResult(
        method=method,
        success=success,
        original_score=original_score,
        mitigated_score=mitigated_score,
        improvement=improvement,
        image=mitigated,
        message=message,
    )


def get_recommended_method(
    result,  # StripeArtifactResult
    has_adjacent_planes: bool = True,
) -> Literal["notch", "directional", "interpolate", "keep_as_is"]:
    """
    Recommend a mitigation method based on artifact characteristics.

    Parameters
    ----------
    result : StripeArtifactResult
        Detection result from detect_stripe_artifacts
    has_adjacent_planes : bool
        Whether the z-plane has valid neighbors for interpolation

    Returns
    -------
    str
        Recommended method
    """
    if not result.has_artifact:
        return "keep_as_is"

    # Periodic stripes with clear frequency peaks → notch filter
    if len(result.frequency_peaks) >= 1 and result.severity in ("moderate", "severe"):
        return "notch"

    # Mild artifacts or no clear peaks → directional filter is gentler
    if result.severity == "mild":
        return "directional"

    # Severe artifact without clear peaks and has neighbors → interpolate
    if result.severity == "severe" and has_adjacent_planes:
        return "interpolate"

    # Default to notch for periodic artifacts
    return "notch"
