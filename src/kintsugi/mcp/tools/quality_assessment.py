"""
Quality Assessment Tools for KINTSUGI MCP Server

Wraps dl_refinement functionality for channel quality assessment.
"""

from __future__ import annotations

import logging
from typing import Any, Dict, List, Optional

import numpy as np

logger = logging.getLogger("kintsugi.mcp.quality_assessment")

# Import session state from signal_isolation
from kintsugi.mcp.tools.signal_isolation import _loaded_images, _add_history


def _get_image(name: str) -> Any:
    """Get a loaded image by name."""
    if name not in _loaded_images:
        raise ValueError(f"Image not loaded: {name}. Use load_channel first.")
    return _loaded_images[name]["data"]


async def assess_quality(
    channel: str,
    use_dl_model: bool = False,
    model_path: Optional[str] = None,
    confidence_threshold: float = 0.85,
    tile_size: int = 512,
) -> Dict[str, Any]:
    """
    Assess channel quality using heuristic and/or deep learning metrics.

    Returns comprehensive quality metrics including SNR, contrast, and
    whether the channel requires manual review.
    """
    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute data if dask array
    if hasattr(data, "compute"):
        # Sample for efficiency
        step = max(1, data.shape[0] // 2000)
        data_sample = data[::step, ::step].compute()
    else:
        data_sample = data

    # Calculate heuristic quality metrics
    metrics = _calculate_heuristic_metrics(data_sample)

    # Calculate overall score
    heuristic_score = _calculate_heuristic_score(metrics)

    # DL model assessment (optional)
    dl_score = None
    if use_dl_model:
        try:
            from kintsugi.dl_refinement import ChannelAssessor

            assessor = ChannelAssessor(
                model_path=model_path,
                tile_size=tile_size,
                confidence_threshold=confidence_threshold,
                use_heuristics=True,
            )

            # Get full assessment
            result = assessor.process_image(data, channel)
            dl_score = result.confidence_score
            metrics["dl_confidence"] = dl_score
            metrics["problem_regions"] = len(result.problem_regions)

        except ImportError:
            logger.warning("DL refinement module not available, using heuristics only")
        except Exception as e:
            logger.warning(f"DL assessment failed: {e}")

    # Combined score
    if dl_score is not None:
        final_score = 0.7 * dl_score + 0.3 * heuristic_score
    else:
        final_score = heuristic_score

    # Determine if review is needed
    requires_review = final_score < confidence_threshold

    # Generate recommendations
    recommendations = _generate_recommendations(metrics, final_score)

    _add_history(channel, "assess_quality", {
        "confidence_threshold": confidence_threshold,
        "use_dl_model": use_dl_model,
    })

    return {
        "status": "success",
        "channel": channel,
        "quality_score": round(final_score, 3),
        "requires_review": requires_review,
        "metrics": {k: round(v, 3) if isinstance(v, float) else v for k, v in metrics.items()},
        "recommendations": recommendations,
        "threshold_used": confidence_threshold,
    }


def _calculate_heuristic_metrics(data: np.ndarray) -> Dict[str, float]:
    """Calculate heuristic quality metrics for an image."""
    metrics = {}

    # Basic intensity statistics
    metrics["min_intensity"] = float(np.min(data))
    metrics["max_intensity"] = float(np.max(data))
    metrics["mean_intensity"] = float(np.mean(data))
    metrics["std_intensity"] = float(np.std(data))

    # Dynamic range
    metrics["dynamic_range"] = metrics["max_intensity"] - metrics["min_intensity"]

    # SNR estimation
    signal = np.percentile(data, 90)
    background = np.percentile(data, 10)
    noise_region = data[data < np.percentile(data, 20)]
    noise_std = np.std(noise_region) if len(noise_region) > 0 else 1.0

    if noise_std > 0:
        metrics["snr"] = float(signal / noise_std)
    else:
        metrics["snr"] = 0.0

    metrics["signal_background_ratio"] = float(signal / max(background, 1))

    # Contrast metrics
    if metrics["mean_intensity"] > 0:
        metrics["coefficient_of_variation"] = metrics["std_intensity"] / metrics["mean_intensity"]
    else:
        metrics["coefficient_of_variation"] = 0.0

    # Histogram spread (entropy proxy)
    hist, _ = np.histogram(data.ravel(), bins=256)
    hist = hist / hist.sum()
    hist = hist[hist > 0]  # Remove zeros for entropy calculation
    metrics["histogram_entropy"] = float(-np.sum(hist * np.log2(hist)))

    # Saturation check
    max_possible = 65535 if data.dtype == np.uint16 else 255
    saturated_high = np.sum(data >= max_possible * 0.99) / data.size
    saturated_low = np.sum(data == 0) / data.size
    metrics["saturation_ratio"] = float(saturated_high + saturated_low)

    # Sharpness (gradient magnitude)
    try:
        gy, gx = np.gradient(data.astype(float))
        metrics["sharpness"] = float(np.sqrt(gx**2 + gy**2).mean())
    except Exception:
        metrics["sharpness"] = 0.0

    # Autofluorescence estimate
    # High background relative to signal suggests AF issues
    if signal > 0:
        metrics["autofluorescence_score"] = float(background / signal)
    else:
        metrics["autofluorescence_score"] = 1.0

    return metrics


def _calculate_heuristic_score(metrics: Dict[str, float]) -> float:
    """Calculate overall heuristic quality score from metrics."""
    scores = []

    # SNR score (0-1, where >10 is excellent)
    snr_score = min(1.0, metrics.get("snr", 0) / 10.0)
    scores.append(snr_score * 0.3)  # 30% weight

    # Dynamic range score
    max_possible = 65535  # Assume 16-bit
    range_score = min(1.0, metrics.get("dynamic_range", 0) / (0.3 * max_possible))
    scores.append(range_score * 0.2)  # 20% weight

    # Low saturation is good
    saturation_score = 1.0 - min(1.0, metrics.get("saturation_ratio", 0) * 10)
    scores.append(saturation_score * 0.15)  # 15% weight

    # Sharpness score (normalized)
    sharpness_score = min(1.0, metrics.get("sharpness", 0) / 100.0)
    scores.append(sharpness_score * 0.15)  # 15% weight

    # Low autofluorescence is good
    af_score = 1.0 - min(1.0, metrics.get("autofluorescence_score", 1.0))
    scores.append(af_score * 0.2)  # 20% weight

    return sum(scores)


def _generate_recommendations(metrics: Dict[str, float], score: float) -> List[str]:
    """Generate processing recommendations based on metrics."""
    recommendations = []

    # SNR issues
    if metrics.get("snr", 0) < 5:
        recommendations.append("Low SNR detected. Consider denoising (median filter) or longer exposure.")

    # High autofluorescence
    if metrics.get("autofluorescence_score", 0) > 0.3:
        recommendations.append("High background detected. Recommend blank subtraction with scale_factor=0.8-1.0.")

    # Saturation issues
    if metrics.get("saturation_ratio", 0) > 0.05:
        recommendations.append("Saturation detected. Consider reducing exposure or checking sample preparation.")

    # Low dynamic range
    if metrics.get("dynamic_range", 0) < 5000:
        recommendations.append("Low dynamic range. CLAHE may help enhance contrast.")

    # Low sharpness
    if metrics.get("sharpness", 0) < 20:
        recommendations.append("Image appears blurry. Check focus or consider deconvolution.")

    # Overall quality
    if score >= 0.85:
        recommendations.insert(0, "✓ Image quality is good. Minimal processing needed.")
    elif score >= 0.7:
        recommendations.insert(0, "⚠ Image quality is acceptable. Some processing recommended.")
    else:
        recommendations.insert(0, "✗ Image quality needs attention. Review recommended.")

    return recommendations


async def compute_snr(channel: str) -> Dict[str, Any]:
    """
    Compute Signal-to-Noise Ratio for a channel.

    Uses the top 20 / bottom 10% method common in microscopy.
    """
    try:
        data = _get_image(channel)
    except ValueError as e:
        return {"error": str(e)}

    # Compute if dask array
    if hasattr(data, "compute"):
        data = data[::5, ::5].compute()  # Sample for speed

    # Flatten and sort
    flat = data.ravel()

    # Top 20 values (signal)
    top_n = 20
    top_values = np.sort(flat)[-top_n:]
    signal_mean = np.mean(top_values)

    # Bottom 10% (background/noise)
    bottom_10_pct = int(len(flat) * 0.1)
    bottom_values = np.sort(flat)[:bottom_10_pct]
    background_mean = np.mean(bottom_values)
    background_std = np.std(bottom_values)

    # Calculate SNR
    if background_std > 0:
        snr = signal_mean / background_std
    else:
        snr = 0.0

    # Signal-to-background ratio
    if background_mean > 0:
        sbr = signal_mean / background_mean
    else:
        sbr = float("inf")

    return {
        "status": "success",
        "channel": channel,
        "snr": round(snr, 2),
        "signal_to_background": round(sbr, 2),
        "signal_mean": round(signal_mean, 2),
        "background_mean": round(background_mean, 2),
        "background_std": round(background_std, 2),
        "interpretation": _interpret_snr(snr),
    }


def _interpret_snr(snr: float) -> str:
    """Interpret SNR value."""
    if snr >= 20:
        return "Excellent - very clean signal"
    elif snr >= 10:
        return "Good - suitable for quantitative analysis"
    elif snr >= 5:
        return "Acceptable - may need denoising"
    elif snr >= 2:
        return "Poor - significant noise, denoising recommended"
    else:
        return "Very poor - consider reacquisition or heavy filtering"


async def batch_assess_quality(
    channels: List[str],
    confidence_threshold: float = 0.85,
) -> Dict[str, Any]:
    """
    Assess quality for multiple channels at once.

    Returns a summary with channels sorted by quality score.
    """
    results = []

    for channel in channels:
        result = await assess_quality(
            channel=channel,
            confidence_threshold=confidence_threshold,
        )
        if result.get("status") == "success":
            results.append({
                "channel": channel,
                "score": result["quality_score"],
                "requires_review": result["requires_review"],
                "snr": result["metrics"].get("snr", 0),
            })
        else:
            results.append({
                "channel": channel,
                "error": result.get("error", "Unknown error"),
            })

    # Sort by score (lowest first - needs attention)
    valid_results = [r for r in results if "score" in r]
    valid_results.sort(key=lambda x: x["score"])

    # Summary statistics
    scores = [r["score"] for r in valid_results]
    needs_review = [r for r in valid_results if r.get("requires_review")]

    return {
        "status": "success",
        "total_channels": len(channels),
        "channels_assessed": len(valid_results),
        "channels_needing_review": len(needs_review),
        "average_score": round(np.mean(scores), 3) if scores else 0,
        "min_score": round(min(scores), 3) if scores else 0,
        "max_score": round(max(scores), 3) if scores else 0,
        "results": results,
        "review_queue": [r["channel"] for r in needs_review],
    }
