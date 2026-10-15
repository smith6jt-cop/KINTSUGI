"""Rule-based parameter suggestions shared by the suggestion tools.

``suggest_parameters`` (workflow) and ``suggest_with_learning`` (learning) used
to carry two diverging copies of the same thresholds. Both now call
:func:`heuristic_suggestion`, so a rule changes in exactly one place and the two
tools cannot disagree about what a "noisy" or "low dynamic range" image is.

Everything here is pure NumPy; nothing imports ``mcp``.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from typing import Any

import numpy as np

logger = logging.getLogger("kintsugi.mcp.heuristics")

#: Operations the heuristics know about, in the order the tools report them.
#: Identical to ``kintsugi.claude.parameter_learning.LEARNING_OPERATIONS``.
HEURISTIC_OPERATIONS: tuple[str, ...] = (
    "blank_subtraction",
    "denoise",
    "clahe",
    "clean_background",
    "gaussian_subtract",
)


@dataclass(frozen=True)
class IntensitySummary:
    """Percentile summary of a (sub-sampled) image used by every rule below."""

    p1: float
    p10: float
    p50: float
    p90: float
    p99: float
    noise_estimate: float

    @property
    def dynamic_range(self) -> float:
        return self.p99 - self.p1

    def to_dict(self) -> dict[str, float]:
        return {
            "p1": self.p1,
            "p10": self.p10,
            "p50": self.p50,
            "p90": self.p90,
            "p99": self.p99,
            "dynamic_range": self.dynamic_range,
            "noise_estimate": self.noise_estimate,
        }


def sample_image(data: Any, target_rows: int = 2000) -> np.ndarray:
    """Sub-sample a loaded image (dask or array) to about ``target_rows`` rows.

    The same stride is used for the signal and its blank so paired statistics
    (means, correlations, weighted-range analysis) are computed on matching grids.
    """
    step = max(1, int(data.shape[0]) // target_rows)
    sample = data[::step, ::step]
    return sample.compute() if hasattr(sample, "compute") else np.asarray(sample)


def summarize_intensity(sample: np.ndarray) -> IntensitySummary:
    """Percentiles and a background-noise estimate (std of the darkest decile)."""
    p1, p10, p50, p90, p99 = (float(v) for v in np.percentile(sample, [1, 10, 50, 90, 99]))
    low = sample[sample < p10]
    noise = float(np.std(low)) if low.size else 0.0
    return IntensitySummary(p1, p10, p50, p90, p99, noise)


def heuristic_suggestion(
    operation: str,
    sample: np.ndarray,
    summary: IntensitySummary,
    blank_sample: np.ndarray | None = None,
    method: str = "global",
) -> dict[str, Any]:
    """Return ``{"parameters": ..., "rationale": ...}`` for ``operation``.

    ``blank_sample`` is the blank/autofluorescence channel sampled with the same
    stride as ``sample``; without it the blank-subtraction rule falls back to
    defaults. ``method="weighted"`` asks for the weighted-subtraction plan.
    """
    p10 = summary.p10
    dynamic_range = summary.dynamic_range
    noise_estimate = summary.noise_estimate

    if operation == "blank_subtraction" and method == "weighted":
        if blank_sample is not None:
            try:
                from kintsugi.signal.autofluorescence import analyze_for_weighted_subtraction

                analysis = analyze_for_weighted_subtraction(sample, blank_sample)
                return {
                    "parameters": {
                        "method": "weighted",
                        "base_scale_factor": analysis["base_scale_factor"],
                        "blank_clip_factor": analysis["blank_clip_factor"],
                        "n_ranges": analysis["n_ranges"],
                        "range_method": analysis["range_method"],
                        "ranges": analysis["ranges"],
                        "transition_width": analysis["transition_width"],
                        "smooth_low": analysis["smooth_low"],
                        "smooth_high": analysis["smooth_high"],
                        "erosion": analysis["erosion"],
                    },
                    "rationale": (
                        f"Weighted: {analysis['n_ranges']} ranges, "
                        f"dim weight={analysis['ranges'][1]['weight']:.2f} "
                        f"bright weight={analysis['ranges'][-1]['weight']:.2f}"
                        if len(analysis["ranges"]) > 1
                        else "Weighted subtraction"
                    ),
                }
            except Exception as e:
                logger.warning(f"Weighted analysis failed, falling back: {e}")
        return {
            "parameters": {
                "method": "weighted",
                "base_scale_factor": 1.0,
                "blank_clip_factor": 0,
                "n_ranges": 5,
                "range_method": "percentile",
                "transition_width": 0.1,
            },
            "rationale": "Default weighted parameters (no blank loaded)",
        }

    if operation == "blank_subtraction":
        if blank_sample is not None:
            blank_mean = float(np.mean(blank_sample))
            signal_mean = float(np.mean(sample))
            if blank_mean > 0:
                suggested_scale = min(1.0, max(0.5, signal_mean / blank_mean * 0.8))
            else:
                suggested_scale = 0.8
            return {
                "parameters": {
                    "blank_scale_factor": round(suggested_scale, 2),
                    "blank_clip_factor": int(np.percentile(blank_sample, 5)),
                    "smooth_low": noise_estimate > 300,
                    "smooth_high": False,
                    "erosion": 1 if suggested_scale > 0.7 else 0,
                },
                "rationale": f"Blank mean={blank_mean:.0f}, Signal mean={signal_mean:.0f}",
            }
        return {
            "parameters": {
                "blank_scale_factor": 0.8,
                "blank_clip_factor": 0,
                "smooth_low": False,
                "smooth_high": False,
                "erosion": 0,
            },
            "rationale": "No blank channel loaded - using defaults",
        }

    if operation == "denoise":
        if noise_estimate > 500:
            return {
                "parameters": {"method": "median", "filter_size": 3},
                "rationale": f"High noise (std={noise_estimate:.0f})",
            }
        if noise_estimate > 200:
            return {
                "parameters": {"method": "median", "filter_size": 2},
                "rationale": f"Moderate noise (std={noise_estimate:.0f})",
            }
        return {
            "parameters": {"method": "median", "filter_size": 1},
            "rationale": f"Low noise (std={noise_estimate:.0f}) - minimal filtering",
        }

    if operation == "clahe":
        if dynamic_range < 10000:
            return {
                "parameters": {"clip_limit": 0.02, "tile_grid_size": 50, "nbins": 128},
                "rationale": f"Low dynamic range ({dynamic_range:.0f}) - aggressive CLAHE",
            }
        if dynamic_range < 30000:
            return {
                "parameters": {"clip_limit": 0.01, "tile_grid_size": 70, "nbins": 128},
                "rationale": f"Moderate dynamic range ({dynamic_range:.0f})",
            }
        return {
            "parameters": {"clip_limit": 0.005, "tile_grid_size": 100, "nbins": 256},
            "rationale": f"Good dynamic range ({dynamic_range:.0f}) - gentle CLAHE",
        }

    if operation == "clean_background":
        return {
            "parameters": {
                "background_threshold": int(p10 * 0.5),
                "remove_small_objects": True,
                "min_object_size": 30,
            },
            "rationale": f"Background level at p10={p10:.0f}",
        }

    if operation == "gaussian_subtract":
        return {
            "parameters": {
                "sigma": 10,
                "scale_factor": 0.1 if dynamic_range > 20000 else 0.2,
            },
            "rationale": "Standard Gaussian subtraction",
        }

    raise ValueError(f"Unknown operation: {operation!r} (expected one of {HEURISTIC_OPERATIONS})")
