"""Realistic tool-selection tasks for the KINTSUGI MCP evaluation harness.

Each task is a natural-language request a scientist might make to Claude while
driving the MCP server. ``expected`` is the set of *acceptable* first tool calls
(more than one is allowed where genuinely ambiguous). Several tasks are
deliberately chosen to probe confusable pairs on the 26-tool surface
(``denoise`` vs ``denoise_advanced``; the three background-removal tools;
``assess_quality`` vs ``compute_snr``; the learning-recall vs record tools).

Keeping these as plain data (no MCP/anthropic import) lets the harness load them
in any environment.
"""

from __future__ import annotations

from typing import TypedDict


class EvalTask(TypedDict):
    id: str
    prompt: str
    expected: set[str]
    rationale: str


TASKS: list[EvalTask] = [
    {
        "id": "load_channel",
        "prompt": "Load the CD8 channel from cycle 2 of my project at /data/proj.",
        "expected": {"load_channel"},
        "rationale": "Direct load request with project/cycle/channel.",
    },
    {
        "id": "list_channels",
        "prompt": "What channels are available in cycle 1 of the project at /data/proj?",
        "expected": {"list_channels"},
        "rationale": "Enumerate available channels for a project/cycle.",
    },
    {
        "id": "quality",
        "prompt": "How good is the CD3 channel I just loaded? Give me a quality score and SNR.",
        "expected": {"assess_quality", "compute_snr"},
        "rationale": "Quality assessment; assess_quality is the richer answer, compute_snr acceptable.",
    },
    {
        "id": "subtract_blank_global",
        "prompt": "Remove the autofluorescence from my loaded CD20 channel using the Blank channel.",
        "expected": {"subtract_blank"},
        "rationale": "Blank-channel AF subtraction with a paired blank.",
    },
    {
        "id": "analyze_weighted_preview",
        "prompt": (
            "Before applying anything, show me the per-intensity-range weights you would use "
            "for weighted autofluorescence subtraction of CD20 against the Blank channel."
        ),
        "expected": {"analyze_weighted_subtraction"},
        "rationale": "Preview-only request (verification-first); must NOT call subtract_blank yet.",
    },
    {
        "id": "denoise_auto",
        "prompt": "This channel is really noisy. Automatically pick the best denoiser and clean it up.",
        "expected": {"denoise_advanced"},
        "rationale": "Auto/adaptive denoising = denoise_advanced, NOT the basic denoise tool.",
    },
    {
        "id": "denoise_basic",
        "prompt": "Just apply a simple median filter to smooth the loaded channel.",
        "expected": {"denoise"},
        "rationale": "Explicit simple median filter = basic denoise, NOT denoise_advanced.",
    },
    {
        "id": "estimate_background_no_blank",
        "prompt": "I don't have a blank channel. Estimate and subtract the background automatically.",
        "expected": {"estimate_background"},
        "rationale": "Parameter-free, blank-free background removal = estimate_background (SMO).",
    },
    {
        "id": "gaussian_subtract",
        "prompt": "Remove the smooth structured background by subtracting a Gaussian-blurred version.",
        "expected": {"gaussian_subtract"},
        "rationale": "Explicit Gaussian-blur subtraction.",
    },
    {
        "id": "thumbnail",
        "prompt": "Show me a small downsampled preview image of the loaded DAPI channel.",
        "expected": {"get_thumbnail"},
        "rationale": "Visual preview = thumbnail (not get_image_stats).",
    },
    {
        "id": "learned_recall",
        "prompt": "What blank-subtraction parameters worked before for CD3 in tonsil tissue?",
        "expected": {"get_learned_parameters", "suggest_with_learning"},
        "rationale": "Recall learned history by tissue/marker/operation.",
    },
    {
        "id": "record_success",
        "prompt": (
            "That result looks great. Record these blank-subtraction parameters so we reuse "
            "them next time for CD3 in tonsil."
        ),
        "expected": {"record_successful_parameters", "approve_and_learn"},
        "rationale": "Persist approved params to the learning DB.",
    },
    {
        "id": "cluster",
        "prompt": "Group my channels so I only have to tune parameters once per similar group.",
        "expected": {"cluster_channels"},
        "rationale": "Feature-similarity clustering of channels.",
    },
    {
        "id": "optimize",
        "prompt": (
            "Automatically search for the best autofluorescence-subtraction parameters for CD8 "
            "vs Blank. I don't want to adjust sliders by hand."
        ),
        "expected": {"optimize_parameters"},
        "rationale": "Automated parameter search (Optuna) = optimize_parameters.",
    },
    {
        "id": "save",
        "prompt": "Save the processed CD8 channel as an OME-TIFF.",
        "expected": {"save_processed"},
        "rationale": "Persist a processed channel to disk.",
    },
]
