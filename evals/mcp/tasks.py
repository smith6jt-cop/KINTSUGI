"""Realistic tool-selection tasks for the KINTSUGI MCP evaluation harness.

Each task is a natural-language request a scientist might make to Claude while
driving the MCP server. ``expected`` is the set of *acceptable* first tool calls
(more than one is allowed where genuinely ambiguous). Several tasks are
deliberately chosen to probe confusable pairs on the tool surface
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
    {
        "id": "clahe",
        "prompt": "The tissue structures are hard to see. Boost the local contrast of the loaded CD31 channel.",
        "expected": {"apply_clahe"},
        "rationale": "Local contrast enhancement = CLAHE.",
    },
    {
        "id": "clean_background",
        "prompt": (
            "Zero out the dim background pixels below 100 counts and drop the tiny specks "
            "smaller than 30 pixels from the loaded channel."
        ),
        "expected": {"clean_background"},
        "rationale": "Threshold + small-object removal = clean_background.",
    },
    {
        "id": "jupyter_tuner",
        "prompt": (
            "Give me a notebook cell with sliders so I can tune the blank-subtraction "
            "parameters for CD3 interactively."
        ),
        "expected": {"generate_jupyter_cell"},
        "rationale": "Interactive widget code for a notebook = generate_jupyter_cell.",
    },
    {
        "id": "image_stats",
        "prompt": "What are the min, max, mean and intensity histogram of the loaded DAPI channel?",
        "expected": {"get_image_stats"},
        "rationale": "Numeric statistics (not a visual preview, not a quality score).",
    },
    {
        "id": "learning_stats",
        "prompt": "How many parameter records does the learning database hold, and for which tissues?",
        "expected": {"get_learning_statistics"},
        "rationale": "Database-level statistics, not a recommendation.",
    },
    {
        "id": "history",
        "prompt": "List every processing step that has been applied to the CD8 channel so far.",
        "expected": {"get_processing_history"},
        "rationale": "Session provenance for one channel.",
    },
    {
        "id": "predict",
        "prompt": (
            "Use the trained predictor model to estimate the signal-isolation parameters for "
            "the loaded CD4 channel from its image features."
        ),
        "expected": {"predict_parameters"},
        "rationale": "Model-based prediction (Random Forest), not heuristic suggestion or search.",
    },
    {
        "id": "propagate",
        "prompt": (
            "Apply the parameters I tuned on the representative channel to every other "
            "channel in cluster 2."
        ),
        "expected": {"propagate_parameters"},
        "rationale": "Apply tuned parameters across a cluster's members.",
    },
    {
        "id": "suggest_heuristic",
        "prompt": (
            "Based only on the image itself (no learned history), what processing parameters "
            "would you start with for the loaded CD68 channel?"
        ),
        "expected": {"suggest_parameters"},
        "rationale": "Image-analysis heuristic suggestion without the learning database.",
    },
    {
        "id": "session_state",
        "prompt": "Which channels are currently loaded in this session, and how much memory do they use?",
        "expected": {"get_session_state"},
        "rationale": "Session introspection (not the project's channel list on disk).",
    },
    {
        "id": "unload",
        "prompt": "I'm done with the DAPI channel; free its memory on the server.",
        "expected": {"unload_channel"},
        "rationale": "Drop a loaded channel from the session.",
    },
    {
        "id": "batch_quality",
        "prompt": (
            "Score all five channels I loaded from cycle 3 in one go and tell me which "
            "ones need a second look."
        ),
        "expected": {"batch_assess_quality"},
        "rationale": "Several channels at once with a review queue; not one assess_quality per channel.",
    },
    {
        "id": "compare_before_after",
        "prompt": (
            "How different is cyc001_CD3_subtracted from the original cyc001_CD3? "
            "Give me the correlation and SSIM between them."
        ),
        "expected": {"compare_channels"},
        "rationale": "Two-channel statistical comparison, not a single-channel quality score.",
    },
    {
        "id": "intensity_profile",
        "prompt": (
            "Is there an illumination gradient across the loaded DAPI channel? Show me the "
            "intensity along a horizontal line through the middle of the image."
        ),
        "expected": {"get_intensity_profile"},
        "rationale": "A line profile answers a gradient question; get_image_stats does not.",
    },
]
