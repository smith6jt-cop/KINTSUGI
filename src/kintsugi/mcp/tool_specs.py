"""
Declarative tool registry for the KINTSUGI MCP server.

Each entry pairs a public MCP tool name with its handler (``module``/``attr`` in
``kintsugi.mcp.tools``), a human-readable ``title``, the top-level ``description``
advertised to clients, and the MCP tool ``annotations`` (behaviour hints). Keeping
this as plain data (no ``mcp`` import) lets ``kintsugi mcp tools``, the docs
generator, and the tests read the registry without the optional ``[mcp]``
dependency installed.

The *input schema* of every tool is derived by the SDK from the handler's own
signature: parameter descriptions, enums and numeric ranges live in
``typing.Annotated[..., pydantic.Field(...)]`` and ``typing.Literal[...]``
annotations on the handler, where they are also enforced. There is no second,
hand-written copy of the schema here.

Annotation keys use the MCP wire spelling (``readOnlyHint`` etc.). ``openWorldHint``
is false for every tool: the server only touches the local project directory
and its own session state.

Naming: tool names are deliberately *not* prefixed with ``kintsugi_``. They are
valid, unique, verb-first snake_case names (SEP-986), and MCP clients namespace
tools by server (Claude Code, for example, as ``mcp__kintsugi__<tool>``); a prefix
would only duplicate that and break existing permission rules, notebooks and eval tasks.
"""

from __future__ import annotations

from typing import Any

# Behaviour-hint presets (MCP ToolAnnotations, camelCase wire keys).
_READ_ONLY: dict[str, bool] = {
    "readOnlyHint": True,
    "destructiveHint": False,
    "idempotentHint": True,
    "openWorldHint": False,
}
# Adds/replaces an image in the in-memory session; re-running yields the same state.
_PROCESSING: dict[str, bool] = {
    "readOnlyHint": False,
    "destructiveHint": False,
    "idempotentHint": True,
    "openWorldHint": False,
}
# Appends a record (learning DB) or runs a stochastic search; repeats accumulate.
_APPENDS: dict[str, bool] = {
    "readOnlyHint": False,
    "destructiveHint": False,
    "idempotentHint": False,
    "openWorldHint": False,
}
# Writes files under the project and may overwrite existing output.
_WRITES_FILES: dict[str, bool] = {
    "readOnlyHint": False,
    "destructiveHint": True,
    "idempotentHint": False,
    "openWorldHint": False,
}

TOOL_SPECS: list[dict[str, Any]] = [
    # ---- Signal isolation ----
    {
        "name": "load_channel",
        "module": "signal_isolation",
        "attr": "load_channel",
        "title": "Load channel",
        "description": (
            "Load a channel image from a KINTSUGI project. Returns image metadata and prepares it "
            "for processing."
        ),
        "annotations": _PROCESSING,
    },
    {
        "name": "subtract_blank",
        "module": "signal_isolation",
        "attr": "subtract_blank",
        "title": "Subtract blank (autofluorescence)",
        "description": (
            "Subtract autofluorescence/blank channel from signal. Uses optimized Dask operations "
            "for large images."
        ),
        "annotations": _PROCESSING,
    },
    {
        "name": "analyze_weighted_subtraction",
        "module": "signal_isolation",
        "attr": "analyze_weighted_subtraction",
        "title": "Preview weighted subtraction",
        "description": (
            "Preview the per-intensity-range weights for weighted autofluorescence subtraction "
            "WITHOUT applying it. Use this to verify the plan before calling subtract_blank with "
            "method='weighted'."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "denoise",
        "module": "signal_isolation",
        "attr": "denoise",
        "title": "Denoise (basic filters)",
        "description": "Apply denoising filters (percentile, uniform, median) to reduce noise.",
        "annotations": _PROCESSING,
    },
    {
        "name": "apply_clahe",
        "module": "signal_isolation",
        "attr": "apply_clahe",
        "title": "Apply CLAHE contrast",
        "description": "Apply Contrast Limited Adaptive Histogram Equalization to enhance local contrast.",
        "annotations": _PROCESSING,
    },
    {
        "name": "clean_background",
        "module": "signal_isolation",
        "attr": "clean_background",
        "title": "Clean background",
        "description": "Remove background and small objects from the image.",
        "annotations": _PROCESSING,
    },
    {
        "name": "gaussian_subtract",
        "module": "signal_isolation",
        "attr": "gaussian_subtract",
        "title": "Gaussian background subtraction",
        "description": "Subtract Gaussian-blurred version to remove structured background/autofluorescence.",
        "annotations": _PROCESSING,
    },
    {
        "name": "denoise_advanced",
        "module": "signal_isolation",
        "attr": "denoise_advanced",
        "title": "Denoise (advanced)",
        "description": (
            "Apply advanced denoising: adaptive (auto-select best filter), bilateral (edge- "
            "preserving), nlm (non-local means), bm3d (block matching 3D), n2v (self-supervised "
            "Noise2Void), or patch_svd."
        ),
        "annotations": _PROCESSING,
    },
    # ---- Quality assessment ----
    {
        "name": "assess_quality",
        "module": "quality_assessment",
        "attr": "assess_quality",
        "title": "Assess quality",
        "description": (
            "Assess channel quality with heuristic metrics (SNR, contrast, dynamic range, "
            "background) and flag channels that need manual review. Returns a 0-1 quality "
            "score, the metrics and recommendations."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "compute_snr",
        "module": "quality_assessment",
        "attr": "compute_snr",
        "title": "Compute SNR",
        "description": "Compute Signal-to-Noise Ratio for a channel.",
        "annotations": _READ_ONLY,
    },
    {
        "name": "batch_assess_quality",
        "module": "quality_assessment",
        "attr": "batch_assess_quality",
        "title": "Batch quality assessment",
        "description": (
            "Assess several loaded channels at once: scores sorted worst-first and a review "
            "queue of the channels below the threshold."
        ),
        "annotations": _READ_ONLY,
    },
    # ---- Visualization ----
    {
        "name": "get_image_stats",
        "module": "visualization",
        "attr": "get_image_stats",
        "title": "Image statistics",
        "description": "Get statistics for a loaded image (min, max, mean, std, histogram).",
        "annotations": _READ_ONLY,
    },
    {
        "name": "get_thumbnail",
        "module": "visualization",
        "attr": "get_thumbnail",
        "title": "Thumbnail preview",
        "description": "Get a downsampled thumbnail of the image for preview.",
        "annotations": _READ_ONLY,
    },
    {
        "name": "compare_channels",
        "module": "visualization",
        "attr": "compare_channels",
        "title": "Compare two channels",
        "description": (
            "Compare two loaded channels statistically (means, difference, correlation, SSIM) "
            "with a one-line interpretation; useful for before/after checks."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "get_intensity_profile",
        "module": "visualization",
        "attr": "get_intensity_profile",
        "title": "Intensity profile",
        "description": (
            "Intensity values along a horizontal or vertical line through a loaded channel, "
            "for spotting gradients or uneven illumination."
        ),
        "annotations": _READ_ONLY,
    },
    # ---- Workflow ----
    {
        "name": "list_channels",
        "module": "workflow",
        "attr": "list_channels",
        "title": "List project channels",
        "description": "List all available channels in the current project/cycle.",
        "annotations": _READ_ONLY,
    },
    {
        "name": "save_processed",
        "module": "workflow",
        "attr": "save_processed",
        "title": "Save processed channel",
        "description": "Save the processed channel to the project output directory.",
        "annotations": _WRITES_FILES,
    },
    {
        "name": "get_processing_history",
        "module": "workflow",
        "attr": "get_processing_history",
        "title": "Processing history",
        "description": "Get the processing history for a channel, showing all operations applied.",
        "annotations": _READ_ONLY,
    },
    {
        "name": "suggest_parameters",
        "module": "workflow",
        "attr": "suggest_parameters",
        "title": "Suggest parameters (heuristic)",
        "description": (
            "Rule-based suggestions for every operation from the channel's intensity "
            "statistics, without the learning database. Prefer suggest_with_learning when "
            "the tissue type and marker are known."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "generate_jupyter_cell",
        "module": "workflow",
        "attr": "generate_jupyter_cell",
        "title": "Generate Jupyter tuning cell",
        "description": (
            "Generate a Jupyter notebook cell for interactive parameter tuning with Kview2 widgets."
        ),
        "annotations": _READ_ONLY,
    },
    # ---- Signal isolation ----
    {
        "name": "cluster_channels",
        "module": "signal_isolation",
        "attr": "cluster_channels_tool",
        "title": "Cluster channels",
        "description": (
            "Cluster channels by image feature similarity. Groups channels that share similar "
            "intensity profiles, SNR, and texture so parameters can be tuned once per cluster."
        ),
        "annotations": _PROCESSING,
    },
    {
        "name": "propagate_parameters",
        "module": "signal_isolation",
        "attr": "propagate_parameters_tool",
        "title": "Propagate cluster parameters",
        "description": (
            "Apply tuned signal isolation parameters from a representative channel to all members "
            "of a cluster."
        ),
        "annotations": _WRITES_FILES,
    },
    {
        "name": "get_session_state",
        "module": "workflow",
        "attr": "get_session_state",
        "title": "Session state",
        "description": (
            "Show what this server session holds: the project root it serves, every loaded "
            "channel with shape/dtype/metadata, per-channel processing step counts, and memory "
            "usage."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "unload_channel",
        "module": "workflow",
        "attr": "unload_channel",
        "title": "Unload channel",
        "description": (
            "Free server memory by dropping a loaded channel (or everything with clear_all). "
            "Files on disk are untouched."
        ),
        "annotations": _PROCESSING,
    },
    {
        "name": "get_workflow_status",
        "module": "workflow",
        "attr": "get_workflow_status",
        "title": "Workflow status",
        "description": (
            "Show the project's persistent workflow state, which survives server restarts: "
            "every channel's status (pending, processing, review, approved, rejected, "
            "needs_adjustment), quality score and step count, plus the processing and review "
            "queues. Call it first to resume earlier work."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "set_review_decision",
        "module": "workflow",
        "attr": "set_review_decision",
        "title": "Record review decision",
        "description": (
            "Record a review decision other than approval (rejected, needs_adjustment, or "
            "pending to reset) for a channel in the project's workflow state; works for "
            "channels that are no longer loaded. Use approve_and_learn to approve."
        ),
        "annotations": _APPENDS,
    },
    # ---- Parameter learning ----
    {
        "name": "get_learned_parameters",
        "module": "learning",
        "attr": "get_learned_parameters",
        "title": "Learned parameters",
        "description": (
            "Get recommended parameters based on learned history from similar tissue/marker "
            "combinations. Returns parameters with confidence scores."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "record_successful_parameters",
        "module": "learning",
        "attr": "record_successful_parameters",
        "title": "Record parameters",
        "description": (
            "Record a successful parameter set for future learning. Call after user approves "
            "processing results."
        ),
        "annotations": _APPENDS,
    },
    {
        "name": "suggest_with_learning",
        "module": "learning",
        "attr": "suggest_with_learning",
        "title": "Suggest parameters (learned + heuristic)",
        "description": (
            "Get comprehensive parameter suggestions combining image analysis AND learned "
            "history. This is the primary recommendation tool."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "approve_and_learn",
        "module": "learning",
        "attr": "approve_and_learn",
        "title": "Approve and learn",
        "description": (
            "Approve processing results and record ALL parameters for learning. Call when user "
            "approves final processed result."
        ),
        "annotations": _APPENDS,
    },
    {
        "name": "get_learning_statistics",
        "module": "learning",
        "attr": "get_learning_statistics",
        "title": "Learning database statistics",
        "description": "Get statistics about the parameter learning database.",
        "annotations": _READ_ONLY,
    },
    # ---- Signal isolation ----
    {
        "name": "optimize_parameters",
        "module": "signal_isolation",
        "attr": "optimize_parameters",
        "title": "Optimize parameters (Optuna)",
        "description": (
            "Find optimal signal isolation parameters via Bayesian optimization (Optuna). "
            "Automatically searches the parameter space — no manual slider adjustment needed."
        ),
        "annotations": _APPENDS,
    },
    {
        "name": "predict_parameters",
        "module": "signal_isolation",
        "attr": "predict_parameters",
        "title": "Predict parameters (Random Forest)",
        "description": (
            "Predict signal isolation parameters from image features using a trained Random "
            "Forest model. Includes confidence scores and uncertainty estimates."
        ),
        "annotations": _READ_ONLY,
    },
    {
        "name": "estimate_background",
        "module": "signal_isolation",
        "attr": "estimate_background",
        "title": "Estimate background (SMO)",
        "description": (
            "Estimate and subtract background using SMO (Silver Mountain Operator). Parameter- "
            "free — no blank channel needed."
        ),
        "annotations": _PROCESSING,
    },
]


def tool_names() -> list[str]:
    """Public tool names in registry order."""
    return [spec["name"] for spec in TOOL_SPECS]
