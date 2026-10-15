"""MCP prompts: reusable, parameterised instructions the client can offer the user.

Kept free of any ``mcp`` import so the text can be rendered and tested without
the optional ``mcp`` extra; :mod:`kintsugi.mcp.server` registers it.
"""

from __future__ import annotations

PROMPT_NAME = "signal_isolation_workflow"
PROMPT_TITLE = "Guided signal isolation"
PROMPT_DESCRIPTION = (
    "Step-by-step autofluorescence removal for one marker channel: load, inspect, "
    "suggest, preview, apply, assess, save, and record what worked."
)


def signal_isolation_workflow(
    channel: str,
    blank_channel: str,
    tissue_type: str = "unknown",
    project_path: str | None = None,
) -> str:
    """Render the guided signal-isolation checklist for a marker/blank pair."""
    project = project_path or "the project the server was started for"
    return f"""Isolate the {channel} signal from autofluorescence using the blank channel \
{blank_channel} in {project} (tissue: {tissue_type}). Work through these steps, \
reporting what you see after each, and stop for confirmation before step 6.

1. list_channels, then load_channel for {channel} and {blank_channel} (same cycle).
2. get_image_stats and get_thumbnail for both; describe intensity range, background \
and any artefacts.
3. suggest_with_learning(channel={channel!r}, tissue_type={tissue_type!r}, \
marker_name={channel!r}, blank_channel={blank_channel!r}) and report the confidence \
and whether the suggestion comes from learned history or heuristics.
4. analyze_weighted_subtraction to preview per-intensity-range weights.
5. subtract_blank with method="weighted" (fall back to "global" only if the preview \
looks wrong); optionally denoise or clean_background.
6. assess_quality and compute_snr on the result; show a get_thumbnail; ask whether the \
result is acceptable.
7. If approved: save_processed, then approve_and_learn with the operations and \
parameters you used so the next session starts from them. If not, adjust one \
parameter at a time and repeat from step 5.
"""
