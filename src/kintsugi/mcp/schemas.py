"""Typed result shapes for selected MCP tools.

The SDK turns a ``TypedDict`` return annotation into the tool's ``outputSchema``
(the generic ``dict[str, Any]`` used elsewhere only advertises "an object"). These
are ``total=False`` with ``NotRequired`` keys so a handler may omit a key
without the SDK serialising it as ``null``; the ``error`` key is the in-band
failure shape every tool shares. Handler bodies keep returning plain dicts.
"""

from __future__ import annotations

from typing import Any

from pydantic import ConfigDict
from typing_extensions import NotRequired, TypedDict


class Histogram(TypedDict, total=False):
    counts: list[int]
    bin_edges: list[float]
    bin_centers: list[float]


class ImageStatistics(TypedDict, total=False):
    shape: list[int]
    dtype: str
    min: float
    max: float
    mean: float
    std: float
    median: float
    percentiles: dict[str, float]
    histogram: Histogram
    __pydantic_config__ = ConfigDict(extra="allow")  # type: ignore[misc]


class ImageStatsResult(TypedDict, total=False):
    """Result of ``get_image_stats``."""

    status: str
    channel: str
    statistics: ImageStatistics
    error: NotRequired[str]


class QualityResult(TypedDict, total=False):
    """Result of ``assess_quality``."""

    status: str
    channel: str
    quality_score: float
    requires_review: bool
    metrics: dict[str, Any]
    recommendations: list[str]
    threshold_used: float
    error: NotRequired[str]


class OperationSuggestion(TypedDict, total=False):
    recommended_parameters: dict[str, Any]
    confidence: float
    source: str
    rationale: str
    learned_records: int


class LearningSuggestion(TypedDict, total=False):
    """Result of ``suggest_with_learning``."""

    status: str
    channel: str
    tissue_type: str
    tissue_type_normalized: str
    marker_name: str
    marker_name_normalized: str
    image_characteristics: dict[str, float]
    operations: dict[str, OperationSuggestion]
    suggested_order: list[str]
    error: NotRequired[str]


class LoadedImage(TypedDict, total=False):
    shape: list[int]
    dtype: str
    metadata: dict[str, Any]


class MemoryUsage(TypedDict, total=False):
    total_mb: float
    total_gb: float


class WorkflowSummary(TypedDict, total=False):
    """Persistent workflow state of a project (``.kintsugi/workflow_state.json``)."""

    session_id: str | None
    sessions: int
    total_channels: int
    status_breakdown: dict[str, int]
    processing_queue: list[str]
    review_queue: list[str]
    approved: list[str]


class SessionStateResult(TypedDict, total=False):
    """Result of ``get_session_state``."""

    status: str
    project_root: str | None
    loaded_images: dict[str, LoadedImage]
    total_loaded: int
    processing_history: dict[str, int]
    memory_usage: MemoryUsage
    workflow: NotRequired[WorkflowSummary]
    error: NotRequired[str]


class WorkflowChannel(TypedDict, total=False):
    status: str
    cycle: str
    quality_score: float | None
    steps: int
    last_operation: str | None
    review_notes: str
    modified: str


class WorkflowStatusResult(WorkflowSummary, total=False):
    """Result of ``get_workflow_status``."""

    status: str
    project_root: str
    channels: dict[str, WorkflowChannel]
    error: NotRequired[str]
