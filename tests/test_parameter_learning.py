"""Tests for the parameter learning engine behind the MCP learning tools.

Covers the alias normalization (token-based, never substring-based), the
visibility and weighting of bootstrap (non-approved) records, and the
per-project engine cache used by the MCP tools.
"""

from __future__ import annotations

import pytest

from kintsugi.claude.parameter_learning import (
    BOOTSTRAP_WEIGHT,
    ParameterLearningEngine,
    normalize_marker_name,
    normalize_tissue_type,
)


@pytest.fixture
def isolated_home(tmp_path, monkeypatch):
    """Keep the user-wide learning database out of the real home directory."""
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("USERPROFILE", str(home))
    return home


@pytest.mark.parametrize(
    ("marker", "expected"),
    [
        ("CD3", "cd3"),
        ("CD3e", "cd3"),
        ("CD31", "cd31"),  # was "cd3" via substring matching
        ("CD34", "cd34"),  # was "cd3"
        ("CD38", "cd38"),  # was "cd3"
        ("CD45", "cd45"),  # was "cd4"
        ("CD45RO", "cd45ro"),
        ("smooth muscle actin", "sma"),  # was "cd4" via the "th" alias
        ("alpha-SMA", "sma"),
        ("PD-L1", "pdl1"),
        ("Ki-67", "ki67"),
        ("DAPI-01", "dapi"),
        ("Hoechst 33342", "dapi"),
        ("anti CD20 B cell", "cd20"),
        ("Blank_2", "blank"),
        ("E-cadherin", "ecad"),
        ("Granzyme B", "granzyme_b"),
        ("SomethingNew", "somethingnew"),
        ("  ", ""),
    ],
)
def test_normalize_marker_name(marker, expected):
    assert normalize_marker_name(marker) == expected


@pytest.mark.parametrize(
    ("tissue", "expected"),
    [
        ("Tonsil", "tonsil"),
        ("palatine tonsil", "tonsil"),
        ("T", "t"),  # was "tonsil" via substring matching
        ("lymph node", "lymph_node"),
        ("LN", "lymph_node"),
        ("human lymph node biopsy", "lymph_node"),
        ("Renal cortex", "kidney"),
        ("adenocarcinoma", "tumor"),
        ("unknown", "unknown"),
    ],
)
def test_normalize_tissue_type(tissue, expected):
    assert normalize_tissue_type(tissue) == expected


def _record(engine, *, approved, scale, improvement=0.3):
    return engine.record_parameters(
        tissue_type="tonsil",
        marker_name="CD3",
        operation="blank_subtraction",
        parameters={"blank_scale_factor": scale, "erosion": 1},
        quality_before=0.5,
        quality_after=0.5 + improvement,
        user_approved=approved,
        user_notes="bootstrap" if not approved else "",
        to_global=False,
    )


def test_bootstrap_records_are_visible_but_capped(isolated_home, tmp_path):
    """Rows written by `kintsugi mcp pretrain` (user_approved=False) inform recommendations."""
    project = tmp_path / "proj"
    project.mkdir()
    engine = ParameterLearningEngine(str(project))
    for scale in (0.8, 0.82, 0.78):
        _record(engine, approved=False, scale=scale)

    rec = engine.recommend_parameters("tonsil", "CD3", "blank_subtraction")
    assert rec["status"] == "success"
    assert rec["bootstrap_records"] == 3
    assert rec["approved_records"] == 0
    assert rec["recommended_parameters"]["blank_scale_factor"] == pytest.approx(0.8, abs=0.02)
    # Bootstrap-only evidence never reaches the "learned" tier (> 0.5).
    assert 0 < rec["confidence"] <= 0.5

    strict = engine.recommend_parameters("tonsil", "CD3", "blank_subtraction", approved_only=True)
    assert strict["status"] == "no_history"


def test_approved_records_outrank_bootstrap(isolated_home, tmp_path):
    project = tmp_path / "proj"
    project.mkdir()
    engine = ParameterLearningEngine(str(project))
    _record(engine, approved=False, scale=0.3)
    _record(engine, approved=True, scale=0.9)

    rec = engine.recommend_parameters("tonsil", "CD3", "blank_subtraction", max_results=5)
    assert rec["approved_records"] == 1 and rec["bootstrap_records"] == 1
    top = rec["similar_records"][0]
    assert top["parameters"]["blank_scale_factor"] == 0.9
    # The weighted aggregate leans towards the approved value.
    assert rec["recommended_parameters"]["blank_scale_factor"] > 0.6
    assert 0 < BOOTSTRAP_WEIGHT < 1


def test_mcp_engine_cache_is_per_project(isolated_home, tmp_path):
    from kintsugi.mcp.tools import learning

    learning._learning_engines.clear()
    a, b = tmp_path / "a", tmp_path / "b"
    a.mkdir()
    b.mkdir()
    engine_a = learning._get_learning_engine(str(a))
    engine_global = learning._get_learning_engine(None)
    engine_b = learning._get_learning_engine(str(b))

    assert engine_a.project_path == a.resolve() or engine_a.project_path == a
    assert engine_b.project_path != engine_a.project_path
    assert engine_global.project_path is None
    assert learning._get_learning_engine(str(a)) is engine_a
    learning._learning_engines.clear()
