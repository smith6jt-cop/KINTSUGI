"""Tests for the parameter learning engine behind the MCP learning tools.

Covers the alias normalization (token-based, never substring-based), the
visibility and weighting of bootstrap (non-approved) records, and the
per-project engine cache used by the MCP tools.
"""

from __future__ import annotations

import pytest

from kintsugi.assist.parameter_learning import (
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


def test_recommend_parameters_many_matches_single_queries(isolated_home, tmp_path):
    engine = ParameterLearningEngine(str(tmp_path / "proj"))
    for op, params in (("denoise", {"filter_size": 2}), ("clahe", {"clip_limit": 0.01})):
        for approved in (True, False):
            engine.record_parameters("tonsil", "CD3", op, params, 0.4, 0.8, user_approved=approved)

    operations = ["denoise", "clahe", "gaussian_subtract"]
    many = engine.recommend_parameters_many("tonsil", "CD3", operations)

    assert list(many) == operations
    for op in operations:
        assert many[op] == engine.recommend_parameters("tonsil", "CD3", op)
    assert many["denoise"]["status"] == "success"
    assert many["denoise"]["recommended_parameters"] == {"filter_size": 2}
    assert many["gaussian_subtract"]["status"] == "no_history"

    # The batched window-function query returns the same rows in the same order as
    # the single-operation query it replaces (recommend_parameters delegates to the
    # batched path, so this is the comparison that catches a divergence).
    rows_many = engine._query_similar_records_many(
        engine.project_db_path, "tonsil", "cd3", operations, 0.0, False
    )
    for op in operations:
        rows_single = engine._query_similar_records(
            engine.project_db_path, "tonsil", "cd3", op, 0.0, False
        )
        assert [r["id"] for r in rows_many[op]] == [r["id"] for r in rows_single]
    approved = engine._query_similar_records_many(
        engine.project_db_path, "tonsil", "cd3", ["denoise"], 0.0, True
    )["denoise"]
    assert len(approved) == 1 and approved[0]["user_approved"] is True


def test_record_parameters_many_writes_every_entry_in_order(isolated_home, tmp_path):
    engine = ParameterLearningEngine(str(tmp_path / "proj"))
    entries = [
        {"operation": "blank_subtraction", "parameters": {"method": "weighted"}},
        {"operation": "denoise", "parameters": {"filter_size": 3}, "algorithm_version": "v2"},
    ]

    ids = engine.record_parameters_many(
        entries,
        tissue_type="tonsil",
        marker_name="CD3",
        quality_before=0.3,
        quality_after=0.9,
        user_approved=True,
        channel_name="cyc001_CD3",
    )

    assert len(ids) == 2 and ids[1] > ids[0] > 0
    assert engine.record_parameters_many([], tissue_type="t", marker_name="m") == []
    stats = engine.get_statistics()
    assert stats["project_records"] == 2 and stats["global_records"] == 2
    denoise = engine.recommend_parameters("tonsil", "CD3", "denoise", algorithm_version="v2")
    assert denoise["status"] == "success" and denoise["recommended_parameters"] == {
        "filter_size": 3
    }
    # Without a project database the ids come from the global database.
    global_only = ParameterLearningEngine(None)
    assert global_only.record_parameters_many(
        entries[:1], tissue_type="tonsil", marker_name="CD3"
    ) == [3]


def test_a_mirrored_record_counts_once(isolated_home, tmp_path):
    """record_parameters writes to the project and the user-wide DB; recall counts it once."""
    engine = ParameterLearningEngine(str(tmp_path / "proj"))
    engine.record_parameters(
        "tonsil", "CD3", "denoise", {"filter_size": 2}, 0.4, 0.8, user_approved=True
    )

    rec = engine.recommend_parameters("tonsil", "CD3", "denoise")

    assert rec["based_on_records"] == 1
    assert rec["approved_records"] == 1
    assert len(rec["similar_records"]) == 1
    # A different recording with identical parameters is still a separate record.
    engine.record_parameters(
        "tonsil", "CD3", "denoise", {"filter_size": 2}, 0.3, 0.9, user_approved=True
    )
    assert engine.recommend_parameters("tonsil", "CD3", "denoise")["based_on_records"] == 2
