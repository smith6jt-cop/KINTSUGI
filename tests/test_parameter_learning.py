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
