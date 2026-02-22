"""
Tests for multi-project batch signal isolation.

Tests cover project discovery, tissue type parsing, cascading recipe resolution,
learned params to recipe conversion, and multi-project orchestration.
"""

from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import numpy as np
import pytest
import tifffile
import yaml

from kintsugi.signal.batch import (
    CleanParams,
    MarkerRecipe,
    SubtractionParams,
)
from kintsugi.signal.batch_multi import (
    MultiProjectResult,
    ProjectInfo,
    ResolvedRecipes,
    _write_recipes_as_param_files,
    discover_projects,
    learned_params_to_recipe,
    parse_tissue_type,
    resolve_recipes_for_project,
)

# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture
def projects_dir(tmp_path):
    """Create a multi-project directory structure for testing."""
    # Project 1: Spleen with config + registered images + own recipes
    p1 = tmp_path / "CX_19-001_SP_CC2-A28"
    _create_project(
        p1,
        tissue="spleen",
        cycles=3,
        markers=["CD3e", "CD20", "CD45"],
        with_recipes=True,
    )

    # Project 2: Lymph node with config + registered images, no recipes
    p2 = tmp_path / "CX_19-002_lymph-node_R1"
    _create_project(
        p2,
        tissue="lymph_node",
        cycles=2,
        markers=["CD3e", "CD8"],
    )

    # Project 3: Thymus with config + registered images, no recipes
    p3 = tmp_path / "CX_20-005_TH_CC1-A"
    _create_project(
        p3,
        tissue="thymus",
        cycles=2,
        markers=["CD3e", "CD4"],
    )

    # Project 4: Already completed (has manifest)
    p4 = tmp_path / "CX_19-003_spleen_CC1-A"
    _create_project(p4, tissue="spleen", cycles=2, markers=["CD3e"])
    manifest_dir = p4 / "data" / "processed" / "signal_isolated"
    manifest_dir.mkdir(parents=True)
    (manifest_dir / "signal_isolation_manifest.json").write_text(
        json.dumps({"timestamp": "2026-01-01", "channels": {}})
    )

    # Project 5: No config (should be skipped)
    p5 = tmp_path / "src_CX_19-004_LN_R2"
    p5.mkdir()

    # Non-directory entry
    (tmp_path / "setup.log").write_text("log")

    return tmp_path


def _create_project(
    path: Path,
    tissue: str,
    cycles: int,
    markers: list[str],
    with_recipes: bool = False,
):
    """Helper to create a minimal project structure."""
    path.mkdir(parents=True, exist_ok=True)

    # Create workflow/config.yaml with channel_names
    config_dir = path / "workflow"
    config_dir.mkdir(parents=True)
    channel_names = {}
    for c in range(1, cycles + 1):
        names = ["DAPI", "Blank_1a", "Blank_1b", "Blank_1c"]
        if c <= len(markers):
            # Replace blanks with markers for this cycle
            names = ["DAPI"] + markers[:3]
        channel_names[c] = names

    config = {"channel_names": channel_names}
    with open(config_dir / "config.yaml", "w") as f:
        yaml.dump(config, f)

    # Create registered images (tiny TIFs)
    registered = path / "data" / "processed" / "registered"
    for c in range(1, cycles + 1):
        cycle_dir = registered / f"cyc{c:02d}"
        cycle_dir.mkdir(parents=True)
        for name in channel_names[c]:
            img = np.random.randint(0, 500, (64, 64), dtype=np.uint16)
            tifffile.imwrite(str(cycle_dir / f"{name}.tif"), img)

    # Create recipes if requested
    if with_recipes:
        recipe_dir = path / "data" / "processed" / "Processing_parameters"
        recipe_dir.mkdir(parents=True)
        for marker in markers:
            (recipe_dir / f"{marker}_param.txt").write_text(
                f"signal_channel: '{marker}'\n"
                f"blankID: 'Blank_1a'\n"
                f"blank_clip_factor: 0\n"
                f"blank_scale_factor: 1.0\n"
                f"smooth_low: False\n"
                f"smooth_high: False\n"
                f"erosion: 0\n"
            )


# =============================================================================
# TestParseTissueType
# =============================================================================


class TestParseTissueType:
    def test_spleen_sp_prefix(self):
        assert parse_tissue_type("CX_19-001_SP_CC2-A28") == "spleen"

    def test_spleen_full_name(self):
        assert parse_tissue_type("CX_19-003_spleen_CC1-A") == "spleen"

    def test_lymph_node_ln_prefix(self):
        assert parse_tissue_type("CX_20-008_LN_n10_reg001") == "lymph_node"

    def test_lymph_node_full_name(self):
        assert parse_tissue_type("CX_19-002_lymph-node_R1") == "lymph_node"

    def test_thymus_th_prefix(self):
        assert parse_tissue_type("CX_20-005_TH_CC1-A") == "thymus"

    def test_thymus_full_name(self):
        assert parse_tissue_type("CX_19-001_thymus_CC1-A") == "thymus"

    def test_unknown_legacy_name(self):
        assert parse_tissue_type("1901CC2A") == "unknown"

    def test_unknown_with_experiment_json(self, tmp_path):
        """Falls back to experiment.json name field."""
        meta = tmp_path / "meta"
        meta.mkdir()
        (meta / "experiment.json").write_text(json.dumps({"name": "Spleen CC2 Sample A"}))
        assert parse_tissue_type("1901CC2A", tmp_path) == "spleen"

    def test_case_insensitive(self):
        assert parse_tissue_type("CX_19-003_SPLEEN_CC1-A") == "spleen"

    def test_src_prefix_spleen(self):
        assert parse_tissue_type("src_CX_19-004_spleen_CC2-2L") == "spleen"

    def test_src_prefix_ln(self):
        assert parse_tissue_type("src_CX_21-010_LN_n3") == "lymph_node"

    def test_src_prefix_sp(self):
        assert parse_tissue_type("src_CX_20-007_SP_CC1-A") == "spleen"


# =============================================================================
# TestDiscoverProjects
# =============================================================================


class TestDiscoverProjects:
    def test_finds_eligible_projects(self, projects_dir):
        projects = discover_projects(projects_dir)
        names = [p.name for p in projects]
        # Should find 3 projects (not CX_19-003_spleen_CC1-A which has manifest,
        # not src_CX_19-004_LN_R2 which has no config)
        assert "CX_19-001_SP_CC2-A28" in names
        assert "CX_19-002_lymph-node_R1" in names
        assert "CX_20-005_TH_CC1-A" in names

    def test_skips_completed_projects(self, projects_dir):
        projects = discover_projects(projects_dir)
        names = [p.name for p in projects]
        assert "CX_19-003_spleen_CC1-A" not in names

    def test_force_includes_completed(self, projects_dir):
        projects = discover_projects(projects_dir, force=True)
        names = [p.name for p in projects]
        assert "CX_19-003_spleen_CC1-A" in names

    def test_skips_no_config(self, projects_dir):
        projects = discover_projects(projects_dir)
        names = [p.name for p in projects]
        assert "src_CX_19-004_LN_R2" not in names

    def test_skips_non_directories(self, projects_dir):
        projects = discover_projects(projects_dir)
        names = [p.name for p in projects]
        assert "setup.log" not in names

    def test_dataset_filter(self, projects_dir):
        projects = discover_projects(projects_dir, dataset="CX_19-002_lymph-node_R1")
        assert len(projects) == 1
        assert projects[0].name == "CX_19-002_lymph-node_R1"

    def test_tissue_type_parsed(self, projects_dir):
        projects = discover_projects(projects_dir)
        by_name = {p.name: p for p in projects}
        assert by_name["CX_19-001_SP_CC2-A28"].tissue_type == "spleen"
        assert by_name["CX_19-002_lymph-node_R1"].tissue_type == "lymph_node"
        assert by_name["CX_20-005_TH_CC1-A"].tissue_type == "thymus"

    def test_detects_own_recipes(self, projects_dir):
        projects = discover_projects(projects_dir)
        by_name = {p.name: p for p in projects}
        assert by_name["CX_19-001_SP_CC2-A28"].has_own_recipes is True
        assert by_name["CX_19-002_lymph-node_R1"].has_own_recipes is False

    def test_counts_registered_channels(self, projects_dir):
        projects = discover_projects(projects_dir)
        for p in projects:
            assert p.n_registered_channels > 0


# =============================================================================
# TestLearnedParamsToRecipe
# =============================================================================


class TestLearnedParamsToRecipe:
    def test_basic_conversion(self):
        params = {
            "recommended_parameters": {
                "blank_name": "Blank_1b",
                "blank_clip_factor": 10,
                "blank_scale_factor": 0.8,
                "smooth_low": True,
                "low_size": 3,
                "erosion": 2,
            },
            "confidence": 0.85,
        }
        recipe = learned_params_to_recipe(params, "CD3e")
        assert recipe.marker_name == "CD3e"
        assert recipe.primary.blank_name == "Blank_1b"
        assert recipe.primary.blank_clip_factor == 10
        assert recipe.primary.blank_scale_factor == 0.8
        assert recipe.primary.smooth_low is True
        assert recipe.primary.erosion == 2
        assert recipe.second is None
        assert recipe.clean is None

    def test_with_second_subtraction(self):
        params = {
            "recommended_parameters": {
                "blank_name": "Blank_1a",
                "blank_scale_factor": 1.0,
                "second_blank": "DAPI",
                "second_clip_factor": 5,
                "second_scale_factor": 0.3,
            },
        }
        recipe = learned_params_to_recipe(params, "CD20")
        assert recipe.second is not None
        assert recipe.second.blank_name == "DAPI"
        assert recipe.second.blank_clip_factor == 5

    def test_empty_recommended_params(self):
        params = {"recommended_parameters": None, "confidence": 0.0}
        recipe = learned_params_to_recipe(params, "CD45")
        assert recipe.primary.blank_name == "Blank_1a"
        assert recipe.primary.blank_scale_factor == 1.0

    def test_defaults_for_missing_keys(self):
        params = {"recommended_parameters": {"blank_name": "Blank_1c"}}
        recipe = learned_params_to_recipe(params, "CD8")
        assert recipe.primary.blank_clip_factor == 0
        assert recipe.primary.smooth_low is False
        assert recipe.primary.erosion == 0


# =============================================================================
# TestResolveRecipes
# =============================================================================


class TestResolveRecipes:
    def test_tier1_own_recipes(self, projects_dir):
        """Project with own recipes should use them."""
        projects = discover_projects(projects_dir)
        p1 = next(p for p in projects if p.name == "CX_19-001_SP_CC2-A28")

        resolved = resolve_recipes_for_project(p1, projects_dir=projects_dir)

        # Markers with own recipes should be tier 1
        for marker in ["CD3e", "CD20", "CD45"]:
            if marker in resolved.sources:
                assert resolved.sources[marker] == "own_recipe"

    def test_tier3_template_recipes(self, projects_dir):
        """Project without own recipes falls to template when learned DB empty."""
        projects = discover_projects(projects_dir)
        p2 = next(p for p in projects if p.name == "CX_19-002_lymph-node_R1")

        # Create template recipes that match p2's markers
        template_dir = projects_dir / "template_recipes"
        template_dir.mkdir()
        (template_dir / "CD3e_param.txt").write_text(
            "signal_channel: 'CD3e'\nblankID: 'Blank_1a'\n"
            "blank_clip_factor: 0\nblank_scale_factor: 1.0\n"
        )

        # Mock learning DB to return no results so template is used
        mock_engine = MagicMock()
        mock_engine.recommend_parameters.return_value = {
            "status": "no_history",
            "confidence": 0.0,
            "recommended_parameters": None,
        }

        with patch(
            "kintsugi.claude.parameter_learning.ParameterLearningEngine",
            return_value=mock_engine,
        ):
            resolved = resolve_recipes_for_project(
                p2,
                template_recipe_dir=template_dir,
                projects_dir=projects_dir,
            )

        assert resolved.sources.get("CD3e") == "template"

    def test_tier4_auto_fallback(self, projects_dir):
        """Markers without any recipe source get 'auto'."""
        projects = discover_projects(projects_dir)
        p2 = next(p for p in projects if p.name == "CX_19-002_lymph-node_R1")

        resolved = resolve_recipes_for_project(p2, projects_dir=projects_dir)

        # Without template or learned DB, markers should fall to auto
        for _marker, source in resolved.sources.items():
            assert source in ("auto", "template", "learned")

    def test_tier2_learned_db(self, projects_dir):
        """Learned DB params override template when confidence is high."""
        projects = discover_projects(projects_dir)
        p2 = next(p for p in projects if p.name == "CX_19-002_lymph-node_R1")

        mock_engine = MagicMock()
        mock_engine.recommend_parameters.return_value = {
            "confidence": 0.9,
            "recommended_parameters": {
                "blank_name": "Blank_1a",
                "blank_scale_factor": 0.75,
            },
        }

        with patch(
            "kintsugi.claude.parameter_learning.ParameterLearningEngine",
            return_value=mock_engine,
        ):
            resolved = resolve_recipes_for_project(
                p2, projects_dir=projects_dir, min_confidence=0.6
            )

        learned_markers = [m for m, s in resolved.sources.items() if s == "learned"]
        assert len(learned_markers) > 0

    def test_summary_counts(self, projects_dir):
        """ResolvedRecipes.summary() returns per-source counts."""
        resolved = ResolvedRecipes(
            recipes={},
            sources={"CD3e": "own_recipe", "CD20": "template", "CD45": "auto"},
        )
        summary = resolved.summary()
        assert summary == {"own_recipe": 1, "template": 1, "auto": 1}


# =============================================================================
# TestWriteRecipesAsParamFiles
# =============================================================================


class TestWriteRecipesAsParamFiles:
    def test_roundtrip(self, tmp_path):
        """Written param files can be read back by load_legacy_recipes."""
        from kintsugi.signal.batch import load_legacy_recipes

        recipes = {
            "CD3e": MarkerRecipe(
                marker_name="CD3e",
                primary=SubtractionParams(
                    blank_name="Blank_1b",
                    blank_clip_factor=10,
                    blank_scale_factor=0.8,
                    smooth_low=True,
                    low_size=3,
                ),
            ),
            "CD20": MarkerRecipe(
                marker_name="CD20",
                primary=SubtractionParams(blank_name="Blank_1a"),
                second=SubtractionParams(
                    blank_name="DAPI",
                    blank_scale_factor=0.3,
                ),
                clean=CleanParams(
                    background_threshold=100,
                    remove_small=True,
                    small_size=50,
                ),
            ),
        }

        _write_recipes_as_param_files(recipes, tmp_path)

        # Read back
        loaded = load_legacy_recipes(tmp_path)
        assert "CD3e" in loaded
        assert "CD20" in loaded
        assert loaded["CD3e"].primary.blank_name == "Blank_1b"
        assert loaded["CD3e"].primary.blank_clip_factor == 10
        assert loaded["CD20"].second is not None
        assert loaded["CD20"].second.blank_name == "DAPI"
        assert loaded["CD20"].clean is not None
        assert loaded["CD20"].clean.background_threshold == 100

    def test_creates_files(self, tmp_path):
        recipes = {
            "CD45": MarkerRecipe(
                marker_name="CD45",
                primary=SubtractionParams(blank_name="Blank_1c"),
            )
        }
        _write_recipes_as_param_files(recipes, tmp_path)
        assert (tmp_path / "CD45_param.txt").exists()


# =============================================================================
# TestProjectInfo
# =============================================================================


class TestProjectInfo:
    def test_to_dict(self, tmp_path):
        test_path = tmp_path / "test"
        info = ProjectInfo(
            name="test",
            path=test_path,
            tissue_type="spleen",
            n_registered_channels=10,
        )
        d = info.to_dict()
        assert d["name"] == "test"
        assert d["path"] == str(test_path)
        assert d["tissue_type"] == "spleen"


# =============================================================================
# TestMultiProjectResult
# =============================================================================


class TestMultiProjectResult:
    def test_to_dict(self):
        result = MultiProjectResult(
            timestamp="2026-01-01",
            total_projects=3,
            completed=2,
            errors=1,
        )
        d = result.to_dict()
        assert d["total_projects"] == 3
        assert d["completed"] == 2
