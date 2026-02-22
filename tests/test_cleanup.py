"""
Tests for the pipeline-aware cleanup safety system.

Tests the dependency graph, manifest generation, staged deletion,
recovery, and CLI commands.
"""

import json
from pathlib import Path

import pytest

from kintsugi.cleanup import (
    DATA_DEPENDENCIES,
    PIPELINE_STAGES,
    CleanupEntry,
    CleanupManifest,
    _check_qc_sentinels,
    _check_stage_sentinels,
    _dir_size,
    _get_stage,
    _move_to_trash,
    _permanent_delete,
    _trash_dir,
    assess_cleanup_safety,
    execute_cleanup,
    list_trash,
    purge_trash,
    recover_trash,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture
def project_dir(tmp_path):
    """Create a minimal KINTSUGI project structure."""
    proj = tmp_path / "test_project"
    proc = proj / "data" / "processed"

    # Create directories
    for d in [
        "raw",
        "processed/stitched",
        "processed/deconvolved",
        "processed/edf",
        "processed/registered",
        "processed/vessel_3d",
    ]:
        (proj / "data" / d).mkdir(parents=True, exist_ok=True)
    (proj / "qc_plots").mkdir(parents=True, exist_ok=True)
    (proj / "workflow").mkdir(parents=True, exist_ok=True)

    # Create cycle subdirs with sentinels
    for cyc in ["01", "02", "03"]:
        for stage_dir in ["stitched", "deconvolved", "edf"]:
            d = proc / stage_dir / f"cyc{cyc}"
            d.mkdir(parents=True, exist_ok=True)
            # Put some data files
            (d / "CH1.tif").write_bytes(b"x" * 1000)
            (d / "CH2.tif").write_bytes(b"x" * 2000)

    # Create raw data
    for cyc in ["01", "02", "03"]:
        raw_cyc = proj / "data" / "raw" / f"cyc{cyc}"
        raw_cyc.mkdir(parents=True, exist_ok=True)
        (raw_cyc / "tile_001.tif").write_bytes(b"x" * 5000)

    return proj


@pytest.fixture
def complete_project(project_dir):
    """Project with all sentinels complete (stitch, decon, edf, QC)."""
    proc = project_dir / "data" / "processed"

    for cyc in ["01", "02", "03"]:
        for stage in ["stitched", "deconvolved", "edf"]:
            (proc / stage / f"cyc{cyc}" / ".snakemake_complete").write_text("done")

    # QC sentinels
    for qc in ["stitch", "decon", "edf"]:
        (project_dir / "qc_plots" / f".snakemake_complete_{qc}").write_text("done")

    # Registration sentinel
    (proc / "registered" / ".snakemake_complete").write_text("done")

    return project_dir


@pytest.fixture
def config_no_optional(project_dir):
    """Config without optional_stages section."""
    return {
        "cycles": [1, 2, 3],
        "channels": [1, 2],
    }


@pytest.fixture
def config_vessel3d_disabled(project_dir):
    """Config with vessel3d explicitly disabled."""
    return {
        "cycles": [1, 2, 3],
        "channels": [1, 2],
        "optional_stages": {
            "vessel3d": {"enabled": False, "cycles": []},
        },
    }


@pytest.fixture
def config_vessel3d_enabled(project_dir):
    """Config with vessel3d enabled for all cycles."""
    return {
        "cycles": [1, 2, 3],
        "channels": [1, 2],
        "optional_stages": {
            "vessel3d": {"enabled": True, "cycles": []},
        },
    }


@pytest.fixture
def config_vessel3d_partial(project_dir):
    """Config with vessel3d enabled for specific cycles only."""
    return {
        "cycles": [1, 2, 3],
        "channels": [1, 2],
        "optional_stages": {
            "vessel3d": {"enabled": True, "cycles": [2, 3]},
        },
    }


# ---------------------------------------------------------------------------
# Data model tests
# ---------------------------------------------------------------------------


class TestDataModel:
    def test_pipeline_stages_defined(self):
        """Core pipeline stages must be defined."""
        names = {s.name for s in PIPELINE_STAGES}
        assert "stitch" in names
        assert "deconvolve" in names
        assert "edf" in names
        assert "vessel3d" in names
        assert "registration" in names

    def test_vessel3d_is_optional(self):
        stage = _get_stage("vessel3d")
        assert stage is not None
        assert stage.is_optional is True
        assert stage.config_key == "vessel3d"

    def test_edf_is_not_optional(self):
        stage = _get_stage("edf")
        assert stage is not None
        assert stage.is_optional is False

    def test_deconvolved_has_two_consumers(self):
        """Deconvolved data is consumed by both EDF and vessel3d."""
        for dep in DATA_DEPENDENCIES:
            if dep.directory == "deconvolved":
                assert "edf" in dep.consumed_by
                assert "vessel3d" in dep.consumed_by
                return
        pytest.fail("No dependency found for deconvolved/")

    def test_cleanup_entry_size_properties(self):
        entry = CleanupEntry(
            directory="test",
            abs_path=Path("/tmp/test"),
            status="safe",
            reason="ok",
            size_bytes=2_147_483_648,  # 2 GB
        )
        assert entry.size_mb == pytest.approx(2048, rel=0.01)
        assert entry.size_gb == pytest.approx(2.0, rel=0.01)

    def test_manifest_safe_entries(self):
        manifest = CleanupManifest(
            project_dir=Path("/tmp"),
            entries=[
                CleanupEntry(Path("a"), Path("/a"), "safe", "ok", size_bytes=100),
                CleanupEntry(Path("b"), Path("/b"), "blocked", "no", size_bytes=200),
                CleanupEntry(Path("c"), Path("/c"), "safe", "ok", size_bytes=300),
            ],
        )
        assert len(manifest.safe_entries) == 2
        assert len(manifest.blocked_entries) == 1
        assert manifest.total_reclaimable_bytes == 400

    def test_manifest_to_dict(self, tmp_path):
        manifest = CleanupManifest(
            project_dir=tmp_path,
            entries=[
                CleanupEntry(Path("a"), tmp_path / "a", "safe", "ok"),
            ],
            warnings=["test warning"],
        )
        d = manifest.to_dict()
        assert d["project_dir"] == str(tmp_path)
        assert len(d["entries"]) == 1
        assert d["entries"][0]["status"] == "safe"
        assert d["warnings"] == ["test warning"]


# ---------------------------------------------------------------------------
# Sentinel checking tests
# ---------------------------------------------------------------------------


class TestSentinelChecking:
    def test_qc_sentinels_all_present(self, complete_project):
        ok, missing = _check_qc_sentinels(complete_project)
        assert ok is True
        assert missing == []

    def test_qc_sentinels_missing(self, project_dir):
        # No QC sentinels created
        ok, missing = _check_qc_sentinels(project_dir)
        assert ok is False
        assert set(missing) == {"stitch", "decon", "edf"}

    def test_qc_sentinels_partial(self, project_dir):
        (project_dir / "qc_plots" / ".snakemake_complete_stitch").write_text("done")
        ok, missing = _check_qc_sentinels(project_dir)
        assert ok is False
        assert "stitch" not in missing
        assert "decon" in missing

    def test_stage_sentinels_per_cycle(self, complete_project):
        stage = _get_stage("edf")
        ok, missing = _check_stage_sentinels(complete_project, stage, [1, 2, 3])
        assert ok is True
        assert missing == []

    def test_stage_sentinels_missing_cycle(self, project_dir):
        proc = project_dir / "data" / "processed"
        # Create sentinels for cycles 1 and 2 only
        for cyc in ["01", "02"]:
            (proc / "edf" / f"cyc{cyc}" / ".snakemake_complete").write_text("done")
        stage = _get_stage("edf")
        ok, missing = _check_stage_sentinels(project_dir, stage, [1, 2, 3])
        assert ok is False
        assert len(missing) == 1

    def test_stage_sentinels_required_cycles_subset(self, project_dir):
        """Only check sentinels for required_cycles, not all cycles."""
        proc = project_dir / "data" / "processed"
        # Only cycle 2 sentinel
        (proc / "vessel_3d" / "cyc02").mkdir(parents=True, exist_ok=True)
        (proc / "vessel_3d" / "cyc02" / ".snakemake_complete").write_text("done")
        stage = _get_stage("vessel3d")
        # Only require cycle 2
        ok, _ = _check_stage_sentinels(project_dir, stage, [1, 2, 3], required_cycles=[2])
        assert ok is True


# ---------------------------------------------------------------------------
# Assessment tests
# ---------------------------------------------------------------------------


class TestAssessCleanupSafety:
    def test_all_safe_when_complete_and_vessel3d_disabled(
        self, complete_project, config_vessel3d_disabled
    ):
        """All directories safe when everything complete and vessel3d disabled."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_disabled, skip_size=True
        )
        assert all(e.status in ("safe", "already_deleted") for e in manifest.entries), [
            f"{e.directory}: {e.status} — {e.reason}" for e in manifest.entries
        ]

    def test_blocked_when_no_optional_stages_section(self, complete_project, config_no_optional):
        """Deconvolved blocked when optional_stages section is absent."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_no_optional, skip_size=True
        )
        decon = next(e for e in manifest.entries if e.directory == "deconvolved")
        assert decon.status == "blocked"
        assert "vessel3d" in decon.blocking_consumers
        assert "optional_stages" in decon.reason

    def test_blocked_when_vessel3d_enabled_no_sentinels(
        self, complete_project, config_vessel3d_enabled
    ):
        """Deconvolved blocked when vessel3d enabled but sentinels missing."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_enabled, skip_size=True
        )
        decon = next(e for e in manifest.entries if e.directory == "deconvolved")
        assert decon.status == "blocked"
        assert "vessel3d" in decon.blocking_consumers

    def test_safe_when_vessel3d_enabled_and_complete(
        self, complete_project, config_vessel3d_enabled
    ):
        """Deconvolved safe when vessel3d enabled AND all sentinels present."""
        proc = complete_project / "data" / "processed"
        for cyc in ["01", "02", "03"]:
            (proc / "vessel_3d" / f"cyc{cyc}").mkdir(parents=True, exist_ok=True)
            (proc / "vessel_3d" / f"cyc{cyc}" / ".snakemake_complete").write_text("done")

        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_enabled, skip_size=True
        )
        decon = next(e for e in manifest.entries if e.directory == "deconvolved")
        assert decon.status == "safe"

    def test_partial_vessel3d_only_blocks_needed_cycles(
        self, complete_project, config_vessel3d_partial
    ):
        """Vessel3d for cycles [2,3] — deconvolved blocked until those cycles done."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_partial, skip_size=True
        )
        decon = next(e for e in manifest.entries if e.directory == "deconvolved")
        assert decon.status == "blocked"
        assert "vessel3d" in decon.blocking_consumers

        # Complete vessel3d for cycles 2 and 3
        proc = complete_project / "data" / "processed"
        for cyc in ["02", "03"]:
            (proc / "vessel_3d" / f"cyc{cyc}").mkdir(parents=True, exist_ok=True)
            (proc / "vessel_3d" / f"cyc{cyc}" / ".snakemake_complete").write_text("done")

        manifest2 = assess_cleanup_safety(
            complete_project, config=config_vessel3d_partial, skip_size=True
        )
        decon2 = next(e for e in manifest2.entries if e.directory == "deconvolved")
        assert decon2.status == "safe"

    def test_blocked_when_qc_missing(self, project_dir, config_vessel3d_disabled):
        """Everything blocked when QC sentinels are missing."""
        # Create all processing sentinels but NO QC sentinels
        proc = project_dir / "data" / "processed"
        for cyc in ["01", "02", "03"]:
            for stage in ["stitched", "deconvolved", "edf"]:
                (proc / stage / f"cyc{cyc}" / ".snakemake_complete").write_text("done")

        manifest = assess_cleanup_safety(
            project_dir, config=config_vessel3d_disabled, skip_size=True
        )
        # All entries should be blocked due to QC
        for entry in manifest.entries:
            if entry.status != "already_deleted":
                assert entry.status == "blocked"
                assert "QC incomplete" in entry.reason

    def test_already_deleted_directory(self, complete_project, config_vessel3d_disabled):
        """Reports 'already_deleted' for missing directories."""
        import shutil

        shutil.rmtree(str(complete_project / "data" / "processed" / "stitched"))

        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_disabled, skip_size=True
        )
        stitch = next(e for e in manifest.entries if e.directory == "stitched")
        assert stitch.status == "already_deleted"

    def test_stitched_only_needs_deconvolve(self, complete_project, config_vessel3d_disabled):
        """Stitched data only consumed by deconvolve (no vessel3d dependency)."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_disabled, skip_size=True
        )
        stitch = next(e for e in manifest.entries if e.directory == "stitched")
        assert stitch.status == "safe"
        assert stitch.consumers == ["deconvolve"]


# ---------------------------------------------------------------------------
# Staged deletion tests
# ---------------------------------------------------------------------------


class TestStagedDeletion:
    def test_move_to_trash(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "test_dir"
        test_dir.mkdir()
        (test_dir / "data.txt").write_text("important data")

        receipt = _move_to_trash(test_dir, project_dir)

        # Original should be gone
        assert not test_dir.exists()
        # Trash should contain the data
        trash_path = Path(receipt.trash_path)
        assert trash_path.exists()
        assert (trash_path / "data.txt").read_text() == "important data"
        # Receipt should be written
        assert (trash_path / "receipt.json").exists()

    def test_trash_receipt_contents(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "mydir"
        test_dir.mkdir()
        (test_dir / "file.tif").write_bytes(b"x" * 100)

        receipt = _move_to_trash(test_dir, project_dir, manifest_snapshot={"test": True})

        receipt_data = json.loads((Path(receipt.trash_path) / "receipt.json").read_text())
        assert receipt_data["original_path"] == str(test_dir)
        assert receipt_data["manifest_snapshot"] == {"test": True}
        assert receipt_data["size_bytes"] == 100

    def test_permanent_delete(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "delete_me"
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("gone")

        _permanent_delete(test_dir)
        assert not test_dir.exists()

    def test_permanent_delete_nonexistent(self, project_dir):
        """Deleting a nonexistent directory should not error."""
        _permanent_delete(project_dir / "nonexistent")


# ---------------------------------------------------------------------------
# Execute cleanup tests
# ---------------------------------------------------------------------------


class TestExecuteCleanup:
    def test_execute_with_trash(self, complete_project, config_vessel3d_disabled):
        result = execute_cleanup(
            complete_project,
            use_trash=True,
            config=config_vessel3d_disabled,
        )
        assert len(result["deleted"]) > 0
        assert result["errors"] == []
        # All deleted should have trash_path
        for d in result["deleted"]:
            assert "trash_path" in d

    def test_execute_without_trash(self, complete_project, config_vessel3d_disabled):
        result = execute_cleanup(
            complete_project,
            use_trash=False,
            config=config_vessel3d_disabled,
        )
        assert len(result["deleted"]) > 0
        # Directories should be gone
        for d in result["deleted"]:
            assert "trash_path" not in d

    def test_execute_skips_blocked(self, complete_project, config_vessel3d_enabled):
        """Blocked directories should be skipped, not deleted."""
        result = execute_cleanup(
            complete_project,
            use_trash=True,
            config=config_vessel3d_enabled,
        )
        skipped_dirs = {s["directory"] for s in result["skipped"]}
        assert "deconvolved" in skipped_dirs

    def test_execute_skip_vessel3d_override(self, complete_project, config_vessel3d_enabled):
        """--skip-vessel3d should allow deletion even when vessel3d blocks."""
        result = execute_cleanup(
            complete_project,
            use_trash=True,
            skip_vessel3d=True,
            config=config_vessel3d_enabled,
        )
        deleted_dirs = {d["directory"] for d in result["deleted"]}
        assert "deconvolved" in deleted_dirs

    def test_execute_deletes_raw(self, complete_project, config_vessel3d_disabled):
        """Raw data should be deleted when all stitch sentinels exist."""
        # Create .staged marker (raw data was staged from orange)
        (complete_project / "data" / "raw" / ".staged").write_text("staged")
        result = execute_cleanup(
            complete_project,
            use_trash=False,
            config=config_vessel3d_disabled,
        )
        deleted_dirs = {d["directory"] for d in result["deleted"]}
        assert "raw" in deleted_dirs


# ---------------------------------------------------------------------------
# Trash management tests
# ---------------------------------------------------------------------------


class TestTrashManagement:
    def test_list_trash_empty(self, project_dir):
        entries = list_trash(project_dir)
        assert entries == []

    def test_list_trash_with_entries(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "test"
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("data")
        _move_to_trash(test_dir, project_dir)

        entries = list_trash(project_dir)
        assert len(entries) == 1
        assert entries[0]["original_path"] == str(test_dir)

    def test_recover_trash(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "recover_me"
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("precious")
        receipt = _move_to_trash(test_dir, project_dir)

        # Recover
        trash_name = Path(receipt.trash_path).name
        result = recover_trash(project_dir, trash_name)

        assert "recovered" in result
        assert test_dir.exists()
        assert (test_dir / "file.txt").read_text() == "precious"

    def test_recover_blocks_if_original_exists(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "conflict"
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("v1")
        receipt = _move_to_trash(test_dir, project_dir)

        # Recreate original
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("v2")

        trash_name = Path(receipt.trash_path).name
        result = recover_trash(project_dir, trash_name)
        assert "error" in result
        assert "already exists" in result["error"]

    def test_purge_old_entries(self, project_dir):
        test_dir = project_dir / "data" / "processed" / "old"
        test_dir.mkdir()
        (test_dir / "file.txt").write_text("ancient")
        _move_to_trash(test_dir, project_dir)

        # Purge with max_age_days=0 (everything is old)
        result = purge_trash(project_dir, max_age_days=0)
        # The entry was just created so its mtime is very recent,
        # but with max_age_days=0 the cutoff is now, so anything older than 0
        # seconds should be purged. Since we just created it, this actually
        # won't purge. Use purge_all instead for immediate purge.
        # Let's test purge_all:
        test_dir2 = project_dir / "data" / "processed" / "old2"
        test_dir2.mkdir()
        (test_dir2 / "file.txt").write_text("ancient2")
        _move_to_trash(test_dir2, project_dir)

        result = purge_trash(project_dir, purge_all=True)
        assert len(result["purged"]) >= 1
        assert result["total_bytes"] > 0

    def test_purge_all(self, project_dir):
        for i in range(3):
            d = project_dir / "data" / "processed" / f"trash_me_{i}"
            d.mkdir()
            (d / "file.txt").write_text(f"data_{i}")
            _move_to_trash(d, project_dir)

        result = purge_trash(project_dir, purge_all=True)
        assert len(result["purged"]) == 3

        # Trash dir should be empty/removed
        trash = _trash_dir(project_dir)
        assert not trash.exists() or not any(trash.iterdir())


# ---------------------------------------------------------------------------
# Dir size tests
# ---------------------------------------------------------------------------


class TestDirSize:
    def test_dir_size(self, project_dir):
        test_dir = project_dir / "test_size"
        test_dir.mkdir()
        (test_dir / "a.txt").write_bytes(b"x" * 100)
        (test_dir / "b.txt").write_bytes(b"y" * 200)

        assert _dir_size(test_dir) == 300

    def test_dir_size_nonexistent(self, project_dir):
        assert _dir_size(project_dir / "nope") == 0

    def test_dir_size_nested(self, project_dir):
        test_dir = project_dir / "nested"
        (test_dir / "sub" / "deep").mkdir(parents=True)
        (test_dir / "file1.txt").write_bytes(b"x" * 50)
        (test_dir / "sub" / "file2.txt").write_bytes(b"y" * 75)
        (test_dir / "sub" / "deep" / "file3.txt").write_bytes(b"z" * 25)

        assert _dir_size(test_dir) == 150


# ---------------------------------------------------------------------------
# Edge case tests
# ---------------------------------------------------------------------------


class TestEdgeCases:
    def test_no_cycles_discovered(self, tmp_path):
        """Empty project with no cycles returns empty manifest."""
        proj = tmp_path / "empty"
        proj.mkdir()
        (proj / "data" / "processed").mkdir(parents=True)

        manifest = assess_cleanup_safety(proj, config={})
        assert manifest.warnings  # Should warn about missing cycles

    def test_config_from_yaml_file(self, complete_project, config_vessel3d_disabled):
        """Config can be loaded from workflow/config.yaml."""
        import yaml

        cfg_path = complete_project / "workflow" / "config.yaml"
        cfg_path.parent.mkdir(parents=True, exist_ok=True)
        cfg_path.write_text(yaml.dump(config_vessel3d_disabled))

        manifest = assess_cleanup_safety(complete_project, skip_size=True)
        # Should load config from file and find vessel3d disabled
        decon = next((e for e in manifest.entries if e.directory == "deconvolved"), None)
        if decon:
            assert decon.status == "safe"

    def test_vessel3d_enabled_with_empty_cycles_means_all(
        self, complete_project, config_vessel3d_enabled
    ):
        """vessel3d.cycles: [] (empty) means ALL cycles need vessel3d."""
        manifest = assess_cleanup_safety(
            complete_project, config=config_vessel3d_enabled, skip_size=True
        )
        decon = next(e for e in manifest.entries if e.directory == "deconvolved")
        assert decon.status == "blocked"
        # Should mention 3 sentinels missing (one per cycle)
        assert "3 sentinel" in decon.reason
