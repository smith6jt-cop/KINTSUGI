"""Tests for the project state directory (``<project>/.kintsugi``) and its migration."""

import logging
from pathlib import Path

import pytest

from kintsugi.assist import state_dir
from kintsugi.assist.state_dir import (
    LEGACY_STATE_DIR_NAME,
    STATE_DIR_NAME,
    project_state_dir,
)


@pytest.fixture(autouse=True)
def _reset_warnings():
    state_dir._warned.clear()
    yield
    state_dir._warned.clear()


def _legacy(project: Path) -> Path:
    legacy = project / LEGACY_STATE_DIR_NAME
    legacy.mkdir(parents=True)
    (legacy / "parameter_learning.db").write_bytes(b"legacy")
    return legacy


def test_read_on_a_fresh_project_creates_nothing(tmp_path):
    assert project_state_dir(tmp_path) == tmp_path / STATE_DIR_NAME
    assert list(tmp_path.iterdir()) == []


def test_read_uses_the_legacy_directory_when_only_it_exists(tmp_path):
    legacy = _legacy(tmp_path)
    assert project_state_dir(tmp_path) == legacy
    assert not (tmp_path / STATE_DIR_NAME).exists()


def test_create_on_a_fresh_project_makes_the_directory(tmp_path):
    current = project_state_dir(tmp_path, create=True)
    assert current == tmp_path / STATE_DIR_NAME and current.is_dir()


def test_first_write_moves_the_legacy_directory(tmp_path):
    _legacy(tmp_path)
    current = project_state_dir(tmp_path, create=True)
    assert current == tmp_path / STATE_DIR_NAME
    assert (current / "parameter_learning.db").read_bytes() == b"legacy"
    assert not (tmp_path / LEGACY_STATE_DIR_NAME).exists()
    assert project_state_dir(tmp_path) == current


def test_both_directories_use_the_new_one_and_warn_once(tmp_path, caplog):
    legacy = _legacy(tmp_path)
    (tmp_path / STATE_DIR_NAME).mkdir()
    with caplog.at_level(logging.WARNING, logger="kintsugi.assist.state_dir"):
        assert project_state_dir(tmp_path) == tmp_path / STATE_DIR_NAME
        assert project_state_dir(tmp_path, create=True) == tmp_path / STATE_DIR_NAME
        project_state_dir(tmp_path, create=True)
    assert len([r for r in caplog.records if "both exist" in r.getMessage()]) == 1
    assert (legacy / "parameter_learning.db").read_bytes() == b"legacy"


def test_a_failed_move_keeps_writing_to_the_legacy_directory(tmp_path, monkeypatch, caplog):
    legacy = _legacy(tmp_path)

    def refuse(self, target):
        raise PermissionError("read-only")

    monkeypatch.setattr(Path, "rename", refuse)
    with caplog.at_level(logging.WARNING, logger="kintsugi.assist.state_dir"):
        assert project_state_dir(tmp_path, create=True) == legacy
    assert not (tmp_path / STATE_DIR_NAME).exists()
    assert "Could not move" in caplog.text


def test_a_move_done_by_another_process_counts_as_success(tmp_path, monkeypatch):
    _legacy(tmp_path)
    real_rename = Path.rename

    def raced(self, target):
        real_rename(self, target)  # the other process wins ...
        raise FileNotFoundError(str(self))  # ... and our rename finds nothing to move

    monkeypatch.setattr(Path, "rename", raced)
    current = project_state_dir(tmp_path, create=True)
    assert current == tmp_path / STATE_DIR_NAME
    assert (current / "parameter_learning.db").exists()


def test_reads_alone_warn_once_when_both_directories_exist(tmp_path, caplog):
    _legacy(tmp_path)
    (tmp_path / STATE_DIR_NAME).mkdir()
    with caplog.at_level(logging.WARNING, logger="kintsugi.assist.state_dir"):
        for _ in range(3):
            assert project_state_dir(tmp_path) == tmp_path / STATE_DIR_NAME
    assert len([r for r in caplog.records if "both exist" in r.getMessage()]) == 1


def test_workflow_state_follows_a_later_migration(tmp_path, monkeypatch):
    """A failed move keeps the legacy directory; a later successful move is followed."""
    from kintsugi.assist.workflow_state import WorkflowState

    legacy = _legacy(tmp_path)
    real_rename = Path.rename

    def refuse(self, target):
        raise PermissionError("busy")

    monkeypatch.setattr(Path, "rename", refuse)
    state = WorkflowState(tmp_path)
    assert state.state_dir == legacy

    monkeypatch.setattr(Path, "rename", real_rename)
    project_state_dir(tmp_path, create=True)  # e.g. a learning write moves it
    state.new_session()
    assert state.json_path == tmp_path / STATE_DIR_NAME / "workflow_state.json"
    assert state.json_path.exists()
    assert not legacy.exists()
