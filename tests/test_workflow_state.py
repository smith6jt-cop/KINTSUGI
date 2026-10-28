"""Tests for the persistent per-project workflow state (``kintsugi.assist.workflow_state``)."""

from __future__ import annotations

import json
import sqlite3
import subprocess
import sys

import numpy as np
import pytest

from kintsugi.assist.workflow_state import (
    LOCK_FILE,
    STATE_FILE,
    STATE_VERSION,
    WorkflowState,
)


def _tree(path):
    return sorted(str(p.relative_to(path)) for p in path.rglob("*"))


def test_construction_and_reads_touch_no_disk(tmp_path):
    state = WorkflowState(tmp_path)
    assert state.channels() == {}
    assert state.get_channel("x") is None
    assert state.summary()["total_channels"] == 0 and state.summary()["session_id"] is None
    assert state.get_review_queue() == [] and state.get_processing_queue() == []
    state.ensure_session()
    assert _tree(tmp_path) == []


def test_first_write_creates_state_and_audit_log(tmp_path):
    state = WorkflowState(tmp_path)
    state.add_channel("cyc001_CD3", cycle="cyc001")

    data = json.loads((tmp_path / ".kintsugi" / STATE_FILE).read_text())
    assert data["version"] == STATE_VERSION
    assert data["channels"]["cyc001_CD3"]["status"] == "pending"
    assert data["current_session"] and len(data["sessions"]) == 1
    rows = sqlite3.connect(tmp_path / ".kintsugi" / "workflow_audit.db").execute(
        "SELECT action, channel FROM audit_log"
    )
    assert list(rows) == [("channel_added", "cyc001_CD3")]
    # The write is atomic: no temporary file is left behind (the lock file stays;
    # only the OS lock on it is released).
    assert sorted(p.name for p in (tmp_path / ".kintsugi").iterdir()) == sorted(
        [STATE_FILE, "workflow_audit.db", LOCK_FILE]
    )


def test_add_channel_is_idempotent(tmp_path):
    state = WorkflowState(tmp_path)
    state.add_channel("ch", cycle="cyc001")
    state.record("ch", "denoise", {"filter_size": 3})
    again = state.add_channel("ch", cycle="other")
    assert again.cycle == "cyc001" and again.status == "processing"
    assert len(again.processing_steps) == 1


def test_record_registers_unknown_channels_and_tracks_quality(tmp_path):
    state = WorkflowState(tmp_path)
    ch = state.record(
        "ch",
        "subtract_blank",
        {"scale": np.float32(0.5), "ranges": np.array([1, 2])},
        image="ch_sub",
        cycle="cyc002",
        quality_after=0.7,
    )
    assert ch.cycle == "cyc002" and ch.status == "processing" and ch.quality_score == 0.7
    step = state.get_channel("ch").processing_steps[0]
    assert step["params"] == {"scale": 0.5, "ranges": [1, 2]}
    assert step["image"] == "ch_sub" and step["session_id"] == state.session_id


def test_decisions_queues_and_unknown_channels(tmp_path):
    state = WorkflowState(tmp_path)
    for name in ("a", "b", "c"):
        state.add_channel(name)
    state.update_channel("b", status="review")
    assert state.set_user_decision(
        "c", "approved", "looks good", final_params={"denoise": {}}, quality_score=0.9
    )

    assert state.get_processing_queue() == ["a"]
    assert state.get_review_queue() == ["b"]
    summary = state.summary()
    assert summary["approved"] == ["c"]
    assert summary["status_breakdown"] == {"pending": 1, "review": 1, "approved": 1}
    approved = state.get_channel("c")
    assert approved.review_notes == "looks good" and approved.quality_score == 0.9

    assert state.update_channel("nope", status="review") is False
    assert state.set_user_decision("nope", "rejected") is False
    assert state.add_processing_step("nope", "denoise", {}) is False
    with pytest.raises(ValueError, match="decision must be one of"):
        state.set_user_decision("a", "maybe")
    with pytest.raises(ValueError, match="unknown ChannelState fields"):
        state.update_channel("a", colour="red")


def test_a_new_instance_resumes_the_session_and_sees_everything(tmp_path):
    first = WorkflowState(tmp_path)
    first.record("ch", "denoise", {})
    second = WorkflowState(tmp_path)
    assert second.ensure_session().session_id == first.session_id
    assert second.get_channel("ch").status == "processing"

    started = WorkflowState(tmp_path).new_session()
    assert started.session_id != first.session_id
    assert json.loads((tmp_path / ".kintsugi" / STATE_FILE).read_text())["current_session"] == (
        started.session_id
    )


def test_two_instances_interleaving_writes_lose_nothing(tmp_path):
    a, b = WorkflowState(tmp_path), WorkflowState(tmp_path)
    a.add_channel("from_a")
    b.add_channel("from_b")
    a.record("from_a", "denoise", {})
    b.set_user_decision("from_b", "rejected")
    channels = WorkflowState(tmp_path).channels()
    assert set(channels) == {"from_a", "from_b"}
    assert channels["from_a"].status == "processing" and channels["from_b"].status == "rejected"


@pytest.mark.skipif(sys.platform == "win32", reason="uses fcntl in the helper process")
def test_a_lock_held_by_another_process_blocks_until_it_exits(tmp_path, monkeypatch):
    from kintsugi.assist import workflow_state as ws

    state = WorkflowState(tmp_path)
    state.add_channel("ch")
    lock = tmp_path / ".kintsugi" / LOCK_FILE
    holder = subprocess.Popen(
        [
            sys.executable,
            "-c",
            "import fcntl, os, sys, time\n"
            f"fd = os.open({str(lock)!r}, os.O_RDWR)\n"
            "fcntl.lockf(fd, fcntl.LOCK_EX)\n"
            "print('locked', flush=True)\n"
            "sys.stdin.read()\n",
        ],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        text=True,
    )
    try:
        assert holder.stdout.readline().strip() == "locked"
        monkeypatch.setattr(ws, "LOCK_TIMEOUT_S", 0.3)
        with pytest.raises(TimeoutError, match="locked by another process"):
            state.record("ch", "denoise", {})
    finally:
        holder.stdin.close()  # the holder exits; the kernel drops its lock
        holder.wait(timeout=10)
    state.record("ch", "denoise", {})  # no stale lock to clean up
    assert lock.exists() and len(state.get_channel("ch").processing_steps) == 1


def test_quality_before_is_kept(tmp_path):
    state = WorkflowState(tmp_path)
    state.add_channel("ch")
    assert state.add_processing_step("ch", "denoise", {}, quality_before=0.3, quality_after=0.8)
    step = state.get_channel("ch").processing_steps[0]
    assert step["quality_before"] == 0.3 and step["quality_after"] == 0.8


def test_an_unreadable_state_file_is_moved_aside(tmp_path, caplog):
    state_dir = tmp_path / ".kintsugi"
    state_dir.mkdir()
    (state_dir / STATE_FILE).write_text("{not json")
    state = WorkflowState(tmp_path)
    assert state.channels() == {}  # reads tolerate it and change nothing
    assert (state_dir / STATE_FILE).read_text() == "{not json"

    state.add_channel("ch")
    assert list(state.channels()) == ["ch"]
    assert [p.name for p in state_dir.glob(f"{STATE_FILE}.corrupt-*")]


def test_legacy_directory_is_read_then_moved_on_write(tmp_path):
    WorkflowState(tmp_path).add_channel("old")
    (tmp_path / ".kintsugi").rename(tmp_path / ".claude_workflow")

    state = WorkflowState(tmp_path)
    assert list(state.channels()) == ["old"]
    assert not (tmp_path / ".kintsugi").exists()
    state.add_channel("new")
    assert not (tmp_path / ".claude_workflow").exists()
    assert set(state.channels()) == {"old", "new"}


def test_export_parameters(tmp_path):
    state = WorkflowState(tmp_path)
    state.record("ch", "denoise", {"filter_size": 3})
    state.set_user_decision("ch", "approved", final_params={"denoise": {"filter_size": 3}})
    out = json.loads(state.export_parameters(tmp_path / "out.json").read_text())
    assert out["channels"]["ch"]["final_params"] == {"denoise": {"filter_size": 3}}


def test_concurrent_processes_lose_no_updates(tmp_path):
    """Four processes register 15 channels each at the same time; all 60 survive."""
    script = (
        "import sys\n"
        "from kintsugi.assist.workflow_state import WorkflowState\n"
        "state = WorkflowState(sys.argv[1])\n"
        "for i in range(15):\n"
        "    state.record(f'{sys.argv[2]}_{i}', 'denoise', {'i': i})\n"
    )
    procs = [
        subprocess.Popen([sys.executable, "-c", script, str(tmp_path), f"p{n}"]) for n in range(4)
    ]
    assert [p.wait(timeout=120) for p in procs] == [0, 0, 0, 0]
    channels = WorkflowState(tmp_path).channels()
    assert len(channels) == 60
    assert all(len(ch.processing_steps) == 1 for ch in channels.values())
