"""Persistent per-project workflow state: channel status, processing steps, review decisions.

State lives in ``<project>/.kintsugi/`` (see :mod:`kintsugi.assist.state_dir`):

* ``workflow_state.json`` - every channel's status, processing steps (with the session
  that recorded them), final parameters and review notes, plus the list of sessions;
* ``workflow_audit.db`` - an append-only SQLite ``audit_log`` of every change.

Channels belong to the project, not to a session, so a channel approved yesterday is
still approved after the MCP server (or notebook) restarts; each recorded step carries
the id of the session that recorded it. Nothing is written until the first change:
constructing a :class:`WorkflowState` or reading from it never touches the disk.

Every change is a read-modify-write of the JSON file under a lock file
(``workflow_state.lock``), written atomically (temporary file + ``os.replace``), so two
processes working on one project (two MCP clients, or a notebook and a client) merge
their changes at channel granularity instead of overwriting each other.

Example::

    state = WorkflowState("/data/projects/tonsil")
    state.add_channel("cyc001_CD3", cycle="cyc001")
    state.record("cyc001_CD3", "subtract_blank", {"blank_scale_factor": 0.8})
    state.set_user_decision("cyc001_CD3", "approved", quality_score=0.92)
    state.get_review_queue()
"""

from __future__ import annotations

import json
import logging
import os
import secrets
import sqlite3
import tempfile
import threading
import time
from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any, TypeVar

from kintsugi.assist._coerce import json_default, plain
from kintsugi.assist.state_dir import project_state_dir

logger = logging.getLogger(__name__)

STATE_VERSION = 2
STATE_FILE = "workflow_state.json"
AUDIT_DB = "workflow_audit.db"
LOCK_FILE = "workflow_state.lock"

STATUSES = ("pending", "processing", "review", "approved", "rejected", "needs_adjustment")
DECISIONS = ("approved", "rejected", "needs_adjustment", "pending", "review")

LOCK_TIMEOUT_S = 20.0

try:  # POSIX record locks (also honoured over NFS, where HPC projects live)
    import fcntl
except ImportError:  # pragma: no cover - Windows
    fcntl = None  # type: ignore[assignment]
    import msvcrt

# POSIX record locks belong to the process, so writers in one process (threads,
# several WorkflowState instances) are serialised by a thread lock per directory.
_thread_locks: dict[Path, threading.Lock] = {}
_thread_locks_guard = threading.Lock()


def _thread_lock(state_dir: Path) -> threading.Lock:
    with _thread_locks_guard:
        return _thread_locks.setdefault(state_dir.resolve(), threading.Lock())


def _try_lock(fd: int) -> bool:
    try:
        if fcntl is not None:
            fcntl.lockf(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        else:  # pragma: no cover - Windows
            os.lseek(fd, 0, os.SEEK_SET)
            msvcrt.locking(fd, msvcrt.LK_NBLCK, 1)
        return True
    except OSError:
        return False


def _unlock(fd: int) -> None:
    if fcntl is not None:
        fcntl.lockf(fd, fcntl.LOCK_UN)
    else:  # pragma: no cover - Windows
        os.lseek(fd, 0, os.SEEK_SET)
        msvcrt.locking(fd, msvcrt.LK_UNLCK, 1)


T = TypeVar("T")


def _now() -> str:
    return datetime.now().isoformat()


@dataclass
class ChannelState:
    """State of one channel in a project."""

    channel_name: str
    cycle: str = "unknown"
    status: str = "pending"  # one of STATUSES
    quality_score: float | None = None
    processing_steps: list[dict[str, Any]] = field(default_factory=list)
    final_params: dict[str, Any] = field(default_factory=dict)
    review_notes: str = ""
    created: str = field(default_factory=_now)
    modified: str = field(default_factory=_now)

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ChannelState:
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class WorkflowSession:
    """One working session (a server run, a notebook run) on a project."""

    session_id: str
    project_path: str
    created: str = field(default_factory=_now)


def _new_session_id() -> str:
    return f"{datetime.now():%Y%m%d_%H%M%S}_{secrets.token_hex(3)}"


def _empty_state() -> dict[str, Any]:
    return {"version": STATE_VERSION, "channels": {}, "sessions": [], "current_session": None}


class WorkflowState:
    """Per-project workflow state backed by ``.kintsugi/workflow_state.json``.

    Args:
        project_path: The KINTSUGI project directory.
        use_sqlite: Also append every change to the ``audit_log`` table of
            ``workflow_audit.db``.
    """

    def __init__(self, project_path: str | Path, use_sqlite: bool = True):
        self.project_path = Path(project_path).resolve()
        self.use_sqlite = use_sqlite
        self._session_id: str | None = None

    # ------------------------------------------------------------------ paths
    # Resolved on every access: reads never create the directory, and an instance
    # follows a migration of a legacy ``.claude_workflow/`` done by another writer.

    @property
    def state_dir(self) -> Path:
        return project_state_dir(self.project_path)

    @property
    def json_path(self) -> Path:
        return self.state_dir / STATE_FILE

    @property
    def db_path(self) -> Path | None:
        return self.state_dir / AUDIT_DB if self.use_sqlite else None

    # ------------------------------------------------------------- persistence

    def _read(self, path: Path | None = None) -> dict[str, Any]:
        """Load the state file; a missing or unreadable file reads as an empty state.

        Reads never modify anything; the next write moves an unreadable file aside
        (see :meth:`_quarantine_if_unreadable`).
        """
        path = path or self.json_path
        if not path.exists():
            return _empty_state()
        try:
            data = json.loads(path.read_text())
            if not isinstance(data, dict) or data.get("version") != STATE_VERSION:
                raise ValueError(f"unsupported workflow state format in {path}")
            data.setdefault("channels", {})
            data.setdefault("sessions", [])
            data.setdefault("current_session", None)
            return data
        except (OSError, ValueError) as exc:
            logger.warning("Ignoring unreadable workflow state %s: %s", path, exc)
            return _empty_state()

    def _quarantine_if_unreadable(self, path: Path) -> None:
        if not path.exists():
            return
        try:
            data = json.loads(path.read_text())
            if isinstance(data, dict) and data.get("version") == STATE_VERSION:
                return
        except (OSError, ValueError):
            pass
        target = path.with_name(f"{path.name}.corrupt-{datetime.now():%Y%m%d_%H%M%S}")
        try:
            path.replace(target)
            logger.warning("Moved unreadable workflow state to %s; starting afresh", target)
        except OSError as exc:  # pragma: no cover - permissions
            logger.warning("Could not move unreadable workflow state %s: %s", path, exc)

    @contextmanager
    def _locked(self) -> Iterator[Path]:
        """Hold an exclusive OS lock on ``workflow_state.lock`` in the state directory.

        The lock is an advisory record lock (``fcntl.lockf``; ``msvcrt.locking`` on
        Windows) on a lock file that is never deleted: the kernel releases it when
        the holder closes it or exits, so a crashed writer cannot leave a stale lock
        and no writer can remove another's.
        """
        state_dir = project_state_dir(self.project_path, create=True)
        lock = state_dir / LOCK_FILE
        deadline = time.monotonic() + LOCK_TIMEOUT_S
        thread_lock = _thread_lock(state_dir)
        if not thread_lock.acquire(timeout=LOCK_TIMEOUT_S):
            raise TimeoutError(f"workflow state is locked by another thread: {lock}")
        try:
            fd = os.open(lock, os.O_CREAT | os.O_RDWR, 0o644)
            try:
                while not _try_lock(fd):
                    if time.monotonic() > deadline:
                        raise TimeoutError(f"workflow state is locked by another process: {lock}")
                    time.sleep(0.05)
                try:
                    yield state_dir
                finally:
                    _unlock(fd)
            finally:
                os.close(fd)
        finally:
            thread_lock.release()

    def _write(self, state_dir: Path, data: dict[str, Any]) -> None:
        path = state_dir / STATE_FILE
        fd, tmp = tempfile.mkstemp(dir=state_dir, prefix=".workflow_state.", suffix=".tmp")
        try:
            with os.fdopen(fd, "w") as f:
                json.dump(data, f, indent=2, default=json_default)
                f.write("\n")
            os.replace(tmp, path)
        except BaseException:
            Path(tmp).unlink(missing_ok=True)
            raise

    def _mutate(
        self,
        action: str,
        change: Callable[[dict[str, Any]], T],
        *,
        channel: str | None = None,
        details: dict[str, Any] | None = None,
        decision: str | None = None,
    ) -> T:
        """Apply ``change`` to the freshly read state under the lock and persist it."""
        with self._locked() as state_dir:
            self._quarantine_if_unreadable(state_dir / STATE_FILE)
            data = self._read(state_dir / STATE_FILE)
            session_id = self._ensure_session_in(data)
            result = change(data)
            self._write(state_dir, data)
            self._audit(state_dir, session_id, action, channel, details, decision)
        return result

    def _ensure_session_in(self, data: dict[str, Any]) -> str:
        if self._session_id is None:
            self._session_id = data.get("current_session") or _new_session_id()
        if not any(s.get("id") == self._session_id for s in data["sessions"]):
            data["sessions"].append({"id": self._session_id, "started": _now()})
        data["current_session"] = self._session_id
        return self._session_id

    def _audit(
        self,
        state_dir: Path,
        session_id: str,
        action: str,
        channel: str | None,
        details: dict[str, Any] | None,
        decision: str | None,
    ) -> None:
        if not self.use_sqlite:
            return
        conn = sqlite3.connect(str(state_dir / AUDIT_DB))
        try:
            conn.execute(
                """
                CREATE TABLE IF NOT EXISTS audit_log (
                    id INTEGER PRIMARY KEY AUTOINCREMENT,
                    timestamp TEXT NOT NULL,
                    session_id TEXT,
                    channel TEXT,
                    action TEXT NOT NULL,
                    details TEXT,
                    user_decision TEXT
                )
                """
            )
            conn.execute(
                "INSERT INTO audit_log (timestamp, session_id, channel, action, details, "
                "user_decision) VALUES (?, ?, ?, ?, ?, ?)",
                (
                    _now(),
                    session_id,
                    channel,
                    action,
                    json.dumps(details, default=json_default) if details else None,
                    decision,
                ),
            )
            conn.commit()
        finally:
            conn.close()

    # ---------------------------------------------------------------- sessions

    def ensure_session(self) -> WorkflowSession:
        """Resume the project's current session, or start one (written on first change)."""
        if self._session_id is None:
            self._session_id = self._read().get("current_session") or _new_session_id()
        return WorkflowSession(self._session_id, str(self.project_path))

    def new_session(self, session_id: str | None = None) -> WorkflowSession:
        """Start a new session and make it the project's current one."""
        self._session_id = session_id or _new_session_id()
        self._mutate("session_created", lambda data: None, details={"session_id": self._session_id})
        return WorkflowSession(self._session_id, str(self.project_path))

    def load_session(self, session_id: str | None = None) -> WorkflowSession | None:
        """Resume ``session_id`` (default: the current one); None if it is unknown."""
        data = self._read()
        session_id = session_id or data.get("current_session")
        if not any(s.get("id") == session_id for s in data["sessions"]):
            return None
        self._session_id = session_id
        return WorkflowSession(session_id, str(self.project_path))

    @property
    def session_id(self) -> str | None:
        return self._session_id

    # ---------------------------------------------------------------- channels

    def add_channel(self, channel_name: str, cycle: str = "unknown", **fields: Any) -> ChannelState:
        """Register a channel; an existing channel keeps its status, steps and parameters."""

        def change(data: dict[str, Any]) -> ChannelState:
            existing = data["channels"].get(channel_name)
            if existing is not None:
                return ChannelState.from_dict(existing)
            state = ChannelState(channel_name=channel_name, cycle=cycle, **fields)
            data["channels"][channel_name] = state.to_dict()
            return state

        if (current := self.get_channel(channel_name)) is not None:
            return current
        return self._mutate("channel_added", change, channel=channel_name)

    def record(
        self,
        channel_name: str,
        operation: str,
        params: dict[str, Any],
        *,
        image: str | None = None,
        cycle: str = "unknown",
        quality_before: float | None = None,
        quality_after: float | None = None,
        status: str | None = "processing",
    ) -> ChannelState:
        """Append a processing step (registering the channel if needed) and set ``status``."""
        step = {
            "operation": operation,
            "params": plain(params),
            "image": image,
            "quality_before": quality_before,
            "quality_after": quality_after,
            "session_id": None,  # filled in under the lock
            "timestamp": _now(),
        }

        def change(data: dict[str, Any]) -> ChannelState:
            raw = data["channels"].get(channel_name)
            state = (
                ChannelState.from_dict(raw)
                if raw is not None
                else ChannelState(channel_name=channel_name, cycle=cycle)
            )
            state.processing_steps.append({**step, "session_id": data["current_session"]})
            if quality_after is not None:
                state.quality_score = quality_after
            if status is not None:
                state.status = status
            state.modified = _now()
            data["channels"][channel_name] = state.to_dict()
            return state

        return self._mutate(
            "step_recorded", change, channel=channel_name, details={"operation": operation}
        )

    def add_processing_step(
        self,
        channel_name: str,
        operation: str,
        params: dict[str, Any],
        quality_before: float | None = None,
        quality_after: float | None = None,
    ) -> bool:
        """Append a step to a known channel; False if the channel is not registered."""
        if self.get_channel(channel_name) is None:
            return False
        self.record(
            channel_name,
            operation,
            params,
            quality_before=quality_before,
            quality_after=quality_after,
            status=None,
        )
        return True

    def update_channel(self, channel_name: str, **fields: Any) -> bool:
        """Update fields of a known channel; False if the channel is not registered."""
        unknown = [k for k in fields if k not in ChannelState.__dataclass_fields__]
        if unknown:
            raise ValueError(f"unknown ChannelState fields: {', '.join(unknown)}")
        if "status" in fields and fields["status"] not in STATUSES:
            raise ValueError(f"status must be one of {STATUSES}, got {fields['status']!r}")

        def change(data: dict[str, Any]) -> bool:
            raw = data["channels"].get(channel_name)
            if raw is None:
                return False
            raw.update(plain(fields))
            raw["modified"] = _now()
            return True

        if self.get_channel(channel_name) is None:
            return False
        return self._mutate("channel_updated", change, channel=channel_name, details=fields)

    def set_user_decision(
        self,
        channel_name: str,
        decision: str,
        notes: str = "",
        *,
        final_params: dict[str, Any] | None = None,
        quality_score: float | None = None,
    ) -> bool:
        """Record a review decision for a known channel; False if it is not registered."""
        if decision not in DECISIONS:
            raise ValueError(f"decision must be one of {DECISIONS}, got {decision!r}")
        if self.get_channel(channel_name) is None:
            return False

        def change(data: dict[str, Any]) -> bool:
            raw = data["channels"].get(channel_name)
            if raw is None:
                return False
            raw["status"] = decision
            raw["review_notes"] = notes
            if final_params is not None:
                raw["final_params"] = plain(final_params)
            if quality_score is not None:
                raw["quality_score"] = float(quality_score)
            raw["modified"] = _now()
            return True

        return self._mutate(
            "user_decision",
            change,
            channel=channel_name,
            details={"notes": notes} if notes else None,
            decision=decision,
        )

    # ------------------------------------------------------------------ queries

    def channels(self) -> dict[str, ChannelState]:
        return {name: ChannelState.from_dict(raw) for name, raw in self._read()["channels"].items()}

    def get_channel(self, channel_name: str) -> ChannelState | None:
        raw = self._read()["channels"].get(channel_name)
        return ChannelState.from_dict(raw) if raw is not None else None

    def get_channels_by_status(self, status: str) -> list[ChannelState]:
        return [ch for ch in self.channels().values() if ch.status == status]

    def get_processing_queue(self) -> list[str]:
        """Channels not yet finished processing (``pending`` or ``processing``)."""
        return [
            ch.channel_name
            for ch in self.channels().values()
            if ch.status in ("pending", "processing", "needs_adjustment")
        ]

    def get_review_queue(self) -> list[str]:
        """Channels awaiting a review decision."""
        return [ch.channel_name for ch in self.get_channels_by_status("review")]

    def summary(self) -> dict[str, Any]:
        """Counts and queues for the whole project (empty when nothing was recorded)."""
        data = self._read()
        channels = {n: ChannelState.from_dict(r) for n, r in data["channels"].items()}
        breakdown: dict[str, int] = {}
        for ch in channels.values():
            breakdown[ch.status] = breakdown.get(ch.status, 0) + 1
        return {
            "session_id": data.get("current_session"),
            "sessions": len(data["sessions"]),
            "total_channels": len(channels),
            "status_breakdown": breakdown,
            "processing_queue": [
                n
                for n, c in channels.items()
                if c.status in ("pending", "processing", "needs_adjustment")
            ],
            "review_queue": [n for n, c in channels.items() if c.status == "review"],
            "approved": [n for n, c in channels.items() if c.status == "approved"],
        }

    def get_session_summary(self) -> dict[str, Any]:
        """Alias of :meth:`summary`."""
        return self.summary()

    def export_parameters(self, output_path: Path | None = None) -> Path:
        """Write the approved channels' parameters and steps to a JSON file."""
        approved = {
            ch.channel_name: {
                "final_params": ch.final_params,
                "processing_steps": ch.processing_steps,
                "quality_score": ch.quality_score,
            }
            for ch in self.get_channels_by_status("approved")
        }
        if output_path is None:
            state_dir = project_state_dir(self.project_path, create=True)
            output_path = state_dir / f"approved_params_{datetime.now():%Y%m%d_%H%M%S}.json"
        Path(output_path).write_text(
            json.dumps({"exported": _now(), "channels": approved}, indent=2, default=json_default)
            + "\n"
        )
        return Path(output_path)
