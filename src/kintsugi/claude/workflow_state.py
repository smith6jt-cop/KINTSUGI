"""
Workflow State Management for Claude-Guided Processing

Manages state persistence across sessions using JSON, SQLite, and Zarr attributes.
"""

from __future__ import annotations

import json
import sqlite3
from dataclasses import asdict, dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Any


@dataclass
class ChannelState:
    """State for a single channel's processing."""

    channel_name: str
    cycle: str
    status: str = "pending"  # pending, processing, review, approved, rejected
    quality_score: float | None = None
    processing_steps: list[dict[str, Any]] = field(default_factory=list)
    final_params: dict[str, Any] = field(default_factory=dict)
    review_notes: str = ""
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    modified: str = field(default_factory=lambda: datetime.now().isoformat())

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> ChannelState:
        return cls(**{k: v for k, v in data.items() if k in cls.__dataclass_fields__})


@dataclass
class WorkflowSession:
    """A workflow session for processing multiple channels."""

    session_id: str
    project_path: str
    created: str = field(default_factory=lambda: datetime.now().isoformat())
    channels: dict[str, ChannelState] = field(default_factory=dict)
    current_channel: str | None = None
    status: str = "active"  # active, paused, completed

    def to_dict(self) -> dict[str, Any]:
        return {
            "session_id": self.session_id,
            "project_path": self.project_path,
            "created": self.created,
            "channels": {k: v.to_dict() for k, v in self.channels.items()},
            "current_channel": self.current_channel,
            "status": self.status,
        }


class WorkflowState:
    """
    Manages workflow state persistence.

    Supports multiple storage backends:
    - JSON: Simple file-based storage for quick access
    - SQLite: Structured storage with audit trail
    - Zarr: Store state alongside image data

    Example:
        state = WorkflowState("C:/Projects/MyExperiment")

        # Start a new session
        session = state.new_session()

        # Add channel to workflow
        state.add_channel("cyc001_CD3")

        # Update channel status
        state.update_channel("cyc001_CD3", status="approved", quality_score=0.92)

        # Get processing queue
        pending = state.get_channels_by_status("pending")
    """

    def __init__(
        self,
        project_path: str | Path,
        use_sqlite: bool = True,
    ):
        """
        Initialize workflow state manager.

        Args:
            project_path: Path to KINTSUGI project
            use_sqlite: Whether to use SQLite for audit trail
        """
        self.project_path = Path(project_path).resolve()
        self.state_dir = self.project_path / ".claude_workflow"
        self.state_dir.mkdir(parents=True, exist_ok=True)

        self.json_path = self.state_dir / "workflow_state.json"
        self.db_path = self.state_dir / "workflow_audit.db" if use_sqlite else None

        # Load or create state
        self._state: dict[str, Any] = self._load_state()

        # Initialize SQLite if needed
        if use_sqlite:
            self._init_db()

        # Current session
        self._session: WorkflowSession | None = None

    def _load_state(self) -> dict[str, Any]:
        """Load state from JSON file."""
        if self.json_path.exists():
            with open(self.json_path) as f:
                return json.load(f)
        return {"sessions": {}, "current_session_id": None}

    def _save_state(self):
        """Save state to JSON file."""
        with open(self.json_path, "w") as f:
            json.dump(self._state, f, indent=2)

    def _init_db(self):
        """Initialize SQLite database for audit trail."""
        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()

        cursor.execute(
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

        cursor.execute(
            """
            CREATE TABLE IF NOT EXISTS channel_history (
                id INTEGER PRIMARY KEY AUTOINCREMENT,
                timestamp TEXT NOT NULL,
                session_id TEXT,
                channel TEXT NOT NULL,
                operation TEXT NOT NULL,
                parameters TEXT,
                quality_before REAL,
                quality_after REAL
            )
        """
        )

        conn.commit()
        conn.close()

    def _log_audit(
        self,
        action: str,
        channel: str | None = None,
        details: dict | None = None,
        user_decision: str | None = None,
    ):
        """Log an action to the audit trail."""
        if self.db_path is None:
            return

        conn = sqlite3.connect(str(self.db_path))
        cursor = conn.cursor()

        cursor.execute(
            """
            INSERT INTO audit_log (timestamp, session_id, channel, action, details, user_decision)
            VALUES (?, ?, ?, ?, ?, ?)
            """,
            (
                datetime.now().isoformat(),
                self._session.session_id if self._session else None,
                channel,
                action,
                json.dumps(details) if details else None,
                user_decision,
            ),
        )

        conn.commit()
        conn.close()

    def new_session(self, session_id: str | None = None) -> WorkflowSession:
        """
        Create a new workflow session.

        Args:
            session_id: Optional custom session ID (auto-generated if not provided)

        Returns:
            New WorkflowSession
        """
        if session_id is None:
            session_id = datetime.now().strftime("%Y%m%d_%H%M%S")

        session = WorkflowSession(
            session_id=session_id,
            project_path=str(self.project_path),
        )

        self._session = session
        self._state["sessions"][session_id] = session.to_dict()
        self._state["current_session_id"] = session_id

        self._save_state()
        self._log_audit("session_created", details={"session_id": session_id})

        return session

    def load_session(self, session_id: str | None = None) -> WorkflowSession | None:
        """
        Load an existing session.

        Args:
            session_id: Session ID to load (loads current if not provided)

        Returns:
            WorkflowSession if found, None otherwise
        """
        if session_id is None:
            session_id = self._state.get("current_session_id")

        if session_id is None or session_id not in self._state["sessions"]:
            return None

        data = self._state["sessions"][session_id]
        self._session = WorkflowSession(
            session_id=data["session_id"],
            project_path=data["project_path"],
            created=data["created"],
            current_channel=data.get("current_channel"),
            status=data.get("status", "active"),
        )

        # Load channel states
        for name, ch_data in data.get("channels", {}).items():
            self._session.channels[name] = ChannelState.from_dict(ch_data)

        return self._session

    def add_channel(
        self,
        channel_name: str,
        cycle: str = "unknown",
        **kwargs,
    ) -> ChannelState:
        """
        Add a channel to the current workflow.

        Args:
            channel_name: Unique channel identifier
            cycle: Cycle name
            **kwargs: Additional ChannelState fields

        Returns:
            Created ChannelState
        """
        if self._session is None:
            self.new_session()

        state = ChannelState(channel_name=channel_name, cycle=cycle, **kwargs)
        self._session.channels[channel_name] = state

        self._state["sessions"][self._session.session_id] = self._session.to_dict()
        self._save_state()
        self._log_audit("channel_added", channel=channel_name)

        return state

    def update_channel(
        self,
        channel_name: str,
        **kwargs,
    ) -> ChannelState | None:
        """
        Update a channel's state.

        Args:
            channel_name: Channel to update
            **kwargs: Fields to update

        Returns:
            Updated ChannelState
        """
        if self._session is None or channel_name not in self._session.channels:
            return None

        state = self._session.channels[channel_name]

        for key, value in kwargs.items():
            if hasattr(state, key):
                setattr(state, key, value)

        state.modified = datetime.now().isoformat()

        self._state["sessions"][self._session.session_id] = self._session.to_dict()
        self._save_state()
        self._log_audit(
            "channel_updated",
            channel=channel_name,
            details=kwargs,
        )

        return state

    def add_processing_step(
        self,
        channel_name: str,
        operation: str,
        params: dict[str, Any],
        quality_before: float | None = None,
        quality_after: float | None = None,
    ):
        """
        Record a processing step for a channel.

        Args:
            channel_name: Channel being processed
            operation: Operation name (e.g., "blank_subtraction")
            params: Parameters used
            quality_before: Quality score before operation
            quality_after: Quality score after operation
        """
        if self._session is None or channel_name not in self._session.channels:
            return

        step = {
            "operation": operation,
            "params": params,
            "timestamp": datetime.now().isoformat(),
            "quality_before": quality_before,
            "quality_after": quality_after,
        }

        self._session.channels[channel_name].processing_steps.append(step)

        self._state["sessions"][self._session.session_id] = self._session.to_dict()
        self._save_state()

        # Log to SQLite
        if self.db_path:
            conn = sqlite3.connect(str(self.db_path))
            cursor = conn.cursor()
            cursor.execute(
                """
                INSERT INTO channel_history
                (timestamp, session_id, channel, operation, parameters, quality_before, quality_after)
                VALUES (?, ?, ?, ?, ?, ?, ?)
                """,
                (
                    datetime.now().isoformat(),
                    self._session.session_id,
                    channel_name,
                    operation,
                    json.dumps(params),
                    quality_before,
                    quality_after,
                ),
            )
            conn.commit()
            conn.close()

    def set_user_decision(
        self,
        channel_name: str,
        decision: str,
        notes: str = "",
    ):
        """
        Record user's decision for a channel.

        Args:
            channel_name: Channel being reviewed
            decision: "approved", "rejected", or "needs_adjustment"
            notes: Optional review notes
        """
        if self._session is None or channel_name not in self._session.channels:
            return

        self._session.channels[channel_name].status = decision
        self._session.channels[channel_name].review_notes = notes

        self._state["sessions"][self._session.session_id] = self._session.to_dict()
        self._save_state()
        self._log_audit(
            "user_decision",
            channel=channel_name,
            user_decision=decision,
            details={"notes": notes},
        )

    def get_channels_by_status(self, status: str) -> list[ChannelState]:
        """Get all channels with a specific status."""
        if self._session is None:
            return []

        return [ch for ch in self._session.channels.values() if ch.status == status]

    def get_processing_queue(self) -> list[str]:
        """Get ordered list of channels pending processing."""
        pending = self.get_channels_by_status("pending")
        return [ch.channel_name for ch in pending]

    def get_review_queue(self) -> list[str]:
        """Get ordered list of channels awaiting review."""
        review = self.get_channels_by_status("review")
        return [ch.channel_name for ch in review]

    def get_session_summary(self) -> dict[str, Any]:
        """Get summary of current session."""
        if self._session is None:
            return {"error": "No active session"}

        channels = self._session.channels
        status_counts = {}
        for ch in channels.values():
            status_counts[ch.status] = status_counts.get(ch.status, 0) + 1

        return {
            "session_id": self._session.session_id,
            "created": self._session.created,
            "status": self._session.status,
            "total_channels": len(channels),
            "status_breakdown": status_counts,
            "current_channel": self._session.current_channel,
            "approved_channels": [
                ch.channel_name for ch in channels.values() if ch.status == "approved"
            ],
        }

    def export_parameters(self, output_path: Path | None = None) -> Path:
        """
        Export all approved parameters to a JSON file.

        Useful for batch processing or documentation.
        """
        if self._session is None:
            raise ValueError("No active session")

        if output_path is None:
            output_path = self.state_dir / f"approved_params_{self._session.session_id}.json"

        approved = self.get_channels_by_status("approved")
        params = {
            ch.channel_name: {
                "final_params": ch.final_params,
                "processing_steps": ch.processing_steps,
                "quality_score": ch.quality_score,
            }
            for ch in approved
        }

        with open(output_path, "w") as f:
            json.dump(
                {
                    "session_id": self._session.session_id,
                    "exported": datetime.now().isoformat(),
                    "channels": params,
                },
                f,
                indent=2,
            )

        return output_path
