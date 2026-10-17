"""Convert numpy values into plain, JSON-serialisable Python values."""

from __future__ import annotations

from typing import Any

import numpy as np


def plain(value: Any) -> Any:
    """Recursively turn numpy scalars/arrays (and tuples) into plain Python values."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {str(k): plain(v) for k, v in value.items()}
    if isinstance(value, list | tuple):
        return [plain(v) for v in value]
    return value


def json_default(value: Any) -> Any:
    """``json.dump(default=...)`` hook: numpy values become Python values, anything else str."""
    converted = plain(value)
    return str(value) if converted is value else converted
