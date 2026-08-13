"""Small, dependency-light I/O helpers shared by Ino-seq modules."""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from pathlib import Path

import pandas as pd


def ensure_parent(path: str | Path) -> Path:
    """Create the parent directory of *path* and return a resolved Path."""
    output = Path(path)
    output.parent.mkdir(parents=True, exist_ok=True)
    return output


def require_columns(frame: pd.DataFrame, columns: Sequence[str], source: str | Path) -> None:
    """Raise a readable error when a tabular input is missing required columns."""
    missing = [column for column in columns if column not in frame.columns]
    if missing:
        raise ValueError(
            f"{source} is missing required column(s): {', '.join(missing)}; "
            f"available columns: {', '.join(map(str, frame.columns))}"
        )


def split_query_names(values: Iterable[object]) -> list[str]:
    """Expand comma-delimited query-name cells while retaining first-seen order."""
    names: list[str] = []
    seen: set[str] = set()
    for value in values:
        if pd.isna(value):
            continue
        for raw_name in str(value).split(","):
            name = raw_name.strip()
            if name and name not in seen:
                names.append(name)
                seen.add(name)
    return names
