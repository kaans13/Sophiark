"""Schema-tolerant normalization of existing biological annotation fields."""

from __future__ import annotations

import json
import re
from collections.abc import Mapping
from typing import Iterable

import pandas as pd


_SPLIT = re.compile(r"\s*(?:\||;|\n|,\s*(?=[A-Z]{2,}|GO:|KEGG:))\s*")


def normalized_text(value: object) -> str:
    """Return one comparison-safe string; never fetch or invent an annotation."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, Mapping):
        return normalized_text(value.get("term") or value.get("name") or "")
    if isinstance(value, (list, tuple, set)):
        return " | ".join(normalized_text(item) for item in value if normalized_text(item))
    text = str(value).strip()
    if not text or text.casefold() in {"nan", "none", "n/a", "na", "veri yok"}:
        return ""
    # A few existing cache fields are JSON arrays; preserve their literal terms.
    if text.startswith("[") and text.endswith("]"):
        try:
            parsed = json.loads(text)
        except (TypeError, ValueError):
            pass
        else:
            if isinstance(parsed, list):
                return normalized_text(parsed)
    return re.sub(r"\s+", " ", text)


def normalized_terms(value: object) -> tuple[str, ...]:
    """Split display/cache values into stable, de-duplicated source terms."""
    text = normalized_text(value)
    if not text:
        return ()
    return tuple(dict.fromkeys(part.strip() for part in _SPLIT.split(text) if part.strip()))


def row_values(row: pd.Series, columns: Iterable[str]) -> tuple[str, ...]:
    """Read all present aliases while keeping source values auditable."""
    values: list[str] = []
    for column in columns:
        if column in row.index:
            values.extend(normalized_terms(row[column]))
    return tuple(dict.fromkeys(values))
