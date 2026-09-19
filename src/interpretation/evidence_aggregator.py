"""Normalize existing directed-evidence rows without inferring direction."""

from __future__ import annotations

import pandas as pd


def aggregate_directed_evidence(evidence: pd.DataFrame | None, symbols: set[str] | None = None) -> list[dict[str, object]]:
    """Return explicit directional records already present in TRRUST/OmniPath output."""
    if evidence is None or evidence.empty or not {"source", "target"}.issubset(evidence.columns):
        return []
    frame = evidence.copy(deep=True)
    if symbols:
        source = frame["source"].fillna("").astype(str).str.upper()
        target = frame["target"].fillna("").astype(str).str.upper()
        frame = frame[source.isin(symbols) | target.isin(symbols)]
    rows: list[dict[str, object]] = []
    for _, row in frame.iterrows():
        directed = bool(row.get("consensus_direction") or row.get("is_directed"))
        if not directed:
            continue
        rows.append({
            "source": str(row["source"]), "target": str(row["target"]),
            "stimulation": bool(row.get("consensus_stimulation") or row.get("is_stimulation")),
            "inhibition": bool(row.get("consensus_inhibition") or row.get("is_inhibition")),
            "provenance": str(row.get("provenance") or row.get("sources") or "belirtilmemiş"),
            "references": str(row.get("references") or ""),
        })
    return rows
