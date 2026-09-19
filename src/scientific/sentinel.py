"""Audit edilebilir local-first sentinel identifier resolution."""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import pandas as pd


SENTINELS = {
    "ENSP00000269305": ("TP53", 90.0), "ENSP00000478887": ("MYC", 85.0),
    "ENSP00000418960": ("BRCA1", 70.0), "ENSP00000275493": ("EGFR", 75.0),
    "ENSP00000335153": ("HSP90AA1", 80.0), "ENSP00000344818": ("UBC", 80.0),
    "ENSP00000344456": ("CTNNB1", 80.0), "ENSP00000382883": ("AKT1", 80.0),
    "ENSP00000263253": ("EP300", 80.0), "ENSP00000346847": ("FN1", 80.0),
    "ENSP00000395213": ("CDK1", 80.0), "ENSP00000253840": ("ESR1", 80.0),
}


def resolve_sentinels(scores: pd.DataFrame, mapping_path: str | Path) -> pd.DataFrame:
    """Legacy ENSP'leri local symbol tablosuyla gÃ¼ncel network node'una Ã§Ã¶zer.

    Bu registry bilinÃ§li olarak manuel ve version-control edilebilir: deprecated
    protein ID'ler iÃ§in canlÄ± API yerine kaynak sembol kullanÄ±lÄ±r.
    """
    mapping = pd.read_csv(mapping_path, usecols=["gene", "Symbol"]).dropna()
    by_symbol = mapping.groupby("Symbol")["gene"].agg(list).to_dict()
    network_nodes = set(scores["gene"].astype(str))
    rows = []
    for original, (symbol, min_score) in SENTINELS.items():
        direct = original if original in network_nodes else None
        candidates = [item for item in by_symbol.get(symbol, []) if item in network_nodes]
        if direct:
            status, node, source = "PASS" if float(scores.loc[scores["gene"] == direct, "Hinterland_Skoru"].iloc[0]) >= min_score else "FAIL", direct, "direct_network_id"
        elif len(candidates) == 1:
            node, source = candidates[0], "local_ensp_with_symbols"
            score = float(scores.loc[scores["gene"] == node, "Hinterland_Skoru"].iloc[0])
            status = "PASS" if score >= min_score else "FAIL"
        elif len(candidates) > 1:
            status, node, source = "AMBIGUOUS_MAPPING", None, "local_symbol_multiple_nodes"
        else:
            status, node, source = "UNMAPPED", None, "local_mapping_not_found"
        score = float(scores.loc[scores["gene"] == node, "Hinterland_Skoru"].iloc[0]) if node else None
        rows.append({
            "original_identifier": original, "resolved_symbol": symbol,
            "resolved_network_node_id": node, "mapping_source": source,
            "species": "Homo sapiens", "status": status,
            "min_hinterland_score": min_score, "observed_hinterland_score": score,
        })
    return pd.DataFrame(rows)


def sentinel_summary(frame: pd.DataFrame) -> dict[str, float | int]:
    evaluated = frame[frame["status"].isin(["PASS", "FAIL"])]
    return {
        "evaluated_success_count": int((evaluated["status"] == "PASS").sum()),
        "evaluated_count": int(len(evaluated)),
        "total_count": int(len(frame)),
        "coverage_rate": float(len(evaluated) / len(frame)) if len(frame) else 0.0,
        "evaluated_success_rate": float((evaluated["status"] == "PASS").mean()) if len(evaluated) else 0.0,
    }
