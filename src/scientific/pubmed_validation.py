"""Snapshot-temelli, leakage yapmayan Literature Concordance validator'Ä±."""
from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable

import pandas as pd


def _hypergeom_sf(overlap: int, universe: int, literature: int, predicted: int) -> float:
    # P(X >= overlap), stdlib ile deterministic ve ek dependency olmadan.
    if overlap <= 0:
        return 1.0
    denominator = math.comb(universe, predicted)
    total = sum(math.comb(literature, i) * math.comb(universe - literature, predicted - i)
                for i in range(overlap, min(literature, predicted) + 1)
                if 0 <= predicted - i <= universe - literature)
    return min(1.0, total / denominator) if denominator else 1.0


def extract_symbols_from_texts(texts: Iterable[str], allowed_symbols: Iterable[str]) -> tuple[set[str], list[str]]:
    """YalnÄ±zca Ã¶nceden tanÄ±mlÄ± symbol evreninde exact-token mention extraction.

    Bu NER veya LLM deÄŸildir; kÄ±sa/ambiguous semboller Ã¶nlemek iÃ§in 3 karakterden
    kÄ±sa sembolleri otomatik Ã§Ä±karmaz ve excluded listesinde raporlar.
    """
    text = "\n".join(str(item) for item in texts).upper()
    found, excluded = set(), []
    for symbol in sorted({str(item).upper() for item in allowed_symbols}):
        if len(symbol) < 3:
            excluded.append(symbol)
            continue
        if re.search(rf"(?<![A-Z0-9]){re.escape(symbol)}(?![A-Z0-9])", text):
            found.add(symbol)
    return found, excluded


def calculate_concordance(predicted: Iterable[str], literature: Iterable[str], universe: Iterable[str]) -> dict:
    predicted_set, literature_set, universe_set = set(predicted), set(literature), set(universe)
    predicted_set &= universe_set
    literature_set &= universe_set
    intersection = predicted_set & literature_set
    denom = len(predicted_set | literature_set)
    return {
        "metric_label": "Literature Concordance",
        "predicted_gene_count": len(predicted_set), "literature_gene_count": len(literature_set),
        "intersection_count": len(intersection), "intersection": sorted(intersection),
        "precision_like_overlap": len(intersection) / len(predicted_set) if predicted_set else 0.0,
        "recall_like_overlap": len(intersection) / len(literature_set) if literature_set else 0.0,
        "jaccard": len(intersection) / denom if denom else 0.0,
        "hypergeometric_p_value": _hypergeom_sf(len(intersection), len(universe_set), len(literature_set), len(predicted_set)),
        "analysis_universe_size": len(universe_set),
    }


def write_validation_snapshot(
    output_dir: str | Path, *, query: dict, pmids: list[str], literature_symbols: Iterable[str],
    predictions: pd.DataFrame, concordance: dict, mapping: pd.DataFrame,
) -> Path:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    payload = {**query, "query_date_utc": query.get("query_date_utc", datetime.now(timezone.utc).isoformat()), "pmid_count": len(pmids)}
    (root / "query.json").write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    (root / "pmids.json").write_text(json.dumps(sorted(pmids), indent=2), encoding="utf-8")
    pd.DataFrame({"literature_symbol": sorted(set(literature_symbols))}).to_csv(root / "literature_genes.csv", index=False)
    predictions.to_csv(root / "sophiark_predictions.csv", index=False)
    mapping.to_csv(root / "mapping.csv", index=False)
    (root / "overlap.json").write_text(json.dumps(concordance, indent=2, ensure_ascii=False), encoding="utf-8")
    return root
