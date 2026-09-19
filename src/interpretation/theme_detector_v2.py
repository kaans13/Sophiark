"""Evidence-backed functional-theme aggregation for deterministic interpretation."""

from __future__ import annotations

from collections import defaultdict
from typing import Iterable

from .theme_ontology import THEME_ONTOLOGY, ThemeDefinition

_BROAD_THEME_IDS = {"membrane_transport"}


def _annotation_text(candidate: dict[str, object]) -> str:
    annotations = candidate.get("annotations", {})
    if not isinstance(annotations, dict):
        return ""
    return " | ".join(
        str(item) for values in annotations.values() if isinstance(values, (list, tuple)) for item in values
    ).casefold()


def detect_themes_from_facts(candidate_facts: Iterable[dict[str, object]], *, min_genes: int = 1,
                             direction: str = "gain") -> list[dict[str, object]]:
    """Report repeated source-term overlap, never a pathway-state inference."""
    candidate_facts = list(candidate_facts)
    found: dict[str, list[str]] = defaultdict(list)
    evidence_terms: dict[str, list[str]] = defaultdict(list)
    for candidate in candidate_facts:
        symbol = str(candidate.get("symbol") or "").strip().upper()
        if not symbol:
            continue
        text = _annotation_text(candidate)
        for definition in THEME_ONTOLOGY:
            matched = [keyword for keyword in definition.keywords if keyword in text]
            if matched:
                if symbol not in found[definition.identifier]:
                    found[definition.identifier].append(symbol)
                evidence_terms[definition.identifier].extend(matched)
    response_magnitude = {
        str(candidate.get("symbol") or "").upper(): max(0.0, -float(candidate.get("delta_pagerank_pct") or 0.0))
        if direction == "loss" else max(0.0, float(candidate.get("delta_pagerank_pct") or 0.0))
        for candidate in candidate_facts
    }
    themes: list[dict[str, object]] = []
    for definition in THEME_ONTOLOGY:
        genes = found.get(definition.identifier, [])
        if len(genes) < min_genes:
            continue
        terms = list(dict.fromkeys(evidence_terms[definition.identifier]))
        response_weight = sum(response_magnitude.get(symbol, 0.0) for symbol in genes)
        specificity_multiplier = .30 if definition.identifier in _BROAD_THEME_IDS else 1.0
        themes.append({
            "theme_id": definition.identifier, "theme": definition.label,
            "mechanism_class": definition.mechanism_class, "genes": genes,
            "gene_count": len(genes), "supporting_keywords": terms,
            "confidence_score": min(0.95, 0.45 + 0.12 * len(genes) + 0.04 * len(terms)),
            "interpretation_priority": round(specificity_multiplier * (response_weight + .08 * len(genes)), 6),
            "evidence": "Mevcut aday anotasyonlarında tema anahtar sözcüğü eşleşmesi.",
            "warning": "Bu tema, yolak aktivasyonu veya inhibisyonu kanıtı değildir.",
        })
    return sorted(themes, key=lambda item: (-float(item["interpretation_priority"]), -int(item["gene_count"]), str(item["theme"])))
