"""Deterministic, source-bounded Turkish candidate interpretation."""

from __future__ import annotations

from typing import Iterable, Mapping

from .biological_rules import fmt_percent
from .confidence import confidence_label
from .sentence_bank import SENTENCES


def validation_suggestions(*, has_annotation: bool, has_localization: bool) -> list[str]:
    """Return recommendations, never an implied result or clinical action."""
    result = [SENTENCES["validation_transcript_protein"], SENTENCES["validation_perturbation"]]
    if has_annotation or has_localization:
        result.append(SENTENCES["validation_localization"])
    return result


def compose_candidate_interpretation(
    fact: Mapping[str, object], *, theme_labels: Iterable[str] = (),
    cluster_names: Iterable[str] = (), theme_peer_symbols: Iterable[str] = (), directed_evidence_count: int = 0,
) -> dict[str, object]:
    """Compose only claims whose inputs are explicitly present in `fact`."""
    symbol = str(fact.get("symbol") or fact.get("gene") or "Bilinmeyen gen")
    delta = fact.get("delta_pagerank_pct")
    delta_value = float(delta) if isinstance(delta, (int, float)) else None
    if delta_value is None:
        first = f"{symbol} için ΔPageRank değeri mevcut değildir; aday yorumu üretilemedi."
    elif delta_value >= 0:
        first = SENTENCES["positive_delta"].format(symbol=symbol, delta=fmt_percent(delta_value, signed=True))
    else:
        first = SENTENCES["negative_delta"].format(symbol=symbol, delta=fmt_percent(delta_value, signed=True))
    hinterland = fact.get("hinterland_baseline")
    bc = fact.get("bc_baseline")
    has_baseline = isinstance(hinterland, (int, float)) or isinstance(bc, (int, float))
    if has_baseline:
        baseline = SENTENCES["baseline_context"].format(
            hinterland=(f"{float(hinterland):.4g}" if isinstance(hinterland, (int, float)) else "veri yok"),
            bc=(f"{float(bc):.4g}" if isinstance(bc, (int, float)) else "veri yok"),
        )
    else:
        baseline = SENTENCES["baseline_missing"]
    themes = tuple(dict.fromkeys(str(item) for item in theme_labels if str(item).strip()))
    clusters = tuple(dict.fromkeys(str(item) for item in cluster_names if str(item).strip()))
    if themes:
        context = SENTENCES["annotation_context"].format(themes=", ".join(themes))
        peers = tuple(
            peer for peer in dict.fromkeys(str(item) for item in theme_peer_symbols)
            if peer and peer.upper() != symbol.upper()
        )
        if peers:
            context += (
                f" Aynı analizde {', '.join(peers[:4])} ile birlikte bu eksende görünmesi, "
                "ortak anotasyon bağlamını güçlendirir; doğrudan moleküler etkileşim kanıtı değildir."
            )
    else:
        context = SENTENCES["low_evidence"]
    cluster_context = SENTENCES["family_context"].format(clusters=", ".join(clusters)) if clusters else ""
    annotations = fact.get("annotations")
    annotation_count = sum(len(items) for items in annotations.values()) if isinstance(annotations, dict) else 0
    score, label = confidence_label(
        has_delta=delta_value is not None, has_baseline=has_baseline,
        annotation_count=annotation_count, directed_evidence_count=directed_evidence_count,
    )
    return {
        "text": " ".join(item for item in (first, baseline, context, cluster_context) if item),
        "confidence_score": score, "confidence_label": label,
        "validation_suggestions": validation_suggestions(
            has_annotation=bool(themes),
            has_localization=bool(isinstance(annotations, dict) and annotations.get("localization")),
        ),
    }
