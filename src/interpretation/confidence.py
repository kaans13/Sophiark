"""Transparent confidence labels for interpretation evidence coverage."""

from __future__ import annotations


def confidence_label(*, has_delta: bool, has_baseline: bool, annotation_count: int,
                     directed_evidence_count: int = 0) -> tuple[float, str]:
    """Score coverage, not biological truth or scientific significance."""
    score = 0.20
    score += 0.30 if has_delta else 0.0
    score += 0.20 if has_baseline else 0.0
    score += min(0.20, annotation_count * 0.05)
    score += min(0.10, directed_evidence_count * 0.05)
    score = min(0.95, score)
    label = "yüksek" if score >= 0.70 else "orta" if score >= 0.45 else "düşük"
    return score, label
