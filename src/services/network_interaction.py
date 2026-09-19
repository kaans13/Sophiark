"""Çoklu hedefler için açıklanabilir predicted network interaction analizi."""

from __future__ import annotations

from pathlib import Path

import pandas as pd


METRIC_COLUMNS = {
    "Global_Efficiency_Change": "Global efficiency mutlak değişimi",
    "Systemic_Network_Shift_Pct": "Sistemik ağ kayması (%)",
}


def _extract_effect(report: pd.DataFrame) -> dict[str, float]:
    first = report.iloc[0]
    return {column: float(first.get(column, 0.0) or 0.0) for column in METRIC_COLUMNS}


def _classification(observed: float, expected: float, tolerance: float) -> str:
    scale = max(abs(expected), 1e-12)
    residual_ratio = (abs(observed) - abs(expected)) / scale
    if residual_ratio > tolerance:
        return "supra-additive network effect"
    if residual_ratio < -tolerance:
        return "sub-additive network effect"
    return "approximately additive"


def predicted_network_interaction(
    *, motor_module, graph, scores: pd.DataFrame, targets: list[str], output_root: Path,
    block_strength: float, damping: float, relative_tolerance: float = 0.10,
) -> pd.DataFrame:
    """Tekli etkilerin toplamı ile kombinasyonu metrik-bazlı karşılaştırır.

    Additive baseline her metrik için ``sum(Effect(single target))`` olarak
    tanımlıdır. Tek bir birleşik synergy skoru yoktur. Sınıf toleransı kullanıcıya
    açık ``relative_tolerance`` parametresidir ve rapora yazılır.
    """
    unique = list(dict.fromkeys(map(str, targets)))
    if len(unique) < 2:
        raise ValueError("Predicted Network Interaction için en az iki hedef gerekir")
    individual: list[dict[str, float]] = []
    for target in unique:
        report = motor_module.run_infection_simulation(
            graph.copy(), scores, output_root / "network_interaction" / target,
            spesifik_hedefler=[target], block_weight_fraction=block_strength, damping=damping,
            redistribution_mode="top_n", exploratory_top_n=250, null_iterations=0,
            efficiency_sample_sources=128, local_efficiency_sample_nodes=64,
        )
        individual.append(_extract_effect(report))
    combined_report = motor_module.run_infection_simulation(
        graph.copy(), scores, output_root / "network_interaction" / "combined",
        spesifik_hedefler=unique, block_weight_fraction=block_strength, damping=damping,
        redistribution_mode="top_n", exploratory_top_n=250, null_iterations=0,
        efficiency_sample_sources=128, local_efficiency_sample_nodes=64,
    )
    combined = _extract_effect(combined_report)
    rows = []
    for column, label in METRIC_COLUMNS.items():
        singles = [effect[column] for effect in individual]
        expected = sum(singles)
        observed = combined[column]
        rows.append({
            "Metrik": label,
            "Tekli Etkiler": " + ".join(f"{value:.6g}" for value in singles),
            "Beklenen Additive Etki": expected,
            "Gözlenen Kombine Etki": observed,
            "Residual": observed - expected,
            "Sınıf": _classification(observed, expected, relative_tolerance),
            "Göreli Tolerans": relative_tolerance,
            "Yorum Sınırı": "Hesaplamalı ağ etkileşimidir; biyolojik synergy iddiası değildir.",
        })
    return pd.DataFrame(rows)
