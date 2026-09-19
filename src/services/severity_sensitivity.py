"""Perturbation attenuation senaryolarını ham metriklerle karşılaştırır."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src.scientific.perturbation import SEVERITY_FRACTIONS


def run_severity_sensitivity(
    *, motor_module, graph, scores: pd.DataFrame, targets: list[str], output_root: Path,
    damping: float,
) -> pd.DataFrame:
    rows = []
    for severity, fraction in SEVERITY_FRACTIONS.items():
        report, raw = motor_module.run_infection_simulation(
            graph.copy(), scores, output_root / "severity_sensitivity" / severity,
            spesifik_hedefler=targets, block_weight_fraction=fraction, damping=damping,
            redistribution_mode="top_n", exploratory_top_n=250, null_iterations=0,
            compute_structural_metrics=True, efficiency_sample_sources=128,
            local_efficiency_sample_nodes=64, return_comparison_metrics=True,
            comparison_genes=targets,
        )
        first = report.iloc[0]
        redistribution = pd.DataFrame(raw.get("redistribution", []))
        rows.append({
            "Senaryo": severity,
            "Kalan Edge Strength": fraction,
            "Global Efficiency Değişimi (%)": float(first["Global_Efficiency_Change_Pct"]),
            "Sistemik Ağ Kayması (%)": float(first["Systemic_Network_Shift_Pct"]),
            "Seçili |ΔPageRank| Burden": float(redistribution.get("Abs_Delta_PageRank", pd.Series(dtype=float)).sum()),
            "Strong Redistribution Gen": int(redistribution.get("Strong_Redistribution", pd.Series(dtype=bool)).fillna(False).sum()),
        })
    return pd.DataFrame(rows)
