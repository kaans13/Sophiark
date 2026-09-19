"""Computational compensation candidates derived from measured post-perturbation change."""

from __future__ import annotations

import pandas as pd

from src.services.perturbation_comparison import run_single_intervention


def predicted_compensation_candidates(*, motor_module, graph, scores: pd.DataFrame, perturbed_gene: str,
                                      output_root, block_strength: float, damping: float, ghost_filter,
                                      gene_to_symbol: dict[str, str] | None = None) -> pd.DataFrame:
    """Return network predictions; this does not claim biological compensation."""
    comparison = run_single_intervention(
        motor_module=motor_module, graph=graph, scores=scores, intervention_gene=perturbed_gene,
        observed_gene=perturbed_gene, output_root=output_root, block_strength=block_strength,
        damping=damping, ghost_filter=ghost_filter,
    )
    gains = pd.DataFrame(comparison.top_pagerank_gains)
    if gains.empty:
        return gains
    lookup = scores.drop_duplicates("gene").set_index("gene", drop=False)
    rows = []
    for _, gain in gains.iterrows():
        gene = str(gain["gene"])
        if gene == perturbed_gene or gene not in lookup.index:
            continue
        base = lookup.loc[gene]
        symbol = str((gene_to_symbol or {}).get(gene) or base.get("Symbol", gene) or gene)
        gate = bool(base.get("Gümrük_Kapisi", False))
        reasons = ["Perturbasyon sonrası PageRank artışı"]
        if gate:
            reasons.append("Compartment Bottleneck / köprü rolü")
        if float(base.get("BC_Skoru", 0.0) or 0.0) > 0:
            reasons.append("Mevcut BC altyapısı")
        rows.append({
            "Gene": symbol, "Protein ID": gene,
            "ΔPageRank (%)": float(gain["delta_pagerank_pct"]),
            "BC Skoru (baseline)": float(base.get("BC_Skoru", 0.0) or 0.0),
            "Compartment Bottleneck": "Evet" if gate else "Hayır",
            "Doku Uygunluğu": "Aktif doku ağında yer alıyor",
            "Açıklama": "; ".join(reasons),
            "Kanıt Türü": "Hesaplamalı ağ tahmini",
        })
    return pd.DataFrame(rows).sort_values("ΔPageRank (%)", ascending=False).head(25).reset_index(drop=True) if rows else pd.DataFrame()
