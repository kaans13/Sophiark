"""Explainable two-stage Target Stress Search using Sophiark's forward engine."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from src.services.perturbation_comparison import full_target_bc_delta, run_single_intervention


def _column(frame: pd.DataFrame, *names: str, default: float = 0.0) -> pd.Series:
    for name in names:
        if name in frame.columns:
            return pd.to_numeric(frame[name], errors="coerce").fillna(default)
    return pd.Series(default, index=frame.index, dtype=float)


def _pareto_layers(frame: pd.DataFrame) -> pd.Series:
    """Return non-dominated sorting layers; no weighted biological score."""
    values = frame[["_distance", "_pagerank", "_bc", "_hinterland", "_gate"]].to_numpy(float, copy=True)
    values[:, 0] *= -1.0  # closer distance is better; the rest are higher-is-better
    remaining = list(range(len(frame)))
    layers = np.zeros(len(frame), dtype=int)
    layer = 1
    while remaining:
        front = []
        for i in remaining:
            dominated = False
            for j in remaining:
                if i == j:
                    continue
                if np.all(values[j] >= values[i]) and np.any(values[j] > values[i]):
                    dominated = True
                    break
            if not dominated:
                front.append(i)
        for i in front:
            layers[i] = layer
        remaining = [i for i in remaining if i not in front]
        layer += 1
    return pd.Series(layers, index=frame.index)


def preselect_candidates(
    *, graph, scores: pd.DataFrame, target_gene: str, regulators: dict,
    gene_to_symbol: dict[str, str], mode: str, limit: int,
) -> pd.DataFrame:
    """Build an evidence-preserving candidate set then Pareto-filter it.

    TRRUST evidence is intentionally only considered direct evidence; a missing
    record never excludes a network-mediated candidate.
    """
    names = set(graph.vs["name"])
    if target_gene not in names:
        raise ValueError("Target gene is not present in the active tissue network.")
    target_row = scores[scores["gene"].astype(str) == str(target_gene)]
    target_symbol = (str(target_row.iloc[0].get("Symbol", "")) if not target_row.empty else "")
    target_symbol = target_symbol or gene_to_symbol.get(target_gene, "")
    direct_symbols = regulators.get(target_symbol, []) if target_symbol else []
    direct_map = {str(item[0]): str(item[1]) if isinstance(item, (tuple, list)) and len(item) > 1 else "reported regulation"
                  for item in direct_symbols}

    name_to_idx = {name: index for index, name in enumerate(graph.vs["name"])}
    target_idx = name_to_idx[target_gene]
    neighbor_ids = {graph.vs[i]["name"] for i in graph.neighbors(target_idx)}
    candidates = set(neighbor_ids)
    for name, symbol in gene_to_symbol.items():
        if symbol in direct_map and name in names:
            candidates.add(name)

    # Same-community high-centrality candidates add non-neighbour pathways
    # without expanding to the entire proteome.
    if not target_row.empty and "Topluluk_ID" in scores.columns:
        community = target_row.iloc[0].get("Topluluk_ID")
        same = scores[scores["Topluluk_ID"] == community].copy()
        same["_priority"] = _column(same, "Hinterland_Skoru")
        candidates.update(same.nlargest(max(30, limit * 3), "_priority")["gene"].astype(str).tolist())

    candidates.discard(target_gene)
    candidates = [gene for gene in candidates if gene in names]
    if not candidates:
        return pd.DataFrame()

    lookup = scores.drop_duplicates("gene").set_index("gene", drop=False)
    distances = graph.distances(source=target_idx, target=[name_to_idx[g] for g in candidates], weights="distance")[0]
    rows = []
    for gene, distance in zip(candidates, distances):
        row = lookup.loc[gene] if gene in lookup.index else pd.Series(dtype=object)
        symbol = str(row.get("Symbol", "") or gene_to_symbol.get(gene, gene))
        directed = direct_map.get(symbol, "")
        if mode in ("Direct Influence", "Doğrudan Etki") and not directed:
            continue
        rows.append({
            "gene": gene, "Candidate": symbol if symbol and symbol != "nan" else gene,
            "Network Distance": float(distance) if np.isfinite(distance) else np.inf,
            "Directed Evidence": directed or "Yüklü TRRUST verisinde yönlü kanıt yok",
            "Tissue Relevance": "Aktif doku ağında yer alıyor",
            "_pagerank": float(row.get("PageRank", row.get("pagerank", 0.0)) or 0.0),
            "_bc": float(row.get("BC_Skoru", 0.0) or 0.0),
            "_hinterland": float(row.get("Hinterland_Skoru", 0.0) or 0.0),
            "_gate": float(bool(row.get("Gümrük_Kapisi", False))),
            "_distance": float(distance) if np.isfinite(distance) else 999999.0,
            "_is_neighbor": gene in neighbor_ids,
            "_direct": bool(directed),
        })
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Pareto Layer"] = _pareto_layers(out)
    # A deterministic tie-breaker only orders candidates within the same
    # non-dominated layer; it is not a new combined biological score.
    out = out.sort_values(["Pareto Layer", "_distance", "_hinterland", "_bc"], ascending=[True, True, False, False])
    return out.head(limit).reset_index(drop=True)


def validate_candidates(
    *, candidates: pd.DataFrame, motor_module, graph, scores: pd.DataFrame,
    target_gene: str, output_root: Path, block_strength: float, damping: float,
    ghost_filter, deep_validation_top_n: int = 0,
) -> pd.DataFrame:
    """Validate each filtered candidate with the existing forward motor."""
    results = []
    for _, candidate in candidates.iterrows():
        comparison = run_single_intervention(
            motor_module=motor_module, graph=graph, scores=scores,
            intervention_gene=str(candidate["gene"]), observed_gene=target_gene,
            output_root=output_root, block_strength=block_strength, damping=damping,
            ghost_filter=ghost_filter,
        )
        observed = comparison.observed_metrics.get(target_gene, {})
        delta_pr = float(observed.get("delta_pagerank_pct", 0.0))
        delta_bc = float(observed.get("delta_local_bc_2hop", float("nan")))
        target_in_stressed = target_gene in set(comparison.report.loc[
            comparison.report.get("Hasar_Tipi", pd.Series(dtype=str)).astype(str).str.contains("\u00dc\u00c7\u00dcNC", na=False), "gene"
        ].astype(str))
        mechanism = []
        if bool(candidate["_is_neighbor"]):
            mechanism.append("PPI komşusu")
        if bool(candidate["_direct"]):
            mechanism.append("TRRUST yönlü kanıtı")
        if float(candidate["Network Distance"]) > 1:
            mechanism.append("aynı topluluk/topolojik aday")
        mechanism.append("forward perturbation ile doğrulandı")
        confidence = "Ağ üzerinde doğrulandı"
        if bool(candidate["_direct"]):
            confidence += " + yönlü kanıt"
        results.append({
            "Aday": candidate["Candidate"], "gen": candidate["gene"],
            "Öngörülen Hedef Stresi": abs(delta_pr),
            "ΔPageRank (%)": delta_pr,
            "Δ Yerel BC (2-hop projeksiyonu)": delta_bc,
            "Ağ Mesafesi": candidate["Network Distance"],
            "Yönlü Kanıt": candidate["Directed Evidence"],
            "Doku Uygunluğu": candidate["Tissue Relevance"],
            "Mekanizma / Açıklama": "; ".join(mechanism),
            "Güven / Kanıt Gücü": confidence,
            "Sistemik Kayma (%)": comparison.system_shift_pct,
            "En Güçlü Pozitif PageRank Yeniden-Dağılımı (%)": comparison.local_stress_pct,
            "Hedef stresli kümede": bool(target_in_stressed),
            "Pareto Katmanı": int(candidate["Pareto Layer"]),
        })
    output = pd.DataFrame(results).sort_values(
        ["Öngörülen Hedef Stresi", "Pareto Katmanı"], ascending=[False, True]
    ).reset_index(drop=True) if results else pd.DataFrame()
    if output.empty:
        return output
    output["Full ΔBC"] = np.nan
    output["BC Doğrulama Modu"] = "2-hop projection"
    for index in output.head(max(0, int(deep_validation_top_n))).index:
        output.loc[index, "Full ΔBC"] = full_target_bc_delta(
            graph, str(output.loc[index, "gen"]), target_gene, block_strength
        )
        output.loc[index, "BC Doğrulama Modu"] = "Full Validation"
    return output
