"""Tissue-by-tissue comparison using unchanged forward simulations."""

from __future__ import annotations

from pathlib import Path
import re

import pandas as pd

from src.models import TissueDifferentialResult
from src.services.analysis_service import execute_simulation, prepare_network


def _safe_directory_name(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_") or "tissue"


def _canonical_targets(targets: list[str], available: set[str]) -> list[str]:
    """Accept the project’s optional taxon prefix without guessing symbols."""
    usable: list[str] = []
    for target in targets:
        raw = str(target)
        candidates = (raw, raw.removeprefix("9606."), raw.removeprefix("10090."))
        match = next((item for item in candidates if item in available), None)
        if match and match not in usable:
            usable.append(match)
    return usable


def _jaccard(left: set[str], right: set[str]) -> float | None:
    union = left | right
    return (len(left & right) / len(union)) if union else None


def run_tissue_differential(*, motor_module, tissues: list[str], targets: list[str], output_root: Path,
                            bc_sample_sources: int | None, block_strength: float, damping: float,
                            ghost_filter, apply_druggability) -> TissueDifferentialResult:
    """Run the existing forward simulation once per requested tissue.

    The function is sequential by design: the cached graph is mutable and the
    existing engine is CPU-bound.  It isolates CSV artifacts per tissue and
    makes no approximation or cross-tissue pooling.
    """
    summaries: list[dict[str, object]] = []
    gene_rows: list[dict[str, object]] = []
    notices: list[str] = []
    critical_sets: dict[str, set[str]] = {}

    for tissue in dict.fromkeys(tissues):
        prepared = prepare_network(
            motor_module, tuple(targets), tissue, bc_sample_sources,
            "cross_tissue_comparable",
        )
        graph, scores = prepared.graph, prepared.scores
        usable = _canonical_targets(targets, set(map(str, graph.vs["name"])))
        if not usable:
            notices.append(f"{tissue}: hedef seçili doku ağında bulunmadı; bu doku karşılaştırmaya alınmadı.")
            continue
        tissue_output = output_root / "tissue_differential" / _safe_directory_name(tissue)
        tissue_output.mkdir(parents=True, exist_ok=True)
        report, _ = execute_simulation(
            motor_module=motor_module, graph=graph, scores=scores, output_dir=tissue_output,
            targets=usable, localization="SNIPER_MODU", block_strength=block_strength, damping=damping,
            apply_druggability=apply_druggability, ghost_filter=ghost_filter,
            scientific_options={
                "redistribution_mode": "top_n", "exploratory_top_n": 250,
                "null_iterations": 0, "tissue": tissue,
                "tissue_normalization_mode": "cross_tissue_comparable",
                "bc_mode": "fast_approximate" if bc_sample_sources else "full_validation",
                "bc_sample_size": bc_sample_sources,
            },
        )
        if report.empty:
            notices.append(f"{tissue}: forward simülasyon raporu üretilemedi.")
            continue
        stress_column = report.get("Hasar_Tipi", pd.Series("", index=report.index)).astype(str)
        stressed = report[stress_column == str(motor_module.HASAR_UCUNCUL)].copy()
        significance = stressed.get("significant_redistribution", pd.Series(False, index=stressed.index)).fillna(False)
        strong = stressed.get("Strong_Redistribution", pd.Series(False, index=stressed.index)).fillna(False)
        critical = stressed[significance | strong].copy()
        if "Abs_Delta_PageRank" in critical:
            critical = critical.sort_values("Abs_Delta_PageRank", ascending=False)
        else:
            critical = critical.sort_values("Hinterland_Skoru", ascending=False)
        critical_sets[tissue] = set(critical["gene"].astype(str))
        summaries.append({
            "Doku": tissue,
            "Ağ Düğümü": graph.vcount(),
            "Ağ Kenarı": graph.ecount(),
            "Stresli Gen": len(stressed),
            "Kritik Gen": len(critical),
            "Global Ağ Verimliliği Değişimi (%)": float(report["Global_Efficiency_Change_Pct"].iloc[0]),
            "Ortalama Local Efficiency Değişimi (%)": float(report["Mean_Local_Efficiency_Change_Pct"].iloc[0]),
            "Sistemik Ağ Kayması (%)": float(report["Systemic_Network_Shift_Pct"].iloc[0]),
            "Hedef Ağda": ", ".join(usable),
        })
        for _, row in critical.iterrows():
            gene_rows.append({
                "Doku": tissue,
                "Gen": row.get("Symbol", row["gene"]),
                "gene": str(row["gene"]),
                "Tier": row.get("Kategori", "—"),
                "Hinterland": float(row.get("Hinterland_Skoru", 0) or 0),
                "Response Direction": row.get("Response_Direction", "—"),
                "|ΔPageRank|": float(row.get("Abs_Delta_PageRank", 0) or 0),
                "q-value": float(row.get("q_value", 1) or 1),
                "Seçim Nedeni": row.get("Selection_Reason", "Strong redistribution"),
            })

    summary = pd.DataFrame(summaries)
    genes = pd.DataFrame(gene_rows)
    if not genes.empty:
        counts = genes.groupby("gene")["Doku"].nunique().rename("Doku Sayısı")
        genes = genes.join(counts, on="gene")
        included_tissue_count = len(summary)
        genes["Yanıt Tipi"] = genes["Doku Sayısı"].map(
            lambda count: "Korunmuş kritik değişim" if count == included_tissue_count else "Dokuya özgü kritik değişim"
        )

    similarity_rows: list[dict[str, object]] = []
    ordered_tissues = list(critical_sets)
    for index, left in enumerate(ordered_tissues):
        for right in ordered_tissues[index + 1:]:
            similarity_rows.append({
                "Doku A": left,
                "Doku B": right,
                "Kritik gen Jaccard": _jaccard(critical_sets[left], critical_sets[right]),
                "Ortak kritik gen": len(critical_sets[left] & critical_sets[right]),
            })
    return TissueDifferentialResult(
        summary=summary, critical_genes=genes, notices=notices, similarity=pd.DataFrame(similarity_rows)
    )
