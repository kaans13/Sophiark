"""Forensic all-node PageRank independence audit for GLP1R in Pancreas.

This script does not alter any Sophiark score, ranking, or scientific engine.
It mirrors the active forward perturbation mechanics on an in-memory graph copy,
then writes audit-only artifacts under ``outputs/forensic_audits``.

Run from the repository root:
    python scripts/glp1r_pancreas_perturbation_independence_audit.py
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

# Support direct execution with ``python scripts/<script>.py`` from any cwd.
REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))

from src import main as motor
from src.biology_logic import _restore_edges, apply_biological_bonus, apply_functional_penalty
from src.config import FORCED_DISTANCE_SCALE
from src.core.analysis_runtime import symbol_options
from src.scientific.efficiency import strength_to_distance


@dataclass(frozen=True)
class SpearmanResult:
    comparison: str
    rho: float | None
    p_value: float | None
    n: int


def _finite_spearman(frame: pd.DataFrame, left: str, right: str, label: str) -> SpearmanResult:
    paired = frame[[left, right]].replace([np.inf, -np.inf], np.nan).dropna()
    if len(paired) < 3 or paired[left].nunique() < 2 or paired[right].nunique() < 2:
        return SpearmanResult(label, None, None, len(paired))
    rho, p_value = spearmanr(paired[left], paired[right])
    return SpearmanResult(label, float(rho), float(p_value), len(paired))


def _resolve_gene(symbol: str) -> str:
    options = symbol_options(False)
    matches = options[options["Symbol"].astype(str).str.upper() == symbol.upper()]
    if len(matches) != 1:
        raise RuntimeError(f"{symbol!r} için benzersiz protein kimliği bulunamadı; eşleşme sayısı={len(matches)}")
    return str(matches.iloc[0]["gene"])


def _apply_forward_perturbation(graph, target: str, *, attenuation: float, damping: float,
                                edge_evidence_weighting_mode: str) -> tuple[dict[str, float], dict[str, float], int]:
    """Run the same in-memory PageRank perturbation steps as run_infection_simulation.

    It intentionally returns every node rather than the engine's selected-report
    subset. Edge state is restored in ``finally`` and no graph/result is persisted.
    """
    names = [str(value) for value in graph.vs["name"]]
    name_to_index = {name: index for index, name in enumerate(names)}
    if target not in name_to_index:
        raise RuntimeError(f"Target {target} aktif Pancreas ağında bulunamadı.")

    bonus_snapshot: dict = {}
    normalized_edges: dict[int, float] = {}
    blocked_edges: dict[int, tuple[float, float]] = {}
    try:
        if edge_evidence_weighting_mode == "legacy_edge_modifiers":
            apply_biological_bonus([target], graph, snapshot=bonus_snapshot)
            apply_functional_penalty([target], graph, snapshot=bonus_snapshot)
        elif edge_evidence_weighting_mode != "structural_only":
            raise ValueError("edge_evidence_weighting_mode structural_only veya legacy_edge_modifiers olmalıdır")

        baseline_values = graph.pagerank(damping=damping, weights="weight", directed=False)
        baseline = dict(zip(names, map(float, baseline_values)))

        # The active engine normalizes unusually long forced-target edge distances
        # before attenuation. This affects only the in-memory forward calculation.
        target_index = name_to_index[target]
        distances = graph.es["distance"]
        mean_distance = float(np.mean(distances)) if distances else 1.0
        for neighbour_index in graph.neighbors(target_index):
            edge_id = graph.get_eid(target_index, neighbour_index, error=False)
            if edge_id == -1 or edge_id in normalized_edges:
                continue
            original_distance = graph.es[edge_id]["distance"]
            if original_distance > mean_distance * FORCED_DISTANCE_SCALE:
                normalized_edges[edge_id] = original_distance
                graph.es[edge_id]["distance"] = mean_distance * FORCED_DISTANCE_SCALE

        for neighbour_index in graph.neighbors(target_index):
            edge_id = graph.get_eid(target_index, neighbour_index, error=False)
            if edge_id == -1 or edge_id in blocked_edges:
                continue
            original_weight = graph.es[edge_id]["weight"]
            original_distance = graph.es[edge_id]["distance"]
            attenuated_weight = original_weight * attenuation
            graph.es[edge_id]["weight"] = attenuated_weight
            graph.es[edge_id]["distance"] = strength_to_distance(attenuated_weight)
            blocked_edges[edge_id] = (original_weight, original_distance)

        perturbed_values = graph.pagerank(damping=damping, weights="weight", directed=False)
        perturbed = dict(zip(names, map(float, perturbed_values)))
        return baseline, perturbed, len(blocked_edges)
    finally:
        for edge_id, (weight, distance) in blocked_edges.items():
            graph.es[edge_id]["weight"] = weight
            graph.es[edge_id]["distance"] = distance
        for edge_id, distance in normalized_edges.items():
            graph.es[edge_id]["distance"] = distance
        _restore_edges(graph, bonus_snapshot)


def _hop_bucket(distance: float) -> str:
    if distance == 1:
        return "direct target neighbor"
    if distance == 2:
        return "2-hop"
    if distance == 3:
        return "3-hop"
    if np.isfinite(distance):
        return "4+ hop"
    return "unreachable"


def run_audit(*, target_symbol: str, tissue: str, attenuation: float, damping: float,
              bc_sample_sources: int, edge_evidence_weighting_mode: str, output_root: Path) -> Path:
    target = _resolve_gene(target_symbol)
    graph, scores = motor.arayuz_icin_motoru_hazirla(
        forced_genes=[target], hedef_doku=tissue, bc_sample_sources=bc_sample_sources,
        tissue_normalization_mode="within_tissue",
    )
    if graph.vcount() == 0 or scores.empty:
        raise RuntimeError(f"{tissue} için aktif ağ hazırlanamadı.")

    active_graph = graph.copy()
    baseline, perturbed, attenuated_edges = _apply_forward_perturbation(
        active_graph, target, attenuation=attenuation, damping=damping,
        edge_evidence_weighting_mode=edge_evidence_weighting_mode,
    )
    ids = [str(value) for value in graph.vs["name"]]
    degree = dict(zip(ids, graph.degree()))
    target_index = graph.vs.find(name=target).index
    hops = dict(zip(ids, graph.distances(source=target_index, weights=None)[0]))
    symbols = dict(zip(symbol_options(False)["gene"].astype(str), symbol_options(False)["Symbol"].astype(str)))
    baseline_context = scores.drop_duplicates("gene").set_index("gene")

    frame = pd.DataFrame({
        "Gene": [symbols.get(gene, gene) for gene in ids],
        "Protein_ID": ids,
        "Baseline_PageRank": [baseline[gene] for gene in ids],
        "Baseline_Degree": [degree[gene] for gene in ids],
        "Hinterland_Skoru": [baseline_context["Hinterland_Skoru"].get(gene, np.nan) for gene in ids],
        "Perturbed_PageRank": [perturbed[gene] for gene in ids],
        "Delta_PageRank": [perturbed[gene] - baseline[gene] for gene in ids],
        "Delta_PageRank_Pct": [((perturbed[gene] - baseline[gene]) / baseline[gene] * 100.0) if baseline[gene] > 0 else np.nan for gene in ids],
        "Hop_Distance_From_GLP1R": [hops[gene] for gene in ids],
    })
    frame["Baseline_PageRank_Rank"] = frame["Baseline_PageRank"].rank(method="min", ascending=False).astype("Int64")
    frame["Perturbed_PageRank_Rank"] = frame["Perturbed_PageRank"].rank(method="min", ascending=False).astype("Int64")
    frame["Rank_Change"] = frame["Baseline_PageRank_Rank"] - frame["Perturbed_PageRank_Rank"]

    correlations = [
        _finite_spearman(frame, "Delta_PageRank_Pct", "Baseline_PageRank", "ΔPageRank_Pct vs Baseline PageRank"),
        _finite_spearman(frame, "Delta_PageRank_Pct", "Baseline_Degree", "ΔPageRank_Pct vs Baseline Degree"),
        _finite_spearman(frame, "Delta_PageRank_Pct", "Hinterland_Skoru", "ΔPageRank_Pct vs Hinterland_Skoru"),
    ]
    ranked = frame.replace([np.inf, -np.inf], np.nan).dropna(subset=["Delta_PageRank_Pct"]).sort_values(
        "Delta_PageRank_Pct", ascending=False, kind="mergesort"
    )
    top100 = ranked.head(100).copy()
    top100["Hop_Bucket"] = top100["Hop_Distance_From_GLP1R"].map(_hop_bucket)
    hop_distribution = top100["Hop_Bucket"].value_counts().reindex(
        ["direct target neighbor", "2-hop", "3-hop", "4+ hop", "unreachable"], fill_value=0
    )
    top25_columns = [
        "Gene", "Delta_PageRank_Pct", "Baseline_PageRank_Rank", "Perturbed_PageRank_Rank",
        "Rank_Change", "Baseline_Degree", "Hinterland_Skoru", "Hop_Distance_From_GLP1R",
    ]
    top25 = ranked.head(25)[top25_columns].copy()
    baseline_cutoff = frame["Hinterland_Skoru"].quantile(.25)
    low_baseline_top100 = int((top100["Hinterland_Skoru"] <= baseline_cutoff).sum())

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    output_dir = output_root / f"{target_symbol.lower()}_{tissue.lower().replace(' ', '_')}_{stamp}"
    output_dir.mkdir(parents=True, exist_ok=False)
    frame.to_csv(output_dir / "all_active_nodes.csv", index=False, encoding="utf-8-sig")
    top25.to_csv(output_dir / "top25_delta_pagerank_pct.csv", index=False, encoding="utf-8-sig")
    hop_distribution.rename_axis("Hop_Bucket").reset_index(name="Top100_Gene_Count").to_csv(
        output_dir / "top100_hop_distribution.csv", index=False, encoding="utf-8-sig"
    )
    summary = {
        "audit": "perturbation_independence",
        "target_symbol": target_symbol,
        "target_protein_id": target,
        "tissue": tissue,
        "active_network_nodes": graph.vcount(),
        "active_network_edges": graph.ecount(),
        "unique_attenuated_edges": attenuated_edges,
        "attenuation_fraction": attenuation,
        "pagerank_damping": damping,
        "edge_evidence_weighting_mode": edge_evidence_weighting_mode,
        "spearman": [asdict(item) for item in correlations],
        "top100_hop_distribution": hop_distribution.to_dict(),
        "top100_lowest_hinterland_quartile_count": low_baseline_top100,
        "baseline_hinterland_first_quartile_cutoff": float(baseline_cutoff),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2), encoding="utf-8")
    lines = [
        "# GLP1R · Pancreas perturbation independence audit", "",
        f"- Target: `{target}` ({target_symbol})", f"- Active network: {graph.vcount():,} nodes / {graph.ecount():,} edges",
        f"- Unique attenuated edges: {attenuated_edges:,}", "", "## Spearman correlations", "",
    ]
    lines.extend(f"- {item.comparison}: rho={item.rho if item.rho is not None else 'N/A'}, p={item.p_value if item.p_value is not None else 'N/A'}, n={item.n}" for item in correlations)
    lines.extend(["", "## Top-100 hop distribution", ""])
    lines.extend(f"- {bucket}: {count}" for bucket, count in hop_distribution.items())
    lines.extend(["", f"Top-100 içinde en düşük baseline Hinterland çeyreğinden kayıt: {low_baseline_top100}.", "",
                  "Yorum: korelasyonlar hub bağımlılığını ölçer; düşük korelasyon ve uzak-hop/low-baseline kayıtlar, sonucu basit hub re-ranking açıklamasından ayıran kanıttır."])
    (output_dir / "audit_summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return output_dir


def main() -> None:
    parser = argparse.ArgumentParser(description="Sophiark all-node perturbation independence audit")
    parser.add_argument("--target", default="GLP1R")
    parser.add_argument("--tissue", default="Pancreas")
    parser.add_argument("--attenuation", type=float, default=0.001)
    parser.add_argument("--damping", type=float, default=0.85)
    parser.add_argument("--bc-sample-sources", type=int, default=1200)
    parser.add_argument("--edge-evidence-weighting-mode", choices=["legacy_edge_modifiers", "structural_only"], default="legacy_edge_modifiers")
    parser.add_argument("--output-root", type=Path, default=Path("outputs") / "forensic_audits")
    args = parser.parse_args()
    output = run_audit(
        target_symbol=args.target, tissue=args.tissue, attenuation=args.attenuation, damping=args.damping,
        bc_sample_sources=args.bc_sample_sources, edge_evidence_weighting_mode=args.edge_evidence_weighting_mode,
        output_root=args.output_root,
    )
    print(f"Audit completed: {output}")


if __name__ == "__main__":
    main()
