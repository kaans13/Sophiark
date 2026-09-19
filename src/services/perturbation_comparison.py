"""Reusable, non-destructive perturbation comparison primitives."""

from __future__ import annotations

from pathlib import Path
from typing import Callable
from uuid import uuid4

import pandas as pd

from src.models import PerturbationComparison
from src.scientific.efficiency import strength_to_distance


def _two_hop_vertices(graph, seed: int, maximum: int = 300) -> list[int]:
    selected, frontier = {seed}, {seed}
    for _ in range(2):
        next_frontier = set()
        for vertex in frontier:
            next_frontier.update(graph.neighbors(vertex))
        new_vertices = next_frontier - selected
        selected.update(new_vertices)
        if len(selected) >= maximum:
            break
        frontier = new_vertices
    return sorted(selected)[:maximum]


def _projected_target_local_bc(graph, intervention_gene: str, observed_gene: str, block_strength: float) -> float:
    """Explicit 2-hop BC projection, not a replacement for the full BC metric."""
    name_to_idx = {name: index for index, name in enumerate(graph.vs["name"])}
    if intervention_gene not in name_to_idx or observed_gene not in name_to_idx:
        return float("nan")
    observed_idx = name_to_idx[observed_gene]
    before_ids = _two_hop_vertices(graph, observed_idx)
    before = graph.induced_subgraph(before_ids)
    before_target = before.vs.find(name=observed_gene).index
    before_bc = float(before.betweenness(vertices=[before_target], weights="distance", directed=False)[0])

    projected = graph.copy()
    intervention_idx = name_to_idx[intervention_gene]
    for neighbor in projected.neighbors(intervention_idx):
        edge_id = projected.get_eid(intervention_idx, neighbor, error=False)
        if edge_id < 0:
            continue
        weight = float(projected.es[edge_id]["weight"])
        new_weight = weight * block_strength
        projected.es[edge_id]["weight"] = new_weight
        projected.es[edge_id]["distance"] = strength_to_distance(new_weight)
    after_ids = _two_hop_vertices(projected, observed_idx)
    after = projected.induced_subgraph(after_ids)
    after_target = after.vs.find(name=observed_gene).index
    after_bc = float(after.betweenness(vertices=[after_target], weights="distance", directed=False)[0])
    return after_bc - before_bc


def full_target_bc_delta(graph, intervention_gene: str, observed_gene: str, block_strength: float) -> float:
    """Tüm kaynak/hedef yollarını kullanan weighted target BC değişimi."""
    name_to_idx = {name: index for index, name in enumerate(graph.vs["name"])}
    if intervention_gene not in name_to_idx or observed_gene not in name_to_idx:
        return float("nan")
    target_index = name_to_idx[observed_gene]
    before = float(graph.betweenness(vertices=[target_index], weights="distance", directed=False)[0])
    projected = graph.copy()
    intervention_index = name_to_idx[intervention_gene]
    for edge_id in projected.incident(intervention_index):
        weight = float(projected.es[edge_id]["weight"])
        new_weight = weight * block_strength
        projected.es[edge_id]["weight"] = new_weight
        projected.es[edge_id]["distance"] = strength_to_distance(new_weight)
    after = float(projected.betweenness(vertices=[target_index], weights="distance", directed=False)[0])
    return after - before


def run_single_intervention(
    *, motor_module, graph, scores: pd.DataFrame, intervention_gene: str,
    observed_gene: str, output_root: Path, block_strength: float, damping: float,
    ghost_filter: Callable[[pd.DataFrame, object], tuple[pd.DataFrame, list]],
) -> PerturbationComparison:
    """Run the unchanged forward engine on a graph copy and retain raw deltas.

    A unique scratch directory prevents a candidate search from overwriting the
    user's normal forward-simulation report.
    """
    scratch = output_root / "target_stress" / uuid4().hex
    scratch.mkdir(parents=True, exist_ok=True)
    result = motor_module.run_infection_simulation(
        graph.copy(), scores, scratch, spesifik_hedefler=[intervention_gene],
        block_weight_fraction=block_strength, damping=damping,
        return_comparison_metrics=True, comparison_genes=[observed_gene],
        # Candidate taramasındaki her aday için pahalı full scientific paneli
        # tekrarlama. Bu çağrı validation hedefi üzerindeki ölçülen ΔPR/yerel
        # değişimi üretir; nihai forward analiz null/FDR ve efficiency'yi ayrıca
        # hesaplar.
        redistribution_mode="top_n", exploratory_top_n=250,
        null_iterations=0, compute_structural_metrics=False,
    )
    report, raw = result
    report, _ = ghost_filter(report, graph)
    projected_delta_bc = _projected_target_local_bc(
        graph, intervention_gene, observed_gene, block_strength
    )
    observed = raw["observed"].setdefault(observed_gene, {})
    observed["delta_local_bc_2hop"] = projected_delta_bc
    return PerturbationComparison(
        intervention_genes=(intervention_gene,), observed_genes=(observed_gene,),
        system_shift_pct=float(raw["system_shift_pct"]),
        local_stress_pct=float(raw["local_stress_pct"]),
        blocked_edges=int(raw["blocked_edges"]),
        observed_metrics=raw["observed"], report=report,
        top_pagerank_gains=tuple(raw.get("top_pagerank_gains", [])),
    )
