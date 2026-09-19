"""Direction-aware PageRank, perturbation, redistribution and BC."""

from __future__ import annotations

import time
from dataclasses import asdict
from typing import Iterable, Mapping

import igraph as ig
import numpy as np
import pandas as pd

from src.scientific.efficiency import strength_to_distance

from .models import (
    DirectedCalculationResult, DirectedEngineStatus, HybridDirectedGraph,
    PerturbationStrategy,
)


def compute_directed_pagerank(
    graph: ig.Graph, *, damping: float = 0.85, max_iter: int = 200,
    tolerance: float = 1e-7,
) -> dict[str, float]:
    """Use igraph's standard dangling-node handling on the directed topology.

    ``max_iter`` and ``tolerance`` are retained in the public contract and run
    metadata.  PRPACK, also used by the classic motor, owns convergence and
    does not expose these knobs through python-igraph.
    """
    if not graph.is_directed():
        raise ValueError("Directed PageRank requires a directed graph")
    if not 0.0 < damping < 1.0:
        raise ValueError("damping must be between 0 and 1")
    if max_iter <= 0 or tolerance <= 0:
        raise ValueError("max_iter and tolerance must be positive")
    np.random.seed(42)
    try:
        values = graph.pagerank(
            damping=damping, weights="weight", directed=True, implementation="prpack",
        )
    except Exception:
        values = graph.pagerank(
            damping=damping, weights="weight", directed=True, implementation="arpack",
        )
    return {str(graph.vs[index]["name"]): float(value) for index, value in enumerate(values)}


def suppress_target_edges(
    graph: ig.Graph, targets: Iterable[str], *, fraction: float,
    strategy: PerturbationStrategy = PerturbationStrategy.INCIDENT_EDGE_ATTENUATION,
) -> tuple[ig.Graph, tuple[int, ...]]:
    if not 0.0 <= fraction <= 1.0:
        raise ValueError("fraction must be between 0 and 1")
    strategy = PerturbationStrategy(strategy)
    perturbed = graph.copy()
    name_to_index = {str(v["name"]): v.index for v in perturbed.vs}
    edge_ids: set[int] = set()
    for target in dict.fromkeys(map(str, targets)):
        index = name_to_index.get(target)
        if index is None:
            continue
        if strategy is PerturbationStrategy.INCIDENT_EDGE_ATTENUATION:
            edge_ids.update(perturbed.incident(index, mode="all"))
    for edge_id in sorted(edge_ids):
        new_weight = float(perturbed.es[edge_id]["weight"]) * fraction
        perturbed.es[edge_id]["weight"] = new_weight
        perturbed.es[edge_id]["distance"] = strength_to_distance(new_weight)
    return perturbed, tuple(sorted(edge_ids))


def compute_directed_betweenness(
    graph: ig.Graph, *, sample_sources: int | None, random_seed: int,
) -> tuple[np.ndarray, dict[str, object]]:
    n = graph.vcount()
    if sample_sources is None or sample_sources >= n:
        values = np.asarray(graph.betweenness(weights="distance", directed=True), dtype=float)
        return values, {"mode": "full", "sample_size": n}
    if sample_sources <= 0:
        return np.zeros(n), {"mode": "not_computed", "sample_size": 0}
    sources = sorted(np.random.default_rng(random_seed).choice(n, size=sample_sources, replace=False).tolist())
    values = np.asarray(
        graph.betweenness(weights="distance", directed=True, sources=sources), dtype=float,
    )
    values *= n / sample_sources
    return values, {"mode": "sampled", "sample_size": sample_sources}


def _relation(
    graph: ig.Graph, targets: set[str], entity: str, names: dict[str, int],
) -> tuple[str, str]:
    entity_index = names.get(entity)
    if entity_index is None:
        return "not_in_graph", "none"
    outgoing = incoming = supported = fallback = False
    for target in targets:
        target_index = names.get(target)
        if target_index is None:
            continue
        for source, dest, label in ((target_index, entity_index, "out"), (entity_index, target_index, "in")):
            edge_id = graph.get_eid(source, dest, directed=True, error=False)
            if edge_id == -1:
                continue
            outgoing |= label == "out"
            incoming |= label == "in"
            is_fallback = bool(graph.es[edge_id]["fallback_bidirectional"])
            fallback |= is_fallback
            supported |= not is_fallback
    relation = "bidirectional" if outgoing and incoming else ("target_to_entity" if outgoing else ("entity_to_target" if incoming else "no_direct_relation"))
    support = "OmniPath" if supported else ("undirected_fallback" if fallback else "none")
    return relation, support


def _entity_direction_context(
    graph: ig.Graph, *, target_indices: list[int], entity_index: int,
    directed_distance: int | None,
) -> tuple[int | None, str, float]:
    incident = set(graph.incident(entity_index, mode="all"))
    if incident:
        supported = sum(not bool(graph.es[edge_id]["fallback_bidirectional"]) for edge_id in incident)
        coverage = supported / len(incident)
    else:
        coverage = 0.0
    signs: set[str] = set()
    for target_index in target_indices:
        for source, target in ((target_index, entity_index), (entity_index, target_index)):
            edge_id = graph.get_eid(source, target, directed=True, error=False)
            if edge_id != -1 and not graph.es[edge_id]["fallback_bidirectional"]:
                signs.update(filter(None, str(graph.es[edge_id]["sign_metadata"]).split(";")))
    return directed_distance, (";".join(sorted(signs)) if signs else "not_available"), coverage


def run_directed_calculation(
    hybrid: HybridDirectedGraph,
    *, targets: Iterable[str], block_weight_fraction: float = 0.001,
    damping: float = 0.85, max_iter: int = 200, tolerance: float = 1e-7,
    bc_sample_sources: int | None = 4, random_seed: int = 42,
    strategy: PerturbationStrategy = PerturbationStrategy.INCIDENT_EDGE_ATTENUATION,
    precomputed_directed_bc: tuple[np.ndarray, dict[str, object]] | None = None,
    source_graph_metadata: Mapping[str, object] | None = None,
) -> DirectedCalculationResult:
    started = time.perf_counter()
    targets = tuple(dict.fromkeys(map(str, targets)))
    present = set(map(str, hybrid.graph.vs["name"]))
    resolved_targets = tuple(target for target in targets if target in present)
    if not resolved_targets:
        return DirectedCalculationResult(
            status=DirectedEngineStatus.UNAVAILABLE,
            directed_graph=hybrid, error="No requested target exists in the directed graph",
        )
    t0 = time.perf_counter()
    baseline = compute_directed_pagerank(
        hybrid.graph, damping=damping, max_iter=max_iter, tolerance=tolerance,
    )
    t1 = time.perf_counter()
    perturbed_graph, suppressed = suppress_target_edges(
        hybrid.graph, resolved_targets, fraction=block_weight_fraction, strategy=strategy,
    )
    perturbed = compute_directed_pagerank(
        perturbed_graph, damping=damping, max_iter=max_iter, tolerance=tolerance,
    )
    t2 = time.perf_counter()
    bc, bc_meta = precomputed_directed_bc or compute_directed_betweenness(
        hybrid.graph, sample_sources=bc_sample_sources, random_seed=random_seed,
    )
    t3 = time.perf_counter()
    target_set = set(resolved_targets)
    name_to_index = {str(v["name"]): v.index for v in hybrid.graph.vs}
    target_indices = [name_to_index[target] for target in resolved_targets]
    target_distances = np.asarray(
        hybrid.graph.distances(source=target_indices, mode="out"), dtype=float,
    )
    minimum_distances = np.min(target_distances, axis=0)
    rows = []
    for index, entity in enumerate(map(str, hybrid.graph.vs["name"])):
        if entity in target_set:
            continue
        before, after = baseline[entity], perturbed[entity]
        delta = after - before
        relation, support = _relation(hybrid.graph, target_set, entity, name_to_index)
        distance_value = minimum_distances[index]
        directed_distance, sign_context, direction_coverage = _entity_direction_context(
            hybrid.graph, target_indices=target_indices, entity_index=index,
            directed_distance=(int(distance_value) if np.isfinite(distance_value) else None),
        )
        rows.append({
            "entity": entity,
            "gene": entity,
            "Directed_PageRank_Baseline": before,
            "Directed_PageRank_Perturbed": after,
            "Directed_Delta_PageRank": delta,
            "Directed_Redistribution_Pct": (delta / before * 100.0 if before > 0 else np.nan),
            "Directed_Abs_Delta_PageRank": abs(delta),
            "Directed_Response_Direction": "Influence Gain" if delta > 0 else ("Influence Loss" if delta < 0 else "No Change"),
            "Directed_BC": float(bc[index]),
            "Direction_Relation": relation,
            "Direction_Support_Type": support,
            "Directed_Hop_Distance": directed_distance,
            # Temporary export compatibility alias; canonical semantics are hops.
            "Directed_Distance": directed_distance,
            "Sign_Context": sign_context,
            "Direction_Source_Coverage": direction_coverage,
            "Alternative_Route_Status": "not_evaluated_in_calculation_engine",
        })
    report = pd.DataFrame(rows).sort_values(
        "Directed_Abs_Delta_PageRank", ascending=False, kind="mergesort",
    ).reset_index(drop=True)
    timings = {
        "baseline_pagerank_seconds": t1 - t0,
        "perturbation_and_pagerank_seconds": t2 - t1,
        "directed_bc_seconds": (0.0 if precomputed_directed_bc is not None else t3 - t2),
        "redistribution_seconds": time.perf_counter() - t3,
        "total_seconds": time.perf_counter() - started,
    }
    return DirectedCalculationResult(
        status=DirectedEngineStatus.AVAILABLE, report=report, directed_graph=hybrid,
        metadata={
            **dict(source_graph_metadata or {}),
            "engine": "directed", "graph_mode": "hybrid_string_omnipath",
            "direction_coverage": hybrid.provenance.direction_coverage,
            "direction_dataset_version": hybrid.provenance.omnipath_fingerprint,
            "fallback_edge_count": hybrid.provenance.fallback_string_edges,
            "perturbation_strategy": PerturbationStrategy(strategy).value,
            "directed_perturbation_strategy": PerturbationStrategy(strategy).value,
            "suppression_factor": block_weight_fraction,
            "suppressed_arc_count": len(suppressed), "targets": resolved_targets,
            "target_pagerank": {
                target: {"baseline": baseline[target], "perturbed": perturbed[target],
                         "delta": perturbed[target] - baseline[target]}
                for target in resolved_targets
            },
            "damping": damping, "max_iter": max_iter, "tolerance": tolerance,
            "bc": bc_meta,
            "directed_distance_semantics": "minimum_unweighted_outgoing_hops",
            "directed_distance_canonical_field": "Directed_Hop_Distance",
            "directed_distance_legacy_alias": "Directed_Distance",
            "unsupported_metrics": (
                "signed_pagerank", "directed_hinterland", "directed_efficiency", "directed_null_model",
            ),
            "hybrid_graph": asdict(hybrid.provenance),
        }, timings=timings,
    )


def safe_run_directed_calculation(*args, **kwargs) -> DirectedCalculationResult:
    try:
        return run_directed_calculation(*args, **kwargs)
    except Exception as error:
        return DirectedCalculationResult(
            status=DirectedEngineStatus.UNAVAILABLE,
            error=f"{type(error).__name__}: {error}",
            metadata={"engine": "directed", "classic_engine_affected": False},
        )
