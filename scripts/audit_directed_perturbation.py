"""Synthetic decision audit for directed target-loss perturbation semantics."""

from __future__ import annotations

import json
from pathlib import Path

import igraph as ig
import numpy as np


GRAPHS = {
    "one_way_chain": [("A", "B"), ("B", "C"), ("C", "D")],
    "competing_route": [("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")],
    "alternative_source": [("A", "C"), ("B", "C"), ("C", "D")],
    "feedback": [("A", "B"), ("B", "C"), ("C", "A")],
    "acceptance_alternatives": [
        ("U", "T"), ("U", "X"), ("T", "D"), ("T", "G"),
        ("X", "E"), ("D", "F"), ("G", "F"), ("E", "F"),
    ],
}
TARGETS = {
    "one_way_chain": ("A", "B", "C", "D"),
    "competing_route": ("A", "B", "C"),
    "alternative_source": ("A",),
    "feedback": ("A", "B", "C"),
    "acceptance_alternatives": ("T",),
}


def _graph(edges) -> ig.Graph:
    names = sorted({name for pair in edges for name in pair})
    index = {name: position for position, name in enumerate(names)}
    graph = ig.Graph(n=len(names), edges=[(index[a], index[b]) for a, b in edges], directed=True)
    graph.vs["name"] = names
    graph.es["weight"] = [1.0] * graph.ecount()
    return graph


def _pagerank(graph: ig.Graph) -> dict[str, float]:
    values = graph.pagerank(damping=0.85, weights="weight", directed=True, implementation="prpack")
    return dict(zip(map(str, graph.vs["name"]), map(float, values)))


def _attenuate(graph: ig.Graph, target: str, mode: str, factor: float = 0.001) -> ig.Graph:
    result = graph.copy()
    index = result.vs.find(name=target).index
    edge_ids = set(result.incident(index, mode=mode))
    for edge_id in edge_ids:
        result.es[edge_id]["weight"] = float(result.es[edge_id]["weight"]) * factor
    return result


def _remove(graph: ig.Graph, target: str) -> ig.Graph:
    result = graph.copy()
    result.delete_vertices(result.vs.find(name=target).index)
    return result


def _measure(graph: ig.Graph, target: str, baseline: dict[str, float], after: dict[str, float]) -> dict:
    names = list(baseline)
    target_index = graph.vs.find(name=target).index
    downstream = {
        graph.vs[index]["name"] for index in graph.subcomponent(target_index, mode="out")
    } - {target}
    upstream = {
        graph.vs[index]["name"] for index in graph.subcomponent(target_index, mode="in")
    } - {target}
    common = [name for name in names if name in after]
    baseline_rank = {name: rank for rank, name in enumerate(sorted(common, key=lambda item: abs(baseline[item]), reverse=True), 1)}
    after_rank = {name: rank for rank, name in enumerate(sorted(common, key=lambda item: abs(after[item]), reverse=True), 1)}
    return {
        "target_baseline": baseline[target],
        "target_perturbed": after.get(target),
        "target_change": None if target not in after else after[target] - baseline[target],
        "downstream_absolute_change": sum(abs(after[name] - baseline[name]) for name in downstream if name in after),
        "upstream_absolute_change": sum(abs(after[name] - baseline[name]) for name in upstream if name in after),
        "total_redistribution_l1": sum(abs(after[name] - baseline[name]) for name in common),
        "median_absolute_rank_shift": float(np.median([abs(after_rank[name] - baseline_rank[name]) for name in common])),
        "source": graph.indegree(target_index) == 0,
        "sink": graph.outdegree(target_index) == 0,
    }


def run_audit() -> dict:
    result = {"suppression_factor": 0.001, "damping": 0.85, "graphs": {}}
    for graph_name, edges in GRAPHS.items():
        graph = _graph(edges)
        baseline = _pagerank(graph)
        target_results = {}
        for target in TARGETS[graph_name]:
            strategies = {
                "incident_edge_attenuation": _pagerank(_attenuate(graph, target, "all")),
                "incoming_only_attenuation": _pagerank(_attenuate(graph, target, "in")),
                "outgoing_only_attenuation": _pagerank(_attenuate(graph, target, "out")),
                "effective_node_knockout": _pagerank(_remove(graph, target)),
            }
            target_results[target] = {
                strategy: _measure(graph, target, baseline, after)
                for strategy, after in strategies.items()
            }
            target_results[target]["incident_equals_incoming_for_pagerank"] = all(
                np.isclose(
                    strategies["incident_edge_attenuation"][name],
                    strategies["incoming_only_attenuation"][name], rtol=1e-12, atol=1e-14,
                ) for name in baseline
            )
            target_results[target]["outgoing_only_is_normalized_away"] = all(
                np.isclose(strategies["outgoing_only_attenuation"][name], baseline[name], rtol=1e-12, atol=1e-14)
                for name in baseline
            )
        result["graphs"][graph_name] = target_results
    result["decision"] = {
        "selected_production_strategy": "incident_edge_attenuation",
        "effective_pagerank_mechanism": "reduced incoming transition probability from upstream nodes",
        "outgoing_scaling": "normalized away when every target outgoing arc is scaled equally",
        "why_selected": "It is the exact directed edge-level analogue of the immutable classic incident-edge attenuation contract.",
        "incoming_only": "PageRank-equivalent, but does not preserve the classic incident-edge contract for future edge-sensitive metrics.",
        "outgoing_only": "Rejected: weighted PageRank row normalization makes uniform outgoing scaling ineffective.",
        "node_removal": "Rejected for v1: changes the node universe and is materially harsher than the classic reference.",
        "transition_influence_suppression": "Deferred: requires an explicitly substochastic/custom propagation model, not standard PageRank.",
        "known_limitation": "A pure directed source with no incoming arcs cannot lose PageRank accessibility under incident attenuation; teleportation remains.",
    }
    return result


if __name__ == "__main__":
    output = Path(__file__).resolve().parents[1] / "docs" / "directed_perturbation_semantics_audit.json"
    output.write_text(json.dumps(run_audit(), indent=2), encoding="utf-8")
    print(output)
