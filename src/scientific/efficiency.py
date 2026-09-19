"""Weighted global/local efficiency hesapları.

STRING ve doku ağırlıkları interaction strength olarak yorumlanır. En kısa yol
algoritmalarında tek merkezi dönüşüm ``distance = 1 / strength`` kullanılır.
Sıfır/negatif/non-finite strength geçilemez kenardır ve ``inf`` distance alır.
Disconnected düğüm çiftleri global efficiency toplamına 0 katkı verir.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Sequence

import igraph as ig
import numpy as np


DISTANCE_EPSILON = 1e-12


def strength_to_distance(strength: float) -> float:
    value = float(strength)
    if not np.isfinite(value) or value <= 0:
        return float("inf")
    return 1.0 / max(value, DISTANCE_EPSILON)


def strengths_to_distances(strengths) -> np.ndarray:
    """Vectorized form of :func:`strength_to_distance`."""
    values = np.asarray(strengths, dtype=float)
    result = np.full(values.shape, np.inf, dtype=float)
    valid = np.isfinite(values) & (values > 0)
    result[valid] = 1.0 / np.maximum(values[valid], DISTANCE_EPSILON)
    return result


def synchronize_distances(graph: ig.Graph) -> None:
    graph.es["distance"] = [strength_to_distance(value) for value in graph.es["weight"]]


def _source_indices(vertex_count: int, sample_size: int | None, seed: int) -> list[int]:
    if sample_size is None or sample_size >= vertex_count:
        return list(range(vertex_count))
    if sample_size <= 0:
        return []
    return sorted(np.random.default_rng(seed).choice(vertex_count, size=sample_size, replace=False).tolist())


def weighted_global_efficiency(
    graph: ig.Graph,
    *,
    sample_sources: int | None = None,
    random_seed: int = 42,
) -> tuple[float, dict[str, object]]:
    n = graph.vcount()
    if n < 2:
        return 0.0, {"mode": "exact", "sample_sources": n, "graph_size": n}
    sources = _source_indices(n, sample_sources, random_seed)
    if not sources:
        return 0.0, {"mode": "invalid_empty_sample", "sample_sources": 0, "graph_size": n}
    matrix = np.asarray(graph.distances(source=sources, weights="distance"), dtype=float)
    efficiencies = np.zeros_like(matrix)
    finite_positive = np.isfinite(matrix) & (matrix > 0)
    efficiencies[finite_positive] = 1.0 / matrix[finite_positive]
    # Her source'un kendisi denominator dışında tutulur. Sampling yalnızca source
    # eksenindedir; her source bütün hedef düğümlere ölçülür.
    value = float(efficiencies.sum() / (len(sources) * (n - 1)))
    return value, {
        "mode": "exact" if len(sources) == n else "deterministic_source_sample",
        "sample_sources": len(sources),
        "graph_size": n,
        "random_seed": random_seed,
        "weighted": True,
        "distance_transform": "1/strength; non-positive edges unreachable",
    }


def _neighborhood_efficiency(graph: ig.Graph, vertex_index: int) -> float:
    neighbors = sorted(set(graph.neighbors(vertex_index)))
    if len(neighbors) < 2:
        return 0.0
    neighborhood = graph.induced_subgraph(neighbors)
    value, _ = weighted_global_efficiency(neighborhood, sample_sources=None)
    return value


def mean_local_efficiency(
    graph: ig.Graph,
    *,
    sample_nodes: int | None = None,
    random_seed: int = 42,
) -> tuple[float, dict[str, object]]:
    indices = _source_indices(graph.vcount(), sample_nodes, random_seed)
    values = [_neighborhood_efficiency(graph, index) for index in indices]
    return (float(np.mean(values)) if values else 0.0), {
        "mode": "exact" if len(indices) == graph.vcount() else "deterministic_node_sample",
        "sample_nodes": len(indices),
        "graph_size": graph.vcount(),
        "random_seed": random_seed,
        "weighted": True,
    }


def target_neighborhood_efficiency(
    graph: ig.Graph,
    target_names: Sequence[str],
    *,
    sample_sources: int | None = None,
    random_seed: int = 42,
) -> tuple[float, dict[str, object]]:
    """Çoklu hedefte hedeflerin kapalı 1-hop ego ağlarının birleşimi.

    Hedefler ve tüm birinci derece komşuları alt-grafa dahildir. Bu, standart
    düğüm-local-efficiency tanımından farklı olarak perturb edilen hedef
    kenarlarını da ölçer; bu nedenle raporda açıkça ``closed_1hop_union`` diye
    etiketlenir. Çoklu hedef birleşimi sıralanmış düğüm indeksleriyle
    deterministiktir.
    """
    name_to_index = {str(v["name"]): v.index for v in graph.vs}
    target_indices = [name_to_index[name] for name in target_names if name in name_to_index]
    neighbors: set[int] = set(target_indices)
    for index in target_indices:
        neighbors.update(graph.neighbors(index))
    ordered = sorted(neighbors)
    if len(ordered) < 2:
        return 0.0, {"mode": "exact", "neighborhood_nodes": len(ordered), "definition": "closed_1hop_union_including_targets"}
    neighborhood = graph.induced_subgraph(ordered)
    value, metadata = weighted_global_efficiency(
        neighborhood, sample_sources=sample_sources, random_seed=random_seed
    )
    metadata.update({
        "neighborhood_nodes": len(ordered),
        "definition": "closed_1hop_union_including_targets",
        "target_count": len(target_indices),
    })
    return value, metadata


@dataclass(frozen=True)
class EfficiencyResult:
    global_efficiency: float
    mean_local_efficiency: float
    target_local_efficiency: float
    metadata: dict[str, object]


def graph_efficiency_summary(
    graph: ig.Graph,
    target_names: Sequence[str],
    *,
    global_sample_sources: int | None = 256,
    local_sample_nodes: int | None = 128,
    random_seed: int = 42,
) -> EfficiencyResult:
    global_value, global_meta = weighted_global_efficiency(
        graph, sample_sources=global_sample_sources, random_seed=random_seed
    )
    local_value, local_meta = mean_local_efficiency(
        graph, sample_nodes=local_sample_nodes, random_seed=random_seed
    )
    target_value, target_meta = target_neighborhood_efficiency(
        graph, target_names, sample_sources=global_sample_sources, random_seed=random_seed
    )
    return EfficiencyResult(
        global_efficiency=global_value,
        mean_local_efficiency=local_value,
        target_local_efficiency=target_value,
        metadata={"global": global_meta, "mean_local": local_meta, "target_local": target_meta},
    )


def percent_change(baseline: float, perturbed: float) -> float:
    return ((perturbed - baseline) / baseline * 100.0) if baseline else 0.0
