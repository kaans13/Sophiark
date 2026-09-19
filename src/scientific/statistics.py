"""Çoklu test düzeltmesi ve perturbasyon null modeli."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import igraph as ig
import numpy as np

from .efficiency import strength_to_distance


def benjamini_hochberg(p_values: Sequence[float]) -> np.ndarray:
    """Benjamini-Hochberg FDR; NaN değerleri 1.0 olarak ele alır."""
    p = np.asarray(p_values, dtype=float)
    p = np.where(np.isfinite(p), np.clip(p, 0.0, 1.0), 1.0)
    if p.size == 0:
        return p
    order = np.argsort(p, kind="mergesort")
    ranked = p[order]
    adjusted = ranked * p.size / np.arange(1, p.size + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.clip(adjusted, 0.0, 1.0)
    return result


def _degree_matched_vertices(
    graph: ig.Graph,
    target_indices: Sequence[int],
    rng: np.random.Generator,
) -> list[int]:
    degrees = np.asarray(graph.degree(), dtype=int)
    chosen: list[int] = []
    forbidden = set(target_indices)
    for target_index in target_indices:
        degree = degrees[target_index]
        # Log2-degree bandı, hub ve düşük dereceli düğümlerin aynı null havuzuna
        # düşmesini önler. Havuz küçükse en yakın dereceye deterministik genişler.
        band = int(np.floor(np.log2(max(1, degree))))
        candidates = np.flatnonzero(
            np.floor(np.log2(np.maximum(1, degrees))).astype(int) == band
        )
        candidates = np.asarray(
            [idx for idx in candidates if idx not in forbidden and idx not in chosen],
            dtype=int,
        )
        if candidates.size == 0:
            order = np.argsort(np.abs(degrees - degree), kind="mergesort")
            candidates = np.asarray(
                [idx for idx in order if idx not in forbidden and idx not in chosen], dtype=int
            )
        if candidates.size:
            chosen.append(int(rng.choice(candidates)))
    return chosen


def empirical_pvalues_degree_matched(
    graph: ig.Graph,
    baseline_pagerank: Mapping[str, float],
    observed_abs_delta: Mapping[str, float],
    target_names: Sequence[str],
    *,
    attenuation_fraction: float,
    damping: float,
    iterations: int,
    random_seed: int,
) -> tuple[dict[str, float], dict[str, object]]:
    """Degree-matched hedef perturbasyonlarından düğüm-bazlı empirical p üret.

    Null, aynı graf üzerinde hedeflerle benzer derece bandındaki düğümlerin
    kenarlarını aynı attenuation fraction ile zayıflatır. Her düğüm için
    ``|ΔPageRank|`` gözlenen değeriyle karşılaştırılır. Eksik/0 iterasyon
    anlamlılık iddiası üretmez (p=1).
    """
    names = list(map(str, graph.vs["name"]))
    name_to_index = {name: index for index, name in enumerate(names)}
    target_indices = [name_to_index[name] for name in target_names if name in name_to_index]
    observed = np.asarray([float(observed_abs_delta.get(name, 0.0)) for name in names])
    if iterations <= 0 or not target_indices:
        return ({name: 1.0 for name in names}, {
            "method": "degree_matched_edge_attenuation",
            "iterations": 0,
            "random_seed": random_seed,
            "degree_matching": "log2 degree bands",
        })

    baseline = np.asarray([float(baseline_pagerank.get(name, 0.0)) for name in names])
    exceedances = np.zeros(len(names), dtype=np.int32)
    rng = np.random.default_rng(random_seed)
    for _ in range(iterations):
        null_targets = _degree_matched_vertices(graph, target_indices, rng)
        snapshots: dict[int, tuple[float, float]] = {}
        for vertex_index in null_targets:
            for edge_id in graph.incident(vertex_index):
                if edge_id in snapshots:
                    continue
                weight = float(graph.es[edge_id]["weight"])
                distance = float(graph.es[edge_id]["distance"])
                new_weight = weight * attenuation_fraction
                snapshots[edge_id] = (weight, distance)
                graph.es[edge_id]["weight"] = new_weight
                graph.es[edge_id]["distance"] = strength_to_distance(new_weight)
        null_pr = np.asarray(graph.pagerank(damping=damping, weights="weight", directed=False))
        null_abs = np.abs(null_pr - baseline)
        exceedances += null_abs >= (observed - np.finfo(float).eps)
        for edge_id, (weight, distance) in snapshots.items():
            graph.es[edge_id]["weight"] = weight
            graph.es[edge_id]["distance"] = distance

    p_values = (exceedances + 1.0) / (iterations + 1.0)
    return dict(zip(names, map(float, p_values))), {
        "method": "degree_matched_edge_attenuation",
        "iterations": iterations,
        "random_seed": random_seed,
        "degree_matching": "log2 degree bands",
    }
