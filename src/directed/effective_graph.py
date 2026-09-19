"""Directed-side replay of the immutable Classic effective-weight context."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import struct
from typing import Iterable

import igraph as ig
import numpy as np

from src.biology_logic import apply_biological_bonus, apply_functional_penalty


@dataclass(frozen=True, slots=True)
class EffectiveGraphSnapshot:
    graph: ig.Graph
    original_fingerprint_before: str
    original_fingerprint_after: str
    effective_fingerprint: str
    modifier_mode: str


def effective_graph_fingerprint(graph: ig.Graph) -> str:
    """Hash the complete node, pair, effective-weight and distance state."""

    if graph.is_directed():
        raise ValueError("Effective Classic source graph must be undirected")
    names = [str(value) for value in graph.vs["name"]]
    digest = hashlib.sha256()
    digest.update(b"sophiark-effective-graph-v1\0")
    for name in sorted(names):
        digest.update(name.encode("utf-8"))
        digest.update(b"\0")
    records = []
    has_distance = "distance" in graph.es.attributes()
    has_rescue = "rescue_bridge" in graph.es.attributes()
    for edge in graph.es:
        left, right = sorted((names[edge.source], names[edge.target]))
        records.append((
            left,
            right,
            float(edge["weight"]),
            float(edge["distance"]) if has_distance else float("nan"),
            bool(edge["rescue_bridge"]) if has_rescue else False,
        ))
    records.sort(key=lambda item: (item[0], item[1]))
    for left, right, weight, distance, rescue in records:
        digest.update(left.encode("utf-8"))
        digest.update(b"\0")
        digest.update(right.encode("utf-8"))
        digest.update(b"\0")
        digest.update(struct.pack("!dd?", weight, distance, rescue))
    return digest.hexdigest()


def prepare_effective_directed_source_graph(
    graph: ig.Graph,
    targets: Iterable[str],
    *,
    edge_evidence_weighting_mode: str = "legacy_edge_modifiers",
) -> EffectiveGraphSnapshot:
    """Copy ``graph`` and reuse Classic's existing modifier functions once.

    Classic calculation code and its graph are not changed.  The function calls
    the same Classic functions in the same order used immediately before the
    Classic baseline PageRank; no biological constants or formulas are copied.
    """

    before = effective_graph_fingerprint(graph)
    effective = graph.copy()
    unique_targets = list(dict.fromkeys(map(str, targets)))
    if edge_evidence_weighting_mode == "legacy_edge_modifiers":
        snapshot: dict[int, tuple[float, float]] = {}
        apply_biological_bonus(unique_targets, effective, snapshot=snapshot)
        apply_functional_penalty(unique_targets, effective, snapshot=snapshot)
    elif edge_evidence_weighting_mode != "structural_only":
        raise ValueError(
            "edge_evidence_weighting_mode structural_only veya legacy_edge_modifiers olmalıdır"
        )
    after = effective_graph_fingerprint(graph)
    if before != after:
        raise RuntimeError("Directed preparation mutated the immutable Classic source graph")
    return EffectiveGraphSnapshot(
        graph=effective,
        original_fingerprint_before=before,
        original_fingerprint_after=after,
        effective_fingerprint=effective_graph_fingerprint(effective),
        modifier_mode=edge_evidence_weighting_mode,
    )


def graph_weight_parity(left: ig.Graph, right: ig.Graph) -> dict[str, float | int]:
    """Graph-wide deterministic equality summary for undirected effective graphs."""

    def weights(graph: ig.Graph) -> tuple[set[str], dict[tuple[str, str], float]]:
        names = [str(value) for value in graph.vs["name"]]
        pairs = {
            tuple(sorted((names[edge.source], names[edge.target]))): float(edge["weight"])
            for edge in graph.es
        }
        return set(names), pairs

    left_nodes, left_pairs = weights(left)
    right_nodes, right_pairs = weights(right)
    shared = sorted(set(left_pairs) & set(right_pairs))
    absolute = np.asarray([abs(left_pairs[pair] - right_pairs[pair]) for pair in shared])
    relative = np.asarray([
        abs(left_pairs[pair] - right_pairs[pair]) / max(abs(left_pairs[pair]), abs(right_pairs[pair]), 1e-300)
        for pair in shared
    ])
    return {
        "left_nodes": len(left_nodes),
        "right_nodes": len(right_nodes),
        "node_mismatch_count": len(left_nodes ^ right_nodes),
        "left_pairs": len(left_pairs),
        "right_pairs": len(right_pairs),
        "pair_mismatch_count": len(set(left_pairs) ^ set(right_pairs)),
        "weight_mismatch_count": int(np.count_nonzero(absolute)),
        "max_absolute_weight_difference": float(absolute.max()) if absolute.size else 0.0,
        "max_relative_weight_difference": float(relative.max()) if relative.size else 0.0,
    }
