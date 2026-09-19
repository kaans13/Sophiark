"""Build a directed calculation topology from STRING edges plus OmniPath direction."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict

import igraph as ig

from src.scientific.efficiency import strength_to_distance
from src.signaling.models import SignalingDataset, SignalingEdge

from .models import HybridDirectedGraph, HybridGraphProvenance


DIRECTION_POLICY = "omnipath_direction_when_available"
FALLBACK_POLICY = "bidirectional_equal_weight"
SYMMETRIC_CONTROL_POLICY = "force_bidirectional_fallback_control"


def directed_graph_cache_key(
    *, tissue_context: str, string_dataset_version: str,
    omnipath_fingerprint: str, direction_policy: str = DIRECTION_POLICY,
    fallback_policy: str = FALLBACK_POLICY, taxon_id: int = 9606,
    effective_graph_fingerprint: str = "unspecified",
    edge_evidence_weighting_mode: str = "legacy_edge_modifiers",
) -> str:
    payload = json.dumps({
        "tissue_context": tissue_context,
        "string_dataset_version": string_dataset_version,
        "omnipath_fingerprint": omnipath_fingerprint,
        "direction_policy": direction_policy,
        "fallback_policy": fallback_policy,
        "taxon_id": int(taxon_id),
        "effective_graph_fingerprint": effective_graph_fingerprint,
        "edge_evidence_weighting_mode": edge_evidence_weighting_mode,
    }, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def _layer1_edge(edge: SignalingEdge, node_layer1: dict[str, str]) -> tuple[str, str] | None:
    source = node_layer1.get(edge.source)
    target = node_layer1.get(edge.target)
    if not source or not target or source == target:
        return None
    return source, target


def build_hybrid_directed_graph(
    string_graph: ig.Graph,
    signaling_dataset: SignalingDataset,
    *,
    tissue_context: str = "None",
    string_dataset_version: str = "local_string",
    force_symmetric_fallback: bool = False,
    effective_graph_fingerprint: str = "unspecified",
    taxon_id: int = 9606,
    edge_evidence_weighting_mode: str = "legacy_edge_modifiers",
) -> HybridDirectedGraph:
    """Orient STRING edges; never add an edge absent from the STRING graph.

    A STRING pair with one or two OmniPath directions receives only those arcs.
    A pair without direction evidence receives equal-weight arcs both ways.
    Sign is retained as metadata and never changes the non-negative weight.
    """

    if string_graph.is_directed():
        raise ValueError("Hybrid builder expects the classic undirected STRING graph")
    if "name" not in string_graph.vs.attributes():
        raise ValueError("STRING graph vertices require a name attribute")
    if "weight" not in string_graph.es.attributes():
        raise ValueError("STRING graph edges require a weight attribute")

    names = [str(value) for value in string_graph.vs["name"]]
    name_set = set(names)
    node_layer1 = {
        node.node_id: str(node.layer1_id)
        for node in signaling_dataset.nodes if node.layer1_id
    }
    evidence: dict[tuple[str, str], list[SignalingEdge]] = defaultdict(list)
    mapped_omnipath_pairs: set[frozenset[str]] = set()
    for edge in signaling_dataset.edges:
        mapped = _layer1_edge(edge, node_layer1)
        if mapped is None:
            continue
        source, target = mapped
        if source not in name_set or target not in name_set:
            continue
        evidence[(source, target)].append(edge)
        mapped_omnipath_pairs.add(frozenset((source, target)))

    vertex_index = {name: index for index, name in enumerate(names)}
    arcs: list[tuple[int, int]] = []
    attrs: dict[str, list] = defaultdict(list)
    string_pairs: set[frozenset[str]] = set()
    resolved = single = both = fallback = 0

    def add_arc(source: str, target: str, weight: float, supported: list[SignalingEdge] | None) -> None:
        arcs.append((vertex_index[source], vertex_index[target]))
        items = supported or []
        resources = sorted({resource for item in items for resource in item.resources})
        references = sorted({reference for item in items for reference in item.references})
        signs = sorted({item.sign_status.value for item in items})
        attrs["weight"].append(weight)
        attrs["distance"].append(strength_to_distance(weight))
        attrs["original_string_weight"].append(weight)
        attrs["direction_source"].append("OmniPath" if items else "undirected_fallback")
        attrs["direction_evidence"].append(";".join(resources))
        attrs["direction_resources"].append(";".join(resources))
        attrs["direction_references"].append(";".join(references))
        attrs["sign_metadata"].append(";".join(signs) if signs else "not_available")
        attrs["fallback_bidirectional"].append(not bool(items))
        attrs["tissue_support"].append(tissue_context)

    active_direction_policy = SYMMETRIC_CONTROL_POLICY if force_symmetric_fallback else DIRECTION_POLICY
    for edge in string_graph.es:
        source, target = names[edge.source], names[edge.target]
        pair = frozenset((source, target))
        string_pairs.add(pair)
        weight = float(edge["weight"])
        forward = [] if force_symmetric_fallback else evidence.get((source, target), [])
        reverse = [] if force_symmetric_fallback else evidence.get((target, source), [])
        if forward or reverse:
            resolved += 1
            if forward and reverse:
                both += 1
            else:
                single += 1
            if forward:
                add_arc(source, target, weight, forward)
            if reverse:
                add_arc(target, source, weight, reverse)
        else:
            fallback += 1
            add_arc(source, target, weight, None)
            add_arc(target, source, weight, None)

    graph = ig.Graph(n=len(names), edges=arcs, directed=True)
    for attribute in string_graph.vs.attributes():
        graph.vs[attribute] = list(string_graph.vs[attribute])
    for attribute, values in attrs.items():
        graph.es[attribute] = values

    fingerprint = signaling_dataset.fingerprint
    cache_key = directed_graph_cache_key(
        tissue_context=tissue_context,
        string_dataset_version=string_dataset_version,
        omnipath_fingerprint=fingerprint,
        direction_policy=active_direction_policy,
        taxon_id=taxon_id,
        effective_graph_fingerprint=effective_graph_fingerprint,
        edge_evidence_weighting_mode=edge_evidence_weighting_mode,
    )
    provenance = HybridGraphProvenance(
        string_edges=string_graph.ecount(), direction_resolved_edges=resolved,
        single_direction_edges=single, bidirectional_evidence_edges=both,
        fallback_string_edges=fallback, directed_arcs=graph.ecount(),
        omnipath_layer1_edges_not_in_string=len(mapped_omnipath_pairs - string_pairs),
        omnipath_mapping_misses=int(signaling_dataset.mapping.unusable_edges),
        affected_unique_nodes=len({node for pair in (mapped_omnipath_pairs & string_pairs) for node in pair}),
        dropped_string_edges=0,
        direction_coverage=(resolved / string_graph.ecount() if string_graph.ecount() else 0.0),
        tissue_context=tissue_context, string_dataset_version=string_dataset_version,
        omnipath_fingerprint=fingerprint, direction_policy=active_direction_policy,
        fallback_policy=FALLBACK_POLICY, cache_key=cache_key,
    )
    return HybridDirectedGraph(graph=graph, provenance=provenance)
