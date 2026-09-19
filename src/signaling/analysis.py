"""Post-process immutable Sophiark responses against the separate signaling graph."""

from __future__ import annotations

import hashlib
import json
import time
from collections import defaultdict
from typing import Iterable

import pandas as pd

from .graph import DirectedSignalingGraph
from .models import (
    ConvergenceRecord, SignalingAnalysisResult, SignalingContextRow, SignalingStatus,
)


LIMITATION = (
    "Directed signaling connectivity does not establish activity, compensation, drug resistance, "
    "causality in the analyzed tissue, or experimental validation. Network redistribution and "
    "biochemical activation sign are separate quantities."
)


def analysis_cache_key(graph: DirectedSignalingGraph, targets: Iterable[str], *, max_depth: int) -> str:
    payload = {
        "dataset": graph.dataset.fingerprint,
        "targets": sorted(map(str, targets)),
        "max_depth": int(max_depth),
        "sign_policy": "uncertain_on_unsigned_conflict_or_mixed-v1",
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode("utf-8")).hexdigest()


def _response_value(row: pd.Series) -> float | None:
    for field in ("Delta_PageRank_Pct", "Delta_PageRank", "PageRank_Degisim_Yuzde"):
        if field in row and pd.notna(row[field]):
            try:
                return float(row[field])
            except (TypeError, ValueError):
                pass
    return None


def _response_identifier(row: pd.Series, graph: DirectedSignalingGraph) -> str | None:
    for field in ("gene", "Symbol", "symbol", "response_entity"):
        if field in row and pd.notna(row[field]):
            resolved = graph.resolve_node(row[field])
            if resolved is not None:
                return resolved
    return None


def _tissue_support(graph: DirectedSignalingGraph, target: str, response: str, supported: set[str] | None) -> str:
    if supported is None:
        return "tissue support unavailable"
    target_supported, response_supported = target in supported, response in supported
    if target_supported and response_supported:
        return "both source and target tissue-supported"
    if target_supported:
        return "source-only"
    if response_supported:
        return "target-only"
    return "neither supported"


def _alternative_route(
    graph: DirectedSignalingGraph, target: str, response: str, *, max_depth: int,
) -> bool:
    target_region = set(graph.bounded_distances(target, max_depth=max_depth))
    response_region = set(graph.bounded_distances(response, max_depth=max_depth))
    shared = (target_region & response_region) - {target, response}
    return bool(shared)


def directed_convergence(
    graph: DirectedSignalingGraph, targets: Iterable[str], *, max_depth: int = 4, limit: int = 100,
) -> tuple[ConvergenceRecord, ...]:
    resolved_targets = tuple(dict.fromkeys(
        node for target in targets if (node := graph.resolve_node(target)) is not None
    ))
    memberships: dict[str, dict[str, int]] = defaultdict(dict)
    for target in resolved_targets:
        for node, distance in graph.bounded_distances(target, max_depth=max_depth).items():
            memberships[node][target] = distance
    records: list[ConvergenceRecord] = []
    for node, distances in memberships.items():
        if len(distances) < 2:
            continue
        signs: dict[str, str] = {}
        resources: set[str] = set()
        references: set[str] = set()
        for target in sorted(distances):
            summary = graph.summarize_paths(target, node, max_depth=max_depth)
            signs[target] = summary.net_sign
            found_resources, found_references = graph.path_evidence(summary.example_path)
            resources.update(found_resources)
            references.update(found_references)
        records.append(ConvergenceRecord(
            convergence_node=node, converging_targets=tuple(sorted(distances)),
            distances=dict(sorted(distances.items())), signed_relations=signs,
            resources=tuple(sorted(resources)), references=tuple(sorted(references)),
        ))
    return tuple(sorted(records, key=lambda item: (sum(item.distances.values()), item.convergence_node))[:limit])


def analyze_signaling_context(
    graph: DirectedSignalingGraph,
    report: pd.DataFrame,
    *,
    targets: Iterable[str],
    max_depth: int = 4,
    tissue_supported_entities: Iterable[str] | None = None,
) -> SignalingAnalysisResult:
    """Annotate only existing response rows; never rerank or alter the report."""

    started = time.perf_counter()
    if not isinstance(report, pd.DataFrame):
        raise TypeError("report must be a pandas DataFrame")
    requested_targets = tuple(map(str, targets))
    resolved_targets = tuple(dict.fromkeys(
        node for target in requested_targets if (node := graph.resolve_node(target)) is not None
    ))
    if not resolved_targets:
        return SignalingAnalysisResult(
            status=SignalingStatus.PARTIAL, error=None,
            metadata={"dataset_fingerprint": graph.dataset.fingerprint, "max_depth": max_depth},
            notices=("No perturbation target resolved in the signaling layer.", LIMITATION),
        )
    supported = None
    if tissue_supported_entities is not None:
        supported = {node for value in tissue_supported_entities if (node := graph.resolve_node(value)) is not None}

    rows: list[SignalingContextRow] = []
    for _, source_row in report.iterrows():
        response = _response_identifier(source_row, graph)
        if response is None:
            continue
        response_type = str(source_row.get("Response_Direction") or source_row.get("Hasar_Tipi") or "network response")
        for target in resolved_targets:
            relation, summary = graph.relation(target, response, max_depth=max_depth)
            resources, references = graph.path_evidence(summary.example_path)
            rows.append(SignalingContextRow(
                perturbation_target=target, response_entity=response,
                response_type=response_type, network_response_value=_response_value(source_row),
                directed_relation=relation, directed_distance=summary.shortest_distance,
                sign=summary.net_sign, example_path=summary.example_path,
                path_count_within_depth=summary.path_count,
                positive_path_count=summary.positive_paths,
                negative_path_count=summary.negative_paths,
                unsigned_path_count=summary.unsigned_paths,
                uncertain_path_count=summary.uncertain_paths,
                alternative_route_candidate=_alternative_route(graph, target, response, max_depth=max_depth),
                tissue_support=_tissue_support(graph, target, response, supported),
                evidence_resources=resources, references=references,
                dataset_version=str(graph.dataset.metadata.get("dataset_version", "local_snapshot")),
            ))
    convergence = directed_convergence(graph, resolved_targets, max_depth=max_depth)
    elapsed_ms = (time.perf_counter() - started) * 1000.0
    status = SignalingStatus.PARTIAL if graph.dataset.mapping.acceptance == "FAILURE" else SignalingStatus.AVAILABLE
    notices = [LIMITATION]
    if graph.dataset.mapping.acceptance == "FAILURE":
        notices.append(
            "Layer-1 edge mapping is below the production acceptance threshold; Layer-2-only "
            "entities are preserved and never mixed into Sophiark scores."
        )
    return SignalingAnalysisResult(
        status=status, rows=tuple(rows), convergence=convergence,
        metadata={
            "dataset_fingerprint": graph.dataset.fingerprint,
            "cache_key": analysis_cache_key(graph, requested_targets, max_depth=max_depth),
            "max_depth": max_depth, "target_count": len(resolved_targets),
            "response_rows": len(report), "annotated_rows": len(rows),
            "annotation_time_ms": round(elapsed_ms, 3), "mapping_acceptance": graph.dataset.mapping.acceptance,
            "dataset_version": graph.dataset.metadata.get("dataset_version", "local_snapshot"),
            "dataset_date": graph.dataset.metadata.get("downloaded_at"),
            "normalized_edge_count": graph.dataset.qa.normalized_interactions,
            "mapped_raw_edge_count": graph.dataset.mapping.fully_layer1_mapped_edges,
            "directed_edge_count": graph.dataset.qa.unique_directed_edges,
            "signed_edge_count": graph.dataset.qa.activating_edges + graph.dataset.qa.inhibitory_edges,
            "unsigned_edge_count": graph.dataset.qa.unsigned_edges,
            "conflicting_edge_count": graph.dataset.qa.conflicting_edges,
            "atomic_entity_mapping_rate": graph.dataset.mapping.atomic_entity_mapping_rate,
            "atomic_edge_mapping_rate": graph.dataset.mapping.atomic_edge_mapping_rate,
            "raw_layer1_edge_mapping_rate": graph.dataset.mapping.edge_mapping_rate,
        },
        notices=tuple(notices),
    )


def safe_analyze_signaling_context(
    graph: DirectedSignalingGraph | None,
    report: pd.DataFrame,
    **kwargs,
) -> SignalingAnalysisResult:
    """Failure boundary: signaling failure never invalidates core output."""

    if graph is None:
        return SignalingAnalysisResult(
            status=SignalingStatus.UNAVAILABLE,
            notices=("OmniPath signaling layer is unavailable; core Sophiark results are unchanged.",),
        )
    try:
        return analyze_signaling_context(graph, report, **kwargs)
    except Exception as exc:
        return SignalingAnalysisResult(
            status=SignalingStatus.UNAVAILABLE,
            notices=("OmniPath signaling analysis failed; core Sophiark results are unchanged.",),
            error=f"{type(exc).__name__}: {exc}",
        )
