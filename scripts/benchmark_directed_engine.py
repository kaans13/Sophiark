"""Reproducible human benchmark for Sophiark's optional directed engine."""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.config import DB_PATH, HIGH_CONF_THRESHOLD  # noqa: E402
from src.directed import (  # noqa: E402
    build_hybrid_directed_graph, compare_classic_and_directed,
    compute_directed_betweenness, run_directed_calculation,
)
from src.graph_engine import build_high_conf_subgraph  # noqa: E402
from src.scientific.efficiency import strength_to_distance  # noqa: E402
from src.signaling.storage import load_signaling_cache  # noqa: E402
from src.signaling.analysis import analyze_signaling_context  # noqa: E402
from src.signaling.graph import DirectedSignalingGraph  # noqa: E402


TARGETS = {
    "EGFR": ("ENSP00000275493",),
    "TP53": ("ENSP00000269305",),
    "MYC": ("ENSP00000478887",),
    "TP53_EGFR_MYC": ("ENSP00000269305", "ENSP00000275493", "ENSP00000478887"),
}


def _rss_mb() -> float | None:
    try:
        import psutil
        return psutil.Process().memory_info().rss / 1_048_576
    except Exception:
        return None


def _classic_reference(graph, targets, fraction: float, damping: float) -> tuple[pd.DataFrame, dict]:
    """The classic engine's structural-only PR attenuation math on a graph copy."""
    baseline_values = graph.pagerank(damping=damping, weights="weight", directed=False)
    baseline = dict(zip(map(str, graph.vs["name"]), map(float, baseline_values)))
    perturbed = graph.copy()
    name_to_index = {str(v["name"]): v.index for v in perturbed.vs}
    edge_ids = set()
    for target in targets:
        if target in name_to_index:
            edge_ids.update(perturbed.incident(name_to_index[target], mode="all"))
    for edge_id in edge_ids:
        weight = float(perturbed.es[edge_id]["weight"]) * fraction
        perturbed.es[edge_id]["weight"] = weight
        perturbed.es[edge_id]["distance"] = strength_to_distance(weight)
    after_values = perturbed.pagerank(damping=damping, weights="weight", directed=False)
    after = dict(zip(map(str, perturbed.vs["name"]), map(float, after_values)))
    rows = []
    for gene in baseline:
        if gene in targets:
            continue
        delta = after[gene] - baseline[gene]
        rows.append({
            "gene": gene, "Delta_PageRank": delta,
            "Delta_PageRank_Pct": delta / baseline[gene] * 100.0 if baseline[gene] else 0.0,
        })
    return pd.DataFrame(rows), {
        target: {"baseline": baseline[target], "perturbed": after[target], "delta": after[target] - baseline[target]}
        for target in targets if target in baseline
    }


def _top(table: pd.DataFrame, response: str, symbols: dict[str, str], count: int = 10) -> list[dict]:
    ranked = table.assign(_abs=pd.to_numeric(table[response], errors="coerce").abs()).nlargest(count, "_abs")
    return [
        {"entity": row["entity"] if "entity" in row else row["gene"],
         "symbol": symbols.get(row["entity"] if "entity" in row else row["gene"], ""),
         "response_pct": float(row[response])}
        for _, row in ranked.iterrows()
    ]


def _region_coverage(graph, node_indices: set[int]) -> dict:
    edge_ids = {edge_id for index in node_indices for edge_id in graph.incident(index, mode="all")}
    supported = {edge_id for edge_id in edge_ids if not graph.es[edge_id]["fallback_bidirectional"]}
    return {
        "arcs": len(edge_ids), "omnipath_arcs": len(supported),
        "resolved_fraction": len(supported) / len(edge_ids) if edge_ids else 0.0,
    }


def _target_coverage(hybrid, targets) -> dict:
    graph = hybrid.graph
    names = {str(v["name"]): v.index for v in graph.vs}
    result = {}
    for target in targets:
        index = names.get(target)
        if index is None:
            result[target] = {"incident_arcs": 0, "omnipath_arcs": 0, "coverage": 0.0}
            continue
        incident = set(graph.incident(index, mode="all"))
        supported = {edge_id for edge_id in incident if not graph.es[edge_id]["fallback_bidirectional"]}
        one_hop = {index, *graph.neighbors(index, mode="all")}
        two_hop = set(one_hop)
        for node in tuple(one_hop):
            two_hop.update(graph.neighbors(node, mode="all"))
        result[target] = {
            "incident_arcs": len(incident), "omnipath_arcs": len(supported),
            "coverage": len(supported) / len(incident) if incident else 0.0,
            "in_degree": graph.indegree(index), "out_degree": graph.outdegree(index),
            "fallback_arc_ratio": 1.0 - (len(supported) / len(incident) if incident else 0.0),
            "one_hop": _region_coverage(graph, one_hop),
            "two_hop": _region_coverage(graph, two_hop),
        }
    return result


def run_benchmark(
    *, tissue: str, bc_sample_sources: int, output: Path,
    scenario_names: tuple[str, ...] = tuple(TARGETS), top_ks: tuple[int, ...] = (20, 50, 100, 300),
    damping: float = 0.85, suppression_factor: float = 0.001,
) -> dict:
    symbols_frame = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbols = dict(zip(symbols_frame["gene"].astype(str), symbols_frame["Symbol"].astype(str)))
    all_targets = tuple(dict.fromkeys(target for group in TARGETS.values() for target in group))
    rss_before = _rss_mb()
    build_started = time.perf_counter()
    graph = build_high_conf_subgraph(DB_PATH, HIGH_CONF_THRESHOLD, forced_genes=list(all_targets), hedef_doku=tissue)
    string_build_seconds = time.perf_counter() - build_started
    rss_after_string = _rss_mb()
    dataset = load_signaling_cache(ROOT / "data" / "processed" / "omnipath_signaling_9606.sqlite")
    rss_after_dataset = _rss_mb()
    directed_build_started = time.perf_counter()
    hybrid = build_hybrid_directed_graph(
        graph, dataset, tissue_context=tissue,
        string_dataset_version="STRING v12.0 local 9606 snapshot",
    )
    directed_build_seconds = time.perf_counter() - directed_build_started
    rss_after_build = _rss_mb()
    signaling_graph = DirectedSignalingGraph(dataset)
    bc_started = time.perf_counter()
    shared_directed_bc = compute_directed_betweenness(
        hybrid.graph, sample_sources=bc_sample_sources, random_seed=42,
    )
    directed_bc_seconds = time.perf_counter() - bc_started
    scenarios = {}
    for label in scenario_names:
        targets = TARGETS[label]
        scenario_started = time.perf_counter()
        classic_started = time.perf_counter()
        classic, classic_target_pr = _classic_reference(graph, targets, suppression_factor, damping)
        classic_seconds = time.perf_counter() - classic_started
        directed = run_directed_calculation(
            hybrid, targets=targets, block_weight_fraction=suppression_factor, damping=damping,
            bc_sample_sources=bc_sample_sources,
            precomputed_directed_bc=shared_directed_bc,
        )
        compare_started = time.perf_counter()
        comparison = compare_classic_and_directed(classic, directed, top_ks=top_ks)
        compare_seconds = time.perf_counter() - compare_started
        table = comparison.table
        amplified = table.nlargest(10, "Response_Delta")
        attenuated = table.nsmallest(10, "Response_Delta")
        layer2_input = directed.report.nlargest(min(20, len(directed.report)), "Directed_Abs_Delta_PageRank").rename(
            columns={"Directed_Redistribution_Pct": "Delta_PageRank_Pct", "Directed_Response_Direction": "Response_Direction"}
        )
        layer2 = analyze_signaling_context(signaling_graph, layer2_input, targets=targets, max_depth=4)
        relation_counts = pd.Series([row.directed_relation for row in layer2.rows]).value_counts().to_dict()
        sign_counts = pd.Series([row.sign for row in layer2.rows]).value_counts().to_dict()
        top_directed_entities = set(directed.report.nlargest(min(100, len(directed.report)), "Directed_Abs_Delta_PageRank")["entity"])
        top_nodes = {hybrid.graph.vs.find(name=entity).index for entity in top_directed_entities}
        classic_positive = int((classic["Delta_PageRank_Pct"] > 0).sum())
        classic_loss = int((classic["Delta_PageRank_Pct"] < 0).sum())
        directed_positive = int((directed.report["Directed_Redistribution_Pct"] > 0).sum())
        directed_loss = int((directed.report["Directed_Redistribution_Pct"] < 0).sum())
        boundary_k = min(100, len(table))
        classic_top = set(table.nsmallest(boundary_k, "Classic_Rank")["entity"])
        directed_top = set(table.nsmallest(boundary_k, "Directed_Rank")["entity"])
        scenarios[label] = {
            "targets": list(targets),
            "target_symbols": [symbols.get(target, target) for target in targets],
            "classic_top_response": _top(classic.rename(columns={"gene": "entity"}), "Delta_PageRank_Pct", symbols, 10),
            "directed_top_response": _top(directed.report, "Directed_Redistribution_Pct", symbols, 10),
            "validation": dict(comparison.validation),
            "effect_size": {
                "mean_absolute_response_delta_pct": float(table["Response_Delta"].abs().mean()),
                "median_absolute_response_delta_pct": float(table["Response_Delta"].abs().median()),
                "mean_absolute_rank_shift": float(table["Rank_Shift"].abs().mean()),
                "redistribution_sign_changes": int(((table["Classic_Redistribution_Pct"] * table["Directed_Redistribution_Pct"]) < 0).sum()),
                "classic_only_top_100": sorted(classic_top - directed_top),
                "directed_only_top_100": sorted(directed_top - classic_top),
            },
            "candidate_counts": {
                "classic_positive": classic_positive, "classic_loss": classic_loss,
                "directed_positive": directed_positive, "directed_loss": directed_loss,
            },
            "target_pagerank": {
                "classic": classic_target_pr,
                "directed": directed.metadata.get("target_pagerank", {}),
            },
            "largest_direction_amplified": [
                {"entity": row.entity, "symbol": symbols.get(row.entity, ""), "response_delta": float(row.Response_Delta)}
                for row in amplified.itertuples()
            ],
            "largest_direction_attenuated": [
                {"entity": row.entity, "symbol": symbols.get(row.entity, ""), "response_delta": float(row.Response_Delta)}
                for row in attenuated.itertuples()
            ],
            "direction_coverage_around_targets": _target_coverage(hybrid, targets),
            "top_response_neighborhood_direction_coverage": _region_coverage(hybrid.graph, top_nodes),
            "layer2_biological_sanity": {
                "status": layer2.status.value, "relation_counts": relation_counts,
                "sign_context_counts": sign_counts,
                "alternative_route_candidates": sum(row.alternative_route_candidate for row in layer2.rows),
                "limitations": list(layer2.notices),
            },
            "runtime_seconds": {
                "classic_reference": classic_seconds,
                **dict(directed.timings),
                "comparison": compare_seconds,
                "compare_total": time.perf_counter() - scenario_started,
            },
        }
    result = {
        "status": "COMPLETE", "tissue": tissue,
        "graph": {
            "nodes": graph.vcount(), "string_edges": graph.ecount(),
            **asdict(hybrid.provenance),
        },
        "performance": {
            "string_graph_build_seconds": string_build_seconds,
            "directed_graph_build_seconds": directed_build_seconds,
            "shared_directed_bc_seconds": directed_bc_seconds,
            "total_memory_impact_mb": (
                rss_after_build - rss_before if rss_before is not None and rss_after_build is not None else None
            ),
            "string_graph_memory_mb": (
                rss_after_string - rss_before if rss_before is not None and rss_after_string is not None else None
            ),
            "omnipath_dataset_memory_mb": (
                rss_after_dataset - rss_after_string
                if rss_after_string is not None and rss_after_dataset is not None else None
            ),
            "hybrid_directed_graph_memory_mb": (
                rss_after_build - rss_after_dataset
                if rss_after_dataset is not None and rss_after_build is not None else None
            ),
        },
        "scenarios": scenarios,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default="Lung")
    parser.add_argument("--bc-sample-sources", type=int, default=4)
    parser.add_argument("--scenario", action="append", choices=tuple(TARGETS), dest="scenarios")
    parser.add_argument("--top-k", type=int, nargs="+", default=[20, 50, 100, 300])
    parser.add_argument("--damping", type=float, default=0.85)
    parser.add_argument("--suppression-factor", type=float, default=0.001)
    parser.add_argument("--perturbation-strategy", choices=["incident_edge_attenuation"], default="incident_edge_attenuation")
    parser.add_argument("--output", type=Path, default=ROOT / "docs" / "directed_engine_real_benchmark.json")
    args = parser.parse_args()
    result = run_benchmark(
        tissue=args.tissue, bc_sample_sources=args.bc_sample_sources, output=args.output,
        scenario_names=tuple(args.scenarios or TARGETS), top_ks=tuple(args.top_k),
        damping=args.damping, suppression_factor=args.suppression_factor,
    )
    print(json.dumps({"output": str(args.output), "performance": result["performance"]}, indent=2))


if __name__ == "__main__":
    main()
