"""Evidence-aware perturbation propagation trace over the existing PPI graph.

Structural PPI paths, directed regulatory annotations and measured simulation
responses are intentionally preserved as separate columns.  No causal or
biochemical signal-flow inference is made here.
"""
from __future__ import annotations

import math
from pathlib import Path

import pandas as pd

from src.models import PropagationTrace
from src.scientific.regulatory import OVERLAY_COLUMNS, load_omnipath_overlay, trrust_records


def _number(value):
    try:
        result = float(value)
        return result if math.isfinite(result) else None
    except (TypeError, ValueError):
        return None


def _truth(value) -> bool:
    return str(value).casefold() in {"true", "1", "yes"} or value is True


def _community_context(value) -> str:
    """Presentation column: preserve IDs while keeping Arrow serialization homogeneous."""
    if value is None or pd.isna(value):
        return "N/A"
    return str(value)


def _response_candidates(report: pd.DataFrame, target_gene: str, max_nodes: int) -> set[str]:
    """KatmanlÄ± preselection: significance > strong |dPR| > critical bridge > top |dPR|.

    Burada weighted composite score yoktur; her kaynak ayrÄ±ca izlenebilir.
    """
    frame = report.drop_duplicates("gene").copy()
    frame["gene"] = frame["gene"].astype(str)
    selected = {target_gene}
    if "significant_redistribution" in frame:
        selected |= set(frame.loc[frame["significant_redistribution"].map(_truth), "gene"])
    if "Strong_Redistribution" in frame:
        selected |= set(frame.loc[frame["Strong_Redistribution"].map(_truth), "gene"])
    if "GÃ¼mrÃ¼k_Kapisi" in frame:
        selected |= set(frame.loc[frame["GÃ¼mrÃ¼k_Kapisi"].map(_truth), "gene"])
    # Top absolute ΔPR is a transparent fallback for finite output size.
    if "Abs_Delta_PageRank" in frame:
        selected |= set(frame.nlargest(max_nodes, "Abs_Delta_PageRank")["gene"])
    if len(selected) <= max_nodes:
        return selected
    # Target + statistically selected nodes are protected; overflow uses abs ΔPR.
    protected = {target_gene}
    for flag in ("significant_redistribution", "Strong_Redistribution", "GÃ¼mrÃ¼k_Kapisi"):
        if flag in frame:
            protected |= set(frame.loc[frame[flag].map(_truth), "gene"])
    room = max(0, max_nodes - len(protected))
    top = frame[~frame["gene"].isin(protected)].nlargest(room, "Abs_Delta_PageRank" if "Abs_Delta_PageRank" in frame else "Hinterland_Skoru")
    return protected | set(top["gene"])


def _edge_evidence(source_symbol: str, target_symbol: str, overlay: pd.DataFrame) -> pd.DataFrame:
    if overlay.empty:
        return overlay
    return overlay[((overlay["source"] == source_symbol) & (overlay["target"] == target_symbol)) |
                   ((overlay["source"] == target_symbol) & (overlay["target"] == source_symbol))]


def _route_type(evidence: pd.DataFrame) -> str:
    if evidence.empty:
        return "Structural Network Route"
    status_column = "evidence_status" if "evidence_status" in evidence.columns else "Evidence_Status"
    statuses = set(evidence[status_column].dropna().astype(str))
    if "conflicting_evidence" in statuses:
        return "Mixed Evidence Route (conflicting directed evidence)"
    return "Mixed Evidence Route"


def filter_trace_nodes(nodes: pd.DataFrame, *, layer_limit: int = 3, critical_only: bool = False,
                       significant_only: bool = False, directed_only: bool = False) -> pd.DataFrame:
    if nodes is None or nodes.empty:
        return pd.DataFrame()
    result = nodes[nodes["Propagation_Layer"] <= layer_limit].copy()
    if critical_only:
        result = result[result["Critical_Context"].astype(bool)]
    if significant_only:
        result = result[result["Significant_Redistribution"].astype(bool)]
    if directed_only:
        result = result[result["Directed_Evidence_Involvement"].astype(bool)]
    return result


def build_propagation_trace(*, graph, report: pd.DataFrame, target_gene: str,
                            baseline_scores: pd.DataFrame, gene_to_symbol: dict[str, str],
                            regulators: dict | None = None, omnipath_path: str | Path | None = None,
                            max_nodes: int = 40, max_routes: int = 12, max_hops: int = 3,
                            trace_context: dict | None = None) -> PropagationTrace:
    """Build bounded, deterministic trace from an already measured forward report."""
    empty = PropagationTrace(target_gene, pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), pd.DataFrame(), {})
    if report is None or report.empty or "gene" not in report.columns:
        return empty
    name_to_idx = {str(item): i for i, item in enumerate(graph.vs["name"])}
    if target_gene not in name_to_idx:
        return empty
    overlay, overlay_meta = load_omnipath_overlay(omnipath_path)
    trrust = trrust_records(regulators or {})
    directed = pd.concat([trrust, overlay], ignore_index=True) if not overlay.empty else trrust
    candidates = _response_candidates(report, target_gene, max_nodes)
    # Always include direct neighbourhood even if it had no visible report row.
    target_idx = name_to_idx[target_gene]
    neighbourhood = {str(graph.vs[index]["name"]) for index in graph.neighbors(target_idx)}
    candidates |= neighbourhood
    candidate_indices = [name_to_idx[item] for item in candidates if item in name_to_idx]
    hops = graph.distances(source=target_idx, target=candidate_indices, weights=None)[0]
    rows, selected_nodes = [], set()
    report_indexed = report.drop_duplicates("gene").set_index("gene", drop=False)
    for index, hop in zip(candidate_indices, hops):
        gene = str(graph.vs[index]["name"])
        if not math.isfinite(hop) or hop > max_hops:
            continue
        source = report_indexed.loc[gene] if gene in report_indexed.index else pd.Series(dtype=object)
        symbol = gene_to_symbol.get(gene, str(source.get("Symbol", gene)))
        # Directed involvement is evaluated only after route selection; target-neighbour evidence is nevertheless visible.
        edge_support = _edge_evidence(gene_to_symbol.get(target_gene, target_gene), symbol, directed)
        critical = _truth(source.get("GÃ¼mrÃ¼k_Kapisi")) or _truth(source.get("Strong_Redistribution"))
        selected_nodes.add(gene)
        rows.append({
            "Gene_Symbol": symbol, "Network_ID": gene, "Propagation_Layer": int(hop),
            "Hop_Distance": int(hop), "Delta_PageRank": _number(source.get("Delta_PageRank")),
            "Response_Direction": source.get("Response_Direction", "N/A"),
            "empirical_p": _number(source.get("empirical_p")), "q_value": _number(source.get("q_value")),
            "Significant_Redistribution": _truth(source.get("significant_redistribution")),
            "Delta_BC": None, "BC_or_Congestion": _number(source.get("BC_Skoru")),
            "WBI_or_Bridge_Change": "N/A", "Community_Context": _community_context(source.get("Topluluk_ID", "N/A")),
            "Tissue_Relevance": source.get("Tissue", (trace_context or {}).get("tissue", "N/A")),
            "Critical_Context": critical, "Directed_Evidence_Involvement": not edge_support.empty,
            "Explanation": "Perturbed target" if hop == 0 else "Measured response node in bounded structural neighbourhood",
        })
    nodes = pd.DataFrame(rows)
    if nodes.empty:
        return empty
    nodes = nodes.sort_values(["Propagation_Layer", "Significant_Redistribution", "Critical_Context", "Delta_PageRank"], ascending=[True, False, False, False], na_position="last").head(max_nodes).reset_index(drop=True)
    edges, routes, route_evidence = [], [], []
    route_targets = nodes[(nodes["Propagation_Layer"] > 0) & (nodes["Network_ID"] != target_gene)].head(max_routes)
    for _, endpoint in route_targets.iterrows():
        path = graph.get_shortest_paths(target_idx, to=name_to_idx[endpoint["Network_ID"]], weights="distance", output="vpath")[0]
        if not path:
            continue
        labels = [gene_to_symbol.get(str(graph.vs[idx]["name"]), str(graph.vs[idx]["name"])) for idx in path]
        evidence_parts = []
        for left, right in zip(path, path[1:]):
            left_id, right_id = str(graph.vs[left]["name"]), str(graph.vs[right]["name"])
            left_symbol, right_symbol = gene_to_symbol.get(left_id, left_id), gene_to_symbol.get(right_id, right_id)
            support = _edge_evidence(left_symbol, right_symbol, directed)
            statuses = ";".join(sorted(set(support.get("evidence_status", pd.Series(dtype=str)).dropna().astype(str)))) if not support.empty else "structural_only"
            edges.append({"Source": left_symbol, "Target": right_symbol, "Structural_Edge": True,
                          "Directed_Evidence": not support.empty, "Evidence_Status": statuses,
                          "Directed_Arrow": bool(not support.empty and ((support["source"] == left_symbol) & (support["target"] == right_symbol)).any())})
            if not support.empty:
                for _, item in support.iterrows():
                    route_evidence.append({"Source": left_symbol, "Target": right_symbol, "Evidence_Layer": item.get("evidence_layer"),
                                           "Evidence_Status": item.get("evidence_status"), "Stimulation": item.get("stimulation"),
                                           "Inhibition": item.get("inhibition"), "Provenance": item.get("provenance"), "References": item.get("references")})
                evidence_parts.append("directed support")
        route_type = _route_type(pd.DataFrame(route_evidence)) if evidence_parts else "Structural Network Route"
        routes.append({"Affected_Gene": endpoint["Gene_Symbol"], "Route": " -> ".join(labels), "Route_Type": route_type,
                       "Path_Length": len(path) - 1, "Affected_Node_Significant": bool(endpoint["Significant_Redistribution"]),
                       "Directed_Edge_Count": len(evidence_parts), "Interpretation": "Predicted Network Propagation Route; not a causal mechanism."})
    edges_df = pd.DataFrame(edges).drop_duplicates() if edges else pd.DataFrame(columns=["Source", "Target", "Structural_Edge", "Directed_Evidence", "Evidence_Status", "Directed_Arrow"])
    evidence_df = pd.DataFrame(route_evidence).drop_duplicates() if route_evidence else pd.DataFrame(columns=["Source", "Target", "Evidence_Layer", "Evidence_Status", "Stimulation", "Inhibition", "Provenance", "References"])
    metadata = {"max_nodes": max_nodes, "max_routes": max_routes, "max_hops": max_hops, "approximation": "bounded candidate subgraph and one weighted shortest route per selected endpoint", "omnipath": overlay_meta, **(trace_context or {})}
    return PropagationTrace(target_gene, nodes, pd.DataFrame(routes), edges_df, evidence_df, metadata)


def export_propagation_trace(trace: PropagationTrace, output_dir: str | Path) -> Path:
    root = Path(output_dir)
    root.mkdir(parents=True, exist_ok=True)
    trace.nodes.to_csv(root / "nodes.csv", index=False)
    trace.edges.to_csv(root / "edges.csv", index=False)
    trace.routes.to_csv(root / "routes.csv", index=False)
    trace.evidence.to_csv(root / "evidence.csv", index=False)
    import json
    (root / "metadata.json").write_text(json.dumps(trace.metadata, ensure_ascii=False, indent=2, default=str), encoding="utf-8")
    return root
