"""Frozen-E1 STRING channel ablation audit for Evidence BETA.

This is a diagnostic script, not an application engine.  It never changes the
Evidence default and it preserves the exact tissue-filtered Classic topology.
"""
from __future__ import annotations

import argparse
import gc
import json
import logging
import pickle
import sys
import tempfile
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic  # noqa: E402
from src.evidence.ablation import AblationCondition, NON_TEXT_CHANNELS, ablation_conditions, apply_e1_ablation  # noqa: E402
from src.evidence.config import EvidenceConfig  # noqa: E402
from src.evidence.engine import _apply_tf_only, _restore, _score_graph, prepare_evidence_network  # noqa: E402
from src.evidence.provenance import build_run_provenance  # noqa: E402
from src.evidence.scores import index_metadata  # noqa: E402
from audit_evidence_emperor_bias import (  # noqa: E402
    MEANINGFUL_RESPONSE_PCT, PANEL, _annotation_correlations, _decile, _domain_frame,
    _matched_controls, _rho, _ranks, _strata_summary, _tail_rows, _top_fraction,
)


TOP_COUNTS = (20, 50, 100)
ABSOLUTE_TOP_COUNTS = (50, 100, 300)


def _edge_records(graph, profiles: pd.DataFrame, confidence: pd.Series) -> pd.DataFrame:
    if len(profiles) != graph.ecount() or len(confidence) != graph.ecount():
        raise ValueError("E1 edge/profile alignment failed")
    rows = pd.DataFrame({
        "protein1": [graph.vs[edge.source]["name"] for edge in graph.es],
        "protein2": [graph.vs[edge.target]["name"] for edge in graph.es],
        "ablated_confidence": confidence.to_numpy(dtype=float),
        "effective_weight": np.asarray(graph.es["weight"], dtype=float),
    })
    selected = ["combined_score", "s_nontext", "text_dependency", "physical_status", *NON_TEXT_CHANNELS]
    for column in selected:
        rows[column] = profiles[column].to_numpy()
    return rows


def _variant_properties(graph, profiles: pd.DataFrame, confidence: pd.Series, base: pd.DataFrame) -> pd.DataFrame:
    edges = _edge_records(graph, profiles, confidence)
    official = pd.to_numeric(edges["combined_score"], errors="coerce") / 1000.0
    edges["ablated_to_official"] = edges["ablated_confidence"] / official.replace(0, np.nan)
    first = edges.rename(columns={"protein1": "gene"})
    second = edges.rename(columns={"protein2": "gene"})
    incident = pd.concat((first, second), ignore_index=True)
    aggregations: dict[str, tuple[str, object]] = {
        "Mean_Ablated_Confidence": ("ablated_confidence", "mean"),
        "Median_Ablated_Confidence": ("ablated_confidence", "median"),
        "Mean_Ablated_to_Official": ("ablated_to_official", "mean"),
        "Median_Ablated_to_Official": ("ablated_to_official", "median"),
        "Median_Effective_Weight": ("effective_weight", "median"),
    }
    for channel in NON_TEXT_CHANNELS:
        aggregations[f"{channel}_Availability"] = (channel, lambda s: float(pd.to_numeric(s, errors="coerce").fillna(0).gt(0).mean()))
        aggregations[f"Median_{channel}"] = (channel, lambda s: float(pd.to_numeric(s, errors="coerce").median() / 1000.0))
    values = incident.groupby("gene", as_index=False).agg(**aggregations)
    return base.merge(values, on="gene", how="left", validate="one_to_one")


def _availability(graph, profiles: pd.DataFrame, base: pd.DataFrame, *, tissue: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    confidence = pd.to_numeric(profiles["s_nontext"], errors="coerce").fillna(0.0)
    edges = _edge_records(graph, profiles, confidence)
    first, second = edges.rename(columns={"protein1": "gene"}), edges.rename(columns={"protein2": "gene"})
    incident = pd.concat((first, second), ignore_index=True)
    node = base[["gene", "Classic_Degree", "Classic_PageRank"]].copy()
    node["Degree_Decile"] = _decile(node["Classic_Degree"])
    node["PageRank_Decile"] = _decile(node["Classic_PageRank"])
    degree_rows: list[dict[str, object]] = []
    pagerank_rows: list[dict[str, object]] = []
    for channel in NON_TEXT_CHANNELS:
        values = incident.groupby("gene", as_index=False).agg(
            Nonzero_Fraction=(channel, lambda s: float(pd.to_numeric(s, errors="coerce").fillna(0).gt(0).mean())),
            Median_Channel_Score=(channel, lambda s: float(pd.to_numeric(s, errors="coerce").median() / 1000.0)),
        )
        merged = node.merge(values, on="gene", how="left", validate="one_to_one")
        for column, destination in (("Degree_Decile", degree_rows), ("PageRank_Decile", pagerank_rows)):
            for stratum, group in merged.groupby(column, observed=True):
                destination.append({"Tissue": tissue, "Channel": channel, "Stratum": str(stratum), "Node_Count": len(group),
                                    "Mean_Node_Nonzero_Fraction": float(group["Nonzero_Fraction"].mean()),
                                    "Median_Node_Nonzero_Fraction": float(group["Nonzero_Fraction"].median()),
                                    "Median_Node_Channel_Score": float(group["Median_Channel_Score"].median())})
    return pd.DataFrame(degree_rows), pd.DataFrame(pagerank_rows)


def _health(edges: pd.DataFrame, *, condition: AblationCondition, tissue: str) -> dict[str, object]:
    confidence = edges["ablated_confidence"].to_numpy(dtype=float)
    weights = edges["effective_weight"].to_numpy(dtype=float)
    return {
        "Condition": condition.identifier, "Condition_Description": condition.description, "Tissue": tissue,
        "Nodes": None, "Edges": len(edges), "Score_Median": float(np.median(confidence)),
        "Score_P5": float(np.quantile(confidence, .05)), "Score_P95": float(np.quantile(confidence, .95)),
        "Score_Near_Zero_Fraction_lt_0_05": float((confidence < .05).mean()),
        "Effective_Weight_Median": float(np.median(weights)), "Effective_Weight_P5": float(np.quantile(weights, .05)),
        "Effective_Weight_P95": float(np.quantile(weights, .95)),
    }


def _simulate_evidence(graph, scores: pd.DataFrame, *, target: str, tissue: str) -> pd.DataFrame:
    """Exact current Evidence perturbation path, without exporting a second engine result."""
    snapshot = _apply_tf_only([target], graph)
    biology_logic.apply_functional_penalty([target], graph, snapshot=snapshot)
    try:
        with tempfile.TemporaryDirectory(prefix="sophiark_channel_ablation_") as directory:
            biology_logic.run_infection_simulation(
                graph, scores, Path(directory), spesifik_hedefler=[target], block_weight_fraction=.001,
                damping=.85, exploratory_top_n=100, compute_structural_metrics=False, null_iterations=0,
                edge_evidence_weighting_mode="structural_only", tissue=tissue,
            )
            result = biology_logic.latest_signed_redistribution()
            return result.copy(deep=True) if isinstance(result, pd.DataFrame) else pd.DataFrame()
    finally:
        _restore(graph, snapshot)


def _simulate_classic(graph, scores: pd.DataFrame, *, target: str, tissue: str) -> pd.DataFrame:
    old_weight = list(graph.es["weight"])
    old_distance = list(graph.es["distance"])
    graph.es["weight"] = graph.es["classic_weight"]
    graph.es["distance"] = (1.0 / np.asarray(graph.es["weight"], dtype=float)).tolist()
    try:
        with tempfile.TemporaryDirectory(prefix="sophiark_channel_classic_") as directory:
            biology_logic.run_infection_simulation(
                graph, scores, Path(directory), spesifik_hedefler=[target], block_weight_fraction=.001,
                damping=.85, exploratory_top_n=100, compute_structural_metrics=False, null_iterations=0,
                edge_evidence_weighting_mode="legacy_edge_modifiers", tissue=tissue,
            )
            result = biology_logic.latest_signed_redistribution()
            return result.copy(deep=True) if isinstance(result, pd.DataFrame) else pd.DataFrame()
    finally:
        graph.es["weight"] = old_weight
        graph.es["distance"] = old_distance


def _benchmark(reference: pd.DataFrame, candidate: pd.DataFrame, *, label: str) -> dict[str, object]:
    left = reference[["gene", "Delta_PageRank_Pct"]].rename(columns={"Delta_PageRank_Pct": "Reference_Response"})
    right = candidate[["gene", "Delta_PageRank_Pct"]].rename(columns={"Delta_PageRank_Pct": "Ablation_Response"})
    values = left.merge(right, on="gene", how="inner", validate="one_to_one")
    ref, ablated = values["Reference_Response"], values["Ablation_Response"]
    meaningful = values[["Reference_Response", "Ablation_Response"]].abs().max(axis=1) >= MEANINGFUL_RESPONSE_PCT
    row: dict[str, object] = {
        "Reference": label, "Nodes_Compared": len(values), "Meaningful_Nodes": int(meaningful.sum()),
        "Signed_Pearson": float(ref.corr(ablated, method="pearson")), "Signed_Spearman": float(ref.corr(ablated, method="spearman")),
        "Sign_Flips_Meaningful": int(((np.sign(ref) != np.sign(ablated)) & meaningful).sum()),
    }
    for sign, cutoffs in (("positive", TOP_COUNTS), ("negative", TOP_COUNTS), ("absolute", ABSOLUTE_TOP_COUNTS)):
        ref_mask = ref > 0 if sign == "positive" else (ref < 0 if sign == "negative" else ref.notna())
        candidate_mask = ablated > 0 if sign == "positive" else (ablated < 0 if sign == "negative" else ablated.notna())
        ref_order = ref if sign == "positive" else (-ref if sign == "negative" else ref.abs())
        candidate_order = ablated if sign == "positive" else (-ablated if sign == "negative" else ablated.abs())
        for count in cutoffs:
            reference_top = set(values.loc[ref_mask].nlargest(count, ref_order.name if ref_order.name in values else "Reference_Response")["gene"])
            # nlargest cannot receive a derived Series name reliably; use an explicit helper column.
            candidate_top = set(values.assign(_order=candidate_order).loc[candidate_mask].nlargest(count, "_order")["gene"])
            if sign == "positive":
                reference_top = set(values.assign(_order=ref_order).loc[ref_mask].nlargest(count, "_order")["gene"])
            elif sign == "negative":
                reference_top = set(values.assign(_order=ref_order).loc[ref_mask].nlargest(count, "_order")["gene"])
            else:
                reference_top = set(values.assign(_order=ref_order).loc[ref_mask].nlargest(count, "_order")["gene"])
            row[f"{sign}_top{count}_overlap"] = len(reference_top & candidate_top)
    return row


def _append_bias_rows(*, ranks: pd.DataFrame, tissue: str, target: str, condition: AblationCondition,
                      degree_rows: list[pd.DataFrame], pagerank_rows: list[pd.DataFrame],
                      annotation_rows: list[pd.DataFrame], matched_rows: list[pd.DataFrame]) -> None:
    for domain in ("absolute", "positive", "negative"):
        for scope in ("all_nodes", "meaningful"):
            current = _domain_frame(ranks, domain, scope)
            if current.empty:
                continue
            degree_rows.extend((_strata_summary(current, "Degree_Decile", tissue=tissue, target=target, rank_domain=domain, scope=scope, audit="Degree"),
                                _tail_rows(current, source="Degree_Decile", tissue=tissue, target=target, rank_domain=domain, scope=scope)))
            pagerank_rows.extend((_strata_summary(current, "PageRank_Decile", tissue=tissue, target=target, rank_domain=domain, scope=scope, audit="PageRank"),
                                  _tail_rows(current, source="PageRank_Decile", tissue=tissue, target=target, rank_domain=domain, scope=scope)))
            for collection in (degree_rows[-2], degree_rows[-1], pagerank_rows[-2], pagerank_rows[-1]):
                collection.insert(0, "Condition", condition.identifier)
                collection.insert(1, "Condition_Description", condition.description)
            if domain == "absolute" and scope == "all_nodes":
                annotation = pd.DataFrame(_annotation_correlations(current, tissue=tissue, target=target))
                annotation.insert(0, "Condition", condition.identifier); annotation.insert(1, "Condition_Description", condition.description)
                annotation_rows.append(annotation)
                matched = _matched_controls(current, tissue=tissue, target=target)
                matched.insert(0, "Condition", condition.identifier); matched.insert(1, "Condition_Description", condition.description)
                matched_rows.append(matched)


def _condition_summary(bias: pd.DataFrame, matched: pd.DataFrame, annotation: pd.DataFrame) -> list[dict[str, object]]:
    rows = []
    for condition in sorted(bias["Condition"].dropna().unique()):
        selected = bias.loc[(bias["Condition"] == condition) & (bias["Rank_Domain"] == "absolute") & (bias["Response_Scope"] == "meaningful")]
        degree = selected.loc[(selected["Audit"] == "Degree") & (selected["Stratum"] == "D10")]
        pagerank = selected.loc[(selected["Audit"] == "PageRank") & (selected["Stratum"] == "D10")]
        annotation_values = annotation.loc[(annotation["Condition"] == condition) & (annotation["Metric"] == "Annotation_Coverage") & (annotation["Response_Scope"] == "meaningful_response")]
        controls = matched.loc[(matched["Condition"] == condition) & (matched["Audit"] == "Coverage_Matched_Control") & (matched["Response_Scope"] == "meaningful_response")]
        per_target = controls.groupby(["Tissue", "Target"], observed=True)["Promotion_Difference_High_Minus_Low"].median()
        rows.append({
            "Condition": condition, "Degree_D10_Median_Promotion_Mean": float(degree["Median_Promotion"].mean()),
            "Degree_D10_Promoted_Targets": int((degree["Median_Promotion"] > 0).sum()),
            "PageRank_D10_Median_Promotion_Mean": float(pagerank["Median_Promotion"].mean()),
            "PageRank_D10_Promoted_Targets": int((pagerank["Median_Promotion"] > 0).sum()),
            "Annotation_Coverage_Rho_Median": float(annotation_values["Spearman_Rho"].median()),
            "Annotation_Coverage_Positive_Targets": int((annotation_values["Spearman_Rho"] > 0).sum()),
            "Matched_Control_Difference_Median": float(per_target.median()) if len(per_target) else None,
            "Matched_Control_Positive_Targets": int((per_target > 0).sum()),
        })
    output = pd.DataFrame(rows)
    current = output.loc[output["Condition"] == "C1"].iloc[0]
    output["Delta_D10_vs_C1"] = current["Degree_D10_Median_Promotion_Mean"] - output["Degree_D10_Median_Promotion_Mean"]
    output["Delta_PR_D10_vs_C1"] = current["PageRank_D10_Median_Promotion_Mean"] - output["PageRank_D10_Median_Promotion_Mean"]
    output["Delta_Matched_vs_C1"] = current["Matched_Control_Difference_Median"] - output["Matched_Control_Difference_Median"]
    return output.to_dict("records")


def run(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.disable(logging.CRITICAL)
    config = EvidenceConfig()
    conditions = ablation_conditions()
    key_counts: dict[tuple[str, ...], int] = defaultdict(int)
    for condition in conditions:
        key_counts[tuple(condition.included_channels)] += 1
    # Retain only C1 (needed as the current-Evidence comparator) and explicit
    # alias profiles (A3/G1, A6/G2, A9/G3).  Every other diagnostic result is
    # summarized before its large graph/property objects are released.
    retained_cache_keys = {tuple(NON_TEXT_CHANNELS), *[key for key, count in key_counts.items() if count > 1]}
    mapping = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbol_map = dict(zip(mapping["gene"].astype(str), mapping["Symbol"].astype(str)))
    symbol_to_id = {symbol.upper(): gene for gene, symbol in symbol_map.items()}
    with (ROOT / "data" / "processed" / "complex_members.pkl").open("rb") as stream:
        complex_members = pickle.load(stream)
    all_bias: list[pd.DataFrame] = []; all_annotation: list[pd.DataFrame] = []; all_matched: list[pd.DataFrame] = []
    benchmarks: list[dict[str, object]] = []; health_rows: list[dict[str, object]] = []; target_local: list[dict[str, object]] = []
    transitions: list[pd.DataFrame] = []; degree_availability: list[pd.DataFrame] = []; pagerank_availability: list[pd.DataFrame] = []
    raw_rank_rows: list[pd.DataFrame] = []
    control_reproduced = True

    for tissue, symbols in PANEL.items():
        print(f"[ablation] preparing frozen E1 {tissue}", flush=True)
        targets = [symbol_to_id[symbol] for symbol in symbols]
        control_graph, control_scores, profiles = prepare_evidence_network(
            config=config, forced_genes=targets, tissue=tissue, bc_sample_sources=4,
        )
        expected_nodes, expected_edges = control_graph.vcount(), control_graph.ecount()
        if any(abs(float(value) - 1.0) > 1e-12 for value in control_graph.es["text_bonus"]) or any(abs(float(value) - 1.0) > 1e-12 for value in control_graph.es["physical_bonus"]):
            raise RuntimeError("Frozen control unexpectedly enabled text or physical weighting")
        from audit_evidence_emperor_bias import _node_properties  # preserve exact existing annotation definitions
        base_properties = _node_properties(control_graph, control_scores, symbol_map, complex_members)
        by_degree, by_pagerank = _availability(control_graph, profiles, base_properties, tissue=tissue)
        degree_availability.append(by_degree); pagerank_availability.append(by_pagerank)
        classic_responses = {symbol: _simulate_classic(control_graph, control_scores, target=target, tissue=tissue) for symbol, target in zip(symbols, targets)}
        cache: dict[tuple[str, ...], dict[str, object]] = {}

        for condition in conditions:
            key = tuple(condition.included_channels)
            if key not in cache:
                graph = control_graph.copy()
                confidence = apply_e1_ablation(graph, profiles, condition, prior=config.string_prior)
                if graph.vcount() != expected_nodes or graph.ecount() != expected_edges:
                    raise RuntimeError(f"E1 topology changed in {condition.identifier}/{tissue}")
                if any(abs(float(value) - 1.0) > 1e-12 for value in graph.es["text_bonus"]) or any(abs(float(value) - 1.0) > 1e-12 for value in graph.es["physical_bonus"]):
                    raise RuntimeError("Optional modifier enabled during ablation")
                scores = _score_graph(graph, 4)
                properties = _variant_properties(graph, profiles, confidence, base_properties)
                edges = _edge_records(graph, profiles, confidence)
                health = _health(edges, condition=condition, tissue=tissue)
                health["Nodes"] = graph.vcount()
                signed = {}
                for target_symbol, target in zip(symbols, targets):
                    print(f"[ablation] {condition.identifier} {tissue}/{target_symbol}", flush=True)
                    signed[target_symbol] = _simulate_evidence(graph, scores, target=target, tissue=tissue)
                cache[key] = {"responses": signed, "properties": properties, "health": health}
                del graph, scores, edges
                gc.collect()
            result = cache[key]
            health = dict(result["health"]); health["Condition"] = condition.identifier; health["Condition_Description"] = condition.description
            health["Included_Channels"] = ";".join(condition.included_channels); health["Excluded_Channels"] = ";".join(condition.excluded_channels)
            health_rows.append(health)
            for target_symbol, target in zip(symbols, targets):
                classic = classic_responses[target_symbol]
                ablated = result["responses"][target_symbol]
                non_target = result["properties"].loc[result["properties"]["gene"] != target].copy()
                ranks = _ranks(classic, ablated).merge(non_target, on="gene", how="left", validate="one_to_one")
                if len(ranks) != len(non_target):
                    raise RuntimeError(f"Response universe mismatch: {condition.identifier}/{tissue}/{target_symbol}")
                raw = ranks[[column for column in ("gene", "Classic_Response", "Evidence_Response", "Classic_Absolute_Rank", "Evidence_Absolute_Rank", "Classic_Degree", "Classic_PageRank", "Degree_Decile", "PageRank_Decile", "Promotion", "Meaningful_Response", "Median_Ablated_Confidence", "Median_Effective_Weight") if column in ranks]].copy()
                raw.insert(0, "Condition", condition.identifier); raw.insert(1, "Tissue", tissue); raw.insert(2, "Target", target_symbol)
                raw["Promotion"] = raw["Classic_Absolute_Rank"] - raw["Evidence_Absolute_Rank"]
                raw_rank_rows.append(raw)
                degree_rows: list[pd.DataFrame] = []; pr_rows: list[pd.DataFrame] = []; annotation_rows: list[pd.DataFrame] = []; matched_rows: list[pd.DataFrame] = []
                _append_bias_rows(ranks=ranks, tissue=tissue, target=target_symbol, condition=condition,
                                  degree_rows=degree_rows, pagerank_rows=pr_rows, annotation_rows=annotation_rows, matched_rows=matched_rows)
                all_bias.extend([*degree_rows, *pr_rows]); all_annotation.extend(annotation_rows); all_matched.extend(matched_rows)
                for label, reference in (("Classic", classic), ("Current_Evidence", cache[tuple(NON_TEXT_CHANNELS)]["responses"][target_symbol])):
                    row = _benchmark(reference, ablated, label=label)
                    row.update({"Condition": condition.identifier, "Condition_Description": condition.description, "Tissue": tissue, "Target": target_symbol})
                    benchmarks.append(row)
                incident = ranks.loc[ranks["gene"] == target]
                # The perturbed target is excluded from signed output; target-local evidence comes from frozen properties.
                local = result["properties"].loc[result["properties"]["gene"] == target].copy()
                if not local.empty:
                    local_row = {"Condition": condition.identifier, "Tissue": tissue, "Target": target_symbol,
                                 "Target_Incident_Edges": int(local.iloc[0]["Incident_Edges"]),
                                 "Target_Mean_Ablated_Confidence": float(local.iloc[0]["Mean_Ablated_Confidence"]),
                                 "Target_Median_Ablated_Confidence": float(local.iloc[0]["Median_Ablated_Confidence"])}
                    for channel in NON_TEXT_CHANNELS:
                        local_row[f"{channel}_Availability"] = float(local.iloc[0][f"{channel}_Availability"])
                    target_local.append(local_row)
                if condition.identifier == "C1":
                    control_reproduced &= bool((ranks["Meaningful_Response"].sum() > 0))
                current_ranks = _ranks(classic, cache[tuple(NON_TEXT_CHANNELS)]["responses"][target_symbol])[["gene", "Evidence_Response", "Evidence_Absolute_Rank"]].rename(columns={"Evidence_Response": "Current_Evidence_Response", "Evidence_Absolute_Rank": "Current_Evidence_Rank"})
                candidates = _domain_frame(ranks, "absolute", "meaningful").merge(current_ranks, on="gene", how="left", validate="one_to_one")
                candidates = candidates.loc[(candidates["Classic_Rank"] > 100) & (candidates["Evidence_Rank"] <= 100)].copy()
                if not candidates.empty:
                    candidates["Condition"] = condition.identifier; candidates["Condition_Description"] = condition.description
                    candidates["Tissue"] = tissue; candidates["Target"] = target_symbol
                    keep = ["Condition", "Condition_Description", "Tissue", "Target", "gene", "Symbol", "Classic_Response", "Classic_Rank",
                            "Current_Evidence_Response", "Current_Evidence_Rank", "Evidence_Response", "Evidence_Rank", "Promotion",
                            "Classic_Degree", "Classic_PageRank", "Median_Official_Combined", "Median_S_nontext", "Median_Ablated_Confidence",
                            "Median_Text_Dependency", "Physical_Support", "Annotation_Coverage", *[f"Median_{channel}" for channel in NON_TEXT_CHANNELS]]
                    transitions.append(candidates[[column for column in keep if column in candidates]].sort_values(["Evidence_Rank", "Promotion"], ascending=[True, False]).head(100))
            if key not in retained_cache_keys:
                cache.pop(key, None)
                gc.collect()

    bias = pd.concat(all_bias, ignore_index=True)
    annotation = pd.concat(all_annotation, ignore_index=True)
    matched = pd.concat(all_matched, ignore_index=True)
    benchmark_frame = pd.DataFrame(benchmarks)
    health = pd.DataFrame(health_rows)
    local = pd.DataFrame(target_local)
    transition_frame = pd.concat(transitions, ignore_index=True) if transitions else pd.DataFrame()
    raw_rank_frame = pd.concat(raw_rank_rows, ignore_index=True)
    summary_conditions = _condition_summary(bias, matched, annotation)
    summary = {
        "status": "COMPLETE", "audit_version": "evidence-channel-ablation-1.0", "control_reproduced": bool(control_reproduced),
        "panel": {tissue: list(symbols) for tissue, symbols in PANEL.items()}, "topology": "E1_WEIGHT_ONLY; fixed nodes and edges for every condition",
        "current_control": {"base": "S_nontext^2", "text_multiplier": "OFF", "physical_multiplier": "OFF", "corum_pairwise": "OFF", "corum_complex_analysis": "post_calculation_only"},
        "conditions": [{"id": item.identifier, "description": item.description, "family": item.family, "included_channels": list(item.included_channels), "excluded_channels": list(item.excluded_channels)} for item in conditions],
        "string_combiner": index_metadata(), "condition_summary": summary_conditions,
        "artifacts": {"summary": "evidence_channel_ablation_summary.json", "bias": "evidence_channel_ablation_bias.csv", "raw_ranks": "evidence_channel_ablation_raw_node_ranks.csv", "benchmarks": "evidence_channel_ablation_benchmarks.csv", "availability_degree": "evidence_channel_availability_by_degree.csv", "availability_pagerank": "evidence_channel_availability_by_pagerank.csv", "matched": "evidence_channel_matched_controls.csv", "transitions": "evidence_channel_candidate_transitions.csv", "health": "evidence_channel_ablation_network_health.csv", "target_local": "evidence_channel_target_local_effect.csv"},
        "motor_mutation": {"classic": False, "directed": False, "evidence_default": False},
        "balanced_variant_created": False,
    }
    bias.to_csv(output_dir / "evidence_channel_ablation_bias.csv", index=False, encoding="utf-8-sig")
    raw_rank_frame.to_csv(output_dir / "evidence_channel_ablation_raw_node_ranks.csv", index=False, encoding="utf-8-sig")
    benchmark_frame.to_csv(output_dir / "evidence_channel_ablation_benchmarks.csv", index=False, encoding="utf-8-sig")
    pd.concat(degree_availability, ignore_index=True).to_csv(output_dir / "evidence_channel_availability_by_degree.csv", index=False, encoding="utf-8-sig")
    pd.concat(pagerank_availability, ignore_index=True).to_csv(output_dir / "evidence_channel_availability_by_pagerank.csv", index=False, encoding="utf-8-sig")
    matched.to_csv(output_dir / "evidence_channel_matched_controls.csv", index=False, encoding="utf-8-sig")
    transition_frame.to_csv(output_dir / "evidence_channel_candidate_transitions.csv", index=False, encoding="utf-8-sig")
    annotation.to_csv(output_dir / "evidence_channel_ablation_annotation.csv", index=False, encoding="utf-8-sig")
    health.to_csv(output_dir / "evidence_channel_ablation_network_health.csv", index=False, encoding="utf-8-sig")
    local.to_csv(output_dir / "evidence_channel_target_local_effect.csv", index=False, encoding="utf-8-sig")
    (output_dir / "evidence_channel_ablation_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("[ablation] artifacts written", flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "evidence_validation" / "channel_ablation")
    args = parser.parse_args()
    print(json.dumps({"status": run(args.output)["status"]}))
