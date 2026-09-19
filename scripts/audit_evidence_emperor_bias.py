"""Cross-target Emperor/Hub-bias validation for the unchanged Evidence BETA."""
from __future__ import annotations

import argparse
import json
import logging
import math
from pathlib import Path
import sys
from typing import Iterable

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic  # noqa: E402
from src.evidence.config import EvidenceConfig  # noqa: E402
from src.evidence.engine import _edge_frame, prepare_evidence_network, run_evidence_simulation  # noqa: E402


PANEL: dict[str, tuple[str, ...]] = {
    "Lung": ("CFTR", "TP53", "EGFR", "GLP1R"),
    "Liver": ("TP53", "MYC", "PTEN", "STAT3"),
}
MEANINGFUL_RESPONSE_PCT = 0.05  # Existing production candidate threshold.
TOP_CUTOFFS = (20, 50, 100)
CHANNEL_COLUMNS = (
    "neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
    "coexpression", "coexpression_transferred", "experiments",
    "experiments_transferred", "database", "database_transferred",
    "textmining", "textmining_transferred",
)


def _rho(left: pd.Series, right: pd.Series) -> float | None:
    values = pd.concat((pd.to_numeric(left, errors="coerce"), pd.to_numeric(right, errors="coerce")), axis=1).dropna()
    if len(values) < 20 or values.iloc[:, 0].nunique() < 2 or values.iloc[:, 1].nunique() < 2:
        return None
    return float(values.iloc[:, 0].corr(values.iloc[:, 1], method="spearman"))


def _decile(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    return pd.qcut(numeric.rank(method="first"), 10, labels=[f"D{i}" for i in range(1, 11)]).astype(str)


def _top_fraction(values: pd.Series, fraction: float) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    rank = numeric.rank(method="first", ascending=False)
    return rank <= math.ceil(numeric.notna().sum() * fraction)


def _node_properties(graph, scores: pd.DataFrame, symbol_map: dict[str, str], complex_members: dict[str, list]) -> pd.DataFrame:
    edges = _edge_frame(graph)
    available = [column for column in CHANNEL_COLUMNS if column in edges]
    numeric = edges[available].apply(pd.to_numeric, errors="coerce").fillna(0)
    edges["Populated_Channel_Count"] = numeric.gt(0).sum(axis=1)
    edges["Has_Experimental"] = numeric[[c for c in ("experiments", "experiments_transferred") if c in numeric]].gt(0).any(axis=1)
    edges["Has_Database"] = numeric[[c for c in ("database", "database_transferred") if c in numeric]].gt(0).any(axis=1)
    transfer = [column for column in available if column.endswith("_transferred")]
    edges["Has_Transferred"] = numeric[transfer].gt(0).any(axis=1)
    official = pd.to_numeric(edges["combined_score"], errors="coerce") / 1000.0
    edges["S_nontext_to_Official"] = pd.to_numeric(edges["s_nontext"], errors="coerce") / official.replace(0, np.nan)
    first = edges.rename(columns={"protein1": "gene"})
    second = edges.rename(columns={"protein2": "gene"})
    incident = pd.concat((first, second), ignore_index=True)
    node = incident.groupby("gene", as_index=False).agg(
        Incident_Edges=("gene", "size"),
        Incident_Edges_With_Populated_STRING=("Populated_Channel_Count", lambda s: int((s > 0).sum())),
        Mean_S_nontext=("s_nontext", "mean"), Median_S_nontext=("s_nontext", "median"),
        Mean_Snontext_to_Official=("S_nontext_to_Official", "mean"), Median_Snontext_to_Official=("S_nontext_to_Official", "median"),
        Mean_Nonzero_Evidence_Channels=("Populated_Channel_Count", "mean"),
        Experimental_Availability=("Has_Experimental", "mean"), Database_Availability=("Has_Database", "mean"),
        Transferred_Availability=("Has_Transferred", "mean"), Median_Text_Dependency=("text_dependency", "median"),
        Physical_Availability=("physical_status", lambda s: float(s.isin(("SUPPORTED", "NOT_SUPPORTED")).mean())),
        Physical_Support=("physical_status", lambda s: float((s == "SUPPORTED").mean())),
        Median_Official_Combined=("combined_score", lambda s: float(pd.to_numeric(s, errors="coerce").median() / 1000.0)),
    )
    node["Symbol"] = node["gene"].map(symbol_map)
    node["CORUM_Membership_Count"] = node["Symbol"].map(lambda symbol: len(complex_members.get(str(symbol), ())))
    base_columns = [column for column in ("gene", "Classic_Degree", "Classic_PageRank", "pagerank_raw", "Hinterland_Skoru") if column in scores]
    node = node.merge(scores[base_columns].drop_duplicates("gene"), on="gene", how="left", validate="one_to_one")
    node["Degree_Decile"] = _decile(node["Classic_Degree"])
    node["PageRank_Decile"] = _decile(node["Classic_PageRank"])
    node["Degree_Top10pct"] = _top_fraction(node["Classic_Degree"], .10)
    node["Degree_Top5pct"] = _top_fraction(node["Classic_Degree"], .05)
    node["Degree_Top1pct"] = _top_fraction(node["Classic_Degree"], .01)
    node["PageRank_Top10pct"] = _top_fraction(node["Classic_PageRank"], .10)
    node["PageRank_Top5pct"] = _top_fraction(node["Classic_PageRank"], .05)
    node["PageRank_Top1pct"] = _top_fraction(node["Classic_PageRank"], .01)
    node["Annotation_Coverage"] = node["Mean_Nonzero_Evidence_Channels"]
    node["Annotation_Top10pct"] = _top_fraction(node["Annotation_Coverage"], .10)
    node["Emperor_Audit_Group"] = node["Degree_Top10pct"] & node["PageRank_Top10pct"] & node["Annotation_Top10pct"]
    return node


def _ranks(classic: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    left = classic[["gene", "Delta_PageRank_Pct"]].rename(columns={"Delta_PageRank_Pct": "Classic_Response"})
    right = evidence[["gene", "Delta_PageRank_Pct"]].rename(columns={"Delta_PageRank_Pct": "Evidence_Response"})
    frame = left.merge(right, on="gene", how="inner", validate="one_to_one")
    for prefix, column in (("Classic", "Classic_Response"), ("Evidence", "Evidence_Response")):
        values = pd.to_numeric(frame[column], errors="coerce")
        frame[f"{prefix}_Absolute_Rank"] = values.abs().rank(method="min", ascending=False)
        frame[f"{prefix}_Positive_Rank"] = values.where(values > 0).rank(method="min", ascending=False)
        frame[f"{prefix}_Negative_Rank"] = values.abs().where(values < 0).rank(method="min", ascending=False)
    frame["Response_Delta"] = frame["Evidence_Response"] - frame["Classic_Response"]
    frame["Meaningful_Response"] = frame[["Classic_Response", "Evidence_Response"]].abs().max(axis=1) >= MEANINGFUL_RESPONSE_PCT
    return frame


def _domain_frame(frame: pd.DataFrame, domain: str, scope: str) -> pd.DataFrame:
    if domain == "absolute":
        result = frame.copy()
        result["Classic_Rank"] = result["Classic_Absolute_Rank"]
        result["Evidence_Rank"] = result["Evidence_Absolute_Rank"]
    elif domain == "positive":
        result = frame.loc[(frame["Classic_Response"] > 0) & (frame["Evidence_Response"] > 0)].copy()
        result["Classic_Rank"] = result["Classic_Positive_Rank"]
        result["Evidence_Rank"] = result["Evidence_Positive_Rank"]
    else:
        result = frame.loc[(frame["Classic_Response"] < 0) & (frame["Evidence_Response"] < 0)].copy()
        result["Classic_Rank"] = result["Classic_Negative_Rank"]
        result["Evidence_Rank"] = result["Evidence_Negative_Rank"]
    if scope == "meaningful":
        result = result.loc[result["Meaningful_Response"]].copy()
    result["Promotion"] = result["Classic_Rank"] - result["Evidence_Rank"]
    for cutoff in TOP_CUTOFFS:
        result[f"Entered_Evidence_Top{cutoff}"] = (result["Evidence_Rank"] <= cutoff) & (result["Classic_Rank"] > cutoff)
        result[f"Left_Classic_Top{cutoff}"] = (result["Classic_Rank"] <= cutoff) & (result["Evidence_Rank"] > cutoff)
    return result


def _strata_summary(frame: pd.DataFrame, group_column: str, *, tissue: str, target: str,
                    rank_domain: str, scope: str, audit: str) -> pd.DataFrame:
    summary = frame.groupby(group_column, observed=True).agg(
        Node_Count=("gene", "size"), Median_Promotion=("Promotion", "median"), Mean_Promotion=("Promotion", "mean"),
        Median_Response_Delta=("Response_Delta", "median"), Fraction_Promoted=("Promotion", lambda s: float((s > 0).mean())),
        **{f"Fraction_Entering_Evidence_Top{cutoff}": (f"Entered_Evidence_Top{cutoff}", "mean") for cutoff in TOP_CUTOFFS},
        **{f"Fraction_Leaving_Classic_Top{cutoff}": (f"Left_Classic_Top{cutoff}", "mean") for cutoff in TOP_CUTOFFS},
    ).reset_index().rename(columns={group_column: "Stratum"})
    summary.insert(0, "Audit", audit); summary.insert(1, "Tissue", tissue); summary.insert(2, "Target", target)
    summary.insert(3, "Rank_Domain", rank_domain); summary.insert(4, "Response_Scope", scope)
    return summary


def _tail_rows(frame: pd.DataFrame, *, source: str, tissue: str, target: str,
               rank_domain: str, scope: str) -> pd.DataFrame:
    definitions = (("Low_D1_D3", frame[source].isin(("D1", "D2", "D3"))), ("Middle_D4_D7", frame[source].isin(("D4", "D5", "D6", "D7"))),
                   ("High_D8_D9", frame[source].isin(("D8", "D9"))), ("D10", frame[source].eq("D10")))
    prefix = "Degree" if source == "Degree_Decile" else "PageRank"
    for fraction, label in ((.10, "Top10pct"), (.05, "Top5pct"), (.01, "Top1pct")):
        values = frame[f"Classic_{prefix}"] if f"Classic_{prefix}" in frame else frame["Classic_Degree" if prefix == "Degree" else "Classic_PageRank"]
        definitions += ((label, _top_fraction(values, fraction)),)
    rows = []
    for label, mask in definitions:
        selected = frame.loc[mask]
        if selected.empty:
            continue
        row = {
            "Audit": f"{prefix}_Tail", "Tissue": tissue, "Target": target, "Rank_Domain": rank_domain,
            "Response_Scope": scope, "Stratum": label, "Node_Count": len(selected),
            "Median_Promotion": float(selected["Promotion"].median()), "Mean_Promotion": float(selected["Promotion"].mean()),
            "Median_Response_Delta": float(selected["Response_Delta"].median()),
            "Fraction_Promoted": float((selected["Promotion"] > 0).mean()),
        }
        for cutoff in TOP_CUTOFFS:
            row[f"Fraction_Entering_Evidence_Top{cutoff}"] = float(selected[f"Entered_Evidence_Top{cutoff}"].mean())
            row[f"Fraction_Leaving_Classic_Top{cutoff}"] = float(selected[f"Left_Classic_Top{cutoff}"].mean())
        rows.append(row)
    return pd.DataFrame(rows)


def _annotation_correlations(frame: pd.DataFrame, *, tissue: str, target: str) -> list[dict[str, object]]:
    metrics = (
        "Incident_Edges_With_Populated_STRING", "Mean_S_nontext", "Median_S_nontext", "Mean_Nonzero_Evidence_Channels",
        "Experimental_Availability", "Database_Availability", "Transferred_Availability", "Median_Text_Dependency",
        "Physical_Availability", "Physical_Support", "CORUM_Membership_Count", "Annotation_Coverage",
    )
    rows = []
    for scope, values in (("all_nodes", frame), ("meaningful_response", frame.loc[frame["Meaningful_Response"]])):
        for metric in metrics:
            rows.append({"Audit": "Promotion_Correlation", "Tissue": tissue, "Target": target, "Response_Scope": scope,
                         "Metric": metric, "Spearman_Rho": _rho(values["Promotion"], values[metric]), "Node_Count": len(values)})
    return rows


def _matched_controls(frame: pd.DataFrame, *, tissue: str, target: str) -> pd.DataFrame:
    rows = []
    for scope, values in (("all_nodes", frame), ("meaningful_response", frame.loc[frame["Meaningful_Response"]])):
        values = values.copy()
        values["Match_Bin"] = values["Degree_Decile"].astype(str) + "|" + values["PageRank_Decile"].astype(str)
        for match_bin, group in values.groupby("Match_Bin", observed=True):
            if len(group) < 10 or group["Annotation_Coverage"].nunique() < 2:
                continue
            median = group["Annotation_Coverage"].median()
            high, low = group.loc[group["Annotation_Coverage"] > median], group.loc[group["Annotation_Coverage"] <= median]
            if high.empty or low.empty:
                continue
            rows.append({"Audit": "Coverage_Matched_Control", "Tissue": tissue, "Target": target, "Response_Scope": scope,
                         "Match_Bin": match_bin, "High_Annotation_N": len(high), "Low_Annotation_N": len(low),
                         "High_Annotation_Median_Promotion": float(high["Promotion"].median()),
                         "Low_Annotation_Median_Promotion": float(low["Promotion"].median()),
                         "Promotion_Difference_High_Minus_Low": float(high["Promotion"].median() - low["Promotion"].median())})
        emperor = values.loc[values["Emperor_Audit_Group"]]
        controls = values.loc[values["Degree_Decile"].eq("D10") & values["PageRank_Decile"].eq("D10") & ~values["Emperor_Audit_Group"]]
        if not emperor.empty and not controls.empty:
            rows.append({"Audit": "Emperor_vs_Matched_Control", "Tissue": tissue, "Target": target, "Response_Scope": scope,
                         "Match_Bin": "D10|D10", "High_Annotation_N": len(emperor), "Low_Annotation_N": len(controls),
                         "High_Annotation_Median_Promotion": float(emperor["Promotion"].median()),
                         "Low_Annotation_Median_Promotion": float(controls["Promotion"].median()),
                         "Promotion_Difference_High_Minus_Low": float(emperor["Promotion"].median() - controls["Promotion"].median())})
    return pd.DataFrame(rows)


def _aggregate_target_rows(frame: pd.DataFrame) -> pd.DataFrame:
    metric_columns = [column for column in frame if column.startswith(("Median_", "Mean_", "Fraction_"))]
    grouped = frame.groupby(["Audit", "Rank_Domain", "Response_Scope", "Stratum"], observed=True)
    aggregate = grouped[metric_columns].mean().reset_index()
    aggregate.insert(1, "Tissue", "ALL"); aggregate.insert(2, "Target", "ALL_TARGET_MEAN")
    aggregate.insert(5, "Target_Conditions", grouped.size().to_numpy())
    return aggregate


def _decision(decile: pd.DataFrame, annotation: pd.DataFrame, matched: pd.DataFrame) -> dict[str, object]:
    selected = decile.loc[(decile["Target"] != "ALL_TARGET_MEAN") & (decile["Rank_Domain"] == "absolute") & (decile["Response_Scope"] == "meaningful")]
    degree = selected.loc[(selected["Audit"] == "Degree") & (selected["Stratum"] == "D10")]
    pagerank = selected.loc[(selected["Audit"] == "PageRank") & (selected["Stratum"] == "D10")]
    def signs(values: pd.Series) -> dict[str, int]:
        return {"promoted": int((values > 0).sum()), "demoted": int((values < 0).sum()), "neutral": int((values == 0).sum()), "conditions": len(values)}
    degree_signs, pr_signs = signs(degree["Median_Promotion"]), signs(pagerank["Median_Promotion"])
    coverage = annotation.loc[(annotation["Metric"] == "Annotation_Coverage") & (annotation["Response_Scope"] == "meaningful_response")]
    coverage_signs = signs(coverage["Spearman_Rho"].dropna())
    matched_selected = matched.loc[
        (matched["Audit"] == "Coverage_Matched_Control")
        & (matched["Response_Scope"] == "meaningful_response")
        & (matched["Target"] != "ALL_TARGET_MEAN")
    ]
    matched_signs = signs(matched_selected.groupby(["Tissue", "Target"], observed=True)["Promotion_Difference_High_Minus_Low"].median())
    systematic = degree_signs["promoted"] >= 6 or pr_signs["promoted"] >= 6 or (coverage_signs["promoted"] >= 6 and matched_signs["promoted"] >= 6)
    return {
        "judgment": "BIAS_WARNING" if systematic else "PARTIAL",
        "rule": "BIAS_WARNING requires a same-direction pattern in at least 6/8 target conditions; otherwise PARTIAL.",
        "degree_d10_median_promotion_direction": degree_signs,
        "pagerank_d10_median_promotion_direction": pr_signs,
        "annotation_promotion_spearman_direction": coverage_signs,
        "matched_coverage_difference_direction": matched_signs,
    }


def run(output_dir: Path) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    logging.disable(logging.CRITICAL)
    mapping = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbol_map = dict(zip(mapping["gene"].astype(str), mapping["Symbol"].astype(str)))
    symbol_to_id = {symbol.upper(): gene for gene, symbol in symbol_map.items()}
    import pickle
    with (ROOT / "data" / "processed" / "complex_members.pkl").open("rb") as stream:
        complex_members = pickle.load(stream)
    config = EvidenceConfig()
    degree_rows: list[pd.DataFrame] = []; pr_rows: list[pd.DataFrame] = []; annotations: list[dict[str, object]] = []
    matched: list[pd.DataFrame] = []; transitions: list[pd.DataFrame] = []; snt_centrality: list[dict[str, object]] = []
    condition_metadata: list[dict[str, object]] = []
    for tissue, symbols in PANEL.items():
        print(f"[audit] preparing {tissue}", flush=True)
        targets = [symbol_to_id[symbol] for symbol in symbols]
        graph, scores, _ = prepare_evidence_network(config=config, forced_genes=targets, tissue=tissue, bc_sample_sources=4)
        properties = _node_properties(graph, scores, symbol_map, complex_members)
        for metric in ("Mean_S_nontext", "Median_S_nontext", "Mean_Snontext_to_Official", "Median_Snontext_to_Official"):
            for centrality in ("Classic_Degree", "Classic_PageRank"):
                snt_centrality.append({"Audit": "S_nontext_Centrality", "Tissue": tissue, "Target": "TISSUE_BASELINE", "Response_Scope": "all_nodes",
                                       "Metric": f"{metric} vs {centrality}", "Spearman_Rho": _rho(properties[metric], properties[centrality]), "Node_Count": len(properties)})
        for stratum, label in (("Degree_Decile", "Classic_Degree"), ("PageRank_Decile", "Classic_PageRank")):
            for bucket, values in properties.groupby(stratum, observed=True):
                for metric in ("Mean_S_nontext", "Median_S_nontext", "Mean_Snontext_to_Official", "Median_Snontext_to_Official"):
                    snt_centrality.append({
                        "Audit": "S_nontext_Decile_Distribution", "Tissue": tissue, "Target": "TISSUE_BASELINE",
                        "Response_Scope": "all_nodes", "Metric": f"{metric} by {label}", "Stratum": str(bucket),
                        "Node_Count": len(values), "Median_Value": float(values[metric].median()), "Mean_Value": float(values[metric].mean()),
                    })
        for target_symbol, target in zip(symbols, targets):
            print(f"[audit] {tissue}/{target_symbol}", flush=True)
            classic_graph = graph.copy()
            classic_graph.es["weight"] = classic_graph.es["classic_weight"]
            classic_graph.es["distance"] = (1.0 / np.asarray(classic_graph.es["weight"], dtype=float)).tolist()
            biology_logic.run_infection_simulation(
                classic_graph, scores, output_dir / "_classic_controls" / tissue / target_symbol,
                spesifik_hedefler=[target], exploratory_top_n=100, compute_structural_metrics=False,
                null_iterations=0, edge_evidence_weighting_mode="legacy_edge_modifiers", tissue=tissue,
            )
            classic = biology_logic.latest_signed_redistribution().copy(deep=True)
            evidence_result = run_evidence_simulation(
                graph=graph.copy(), scores=scores, out_dir=output_dir / "_evidence_runs" / tissue / target_symbol,
                targets=[target], config=config, tissue=tissue, block_weight_fraction=.001, damping=.85,
                scientific_options={"exploratory_top_n": 100, "compute_structural_metrics": False, "null_iterations": 0},
            )
            # Signed redistribution deliberately excludes the directly perturbed target.
            # Baseline strata must use that exact non-target universe as well.
            non_target_properties = properties.loc[properties["gene"] != target].copy()
            ranks = _ranks(classic, evidence_result.signed_response).merge(non_target_properties, on="gene", how="left", validate="one_to_one")
            if len(ranks) != len(non_target_properties):
                raise RuntimeError(f"Node-universe mismatch for {target_symbol}/{tissue}: {len(ranks)} != {len(non_target_properties)}")
            condition_metadata.append({"tissue": tissue, "target": target_symbol, "nodes": len(ranks), "meaningful_nodes": int(ranks["Meaningful_Response"].sum()),
                                       "text_multiplier_enabled": evidence_result.provenance["perturbation"]["actual_text_multiplier_enabled"],
                                       "physical_multiplier_enabled": evidence_result.provenance["perturbation"]["actual_physical_multiplier_enabled"]})
            for domain in ("absolute", "positive", "negative"):
                for scope in ("all_nodes", "meaningful"):
                    current = _domain_frame(ranks, domain, scope)
                    if current.empty:
                        continue
                    degree_rows.extend((_strata_summary(current, "Degree_Decile", tissue=tissue, target=target_symbol, rank_domain=domain, scope=scope, audit="Degree"),
                                        _tail_rows(current, source="Degree_Decile", tissue=tissue, target=target_symbol, rank_domain=domain, scope=scope)))
                    pr_rows.extend((_strata_summary(current, "PageRank_Decile", tissue=tissue, target=target_symbol, rank_domain=domain, scope=scope, audit="PageRank"),
                                    _tail_rows(current, source="PageRank_Decile", tissue=tissue, target=target_symbol, rank_domain=domain, scope=scope)))
                    if domain == "absolute" and scope == "all_nodes":
                        annotations.extend(_annotation_correlations(current, tissue=tissue, target=target_symbol))
                        matched.append(_matched_controls(current, tissue=tissue, target=target_symbol))
            absolute = _domain_frame(ranks, "absolute", "all_nodes")
            transition = absolute.loc[(absolute["Classic_Rank"] > 100) & (absolute["Evidence_Rank"] <= 100)].copy()
            transition["Tissue"] = tissue; transition["Target"] = target_symbol
            transition["Transition"] = "Classic_low_to_Evidence_high_absolute_top100"
            transitions.append(transition[[column for column in (
                "Tissue", "Target", "Transition", "gene", "Symbol", "Classic_Rank", "Evidence_Rank", "Promotion", "Classic_Response", "Evidence_Response",
                "Classic_Degree", "Classic_PageRank", "Mean_S_nontext", "Median_S_nontext", "Median_Official_Combined", "Median_Text_Dependency",
                "Experimental_Availability", "Database_Availability", "Physical_Support", "Annotation_Coverage", "CORUM_Membership_Count",
                "Degree_Decile", "PageRank_Decile", "Degree_Top10pct", "PageRank_Top10pct", "Annotation_Top10pct", "Emperor_Audit_Group",
            ) if column in transition]])
    degree = pd.concat(degree_rows, ignore_index=True); pagerank = pd.concat(pr_rows, ignore_index=True)
    degree = pd.concat((degree, _aggregate_target_rows(degree)), ignore_index=True)
    pagerank = pd.concat((pagerank, _aggregate_target_rows(pagerank)), ignore_index=True)
    annotation = pd.DataFrame([*annotations, *snt_centrality])
    annotation_summary = annotation.groupby(["Audit", "Response_Scope", "Metric"], observed=True)["Spearman_Rho"].agg(["median", "mean", "count", lambda s: int((s > 0).sum()), lambda s: int((s < 0).sum())]).reset_index()
    annotation_summary.columns = ["Audit", "Response_Scope", "Metric", "Median_Target_Rho", "Mean_Target_Rho", "Target_Conditions", "Positive_Rho_Targets", "Negative_Rho_Targets"]
    annotation_summary.insert(1, "Tissue", "ALL"); annotation_summary.insert(2, "Target", "ALL_TARGET_MEAN")
    annotation = pd.concat((annotation, annotation_summary), ignore_index=True, sort=False)
    matched_frame = pd.concat(matched, ignore_index=True)
    matched_aggregate = matched_frame.groupby(["Audit", "Response_Scope"], observed=True)["Promotion_Difference_High_Minus_Low"].agg(["median", "mean", "count", lambda s: int((s > 0).sum()), lambda s: int((s < 0).sum())]).reset_index()
    matched_aggregate.columns = ["Audit", "Response_Scope", "Median_Target_Difference", "Mean_Target_Difference", "Target_Conditions", "Positive_Difference_Targets", "Negative_Difference_Targets"]
    matched_aggregate.insert(1, "Tissue", "ALL"); matched_aggregate.insert(2, "Target", "ALL_TARGET_MEAN")
    matched_frame = pd.concat((matched_frame, matched_aggregate), ignore_index=True, sort=False)
    transition_frame = pd.concat(transitions, ignore_index=True).sort_values(["Tissue", "Target", "Evidence_Rank", "Promotion"], ascending=[True, True, True, False])
    degree.to_csv(output_dir / "evidence_degree_decile_audit.csv", index=False, encoding="utf-8-sig")
    pagerank.to_csv(output_dir / "evidence_pagerank_decile_audit.csv", index=False, encoding="utf-8-sig")
    annotation.to_csv(output_dir / "evidence_annotation_bias.csv", index=False, encoding="utf-8-sig")
    matched_frame.to_csv(output_dir / "evidence_matched_controls.csv", index=False, encoding="utf-8-sig")
    transition_frame.to_csv(output_dir / "evidence_candidate_transitions.csv", index=False, encoding="utf-8-sig")
    decision = _decision(pd.concat((degree, pagerank), ignore_index=True, sort=False), annotation, matched_frame)
    summary = {
        "status": decision["judgment"], "audit_version": "emperor-bias-validation-1.0", "config": config.to_dict(),
        "rank_convention": "Promotion = ClassicRank - EvidenceRank; positive values are Evidence promotions.",
        "rank_domains": {"primary": "absolute |response|", "secondary": ["positive", "negative"]},
        "meaningful_response_threshold_pct": MEANINGFUL_RESPONSE_PCT,
        "panel": condition_metadata, "decision": decision,
        "artifacts": {"degree": "evidence_degree_decile_audit.csv", "pagerank": "evidence_pagerank_decile_audit.csv", "annotation": "evidence_annotation_bias.csv", "matched": "evidence_matched_controls.csv", "transitions": "evidence_candidate_transitions.csv"},
        "motor_mutation": {"classic": False, "directed": False, "evidence_default": False},
    }
    (output_dir / "evidence_emperor_bias_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print("[audit] artifacts written", flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "evidence_validation" / "emperor_bias")
    args = parser.parse_args()
    print(json.dumps({"status": run(args.output)["status"]}))
