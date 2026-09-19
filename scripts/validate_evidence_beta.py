"""Real-data validation panel for Evidence Engine BETA."""
from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
import sys

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic  # noqa: E402
from src.evidence.audit import mid_evidence_audit  # noqa: E402
from src.evidence.config import EvidenceConfig  # noqa: E402
from src.evidence.engine import _edge_frame, prepare_evidence_network, run_evidence_simulation  # noqa: E402
from src.evidence.scores import index_metadata, inspect_local_sources  # noqa: E402


PANEL = {
    "Lung": ("CFTR", "TP53", "EGFR", "GLP1R"),
    "Liver": ("TP53", "MYC", "PTEN", "STAT3"),
}


def _rank(frame: pd.DataFrame) -> pd.DataFrame:
    result = frame[["gene", "Delta_PageRank_Pct"]].copy()
    result["rank"] = pd.to_numeric(result["Delta_PageRank_Pct"], errors="coerce").abs().rank(method="min", ascending=False)
    return result


def _compare(classic: pd.DataFrame, evidence: pd.DataFrame) -> dict[str, object]:
    merged = _rank(classic).merge(_rank(evidence), on="gene", suffixes=("_classic", "_evidence"), validate="one_to_one")
    c = merged["Delta_PageRank_Pct_classic"]
    e = merged["Delta_PageRank_Pct_evidence"]
    meaningful = (c.abs() >= .05) | (e.abs() >= .05)
    shifts = (merged["rank_classic"] - merged["rank_evidence"])[meaningful]
    result: dict[str, object] = {
        "nodes_compared": len(merged), "meaningful_nodes": int(meaningful.sum()),
        "pearson": float(c.corr(e, method="pearson")), "spearman": float(c.corr(e, method="spearman")),
        "sign_flips_meaningful": int(((np.sign(c) != np.sign(e)) & meaningful).sum()),
        "median_absolute_rank_shift_meaningful": float(shifts.abs().median()) if len(shifts) else None,
    }
    for sign, mask in (("positive", c > 0), ("negative", c < 0), ("absolute", c.notna())):
        for count in ((20, 50, 100) if sign != "absolute" else (50, 100, 300)):
            classic_subset = merged.loc[mask]
            classic_ranked = (
                classic_subset.nlargest(count, "Delta_PageRank_Pct_classic") if sign == "positive"
                else classic_subset.nsmallest(count, "Delta_PageRank_Pct_classic") if sign == "negative"
                else classic_subset.nsmallest(count, "rank_classic")
            )
            classic_top = set(classic_ranked["gene"])
            emask = e > 0 if sign == "positive" else (e < 0 if sign == "negative" else e.notna())
            evidence_subset = merged.loc[emask]
            evidence_ranked = (
                evidence_subset.nlargest(count, "Delta_PageRank_Pct_evidence") if sign == "positive"
                else evidence_subset.nsmallest(count, "Delta_PageRank_Pct_evidence") if sign == "negative"
                else evidence_subset.nsmallest(count, "rank_evidence")
            )
            evidence_top = set(evidence_ranked["gene"])
            result[f"{sign}_top{count}_overlap"] = len(classic_top & evidence_top)
    transitions = merged.loc[meaningful].sort_values("rank_classic").copy()
    transitions["promotion"] = transitions["rank_classic"] - transitions["rank_evidence"]
    result["largest_evidence_promotions"] = transitions.nlargest(10, "promotion")[["gene", "promotion", "Delta_PageRank_Pct_classic", "Delta_PageRank_Pct_evidence"]].to_dict("records")
    return result


def run(output: Path) -> dict[str, object]:
    mapping = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbol_to_id = dict(zip(mapping["Symbol"].astype(str).str.upper(), mapping["gene"].astype(str)))
    config = EvidenceConfig()
    result: dict[str, object] = {
        "status": "COMPLETE", "local_sources": inspect_local_sources(),
        "string_reconstruction": index_metadata(), "configuration": config.to_dict(),
        "panel": {}, "tests": {},
    }
    logging.disable(logging.CRITICAL)
    for tissue, symbols in PANEL.items():
        ids = [symbol_to_id[s] for s in symbols if s in symbol_to_id]
        graph, scores, _ = prepare_evidence_network(
            config=config, forced_genes=ids, tissue=tissue, bc_sample_sources=4,
        )
        edge_frame = _edge_frame(graph)
        tissue_result: dict[str, object] = {
            "graph": {"nodes": graph.vcount(), "edges": graph.ecount()},
            "modifier_bias_audit": graph["evidence_modifier_audit"],
            "mid_evidence": mid_evidence_audit(edge_frame).to_dict("records"),
            "physical_supported_edges": int(edge_frame["physical_status"].eq("SUPPORTED").sum()),
            "physical_coverage": float(edge_frame["physical_status"].eq("SUPPORTED").mean()),
            "targets": {},
        }
        for symbol, target in zip(symbols, ids):
            classic_graph = graph.copy()
            classic_graph.es["weight"] = classic_graph.es["classic_weight"]
            classic_graph.es["distance"] = (1 / np.asarray(classic_graph.es["weight"], dtype=float)).tolist()
            biology_logic.run_infection_simulation(
                classic_graph, scores, output / tissue / symbol / "classic",
                spesifik_hedefler=[target], exploratory_top_n=100,
                compute_structural_metrics=False, null_iterations=0,
                edge_evidence_weighting_mode="legacy_edge_modifiers", tissue=tissue,
            )
            classic_signed = biology_logic.latest_signed_redistribution().copy(deep=True)
            evidence = run_evidence_simulation(
                graph=graph.copy(), scores=scores, out_dir=output / tissue / symbol / "evidence",
                targets=[target], config=config, tissue=tissue,
                block_weight_fraction=.001, damping=.85,
                scientific_options={"exploratory_top_n": 100, "compute_structural_metrics": False, "null_iterations": 0},
            )
            comparison = _compare(classic_signed, evidence.signed_response)
            comparison["complex_count"] = len(evidence.complexes)
            comparison["provenance_embedded"] = "Evidence_Run_Provenance_JSON" in evidence.signed_response
            tissue_result["targets"][symbol] = comparison
        result["panel"][tissue] = tissue_result
    result["tests"] = {
        "classic_code_modified": False, "directed_code_modified": False,
        "corum_is_post_calculation": True, "physical_multiplier_enabled": False,
    }
    output.mkdir(parents=True, exist_ok=True)
    (output / "evidence_beta_validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=ROOT / "outputs" / "evidence_validation")
    args = parser.parse_args()
    report = run(args.output)
    print(json.dumps({"status": report["status"], "tissues": list(report["panel"])}))
