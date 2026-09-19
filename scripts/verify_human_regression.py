from __future__ import annotations

import argparse
import json
import logging
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic, main  # noqa: E402


TARGETS = ("CFTR", "TP53", "EGFR", "MYC", "BAD")
EFFICIENCY_COLUMNS = (
    "Global_Efficiency_Baseline", "Global_Efficiency_Perturbed",
    "Global_Efficiency_Change_Pct", "Mean_Local_Efficiency_Baseline",
    "Mean_Local_Efficiency_Perturbed", "Mean_Local_Efficiency_Change_Pct",
    "Target_Local_Efficiency_Baseline", "Target_Local_Efficiency_Perturbed",
    "Target_Local_Efficiency_Change_Pct",
)


def snapshot() -> dict:
    symbols = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbol_to_id = {str(symbol).upper(): str(gene) for gene, symbol in zip(symbols["gene"], symbols["Symbol"])}
    identifiers = [symbol_to_id[symbol] for symbol in TARGETS]
    graph, scores = main.arayuz_icin_motoru_hazirla(
        forced_genes=identifiers, hedef_doku="Lung", bc_sample_sources=4,
    )
    results = {"graph": {"nodes": graph.vcount(), "edges": graph.ecount()}, "targets": {}}
    with tempfile.TemporaryDirectory() as directory:
        for symbol, identifier in zip(TARGETS, identifiers):
            report, metrics = biology_logic.run_infection_simulation(
                graph, scores, Path(directory) / symbol, spesifik_hedefler=[identifier],
                redistribution_mode="top_n", exploratory_top_n=20, null_iterations=0,
                compute_structural_metrics=True, efficiency_sample_sources=256,
                local_efficiency_sample_nodes=128, tissue="Lung",
                return_comparison_metrics=True, comparison_genes=[identifier],
            )
            signed = biology_logic.latest_signed_redistribution()
            candidates = signed[signed["gene"] != identifier].copy()
            positive = candidates[candidates["Delta_PageRank_Pct"] > 0].nlargest(20, "Delta_PageRank_Pct")
            negative = candidates[candidates["Delta_PageRank_Pct"] < 0].nsmallest(20, "Delta_PageRank_Pct")
            target_row = report[report["gene"] == identifier].iloc[0]
            score_row = scores[scores["gene"] == identifier].iloc[0]
            results["targets"][symbol] = {
                "id": identifier,
                "observed": metrics["observed"][identifier],
                "strongest_positive": {"id": str(positive.iloc[0]["gene"]), "pct": float(positive.iloc[0]["Delta_PageRank_Pct"])},
                "strongest_negative": {"id": str(negative.iloc[0]["gene"]), "pct": float(negative.iloc[0]["Delta_PageRank_Pct"])},
                "positive_top20": positive["gene"].astype(str).tolist(),
                "negative_top20": negative["gene"].astype(str).tolist(),
                "report_rows": len(report),
                "bc": float(score_row["BC_Skoru"]),
                "hinterland": float(score_row["Hinterland_Skoru"]),
                "gateway": bool(score_row["Gümrük_Kapisi"]),
                "efficiency": {column: float(target_row[column]) for column in EFFICIENCY_COLUMNS},
            }
    return results


def compare(reference: dict, current: dict) -> list[str]:
    differences: list[str] = []

    def walk(left, right, path="root"):
        if isinstance(left, dict) and isinstance(right, dict):
            if set(left) != set(right):
                differences.append(f"{path}: keys differ")
                return
            for key in left:
                walk(left[key], right[key], f"{path}.{key}")
        elif isinstance(left, list) and isinstance(right, list):
            if left != right:
                differences.append(f"{path}: list differs")
        elif isinstance(left, (int, float)) and isinstance(right, (int, float)):
            if not np.isclose(left, right, rtol=1e-10, atol=1e-14, equal_nan=True):
                differences.append(f"{path}: {left!r} != {right!r}")
        elif left != right:
            differences.append(f"{path}: {left!r} != {right!r}")

    walk(reference, current)
    return differences


def main_cli() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--reference", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    logging.disable(logging.CRITICAL)
    current = snapshot()
    if args.output:
        args.output.write_text(json.dumps(current, indent=2), encoding="utf-8")
    differences = compare(json.loads(args.reference.read_text(encoding="utf-8")), current)
    print(json.dumps({"status": "PASS" if not differences else "FAIL", "difference_count": len(differences), "differences": differences[:20]}))
    raise SystemExit(0 if not differences else 1)


if __name__ == "__main__":
    main_cli()
