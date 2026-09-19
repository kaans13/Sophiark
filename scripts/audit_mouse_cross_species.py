from __future__ import annotations

import json
import logging
import sys
import tempfile
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic as human_biology
from src import main as human_main
from src.cross_species import ComparisonConfig, OrthologyIndex, compare_ortholog_responses
from src.mouse import biology_logic as mouse_biology
from src.mouse import main as mouse_main
from src.mouse.config import MOUSE_SYMBOL_MAP


OUT = ROOT / "outputs_mouse" / "audits" / "mouse_cross_species_benchmark.json"
TARGETS = ("Cftr", "Trp53", "Egfr", "Myc", "Pten", "Stat3", "Brca1", "Nfkb1")


def run_target(module, graph, scores, target, out_dir, symbol_map):
    started = time.perf_counter()
    report, metrics = module.run_infection_simulation(
        graph, scores, out_dir, spesifik_hedefler=[target],
        redistribution_mode="top_n", exploratory_top_n=20, null_iterations=0,
        compute_structural_metrics=True, efficiency_sample_sources=32,
        local_efficiency_sample_nodes=32, tissue="lung",
        return_comparison_metrics=True, comparison_genes=[target],
    )
    signed = module.latest_signed_redistribution()
    candidates = signed[signed["gene"] != target].copy()
    values = pd.to_numeric(candidates["Delta_PageRank_Pct"], errors="coerce")
    positive = candidates.loc[values.idxmax()]
    negative = candidates.loc[values.idxmin()]
    primary = report.loc[report["gene"] == target].iloc[0]
    observed = metrics["observed"][target]
    return {
        "protein_id": target,
        "degree": graph.degree(target),
        "baseline_pagerank": observed["pagerank_before"],
        "perturbed_pagerank": observed["pagerank_after"],
        "target_delta_pct": observed["delta_pagerank_pct"],
        "strongest_positive": symbol_map.get(str(positive["gene"]), str(positive["gene"])),
        "strongest_positive_pct": float(positive["Delta_PageRank_Pct"]),
        "strongest_negative": symbol_map.get(str(negative["gene"]), str(negative["gene"])),
        "strongest_negative_pct": float(negative["Delta_PageRank_Pct"]),
        "target_local_efficiency_change_pct": float(primary["Target_Local_Efficiency_Change_Pct"]),
        "runtime_seconds": time.perf_counter() - started,
    }


def main() -> None:
    logging.disable(logging.CRITICAL)
    mouse_symbol_to_id = {str(symbol).casefold(): protein for protein, symbol in MOUSE_SYMBOL_MAP.items()}
    human_symbols = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    human_id_to_symbol = dict(zip(human_symbols["gene"].astype(str), human_symbols["Symbol"].astype(str)))
    human_symbol_to_id = {symbol.casefold(): protein for protein, symbol in human_id_to_symbol.items()}

    prepare_started = time.perf_counter()
    mouse_ids = [mouse_symbol_to_id[symbol.casefold()] for symbol in TARGETS]
    mouse_graph, mouse_scores = mouse_main.arayuz_icin_motoru_hazirla(
        forced_genes=mouse_ids, hedef_doku="lung", bc_sample_sources=16,
    )
    human_target = human_symbol_to_id["tp53"]
    human_graph, human_scores = human_main.arayuz_icin_motoru_hazirla(
        forced_genes=[human_target], hedef_doku="Lung", bc_sample_sources=16,
    )
    preparation_seconds = time.perf_counter() - prepare_started

    with tempfile.TemporaryDirectory() as directory:
        base = Path(directory)
        smoke = {
            symbol: run_target(
                mouse_biology, mouse_graph, mouse_scores, mouse_symbol_to_id[symbol.casefold()],
                base / symbol, MOUSE_SYMBOL_MAP,
            )
            for symbol in TARGETS
        }
        human_tp53 = run_target(
            human_biology, human_graph, human_scores, human_target, base / "human_tp53",
            human_id_to_symbol,
        )
        human_signed = human_biology.latest_signed_redistribution().copy()
        human_signed["Symbol"] = human_signed["gene"].map(human_id_to_symbol)

        trp53_id = mouse_symbol_to_id["trp53"]
        mouse_biology.run_infection_simulation(
            mouse_graph, mouse_scores, base / "mouse_trp53_compare",
            spesifik_hedefler=[trp53_id], redistribution_mode="top_n",
            exploratory_top_n=20, null_iterations=0, compute_structural_metrics=False,
            tissue="lung",
        )
        mouse_signed = mouse_biology.latest_signed_redistribution().copy()
        mouse_signed["Symbol"] = mouse_signed["gene"].map(MOUSE_SYMBOL_MAP)

        orthology = OrthologyIndex.from_snapshot(
            ROOT / "data" / "orthology" / "ensembl_human_mouse_orthology.tsv",
            ROOT / "data" / "orthology" / "ensembl_human_mouse_orthology.metadata.json",
        )
        config = ComparisonConfig(500, .85, .001, "top_n", "sophiark-classic-v5.3")
        comparison = compare_ortholog_responses(
            human_signed, mouse_signed, orthology, human_config=config, mouse_config=config,
        )

        multi_ids = [mouse_symbol_to_id[name] for name in ("trp53", "egfr", "myc")]
        multi_report, multi_metrics = mouse_biology.run_infection_simulation(
            mouse_graph, mouse_scores, base / "mouse_multi", spesifik_hedefler=multi_ids,
            redistribution_mode="top_n", exploratory_top_n=20, null_iterations=0,
            compute_structural_metrics=False, return_comparison_metrics=True, tissue="lung",
        )

    payload = {
        "mouse_graph": {"nodes": mouse_graph.vcount(), "edges": mouse_graph.ecount(), "tissue": "lung"},
        "human_graph": {"nodes": human_graph.vcount(), "edges": human_graph.ecount(), "tissue": "Lung"},
        "preparation_seconds": preparation_seconds,
        "mouse_classic_smoke": smoke,
        "human_tp53": human_tp53,
        "human_mouse_tp53_comparison": {
            "human_target": "TP53", "mouse_target": "Trp53",
            "metrics": comparison["metrics"],
        },
        "mouse_multi_target": {
            "requested": ["Trp53", "Egfr", "Myc"], "mapped": len(multi_ids),
            "target_ids": multi_ids, "report_rows": len(multi_report),
            "blocked_unique_edges": multi_metrics["blocked_edges"],
        },
        "orthology": orthology.statistics(),
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False), encoding="utf-8")
    print(json.dumps(payload, ensure_ascii=False))


if __name__ == "__main__":
    main()
