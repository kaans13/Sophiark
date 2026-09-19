"""Run the real CFTR/Lung Classic and Directed presentation acceptance check."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import tempfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import main as motor
from src.biology_logic import latest_signed_redistribution
from src.interpretation import interpret_forward_results, select_redistribution_losses
from src.services.directed_analysis import run_optional_directed_engine
from src.services.enrichment import run_post_simulation_enrichment
from src.services.engine_presentation import compare_ui_frame, directed_ui_frame


CFTR = "ENSP00000003084"
ANO2 = "ENSP00000348453"
RFLNB = "ENSP00000331915"


def _terms(frame: pd.DataFrame) -> list[str]:
    return frame.head(5).get("Term", pd.Series(dtype=str)).astype(str).tolist()


def _json_record(frame: pd.DataFrame, gene: str) -> dict[str, object]:
    symbol_column = next((column for column in ("Gen", "Gene", "Symbol") if column in frame), None)
    if symbol_column is None:
        raise RuntimeError(f"Comparison display has no symbol column: {frame.columns.tolist()}")
    row = frame.loc[frame[symbol_column] == gene].head(1)
    if row.empty:
        raise RuntimeError(f"Comparison display has no record for {gene}")
    return json.loads(row.to_json(orient="records"))[0]


def main() -> None:
    symbols = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    symbol_map = dict(zip(symbols["gene"].astype(str), symbols["Symbol"].astype(str)))
    graph, scores = motor.arayuz_icin_motoru_hazirla(
        forced_genes=[CFTR], hedef_doku="Lung", bc_sample_sources=4,
    )
    with tempfile.TemporaryDirectory(prefix="sophiark_cftr_verify_") as output:
        classic_report = motor.run_infection_simulation(
            graph.copy(), scores, Path(output), spesifik_hedefler=[CFTR],
            exploratory_top_n=20, tissue="Lung", compute_structural_metrics=False,
        )
    classic_signed = latest_signed_redistribution()
    bundle = run_optional_directed_engine(
        graph=graph, targets=[CFTR], tissue="Lung", classic_report=classic_signed,
        project_root=ROOT, mode="Compare", block_weight_fraction=.001,
        damping=.85, bc_sample_sources=4, scores=scores,
        gene_to_symbol=symbol_map, candidate_limit=20,
    )
    if bundle.presentation is None:
        raise RuntimeError(bundle.calculation.error or "Directed presentation unavailable")
    directed = bundle.presentation
    if bundle.comparison is None:
        raise RuntimeError("Compare result unavailable")
    directed_display = directed_ui_frame(directed)
    compare_display = compare_ui_frame(bundle.comparison.table, directed)
    background = [symbol_map.get(str(gene), "") for gene in graph.vs["name"]]
    classic_candidates = classic_report[classic_report["Hasar_Tipi"] == motor.HASAR_UCUNCUL].copy()
    classic_candidates["Symbol"] = classic_candidates["gene"].map(symbol_map)
    classic_enrichment, classic_notices = run_post_simulation_enrichment(
        classic_candidates["Symbol"].dropna().tolist(), is_mouse=False,
        background_symbols=background,
    )
    directed_enrichment, directed_notices = run_post_simulation_enrichment(
        directed.candidates["gene_symbol"].replace("Unresolved", pd.NA).dropna().tolist(),
        is_mouse=False, background_symbols=background,
    )
    classic_interp = interpret_forward_results(
        report=classic_report, candidates=classic_candidates, targets=[CFTR],
        gene_to_symbol=symbol_map, tissue="Lung", organism="Homo sapiens",
        enrichment=classic_enrichment,
        loss_candidates=select_redistribution_losses(classic_signed),
    )
    directed_interp = interpret_forward_results(
        report=directed.report, candidates=directed.candidates, targets=[CFTR],
        gene_to_symbol=symbol_map, tissue="Lung", organism="Homo sapiens",
        enrichment=directed_enrichment, loss_candidates=directed.losses,
    )
    classic_values = classic_signed.set_index("gene")["Delta_PageRank_Pct"]
    directed_values = directed.full_response.set_index("gene")["Delta_PageRank_Pct"]
    classic_top = classic_candidates.sort_values("Delta_PageRank_Pct", ascending=False)["Symbol"].tolist()
    directed_top = [
        symbol if symbol != "Unresolved" else f"Unresolved ({protein_id})"
        for symbol, protein_id in zip(directed.candidates["gene_symbol"], directed.candidates["protein_id"])
    ]
    result = {
        "classic": {"ANO2": float(classic_values[ANO2]), "RFLNB": float(classic_values[RFLNB]), "top20": classic_top,
                    "go_kegg_top_terms": _terms(classic_enrichment), "notices": classic_notices,
                    "interpretation_top_candidate": classic_interp.top_candidates[0]["gene"]},
        "directed": {"ANO2": float(directed_values[ANO2]), "RFLNB": float(directed_values[RFLNB]), "top20": directed_top,
                     "go_kegg_top_terms": _terms(directed_enrichment), "notices": directed_notices,
                     "interpretation_top_candidate": directed_interp.top_candidates[0]["gene"],
                     "candidate_fingerprint": directed.candidate_fingerprint},
        "top20_overlap": sorted(set(classic_top) & {item.split(" (")[0] for item in directed_top}),
        "directed_display_columns": directed_display.columns.tolist(),
        "compare_display_columns": compare_display.columns.tolist(),
        "compare_ano2": _json_record(compare_display, "ANO2"),
        "compare_rflnb": _json_record(compare_display, "RFLNB"),
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
