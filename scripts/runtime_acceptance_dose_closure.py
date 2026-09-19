"""Real active dose-response closure with ordinary Classic state isolation."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import biology_logic, config, main

TARGET = "ENSP00000003084"
OUT = ROOT / "outputs" / "runtime_acceptance_closure" / "dose_response"


def _weights_digest(graph) -> str:
    return hashlib.sha256(np.asarray(graph.es["weight"], dtype=np.float64).tobytes()).hexdigest()


def _ordinary(graph, scores, label: str) -> pd.DataFrame:
    biology_logic.run_infection_simulation(
        graph, scores, OUT / label, spesifik_hedefler=[TARGET], block_weight_fraction=0.001, damping=0.85,
        exploratory_top_n=20, null_iterations=0, compute_structural_metrics=False, tissue="Lung",
        edge_evidence_weighting_mode="legacy_edge_modifiers",
    )
    return biology_logic.latest_signed_redistribution().copy(deep=True)


def _max_deltas(before: pd.DataFrame, after: pd.DataFrame) -> dict[str, float]:
    columns = ["PageRank_Baseline", "PageRank_Perturbed", "Delta_PageRank_Pct"]
    left = before.set_index("gene").sort_index()
    right = after.set_index("gene").sort_index()
    if not left.index.equals(right.index):
        raise AssertionError("ordinary Classic universe changed after dose-response")
    return {column: float(np.nanmax(np.abs(pd.to_numeric(left[column], errors="coerce") - pd.to_numeric(right[column], errors="coerce")))) for column in columns}


def main_run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    graph, scores = main.arayuz_icin_motoru_hazirla(forced_genes=[TARGET], hedef_doku="Lung", bc_sample_sources=4)
    initial_digest = _weights_digest(graph)
    ordinary_before = _ordinary(graph, scores, "ordinary_before")
    before_digest = _weights_digest(graph)
    dose_started = time.perf_counter()
    dose = biology_logic.run_pharmacological_dose_response(
        graph, scores, OUT, spesifik_hedefler=[TARGET], dozlar=list(config.DOZLAR),
        hill_n=config.HILL_COEFFICIENT,
    )
    dose_seconds = time.perf_counter() - dose_started
    after_dose_digest = _weights_digest(graph)
    ordinary_after = _ordinary(graph, scores, "ordinary_after")
    final_digest = _weights_digest(graph)
    required = ["Survival_Fraction", "Inhibition_Pct", "Hill_n", "RWR_Baseline", "RWR_Signal", "Signal_Kayip_Pct", "Max_Congestion_Score"]
    mandatory_finite = {column: bool(np.isfinite(pd.to_numeric(dose[column], errors="coerce")).all()) for column in required}
    values = dose["Survival_Fraction"].astype(float).tolist()
    result = {
        "target": TARGET, "tissue": "Lung", "configured_doses": list(config.DOZLAR), "configured_hill_n": config.HILL_COEFFICIENT,
        "seconds_total": time.perf_counter() - started, "seconds_dose": dose_seconds,
        "dose_rows": len(dose), "dose_values": values, "dose_values_unique": len(values) == len(set(values)),
        "dose_order_matches_config_descending": values == sorted(config.DOZLAR, reverse=True),
        "hill_values": dose["Hill_n"].astype(int).tolist(), "required_columns": required,
        "required_columns_present": all(column in dose for column in required), "mandatory_finite": mandatory_finite,
        "output_csv_exists": (OUT / "farmakolojik_doz_yanit_raporu.csv").exists(),
        "output_csv_rows": len(pd.read_csv(OUT / "farmakolojik_doz_yanit_raporu.csv")),
        "graph_digests": {"initial": initial_digest, "before_dose": before_digest, "after_dose": after_dose_digest, "final": final_digest},
        "ordinary_classic_max_deltas": _max_deltas(ordinary_before, ordinary_after),
        "ordinary_universe_rows": len(ordinary_before),
    }
    (OUT / "validation.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main_run()
