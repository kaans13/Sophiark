"""Bounded CFTR standalone Classic/Directed and Unified parity closure."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.unified.service import _default_engine_runs, run_unified_analysis

TARGET = "ENSP00000003084"
TISSUE = "Lung"
OUT = ROOT / "outputs" / "runtime_acceptance_closure"


def _deltas(left: pd.DataFrame, right: pd.DataFrame, key: str, columns: list[str]) -> dict[str, float]:
    a = left.set_index(key).sort_index()
    b = right.set_index(key).sort_index()
    if not a.index.equals(b.index):
        raise AssertionError(f"universe mismatch for {key}")
    return {column: float(np.nanmax(np.abs(pd.to_numeric(a[column], errors="coerce") - pd.to_numeric(b[column], errors="coerce")))) for column in columns}


def main() -> None:
    started = time.perf_counter()
    classic, directed, metadata, symbols, errors = _default_engine_runs(
        target=TARGET, tissue=TISSUE, block_strength=.001, damping=.85, bc_sample_sources=4,
    )
    standalone_seconds = time.perf_counter() - started
    old = ROOT / "outputs" / "runtime_acceptance" / "unified_panel" / "CFTR"
    old_classic = pd.read_csv(old / "unified_classic_raw.csv")
    old_directed = pd.read_csv(old / "unified_directed_raw.csv")
    classic_columns = ["PageRank_Baseline", "PageRank_Perturbed", "Delta_PageRank_Pct"]
    directed_columns = ["Directed_PageRank_Baseline", "Directed_PageRank_Perturbed", "Directed_Redistribution_Pct", "Directed_BC"]
    pre_fix = {
        "classic": _deltas(classic, old_classic, "gene", classic_columns),
        "directed": _deltas(directed, old_directed, "entity", directed_columns),
    }
    report = run_unified_analysis(
        target=TARGET, tissue=TISSUE, block_strength=.001, damping=.85, bc_sample_sources=4, use_cache=False,
        engine_runner=lambda **_kwargs: (classic.copy(deep=True), directed.copy(deep=True), metadata, symbols, errors),
    )
    unified = {
        "classic": _deltas(classic, report.classic_raw, "gene", classic_columns),
        "directed": _deltas(directed, report.directed_raw, "entity", directed_columns),
    }
    result = {
        "target": TARGET, "tissue": TISSUE, "standalone_seconds": standalone_seconds,
        "classic_rows": len(classic), "directed_rows": len(directed), "directed_status": metadata.get("status"), "errors": errors,
        "pre_fix_artifact_max_deltas": pre_fix, "unified_standalone_max_deltas": unified,
        "unified_status": report.status.value, "unified_classic_status": report.classic_status, "unified_directed_status": report.directed_status,
        "evidence_context": report.context_availability.get("evidence_provenance"), "corum_context": report.context_availability.get("corum_complex_context"),
    }
    (OUT / "bounded_parity.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
