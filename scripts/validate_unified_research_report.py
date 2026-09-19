"""Run a small real-data Unified Research Report validation panel."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.unified import export_unified_report, run_unified_analysis
from src.unified.service import _default_engine_runs


DEFAULT_TARGETS = {
    "TGFBR1": "ENSP00000447297",
    "BAD": "ENSP00000378040",
    "ANO2": "ENSP00000348453",
    "EGFR": "ENSP00000275493",
    "CFTR": "ENSP00000003084",
}


def _max_delta(left, right, entity: str, columns: tuple[str, ...]) -> dict[str, float]:
    if left.empty or right.empty:
        return {column: float("nan") for column in columns}
    a = left.set_index(entity).sort_index(); b = right.set_index(entity).sort_index()
    if not a.index.equals(b.index):
        raise AssertionError("Standalone and unified response universes differ")
    return {
        column: float(np.nanmax(np.abs(a[column].astype(float).to_numpy() - b[column].astype(float).to_numpy())))
        for column in columns
    }


def run(*, targets: dict[str, str], tissue: str, output: Path, verify_targets: set[str]) -> dict[str, object]:
    rows = []
    for symbol, target in targets.items():
        report = run_unified_analysis(target=target, tissue=tissue, bc_sample_sources=4, use_cache=False)
        item: dict[str, object] = {
            "symbol": symbol, "target": target, "status": report.status.value,
            "classic_status": report.classic_status, "directed_status": report.directed_status,
            "classic_rows": len(report.classic_raw), "directed_rows": len(report.directed_raw),
            "candidate_rows": len(report.candidates), "agreement": dict(report.agreement),
            "context": dict(report.context_availability), "errors": dict(report.errors),
        }
        if symbol in verify_targets:
            standalone_classic, standalone_directed, _, _, _ = _default_engine_runs(
                target=target, tissue=tissue, block_strength=.001, damping=.85, bc_sample_sources=4,
            )
            item["classic_max_deltas"] = _max_delta(
                standalone_classic, report.classic_raw, "gene",
                ("PageRank_Baseline", "PageRank_Perturbed", "Delta_PageRank_Pct"),
            )
            item["directed_max_deltas"] = _max_delta(
                standalone_directed, report.directed_raw, "entity",
                ("Directed_PageRank_Baseline", "Directed_PageRank_Perturbed", "Directed_Redistribution_Pct", "Directed_BC"),
            )
        export_unified_report(report, output / symbol)
        rows.append(item)
        print(json.dumps(item, default=str), flush=True)
    result = {"status": "COMPLETE", "tissue": tissue, "runs": rows}
    output.mkdir(parents=True, exist_ok=True)
    (output / "unified_validation_summary.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--tissue", default="Lung")
    parser.add_argument("--output", type=Path, default=Path("outputs/evidence_validation/unified_report_panel"))
    parser.add_argument("--targets", nargs="*", choices=tuple(DEFAULT_TARGETS), default=list(DEFAULT_TARGETS))
    parser.add_argument("--verify-targets", nargs="*", choices=tuple(DEFAULT_TARGETS), default=("TGFBR1", "ANO2"))
    args = parser.parse_args()
    run(targets={symbol: DEFAULT_TARGETS[symbol] for symbol in args.targets}, tissue=args.tissue, output=args.output, verify_targets=set(args.verify_targets))
