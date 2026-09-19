"""Inject reconstructed CORUM in memory and compare it with frozen production data."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import pickle
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import biology_logic  # noqa: E402
from src.unified.service import _default_engine_runs  # noqa: E402


TARGETS = {
    "CFTR": "ENSP00000003084",
    "TGFBR1": "ENSP00000447297",
    "ANO2": "ENSP00000348453",
}


def _compare(left: pd.DataFrame, right: pd.DataFrame) -> dict[str, object]:
    key = "gene" if "gene" in left.columns else "entity_id"
    a = left.sort_values(key).reset_index(drop=True)
    b = right.sort_values(key).reset_index(drop=True)
    numeric = sorted(set(a.select_dtypes(include="number").columns) & set(b.select_dtypes(include="number").columns))
    deltas = {
        column: float((a[column].astype(float) - b[column].astype(float)).abs().max())
        for column in numeric
    }
    return {
        "rows": len(a),
        "keys_equal": a[key].astype(str).tolist() == b[key].astype(str).tolist(),
        "max_numeric_delta": max(deltas.values(), default=0.0),
        "numeric_deltas": deltas,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    artifacts = sorted((root / "data/derived/staging/legacy-provenance/corum").glob("*/classic_complex_members.pkl"))
    if not artifacts:
        raise SystemExit("run scripts/audit_legacy_provenance.py first")
    with (root / "data/processed/complex_members.pkl").open("rb") as stream:
        production = pickle.load(stream)
    with artifacts[-1].open("rb") as stream:
        reconstructed = pickle.load(stream)
    semantic_equal = set(production) == set(reconstructed) and all(
        set(production[key]) == set(reconstructed[key]) for key in production
    )
    if not semantic_equal:
        raise SystemExit("reconstructed CORUM is not semantically equivalent")
    runs = []
    original = biology_logic.COMPLEX_MEMBERS
    try:
        for symbol, target in TARGETS.items():
            biology_logic.COMPLEX_MEMBERS = production
            classic_a, directed_a, _, _, errors_a = _default_engine_runs(
                target=target, tissue="Lung", block_strength=0.1, damping=0.85, bc_sample_sources=4,
            )
            biology_logic.COMPLEX_MEMBERS = reconstructed
            classic_b, directed_b, _, _, errors_b = _default_engine_runs(
                target=target, tissue="Lung", block_strength=0.1, damping=0.85, bc_sample_sources=4,
            )
            runs.append({
                "symbol": symbol,
                "target": target,
                "classic": _compare(classic_a, classic_b),
                "directed": _compare(directed_a, directed_b),
                "errors": {"production": errors_a, "reconstructed": errors_b},
            })
    finally:
        biology_logic.COMPLEX_MEMBERS = original
    report = {
        "status": "COMPLETE",
        "semantic_corum_equivalence": semantic_equal,
        "production_mutated": False,
        "runs": runs,
    }
    output = (args.output or root / "data/derived/staging/legacy-provenance/corum-parity.json").resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "max_classic_delta": max(run["classic"]["max_numeric_delta"] for run in runs),
        "max_directed_delta": max(run["directed"]["max_numeric_delta"] for run in runs),
        "output": str(output),
    }))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
