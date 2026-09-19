"""Generate a staging-only HPA and Classic CORUM provenance report."""

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

from src.product.legacy_provenance import (  # noqa: E402
    audit_hpa_mixed_artifact,
    inspect_hpa_git_evidence,
    reconstruct_classic_corum,
)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--project-root", type=Path, default=ROOT)
    parser.add_argument("--output-root", type=Path)
    args = parser.parse_args()
    root = args.project_root.resolve()
    output = (args.output_root or root / "data/derived/staging/legacy-provenance").resolve()
    hpa_builds = sorted((root / "data/derived/staging/hpa_human_tissue").glob("*/tissue_expression.sqlite"))
    if not hpa_builds:
        raise SystemExit("HPA clean shadow is required; run scripts/data_build.py --dataset hpa")
    hpa = audit_hpa_mixed_artifact(root / "data/raw/hinterland_core.db", hpa_builds[-1])
    hpa["source_evidence"] = inspect_hpa_git_evidence(root)
    raw = pd.read_csv(root / "data/raw/rna_tissue_hpa.tsv", sep="\t")
    with (root / "data/processed/hpa_pivot.pkl").open("rb") as stream:
        pivot = pickle.load(stream)
    exact = sum(
        pivot.get(str(row.Gene), {}).get(str(row.Tissue)) == float(row.TPM)
        for row in raw.itertuples(index=False)
    )
    hpa["tpm_source_chain"] = {
        "raw_rows": len(raw), "pivot_rows": sum(len(value) for value in pivot.values()),
        "raw_to_pivot_exact": exact, "semantic_exact": exact == len(raw),
    }
    corum = reconstruct_classic_corum(
        root / "data/raw/coreComplexes.txt",
        root / "data/processed/complex_members.pkl",
        output,
        root / "data/processed/ensp_with_symbols.csv",
        root / "data/processed/mapping.parquet",
    )
    output.mkdir(parents=True, exist_ok=True)
    aggregate = {
        "status": "PARTIAL",
        "hpa": hpa,
        "corum": corum,
        "production_promotion": False,
    }
    report_path = output / "legacy-provenance-report.json"
    report_path.write_text(json.dumps(aggregate, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps({"report": str(report_path), "hpa": hpa["classification"], "corum": corum["classification"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
