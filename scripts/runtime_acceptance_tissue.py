"""Bounded real active tissue-differential acceptance run."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src import main
from src.services.tissue_differential import run_tissue_differential


def main_run() -> None:
    output = ROOT / "outputs" / "runtime_acceptance" / "tissue_smoke"
    started = time.perf_counter()
    result = run_tissue_differential(
        motor_module=main, tissues=["Lung", "Liver"], targets=["ENSP00000003084"], output_root=output,
        bc_sample_sources=4, block_strength=0.001, damping=0.85,
        ghost_filter=lambda frame, _graph: (frame, []), apply_druggability=lambda frame: frame,
    )
    payload = {
        "seconds": time.perf_counter() - started,
        "summary_rows": len(result.summary), "summary": result.summary.to_dict("records"),
        "critical_rows": len(result.critical_genes), "similarity": result.similarity.to_dict("records"),
        "notices": list(result.notices),
    }
    (output / "validation.json").write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2, default=str))


if __name__ == "__main__":
    main_run()
