"""One real bounded Evidence BETA execution using the active engine entry point."""
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

from src.evidence.config import EvidenceConfig
from src.evidence.engine import prepare_evidence_network, run_evidence_simulation


def main() -> None:
    mapping = pd.read_csv(ROOT / "data" / "processed" / "ensp_with_symbols.csv")
    target = str(mapping.loc[mapping["Symbol"].eq("CFTR"), "gene"].iloc[0])
    config = EvidenceConfig()
    started = time.perf_counter()
    graph, scores, metadata = prepare_evidence_network(config=config, forced_genes=[target], tissue="Lung", bc_sample_sources=4)
    before = hashlib.sha256(np.asarray(graph.es["weight"], dtype=np.float64).tobytes()).hexdigest()
    result = run_evidence_simulation(
        graph=graph.copy(), scores=scores, out_dir=ROOT / "outputs" / "runtime_acceptance" / "evidence_smoke",
        targets=[target], config=config, tissue="Lung", block_weight_fraction=.001, damping=.85,
        scientific_options={"exploratory_top_n": 20, "compute_structural_metrics": False, "null_iterations": 0},
    )
    payload = {
        "seconds": time.perf_counter() - started, "target": target,
        "graph": {"nodes": graph.vcount(), "edges": graph.ecount(), "base_unchanged": before == hashlib.sha256(np.asarray(graph.es["weight"], dtype=np.float64).tobytes()).hexdigest()},
        "config": config.to_dict(), "metadata": metadata,
        "status": "COMPLETE", "signed_rows": len(result.signed_response), "complex_rows": len(result.complexes),
        "comparison_rows": len(result.comparison), "provenance_embedded": "Evidence_Run_Provenance_JSON" in result.signed_response,
    }
    target_path = ROOT / "outputs" / "runtime_acceptance" / "evidence_validation.json"
    target_path.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2, default=str))


if __name__ == "__main__":
    main()
