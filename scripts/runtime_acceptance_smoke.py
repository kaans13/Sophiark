"""Bounded real runtime checks for non-Unified active simulation services."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src import main
from src.core.analysis_runtime import RuntimeContext, run_compensation_analysis, run_target_stress_search
from src.services.analysis_service import execute_simulation, prepare_network


OUT = Path("outputs/runtime_acceptance/secondary_smoke")
TARGETS = ["ENSP00000003084", "ENSP00000447297"]


def _identity(frame, _graph=None):
    return frame, []


def _graph_digest(graph) -> str:
    return hashlib.sha256(np.asarray(graph.es["weight"], dtype=np.float64).tobytes()).hexdigest()


def _run(*, graph, scores, targets, mode, null_iterations=0):
    started = time.perf_counter()
    report, ghosts = execute_simulation(
        motor_module=main, graph=graph, scores=scores, output_dir=OUT / mode,
        targets=targets, localization="SNIPER_MODU", block_strength=0.001, damping=0.85,
        apply_druggability=lambda frame: frame, ghost_filter=_identity,
        scientific_options={
            "redistribution_mode": "null_fdr" if null_iterations else "top_n",
            "exploratory_top_n": 20,
            "null_iterations": null_iterations,
            "random_seed": 42,
            "compute_structural_metrics": False,
            "tissue": "Lung",
            "tissue_normalization_mode": "within_tissue",
            "bc_mode": "fast_approximate",
            "bc_sample_size": 4,
            "edge_evidence_weighting_mode": "structural_only",
        },
    )
    return report, ghosts, time.perf_counter() - started


def main_run() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    started = time.perf_counter()
    prepared = prepare_network(main, tuple(TARGETS), "Lung", 4, "within_tissue")
    graph, scores = prepared.graph, prepared.scores
    graph_before = _graph_digest(graph)
    result: dict[str, object] = {
        "acceptance_configuration": {"tissue": "Lung", "bc_sample_sources": 4, "top_n": 20, "null_iterations": 3},
        "prepare_seconds": time.perf_counter() - started,
        "network": {"nodes": graph.vcount(), "edges": graph.ecount(), "base_weight_digest": graph_before},
    }

    multi, ghosts, elapsed = _run(graph=graph, scores=scores, targets=TARGETS, mode="multi_target")
    result["multi_target"] = {
        "seconds": elapsed, "rows": len(multi), "ghosts": len(ghosts),
        "requested_targets": TARGETS,
        "resolved_targets": sorted(set(multi.get("Hasar_Hedefleri", []).astype(str).iloc[0].split(" | "))) if not multi.empty else [],
        "primary_rows": int((multi.get("Hasar_Tipi", "") == main.HASAR_BIRINCIL).sum()) if not multi.empty else 0,
        "combined_damage_values": sorted(set(map(str, multi.get("Kombine_Hasar", [])))),
        "metadata_modes": sorted(set(map(str, multi.get("Redistribution_Selection_Mode", [])))),
    }

    null, ghosts, elapsed = _run(graph=graph, scores=scores, targets=[TARGETS[0]], mode="null_smoke", null_iterations=3)
    result["null_smoke"] = {
        "seconds": elapsed, "rows": len(null), "ghosts": len(ghosts),
        "configured_iterations": 3,
        "reported_iterations": sorted(set(map(int, null.get("Null_Model_Iterations", [])))) if not null.empty else [],
        "selection_modes": sorted(set(map(str, null.get("Redistribution_Selection_Mode", [])))),
        "base_graph_unchanged": graph_before == _graph_digest(graph),
    }

    context = RuntimeContext(motor_module=main, output_dir=OUT, is_mouse=False, tissue="Lung")
    mapping = dict(zip(scores.get("gene", []), scores.get("Symbol", [])))
    stress = run_target_stress_search(
        context=context, graph=graph, scores=scores, target_gene=TARGETS[0], regulators={}, gene_to_symbol=mapping,
        mode="Network-mediated", limit=1, block_strength=0.001, damping=0.85, ghost_filter=_identity,
    )
    result["target_stress_search"] = {"rows": len(stress), "columns": list(stress.columns), "ids": stress.get("gen", []).astype(str).tolist() if not stress.empty else []}
    compensation = run_compensation_analysis(
        context=context, graph=graph, scores=scores, perturbed_gene=TARGETS[0], block_strength=0.001,
        damping=0.85, ghost_filter=_identity,
    )
    result["compensation"] = {"rows": len(compensation), "columns": list(compensation.columns), "ids": compensation.get("Protein ID", []).astype(str).tolist() if not compensation.empty else []}

    invalid, ghosts, elapsed = _run(graph=graph, scores=scores, targets=["ENSP_INVALID_RUNTIME_ACCEPTANCE"], mode="invalid_target")
    result["invalid_target"] = {"seconds": elapsed, "rows": len(invalid), "ghosts": len(ghosts), "is_empty": bool(invalid.empty)}
    result["base_graph_unchanged_after_all"] = graph_before == _graph_digest(graph)
    result["total_seconds"] = time.perf_counter() - started
    (OUT / "secondary_smoke.json").write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main_run()
