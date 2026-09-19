"""File-backed integration validation for the scientific hardening sprint."""

from __future__ import annotations

import json
from pathlib import Path
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))


OUTPUT = Path("outputs") / "scientific_validation" / "hardening_integration.json"
OUTPUT.parent.mkdir(parents=True, exist_ok=True)
state: dict[str, object] = {"status": "started", "steps": {}}


def save() -> None:
    OUTPUT.write_text(json.dumps(state, ensure_ascii=False, indent=2, default=str), encoding="utf-8")


def metric(report, name):
    return float(report[name].iloc[0])


try:
    from src import main as human
    from src.biology_logic import validate_against_known_biology
    from src.metrics import validate_jaccard_stability
    from src.scientific.validation import pagerank_robustness

    state["steps"]["human_prepare"] = "running"
    save()
    human_graph, human_scores = human.arayuz_icin_motoru_hazirla(
        forced_genes=["ENSP00000003084"], hedef_doku=None, bc_sample_sources=300,
        tissue_normalization_mode="within_tissue",
    )
    state["steps"]["human_prepare"] = {
        "nodes": human_graph.vcount(), "edges": human_graph.ecount(), "rows": len(human_scores)
    }
    save()

    report = human.run_infection_simulation(
        human_graph.copy(), human_scores,
        OUTPUT.parent / "human_cftr", spesifik_hedefler=["ENSP00000003084"],
        block_weight_fraction=.001, damping=.85,
        redistribution_mode="null_fdr", null_iterations=19,
        exploratory_top_n=100, efficiency_sample_sources=64,
        local_efficiency_sample_nodes=32, tissue="None",
        bc_mode="fast_approximate", bc_sample_size=300,
    )
    directions = report.get("Response_Direction")
    state["steps"]["human_forward"] = {
        "rows": len(report),
        "influence_gain": int((directions == "Influence Gain").sum()),
        "influence_loss": int((directions == "Influence Loss").sum()),
        "significant_fdr": int(report.get("significant_redistribution", False).fillna(False).sum()),
        "systemic_network_shift_pct": metric(report, "Systemic_Network_Shift_Pct"),
        "global_efficiency_baseline": metric(report, "Global_Efficiency_Baseline"),
        "global_efficiency_perturbed": metric(report, "Global_Efficiency_Perturbed"),
        "global_efficiency_change_pct": metric(report, "Global_Efficiency_Change_Pct"),
        "mean_local_efficiency_change_pct": metric(report, "Mean_Local_Efficiency_Change_Pct"),
        "target_local_efficiency_change_pct": metric(report, "Target_Local_Efficiency_Change_Pct"),
        "deprecated_alias_preserved": metric(report, "Efficiency_Kayip_Pct") == metric(report, "Systemic_Network_Shift_Pct"),
    }
    save()

    robustness, robustness_meta = pagerank_robustness(human_graph, top_k=100)
    state["steps"]["pagerank_robustness"] = {
        "metadata": robustness_meta,
        "rows": robustness.to_dict("records"),
    }
    state["steps"]["sentinel"] = validate_against_known_biology(human_scores, human_graph)
    state["steps"]["jaccard_stability"] = validate_jaccard_stability(human_graph, human_scores)
    save()

    from src.graph_engine import build_high_conf_subgraph
    from src.config import DB_PATH, HIGH_CONF_THRESHOLD
    liver_graph = build_high_conf_subgraph(
        DB_PATH, HIGH_CONF_THRESHOLD, hedef_doku="Liver",
        tissue_normalization_mode="cross_tissue_comparable",
    )
    state["steps"]["cross_tissue_graph"] = {
        "nodes": liver_graph.vcount(), "edges": liver_graph.ecount(),
        "normalization": liver_graph["tissue_normalization"],
    }
    save()

    from src.mouse import main as mouse
    mouse_graph, mouse_scores = mouse.arayuz_icin_motoru_hazirla(
        forced_genes=["ENSMUSP00000049228"], hedef_doku=None, bc_sample_sources=300,
        tissue_normalization_mode="within_tissue",
    )
    mouse_report = mouse.run_infection_simulation(
        mouse_graph.copy(), mouse_scores,
        Path("outputs_mouse") / "scientific_validation" / "mouse_cftr",
        spesifik_hedefler=["ENSMUSP00000049228"], block_weight_fraction=.001,
        damping=.85, redistribution_mode="top_n", null_iterations=0,
        exploratory_top_n=100, efficiency_sample_sources=32,
        local_efficiency_sample_nodes=16, tissue="None",
        bc_mode="fast_approximate", bc_sample_size=300,
    )
    mouse_directions = mouse_report.get("Response_Direction")
    state["steps"]["mouse_forward"] = {
        "nodes": mouse_graph.vcount(), "edges": mouse_graph.ecount(),
        "rows": len(mouse_report),
        "target_format_ok": str(mouse_report.iloc[0]["gene"]).startswith("ENSMUSP"),
        "influence_gain": int((mouse_directions == "Influence Gain").sum()),
        "influence_loss": int((mouse_directions == "Influence Loss").sum()),
        "global_efficiency_change_pct": metric(mouse_report, "Global_Efficiency_Change_Pct"),
    }
    save()

    state["status"] = "complete"
except Exception as error:
    state["status"] = "failed"
    state["error"] = {"type": type(error).__name__, "message": str(error), "traceback": traceback.format_exc()}
finally:
    save()
