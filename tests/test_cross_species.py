from __future__ import annotations

import hashlib
import json
import sqlite3
import tempfile
from pathlib import Path

import igraph as ig
import pandas as pd
import pytest

from src.cross_species import (
    ComparisonConfig,
    OrthologyClass,
    OrthologyIndex,
    SuitabilityStatus,
    check_mouse_suitability,
    compare_ortholog_responses,
    discover_functional_proxies,
    resolve_mouse_entity,
)
from src.graph_engine import _edge_attrs as human_edge_attrs
from src.mouse.graph_engine import _edge_attrs as mouse_edge_attrs
from src.mouse.data_loader import fetch_degree_summary
from src.services.directed_analysis import run_optional_directed_engine


HEADER = (
    "Gene stable ID\tGene name\tMouse gene stable ID\tMouse gene name\t"
    "Mouse homology type\tMouse orthology confidence [0 low, 1 high]\n"
)


def orthology_index(tmp_path: Path) -> OrthologyIndex:
    content = (
        HEADER
        + "ENSG1\tCFTR\tENSMUSG1\tCftr\tortholog_one2one\t1\n"
        + "ENSG2\tTP53\tENSMUSG2\tTrp53\tortholog_one2one\t1\n"
        + "ENSG3\tMULTI\tENSMUSG3\tMulti1\tortholog_one2many\t1\n"
        + "ENSG3\tMULTI\tENSMUSG4\tMulti2\tortholog_one2many\t1\n"
        + "ENSG4\tNONE\t\t\t\t\n"
    ).encode()
    tsv = tmp_path / "orthology.tsv"
    metadata = tmp_path / "orthology.json"
    tsv.write_bytes(content)
    metadata.write_text(json.dumps({
        "source": "Ensembl Compara test fixture", "human_taxon_id": 9606,
        "mouse_taxon_id": 10090, "sha256": hashlib.sha256(content).hexdigest(),
    }), encoding="utf-8")
    return OrthologyIndex.from_snapshot(tsv, metadata)


def test_authoritative_orthology_supports_one_to_one_one_to_many_and_no_ortholog(tmp_path):
    index = orthology_index(tmp_path)
    assert index.one_to_one("cftr").mouse_gene == "Cftr"
    assert {row.orthology_class for row in index.lookup("MULTI")} == {OrthologyClass.ONE_TO_MANY}
    assert {row.mouse_gene for row in index.lookup("MULTI")} == {"Multi1", "Multi2"}
    assert index.lookup("NONE")[0].orthology_class is OrthologyClass.NO_ORTHOLOG
    assert index.lookup("missing") == ()


def test_orthology_snapshot_is_checksum_and_species_safe(tmp_path):
    index = orthology_index(tmp_path)
    assert index.metadata["human_taxon_id"] == 9606
    tsv = tmp_path / "orthology.tsv"
    tsv.write_text(tsv.read_text() + "corrupt", encoding="utf-8")
    with pytest.raises(ValueError, match="CHECKSUM"):
        OrthologyIndex.from_snapshot(tsv, tmp_path / "orthology.json")


def test_suitability_distinguishes_available_graph_context_multiple_and_unresolved(tmp_path):
    index = orthology_index(tmp_path)
    mapping = {"cftr": "ENSMUSP1", "trp53": "ENSMUSP2"}
    available = check_mouse_suitability(
        "CFTR", index, mapping, mouse_graph_nodes={"ENSMUSP1"},
        context_nodes={"ENSMUSP1"}, context="lung",
    )
    assert available.status is SuitabilityStatus.DIRECT_ORTHOLOG_AVAILABLE
    missing_context = check_mouse_suitability(
        "CFTR", index, mapping, mouse_graph_nodes={"ENSMUSP1"}, context_nodes=set(),
    )
    assert missing_context.status is SuitabilityStatus.ORTHOLOG_EXISTS_CONTEXT_UNAVAILABLE
    assert check_mouse_suitability("MULTI", index, mapping).status is SuitabilityStatus.MULTIPLE_ORTHOLOGS
    assert check_mouse_suitability("UNKNOWN", index, mapping).status is SuitabilityStatus.UNRESOLVED


def test_mouse_entity_resolution_is_species_safe_and_uses_authoritative_orthology_alias(tmp_path):
    index = orthology_index(tmp_path)
    mapping = {"Cftr": "ENSMUSP1", "Trp53": "ENSMUSP2"}
    assert resolve_mouse_entity("CFTR", mapping, index).protein_id == "ENSMUSP1"
    tp53 = resolve_mouse_entity("TP53", mapping, index)
    assert tp53.official_symbol == "Trp53"
    assert tp53.mapping_type == "AUTHORITATIVE_ONE_TO_ONE_ORTHOLOG_ALIAS"
    assert resolve_mouse_entity("Trp53", mapping, index).protein_id == "ENSMUSP2"
    assert resolve_mouse_entity("ENSP00000269305", mapping, index).status == "SPECIES_MISMATCH"
    assert resolve_mouse_entity("MULTI", mapping, index).status == "AMBIGUOUS"


def test_functional_proxy_is_explanatory_and_never_automatically_selected():
    candidates = discover_functional_proxies(
        {"family": ["ABC"], "pathway": ["P1"], "go_mf": ["M1"], "go_bp": ["B1"]},
        {
            "CandidateA": {"family": ["ABC"], "pathway": ["P1"], "go_mf": ["M1"], "evidence_source": "GO fixture"},
            "CandidateB": {"go_bp": ["B1"], "evidence_source": "GO fixture"},
        },
    )
    assert candidates[0].mouse_gene == "CandidateA"
    assert candidates[0].same_family is True
    assert all(candidate.automatically_selected is False for candidate in candidates)


def test_cross_species_comparison_uses_one_to_one_rank_percentiles_and_direction_classes(tmp_path):
    index = orthology_index(tmp_path)
    human = pd.DataFrame({"Symbol": ["CFTR", "TP53", "MULTI"], "Delta_PageRank_Pct": [4.0, -2.0, 9.0]})
    mouse = pd.DataFrame({"Symbol": ["Cftr", "Trp53", "Multi1"], "Delta_PageRank_Pct": [3.0, 1.0, 8.0]})
    config = ComparisonConfig(500, 0.85, 0.001, "top_n", "v1")
    result = compare_ortholog_responses(human, mouse, index, human_config=config, mouse_config=config)
    assert result["metrics"]["comparable_one_to_one_count"] == 2
    assert result["metrics"]["concordant_positive"] == 1
    assert result["metrics"]["discordant"] == 1
    assert "rank_percentile_Human" in result["matched_responses"]
    assert "rank_percentile_Mouse" in result["matched_responses"]


def test_cross_species_comparison_refuses_silent_config_mismatch(tmp_path):
    index = orthology_index(tmp_path)
    frame = pd.DataFrame({"Symbol": ["CFTR"], "Delta_PageRank_Pct": [1.0]})
    with pytest.raises(ValueError, match="CONFIG_MISMATCH"):
        compare_ortholog_responses(
            frame, frame.rename(columns={"Symbol": "Symbol"}), index,
            human_config=ComparisonConfig(500, .85, .001, "top_n", "v1"),
            mouse_config=ComparisonConfig(500, .75, .001, "top_n", "v1"),
        )


def test_mouse_combined_score_base_weight_and_rescue_distance_match_human_contract():
    for score in (500, 777, 999):
        assert mouse_edge_attrs(score) == human_edge_attrs(score)
        assert mouse_edge_attrs(score)["weight"] == pytest.approx((score / 1000.0) ** 2)


def test_human_mouse_synthetic_pagerank_perturbation_and_selection_contracts_match(monkeypatch):
    from src import biology_logic as human
    from src.mouse import biology_logic as mouse

    def fixture(prefix: str):
        graph = ig.Graph(edges=[(0, 1), (1, 2), (1, 3), (2, 3)], directed=False)
        names = [f"{prefix}{i}" for i in range(4)]
        graph.vs["name"] = names
        graph.vs["lokalizasyon"] = ["Unknown"] * 4
        graph.es["weight"] = [1.0, .8, .6, .4]
        graph.es["distance"] = [1.0 / value for value in graph.es["weight"]]
        frame = pd.DataFrame({
            "gene": names, "Hinterland_Skoru": [80, 70, 60, 50],
            "Lokalizasyon": ["Unknown"] * 4, "BC_Skoru": [0.0] * 4,
        })
        return graph, frame, names[1]

    for module in (human, mouse):
        monkeypatch.setattr(module, "apply_biological_bonus", lambda *args, **kwargs: None)
        monkeypatch.setattr(module, "apply_functional_penalty", lambda *args, **kwargs: None)
    outputs = []
    with tempfile.TemporaryDirectory() as directory:
        for module, prefix in ((human, "ENSP"), (mouse, "ENSMUSP")):
            graph, frame, target = fixture(prefix)
            module.run_infection_simulation(
                graph, frame, Path(directory) / prefix, spesifik_hedefler=[target],
                redistribution_mode="top_n", exploratory_top_n=4, null_iterations=0,
                compute_structural_metrics=False, damping=.85, block_weight_fraction=.001,
            )
            signed = module.latest_signed_redistribution().sort_values("gene").reset_index(drop=True)
            outputs.append(signed[["Delta_PageRank", "Delta_PageRank_Pct", "Response_Direction"]])
    pd.testing.assert_frame_equal(outputs[0], outputs[1], check_exact=False, rtol=1e-12, atol=1e-14)


def test_human_mouse_hinterland_and_bc_contracts_match_on_same_graph():
    from src.metrics import build_hinterland_score as human_hinterland, detect_customs_gates as human_bc
    from src.mouse.metrics import build_hinterland_score as mouse_hinterland, detect_customs_gates as mouse_bc

    graph = ig.Graph(edges=[(0, 1), (1, 2), (1, 3)], directed=False)
    graph.vs["name"] = ["A", "B", "C", "D"]
    graph.vs["lokalizasyon"] = ["Nucleus", "Nucleus", "Cytoplasm", "Cytoplasm"]
    graph.es["weight"] = [1.0, .8, .6]
    graph.es["distance"] = [1.0, 1.25, 1.0 / .6]
    base = pd.DataFrame({"gene": graph.vs["name"], "k_i": [1, 3, 1, 1], "max_cs": [800] * 4, "Lokalizasyon": graph.vs["lokalizasyon"]})
    pr = dict(zip(graph.vs["name"], graph.pagerank(weights="weight", damping=.85)))
    human_score = human_hinterland(base, pr)
    mouse_score = mouse_hinterland(base, pr)
    assert human_score["Hinterland_Skoru"].tolist() == mouse_score["Hinterland_Skoru"].tolist()
    human_result = human_bc(graph, human_score, bc_sample_sources=None)
    mouse_result = mouse_bc(graph, mouse_score, bc_sample_sources=None)
    assert human_result["BC_Skoru"].tolist() == mouse_result["BC_Skoru"].tolist()


def test_mouse_degree_summary_deduplicates_symmetric_string_rows(tmp_path):
    db = tmp_path / "mouse.sqlite"
    con = sqlite3.connect(db)
    con.execute("CREATE TABLE mouse_interactions(protein1 TEXT, protein2 TEXT, combined_score INTEGER)")
    con.executemany("INSERT INTO mouse_interactions VALUES (?,?,?)", [
        ("ENSMUSP1", "ENSMUSP2", 800), ("ENSMUSP2", "ENSMUSP1", 800),
        ("ENSMUSP1", "ENSMUSP3", 700), ("ENSMUSP3", "ENSMUSP1", 700),
    ])
    con.commit(); con.close()
    result = fetch_degree_summary(str(db), 500).set_index("gene")
    assert result.loc["ENSMUSP1", "k_i"] == 2
    assert result.loc["ENSMUSP2", "k_i"] == 1


def test_mouse_directed_is_fail_closed_below_mapping_acceptance_threshold(tmp_path):
    bundle = run_optional_directed_engine(
        graph=None, targets=("ENSMUSP1",), tissue="lung", classic_report=None,
        project_root=tmp_path, mode="Directed", block_weight_fraction=.001,
        damping=.85, taxon_id=10090,
    )
    assert bundle.calculation.status.value == "UNAVAILABLE"
    assert "mapping" in bundle.calculation.error.lower()


def test_mouse_ui_disables_research_interpretation_and_directed_production():
    source = Path("app.py").read_text(encoding="utf-8")
    assert "MOUSE_RESEARCH_EXPLORER_DISABLED" in source
    assert "MOUSE_BIOLOGICAL_INTERPRETATION_DISABLED" in source
    assert '_engine_options = ["Classic"] if IS_MOUSE' in source
    assert "View Functional Proxy Candidates" in source
    assert "No candidate is automatically used as a perturbation target" in source
