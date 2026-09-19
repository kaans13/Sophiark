from __future__ import annotations

from dataclasses import replace
from pathlib import Path
import tempfile

import igraph as ig
import numpy as np
import pandas as pd
from unittest.mock import patch

from src.directed import (
    build_hybrid_directed_graph, compare_classic_and_directed,
    compute_directed_pagerank, directed_graph_cache_key,
    effective_graph_fingerprint, graph_weight_parity,
    prepare_effective_directed_source_graph, run_directed_calculation,
    suppress_target_edges,
)
from src.signaling.models import (
    DatasetQA, LayerMembership, MappingReport, SignalingDataset,
    SignalingEdge, SignalingNode, SignStatus,
)
from src.services.directed_analysis import directed_export_tables, run_optional_directed_engine
from scripts.audit_directed_perturbation import run_audit


def _string_graph(edges, weights=None):
    names = sorted({item for edge in edges for item in edge})
    index = {name: i for i, name in enumerate(names)}
    graph = ig.Graph(n=len(names), edges=[(index[a], index[b]) for a, b in edges], directed=False)
    graph.vs["name"] = names
    graph.es["weight"] = list(weights or [1.0] * len(edges))
    graph.es["distance"] = [1.0 / value for value in graph.es["weight"]]
    graph.es["rescue_bridge"] = [False] * graph.ecount()
    return graph


def _dataset(directions, fingerprint="fixture"):
    names = sorted({item for edge in directions for item in edge})
    nodes = tuple(SignalingNode(
        node_id=name, source_identifier=name, symbol=name,
        layer_membership=LayerMembership.LAYER_1_AND_2,
        layer1_id=name, mapping_status="EXACT",
    ) for name in names)
    edges = tuple(SignalingEdge(
        source=source, target=target, source_identifier=source,
        source_symbol=source, target_identifier=target, target_symbol=target,
        directed=True, stimulation=True, inhibition=False, canonical_sign=1,
        sign_status=SignStatus.ACTIVATION, interaction="fixture",
        resources=("fixture",), references=(f"{source}>{target}",),
    ) for source, target in directions)
    qa = DatasetQA(len(edges), len(edges), len(nodes), len(edges), len(edges), 0, 0, 0, 0, len(edges), 0)
    mapping = MappingReport(
        len(nodes), len(nodes), 0, 0, 0, len(edges), len(edges), 0, 0, 0,
        1.0, 1.0, len(nodes), len(nodes), len(edges), len(edges), 1.0, 1.0, "PREFERRED",
    )
    return SignalingDataset(nodes, edges, qa, mapping, {"dataset_version": "fixture-v1"}, fingerprint)


def test_single_omnipath_direction_removes_classic_reverse_arc():
    hybrid = build_hybrid_directed_graph(_string_graph([("A", "B")]), _dataset([("A", "B")]))
    graph = hybrid.graph
    assert graph.is_directed()
    assert graph.get_eid("A", "B", directed=True, error=False) != -1
    assert graph.get_eid("B", "A", directed=True, error=False) == -1
    assert hybrid.provenance.direction_resolved_edges == 1
    assert hybrid.provenance.fallback_string_edges == 0


def test_fallback_is_bidirectional_with_equal_original_string_weight():
    graph = build_hybrid_directed_graph(
        _string_graph([("A", "B")], [0.37]), _dataset([]),
    ).graph
    assert graph.ecount() == 2
    assert graph.es["weight"] == [0.37, 0.37]
    assert all(graph.es["fallback_bidirectional"])


def test_omnipath_only_edge_is_never_added_to_calculation_graph():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B")]), _dataset([("A", "C"), ("A", "B")]),
    )
    assert "C" not in hybrid.graph.vs["name"]
    assert hybrid.graph.ecount() == 1


def test_symmetric_directed_pagerank_matches_classic_undirected():
    classic = _string_graph([("A", "B"), ("B", "C"), ("C", "A")], [0.5, 0.8, 1.0])
    hybrid = build_hybrid_directed_graph(classic, _dataset([]))
    classic_values = np.asarray(classic.pagerank(damping=0.85, weights="weight", directed=False))
    directed = compute_directed_pagerank(hybrid.graph, damping=0.85)
    directed_values = np.asarray([directed[name] for name in classic.vs["name"]])
    np.testing.assert_allclose(classic_values, directed_values, rtol=1e-10, atol=1e-12)


def test_effective_modifier_context_and_full_symmetric_equivalence_gate():
    import src.biology_logic as classic_logic

    graph = _string_graph(
        [("T", "A"), ("T", "B"), ("A", "C"), ("B", "C"), ("A", "B")],
        [1.0, 0.8, 0.6, 0.4, 0.3],
    )
    scores = pd.DataFrame({
        "gene": graph.vs["name"],
        "Hinterland_Skoru": [40.0, 30.0, 20.0, 10.0],
        "Lokalizasyon": ["fixture"] * 4,
    })
    symbols = {"T": "TARGET", "A": "ALPHA", "B": "BETA", "C": "GAMMA"}
    original_fingerprint = effective_graph_fingerprint(graph)
    with (
        patch.object(classic_logic, "BONUS_ENABLED", True),
        patch.object(classic_logic, "TF_TARGETS", {"TARGET": {"ALPHA"}}),
        patch.object(classic_logic, "COMPLEX_MEMBERS", {
            "TARGET": ["complex-1"], "BETA": ["complex-1"],
        }),
        patch.object(classic_logic, "_ensp_to_symbol", side_effect=lambda value: symbols.get(value)),
        patch.object(classic_logic, "_save_symbol_cache"),
        patch.object(
            classic_logic.fetch_localizations_mygene,
            "_mf_cache",
            {"T": "kinase", "A": "kinase", "B": "other"},
            create=True,
        ),
    ):
        source = prepare_effective_directed_source_graph(graph, ["T"])
        assert source.original_fingerprint_before == source.original_fingerprint_after
        assert effective_graph_fingerprint(graph) == original_fingerprint
        parity = graph_weight_parity(source.graph, source.graph.copy())
        assert parity["node_mismatch_count"] == 0
        assert parity["pair_mismatch_count"] == 0
        assert parity["weight_mismatch_count"] == 0

        pair_weights = {
            frozenset((source.graph.vs[e.source]["name"], source.graph.vs[e.target]["name"])): e["weight"]
            for e in source.graph.es
        }
        assert pair_weights[frozenset(("T", "A"))] == 0.84
        assert pair_weights[frozenset(("T", "B"))] == 0.96

        with tempfile.TemporaryDirectory() as output:
            _, classic_metrics = classic_logic.run_infection_simulation(
                graph.copy(), scores, Path(output), spesifik_hedefler=["T"],
                exploratory_top_n=20, compute_structural_metrics=False,
                return_comparison_metrics=True, comparison_genes=["T"],
            )
        classic_signed = classic_logic.latest_signed_redistribution().set_index("gene")

    symmetric = build_hybrid_directed_graph(
        source.graph, _dataset([("T", "A")]), force_symmetric_fallback=True,
        effective_graph_fingerprint=source.effective_fingerprint,
    )
    result = run_directed_calculation(
        symmetric, targets=["T"], block_weight_fraction=0.001,
        damping=0.85, bc_sample_sources=0,
    )
    directed = result.report.set_index("entity")
    np.testing.assert_allclose(
        classic_signed.loc[directed.index, "PageRank_Baseline"],
        directed["Directed_PageRank_Baseline"], rtol=1e-10, atol=1e-12,
    )
    np.testing.assert_allclose(
        classic_signed.loc[directed.index, "PageRank_Perturbed"],
        directed["Directed_PageRank_Perturbed"], rtol=1e-10, atol=1e-12,
    )
    np.testing.assert_allclose(
        classic_signed.loc[directed.index, "Delta_PageRank_Pct"],
        directed["Directed_Redistribution_Pct"], rtol=1e-8, atol=1e-10,
    )
    classic_target = classic_metrics["observed"]["T"]
    directed_target = result.metadata["target_pagerank"]["T"]
    np.testing.assert_allclose(
        [classic_target["pagerank_before"], classic_target["pagerank_after"]],
        [directed_target["baseline"], directed_target["perturbed"]],
        rtol=1e-10, atol=1e-12,
    )


def test_one_way_chain_has_no_implicit_reverse_and_target_suppression_uses_directed_arcs():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C")]),
        _dataset([("A", "B"), ("B", "C")]),
    )
    graph = hybrid.graph
    assert graph.get_eid("B", "A", directed=True, error=False) == -1
    assert graph.get_eid("C", "B", directed=True, error=False) == -1
    perturbed_a, edges_a = suppress_target_edges(graph, ["A"], fraction=0.001)
    perturbed_c, edges_c = suppress_target_edges(graph, ["C"], fraction=0.001)
    assert edges_a == (graph.get_eid("A", "B", directed=True),)
    assert edges_c == (graph.get_eid("B", "C", directed=True),)
    assert perturbed_a.es[edges_a[0]]["weight"] == 0.001
    assert perturbed_c.es[edges_c[0]]["weight"] == 0.001


def test_fork_runs_combined_directed_perturbation_and_redistribution():
    string = _string_graph([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")])
    hybrid = build_hybrid_directed_graph(
        string, _dataset([("A", "B"), ("A", "C"), ("B", "D"), ("C", "D")]),
    )
    result = run_directed_calculation(hybrid, targets=["B", "C"], bc_sample_sources=None)
    assert result.status.value == "AVAILABLE"
    assert result.metadata["targets"] == ("B", "C")
    assert result.metadata["suppressed_arc_count"] == 4
    assert set(["Directed_PageRank_Baseline", "Directed_PageRank_Perturbed", "Directed_Redistribution_Pct"]).issubset(result.report)


def test_reversing_an_edge_changes_directed_pagerank():
    string = _string_graph([("A", "B"), ("B", "C")])
    forward = build_hybrid_directed_graph(string, _dataset([("A", "B"), ("B", "C")])).graph
    reverse = build_hybrid_directed_graph(string, _dataset([("B", "A"), ("B", "C")])).graph
    assert compute_directed_pagerank(forward) != compute_directed_pagerank(reverse)


def test_final_acceptance_direction_enters_pagerank_perturbation_and_redistribution():
    string = _string_graph([("A", "B"), ("B", "C"), ("A", "C")], [1.0, 0.7, 0.2])
    hybrid = build_hybrid_directed_graph(string, _dataset([("A", "B"), ("B", "C"), ("A", "C")]))
    assert hybrid.graph.get_eid("A", "B", directed=True, error=False) != -1
    assert hybrid.graph.get_eid("B", "A", directed=True, error=False) == -1
    baseline = compute_directed_pagerank(hybrid.graph)
    result = run_directed_calculation(hybrid, targets=["B"], bc_sample_sources=None)
    rows = result.report.set_index("entity")
    assert rows.loc["A", "Directed_PageRank_Baseline"] == baseline["A"]
    assert result.metadata["suppressed_arc_count"] == 2
    assert result.metadata["directed_perturbation_strategy"] == "incident_edge_attenuation"
    assert result.metadata["suppression_factor"] == 0.001
    assert result.metadata["target_pagerank"]["B"]["perturbed"] < result.metadata["target_pagerank"]["B"]["baseline"]
    assert np.isfinite(rows["Directed_Redistribution_Pct"]).all()
    assert (rows["Directed_Delta_PageRank"].abs() > 0).any()


def test_comparison_reports_validation_and_rank_shift():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C"), ("A", "C")]),
        _dataset([("A", "B"), ("B", "C"), ("C", "A")]),
    )
    directed = run_directed_calculation(hybrid, targets=["B"], bc_sample_sources=None)
    classic = pd.DataFrame({"gene": ["A", "C"], "Delta_PageRank_Pct": [2.0, -1.0]})
    comparison = compare_classic_and_directed(classic, directed)
    assert {"Response_Delta", "Classic_Rank", "Directed_Rank", "Rank_Shift", "Directed_BC"}.issubset(comparison.table)
    assert "spearman_rank_correlation" in comparison.validation
    assert "overlap_at_50" in comparison.validation


def test_cache_key_changes_for_every_scientific_context_input():
    base = dict(tissue_context="Lung", string_dataset_version="STRING-v12", omnipath_fingerprint="abc")
    key = directed_graph_cache_key(**base)
    assert key != directed_graph_cache_key(**{**base, "tissue_context": "Brain"})
    assert key != directed_graph_cache_key(**{**base, "string_dataset_version": "STRING-v13"})
    assert key != directed_graph_cache_key(**{**base, "omnipath_fingerprint": "def"})
    assert key != directed_graph_cache_key(**base, direction_policy="other")
    assert key != directed_graph_cache_key(**base, fallback_policy="other")
    assert key != directed_graph_cache_key(**base, taxon_id=10090)
    assert key != directed_graph_cache_key(**base, effective_graph_fingerprint="changed")
    assert key != directed_graph_cache_key(**base, edge_evidence_weighting_mode="structural_only")


def test_effective_graph_fingerprint_changes_when_weight_changes_without_count_change():
    graph = _string_graph([("A", "B"), ("B", "C")], [0.5, 0.8])
    before = effective_graph_fingerprint(graph)
    graph.es[0]["weight"] = 0.51
    assert effective_graph_fingerprint(graph) != before


def test_hybrid_cache_uses_content_fingerprint_and_context_not_object_identity():
    import src.services.directed_analysis as service

    service._HYBRID_GRAPH_CACHE.clear()
    graph = _string_graph([("A", "B"), ("B", "C")], [0.5, 0.8])
    dataset = _dataset([("A", "B")])
    first_fingerprint = effective_graph_fingerprint(graph)
    first = service._cached_hybrid_graph(
        graph=graph, dataset=dataset, tissue="Lung", string_dataset_version="fixture",
        effective_graph_fingerprint=first_fingerprint, taxon_id=9606,
        edge_evidence_weighting_mode="legacy_edge_modifiers",
    )
    graph.es[0]["weight"] = 0.51
    graph.es[0]["distance"] = 1.0 / 0.51
    second_fingerprint = effective_graph_fingerprint(graph)
    second = service._cached_hybrid_graph(
        graph=graph, dataset=dataset, tissue="Lung", string_dataset_version="fixture",
        effective_graph_fingerprint=second_fingerprint, taxon_id=9606,
        edge_evidence_weighting_mode="legacy_edge_modifiers",
    )
    other_context = service._cached_hybrid_graph(
        graph=graph, dataset=dataset, tissue="Brain", string_dataset_version="fixture",
        effective_graph_fingerprint=second_fingerprint, taxon_id=9606,
        edge_evidence_weighting_mode="legacy_edge_modifiers",
    )
    assert first is not second
    assert second is not other_context
    assert first.provenance.cache_key != second.provenance.cache_key
    assert second.provenance.cache_key != other_context.provenance.cache_key


def test_comparison_requires_full_response_and_never_fills_missing_with_zero():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C")]), _dataset([]),
    )
    directed = run_directed_calculation(hybrid, targets=["B"], bc_sample_sources=0)
    incomplete = pd.DataFrame({"gene": ["A"], "Delta_PageRank_Pct": [1.0]})
    with np.testing.assert_raises_regex(ValueError, "COMPARISON_INCOMPLETE"):
        compare_classic_and_directed(incomplete, directed)


def test_incident_attenuation_changes_downstream_response_when_incoming_alternatives_exist():
    string = _string_graph([
        ("U", "T"), ("U", "X"), ("T", "D"), ("X", "D"), ("D", "E"),
    ])
    dataset = _dataset([
        ("U", "T"), ("U", "X"), ("T", "D"), ("X", "D"), ("D", "E"),
    ])
    hybrid = build_hybrid_directed_graph(string, dataset)
    result = run_directed_calculation(hybrid, targets=["T"], bc_sample_sources=None)
    target = result.metadata["target_pagerank"]["T"]
    assert target["perturbed"] < target["baseline"]
    downstream = result.report.set_index("entity").loc[["D", "E"], "Directed_Delta_PageRank"]
    assert (downstream.abs() > 0).any()


def test_directed_failure_is_isolated_from_classic_report(tmp_path):
    classic = pd.DataFrame({"gene": ["A"], "Delta_PageRank_Pct": [2.5]})
    before = classic.copy(deep=True)
    bundle = run_optional_directed_engine(
        graph=_string_graph([("A", "B")]), targets=["A"], tissue="Lung",
        classic_report=classic, project_root=tmp_path, mode="Compare",
        block_weight_fraction=0.001, damping=0.85,
    )
    assert bundle.calculation.status.value == "UNAVAILABLE"
    pd.testing.assert_frame_equal(classic, before, check_exact=True)


def test_strategy_audit_exposes_normalization_and_competing_route_effect():
    audit = run_audit()
    chain = audit["graphs"]["one_way_chain"]
    assert all(item["outgoing_only_is_normalized_away"] for item in chain.values())
    competing_b = audit["graphs"]["competing_route"]["B"]
    assert competing_b["incident_equals_incoming_for_pagerank"]
    assert competing_b["incident_edge_attenuation"]["target_change"] < 0
    assert competing_b["incident_edge_attenuation"]["total_redistribution_l1"] > 0
    acceptance = audit["graphs"]["acceptance_alternatives"]["T"]
    assert acceptance["incident_edge_attenuation"]["target_change"] < 0
    assert acceptance["incident_edge_attenuation"]["downstream_absolute_change"] > 0


def test_mixed_hybrid_policy_keeps_resolved_and_fallback_arcs_distinct():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C"), ("C", "D")]),
        _dataset([("A", "B"), ("C", "D")]),
    )
    graph = hybrid.graph
    assert graph.get_eid("B", "A", directed=True, error=False) == -1
    assert graph.get_eid("B", "C", directed=True, error=False) != -1
    assert graph.get_eid("C", "B", directed=True, error=False) != -1
    assert graph.get_eid("D", "C", directed=True, error=False) == -1
    middle = graph.es[graph.get_eid("B", "C", directed=True)]
    assert middle["direction_source"] == "undirected_fallback"


def test_sign_is_metadata_only_and_does_not_change_calculation():
    string = _string_graph([("A", "B"), ("B", "C")], [0.7, 0.4])
    activating = _dataset([("A", "B"), ("B", "C")])
    inhibiting = replace(
        activating,
        edges=tuple(replace(
            edge, stimulation=False, inhibition=True, canonical_sign=-1,
            sign_status=SignStatus.INHIBITION,
        ) for edge in activating.edges),
        fingerprint="inhibiting-fixture",
    )
    activation_graph = build_hybrid_directed_graph(string, activating).graph
    inhibition_graph = build_hybrid_directed_graph(string, inhibiting).graph
    assert activation_graph.es["weight"] == inhibition_graph.es["weight"]
    assert compute_directed_pagerank(activation_graph) == compute_directed_pagerank(inhibition_graph)


def test_directed_hop_distance_is_canonical_and_legacy_alias_matches():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C")]), _dataset([("A", "B"), ("B", "C")]),
    )
    result = run_directed_calculation(hybrid, targets=["A"], bc_sample_sources=0)
    rows = result.report.set_index("entity")
    assert rows.loc["B", "Directed_Hop_Distance"] == 1
    assert rows.loc["C", "Directed_Hop_Distance"] == 2
    assert rows["Directed_Hop_Distance"].equals(rows["Directed_Distance"])
    assert result.metadata["directed_distance_semantics"] == "minimum_unweighted_outgoing_hops"


def test_compare_export_retains_engine_strategy_coverage_and_rank_data():
    hybrid = build_hybrid_directed_graph(
        _string_graph([("A", "B"), ("B", "C"), ("A", "C")]),
        _dataset([("A", "B"), ("B", "C"), ("A", "C")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["B"], bc_sample_sources=None)
    classic = pd.DataFrame({"gene": ["A", "C"], "Delta_PageRank_Pct": [1.0, -1.0]})
    from src.services.directed_analysis import DirectedRunBundle
    bundle = DirectedRunBundle(calculation, compare_classic_and_directed(classic, calculation))
    tables = directed_export_tables(bundle)
    assert {"Classic vs Directed", "Directed Engine Provenance"}.issubset(tables)
    assert {"Classic_Rank", "Directed_Rank", "Rank_Shift"}.issubset(tables["Classic vs Directed"])
    provenance = tables["Directed Engine Provenance"].iloc[0]
    assert provenance["directed_perturbation_strategy"] == "incident_edge_attenuation"
    assert float(provenance["direction_coverage"]) >= 0.0


def test_directed_presentation_is_single_source_for_ui_aliases_and_export():
    from src.services.engine_presentation import build_directed_presentation

    hybrid = build_hybrid_directed_graph(
        _string_graph([("TARGET", "ANO2"), ("ANO2", "RFLNB"), ("TARGET", "RFLNB")]),
        _dataset([("TARGET", "ANO2"), ("ANO2", "RFLNB"), ("TARGET", "RFLNB")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["TARGET"], bc_sample_sources=None)
    scores = pd.DataFrame({
        "gene": ["TARGET", "ANO2", "RFLNB"],
        "Lokalizasyon": ["membrane", "membrane", "cytoplasm"],
    })
    presentation = build_directed_presentation(
        calculation, targets=["TARGET"], scores=scores,
        gene_to_symbol={"TARGET": "CFTR", "ANO2": "ANO2"}, candidate_limit=20,
        minimum_abs_response_pct=0.0,
    )
    assert set(presentation.report["Engine"]) == {"directed"}
    assert set(presentation.report["Redistribution_Selection_Mode"]) == {"top_n"}
    assert presentation.metadata["candidate_count"] == len(presentation.candidates)
    assert "protein_id" in presentation.export
    assert presentation.export.set_index("protein_id").loc["ANO2", "gene"] == "ANO2"
    assert presentation.export.set_index("protein_id").loc["RFLNB", "gene"] == "Unresolved"
    assert presentation.export.set_index("protein_id").loc["RFLNB", "Symbol_Mapping_Status"] == "UNRESOLVED"


def test_active_result_flow_never_mixes_classic_and_directed_candidates():
    from src.services.directed_analysis import DirectedRunBundle
    from src.services.engine_presentation import build_directed_presentation, select_active_result_flow

    hybrid = build_hybrid_directed_graph(
        _string_graph([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
        _dataset([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["T"], bc_sample_sources=None)
    presentation = build_directed_presentation(
        calculation, targets=["T"], scores=pd.DataFrame({"gene": ["T", "D1", "D2"]}),
        gene_to_symbol={"T": "TARGET", "D1": "DIRECTED_ONE", "D2": "DIRECTED_TWO"},
        candidate_limit=20, minimum_abs_response_pct=0.0,
    )
    bundle = DirectedRunBundle(calculation=calculation, presentation=presentation)
    classic = pd.DataFrame({"gene": ["T", "CLASSIC_ONLY"], "Hasar_Tipi": ["Birincil_Hasar_Hedef", "Üçüncül_Hasar_Stresli"]})

    directed = select_active_result_flow(
        classic_report=classic, classic_signed_response=classic,
        requested_mode="Directed", directed_bundle=bundle,
    )
    assert directed.engine == "directed"
    assert "CLASSIC_ONLY" not in set(directed.report["gene"])
    assert directed.signed_response is not None
    assert set(directed.signed_response["Engine"]) == {"directed"}

    for mode in ("Classic", "Compare"):
        selected = select_active_result_flow(
            classic_report=classic, classic_signed_response=classic,
            requested_mode=mode, directed_bundle=bundle,
        )
        assert selected.engine == "classic"
        pd.testing.assert_frame_equal(selected.report, classic)


def test_candidate_cache_identity_is_engine_sensitive_even_for_same_gene_set():
    from src.services.engine_presentation import candidate_set_fingerprint

    genes = ["ENSP_A", "ENSP_B"]
    assert candidate_set_fingerprint("classic", genes) != candidate_set_fingerprint("directed", genes)
    assert candidate_set_fingerprint("directed", genes) == candidate_set_fingerprint("Directed", reversed(genes))


def test_ui_modes_expose_their_actual_base_scientific_motor():
    from src.services.engine_presentation import base_engine_for_mode

    assert base_engine_for_mode("Classic") == "classic"
    assert base_engine_for_mode("Directed") == "classic"
    assert base_engine_for_mode("Compare") == "classic"
    assert base_engine_for_mode("Evidence") == "evidence"


def test_directed_candidates_feed_interpretation_and_enrichment_without_classic_rows():
    from src.interpretation import interpret_forward_results
    from src.services.enrichment import local_over_representation
    from src.services.engine_presentation import build_directed_presentation

    hybrid = build_hybrid_directed_graph(
        _string_graph([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
        _dataset([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["D1"], bc_sample_sources=None)
    mapping = {"T": "TARGET", "D1": "DIRECTED_ONE", "D2": "DIRECTED_TWO"}
    presentation = build_directed_presentation(
        calculation, targets=["D1"], scores=pd.DataFrame({"gene": list(mapping)}),
        gene_to_symbol=mapping, candidate_limit=20, minimum_abs_response_pct=0.0,
    )
    interpretation = interpret_forward_results(
        report=presentation.report, candidates=presentation.candidates,
        loss_candidates=presentation.losses, targets=["D1"], gene_to_symbol=mapping,
        tissue="fixture", organism="Homo sapiens",
    )
    assert interpretation.top_candidates[0]["gene"] in set(presentation.candidates["gene_symbol"])
    assert all(item["gene"] != "CLASSIC_ONLY" for item in interpretation.top_candidates)

    selected = presentation.candidates["gene_symbol"].tolist()
    enrichment = local_over_representation(
        selected, [*selected, "CLASSIC_ONLY"],
        {"Directed term": selected, "Classic term": ["CLASSIC_ONLY"]},
    )
    assert enrichment["Term"].tolist() == ["Directed term"]


def test_directed_ui_projection_has_exact_contract_and_real_directed_values():
    from src.services.engine_presentation import DIRECTED_UI_COLUMNS, build_directed_presentation, directed_ui_frame

    hybrid = build_hybrid_directed_graph(
        _string_graph([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
        _dataset([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["D2"], bc_sample_sources=None)
    presentation = build_directed_presentation(
        calculation, targets=["D2"],
        scores=pd.DataFrame({"gene": ["T", "D1", "D2"], "Lokalizasyon": ["nucleus", "membrane", "cytoplasm"]}),
        gene_to_symbol={"T": "TARGET", "D1": "DIRECTED_ONE", "D2": "DIRECTED_TWO"},
        candidate_limit=20, minimum_abs_response_pct=0.0,
    )
    shown = directed_ui_frame(presentation)
    assert tuple(shown.columns) == DIRECTED_UI_COLUMNS
    source = presentation.full_response.set_index("gene_symbol")
    row = shown.set_index("Gen").loc["DIRECTED_ONE"]
    assert row["ENSP"] == "D1"
    assert row["Directed % değişim"] == source.loc["DIRECTED_ONE", "Delta_PageRank_Pct"]
    assert row["Direction"] == source.loc["DIRECTED_ONE", "Response_Direction"]
    assert row["Rank"] == source.loc["DIRECTED_ONE", "Directed_Rank"]
    assert row["Directed BC"] == source.loc["DIRECTED_ONE", "BC_Skoru"]
    assert row["Localization"] == "membrane"


def test_compare_ui_projection_visibly_shows_classic_directed_difference():
    from src.services.engine_presentation import COMPARE_UI_COLUMNS, build_directed_presentation, compare_ui_frame

    hybrid = build_hybrid_directed_graph(
        _string_graph([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
        _dataset([("T", "D1"), ("D1", "D2"), ("T", "D2")]),
    )
    calculation = run_directed_calculation(hybrid, targets=["D2"], bc_sample_sources=None)
    presentation = build_directed_presentation(
        calculation, targets=["D2"],
        scores=pd.DataFrame({"gene": ["T", "D1", "D2"], "Lokalizasyon": ["nucleus", "membrane", "cytoplasm"]}),
        gene_to_symbol={"T": "TARGET", "D1": "DIRECTED_ONE", "D2": "DIRECTED_TWO"},
        candidate_limit=20, minimum_abs_response_pct=0.0,
    )
    classic = pd.DataFrame({"gene": ["T", "D1"], "Delta_PageRank_Pct": [1.0, -2.0]})
    comparison = compare_classic_and_directed(classic, calculation)
    shown = compare_ui_frame(comparison.table, presentation)
    assert tuple(shown.columns) == COMPARE_UI_COLUMNS
    assert (shown["Fark"].abs() > 0).any()
    row = shown.set_index("Gen").loc["DIRECTED_ONE"]
    assert row["ENSP"] == "D1"
    assert row["Classic % değişim"] == -2.0
    assert row["Directed % değişim"] != row["Classic % değişim"]
    assert row["Localization"] == "membrane"
