"""Read-only presentation boundaries; local browser acceptance is separate."""
import pandas as pd
from pandas.testing import assert_frame_equal

from src.ui.research_workspace import (
    _candidate_rows_html,
    evidence_html,
    matching_enrichment,
    search_results,
)
from src.state import build_simulation_result_context, simulation_result_context_matches


def test_search_is_literal_order_preserving_and_does_not_reduce_source_schema():
    source = pd.DataFrame({"gene": ["A.1", "AX1", "B"], "precise": [1.23456789123, -0.0, 3.0], "context": ["x", "y", "A.1"]})
    before = source.copy(deep=True)
    result = search_results(source, "A.1")
    assert list(result.index) == [0, 2]
    assert list(result.columns) == list(source.columns)
    assert result.iloc[0]["precise"] == source.iloc[0]["precise"]
    assert_frame_equal(source, before)
    assert search_results(source, "no match").empty
    assert_frame_equal(search_results(source, ""), source)


def test_pathway_evidence_requires_exact_recorded_membership():
    source = pd.DataFrame({"Genes": ["CFTR;ANO2", "CFTR2;ANO2", "CFTR-AS1", ["CFTR", "B"]], "Adjusted P-value": [.01, .02, .03, .04]})
    before = source.copy(deep=True)
    matches = matching_enrichment(source, "CFTR")
    assert list(matches.index) == [0, 3]
    assert_frame_equal(source, before)
    assert matching_enrichment(None, "CFTR").empty


def test_untrusted_annotations_are_escaped_and_missing_evidence_is_not_invented():
    html = evidence_html(pd.Series({"Lokalizasyon": "<script>alert(1)</script>", "q_value": None}), ("Lokalizasyon", "q_value"))
    assert "<script>" not in html
    assert "&lt;script&gt;" in html
    assert "q_value" not in html
    assert "kanıt yok" in evidence_html(pd.Series(dtype=object), ("q_value",))


def test_candidate_object_rows_escape_annotations_and_keep_signed_values():
    source = pd.DataFrame({
        "gene": ["ENSP1", "ENSP2"],
        "Symbol": ["<script>A</script>", "SAFE"],
        "Delta_PageRank_Pct": [-2.5, 1.25],
        "Hinterland_Skoru": [10.0, 20.0],
        "Lokalizasyon": ["Cell_Membrane", "Nucleus"],
    })
    before = source.copy(deep=True)
    html = _candidate_rows_html(source, identity="gene", symbol_field="Symbol")
    assert "<script>" not in html
    assert "&lt;script&gt;A&lt;/script&gt;" in html
    assert "-2.50%" in html and "+1.25%" in html
    assert html.index("ENSP1") < html.index("ENSP2")
    assert_frame_equal(source, before)


def test_betweenness_accuracy_is_part_of_completed_result_identity():
    inputs = dict(species="İnsan", tissue="Lung", engine="Classic", targets=("CFTR",), localization="SNIPER_MODU", block_strength=.001, damping=.85, candidate_limit=20, tissue_normalization_mode="within_tissue", bc_sample_sources=1200)
    context = build_simulation_result_context(**inputs)
    assert simulation_result_context_matches(context, **inputs)
    assert not simulation_result_context_matches(context, **{**inputs, "bc_sample_sources": 3000})
