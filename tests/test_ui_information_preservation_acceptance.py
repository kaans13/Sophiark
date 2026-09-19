"""Acceptance checks for the UI information-preservation audit.

These checks deliberately exercise presentation and bundle boundaries only.  They
must never invoke either scientific engine.
"""

from __future__ import annotations

import io
from math import ceil
from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.product.analysis_bundle import AnalysisBundleReader, BundleError
from src import state
from src.ui.metric_presentation import prepare_display_frame
from src.ui.tables import dataframe_page, export_bytes
from src.ui.unified_research_report import render_unified_research_report
from src.unified.service import UnifiedResearchReport, UnifiedStatus
from streamlit.testing.v1 import AppTest


def test_display_projection_never_changes_canonical_values_or_columns() -> None:
    source = pd.DataFrame({
        "entity_id": ["ENSP_A"],
        "Delta_PageRank_Pct": [-12.3456789],
        "q_value": [0.004321],
        "significant_redistribution": [True],
    })
    before = source.copy(deep=True)

    display = prepare_display_frame(source)

    assert set(display.columns) == set(source.columns)
    assert source.loc[0, "entity_id"] == "ENSP_A"
    assert source.loc[0, "Delta_PageRank_Pct"] == -12.3456789
    assert source.loc[0, "q_value"] == 0.004321
    assert bool(source.loc[0, "significant_redistribution"])
    assert_frame_equal(source, before)


def test_unified_surface_renders_raw_tables_and_does_not_mutate_them() -> None:
    candidates = pd.DataFrame({
        "entity_id": ["ENSP_A", "ENSP_B", "ENSP_LOW"],
        "finding_class": ["ROBUST_CROSS_MODEL", "DIRECTION_SENSITIVE", "LOW_AGREEMENT"],
        "classic_response": [-0.25, 0.0, 0.75],
        "directed_response": [-0.5, 0.0, -0.25],
        "Evidence_Provenance_Status": ["AVAILABLE", "UNAVAILABLE", "AVAILABLE"],
    })
    classic = pd.DataFrame({"gene": ["ENSP_A"], "Delta_PageRank_Pct": [-12.3456789], "q_value": [0.004321]})
    directed = pd.DataFrame({"entity": ["ENSP_A"], "Directed_Redistribution_Pct": [-2.5], "Directed_BC": [0.0]})
    report = UnifiedResearchReport(
        target="ENSP_TARGET", tissue="Lung", status=UnifiedStatus.COMPLETE,
        classic_status="COMPLETE", directed_status="AVAILABLE", candidates=candidates,
        classic_raw=classic, directed_raw=directed, complex_context=pd.DataFrame(),
        agreement={"sign_agreement_pct": 98.1}, provenance={}, context_availability={}, errors={},
    )
    before = (candidates.copy(deep=True), classic.copy(deep=True), directed.copy(deep=True))
    app = AppTest.from_string('''\
import pandas as pd
from src.unified.service import UnifiedResearchReport, UnifiedStatus
from src.ui.unified_research_report import render_unified_research_report
report = UnifiedResearchReport(
    target="ENSP_TARGET", tissue="Lung", status=UnifiedStatus.COMPLETE,
    classic_status="COMPLETE", directed_status="AVAILABLE",
    candidates=pd.DataFrame({"entity_id":["ENSP_A","ENSP_B","ENSP_LOW"], "finding_class":["ROBUST_CROSS_MODEL","DIRECTION_SENSITIVE","LOW_AGREEMENT"], "classic_response":[-0.25,0.0,0.75], "directed_response":[-0.5,0.0,-0.25], "Evidence_Provenance_Status":["AVAILABLE","UNAVAILABLE","AVAILABLE"]}),
    classic_raw=pd.DataFrame({"gene":["ENSP_A"], "Delta_PageRank_Pct":[-12.3456789], "q_value":[0.004321]}),
    directed_raw=pd.DataFrame({"entity":["ENSP_A"], "Directed_Redistribution_Pct":[-2.5], "Directed_BC":[0.0]}),
    complex_context=pd.DataFrame(), agreement={"sign_agreement_pct":98.1}, provenance={}, context_availability={}, errors={},
)
render_unified_research_report(report)
''')
    app.run(timeout=20)

    assert not app.exception
    assert {tab.label for tab in app.tabs} >= {"Düşük uyum", "Tüm adaylar", "Evidence provenance", "Ham motor sonuçları"}
    assert len(app.number_input) == 4
    assert_frame_equal(candidates, before[0])
    assert_frame_equal(classic, before[1])
    assert_frame_equal(directed, before[2])


def test_real_bundle_keeps_canonical_raw_tables_and_rejects_invalid_fixture() -> None:
    root = Path(__file__).resolve().parents[1]
    saved = AnalysisBundleReader().load(root / "outputs" / "analysis_bundle_validation" / "CFTR_Lung.sophiark")

    assert len(saved.report.candidates) == 100
    assert len(saved.report.classic_raw) == len(saved.report.directed_raw) == 16_954
    assert {"entity_id", "classic_response", "directed_response", "finding_class"} <= set(saved.report.candidates)
    assert {"gene", "Delta_PageRank_Pct", "q_value"} <= set(saved.report.classic_raw)
    assert {"entity", "Directed_Redistribution_Pct", "Directed_BC"} <= set(saved.report.directed_raw)

    with pytest.raises(BundleError):
        AnalysisBundleReader().load(root / "outputs" / "ui_acceptance" / "fixtures" / "invalid_bundle.sophiark")


@pytest.mark.parametrize("raw_name,identity_column,value_column", [
    ("classic_raw", "gene", "Delta_PageRank_Pct"),
    ("directed_raw", "entity", "Directed_Redistribution_Pct"),
])
def test_full_raw_inspector_reaches_first_boundary_beyond_500_and_final_rows(
    raw_name: str, identity_column: str, value_column: str,
) -> None:
    root = Path(__file__).resolve().parents[1]
    raw = getattr(AnalysisBundleReader().load(root / "outputs" / "analysis_bundle_validation" / "CFTR_Lung.sophiark").report, raw_name)
    page_size = 500
    pages = ceil(len(raw) / page_size)

    first = dataframe_page(raw, page=1, page_size=page_size)
    boundary = dataframe_page(raw, page=2, page_size=page_size)
    final = dataframe_page(raw, page=pages, page_size=page_size)

    assert len(raw) == 16_954
    assert first.start == 0 and first.frame.iloc[0][identity_column] == raw.iloc[0][identity_column]
    assert first.frame.iloc[-1][identity_column] == raw.iloc[499][identity_column]
    assert boundary.start == 500 and boundary.frame.iloc[0][identity_column] == raw.iloc[500][identity_column]
    assert boundary.frame.iloc[0][value_column] == raw.iloc[500][value_column]
    assert final.end == len(raw) and final.frame.iloc[-1][identity_column] == raw.iloc[-1][identity_column]
    assert final.frame.iloc[-1][value_column] == raw.iloc[-1][value_column]

    inspected = pd.concat([dataframe_page(raw, page=page, page_size=page_size).frame for page in range(1, pages + 1)], ignore_index=True)
    assert_frame_equal(inspected, raw.reset_index(drop=True))


def test_low_agreement_candidates_have_full_inspection_and_canonical_export() -> None:
    candidates = pd.DataFrame({
        "entity_id": ["ENSP_LOW_1", "ENSP_LOW_2", "ENSP_OTHER"],
        "symbol": ["LOW1", "LOW2", "OTHER"],
        "classic_response": [-0.125, 0.75, 0.25],
        "directed_response": [0.5, -0.25, 0.1],
        "classic_absolute_response_rank": [501, 502, 1],
        "directed_absolute_response_rank": [601, 602, 2],
        "finding_class": ["LOW_AGREEMENT", "LOW_AGREEMENT", "ROBUST_CROSS_MODEL"],
    })
    low = candidates.loc[candidates["finding_class"].eq("LOW_AGREEMENT")].copy()
    page = dataframe_page(low, page=1, page_size=100)
    payload, _, extension = export_bytes(low, "CSV (.csv)")
    exported = pd.read_csv(io.BytesIO(payload))

    assert len(low) == len(page.frame) == len(exported) == 2
    assert list(page.frame["entity_id"]) == ["ENSP_LOW_1", "ENSP_LOW_2"]
    assert list(exported["entity_id"]) == ["ENSP_LOW_1", "ENSP_LOW_2"]
    assert list(exported["classic_response"]) == [-0.125, 0.75]
    assert list(exported["directed_absolute_response_rank"]) == [601, 602]
    assert extension == "csv"


def test_paginated_ui_control_reaches_a_final_row_without_reordering() -> None:
    app = AppTest.from_string('''\
import pandas as pd
from src.ui.tables import render_dataframe
frame = pd.DataFrame({"entity_id":[f"ENSP_{i:04d}" for i in range(1001)], "rank":range(1, 1002), "response":[i / 10 for i in range(1001)]})
render_dataframe(frame, label="pagination_acceptance", key="pagination_acceptance", row_limit=500, paginate=True)
''')
    app.run(timeout=20)
    assert app.number_input[0].value == 1
    assert app.dataframe[0].value.iloc[0]["entity_id"] == "ENSP_0000"

    app.number_input[0].set_value(3).run(timeout=20)
    final = app.dataframe[0].value
    assert len(final) == 1
    assert final.iloc[0]["entity_id"] == "ENSP_1000"
    assert final.iloc[0]["rank"] == 1001
    assert final.iloc[0]["response"] == 100.0


def test_invalid_bundle_load_clears_transient_display_and_valid_load_recovers() -> None:
    session = {"saved_unified_result": "previous-valid-report", "durable_history": ["A.sophiark"]}

    with patch.object(state.st, "session_state", session):
        with pytest.raises(BundleError):
            state.load_saved_unified_result(lambda: AnalysisBundleReader().load(b"not a bundle"))
        assert session["saved_unified_result"] is None
        assert session["durable_history"] == ["A.sophiark"]

        result = state.load_saved_unified_result(lambda: "valid-C-report")
        assert result == session["saved_unified_result"] == "valid-C-report"


def test_supported_ui_entrypoint_is_root_app_not_internal_page_script() -> None:
    root = Path(__file__).resolve().parents[1]
    readme = (root / "readme.md").read_text(encoding="utf-8")
    page_source = (root / "src" / "ui" / "pages" / "01_analiz.py").read_text(encoding="utf-8")

    assert "streamlit run app.py" in readme
    assert "NOT_SUPPORTED_ENTRYPOINT" in page_source


def test_stable_output_summary_and_candidate_cards_remain_directly_accessible() -> None:
    root = Path(__file__).resolve().parents[1]
    app_source = (root / "app.py").read_text(encoding="utf-8")

    # The stable release exposed these scientific readouts in its primary result
    # summary.  A compact overview may coexist, but it must not replace them.
    for readout in (
        'ui_t("mean_local_efficiency_change")',
        'ui_t("global_efficiency_change")',
        'ui_t("target_local_efficiency_change")',
        'summary_cell_html("Pozitif yanıt"',
        'summary_cell_html("Negatif yanıt"',
        'summary_cell_html("En güçlü kayıp"',
        'summary_cell_html("Motor / anlamlılık"',
        "target_html=_hedef_ozet_html",
        'context_label=f"{doku_label} · {_summary_mode_label}"',
    ):
        assert readout in app_source

    assert 'st.markdown("### Ayrıntılı çıktı özeti")' in app_source
    assert 'title="Kararlı sonuç ölçümleri"' in app_source

    network_block = app_source[
        app_source.index("with _network_tab:") : app_source.index("with _biology_tab:")
    ]
    assert "_candidate_surface = st.container()" in network_block
    assert "_cards_surface = st.container()" in network_block
    assert "_system_response_surface = st.container()" in network_block
    assert 'st.expander("Hedef ve öne çıkan aday kartları"' not in app_source

    interpretation_source = (root / "src" / "ui" / "biological_interpretation.py").read_text(encoding="utf-8")
    for stable_visible_output in (
        'label="biyolojik_negatif_ozet"',
        'st.markdown("### Sistem düzeyi ağ yanıtı")',
        "st.write(interpretation.system_response)",
        "render_scientific_note(interpretation.limitations)",
    ):
        assert stable_visible_output in interpretation_source

    # Keep the redesign as a clearer projection of the stable release, not as a
    # place to expose additional scientific result datasets.
    assert "Tam işaretli ağ yanıtı · seçilmemiş kayıtlar dahil" not in app_source
    assert 'label="full_signed_response"' not in app_source


def test_primary_workspace_has_no_permanent_sidebar_and_keeps_compact_navigation() -> None:
    root = Path(__file__).resolve().parents[1]
    app_source = (root / "app.py").read_text(encoding="utf-8")
    workspace_source = (root / "src" / "ui" / "research_workspace.py").read_text(encoding="utf-8")

    assert "st.sidebar" not in app_source
    assert "_config_panel = st.popover(" in app_source
    assert '"Analizi değiştir" if st.session_state.get("rapor_df") is not None else "Yeni analiz"' in app_source
    assert "_shell_actions = st.columns" not in app_source
    assert 'st.session_state.get("rapor_df") is None and SOK_RAPORU.exists()' not in app_source
    assert app_source.index("if _result_is_stale:") < app_source.index('_network_status_slot.status("Graf hazırlanıyor..."')
    assert "if not ates_buton:\n        st.stop()" in app_source
    for destination in (
        "1 · Sonucu Anla", "2 · Ağı Keşfet", "3 · Kanıtı Doğrula",
        "Ağ Yanıtı", "Biyolojik Bağlam", "Kanıt & Kaynaklar",
        "Ham Veri & İndirme", "İleri Araçlar",
    ):
        assert destination in workspace_source
