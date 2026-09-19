"""Regression tests for the read-only network response presentation layer."""

from __future__ import annotations

import io
import unittest
from pathlib import Path
from unittest.mock import Mock, patch

import numpy as np
import pandas as pd
from streamlit.testing.v1 import AppTest

from src import state as app_state
from src.core import cache_manager
from src.ui.metric_presentation import format_metric, metric_label, metric_printf
from src.ui.network_response import (
    _analytical_rows_html,
    _display_table,
    _nonzero_responses,
    build_network_response_matrix,
)
from src.ui.tables import export_bytes


class NetworkResponsePresentationTests(unittest.TestCase):
    @staticmethod
    def _redistribution(size: int) -> pd.DataFrame:
        return pd.DataFrame({
            "gene": [f"GENE_{index}" for index in range(size)],
            "Symbol": [f"G{index}" for index in range(size)],
            "Delta_PageRank_Pct": np.linspace(.1, 2.0, size),
            "Hinterland_Skoru": np.linspace(1.0, 95.0, size),
            "BC_Skoru": np.linspace(0.0, 1.0, size),
            "Gümrük_Kapisi": [index % 5 == 0 for index in range(size)],
            "Lokalizasyon": ["Nucleus"] * size,
        })

    def test_20_100_200_rows_are_presented_without_reselection(self) -> None:
        for size in (20, 100, 200):
            with self.subTest(size=size):
                source = self._redistribution(size)
                matrix = build_network_response_matrix(None, None, source)
                self.assertEqual(len(matrix), size)
                self.assertEqual(list(matrix["Protein kimliği"]), list(source["gene"]))

    def test_sources_are_not_mutated(self) -> None:
        target = self._redistribution(2)
        loss = self._redistribution(3).assign(Delta_PageRank_Pct=[-1.0, -2.0, -3.0])
        redistribution = self._redistribution(5)
        snapshots = [frame.copy(deep=True) for frame in (target, loss, redistribution)]
        build_network_response_matrix(target, loss, redistribution)
        for source, before in zip((target, loss, redistribution), snapshots):
            pd.testing.assert_frame_equal(source, before)

    def test_same_gene_gets_linked_roles_without_duplicate_source_changes(self) -> None:
        target = pd.DataFrame({"gene": ["A"], "Symbol": ["GENEA"]})
        loss = pd.DataFrame({"gene": ["A"], "Delta_PageRank_Pct": [-2.5]})
        redistribution = pd.DataFrame({"gene": ["A", "B"], "Delta_PageRank_Pct": [1.0, 2.0]})
        matrix = build_network_response_matrix(target, loss, redistribution)
        self.assertEqual(list(matrix["Protein kimliği"]), ["A", "B"])
        self.assertIn("Pertürbasyon hedefi", matrix.iloc[0]["Ağ yanıtı"])
        self.assertIn("Ağ önemi kaybı", matrix.iloc[0]["Ağ yanıtı"])
        self.assertIn("Yeniden dağılım", matrix.iloc[0]["Ağ yanıtı"])

    def test_missing_columns_nan_and_empty_frames_are_safe(self) -> None:
        minimal = pd.DataFrame({"gene": ["A", "B"], "Symbol": [np.nan, "GENEB"]})
        matrix = build_network_response_matrix(pd.DataFrame(), None, minimal)
        self.assertEqual(len(matrix), 2)
        self.assertEqual(matrix.iloc[0]["Gen / protein"], "A")
        self.assertTrue(pd.isna(matrix.iloc[0]["BC"]))
        self.assertTrue(build_network_response_matrix(None, None, None).empty)

    def test_empty_loss_and_empty_redistribution_do_not_remove_targets(self) -> None:
        target = pd.DataFrame({"gene": ["TARGET"], "Symbol": ["T"]})
        matrix = build_network_response_matrix(target, pd.DataFrame(), pd.DataFrame())
        self.assertEqual(list(matrix["Protein kimliği"]), ["TARGET"])
        self.assertEqual(matrix.iloc[0]["Ağ yanıtı"], "Pertürbasyon hedefi")

    def test_analytical_rows_balance_both_response_directions_without_recalculation(self) -> None:
        source = pd.DataFrame({
            "Gen / protein": ["TARGET", *[f"LOSS_{i}" for i in range(20)], *[f"GAIN_{i}" for i in range(20)]],
            "Protein kimliği": ["T", *[f"L{i}" for i in range(20)], *[f"G{i}" for i in range(20)]],
            "Ağ yanıtı": ["Pertürbasyon hedefi", *(["Ağ önemi kaybı"] * 20), *(["Yeniden dağılım"] * 20)],
            "Yanıt yönü": ["—", *(["↓ Ağ önemi azaldı"] * 20), *(["↑ Ağ önemi arttı"] * 20)],
            "Düzenleyici ilişki": ["Yön kanıtı yok"] * 41,
            "Ağ önemi değişimi (%)": [float("nan"), *[-float(i + 1) for i in range(20)], *[float(i + 1) for i in range(20)]],
            "Hinterland Skoru": [1.0] * 41,
            "Hücresel konum": ["Nucleus"] * 41,
        })
        before = source.copy(deep=True)
        html = _analytical_rows_html(source, limit=9)
        self.assertIn("LOSS_0", html)
        self.assertIn("GAIN_0", html)
        self.assertIn("4 azalan + 4 artan", html)
        self.assertNotIn("LOSS_10", html)
        pd.testing.assert_frame_equal(source, before)

    def test_a_second_simulation_view_has_no_rows_from_the_first(self) -> None:
        first = build_network_response_matrix(None, None, pd.DataFrame({"gene": ["OLD"]}))
        second = build_network_response_matrix(None, None, pd.DataFrame({"gene": ["NEW"]}))
        self.assertEqual(list(first["Protein kimliği"]), ["OLD"])
        self.assertEqual(list(second["Protein kimliği"]), ["NEW"])

    def test_scientific_values_pass_through_exactly(self) -> None:
        source = pd.DataFrame({
            "gene": ["A"],
            "Delta_PageRank_Pct": [-12.3456789],
            "Hinterland_Skoru": [64.86],
            "BC_Skoru": [26252.0],
            "Gümrük_Kapisi": [True],
            "significant_redistribution": [True],
        })
        matrix = build_network_response_matrix(None, source, None)
        self.assertEqual(matrix.iloc[0]["Ağ önemi değişimi (%)"], -12.3456789)
        self.assertEqual(matrix.iloc[0]["Hinterland Skoru"], 64.86)
        self.assertEqual(matrix.iloc[0]["BC"], 26252.0)
        self.assertEqual(matrix.iloc[0]["Compartment Bottleneck"], "✓")
        self.assertEqual(matrix.iloc[0]["İstatistiksel destek"], "✓")
        self.assertIn("İstatistiksel destek", matrix.iloc[0]["Ağ yanıtı"])

    def test_positive_and_negative_nonzero_responses_are_combined_without_mutation(self) -> None:
        losses = pd.DataFrame({
            "gene": ["LOSS", "ZERO_LOSS"],
            "Delta_PageRank_Pct": [-2.5, 0.0],
            "Kategori": ["Mavi", "Açık Kırmızı"],
        })
        gains = pd.DataFrame({
            "gene": ["GAIN", "ZERO_GAIN"],
            "Delta_PageRank_Pct": [3.5, 0.0],
            "Kategori": ["Sarı", "Mavi"],
        })
        losses_before = losses.copy(deep=True)
        gains_before = gains.copy(deep=True)
        combined = _nonzero_responses(losses, gains)
        self.assertEqual(list(combined["gene"]), ["LOSS", "GAIN"])
        self.assertEqual(list(combined["Delta_PageRank_Pct"]), [-2.5, 3.5])
        pd.testing.assert_frame_equal(losses, losses_before)
        pd.testing.assert_frame_equal(gains, gains_before)

    def test_location_separators_are_removed_only_in_the_matrix(self) -> None:
        source = pd.DataFrame({
            "gene": ["A"], "Lokalizasyon": ["Cell_Membrane|Endoplasmic_Reticulum"],
        })
        before = source.copy(deep=True)
        matrix = build_network_response_matrix(None, None, source)
        self.assertEqual(matrix.iloc[0]["Hücresel konum"], "Cell Membrane|Endoplasmic Reticulum")
        pd.testing.assert_frame_equal(source, before)

    def test_matrix_keeps_recorded_mygene_go_context_and_gateway_is_last(self) -> None:
        source = pd.DataFrame({
            "gene": ["ENSP00000423820"],
            "Symbol": ["ZFP62"],
            "Protein_Adi": ["Zinc finger protein 62"],
            "MyGene Adı": ["ZFP62"],
            "GO Biyolojik Süreç": [["DNA-templated transcription"]],
            "GO Moleküler İşlev": [["DNA binding"]],
            "GO Hücresel Bileşen": [["nucleus"]],
            "Gümrük_Kapisi": [True],
        })
        before = source.copy(deep=True)
        matrix = build_network_response_matrix(source, None, None)

        self.assertEqual(matrix.iloc[0]["Protein adı"], "Zinc finger protein 62")
        self.assertEqual(matrix.iloc[0]["MyGene Adı"], "ZFP62")
        self.assertEqual(matrix.iloc[0]["GO Biyolojik Süreç"], "DNA-templated transcription")
        self.assertEqual(matrix.columns[-1], "Compartment Bottleneck")
        pd.testing.assert_frame_equal(source, before)

    def test_source_exports_append_context_and_put_gateway_last_without_mutation(self) -> None:
        source = pd.DataFrame({
            "gene": ["A"],
            "Gümrük_Kapisi": [True],
            "Protein_Adi": ["Protein A"],
            "MyGene Adı": ["GENEA"],
            "GO Biyolojik Süreç": [["signal transduction"]],
        })
        before = source.copy(deep=True)
        display, export = _display_table(source, ["gene", "Gümrük_Kapisi"])

        self.assertEqual(list(export.columns), [
            "gene", "Protein_Adi", "MyGene Adı", "GO Biyolojik Süreç", "Compartment Bottleneck",
        ])
        self.assertEqual(export.columns[-1], "Compartment Bottleneck")
        self.assertEqual(display.iloc[0]["GO Biyolojik Süreç"], ["signal transduction"])
        pd.testing.assert_frame_equal(source, before)

    def test_matrix_csv_and_excel_exports_keep_context_and_gateway_order(self) -> None:
        source = pd.DataFrame({
            "gene": ["A"], "Symbol": ["GENEA"], "MyGene Adı": ["Gene A"],
            "GO Biyolojik Süreç": [["signal transduction"]], "Gümrük_Kapisi": [True],
        })
        matrix = build_network_response_matrix(source, None, None)
        for fmt, reader in (
            ("CSV (.csv)", lambda data: pd.read_csv(io.BytesIO(data))),
            ("Excel (.xlsx)", lambda data: pd.read_excel(io.BytesIO(data))),
        ):
            with self.subTest(fmt=fmt):
                payload, _, _ = export_bytes(matrix, fmt)
                restored = reader(payload)
                self.assertEqual(restored.columns[-1], "Compartment Bottleneck")
                self.assertEqual(restored.iloc[0]["MyGene Adı"], "Gene A")
                self.assertEqual(restored.iloc[0]["GO Biyolojik Süreç"], "signal transduction")

    def test_matrix_exposes_recorded_target_stressed_gene_pubmed_search(self) -> None:
        pubmed_url = "https://pubmed.ncbi.nlm.nih.gov/?term=%28CFTR%29+AND+SCNN1B"
        source = pd.DataFrame({
            "gene": ["STRESSED"], "Symbol": ["SCNN1B"],
            "PubMed_Makale": [pubmed_url], "Gümrük_Kapisi": [False],
        })
        matrix = build_network_response_matrix(None, None, source)

        self.assertEqual(matrix.iloc[0]["PubMed araması"], pubmed_url)
        self.assertEqual(matrix.columns[-1], "Compartment Bottleneck")


class ConditionalPresentationTests(unittest.TestCase):
    @staticmethod
    def _run(script: str) -> AppTest:
        app = AppTest.from_string(script)
        app.run(timeout=20)
        return app

    def test_go_present_kegg_absent_and_reverse_are_safe(self) -> None:
        for source in ("GO_Biological_Process_2021", "KEGG_2021_Human"):
            with self.subTest(source=source):
                app = self._run(f'''\
import pandas as pd
from src.ui.network_response import render_functional_context
frame = pd.DataFrame({{
    "Kaynak": ["{source}"], "Term": ["term"], "Overlap": ["2/10"],
    "Adjusted P-value": [0.01], "Genes": ["A;B"],
}})
render_functional_context(frame, analysis_ran=True)
''')
                self.assertEqual(len(app.exception), 0)
                self.assertEqual([tab.label for tab in app.tabs], ["GO", "KEGG"])
                self.assertGreaterEqual(len(app.dataframe), 1)

    def test_no_significant_fdr_candidates_is_an_explicit_empty_result(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_statistical_support
frame = pd.DataFrame(columns=["gene", "q_value", "significant_redistribution"])
render_statistical_support(frame, columns=list(frame.columns))
''')
        self.assertEqual(len(app.exception), 0)
        self.assertTrue(any("aday bulunamadı" in item.value for item in app.info))

    def test_functional_not_run_no_result_and_error_states_are_safe(self) -> None:
        cases = (
            ("render_functional_context(None, analysis_ran=False)", "info"),
            ("render_functional_context(None, analysis_ran=True)", "info"),
            (
                "render_functional_context(None, analysis_ran=True, notices=['Yolak analizi tamamlanamadı.'])",
                "warning",
            ),
        )
        for call, expected in cases:
            with self.subTest(call=call):
                app = self._run(f'''\
from src.ui.network_response import render_functional_context
{call}
''')
                self.assertEqual(len(app.exception), 0)
                self.assertGreaterEqual(len(getattr(app, expected)), 1)

    def test_go_and_kegg_are_both_kept_as_complete_tables(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_functional_context
frame = pd.DataFrame({
    "Kaynak": ["GO_Biological_Process_2021", "KEGG_2021_Human"],
    "Term": ["GO term", "KEGG term"],
    "Overlap": ["2/10", "3/10"],
    "Adjusted P-value": [0.01, 0.02],
    "Genes": ["A;B", "A;B;C"],
})
render_functional_context(frame, analysis_ran=True)
''')
        self.assertEqual(len(app.exception), 0)
        self.assertEqual([tab.label for tab in app.tabs], ["GO", "KEGG"])
        self.assertEqual(len(app.dataframe), 2)
        self.assertEqual(app.dataframe[0].value.iloc[0]["Terim"], "GO term")
        self.assertEqual(app.dataframe[1].value.iloc[0]["Terim"], "KEGG term")

    def test_significant_fdr_result_keeps_source_values(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_statistical_support
frame = pd.DataFrame({
    "gene": ["A"], "empirical_p": [0.0012345], "q_value": [0.004321],
    "significant_redistribution": [True], "Selection_Reason": ["null_fdr"],
})
render_statistical_support(frame, columns=list(frame.columns))
''')
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(len(app.dataframe), 1)
        shown = app.dataframe[0].value
        self.assertEqual(shown.iloc[0]["empirical_p"], 0.0012345)
        self.assertEqual(shown.iloc[0]["q_value"], 0.004321)
        self.assertTrue(bool(shown.iloc[0]["significant_redistribution"]))

    def test_old_report_with_only_identifier_is_safe(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_network_response
old_report = pd.DataFrame({"gene": ["OLD_ID"]})
render_network_response(targets=old_report, losses=None, redistribution=None)
''')
        self.assertEqual(len(app.exception), 0)
        self.assertGreaterEqual(len(app.dataframe), 2)

    def test_loss_view_explains_missing_signed_result_without_claiming_missing_data_file(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_network_response
render_network_response(
    targets=pd.DataFrame({"gene": ["TARGET"]}), losses=None, redistribution=None,
    loss_empty_message="Bu oturumda iki-yönlü signed sonuç bellekte yok.",
)
''')
        self.assertEqual(len(app.exception), 0)
        self.assertTrue(any("signed sonuç bellekte yok" in item.value for item in app.info))

    def test_network_response_matrices_and_source_tables_offer_independent_exports(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_network_response
frame = pd.DataFrame({
    "gene": ["A"], "Symbol": ["GENEA"], "Delta_PageRank_Pct": [-2.0],
    "Hinterland_Skoru": [10.0], "BC_Skoru": [0.5], "Gümrük_Kapisi": [True],
    "MyGene Adı": ["GENEA"], "GO Biyolojik Süreç": [["signal transduction"]],
})
render_network_response(targets=frame, losses=frame, redistribution=frame)
''')
        self.assertEqual(len(app.exception), 0)
        # Overview, targets, loss matrix/source, and redistribution matrix/source
        # each provide a local CSV/Excel choice and an independent download action.
        self.assertGreaterEqual(len(app.radio), 6)
        self.assertGreaterEqual(len(app.get("download_button")), 6)

    def test_empty_gene_detail_groups_are_not_rendered_as_errors(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.network_response import render_network_response
target = pd.DataFrame({"gene": ["TARGET"], "Hinterland_Skoru": [10.0]})
render_network_response(targets=target, losses=None, redistribution=None)
''')
        self.assertEqual(len(app.exception), 0)
        messages = [item.value for item in app.info]
        self.assertFalse(any("bu sonuçta gösterilecek mevcut alan yok" in message for message in messages))

    def test_multiple_unnamed_downloadable_tables_have_unique_keys(self) -> None:
        app = self._run('''\
import pandas as pd
from src.ui.tables import render_dataframe
render_dataframe(pd.DataFrame({"A": [1]}))
render_dataframe(pd.DataFrame({"B": [2]}))
render_dataframe(pd.DataFrame({"C": [3]}))
''')
        self.assertEqual(len(app.exception), 0)
        self.assertEqual(len(app.radio), 3)
        self.assertEqual(len(app.dataframe), 3)


class MetricPresentationTests(unittest.TestCase):
    def test_shared_labels_formats_sign_and_missing_behavior(self) -> None:
        self.assertEqual(metric_label("Delta_PageRank_Pct"), "Ağ önemi değişimi (%)")
        self.assertEqual(metric_label("Hinterland_Skoru"), "Ağdaki Önemi")
        self.assertEqual(metric_label("BC_Skoru"), "Geçiş Merkeziliği")
        self.assertEqual(metric_printf("Delta_PageRank_Pct"), "%+.4f%%")
        self.assertEqual(format_metric("Delta_PageRank_Pct", -2.84), "-2.8400%")
        self.assertEqual(format_metric("Delta_PageRank_Pct", 2.84), "+2.8400%")
        self.assertEqual(format_metric("Hinterland_Skoru", 64.86), "64.86")
        self.assertEqual(format_metric("BC_Skoru", 26252), "26252.000")
        self.assertEqual(format_metric("empirical_p", 0.001), "0.001000")
        self.assertEqual(format_metric("q_value", 0.02), "0.020000")
        self.assertEqual(format_metric("Gümrük_Kapisi", True), "✓")
        self.assertEqual(format_metric("Drug_Score", np.nan), "—")


class ExportAndHierarchyRegressionTests(unittest.TestCase):
    def test_csv_and_xlsx_roundtrip_keep_all_rows_columns_and_values(self) -> None:
        source = pd.DataFrame({
            "gene": ["A", "B"],
            "Delta_PageRank_Pct": [-2.84, 1.25],
            "Hinterland_Skoru": [64.86, 11.5],
            "BC_Skoru": [26252.0, 3.0],
            "Gümrük_Kapisi": [True, False],
            "annotation": [["x", "y"], np.nan],
        })
        before = source.copy(deep=True)
        csv_payload, _, csv_ext = export_bytes(source, "CSV (.csv)")
        csv_roundtrip = pd.read_csv(io.BytesIO(csv_payload))
        self.assertEqual(csv_ext, "csv")
        self.assertEqual(csv_roundtrip.shape, source.shape)
        self.assertEqual(list(csv_roundtrip.columns), list(source.columns))

        xlsx_payload, _, xlsx_ext = export_bytes(source, "Excel (.xlsx)")
        xlsx_roundtrip = pd.read_excel(io.BytesIO(xlsx_payload))
        self.assertEqual(xlsx_ext, "xlsx")
        self.assertEqual(xlsx_roundtrip.shape, source.shape)
        self.assertEqual(list(xlsx_roundtrip.columns), list(source.columns))
        pd.testing.assert_frame_equal(source, before)

    def test_result_surface_order_matches_the_desktop_reading_hierarchy(self) -> None:
        app_source = Path("app.py").read_text(encoding="utf-8")
        declarations = [
            "_cards_surface = st.container()",
            "_system_response_surface = st.container()",
            "_biological_context_surface = st.container()",
            "_evidence_surface = st.container()",
            "_functional_surface = st.container()",
            "_deep_interpretation_surface = st.container()",
            "_data_audit_surface = st.container()",
            "_advanced_surface = st.container()",
        ]
        positions = [app_source.index(declaration) for declaration in declarations]
        self.assertEqual(positions, sorted(positions))
        for label in (
            "Ağ yanıtı", "Biyolojik bağlam", "İSTATİSTİKSEL DESTEK",
            "Biyolojik yorum", "Yorum", "DATA & AUDIT", "İleri araçlar",
        ):
            self.assertIn(label, app_source + Path("src/ui/network_response.py").read_text(encoding="utf-8"))


class ResultStateSafetyTests(unittest.TestCase):
    def test_result_context_requires_exact_target_tissue_engine_and_parameters(self) -> None:
        context = app_state.build_simulation_result_context(
            species="İnsan", tissue="Pancreas", engine="Classic",
            targets=("ENSP_GLP1R",), localization="Mitochondrion",
            block_strength=.001, damping=.85, candidate_limit=20,
            tissue_normalization_mode="within_tissue",
        )
        current = dict(
            species="İnsan", tissue="Pancreas", engine="Classic",
            targets=("ENSP_GLP1R",), localization="Mitochondrion",
            block_strength=.001, damping=.85, candidate_limit=20,
            tissue_normalization_mode="within_tissue",
        )
        self.assertTrue(app_state.simulation_result_context_matches(context, **current))
        for field, value in (
            ("tissue", "Liver"),
            ("engine", "Evidence"),
            ("targets", ("ENSP_OTHER",)),
            ("damping", .90),
        ):
            changed = {**current, field: value}
            self.assertFalse(app_state.simulation_result_context_matches(context, **changed))

    def test_disk_fallback_result_is_bound_to_the_same_run_context(self) -> None:
        context_inputs = dict(
            species="İnsan", tissue="Lung", engine="Classic",
            targets=("ENSP_CFTR",), localization="Mitochondrion",
            block_strength=.001, damping=.85, candidate_limit=20,
            tissue_normalization_mode="within_tissue", bc_sample_sources=1200,
        )
        report = pd.DataFrame({"gene": ["ENSP_CFTR"], "Hasar_Tipi": ["Birincil"]})
        session = {}

        with patch.object(app_state.st, "session_state", session):
            context = app_state.build_simulation_result_context(**context_inputs)
            app_state.bind_simulation_result(
                report=report, context=context, ghosts=["ENSP_GHOST"], target_label="CFTR",
            )

        self.assertIs(session["rapor_df"], report)
        self.assertTrue(session["sim_tamam"])
        self.assertEqual(session["sim_hedef"], "CFTR")
        self.assertEqual(session["hayalet_genler"], ["ENSP_GHOST"])
        self.assertTrue(app_state.simulation_result_context_matches(
            session["simulation_result_context"], **context_inputs,
        ))

    def test_new_simulation_clears_presentation_result_objects_only(self) -> None:
        session = {
            "rapor_df": pd.DataFrame({"gene": ["OLD"]}),
            "fonksiyon_sonuc": object(),
            "hayalet_genler": ["OLD"],
            "_mygene_query_hash": "old",
            "enrichment_results": pd.DataFrame({"Term": ["old"]}),
            "enrichment_notices": ["old"],
            "signed_redistribution_result": pd.DataFrame({"gene": ["OLD"]}),
            "dose_df": pd.DataFrame({"old": [1]}),
            "sweep_df": pd.DataFrame({"old": [1]}),
            "block_strength": .001,
        }
        with patch.object(app_state.st, "session_state", session):
            app_state.reset_simulation_result()
        for key in (
            "rapor_df", "fonksiyon_sonuc", "hayalet_genler", "_mygene_query_hash",
            "enrichment_results", "signed_redistribution_result",
        ):
            self.assertIsNone(session[key])
        self.assertEqual(session["enrichment_notices"], [])
        self.assertNotIn("dose_df", session)
        self.assertNotIn("sweep_df", session)
        self.assertEqual(session["block_strength"], .001)

    def test_hard_reset_purges_network_presentation_state(self) -> None:
        fake_streamlit = Mock()
        fake_streamlit.session_state = {
            "rapor_df": pd.DataFrame({"gene": ["OLD"]}),
            "signed_redistribution_result": pd.DataFrame({"gene": ["OLD"]}),
            "enrichment_results": pd.DataFrame({"Term": ["old"]}),
            "enrichment_notices": ["old"],
            "sidebar_choice": "preserved by existing reset contract",
        }
        with patch.object(cache_manager, "st", fake_streamlit), patch.object(cache_manager, "force_gc", return_value=0):
            purged, collected = cache_manager.scorched_earth_reset()
        self.assertEqual(purged, 4)
        self.assertEqual(collected, 0)
        self.assertEqual(fake_streamlit.session_state, {"sidebar_choice": "preserved by existing reset contract"})


if __name__ == "__main__":
    unittest.main()
