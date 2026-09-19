"""Focused tests for the pure Research Context semantic adapter."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
import unittest

import numpy as np
import pandas as pd

from src.research.semantic_adapter import (
    FunctionalContextKind,
    ResponseDirection,
    ResultRole,
    ResultSource,
    adapt_snapshot,
)
from src.research.snapshot import build_snapshot


class SemanticAdapterTests(unittest.TestCase):
    @staticmethod
    def _snapshot(*, selection_mode: str = "null_fdr", signed=True, enrichment=True):
        report = pd.DataFrame({
            "gene": ["TARGET", "GAIN_B", "LOSS_REPORT", "GAIN_A", "EXPLORATORY"],
            "Symbol": ["T", "B", "L", "A", "E"],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", *("Üçüncül_Hasar_Stresli",) * 4],
            "Response_Direction": [None, "Influence Gain", "Influence Loss", "Influence Gain", "Influence Gain"],
            "Delta_PageRank_Pct": [None, 1.0, -8.0, 9.0, 0.5],
            "PageRank_Baseline": [0.2, 0.1, 0.3, 0.4, 0.5],
            "Hinterland_Skoru": [50.0, 20.0, 30.0, 40.0, 10.0],
            "BC_Skoru": [5.0, 2.0, 3.0, 4.0, 1.0],
            "Gümrük_Kapisi": [True, False, True, False, False],
            "q_value": [None, 0.001, 0.002, 0.9, 0.0001],
            "significant_redistribution": [False, True, True, False, "False"],
            "Systemic_Network_Shift_Pct": [2.5] * 5,
            "Lokalizasyon": [{"term": ["membrane"]}, "nucleus", "cytosol", "membrane", None],
        })
        signed_frame = pd.DataFrame({
            "gene": ["LOSS_SHALLOW", "GAIN_SIGNED", "LOSS_DEEP", "ZERO"],
            "Response_Direction": ["Influence Loss", "Influence Gain", "Influence Loss", "No Change"],
            "Delta_PageRank_Pct": [-1.0, 7.0, -20.0, 0.0],
            "q_value": [0.9, 0.01, 0.001, 0.001],
            "significant_redistribution": [False, True, True, True],
        }) if signed else None
        enrichment_frame = pd.DataFrame({
            "Pathway_Family": ["KEGG", "GO:BP", "CUSTOM", "GO"],
            "Kaynak": ["KEGG_2021_Human", "GO_Biological_Process_2021", "KEGG_fake", "GO_Legacy"],
            "Term": ["Pathway one", "GO later", "GO misleading term", "GO last"],
            "Adjusted P-value": [0.01, 0.02, 0.03, 0.04],
            "Genes": ["A;B", "B;C", "C", "D;E"],
        }) if enrichment else None
        return build_snapshot(
            report, signed_frame, enrichment_frame,
            species="Homo sapiens", taxon_id=9606, tissue="Pancreas",
            targets=("TARGET", "MISSING_TARGET"), attenuation=0.01,
            selection_mode=selection_mode, test_limit=200, tested_count=200,
            returned_count=5, top_n=20,
            threshold_parameters={"threshold": 0.05, "fdr_alpha": 0.05},
            created_at="2026-01-01T00:00:00Z",
        )

    def test_roles_keep_source_order_and_do_not_rerank(self) -> None:
        view = adapt_snapshot(self._snapshot())
        self.assertEqual(view.target_ids, ("TARGET", "MISSING_TARGET"))
        self.assertEqual([row.canonical_id for row in view.target_records], ["TARGET"])
        self.assertEqual(
            [row.canonical_id for row in view.positive_redistribution],
            ["GAIN_B", "GAIN_A", "EXPLORATORY"],
        )
        self.assertEqual(
            [row.canonical_id for row in view.network_losses],
            ["LOSS_SHALLOW", "LOSS_DEEP"],
        )
        self.assertEqual(
            [row.source_position for row in view.network_losses], [0, 2],
        )

    def test_fdr_uses_stored_flag_report_order_and_null_fdr_scope_only(self) -> None:
        view = adapt_snapshot(self._snapshot())
        self.assertEqual([row.canonical_id for row in view.fdr_supported], ["GAIN_B", "LOSS_REPORT"])
        self.assertEqual(
            view.fdr_supported[0].roles,
            (ResultRole.POSITIVE_REDISTRIBUTION, ResultRole.FDR_SUPPORTED),
        )
        self.assertEqual(
            view.fdr_supported[1].roles,
            (ResultRole.NETWORK_LOSS, ResultRole.FDR_SUPPORTED),
        )
        top_n = adapt_snapshot(self._snapshot(selection_mode="top_n"))
        self.assertEqual(top_n.fdr_supported, ())
        self.assertTrue(top_n.network_losses)
        self.assertTrue(all(
            ResultRole.FDR_SUPPORTED not in row.roles
            for row in top_n.network_losses
        ))

    def test_fdr_is_never_recomputed_from_q_value(self) -> None:
        report = pd.DataFrame({
            "gene": ["TARGET", "LOW_Q"],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", "Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [None, 2.0],
            "q_value": [None, 0.000001],
        })
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=("TARGET",), attenuation=None, selection_mode="null_fdr",
            test_limit=None, tested_count=2, returned_count=2,
        )
        view = adapt_snapshot(snapshot)
        self.assertEqual(view.fdr_supported, ())
        self.assertTrue(any("q-values were not reinterpreted" in notice for notice in view.notices))

    def test_functional_rows_are_authoritatively_partitioned_in_source_order(self) -> None:
        view = adapt_snapshot(self._snapshot())
        self.assertEqual([row.term for row in view.go_terms], ["GO later", "GO last"])
        self.assertEqual([row.term for row in view.kegg_pathways], ["Pathway one"])
        self.assertEqual([row.term for row in view.unclassified_enrichment], ["GO misleading term"])
        self.assertEqual(
            [row.term for row in view.functional_rows],
            ["Pathway one", "GO later", "GO misleading term", "GO last"],
        )
        self.assertEqual(view.go_terms[0].contract.metric_id, "go_enrichment")
        self.assertEqual(view.kegg_pathways[0].contract.metric_id, "kegg_enrichment")
        self.assertEqual(view.go_terms[0].members, ("B", "C"))
        self.assertIs(view.unclassified_enrichment[0].kind, FunctionalContextKind.UNCLASSIFIED)

    def test_metrics_are_raw_values_bound_to_existing_contracts(self) -> None:
        row = adapt_snapshot(self._snapshot()).positive_redistribution[0]
        self.assertEqual(row.metric("delta_pagerank_pct").raw_value, 1.0)
        self.assertEqual(row.metric("delta_pagerank_pct").contract.display_name, "Ağ önemi değişimi (%)")
        self.assertEqual(row.metric("hinterland").source_field, "Hinterland_Skoru")
        self.assertEqual(row.metric("systemic_network_shift_pct").raw_value, 2.5)
        self.assertIsNone(row.metric("global_efficiency_change_pct"))

    def test_snapshot_scope_is_preserved_without_row_count_repair(self) -> None:
        snapshot = self._snapshot()
        view = adapt_snapshot(snapshot)
        self.assertEqual(view.scope, snapshot.scope)
        self.assertEqual(view.scope.tested_count, 200)
        self.assertEqual(view.scope.returned_count, 5)
        self.assertEqual(view.scope.top_n, 20)
        self.assertEqual(view.scope.threshold, 0.05)

    def test_source_and_nested_values_are_immutable_and_detached(self) -> None:
        snapshot = self._snapshot()
        before = snapshot.report.records()
        view = adapt_snapshot(snapshot)
        row = view.target_records[0]
        with self.assertRaises(TypeError):
            row.values["gene"] = "MUTATED"  # type: ignore[index]
        with self.assertRaises(TypeError):
            row.values["Lokalizasyon"]["term"] = ()  # type: ignore[index]
        with self.assertRaises(FrozenInstanceError):
            row.direction = ResponseDirection.LOSS  # type: ignore[misc]
        self.assertEqual(snapshot.report.records(), before)
        self.assertEqual(adapt_snapshot(snapshot), view)

    def test_missing_sources_and_empty_legacy_results_are_safe(self) -> None:
        report = pd.DataFrame({"gene": ["TARGET", "LEGACY_B", "LEGACY_A"], "Delta_PageRank_Pct": [np.nan, 0.5, 2.0]})
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=("TARGET",), attenuation=None, selection_mode="top_n",
            test_limit=None, tested_count=None, returned_count=3,
        )
        view = adapt_snapshot(snapshot)
        self.assertEqual([row.canonical_id for row in view.positive_redistribution], ["LEGACY_B", "LEGACY_A"])
        self.assertEqual(view.network_losses, ())
        self.assertEqual(view.go_terms, ())
        self.assertEqual(view.kegg_pathways, ())
        self.assertTrue(any("not inferred from the report" in notice for notice in view.notices))

    def test_explicit_direction_wins_over_contradictory_legacy_delta(self) -> None:
        signed = pd.DataFrame({
            "gene": ["NOT_A_LOSS", "REAL_LOSS"],
            "Response_Direction": ["Influence Gain", "Influence Loss"],
            "Delta_PageRank_Pct": [-99.0, 99.0],
        })
        snapshot = build_snapshot(
            pd.DataFrame({"gene": ["TARGET"]}), signed, pd.DataFrame(),
            species="human", taxon_id=9606, tissue="None", targets=("TARGET",),
            attenuation=None, selection_mode="top_n", test_limit=None,
            tested_count=2, returned_count=1,
        )
        view = adapt_snapshot(snapshot)
        self.assertEqual([row.canonical_id for row in view.network_losses], ["REAL_LOSS"])
        self.assertIs(view.network_losses[0].direction, ResponseDirection.LOSS)
        self.assertIs(view.network_losses[0].source, ResultSource.SIGNED_REDISTRIBUTION)

    def test_blank_identity_is_not_guessed_from_symbol(self) -> None:
        report = pd.DataFrame({
            "gene": ["TARGET", ""], "Symbol": ["T", "SYMBOL_ONLY"],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", "Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [None, 1.0],
        })
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=("TARGET",), attenuation=None, selection_mode="top_n",
            test_limit=None, tested_count=2, returned_count=2,
        )
        row = adapt_snapshot(snapshot).positive_redistribution[0]
        self.assertIsNone(row.canonical_id)
        self.assertEqual(row.symbol, "SYMBOL_ONLY")


if __name__ == "__main__":
    unittest.main()
