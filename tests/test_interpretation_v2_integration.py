"""End-to-end read-only contract for deterministic interpretation V2."""

from __future__ import annotations

import unittest

import pandas as pd

from src.interpretation import interpret_forward_results


class InterpretationV2IntegrationTests(unittest.TestCase):
    def test_cftr_lung_reference_fields_pass_through_unchanged(self) -> None:
        report = pd.DataFrame({
            "Systemic_Network_Shift_Pct": [0.1275],
            "Top_Positive_PageRank_Mean_Pct": [0.5586],
            "BC_Graf_Dugum_Sayisi": [16955], "BC_Graf_Kenar_Sayisi": [513176],
            "Benzersiz_Bloklanan_Kenar": [209], "Redistribution_Selection_Mode": ["top_n"],
        })
        candidates = pd.DataFrame({
            "Symbol": ["ANO2", "SLC4A3", "SLC4A5", "LRRC8B", "LRRC8C", "AQP8"],
            "Delta_PageRank_Pct": [2.216, 1.605, 1.585, 1.20, 1.10, 1.00],
            "Hinterland_Skoru": [33.0, 22.0, 21.0, 12.0, 16.0, 8.0],
            "BC_Skoru": [2.0, 1.0, 1.0, 1.5, 1.2, 0.8],
            "GO Biyolojik Süreç": ["chloride transport", "bicarbonate transport", "bicarbonate transport", "cell volume regulation", "anion transport", "water transport"],
        })
        report_before, candidates_before = report.copy(deep=True), candidates.copy(deep=True)
        result = interpret_forward_results(
            report=report, candidates=candidates, targets=["CFTR"], gene_to_symbol={"CFTR": "CFTR"},
            tissue="Lung", organism="Homo sapiens", edge_attenuation=.001,
            network_nodes=16955, network_edges=513176,
        )
        pd.testing.assert_frame_equal(report, report_before)
        pd.testing.assert_frame_equal(candidates, candidates_before)
        self.assertEqual(result.top_candidates[0]["gene"], "ANO2")
        self.assertEqual(result.statistical_validation, None)
        self.assertTrue(any(item["theme"] == "Klorür / anyon taşınması" for item in result.functional_themes))
        self.assertTrue(any({"LRRC8B", "LRRC8C"}.issubset(set(item["members"])) for item in result.family_clusters))
        self.assertIn("16,955", result.overview)
        self.assertIn("0.1275%", result.system_response)

    def test_mouse_uses_deterministic_degraded_family_context(self) -> None:
        report = pd.DataFrame({"Redistribution_Selection_Mode": ["top_n"]})
        candidates = pd.DataFrame({"Symbol": ["Cftr", "Ano2"], "Delta_PageRank_Pct": [1.0, .5]})
        result = interpret_forward_results(
            report=report, candidates=candidates, targets=["Cftr"], gene_to_symbol={"Cftr": "Cftr"},
            tissue="Lung", organism="Mus musculus",
        )
        self.assertEqual(result.family_clusters, [])
        self.assertTrue(result.top_candidates[0]["interpretation"])


if __name__ == "__main__":
    unittest.main()
