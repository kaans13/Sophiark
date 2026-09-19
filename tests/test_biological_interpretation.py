from __future__ import annotations

import unittest

import pandas as pd

from src.interpretation import interpret_forward_results


class BiologicalInterpretationTests(unittest.TestCase):
    def _report(self) -> pd.DataFrame:
        return pd.DataFrame({
            "Systemic_Network_Shift_Pct": [0.1275] * 4,
            "Global_Efficiency_Change_Pct": [-0.02] * 4,
            "Mean_Local_Efficiency_Change_Pct": [-0.10] * 4,
            "Target_Local_Efficiency_Change_Pct": [-0.70] * 4,
            "Redistribution_Selection_Mode": ["null_fdr"] * 4,
        })

    def _candidates(self) -> pd.DataFrame:
        return pd.DataFrame({
            "Symbol": ["ANO2", "SLC4A3", "LRRC8B", "LRRC8C"],
            "Delta_PageRank_Pct": [2.216, 1.605, 1.2, 1.1],
            "Hinterland_Skoru": [33.0, 22.0, 12.0, 16.0],
            "BC_Skoru": [2.0, 1.0, 1.5, 1.2],
            "Lokalizasyon": ["Membrane"] * 4,
            "GO Biyolojik Süreç": [
                "chloride transport", "bicarbonate transport", "cell volume regulation", "anion transport",
            ],
            "q_value": [.01, .02, .03, .04],
            "empirical_p": [.01, .02, .03, .04],
            "significant_redistribution": [True] * 4,
        })

    def test_interpretation_is_read_only_and_uses_safe_language(self):
        report = self._report()
        candidates = self._candidates()
        report_before = report.copy(deep=True)
        candidates_before = candidates.copy(deep=True)
        result = interpret_forward_results(
            report=report, candidates=candidates, targets=["CFTR_ID"], gene_to_symbol={"CFTR_ID": "CFTR"},
            tissue="Lung", organism="Homo sapiens", edge_attenuation=.001, network_nodes=16955, network_edges=513176,
            directed_evidence=pd.DataFrame([{
                "source": "NFKB1", "target": "ANO2", "consensus_direction": True,
                "consensus_stimulation": True, "consensus_inhibition": False, "provenance": "TRRUST",
            }]),
            enrichment=pd.DataFrame([{"Term": "Ion transport", "Adjusted P-value": .01, "Genes": "ANO2;LRRC8B"}]),
        )
        pd.testing.assert_frame_equal(report, report_before)
        pd.testing.assert_frame_equal(candidates, candidates_before)
        self.assertEqual(result.top_candidates[0]["gene"], "ANO2")
        self.assertIn("göreli", result.top_candidates[0]["interpretation"].casefold())
        self.assertTrue(any(item["family"] == "LRRC8" for item in result.family_patterns))
        self.assertTrue(any(item["theme"] == "Klorür / anyon taşınması" for item in result.functional_themes))
        self.assertIsNotNone(result.statistical_validation)
        self.assertIn("Grafik verimliliği", result.system_response)
        self.assertNotIn("activated", " ".join(item["interpretation"] for item in result.top_candidates).casefold())

    def test_top_n_suppresses_statistical_validation(self):
        report = self._report().assign(Redistribution_Selection_Mode="top_n")
        result = interpret_forward_results(
            report=report, candidates=self._candidates(), targets=[], gene_to_symbol={}, tissue="Lung", organism="Homo sapiens",
        )
        self.assertIsNone(result.statistical_validation)

    def test_candidate_identity_keeps_canonical_id_separate_from_symbol(self):
        candidates = self._candidates().copy()
        candidates.insert(0, "gene", ["ENSP00000355629", "ENSP00000497669", "ENSP00000300000", "ENSP00000300001"])
        candidates.loc[1, "Symbol"] = None
        original = candidates.copy(deep=True)

        result = interpret_forward_results(
            report=self._report(), candidates=candidates, targets=[], gene_to_symbol={}, tissue="Lung", organism="Homo sapiens",
        )

        self.assertEqual(result.top_candidates[0]["gene"], "ANO2")
        self.assertEqual(result.top_candidates[0]["entity_id"], "ENSP00000355629")
        missing_symbol = next(item for item in result.top_candidates if item["entity_id"] == "ENSP00000497669")
        self.assertEqual(missing_symbol["gene"], "—")
        pd.testing.assert_frame_equal(candidates, original)


if __name__ == "__main__":
    unittest.main()
