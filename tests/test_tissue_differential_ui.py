from __future__ import annotations

import unittest

import pandas as pd

from src.models import TissueDifferentialResult
from src.services.propagation_trace import _community_context
from src.ui.pages.tissue_differential import _narrative


class TissueDifferentialPresentationTests(unittest.TestCase):
    def _result(self, summary: pd.DataFrame) -> TissueDifferentialResult:
        return TissueDifferentialResult(summary=summary, critical_genes=pd.DataFrame(), notices=[], similarity=pd.DataFrame())

    def test_narrative_uses_canonical_local_efficiency_column(self):
        result = self._result(pd.DataFrame({
            "Doku": ["Lung", "Pancreas"],
            "Ortalama Local Efficiency Değişimi (%)": [0.12, 0.45],
        }))
        self.assertIn("Pancreas", _narrative(result))
        self.assertIn("0.450", _narrative(result))

    def test_narrative_handles_missing_canonical_column_without_keyerror(self):
        result = self._result(pd.DataFrame({"Doku": ["Lung"], "Sistemik Ağ Kayması (%)": [0.1275]}))
        self.assertIn("canonical yerel verimlilik metriği", _narrative(result))

    def test_community_context_is_always_string_for_arrow_serialization(self):
        self.assertEqual(_community_context(12), "12")
        self.assertEqual(_community_context("N/A"), "N/A")
        self.assertEqual(_community_context(None), "N/A")


if __name__ == "__main__":
    unittest.main()
