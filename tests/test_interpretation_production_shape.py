"""Production-shaped regression: raw report plus existing local annotation cache."""

from __future__ import annotations

from pathlib import Path
import unittest

import pandas as pd

from src.interpretation import interpret_forward_results


ROOT = Path(__file__).resolve().parents[1]
REPORT_PATH = ROOT / "outputs" / "enfeksiyon_sok_dalgasi_raporu.csv"
SYMBOL_PATH = ROOT / "data" / "processed" / "ensp_with_symbols.csv"
CACHE_PATH = ROOT / "data" / "processed" / "mygene_cache.pkl"


@unittest.skipUnless(REPORT_PATH.exists() and SYMBOL_PATH.exists() and CACHE_PATH.exists(), "Yerel production-shape fixture mevcut değil")
class ProductionShapeInterpretationTests(unittest.TestCase):
    def test_raw_forward_report_uses_local_annotation_fallback(self) -> None:
        report = pd.read_csv(REPORT_PATH)
        symbol_frame = pd.read_csv(SYMBOL_PATH, encoding="utf-8-sig")
        mapping = dict(zip(symbol_frame["gene"].astype(str), symbol_frame["Symbol"].astype(str)))
        candidates = report.iloc[1:].copy(deep=True)
        candidates["Symbol"] = candidates["gene"].astype(str).map(mapping).fillna("")
        candidates = candidates.sort_values("Delta_PageRank_Pct", ascending=False, kind="mergesort")
        self.assertNotIn("GO Biyolojik Süreç", candidates.columns)  # production-shaped input is raw
        result = interpret_forward_results(
            report=report, candidates=candidates, targets=[str(report["Hasar_Hedefleri"].dropna().iloc[0])],
            gene_to_symbol=mapping, tissue="Lung", organism="Homo sapiens",
        )
        self.assertGreater(result.debug["annotation_source_counts"]["yerel_mygene_cache"], 0)
        self.assertTrue(result.functional_themes)
        self.assertIsNotNone(result.analysis_summary["primary_theme"])
        self.assertIn("pozitif topolojik yeniden dağılım", result.analysis_summary["primary"])
        self.assertTrue(result.top_candidates[0]["interpretation"])


if __name__ == "__main__":
    unittest.main()
