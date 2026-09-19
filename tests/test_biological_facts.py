"""Tests for Phase 2 facts and data-driven theme ontology."""

from __future__ import annotations

import unittest

import pandas as pd

from src.interpretation.biological_facts import extract_biological_facts, extract_candidate_facts
from src.interpretation.theme_detector_v2 import detect_themes_from_facts


class BiologicalFactsTests(unittest.TestCase):
    def _facts_for(self, rows: list[dict[str, object]]) -> list[dict[str, object]]:
        frame = pd.DataFrame(rows)
        before = frame.copy(deep=True)
        facts = extract_biological_facts(frame)
        pd.testing.assert_frame_equal(frame, before)
        return facts["candidates"]

    def assert_theme(self, rows: list[dict[str, object]], label: str) -> None:
        themes = detect_themes_from_facts(self._facts_for(rows))
        self.assertTrue(any(item["theme"] == label for item in themes), themes)

    def test_cftr_context_uses_annotation_evidence(self) -> None:
        rows = [
            {"Symbol": "ANO2", "Delta_PageRank_Pct": 2.2, "GO Biyolojik Süreç": "chloride transport; anion channel activity"},
            {"Symbol": "SLC4A3", "Delta_PageRank_Pct": 1.6, "GO Biyolojik Süreç": "bicarbonate transport; ion homeostasis"},
            {"Symbol": "AQP8", "Delta_PageRank_Pct": 1.1, "GO Biyolojik Süreç": "water transport"},
        ]
        self.assert_theme(rows, "Klorür / anyon taşınması")
        self.assert_theme(rows, "Bikarbonat taşınması")
        self.assert_theme(rows, "Su taşınması")

    def test_ryr2_context_requires_calcium_annotations(self) -> None:
        rows = [
            {"Symbol": "RYR2", "GO Biyolojik Süreç": "calcium ion release from sarcoplasmic reticulum; cardiac muscle contraction"},
            {"Symbol": "CASQ2", "GO Biyolojik Süreç": "calcium ion homeostasis; sarcoplasmic reticulum"},
        ]
        self.assert_theme(rows, "Kalsiyum işlenmesi")
        self.assert_theme(rows, "Uyarılma-kasılma eşleşmesi")
        self.assert_theme(rows, "Sarkoplazmik retikulum kalsiyum bağlamı")

    def test_ldlr_context_requires_lipid_annotations(self) -> None:
        rows = [
            {"Symbol": "LDLR", "GO Biyolojik Süreç": "cholesterol homeostasis; low-density lipoprotein receptor activity"},
            {"Symbol": "APOB", "GO Biyolojik Süreç": "lipoprotein particle clearance; cholesterol transport"},
        ]
        self.assert_theme(rows, "Kolesterol işlenmesi")
        self.assert_theme(rows, "Lipoprotein işlenmesi")

    def test_empty_or_missing_annotation_is_safe(self) -> None:
        facts = extract_candidate_facts(pd.DataFrame({"Symbol": ["GENE1"], "GO Biyolojik Süreç": [None]}))
        self.assertEqual(facts[0].annotations["go_biological_process"], ())
        self.assertEqual(detect_themes_from_facts([fact.to_dict() for fact in facts]), [])

    def test_mouse_fallback_is_explicit(self) -> None:
        result = extract_biological_facts(pd.DataFrame({"Symbol": ["Cftr"]}), organism="Mus musculus")
        self.assertFalse(result["family_reference_available"])
        self.assertIn("fare", str(result["family_reference_warning"]).casefold())


if __name__ == "__main__":
    unittest.main()
