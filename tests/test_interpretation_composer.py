"""Safety and determinism tests for Turkish sentence-bank interpretation."""

from __future__ import annotations

import unittest

from src.interpretation.candidate_interpreter import compose_candidate_interpretation
from src.interpretation.sentence_bank import FORBIDDEN_ASSERTION_TOKENS


class InterpretationComposerTests(unittest.TestCase):
    def test_composer_is_deterministic_and_source_bounded(self) -> None:
        fact = {
            "symbol": "ANO2", "delta_pagerank_pct": 2.216,
            "hinterland_baseline": 33.0, "bc_baseline": 2.0,
            "annotations": {"go_biological_process": ("chloride transport",), "localization": ("membrane",)},
        }
        first = compose_candidate_interpretation(fact, theme_labels=["Klorür / anyon taşınması"])
        second = compose_candidate_interpretation(fact, theme_labels=["Klorür / anyon taşınması"])
        self.assertEqual(first, second)
        self.assertIn("+2.2160%", first["text"])
        self.assertIn("başlangıç topolojik bağlam", first["text"])
        self.assertEqual(len(first["validation_suggestions"]), 3)

    def test_forbidden_positive_assertions_are_absent(self) -> None:
        output = compose_candidate_interpretation({"symbol": "GENE1", "delta_pagerank_pct": 1.0, "annotations": {}})
        lowered = str(output).casefold()
        for token in FORBIDDEN_ASSERTION_TOKENS:
            self.assertNotIn(token.casefold(), lowered)

    def test_missing_schema_does_not_substitute_another_metric(self) -> None:
        output = compose_candidate_interpretation({"symbol": "GENE1", "annotations": {}})
        self.assertIn("mevcut değildir", output["text"])
        self.assertEqual(output["confidence_label"], "düşük")


if __name__ == "__main__":
    unittest.main()
