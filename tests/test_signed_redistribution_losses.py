"""Read-only tests for signed redistribution presentation; no graph simulation."""

from __future__ import annotations

import unittest
from types import SimpleNamespace

import pandas as pd
from src import biology_logic
from src.interpretation import interpret_forward_results, select_redistribution_losses
from src.services.analysis_service import execute_simulation
from src.services.engine_presentation import select_active_result_flow
from src.ui.tables import export_bytes


class SignedRedistributionLossTests(unittest.TestCase):
    def setUp(self) -> None:
        self.signed = pd.DataFrame({
            "gene": ["GAIN", "ATP2A2", "MID", "ZERO"],
            "Symbol": ["GAIN", "ATP2A2", "MID", "ZERO"],
            "Delta_PageRank_Pct": [2.2, -5.0, -0.1, 0.0],
            "Hinterland_Skoru": [10.0, 80.0, 12.0, 5.0],
            "BC_Skoru": [.1, .7, .2, .0],
            "GO Biyolojik Süreç": ["calcium transport", "calcium transport", "membrane transport", ""],
        })

    def test_sign_selection_excludes_positive_and_zero(self) -> None:
        losses = select_redistribution_losses(self.signed, min_abs_delta_pct=.05, limit=None)
        self.assertEqual(list(losses["gene"]), ["ATP2A2", "MID"])
        self.assertTrue((losses["Delta_PageRank_Pct"] < 0).all())

    def test_losses_are_sorted_most_negative_first(self) -> None:
        source = self.signed.assign(Delta_PageRank_Pct=[1.0, -.1, -5.0, -1.0])
        losses = select_redistribution_losses(source, min_abs_delta_pct=.05, limit=None)
        self.assertEqual(list(losses["Delta_PageRank_Pct"]), [-5.0, -1.0, -.1])

    def test_positive_input_is_not_mutated(self) -> None:
        positives = self.signed[self.signed["Delta_PageRank_Pct"] > .05].copy(deep=True)
        before = positives.copy(deep=True)
        select_redistribution_losses(self.signed, min_abs_delta_pct=.05)
        pd.testing.assert_frame_equal(positives, before)

    def test_negative_interpretation_and_canonical_lookup_are_safe(self) -> None:
        gains = self.signed[self.signed["Delta_PageRank_Pct"] > .05]
        losses = select_redistribution_losses(self.signed, min_abs_delta_pct=.05)
        interpretation = interpret_forward_results(
            report=pd.DataFrame({"Redistribution_Selection_Mode": ["top_n"]}), candidates=gains,
            loss_candidates=losses, targets=["RYR2"], gene_to_symbol={"RYR2": "RYR2"},
            tissue="Heart Muscle", organism="Homo sapiens", canonical_components=["ATP2A2"],
        )
        self.assertEqual(interpretation.redistribution_balance["strongest_loss"]["gene"], "ATP2A2")
        self.assertEqual(interpretation.canonical_component_check[0]["side"], "negatif")
        text = interpretation.loss_candidates[0]["interpretation"].casefold()
        for forbidden in ("downregulated", "inhibited", "suppressed", "expression decreased"):
            self.assertNotIn(forbidden, text)

    def test_separate_loss_export_does_not_replace_existing_export(self) -> None:
        gains = self.signed[self.signed["Delta_PageRank_Pct"] > .05]
        losses = select_redistribution_losses(self.signed, min_abs_delta_pct=.05)
        gain_payload, _, gain_ext = export_bytes(gains, "CSV (.csv)")
        loss_payload, _, loss_ext = export_bytes(losses, "CSV (.csv)")
        self.assertEqual(gain_ext, "csv")
        self.assertEqual(loss_ext, "csv")
        self.assertIn(b"GAIN", gain_payload)
        self.assertIn(b"ATP2A2", loss_payload)

    def test_ui_handoff_keeps_signed_frame_when_presentation_transforms_copy_report(self) -> None:
        signed = self.signed
        class Graph:
            def copy(self):
                return self

        class Motor:
            def run_infection_simulation(self, *_args, **_kwargs):
                return pd.DataFrame({"gene": ["GAIN"]})

        previous = biology_logic._LATEST_SIGNED_REDISTRIBUTION
        handoff: dict[str, pd.DataFrame] = {}
        try:
            biology_logic._LATEST_SIGNED_REDISTRIBUTION = signed
            execute_simulation(
                motor_module=Motor(), graph=Graph(), scores=pd.DataFrame(), output_dir="outputs",
                targets=["TARGET"], localization="", block_strength=.001, damping=.85,
                apply_druggability=lambda frame: frame.copy(deep=True),
                ghost_filter=lambda frame, _graph: (frame.copy(deep=True), []), scientific_options={},
                signed_result_sink=lambda frame: handoff.__setitem__("signed", frame),
            )
            losses = select_redistribution_losses(handoff.get("signed"), limit=None)
            self.assertEqual(list(losses["gene"]), ["ATP2A2", "MID"])
        finally:
            biology_logic._LATEST_SIGNED_REDISTRIBUTION = previous

    def test_evidence_handoff_keeps_evidence_signed_frame_and_active_flow(self) -> None:
        """Evidence presentation must not fall back to a Classic signed result."""
        from src.evidence import engine as evidence_engine

        evidence_signed = self.signed.assign(
            Non_Text_Evidence=[.9, .8, .7, .6],
            Text_Dependency=[.1, .2, .3, .4],
            Applied_Text_Modifier=[1.0] * 4,
            Physical_Supported_Edges=[1, 2, 0, 0],
        )

        class Graph:
            def copy(self):
                return self

        class EvidenceMotor:
            __name__ = "src.evidence.main"

            def run_infection_simulation(self, *_args, **_kwargs):
                return pd.DataFrame({"gene": ["GAIN"], "Hasar_Tipi": ["Üçüncül_Hasar_Stresli"]})

        previous = evidence_engine._LATEST_RESULT
        handoff: dict[str, pd.DataFrame] = {}
        try:
            evidence_engine._LATEST_RESULT = SimpleNamespace(signed_response=evidence_signed)
            report, _ = execute_simulation(
                motor_module=EvidenceMotor(), graph=Graph(), scores=pd.DataFrame(), output_dir="outputs",
                targets=["TARGET"], localization="", block_strength=.001, damping=.85,
                apply_druggability=lambda frame: frame.copy(deep=True),
                ghost_filter=lambda frame, _graph: (frame.copy(deep=True), []), scientific_options={},
                signed_result_sink=lambda frame: handoff.__setitem__("signed", frame),
            )
            active = select_active_result_flow(
                classic_report=report, classic_signed_response=handoff.get("signed"),
                requested_mode="Evidence", directed_bundle=None,
            )
            self.assertEqual(active.engine, "evidence")
            self.assertTrue(active.signed_response.equals(evidence_signed))
            self.assertIn("Non_Text_Evidence", active.signed_response.columns)
            self.assertEqual(list(select_redistribution_losses(active.signed_response)["gene"]), ["ATP2A2", "MID"])
        finally:
            evidence_engine._LATEST_RESULT = previous


if __name__ == "__main__":
    unittest.main()
