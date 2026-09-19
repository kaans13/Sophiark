from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import igraph as ig
import pandas as pd


def tiny_inputs(prefix: str):
    names = [f"{prefix}P{i}" for i in range(5)]
    graph = ig.Graph(edges=[(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)], directed=False)
    graph.vs["name"] = names
    graph.es["weight"] = [1.0] * graph.ecount()
    graph.es["distance"] = [1.0] * graph.ecount()
    frame = pd.DataFrame({
        "gene": names,
        "Symbol": [f"S{i}" for i in range(5)],
        "PageRank": [0.25, 0.24, 0.20, 0.17, 0.14],
        "BC_Skoru": [0.0, 3.0, 4.0, 3.0, 0.0],
        "Hinterland_Skoru": [100, 80, 60, 40, 20],
        "Topluluk_ID": [1] * 5,
        "Lokalizasyon": ["Nucleus"] * 5,
        "Gümrük_Kapisi": [False] * 5,
    })
    return graph, frame, names[1]


class ForwardSmokeTests(unittest.TestCase):
    def _run(self, module, prefix):
        item, frame, target = tiny_inputs(prefix)
        old_bonus, old_penalty = module.apply_biological_bonus, module.apply_functional_penalty
        module.apply_biological_bonus = lambda *args, **kwargs: None
        module.apply_functional_penalty = lambda *args, **kwargs: None
        try:
            with tempfile.TemporaryDirectory() as directory:
                result = module.run_infection_simulation(
                    item, frame, Path(directory), spesifik_hedefler=[target],
                    redistribution_mode="top_n", exploratory_top_n=4,
                    null_iterations=0, efficiency_sample_sources=None,
                    local_efficiency_sample_nodes=None, tissue="Lung",
                )
                self.assertFalse(result.empty)
                for column in (
                    "Systemic_Network_Shift_Pct", "Global_Efficiency_Baseline",
                    "Global_Efficiency_Perturbed", "Mean_Local_Efficiency_Baseline",
                    "Target_Local_Efficiency_Change_Pct", "Response_Direction",
                    "Tissue", "Tissue_Normalization_Mode",
                ):
                    self.assertIn(column, result.columns)
                self.assertEqual(set(result["Tissue"]), {"Lung"})
        finally:
            module.apply_biological_bonus, module.apply_functional_penalty = old_bonus, old_penalty

    def test_human_forward(self):
        from src import biology_logic
        self._run(biology_logic, "ENSP")

    def test_mouse_forward(self):
        from src.mouse import biology_logic
        self._run(biology_logic, "ENSMUSP")

    def test_multi_target_predicted_network_interaction(self):
        from src import biology_logic
        from src.services.network_interaction import predicted_network_interaction
        item, frame, _ = tiny_inputs("ENSP")
        targets = ["ENSPP1", "ENSPP3"]
        old_bonus, old_penalty = biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty
        biology_logic.apply_biological_bonus = lambda *args, **kwargs: None
        biology_logic.apply_functional_penalty = lambda *args, **kwargs: None
        try:
            with tempfile.TemporaryDirectory() as directory:
                result = predicted_network_interaction(
                    motor_module=biology_logic, graph=item, scores=frame,
                    targets=targets, output_root=Path(directory),
                    block_strength=.1, damping=.85, relative_tolerance=.1,
                )
            self.assertEqual(len(result), 2)
            self.assertIn("Beklenen Additive Etki", result)
            self.assertIn("Gözlenen Kombine Etki", result)
            self.assertTrue(result["Yorum Sınırı"].str.contains("biyolojik synergy iddiası değildir").all())
        finally:
            biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty = old_bonus, old_penalty

    def test_severity_sensitivity_runs_every_declared_scenario(self):
        from src import biology_logic
        from src.scientific.perturbation import SEVERITY_FRACTIONS
        from src.services.severity_sensitivity import run_severity_sensitivity

        item, frame, target = tiny_inputs("ENSP")
        old_bonus, old_penalty = biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty
        biology_logic.apply_biological_bonus = lambda *args, **kwargs: None
        biology_logic.apply_functional_penalty = lambda *args, **kwargs: None
        try:
            with tempfile.TemporaryDirectory() as directory:
                result = run_severity_sensitivity(
                    motor_module=biology_logic, graph=item, scores=frame,
                    targets=[target], output_root=Path(directory), damping=.85,
                )
            self.assertEqual(result["Senaryo"].tolist(), list(SEVERITY_FRACTIONS))
            self.assertEqual(result["Kalan Edge Strength"].tolist(), list(SEVERITY_FRACTIONS.values()))
            self.assertTrue(result["Global Efficiency Değişimi (%)"].notna().all())
            self.assertTrue(result["Sistemik Ağ Kayması (%)"].notna().all())
        finally:
            biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty = old_bonus, old_penalty

    def test_target_stress_and_compensation_regression(self):
        from src import biology_logic
        from src.services.compensation import predicted_compensation_candidates
        from src.services.target_stress import preselect_candidates, validate_candidates
        item, frame, target = tiny_inputs("ENSP")
        mapping = dict(zip(frame["gene"], frame["Symbol"]))
        candidates = preselect_candidates(
            graph=item, scores=frame, target_gene=target,
            regulators={"S1": [("S0", "Activation")]}, gene_to_symbol=mapping,
            mode="Birleşik", limit=3,
        )
        ghost = lambda report, graph: (report, [])
        with tempfile.TemporaryDirectory() as directory:
            validated = validate_candidates(
                candidates=candidates, motor_module=biology_logic, graph=item,
                scores=frame, target_gene=target, output_root=Path(directory),
                block_strength=.1, damping=.85, ghost_filter=ghost,
                deep_validation_top_n=1,
            )
            compensation = predicted_compensation_candidates(
                motor_module=biology_logic, graph=item, scores=frame,
                perturbed_gene=target, output_root=Path(directory),
                block_strength=.1, damping=.85, ghost_filter=ghost,
            )
        self.assertFalse(validated.empty)
        self.assertEqual(validated.iloc[0]["BC Doğrulama Modu"], "Full Validation")
        self.assertIn("Full ΔBC", validated)
        self.assertFalse(compensation.empty)


if __name__ == "__main__":
    unittest.main()
