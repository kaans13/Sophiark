from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory
from unittest import TestCase
from unittest.mock import patch

import igraph as ig
import numpy as np
import pandas as pd

from src import biology_logic


class DoseResponseControlFlowTests(TestCase):
    def _graph(self) -> ig.Graph:
        graph = ig.Graph(n=2, edges=[(0, 1)], directed=False)
        graph.vs["name"] = ["ENSP00000003084", "ENSP00000447297"]
        graph.es["weight"] = [0.8]
        graph.es["distance"] = [1.25]
        return graph

    @staticmethod
    def _scores() -> pd.DataFrame:
        return pd.DataFrame({
            "gene": ["ENSP00000003084", "ENSP00000447297"],
            "Gümrük_Kapisi": [True, False], "Lokalizasyon": ["Membrane", "Cytosol"],
        })

    def test_existing_dose_loop_returns_one_structured_row_per_dose_without_graph_leak(self):
        graph = self._graph()
        original_weight = list(graph.es["weight"])
        original_distance = list(graph.es["distance"])
        observed_survival_fractions: list[float] = []
        original_hill = biology_logic._hill_inhibition

        def tracked_hill(weight, survival_fraction, hill_n):
            observed_survival_fractions.append(float(survival_fraction))
            return original_hill(weight, survival_fraction, hill_n)

        with TemporaryDirectory() as directory, \
             patch.object(biology_logic, "resolve_sentinels", return_value=pd.DataFrame()), \
             patch.object(biology_logic, "sentinel_summary", return_value={"evaluated_success_count": 0, "evaluated_count": 0, "total_count": 0}), \
             patch.object(biology_logic, "apply_biological_bonus"), \
             patch.object(biology_logic, "apply_functional_penalty"), \
             patch.object(biology_logic, "_calculate_congestion", return_value=np.zeros(2)), \
             patch("src.metrics.plot_dose_response"), \
             patch("src.metrics.plot_congestion"), \
             patch.object(biology_logic, "_hill_inhibition", side_effect=tracked_hill):
            result = biology_logic.run_pharmacological_dose_response(
                graph, self._scores(), Path(directory), spesifik_hedefler=["ENSP00000003084"],
                dozlar=[0.9, 0.5, 0.1], hill_n=2,
            )

        self.assertIsInstance(result, pd.DataFrame)
        self.assertEqual(result["Survival_Fraction"].tolist(), [0.9, 0.5, 0.1])
        self.assertEqual(result["Survival_Fraction"].nunique(), 3)
        self.assertEqual(result["Hill_n"].tolist(), [2, 2, 2])
        self.assertEqual(observed_survival_fractions, [0.9, 0.5, 0.1])
        self.assertEqual(list(graph.es["weight"]), original_weight)
        self.assertEqual(list(graph.es["distance"]), original_distance)
