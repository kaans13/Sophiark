from __future__ import annotations

import unittest

import igraph as ig
import pandas as pd

from src.services.report_context import critical_rows, symbol_first_report
from src.ui.network_visualization import build_network_figure


class UIRecoveryTests(unittest.TestCase):
    def setUp(self):
        self.graph = ig.Graph.Ring(8)
        self.graph.vs["name"] = [f"ENSP{i}" for i in range(8)]
        self.report = pd.DataFrame({
            "gene": self.graph.vs["name"],
            "Hinterland_Skoru": range(8),
            "Delta_PageRank": [value / 100 for value in range(8)],
            "Abs_Delta_PageRank": [value / 100 for value in range(8)],
            "Strong_Redistribution": [False, True, False, True, False, True, False, True],
            "significant_redistribution": [False, False, False, True, False, False, False, True],
        })
        self.mapping = {f"ENSP{i}": f"GENE{i}" for i in range(8)}

    def test_symbol_first_table_keeps_identifier_secondary(self):
        shown = symbol_first_report(self.report, self.mapping)
        self.assertEqual(list(shown.columns[:2]), ["Gen Sembolü", "Protein Kimliği"])
        self.assertEqual(shown.iloc[0]["Gen Sembolü"], "GENE0")
        self.assertEqual(shown.iloc[0]["Protein Kimliği"], "ENSP0")

    def test_critical_rows_uses_existing_redistribution_flag(self):
        selected = critical_rows(self.report)
        self.assertEqual(set(selected["gene"]), {"ENSP1", "ENSP3", "ENSP5", "ENSP7"})

    def test_2d_and_3d_networks_render_alias_hover_without_id_labels(self):
        for dimensions in (2, 3):
            figure = build_network_figure(
                self.graph, self.report, ["ENSP0"], self.mapping,
                dimensions=dimensions, max_nodes=8, show_labels=False,
            )
            self.assertIsNotNone(figure)
            self.assertEqual(len(figure.data), 2)
            self.assertTrue(all(label == "" for label in figure.data[1].text))
            self.assertIn("GENE0", figure.data[1].hovertext[0])
            self.assertIn("Protein kimliği", figure.data[1].hovertext[0])


if __name__ == "__main__":
    unittest.main()
