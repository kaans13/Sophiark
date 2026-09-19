from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import igraph as ig
import pandas as pd

from src.services.propagation_trace import build_propagation_trace, export_propagation_trace, filter_trace_nodes


def fixture_graph():
    graph = ig.Graph(edges=[(0, 1), (1, 2), (2, 3), (0, 4)], directed=False)
    graph.vs["name"] = ["T", "A", "B", "C", "D"]
    graph.es["weight"] = [1.0] * 4
    graph.es["distance"] = [1.0] * 4
    return graph


def fixture_report():
    return pd.DataFrame({
        "gene": ["T", "A", "B", "C", "D"], "Delta_PageRank": [None, .1, -.2, .05, .01],
        "Abs_Delta_PageRank": [None, .1, .2, .05, .01], "Response_Direction": ["N/A", "Influence Gain", "Influence Loss", "Influence Gain", "Influence Gain"],
        "empirical_p": [None, .01, .02, .3, .4], "q_value": [None, .02, .03, .5, .6],
        "significant_redistribution": [False, True, True, False, False], "Strong_Redistribution": [False, True, True, False, False],
        "GÃ¼mrÃ¼k_Kapisi": [False, False, True, False, False], "BC_Skoru": [1, 2, 3, 2, 1], "Topluluk_ID": [1, 1, 1, 2, 1],
    })


class PropagationTraceTests(unittest.TestCase):
    def _snapshot(self, root: Path):
        path = root / "omnipath_9606.tsv"
        pd.DataFrame([
            {"source_genesymbol": "T", "target_genesymbol": "A", "is_directed": True, "is_stimulation": True, "is_inhibition": False, "consensus_direction": True, "consensus_stimulation": True, "consensus_inhibition": False, "sources": "DB", "references": "PMID:1"},
            {"source_genesymbol": "A", "target_genesymbol": "B", "is_directed": True, "is_stimulation": True, "is_inhibition": True, "consensus_direction": True, "consensus_stimulation": True, "consensus_inhibition": True, "sources": "DB", "references": "PMID:2"},
        ]).to_csv(path, sep="\t", index=False)
        (root / "omnipath_9606_metadata.json").write_text(json.dumps({"release": "test"}), encoding="utf-8")
        return path

    def test_layers_routes_evidence_and_determinism(self):
        with tempfile.TemporaryDirectory() as directory:
            snapshot = self._snapshot(Path(directory))
            kwargs = dict(graph=fixture_graph(), report=fixture_report(), target_gene="T", baseline_scores=pd.DataFrame(), gene_to_symbol={x: x for x in "TABCD"}, omnipath_path=snapshot, max_nodes=10, max_routes=4)
            first, second = build_propagation_trace(**kwargs), build_propagation_trace(**kwargs)
            self.assertEqual(first.nodes["Propagation_Layer"].tolist(), second.nodes["Propagation_Layer"].tolist())
            self.assertEqual(set(first.nodes["Propagation_Layer"]), {0, 1, 2, 3})
            self.assertTrue((first.evidence["Evidence_Status"] == "conflicting_evidence").any())
            self.assertTrue((first.edges["Structural_Edge"] == True).all())
            self.assertGreater(len(first.routes), 0)

    def test_filters_no_route_and_export(self):
        trace = build_propagation_trace(graph=fixture_graph(), report=fixture_report(), target_gene="T", baseline_scores=pd.DataFrame(), gene_to_symbol={x: x for x in "TABCD"}, max_nodes=10, max_routes=1)
        self.assertEqual(len(filter_trace_nodes(trace.nodes, layer_limit=1)), 3)
        self.assertLessEqual(len(trace.routes), 1)
        with tempfile.TemporaryDirectory() as directory:
            root = export_propagation_trace(trace, Path(directory) / "trace")
            self.assertTrue((root / "nodes.csv").exists())
            self.assertTrue((root / "metadata.json").exists())

    def test_missing_target_returns_empty_trace(self):
        trace = build_propagation_trace(graph=fixture_graph(), report=fixture_report(), target_gene="ABSENT", baseline_scores=pd.DataFrame(), gene_to_symbol={})
        self.assertTrue(trace.nodes.empty)

    def test_human_and_mouse_identifier_traces(self):
        for target, neighbour in [("ENSP000001", "ENSP000002"), ("ENSMUSP000001", "ENSMUSP000002")]:
            graph = ig.Graph(edges=[(0, 1)], directed=False)
            graph.vs["name"] = [target, neighbour]
            graph.es["weight"] = [1.0]
            graph.es["distance"] = [1.0]
            report = pd.DataFrame({"gene": [target, neighbour], "Abs_Delta_PageRank": [0.0, 0.1], "Delta_PageRank": [None, 0.1], "significant_redistribution": [False, True], "Strong_Redistribution": [False, True]})
            trace = build_propagation_trace(graph=graph, report=report, target_gene=target, baseline_scores=pd.DataFrame(), gene_to_symbol={target: "Target", neighbour: "Neighbor"})
            self.assertEqual(set(trace.nodes["Network_ID"]), {target, neighbour})


if __name__ == "__main__":
    unittest.main()
