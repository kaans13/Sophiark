from __future__ import annotations

import tempfile
import unittest
from pathlib import Path
import sqlite3

import pandas as pd

from src.research.entity_resolution import EntityResolver, load_human_project_index
from src.signaling.analysis import analyze_signaling_context, directed_convergence, safe_analyze_signaling_context
from src.signaling.graph import DirectedSignalingGraph, path_sign
from src.signaling.ingestion import load_omnipath_signaling_dataset
from src.signaling.models import (
    DatasetQA, LayerMembership, MappingReport, SignalingDataset, SignalingEdge,
    SignalingNode, SignalingStatus, SignStatus,
)
from src.signaling.storage import load_signaling_cache, write_signaling_cache


def edge(source: str, target: str, sign: SignStatus) -> SignalingEdge:
    canonical = 1 if sign is SignStatus.ACTIVATION else -1 if sign is SignStatus.INHIBITION else 0
    return SignalingEdge(
        source=source, target=target, source_identifier=source, source_symbol=source,
        target_identifier=target, target_symbol=target, directed=True,
        stimulation=sign in {SignStatus.ACTIVATION, SignStatus.CONFLICTING},
        inhibition=sign in {SignStatus.INHIBITION, SignStatus.CONFLICTING},
        canonical_sign=canonical, sign_status=sign, interaction="fixture",
        resources=("fixture",), references=(f"PMID:{source}{target}",),
    )


def graph_fixture() -> DirectedSignalingGraph:
    names = ("A", "B", "C", "D", "E", "X")
    nodes = tuple(SignalingNode(
        node_id=name, source_identifier=name, symbol=name,
        layer_membership=LayerMembership.LAYER_1_AND_2, layer1_id=name,
        mapping_status="EXACT",
    ) for name in names)
    edges = (
        edge("A", "B", SignStatus.ACTIVATION),
        edge("B", "C", SignStatus.ACTIVATION),
        edge("C", "D", SignStatus.INHIBITION),
        edge("A", "E", SignStatus.INHIBITION),
        edge("E", "D", SignStatus.INHIBITION),
        edge("D", "B", SignStatus.ACTIVATION),  # cycle
        edge("X", "C", SignStatus.UNSIGNED),
    )
    qa = DatasetQA(7, 7, 6, 7, 4, 2, 1, 0, 0, 7, 0)
    mapping = MappingReport(6, 6, 0, 0, 0, 7, 7, 0, 0, 0, 1.0, 1.0, 6, 6, 7, 7, 1.0, 1.0, "PREFERRED")
    return DirectedSignalingGraph(SignalingDataset(nodes, edges, qa, mapping, {}, "fixture-hash"))


class DirectedGraphTests(unittest.TestCase):
    def test_direction_is_not_implicitly_reversed(self):
        graph = graph_fixture()
        self.assertEqual(graph.bounded_distances("A", max_depth=1), {"B": 1, "E": 1})
        self.assertNotIn("A", graph.bounded_distances("B", max_depth=1))

    def test_sign_products_and_uncertainty(self):
        self.assertEqual(path_sign((edge("A", "B", SignStatus.ACTIVATION), edge("B", "C", SignStatus.ACTIVATION))), "activation")
        self.assertEqual(path_sign((edge("A", "B", SignStatus.ACTIVATION), edge("B", "C", SignStatus.INHIBITION))), "inhibition")
        self.assertEqual(path_sign((edge("A", "B", SignStatus.INHIBITION), edge("B", "C", SignStatus.INHIBITION))), "activation")
        self.assertEqual(path_sign((edge("A", "B", SignStatus.UNSIGNED),)), "uncertain")
        self.assertEqual(path_sign((edge("A", "B", SignStatus.CONFLICTING),)), "uncertain")

    def test_shortest_depth_cycle_protection_and_bounded_paths(self):
        graph = graph_fixture()
        summary = graph.summarize_paths("A", "D", max_depth=4)
        self.assertEqual(summary.shortest_distance, 2)
        self.assertEqual(summary.path_count, 2)
        self.assertEqual(summary.net_sign, "uncertain")  # one positive and one negative route
        self.assertEqual(graph.summarize_paths("A", "D", max_depth=1).path_count, 0)
        self.assertLessEqual(len(graph.bounded_paths("A", "D", max_depth=10)), 10)

    def test_multi_target_relations_and_convergence(self):
        graph = graph_fixture()
        convergence = directed_convergence(graph, ("A", "X"), max_depth=3)
        record = next(item for item in convergence if item.convergence_node == "C")
        self.assertEqual(record.distances, {"A": 2, "X": 1})
        report = pd.DataFrame({"gene": ["D"], "Delta_PageRank_Pct": [2.71], "Response_Direction": ["gain"]})
        result = analyze_signaling_context(graph, report, targets=("A", "X"), max_depth=4)
        self.assertEqual(len(result.rows), 2)
        self.assertEqual({item.perturbation_target for item in result.rows}, {"A", "X"})

    def test_failure_isolated_from_core_report(self):
        report = pd.DataFrame({"gene": ["A"], "Delta_PageRank": [0.1]})
        before = report.copy(deep=True)
        result = safe_analyze_signaling_context(None, report, targets=("A",))
        self.assertEqual(result.status, SignalingStatus.UNAVAILABLE)
        pd.testing.assert_frame_equal(report, before)

    def test_enabled_and_disabled_layer_preserve_all_core_regression_fields_exactly(self):
        core_fields = [
            "PageRank", "BC_Skoru", "Hinterland_Skoru", "Delta_PageRank",
            "Systemic_Network_Shift_Pct", "Global_Efficiency_Baseline",
            "Global_Efficiency_Perturbed", "Mean_Local_Efficiency_Baseline",
            "Mean_Local_Efficiency_Perturbed", "Compartment_Bottleneck",
        ]
        report = pd.DataFrame([dict(zip(core_fields, [
            .25, 4.0, 80.0, .002, .7, .9, .89, .8, .79, "gateway",
        ]), gene="D")])
        golden = report.copy(deep=True)
        enabled = safe_analyze_signaling_context(graph_fixture(), report, targets=("A",))
        disabled = safe_analyze_signaling_context(None, report, targets=("A",))
        self.assertEqual(enabled.status, SignalingStatus.AVAILABLE)
        self.assertEqual(disabled.status, SignalingStatus.UNAVAILABLE)
        pd.testing.assert_frame_equal(report[core_fields], golden[core_fields], check_exact=True)

    def test_research_explorer_export_contains_signaling_and_provenance(self):
        from src.research.ui import build_research_export_tables, render_research_explorer
        from tests.test_research_ui import FakeUI, _bundle

        bundle = _bundle()
        report = pd.DataFrame({"gene": ["D"], "Delta_PageRank_Pct": [2.71]})
        bundle.signaling_result = analyze_signaling_context(
            graph_fixture(), report, targets=("A",), max_depth=4,
        )
        tables = build_research_export_tables(bundle)
        row = tables["Directed Signaling Context"][0]
        self.assertEqual(row["network_response_value"], 2.71)
        self.assertTrue(row["example_path"])
        self.assertIn("fixture", row["evidence_resources"])
        self.assertTrue(row["references"])
        ui = FakeUI()
        render_research_explorer(bundle, ui=ui, selected_section="Overview")
        self.assertIn("Yönlü sinyalleme bağlamı", ui.text)


class MappingIngestionTests(unittest.TestCase):
    def test_sqlite_cache_round_trip_indexes_and_stale_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.tsv"
            source.write_text("fixture-v1", encoding="utf-8")
            cache = root / "signaling.sqlite"
            dataset = graph_fixture().dataset
            write_signaling_cache(dataset, cache, source_snapshot=source)
            restored = load_signaling_cache(cache, source_snapshot=source)
            self.assertEqual(restored.fingerprint, dataset.fingerprint)
            key = lambda item: (item.source, item.target)
            self.assertEqual(tuple(sorted(restored.edges, key=key)), tuple(sorted(dataset.edges, key=key)))
            connection = sqlite3.connect(cache)
            try:
                indexes = {row[1] for row in connection.execute("PRAGMA index_list(omnipath_signaling_edges)")}
            finally:
                connection.close()
            self.assertTrue({"idx_omnipath_signaling_source", "idx_omnipath_signaling_target", "idx_omnipath_signaling_pair"} <= indexes)
            source.write_text("fixture-v2", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "Stale"):
                load_signaling_cache(cache, source_snapshot=source)

    def test_valid_ambiguous_unmapped_and_layer2_only(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            processed = root / "data" / "processed"
            processed.mkdir(parents=True)
            pd.DataFrame({
                "gene": ["ENSP_A", "ENSP_B1", "ENSP_B2"],
                "Symbol": ["A", "B", "B"],
            }).to_csv(processed / "ensp_with_symbols.csv", index=False)
            index = load_human_project_index(root)
            resolver = EntityResolver((index,))
            snapshot = root / "omnipath_9606.tsv"
            pd.DataFrame([
                {
                    "source": "P-A", "target": "P-B", "source_genesymbol": "A", "target_genesymbol": "B",
                    "is_directed": True, "is_stimulation": True, "is_inhibition": False,
                    "consensus_direction": True, "consensus_stimulation": True,
                    "consensus_inhibition": False, "sources": "DB1", "references": "PMID:1",
                },
                {
                    "source": "P-A", "target": "COMPLEX:X_Y", "source_genesymbol": "A", "target_genesymbol": "X_Y",
                    "is_directed": True, "is_stimulation": False, "is_inhibition": True,
                    "consensus_direction": True, "consensus_stimulation": False,
                    "consensus_inhibition": True, "sources": "DB2", "references": "PMID:2",
                },
                {
                    "source": "P-A", "target": "P-B", "source_genesymbol": "A", "target_genesymbol": "B",
                    "is_directed": True, "is_stimulation": False, "is_inhibition": True,
                    "consensus_direction": True, "consensus_stimulation": False,
                    "consensus_inhibition": True, "sources": "DB3", "references": "PMID:3",
                },
            ]).to_csv(snapshot, sep="\t", index=False)
            dataset = load_omnipath_signaling_dataset(
                snapshot, resolver=resolver, taxon_id=9606, minimum_rows=1,
            )
        self.assertEqual(dataset.mapping.mapped_layer1_entities, 1)
        self.assertEqual(dataset.mapping.ambiguous_entities, 1)
        self.assertEqual(dataset.mapping.signaling_only_entities, 2)
        self.assertTrue(any(node.layer_membership is LayerMembership.SIGNALING_ONLY for node in dataset.nodes))
        self.assertEqual(dataset.qa.unique_directed_edges, 2)
        self.assertEqual(dataset.qa.conflicting_edges, 1)
        self.assertEqual(dataset.qa.duplicates_removed, 1)

    def test_abnormally_small_production_snapshot_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            snapshot = Path(directory) / "omnipath.tsv"
            pd.DataFrame([{
                "source": "P-A", "target": "P-B", "source_genesymbol": "A", "target_genesymbol": "B",
                "is_directed": True, "is_stimulation": True, "is_inhibition": False,
                "sources": "DB", "references": "PMID:1",
            }]).to_csv(snapshot, sep="\t", index=False)
            with self.assertRaisesRegex(ValueError, "abnormally small"):
                load_omnipath_signaling_dataset(
                    snapshot, resolver=EntityResolver(()), taxon_id=9606,
                )


if __name__ == "__main__":
    unittest.main()
