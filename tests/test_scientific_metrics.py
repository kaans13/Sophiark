from __future__ import annotations

import sqlite3
from contextlib import closing
import tempfile
import unittest
from pathlib import Path

import igraph as ig
import numpy as np

from src.scientific.efficiency import (
    mean_local_efficiency,
    strength_to_distance,
    weighted_global_efficiency,
)
from src.scientific.perturbation import SEVERITY_FRACTIONS
from src.scientific.redistribution import classify_redistribution
from src.scientific.statistics import benjamini_hochberg
from src.scientific.statistics import empirical_pvalues_degree_matched
from src.scientific.regulatory import load_omnipath_overlay, trrust_records
from src.scientific.tissue import pooled_expression_reference, tissue_multiplier
from src.scientific.validation import sampled_betweenness
from src.services.enrichment import LIBRARIES, local_over_representation


def graph(edges, weights, n=None):
    item = ig.Graph(n=n, edges=edges, directed=False)
    item.vs["name"] = [f"g{i}" for i in range(item.vcount())]
    item.es["weight"] = weights
    item.es["distance"] = [strength_to_distance(value) for value in weights]
    return item


class EfficiencySanityTests(unittest.TestCase):
    def test_complete_graph(self):
        item = graph([(i, j) for i in range(4) for j in range(i + 1, 4)], [1.0] * 6, n=4)
        self.assertAlmostEqual(weighted_global_efficiency(item)[0], 1.0)
        self.assertAlmostEqual(mean_local_efficiency(item)[0], 1.0)

    def test_disconnected_graph(self):
        item = graph([(0, 1), (2, 3)], [1.0, 1.0], n=4)
        self.assertAlmostEqual(weighted_global_efficiency(item)[0], 1.0 / 3.0)

    def test_star_and_path_graphs(self):
        star = graph([(0, 1), (0, 2), (0, 3)], [1.0] * 3, n=4)
        path = graph([(0, 1), (1, 2)], [1.0, 1.0], n=3)
        self.assertAlmostEqual(weighted_global_efficiency(star)[0], 0.75)
        self.assertAlmostEqual(mean_local_efficiency(star)[0], 0.0)
        self.assertAlmostEqual(weighted_global_efficiency(path)[0], 5.0 / 6.0)

    def test_weighted_triangle_uses_inverse_strength(self):
        item = graph([(0, 1), (1, 2), (0, 2)], [1.0, 1.0, 0.25], n=3)
        self.assertAlmostEqual(weighted_global_efficiency(item)[0], 5.0 / 6.0)
        self.assertTrue(np.isinf(strength_to_distance(0.0)))
        self.assertTrue(np.isinf(strength_to_distance(-1.0)))

    def test_stronger_attenuation_reduces_efficiency_on_star(self):
        baseline = graph([(0, 1), (0, 2), (0, 3)], [1.0] * 3, n=4)
        values = []
        for name in ("Mild", "Moderate", "Strong", "Near-complete"):
            perturbed = baseline.copy()
            fraction = SEVERITY_FRACTIONS[name]
            perturbed.es["weight"] = [fraction] * 3
            perturbed.es["distance"] = [strength_to_distance(fraction)] * 3
            values.append(weighted_global_efficiency(perturbed)[0])
        self.assertEqual(values, sorted(values, reverse=True))


class StatisticalTests(unittest.TestCase):
    def test_bh_is_monotone_in_ranked_order(self):
        p = np.array([0.01, 0.04, 0.03, 0.002])
        q = benjamini_hochberg(p)
        order = np.argsort(p)
        self.assertTrue(np.all(np.diff(q[order]) >= -1e-15))
        self.assertTrue(np.all((0 <= q) & (q <= 1)))

    def test_two_sided_redistribution(self):
        frame = classify_redistribution(
            {"a": .5, "b": .3, "c": .2}, {"a": .4, "b": .4, "c": .2},
            empirical_p={"a": .01, "b": .01, "c": 1}, target_names=[],
            selection_mode="top_n", exploratory_top_n=3, percentile=90, fdr_alpha=.05,
        )
        directions = dict(zip(frame["gene"], frame["Response_Direction"]))
        self.assertEqual(directions["a"], "Influence Loss")
        self.assertEqual(directions["b"], "Influence Gain")

    def test_degree_matched_null_is_reproducible(self):
        item = graph([(0, 1), (1, 2), (2, 3), (3, 4), (1, 3)], [1.0] * 5, n=5)
        baseline_values = item.pagerank(weights="weight", directed=False)
        baseline = dict(zip(item.vs["name"], baseline_values))
        observed = {name: 1e-5 for name in item.vs["name"]}
        first, meta = empirical_pvalues_degree_matched(
            item, baseline, observed, ["g1"], attenuation_fraction=.1,
            damping=.85, iterations=5, random_seed=42,
        )
        second, _ = empirical_pvalues_degree_matched(
            item, baseline, observed, ["g1"], attenuation_fraction=.1,
            damping=.85, iterations=5, random_seed=42,
        )
        self.assertEqual(first, second)
        self.assertEqual(meta["iterations"], 5)


class RegulatoryLayerTests(unittest.TestCase):
    def test_trrust_kept_as_directed_overlay(self):
        frame = trrust_records({"TP53": [("MDM2", "Repression")]})
        self.assertEqual(frame.iloc[0]["source"], "MDM2")
        self.assertEqual(frame.iloc[0]["target"], "TP53")
        self.assertTrue(bool(frame.iloc[0]["inhibition"]))
        self.assertEqual(frame.iloc[0]["evidence_layer"], "directed_regulatory")

    def test_missing_omnipath_is_not_negative_evidence(self):
        frame, metadata = load_omnipath_overlay(None)
        self.assertTrue(frame.empty)
        self.assertFalse(metadata["available"])
        self.assertIn("negatif evidence değildir", metadata["reason"])


class ReproducibilityTests(unittest.TestCase):
    def test_sampled_bc_same_seed(self):
        item = graph([(0, 1), (1, 2), (2, 3), (1, 3)], [1.0] * 4, n=4)
        first, _ = sampled_betweenness(item, sample_sources=2, random_seed=42)
        second, _ = sampled_betweenness(item, sample_sources=2, random_seed=42)
        np.testing.assert_array_equal(first, second)

    def test_cross_tissue_reference_is_deterministic(self):
        with tempfile.TemporaryDirectory() as directory:
            db = Path(directory) / "expression.db"
            with closing(sqlite3.connect(db)) as connection:
                connection.execute("CREATE TABLE tissue_expression(expression_level REAL)")
                connection.executemany("INSERT INTO tissue_expression VALUES (?)", [(1,), (2,), (3,), (100,)])
                connection.commit()
            pooled_expression_reference.cache_clear()
            import pandas as pd
            a, meta_a = tissue_multiplier(pd.Series([1.0, 3.0]), pd.Series([2.0, 3.0]), mode="cross_tissue_comparable", db_path=str(db), expression_table="tissue_expression")
            b, meta_b = tissue_multiplier(pd.Series([1.0, 3.0]), pd.Series([2.0, 3.0]), mode="cross_tissue_comparable", db_path=str(db), expression_table="tissue_expression")
            np.testing.assert_array_equal(a.to_numpy(), b.to_numpy())
            self.assertEqual(meta_a, meta_b)
            pooled_expression_reference.cache_clear()


class EnrichmentTests(unittest.TestCase):
    def test_human_and_mouse_have_explicit_species_libraries(self):
        self.assertTrue(any("Human" in item.name for item in LIBRARIES["human"]))
        self.assertTrue(any("Mouse" in item.name for item in LIBRARIES["mouse"]))

    def test_local_ora_uses_custom_background(self):
        result = local_over_representation(
            ["A", "B"], ["A", "B", "C", "D", "E"],
            {"pathway": ["A", "B", "X"]},
        )
        self.assertEqual(int(result.iloc[0]["Background_Universe_Size"]), 5)
        self.assertEqual(int(result.iloc[0]["Matched_Genes"]), 2)
        self.assertIn("Adjusted P-value", result)


class ExportTests(unittest.TestCase):
    def test_excel_export_contains_payload(self):
        from src.ui.tables import export_bytes
        payload, mime, extension = export_bytes(
            __import__("pandas").DataFrame({"Sophiark_Version": ["test"], "value": [1.0]}),
            "Excel (.xlsx)",
        )
        self.assertEqual(extension, "xlsx")
        self.assertGreater(len(payload), 1000)
        self.assertIn("spreadsheetml", mime)

    def test_metadata_dict_is_serialized_into_export_column(self):
        import pandas as pd
        from src.scientific.metadata import attach_metadata_columns
        result = attach_metadata_columns(
            pd.DataFrame({"gene": ["A"]}),
            {"cache_snapshot_identifiers": {"symbol_cache": "abc"}},
        )
        self.assertIn('"symbol_cache": "abc"', result.iloc[0]["Cache_Snapshot_IDs"])


if __name__ == "__main__":
    unittest.main()
