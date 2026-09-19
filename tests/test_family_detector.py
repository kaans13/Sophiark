"""Regression tests for the read-only deterministic family/complex detector."""

from __future__ import annotations

import unittest

import pandas as pd

from src.interpretation.family_detector import detect_family_complex_clusters


def _clusters(*symbols: str) -> list[dict[str, object]]:
    frame = pd.DataFrame({
        "Symbol": list(symbols),
        "Delta_PageRank_Pct": [float(index + 1) / 10 for index in range(len(symbols))],
        "Hinterland_Skoru": [10 + index for index in range(len(symbols))],
        "Community_Context": [1] * len(symbols),
    })
    return detect_family_complex_clusters(symbols, candidates=frame, direct_neighbors=symbols[:1])


class FamilyComplexDetectorTests(unittest.TestCase):
    def assert_named_cluster(self, clusters: list[dict[str, object]], needle: str, cluster_type: str | None = None) -> None:
        matching = [item for item in clusters if needle.casefold() in str(item["name"]).casefold()]
        self.assertTrue(matching, f"{needle!r} kümesi bulunamadı: {clusters}")
        if cluster_type:
            self.assertTrue(any(item["cluster_type"] == cluster_type for item in matching))

    def test_curated_paralog_and_family_examples(self) -> None:
        cases = {
            "PRELI": ("PRELID1", "PRELID2", "PRELID3A", "PRELID3B"),
            "Metallothioneins": ("MT1A", "MT1B", "MT1E"),
            "NBPF": ("NBPF1", "NBPF3", "NBPF10"),
            "Taste 2": ("TAS2R1", "TAS2R3", "TAS2R4"),
            "Receptor (G protein": ("RAMP1", "RAMP2", "RAMP3"),
            "Solute carrier family 25": ("SLC25A1", "SLC25A2", "SLC25A3"),
            "Potassium voltage-gated": ("KCNE1", "KCNE2", "KCNE3"),
            "3-hydroxyacyl": ("HACD1", "HACD2", "HACD3"),
            "CatSpermasome": ("CATSPER1", "CATSPER2", "CATSPERB"),
        }
        for expected, symbols in cases.items():
            with self.subTest(expected=expected):
                self.assert_named_cluster(_clusters(*symbols), expected)

    def test_notch2nl_is_explicit_low_confidence_heuristic(self) -> None:
        clusters = _clusters("NOTCH2NLA", "NOTCH2NLB", "NOTCH2NLC")
        item = next(cluster for cluster in clusters if "NOTCH2NL" in str(cluster["name"]))
        self.assertTrue(item["heuristic_used"])
        self.assertEqual(item["confidence_label"], "düşük")
        self.assertTrue(item["warning"])

    def test_srgap2_uses_curated_shared_domain_group(self) -> None:
        clusters = _clusters("SRGAP2", "SRGAP2C")
        self.assert_named_cluster(clusters, "F-BAR")

    def test_calcium_subgroups_remain_distinct_but_share_system(self) -> None:
        clusters = _clusters("CACNA1A", "CACNA1B", "CACNB1", "CACNB2", "CACNG1", "CACNG2", "CACNA2D1", "CACNA2D2")
        self.assert_named_cluster(clusters, "alpha1")
        self.assert_named_cluster(clusters, "beta")
        self.assert_named_cluster(clusters, "gamma")
        self.assert_named_cluster(clusters, "alpha2delta")
        self.assert_named_cluster(clusters, "kalsiyum kanalı alt birimleri", "functional_system")

    def test_adaptor_groups_are_complexes(self) -> None:
        self.assert_named_cluster(_clusters("AP2A2", "AP2B1"), "Adaptor related protein complex 2", "protein_complex")
        self.assert_named_cluster(_clusters("AP1G1", "AP1B1"), "Adaptor related protein complex 1", "protein_complex")

    def test_false_positive_prefixes_do_not_merge(self) -> None:
        self.assertEqual(_clusters("SLC2A1", "SLC25A1"), [])
        clusters = _clusters("CACNA1A", "CACNB1")
        self.assertFalse(any(cluster["cluster_type"] == "family" for cluster in clusters))
        self.assertEqual(_clusters("AP2A2", "AP3B1"), [])
        self.assertEqual(_clusters("RANDOM1", "RANDOM2"), [])

    def test_mettl_is_conservative(self) -> None:
        self.assertEqual(_clusters("METTL1", "METTL3"), [])
        clusters = _clusters("METTL3", "METTL16")
        self.assertTrue(all(not item["heuristic_used"] for item in clusters))

    def test_metrics_and_topology_warning_are_explicit(self) -> None:
        cluster = next(item for item in _clusters("PRELID1", "PRELID2", "PRELID3A", "PRELID3B") if "PRELI" in str(item["name"]))
        self.assertEqual(cluster["member_count"], 4)
        self.assertIn("mean_delta_pagerank_pct", cluster["metrics"])
        self.assertIn(cluster["topology_amplification_risk"], {"düşük", "orta", "yüksek"})
        self.assertTrue(cluster["possible_topology_amplification"])


if __name__ == "__main__":
    unittest.main()
