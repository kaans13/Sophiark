from __future__ import annotations

import json
import tempfile
import unittest
from pathlib import Path

import pandas as pd

from src.scientific.mygene_cache import PersistentMyGeneCache
from src.scientific.pubmed_validation import calculate_concordance, extract_symbols_from_texts
from src.scientific.regulatory import load_omnipath_overlay
from src.scientific.sentinel import resolve_sentinels, sentinel_summary


class OmniPathTests(unittest.TestCase):
    def test_snapshot_preserves_conflict_and_provenance(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            snapshot = root / "omnipath_9606.tsv"
            pd.DataFrame([{
                "source_genesymbol": "A", "target_genesymbol": "B", "is_directed": True,
                "is_stimulation": True, "is_inhibition": True, "consensus_direction": True,
                "consensus_stimulation": True, "consensus_inhibition": True,
                "sources": "DB1;DB2", "references": "PMID:1",
            }]).to_csv(snapshot, sep="\t", index=False)
            (root / "omnipath_metadata.json").write_text(json.dumps({"release": "test-v1"}), encoding="utf-8")
            frame, metadata = load_omnipath_overlay(snapshot)
            self.assertEqual(frame.iloc[0]["evidence_status"], "conflicting_evidence")
            self.assertEqual(frame.iloc[0]["provenance"], "DB1;DB2")
            self.assertEqual(metadata["records"], 1)

    def test_missing_snapshot_is_not_negative_evidence(self):
        frame, metadata = load_omnipath_overlay(None)
        self.assertTrue(frame.empty)
        self.assertEqual(metadata["status"], "OmniPath unavailable")


class PubMedValidationTests(unittest.TestCase):
    def test_deterministic_exact_symbol_extraction_and_overlap(self):
        literature, excluded = extract_symbols_from_texts(["CFTR and SCNN1B interact; AR is not resolved."], ["CFTR", "SCNN1B", "AR"])
        self.assertEqual(literature, {"CFTR", "SCNN1B"})
        self.assertEqual(excluded, ["AR"])
        first = calculate_concordance(["CFTR", "TP53"], literature, ["CFTR", "SCNN1B", "TP53", "EGFR"])
        second = calculate_concordance(["CFTR", "TP53"], literature, ["CFTR", "SCNN1B", "TP53", "EGFR"])
        self.assertEqual(first, second)
        self.assertEqual(first["intersection_count"], 1)


class SentinelTests(unittest.TestCase):
    def test_legacy_identifier_resolves_through_local_symbol(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            mapping = root / "mapping.csv"
            pd.DataFrame({"gene": ["ENSP00000441543"], "Symbol": ["UBC"]}).to_csv(mapping, index=False)
            scores = pd.DataFrame({"gene": ["ENSP00000441543"], "Hinterland_Skoru": [99.0]})
            result = resolve_sentinels(scores, mapping)
            ubc = result[result["resolved_symbol"] == "UBC"].iloc[0]
            self.assertEqual(ubc["status"], "PASS")
            summary = sentinel_summary(result)
            self.assertLess(summary["coverage_rate"], 1.0)


class MyGeneCacheTests(unittest.TestCase):
    def test_positive_negative_and_species_separated_cache(self):
        with tempfile.TemporaryDirectory() as directory:
            cache = PersistentMyGeneCache(Path(directory) / "mygene.sqlite")
            calls = []
            def remote(ids):
                calls.append(ids)
                return [{"query": "A", "symbol": "A"}]
            first = cache.resolve_batch("mouse", ["A", "MISSING"], "go.CC", remote)
            second = cache.resolve_batch("mouse", ["A", "MISSING"], "go.CC", remote)
            other_species = cache.resolve_batch("human", ["A"], "go.CC", remote)
            self.assertEqual(first["A"]["symbol"], "A")
            self.assertIsNone(second["MISSING"])
            self.assertEqual(len(calls), 2)  # first mouse + independent human
            self.assertEqual(other_species["A"]["symbol"], "A")
            self.assertGreaterEqual(cache.stats.negative_hits, 1)


if __name__ == "__main__":
    unittest.main()
