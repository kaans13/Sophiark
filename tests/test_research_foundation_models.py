"""Focused tests for immutable Research Context foundation models."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, fields
from datetime import datetime, timezone
import unittest

import numpy as np
import pandas as pd

from src.research.contracts import (
    CONTRACTS,
    contracts_for_source_field,
    format_metric,
    get_contract,
    get_metric_definition,
    get_metric_label,
    get_metric_limitations,
)
from src.research.models import (
    AliasType,
    Entity,
    EntityAlias,
    EntityResolution,
    EntityType,
    Observation,
    ObservationScope,
    ObservationType,
    ResolutionStatus,
)
from src.research.snapshot import (
    build_simulation_result_snapshot,
    freeze_result_table,
    thaw_result_table,
)


class SnapshotTests(unittest.TestCase):
    @staticmethod
    def _frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        report = pd.DataFrame({
            "gene": ["ENSP_TARGET", "ENSP_GAIN"],
            "Hasar_Tipi": ["Birincil_Hasar", "Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [np.nan, 2.84],
            "metadata": [{"source": "core"}, ["local", "computed"]],
        })
        signed = pd.DataFrame({
            "gene": ["ENSP_GAIN", "ENSP_LOSS"],
            "Delta_PageRank_Pct": [2.84, -1.25],
            "Response_Direction": ["Influence Gain", "Influence Loss"],
        })
        enrichment = pd.DataFrame({
            "Term": ["Example term"], "Adjusted P-value": [0.02],
            "Genes": ["GAIN;LOSS"], "Pathway_Family": ["GO:BP"],
        })
        return report, signed, enrichment

    def _build(self, *, created_at: str, threshold_parameters=None):
        report, signed, enrichment = self._frames()
        return build_simulation_result_snapshot(
            report, signed, enrichment,
            species="Homo sapiens", taxon_id=9606, tissue="Pancreas",
            targets=["ENSP_TARGET"], attenuation=0.001,
            selection_mode="top_n", test_limit=200, tested_count=200,
            returned_count=20, top_n=20,
            threshold_parameters=threshold_parameters or {"minimum_delta_pct": 0.05, "threshold": 0.05},
            created_at=created_at,
        )

    def test_snapshot_is_content_addressed_and_time_independent(self) -> None:
        first = self._build(created_at="2026-01-01T00:00:00+00:00")
        second = self._build(created_at="2026-02-01T00:00:00+00:00")
        self.assertEqual(first.scientific_fingerprint, second.scientific_fingerprint)
        self.assertEqual(first.snapshot_fingerprint, second.snapshot_fingerprint)
        self.assertEqual(first.simulation_id, second.simulation_id)
        self.assertEqual(first.snapshot_id, second.snapshot_id)
        self.assertNotEqual(first.created_at, second.created_at)

    def test_scope_and_provenance_are_explicit(self) -> None:
        snapshot = self._build(created_at="2026-01-01T00:00:00Z")
        self.assertEqual(snapshot.scope.selection_mode, "top_n")
        self.assertEqual(snapshot.scope.tested_count, 200)
        self.assertEqual(snapshot.scope.returned_count, 20)
        self.assertEqual(snapshot.scope.top_n, 20)
        self.assertEqual(snapshot.scope.species, "Homo sapiens")
        self.assertEqual(snapshot.scope.taxon_id, 9606)
        self.assertEqual(snapshot.scope.threshold, 0.05)
        self.assertEqual(len(snapshot.provenance), 3)
        self.assertTrue(all(item.snapshot_id for item in snapshot.provenance))

    def test_source_frames_and_nested_values_are_detached(self) -> None:
        report, signed, enrichment = self._frames()
        report_before = report.copy(deep=True)
        snapshot = build_simulation_result_snapshot(
            report, signed, enrichment,
            species="human", taxon_id=9606, tissue="Pancreas",
            targets=["ENSP_TARGET"], attenuation=0.001,
            selection_mode="top_n", test_limit=20, tested_count=100,
            returned_count=20, top_n=20, created_at=datetime.now(timezone.utc),
        )
        report.loc[0, "gene"] = "MUTATED"
        report.loc[0, "metadata"]["source"] = "mutated"
        records = snapshot.report.records()
        self.assertEqual(records[0]["gene"], "ENSP_TARGET")
        self.assertEqual(records[0]["metadata"], {"source": "core"})
        records[0]["metadata"]["source"] = "caller mutation"
        self.assertEqual(snapshot.report.records()[0]["metadata"], {"source": "core"})
        self.assertEqual(report_before.loc[1, "gene"], "ENSP_GAIN")

    def test_scientific_or_scope_change_changes_fingerprint(self) -> None:
        baseline = self._build(created_at="2026-01-01T00:00:00Z")
        changed = self._build(
            created_at="2026-01-01T00:00:00Z",
            threshold_parameters={"minimum_delta_pct": 0.10, "threshold": 0.10},
        )
        self.assertNotEqual(baseline.scientific_fingerprint, changed.scientific_fingerprint)

    def test_freeze_thaw_handles_missing_nested_and_nonfinite_without_source_changes(self) -> None:
        source = pd.DataFrame({
            "gene": ["A", "B", "C"],
            "value": [np.nan, np.inf, -np.inf],
            "terms": [["x", "y"], {"b": 2, "a": 1}, None],
        })
        before = source.copy(deep=True)
        frozen = freeze_result_table(source, name="fixture")
        thawed = thaw_result_table(frozen)
        self.assertEqual(frozen.row_count, 3)
        self.assertTrue(pd.isna(thawed.loc[0, "value"]))
        self.assertEqual(thawed.loc[1, "value"], np.inf)
        self.assertEqual(thawed.loc[2, "value"], -np.inf)
        self.assertEqual(thawed.loc[0, "terms"], ["x", "y"])
        self.assertEqual(thawed.loc[1, "terms"], {"a": 1, "b": 2})
        pd.testing.assert_frame_equal(source, before)

    def test_species_taxon_mismatch_is_rejected(self) -> None:
        with self.assertRaisesRegex(ValueError, "species/taxon mismatch"):
            build_simulation_result_snapshot(
                pd.DataFrame(), None, None,
                species="Mus musculus", taxon_id=9606, tissue="None", targets=[],
                attenuation=None, selection_mode="top_n", test_limit=0,
                tested_count=0, returned_count=0,
            )

    def test_snapshot_fields_and_scope_mappings_are_immutable(self) -> None:
        snapshot = self._build(created_at="2026-01-01T00:00:00Z")
        with self.assertRaises(FrozenInstanceError):
            snapshot.species = "Mouse"  # type: ignore[misc]
        with self.assertRaises(TypeError):
            snapshot.threshold_parameters["threshold"] = 1.0  # type: ignore[index]


class EntityAndObservationModelTests(unittest.TestCase):
    def test_entity_key_is_species_and_type_safe(self) -> None:
        human = Entity(9606, EntityType.PROTEIN, "ENSP0001", symbol="SAME")
        mouse = Entity(10090, EntityType.PROTEIN, "ENSMUSP0001", symbol="SAME")
        self.assertEqual(human.logical_key, "9606:protein:ENSP0001")
        self.assertNotEqual(human.logical_key, mouse.logical_key)

    def test_resolution_requires_alias_provenance_and_never_selects_ambiguous(self) -> None:
        entity = Entity(9606, EntityType.PROTEIN, "ENSP0001", symbol="ABC")
        alias = EntityAlias(AliasType.GENE_SYMBOL, "ABC", "local_symbol_map")
        resolved = EntityResolution(
            query="ABC", taxon_id=9606, entity_type=EntityType.PROTEIN,
            status=ResolutionStatus.ALIAS, entity=entity, matched_alias=alias,
            source="local_symbol_map",
        )
        self.assertIs(resolved.entity, entity)
        with self.assertRaises(ValueError):
            EntityResolution(
                query="ABC", taxon_id=9606, entity_type=EntityType.PROTEIN,
                status=ResolutionStatus.AMBIGUOUS, entity=entity, candidates=(entity,),
            )

    def test_observation_has_no_score_and_preserves_scope(self) -> None:
        entity = Entity(9606, EntityType.PROTEIN, "ENSP0001", symbol="ABC")
        scope = ObservationScope(
            selection_mode="top_n", tested_count=200, returned_count=20,
            threshold=None, top_n=20, tissue="Pancreas",
            species="Homo sapiens", taxon_id=9606,
        )
        observation = Observation(
            observation_id="obs-1", simulation_id="sim-1",
            type=ObservationType.FAMILY_COOCCURRENCE, members=(entity,),
            basis="Returned members share an authoritative annotation.",
            scope=scope, source_fields=("family",),
            created_at="2026-01-01T00:00:00+00:00",
            limitations="This does not establish enrichment or causality.",
        )
        self.assertEqual(observation.scope.returned_count, 20)
        self.assertNotIn("score", {item.name for item in fields(observation)})


class InterpretationContractTests(unittest.TestCase):
    def test_master_spec_contracts_are_covered(self) -> None:
        required = {
            "delta_pagerank_pct", "pagerank_baseline", "pagerank_perturbed",
            "hinterland", "betweenness_centrality", "gateway", "degree", "community",
            "drug_score", "essentiality", "empirical_p", "q_value",
            "significant_redistribution", "global_efficiency_change_pct",
            "mean_local_efficiency_change_pct", "target_local_efficiency_change_pct",
            "go_enrichment", "kegg_enrichment", "localization", "tf_annotations",
        }
        self.assertTrue(required.issubset(CONTRACTS))

    def test_existing_ui_semantics_are_reused_for_verified_fields(self) -> None:
        expected = {
            "Delta_PageRank_Pct": ("Ağ önemi değişimi (%)", "+2.8400%"),
            "Hinterland_Skoru": ("Ağdaki Önemi", "64.86"),
            "BC_Skoru": ("Geçiş Merkeziliği", "26252.000"),
            "PageRank_Baseline": ("Başlangıç ağ önemi", "0.123457"),
            "Gümrük_Kapisi": ("Compartment Bottleneck", "✓"),
            "q_value": ("q-değeri (FDR)", "0.020000"),
            "Efficiency_Kayip_Pct": ("Sistemik ağ kayması (%) · legacy", "+2.00%"),
            "Local_Efficiency_Kayip_Pct": ("Pozitif ağ önemi yanıtı ortalaması (%) · legacy", "+2.00%"),
        }
        values = {
            "Delta_PageRank_Pct": 2.84, "Hinterland_Skoru": 64.86,
            "BC_Skoru": 26252, "PageRank_Baseline": 0.1234567,
            "Gümrük_Kapisi": True, "q_value": 0.02,
            "Efficiency_Kayip_Pct": 2.0, "Local_Efficiency_Kayip_Pct": 2.0,
        }
        for field, (label, rendered) in expected.items():
            with self.subTest(field=field):
                self.assertEqual(get_metric_label(field), label)
                self.assertEqual(format_metric(field, values[field]), rendered)

    def test_contract_boundaries_and_missing_values_are_deterministic(self) -> None:
        self.assertIn("gene expression changed", get_metric_limitations("Delta_PageRank_Pct"))
        self.assertIn("it is not global efficiency", get_metric_definition("legacy_efficiency_kayip_pct"))
        self.assertEqual(format_metric("Drug_Score", np.nan), "—")
        self.assertEqual(format_metric("unknown", None), "—")
        self.assertIsNone(get_contract("unknown"))

    def test_shared_enrichment_fields_remain_explicitly_ambiguous(self) -> None:
        matches = contracts_for_source_field("Term")
        self.assertEqual({item.metric_id for item in matches}, {"go_enrichment", "kegg_enrichment"})
        self.assertIsNone(get_contract("Term"))


if __name__ == "__main__":
    unittest.main()
