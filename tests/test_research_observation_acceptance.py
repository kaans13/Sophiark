"""Independent acceptance tests for master clauses 26--34.

These tests deliberately consume only frozen Phase 3 context and typed
relationship evidence.  They do not exercise or modify the scientific core.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import FrozenInstanceError, fields
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import time
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal

from src.research.config import DEFAULT_RESEARCH_CONFIG
from src.research.entity_resolution import EntityResolver, load_human_project_index
from src.research.local_context import LocalResearchContext, build_local_research_context
from src.research.models import (
    EvidenceAssertion,
    Observation,
    ObservationType,
    ProvenanceKind,
    ProvenanceRecord,
    RelationshipType,
)
from src.research.observations import detect_observations
from src.research.providers import (
    FamilyMembership,
    LocalEvidenceAdapters,
    LocalProviderResult,
    LocalProviderStatus,
)
from src.research.relationships import build_target_response_relationships
from src.research.snapshot import build_snapshot


_CREATED_AT = "2026-01-01T00:00:00Z"


def _provenance(source: str, kind: ProvenanceKind = ProvenanceKind.LOCAL_ANNOTATION) -> ProvenanceRecord:
    return ProvenanceRecord(
        source=source,
        kind=kind,
        version="acceptance-v1",
        locator=f"fixture://{source}",
    )


class _NoMatchProvider:
    def __init__(self, provider_id: str, method_name: str) -> None:
        self.provider_id = provider_id
        setattr(self, method_name, self._lookup)

    def _lookup(self, values, *, taxon_id: int) -> LocalProviderResult:
        requested = len(tuple(values))
        return LocalProviderResult(
            provider_id=self.provider_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.NO_MATCH,
            requested_count=requested,
        )


def _providers(*, with_families: bool = False) -> LocalEvidenceAdapters:
    if with_families:
        provenance = _provenance("authoritative_family_fixture")

        def get_families(values, *, taxon_id: int) -> LocalProviderResult:
            records = (
                FamilyMembership(
                    taxon_id, "P1", "P1", "FAMILY-A", "Acceptance family",
                    "approved_symbol", (), None, None, provenance,
                ),
                FamilyMembership(
                    taxon_id, "P2", "P2", "FAMILY-A", "Acceptance family",
                    "approved_symbol", (), None, None, provenance,
                ),
                # This singleton is present specifically to prove that the
                # configured minimum is presentation-only suppression.
                FamilyMembership(
                    taxon_id, "N1", "N1", "FAMILY-SINGLETON", "Singleton family",
                    "approved_symbol", (), None, None, provenance,
                ),
            )
            return LocalProviderResult(
                provider_id="local_hgnc_family",
                taxon_id=taxon_id,
                status=LocalProviderStatus.AVAILABLE,
                records=records,
                provenance=(provenance,),
                requested_count=len(tuple(values)),
                matched_count=3,
            )

        family_provider = SimpleNamespace(get_families=get_families)
    else:
        family_provider = _NoMatchProvider("local_hgnc_family", "get_families")

    return LocalEvidenceAdapters(
        mygene=_NoMatchProvider("local_mygene", "get_annotations"),
        families=family_provider,
        complexes=_NoMatchProvider("local_complex_membership", "get_memberships"),
        uniprot=_NoMatchProvider("local_uniprot_partial", "get_context"),
        regulatory=_NoMatchProvider("local_regulatory", "get_relations"),
    )


def _nested_keys(value) -> set[str]:
    keys: set[str] = set()
    if isinstance(value, Mapping):
        for key, item in value.items():
            keys.add(str(key).casefold())
            keys.update(_nested_keys(item))
    elif isinstance(value, (tuple, list, set, frozenset)):
        for item in value:
            keys.update(_nested_keys(item))
    return keys


class ObservationAcceptanceTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def _resolver(self, identifiers: list[tuple[str, str]]) -> EntityResolver:
        path = self.root / "data" / "processed" / "ensp_with_symbols.csv"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(
            "gene,Symbol\n" + "".join(f"{canonical},{symbol}\n" for canonical, symbol in identifiers),
            encoding="utf-8",
        )
        return EntityResolver((load_human_project_index(self.root),))

    def _context(
        self,
        report: pd.DataFrame,
        signed: pd.DataFrame | None = None,
        enrichment: pd.DataFrame | None = None,
        *,
        identifiers: list[tuple[str, str]],
        targets: tuple[str, ...] = (),
        selection_mode: str = "top_n",
        top_n: int | None = None,
        with_families: bool = False,
    ) -> LocalResearchContext:
        snapshot = build_snapshot(
            report,
            signed,
            enrichment,
            species="Homo sapiens",
            taxon_id=9606,
            tissue="Pancreas",
            targets=targets,
            attenuation=0.01,
            selection_mode=selection_mode,
            test_limit=200,
            tested_count=200,
            returned_count=len(report.index),
            top_n=top_n,
            threshold_parameters={"fdr_alpha": 0.05},
            created_at=_CREATED_AT,
        )
        return build_local_research_context(
            snapshot,
            resolver=self._resolver(identifiers),
            providers=_providers(with_families=with_families),
        )

    def _rich_fixture(self):
        report = pd.DataFrame({
            "gene": ["TARGET", "P1", "P2", "N1"],
            "Symbol": ["T", "P1", "P2", "N1"],
            "Hasar_Tipi": [
                "Birincil_Hasar_Hedef",
                "Üçüncül_Hasar_Stresli",
                "Üçüncül_Hasar_Stresli",
                "Üçüncül_Hasar_Stresli",
            ],
            "Response_Direction": [None, "Influence Gain", "Influence Gain", "Influence Loss"],
            "Delta_PageRank_Pct": [None, 8.0, 4.0, -3.0],
            "Lokalizasyon": ["Cytosol", "Nucleus", "Nucleus", "Membrane"],
            "Gümrük_Kapisi": [False, True, False, True],
            "significant_redistribution": [False, True, False, True],
            # P2 deliberately has the smallest q-value but no stored FDR flag.
            "q_value": [None, 0.02, 0.000001, 0.03],
            "empirical_p": [None, 0.01, 0.000001, 0.02],
            "Düzenleyici_TFler": ["TF-T", "TF-P", "TF-P", "TF-N"],
        })
        signed = pd.DataFrame({
            "gene": ["N1", "N2"],
            "Symbol": ["N1", "N2"],
            "Response_Direction": ["Influence Loss", "Influence Loss"],
            "Delta_PageRank_Pct": [-3.0, -9.0],
            "Lokalizasyon": ["Membrane", "Membrane"],
            "Gümrük_Kapisi": [True, False],
            "significant_redistribution": [True, False],
        })
        enrichment = pd.DataFrame({
            "Pathway_Family": ["GO", "KEGG"],
            "Kaynak": ["GO_Biological_Process_2021", "KEGG_2021_Human"],
            "Term": ["shared response process", "Acceptance pathway"],
            "Adjusted P-value": [0.01, 0.02],
            "Genes": ["T;P1;P2", "P1;P2"],
        })
        originals = tuple(frame.copy(deep=True) for frame in (report, signed, enrichment))
        context = self._context(
            report,
            signed,
            enrichment,
            identifiers=[("TARGET", "T"), ("P1", "P1"), ("P2", "P2"), ("N1", "N1"), ("N2", "N2")],
            targets=("TARGET",),
            selection_mode="null_fdr",
            top_n=20,
            with_families=True,
        )
        by_id = {item.entity.canonical_id: item.entity for item in context.entities}
        explicit = (
            EvidenceAssertion(
                assertion_id="fixture-ppi",
                subject=by_id["TARGET"],
                relationship_type=RelationshipType.PPI_EDGE,
                object=by_id["P1"],
                statement="A curated direct protein-protein interaction edge is recorded.",
                provenance=(_provenance("curated_ppi", ProvenanceKind.CURATED_DATABASE),),
            ),
            EvidenceAssertion(
                assertion_id="fixture-comention",
                subject=by_id["TARGET"],
                relationship_type=RelationshipType.CO_MENTIONED_IN_PUBLICATION,
                object=by_id["P1"],
                statement="The entities are co-mentioned in publication PMID:12345.",
                provenance=(_provenance("publication_12345", ProvenanceKind.LITERATURE),),
                qualifiers={"publication_id": "PMID:12345"},
            ),
        )
        relationships = build_target_response_relationships(context, explicit_evidence=explicit)
        return context, relationships, originals, (report, signed, enrichment)

    def test_all_and_only_seven_types_scope_family_function_gateway_and_fdr(self) -> None:
        context, relationships, *_ = self._rich_fixture()
        result = detect_observations(context, relationships, created_at=_CREATED_AT)

        self.assertEqual(len(ObservationType), 7)
        self.assertEqual({item.type for item in result.observations}, set(ObservationType))
        self.assertTrue(all(item.observation.scope.top_n == 20 for item in result.observations))
        self.assertTrue(all(item.observation.scope.selection_mode == "null_fdr" for item in result.observations))
        self.assertEqual(result.base_observations, tuple(item.observation for item in result.observations))

        family = result.of_type(ObservationType.FAMILY_COOCCURRENCE)
        self.assertEqual(len(family), 1)
        self.assertEqual(family[0].details["family_id"], "FAMILY-A")
        self.assertEqual(family[0].details["member_count"], 2)
        self.assertEqual(family[0].details["presentation_min_members"], 2)
        self.assertIs(family[0].details["presentation_only"], True)
        self.assertTrue(DEFAULT_RESEARCH_CONFIG.presentation_limits.family_threshold_presentation_only)
        self.assertNotIn("FAMILY-SINGLETON", repr(family))

        functional = result.of_type(ObservationType.FUNCTIONAL_COOCCURRENCE)
        self.assertGreaterEqual(len(functional), 2)
        self.assertEqual(
            {item.details["context_kind"] for item in functional},
            {"GO", "KEGG"},
        )
        for item in functional:
            language = f"{item.observation.basis} {item.observation.limitations}".casefold()
            self.assertTrue(any(phrase in language for phrase in (
                "observed among returned", "shared annotation", "frequently represented within the returned set",
            )))
            self.assertTrue(all(term not in language for term in (
                "enriched", "dominant biological mechanism", "activated", "suppressed",
            )))

        gateway = result.of_type(ObservationType.GATEWAY_RESPONSE_OVERLAP)
        self.assertEqual(len(gateway), 1)
        self.assertGreaterEqual(gateway[0].details["member_count"], 2)

        fdr = result.of_type(ObservationType.FDR_SUPPORTED_SUBSET)
        self.assertEqual(len(fdr), 1)
        self.assertEqual(fdr[0].details["selection_mode"], "null_fdr")
        self.assertEqual({member.canonical_id for member in fdr[0].members}, {"P1", "N1"})
        self.assertNotIn("P2", {member.canonical_id for member in fdr[0].members})

    def test_positive_negative_contrast_is_descriptive_only(self) -> None:
        context, relationships, *_ = self._rich_fixture()
        result = detect_observations(context, relationships, created_at=_CREATED_AT)
        contrast = result.of_type(ObservationType.POSITIVE_NEGATIVE_CONTRAST)

        self.assertEqual(len(contrast), 1)
        details = contrast[0].details
        self.assertEqual(details["positive_total"], 2)
        self.assertEqual(details["negative_total"], 2)
        self.assertTrue(details["annotations"])
        for row in details["annotations"]:
            self.assertTrue({
                "annotation", "positive_count", "negative_count",
                "positive_fraction", "negative_fraction",
            } <= set(row))
            self.assertAlmostEqual(row["positive_fraction"], row["positive_count"] / 2)
            self.assertAlmostEqual(row["negative_fraction"], row["negative_count"] / 2)
        nucleus = next(row for row in details["annotations"] if "nucleus" in row["annotation"].casefold())
        membrane = next(row for row in details["annotations"] if "membrane" in row["annotation"].casefold())
        self.assertEqual((nucleus["positive_count"], nucleus["negative_count"]), (2, 0))
        self.assertEqual((membrane["positive_count"], membrane["negative_count"]), (0, 2))
        forbidden_keys = {"p", "p_value", "q_value", "enrichment", "significance", "significant"}
        self.assertTrue(_nested_keys(details).isdisjoint(forbidden_keys))

    def test_relationship_types_remain_distinct_and_comention_is_not_interaction(self) -> None:
        context, relationships, *_ = self._rich_fixture()
        relationship_types = {item.relationship_type for item in relationships.assertions}
        self.assertTrue({
            RelationshipType.PPI_EDGE,
            RelationshipType.SHARED_GO,
            RelationshipType.CO_MENTIONED_IN_PUBLICATION,
        } <= relationship_types)

        ppi = next(item for item in relationships.assertions if "protein-protein" in item.statement)
        comention = next(item for item in relationships.assertions if "PMID:12345" in item.statement)
        self.assertIs(ppi.relationship_type, RelationshipType.PPI_EDGE)
        self.assertIs(comention.relationship_type, RelationshipType.CO_MENTIONED_IN_PUBLICATION)
        self.assertNotIn("interact", comention.statement.casefold())

        detected = detect_observations(context, relationships, created_at=_CREATED_AT)
        target_relationships = detected.of_type(ObservationType.TARGET_RESPONSE_RELATIONSHIP)
        self.assertTrue(target_relationships)
        detected_types = {
            assertion.relationship_type
            for item in target_relationships
            for assertion in item.assertions
        }
        self.assertTrue({
            RelationshipType.PPI_EDGE,
            RelationshipType.SHARED_GO,
            RelationshipType.CO_MENTIONED_IN_PUBLICATION,
        } <= detected_types)
        self.assertNotIn("related_to", {item.value for item in detected_types})

    def test_repeatability_immutability_required_fields_and_no_scores(self) -> None:
        context, relationships, originals, frames_now = self._rich_fixture()
        before = (
            context.entities,
            context.unresolved_entities,
            context.functional_context,
            context.roles,
            context.annotations,
            relationships.assertions,
        )
        first = detect_observations(context, relationships, created_at=_CREATED_AT)
        second = detect_observations(context, relationships, created_at=_CREATED_AT)

        self.assertEqual(first, second)
        self.assertEqual(
            [item.observation_id for item in first.observations],
            [item.observation_id for item in second.observations],
        )
        self.assertEqual(before, (
            context.entities,
            context.unresolved_entities,
            context.functional_context,
            context.roles,
            context.annotations,
            relationships.assertions,
        ))
        for original, current in zip(originals, frames_now):
            assert_frame_equal(original, current)

        required = {
            "observation_id", "simulation_id", "type", "members", "basis",
            "scope", "source_fields", "created_at", "limitations",
        }
        self.assertTrue(required <= {field.name for field in fields(Observation)})
        forbidden_field_fragments = ("score", "importance", "novelty", "evidence_strength", "priority")
        for item in first.observations:
            names = {field.name.casefold() for field in fields(item)} | _nested_keys(item.details)
            self.assertFalse(any(fragment in name for name in names for fragment in forbidden_field_fragments))
            language = f"{item.observation.basis} {item.observation.limitations}".casefold()
            self.assertFalse(any(fragment in language for fragment in ("score", "importance", "novelty")))
        with self.assertRaises(TypeError):
            first.observations[0].details["score"] = 1
        with self.assertRaises(FrozenInstanceError):
            first.observations[0].observation.basis = "changed"

    def test_fdr_requires_existing_null_fdr_role_not_stale_top_n_flag(self) -> None:
        report = pd.DataFrame({
            "gene": ["P1"],
            "Symbol": ["P1"],
            "Hasar_Tipi": ["Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [2.0],
            "significant_redistribution": [True],
            "q_value": [0.000001],
        })
        context = self._context(
            report,
            identifiers=[("P1", "P1")],
            selection_mode="top_n",
            top_n=20,
        )
        detected = detect_observations(context, created_at=_CREATED_AT)

        self.assertEqual(detected.of_type(ObservationType.FDR_SUPPORTED_SUBSET), ())
        self.assertNotIn(
            "fdr_supported",
            {annotation.annotation_type for annotation in context.entities[0].annotations},
        )

    def test_empty_and_old_report_shapes_are_safe(self) -> None:
        empty = self._context(
            pd.DataFrame(),
            identifiers=[],
            selection_mode="top_n",
            top_n=20,
        )
        empty_result = detect_observations(empty, created_at=_CREATED_AT)
        self.assertEqual(empty_result.observations, ())

        old_report = pd.DataFrame({
            "gene": ["OLD_TARGET", "OLD_RESPONSE"],
            "Delta_PageRank_Pct": [None, 1.0],
        })
        old = self._context(
            old_report,
            identifiers=[("OLD_TARGET", "OT"), ("OLD_RESPONSE", "OR")],
            targets=("OLD_TARGET",),
            selection_mode="top_n",
            top_n=20,
        )
        old_result = detect_observations(old, created_at=_CREATED_AT)
        self.assertIsInstance(old_result.observations, tuple)
        self.assertEqual(old_result.of_type(ObservationType.FDR_SUPPORTED_SUBSET), ())

    def test_200_entity_detector_budget(self) -> None:
        size = 200
        identifiers = [(f"E{index:03d}", f"S{index:03d}") for index in range(size)]
        report = pd.DataFrame({
            "gene": [canonical for canonical, _ in identifiers],
            "Symbol": [symbol for _, symbol in identifiers],
            "Hasar_Tipi": ["Üçüncül_Hasar_Stresli"] * size,
            "Delta_PageRank_Pct": [1.0] * size,
            "Lokalizasyon": ["Nucleus" if index % 2 else "Membrane" for index in range(size)],
            "Gümrük_Kapisi": [False] * size,
        })
        context = self._context(
            report,
            identifiers=identifiers,
            selection_mode="top_n",
            top_n=200,
        )

        started = time.perf_counter()
        detected = detect_observations(context, created_at=_CREATED_AT)
        elapsed = time.perf_counter() - started

        self.assertTrue(detected.of_type(ObservationType.LOCALIZATION_PATTERN))
        self.assertLess(elapsed, 2.0, f"200-entity detector took {elapsed:.3f}s")


if __name__ == "__main__":
    unittest.main()
