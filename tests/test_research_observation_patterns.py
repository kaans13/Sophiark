"""Focused tests for all seven deterministic Research Context observations."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, fields, replace
import inspect
from types import SimpleNamespace
import unittest

from src.research.config import PresentationLimits, ResearchConfig
from src.research.contracts import CONTRACTS
from src.research.local_context import LocalAnnotation, LocalEntityContext, LocalResearchContext
from src.research.models import (
    Entity,
    EntityResolution,
    EntityType,
    EvidenceAssertion,
    Observation,
    ObservationScope,
    ObservationType,
    ProvenanceKind,
    ProvenanceRecord,
    RelationshipType,
    ResolutionStatus,
)
from src.research.observations import detect_observations
from src.research.providers import FamilyMembership
from src.research.semantic_adapter import (
    FunctionalContextKind,
    FunctionalContextRow,
    ResponseDirection,
    ResultRole,
    SemanticResultView,
)


class ObservationPatternTests(unittest.TestCase):
    @staticmethod
    def _provenance(source: str = "fixture") -> ProvenanceRecord:
        return ProvenanceRecord(
            source=source,
            kind=ProvenanceKind.LOCAL_ANNOTATION,
            retrieved_at="2026-01-02T03:04:05Z",
            locator=f"fixture://{source}",
        )

    @classmethod
    def _family(cls, symbol: str, group_id: str = "HGNC:42") -> FamilyMembership:
        return FamilyMembership(
            taxon_id=9606,
            query_identifier=symbol,
            approved_symbol=symbol,
            group_id=group_id,
            group_name="Fixture receptor family",
            matched_via="approved_symbol",
            aliases=(),
            ensembl_gene_id=None,
            ncbi_gene_id=None,
            provenance=cls._provenance("hgnc"),
        )

    @classmethod
    def _entity_context(
        cls,
        canonical_id: str,
        symbol: str,
        roles: tuple[ResultRole, ...],
        annotations: tuple[LocalAnnotation, ...],
    ) -> LocalEntityContext:
        entity = Entity(9606, EntityType.PROTEIN, canonical_id, symbol=symbol)
        resolution = EntityResolution(
            query=canonical_id,
            taxon_id=9606,
            entity_type=EntityType.PROTEIN,
            status=ResolutionStatus.EXACT,
            entity=entity,
            source="fixture_resolver",
        )
        return LocalEntityContext(
            entity=entity,
            resolution=resolution,
            roles=roles,
            directions=(
                ResponseDirection.GAIN
                if ResultRole.POSITIVE_REDISTRIBUTION in roles
                else ResponseDirection.LOSS
                if ResultRole.NETWORK_LOSS in roles
                else ResponseDirection.UNAVAILABLE,
            ),
            metrics=(),
            annotations=annotations,
            source_records=(),
            provenance=(cls._provenance(),),
        )

    @classmethod
    def _context(cls, *, selection_mode: str = "null_fdr") -> LocalResearchContext:
        provenance = cls._provenance("snapshot")
        family_p1 = LocalAnnotation(
            "authoritative_family_membership", cls._family("P1"),
            "local_hgnc_family", "group_id/group_name", (provenance,),
        )
        family_p2 = LocalAnnotation(
            "authoritative_family_membership", cls._family("P2"),
            "local_hgnc_family", "group_id/group_name", (provenance,),
        )
        # The same annotation name with a non-authoritative value must not be
        # accepted as family evidence.
        guessed_family = LocalAnnotation(
            "authoritative_family_membership", "symbol-prefix guess",
            "fixture", "guess", (provenance,),
        )
        go_p1 = LocalAnnotation(
            "go_biological_process", ("calcium signaling",),
            "local_mygene", "go_bp", (provenance,),
        )
        go_p2 = LocalAnnotation(
            "go_biological_process", ("calcium signaling",),
            "local_mygene", "go_bp", (provenance,),
        )
        kegg_row = FunctionalContextRow(
            kind=FunctionalContextKind.KEGG,
            source_position=0,
            contract=CONTRACTS["kegg_enrichment"],
            term="Fixture signaling pathway",
            source_library="KEGG_2021_Human",
            pathway_family="KEGG",
            members=("P1", "N1"),
            values={"Term": "Fixture signaling pathway", "Genes": "P1;N1"},
        )
        kegg_p1 = LocalAnnotation(
            "kegg_enrichment", kegg_row, "sophiark_enrichment", "Genes", (provenance,),
        )
        kegg_n1 = LocalAnnotation(
            "kegg_enrichment", kegg_row, "sophiark_enrichment", "Genes", (provenance,),
        )
        reactome_p2 = LocalAnnotation(
            "reactome_pathway", ("R-HSA-123 Fixture pathway",),
            "local_reactome", "pathway", (provenance,),
        )
        reactome_n1 = LocalAnnotation(
            "reactome_pathway", ("R-HSA-123 Fixture pathway",),
            "local_reactome", "pathway", (provenance,),
        )

        def local(value: str) -> LocalAnnotation:
            return LocalAnnotation("localization", value, "sophiark_report", "Lokalizasyon", (provenance,))

        def gateway(value: bool) -> LocalAnnotation:
            return LocalAnnotation("gateway", value, "sophiark_report", "Gümrük_Kapisi", (provenance,))

        target = cls._entity_context(
            "ENSP_T", "T", (ResultRole.PERTURBATION_TARGET,),
            (local("Membrane"), gateway(True)),
        )
        p1 = cls._entity_context(
            "ENSP_P1", "P1",
            (ResultRole.POSITIVE_REDISTRIBUTION, ResultRole.FDR_SUPPORTED),
            (family_p1, go_p1, kegg_p1, local("Nucleus"), gateway(True)),
        )
        p2 = cls._entity_context(
            "ENSP_P2", "P2", (ResultRole.POSITIVE_REDISTRIBUTION,),
            (family_p2, go_p2, reactome_p2, local("Nucleus"), gateway(False)),
        )
        n1 = cls._entity_context(
            "ENSP_N1", "N1", (ResultRole.NETWORK_LOSS,),
            (kegg_n1, reactome_n1, local("Cytosol"), gateway(False)),
        )
        n2 = cls._entity_context(
            "ENSP_N2", "N2", (ResultRole.NETWORK_LOSS, ResultRole.FDR_SUPPORTED),
            (guessed_family, local("Cytosol"), gateway(True)),
        )
        scope = ObservationScope(
            selection_mode=selection_mode,
            tested_count=200,
            returned_count=4,
            threshold=0.05,
            top_n=20,
            tissue="Pancreas",
            species="Homo sapiens",
            taxon_id=9606,
            threshold_parameters={"fdr_alpha": 0.05},
        )
        semantic = SemanticResultView(
            snapshot_id="snapshot-1",
            simulation_id="simulation-1",
            scope=scope,
            target_ids=("ENSP_T",),
            target_records=(),
            positive_redistribution=(),
            network_losses=(),
            fdr_supported=(),
            go_terms=(),
            kegg_pathways=(),
            unclassified_enrichment=(),
        )
        entities = (target, p1, p2, n1, n2)
        return LocalResearchContext(
            snapshot_id="snapshot-1",
            simulation_id="simulation-1",
            semantic_result=semantic,
            entities=entities,
            unresolved_entities=(),
            roles={item.entity.logical_key: item.roles for item in entities},
            annotations={item.entity.logical_key: item.annotations for item in entities},
            functional_context=(kegg_row,),
            source_statuses=(),
            provenance=(provenance,),
        )

    @classmethod
    def _assertions(cls, context: LocalResearchContext):
        by_id = {item.entity.canonical_id: item.entity for item in context.entities}
        provenance = (cls._provenance("relationships"),)
        return (
            EvidenceAssertion(
                "co-mention", by_id["ENSP_T"], RelationshipType.CO_MENTIONED_IN_PUBLICATION,
                by_id["ENSP_P1"], "Co-mentioned in one publication.", provenance,
            ),
            EvidenceAssertion(
                "shared-go", by_id["ENSP_T"], RelationshipType.SHARED_GO,
                by_id["ENSP_P1"], "Shared GO annotation.", provenance,
            ),
            EvidenceAssertion(
                "ppi", by_id["ENSP_T"], RelationshipType.PPI_EDGE,
                by_id["ENSP_P1"], "Existing PPI edge.", provenance,
            ),
            EvidenceAssertion(
                "directed", by_id["ENSP_N1"], RelationshipType.DIRECTED_REGULATION,
                by_id["ENSP_T"], "Directed response-to-target record.", provenance,
            ),
            EvidenceAssertion(
                "not-target-response", by_id["ENSP_P1"], RelationshipType.SHARED_FAMILY,
                by_id["ENSP_P2"], "Response-response assertion.", provenance,
            ),
        )

    def test_all_seven_v1_types_are_detected_without_scores(self) -> None:
        context = self._context()
        result = detect_observations(context, self._assertions(context))
        self.assertEqual(
            {item.type for item in result.observations},
            set(ObservationType),
        )
        self.assertTrue(all(item.observation.scope is context.semantic_result.scope for item in result.observations))
        self.assertNotIn("score", {item.name for item in fields(Observation)})
        self.assertFalse(any("score" in item.details for item in result.observations))

    def test_family_requires_authoritative_membership_and_presentation_threshold(self) -> None:
        context = self._context()
        detected = detect_observations(context)
        families = detected.of_type(ObservationType.FAMILY_COOCCURRENCE)
        self.assertEqual(len(families), 1)
        self.assertEqual([item.canonical_id for item in families[0].members], ["ENSP_P1", "ENSP_P2"])
        self.assertEqual(families[0].details["family_id"], "HGNC:42")
        self.assertTrue(families[0].details["presentation_only"])
        self.assertNotIn("significant", families[0].basis.casefold())

        config = ResearchConfig(presentation_limits=PresentationLimits(family_min_members=3))
        self.assertEqual(
            detect_observations(context, config=config).of_type(ObservationType.FAMILY_COOCCURRENCE),
            (),
        )

    def test_functional_cooccurrence_keeps_go_kegg_reactome_distinct(self) -> None:
        functional = detect_observations(self._context()).of_type(
            ObservationType.FUNCTIONAL_COOCCURRENCE
        )
        self.assertEqual(
            [item.details["context_kind"] for item in functional],
            ["GO", "KEGG", "Reactome"],
        )
        self.assertEqual(
            [item.details["term"] for item in functional],
            ["calcium signaling", "Fixture signaling pathway", "R-HSA-123 Fixture pathway"],
        )
        self.assertTrue(all("activated" not in item.basis.casefold() for item in functional))
        self.assertTrue(all("dominant" not in item.basis.casefold() for item in functional))

    def test_localization_is_a_descriptive_distribution(self) -> None:
        observation = detect_observations(self._context()).of_type(
            ObservationType.LOCALIZATION_PATTERN
        )[0]
        self.assertEqual(dict(observation.details["distribution"]), {"Nucleus": 2, "Cytosol": 2})
        self.assertEqual(observation.details["total_count"], 4)
        self.assertEqual(observation.details["returned_response_count"], 4)
        self.assertEqual(dict(observation.details["fractions_of_annotated_members"]), {
            "Nucleus": 0.5, "Cytosol": 0.5,
        })

    def test_positive_negative_contrast_has_counts_and_fractions_only(self) -> None:
        observation = detect_observations(self._context()).of_type(
            ObservationType.POSITIVE_NEGATIVE_CONTRAST
        )[0]
        rows = tuple(observation.details["annotations"])
        nucleus = next(row for row in rows if row["annotation"] == "Nucleus")
        cytosol = next(row for row in rows if row["annotation"] == "Cytosol")
        self.assertEqual(
            (nucleus["positive_count"], nucleus["negative_count"], nucleus["positive_fraction"], nucleus["negative_fraction"]),
            (2, 0, 1.0, 0.0),
        )
        self.assertEqual(
            (cytosol["positive_count"], cytosol["negative_count"], cytosol["positive_fraction"], cytosol["negative_fraction"]),
            (0, 2, 0.0, 1.0),
        )
        forbidden = {"p_value", "q_value", "significance", "score", "rank"}
        self.assertTrue(all(forbidden.isdisjoint(row) for row in rows))
        self.assertIn("No p-value", observation.limitations)

    def test_gateway_and_fdr_use_only_preserved_roles(self) -> None:
        result = detect_observations(self._context())
        gateway = result.of_type(ObservationType.GATEWAY_RESPONSE_OVERLAP)[0]
        overlaps = gateway.details["overlap_roles"]
        self.assertEqual(overlaps[ResultRole.POSITIVE_REDISTRIBUTION.value], ("9606:protein:ENSP_P1",))
        self.assertEqual(overlaps[ResultRole.NETWORK_LOSS.value], ("9606:protein:ENSP_N2",))
        self.assertEqual(overlaps[ResultRole.FDR_SUPPORTED.value], (
            "9606:protein:ENSP_P1", "9606:protein:ENSP_N2",
        ))
        fdr = result.of_type(ObservationType.FDR_SUPPORTED_SUBSET)[0]
        self.assertEqual([item.canonical_id for item in fdr.members], ["ENSP_P1", "ENSP_N2"])
        self.assertEqual(fdr.details["selection_mode"], "null_fdr")

        top_n = detect_observations(self._context(selection_mode="top_n"))
        self.assertEqual(top_n.of_type(ObservationType.FDR_SUPPORTED_SUBSET), ())
        top_gateway = top_n.of_type(ObservationType.GATEWAY_RESPONSE_OVERLAP)[0]
        self.assertNotIn(ResultRole.FDR_SUPPORTED.value, top_gateway.details["overlap_roles"])

    def test_relationship_assertion_types_and_orientation_are_preserved(self) -> None:
        context = self._context()
        assertions = self._assertions(context)
        wrapper = SimpleNamespace(
            assertions=assertions,
            notices=("relationship fixture notice",),
            rejected_explicit_count=1,
        )
        result = detect_observations(context, wrapper)
        relationships = result.of_type(ObservationType.TARGET_RESPONSE_RELATIONSHIP)
        self.assertEqual(len(relationships), 2)
        p1 = next(item for item in relationships if item.details["response"].endswith("ENSP_P1"))
        self.assertEqual(
            tuple(item.relationship_type for item in p1.assertions),
            (
                RelationshipType.PPI_EDGE,
                RelationshipType.SHARED_GO,
                RelationshipType.CO_MENTIONED_IN_PUBLICATION,
            ),
        )
        self.assertEqual(p1.details["relationship_types"], (
            "PPI_EDGE", "SHARED_GO", "CO_MENTIONED_IN_PUBLICATION",
        ))
        directed = next(item for item in relationships if item.details["response"].endswith("ENSP_N1"))
        self.assertIs(directed.assertions[0].relationship_type, RelationshipType.DIRECTED_REGULATION)
        self.assertEqual(directed.assertions[0].subject.canonical_id, "ENSP_N1")
        self.assertEqual(result.rejected_relationship_count, 2)
        self.assertIn("relationship fixture notice", result.notices)

    def test_ids_are_stable_but_created_at_remains_explicit(self) -> None:
        context = self._context()
        assertions = self._assertions(context)
        first = detect_observations(
            context, assertions, created_at="2026-02-01T00:00:00Z",
        )
        second = detect_observations(
            context, reversed(assertions), created_at="2026-03-01T00:00:00Z",
        )
        first_ids = {item.observation_id for item in first.observations}
        second_ids = {item.observation_id for item in second.observations}
        self.assertEqual(first_ids, second_ids)
        self.assertTrue(all(item.created_at == "2026-02-01T00:00:00+00:00" for item in first.observations))
        from_provenance = detect_observations(context)
        self.assertTrue(all(item.created_at == "2026-01-02T03:04:05+00:00" for item in from_provenance.observations))

    def test_result_and_nested_details_are_immutable(self) -> None:
        result = detect_observations(self._context())
        family = result.of_type(ObservationType.FAMILY_COOCCURRENCE)[0]
        with self.assertRaises(FrozenInstanceError):
            result.observations = ()  # type: ignore[misc]
        with self.assertRaises(TypeError):
            family.details["member_count"] = 99  # type: ignore[index]
        contrast = result.of_type(ObservationType.POSITIVE_NEGATIVE_CONTRAST)[0]
        with self.assertRaises(TypeError):
            contrast.details["annotations"][0]["positive_count"] = 99  # type: ignore[index]

    def test_empty_context_is_safe_and_does_not_fabricate_patterns(self) -> None:
        context = self._context()
        empty = replace(
            context,
            entities=(), roles={}, annotations={}, functional_context=(),
        )
        result = detect_observations(empty)
        self.assertEqual(result.observations, ())
        self.assertTrue(any("No resolved returned" in notice for notice in result.notices))

    def test_module_has_no_core_session_file_or_network_dependency(self) -> None:
        import src.research.observations as module

        tree = ast.parse(inspect.getsource(module))
        imported_roots = {
            alias.name.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertTrue(imported_roots.isdisjoint({
            "streamlit", "requests", "httpx", "urllib", "socket",
            "pandas", "numpy", "igraph", "networkx",
        }))
        source = inspect.getsource(module)
        self.assertNotIn("open(", source)
        self.assertNotIn("read_csv", source)
        self.assertNotIn("pagerank", source.casefold())
        self.assertNotIn("benjamini", source.casefold())


if __name__ == "__main__":
    unittest.main()
