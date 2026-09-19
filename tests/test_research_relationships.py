"""Focused tests for evidence-only target/response relationships."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
import unittest

from src.research.contracts import CONTRACTS
from src.research.local_context import (
    LocalAnnotation,
    LocalEntityContext,
    LocalResearchContext,
)
from src.research.models import (
    Entity,
    EntityResolution,
    EntityType,
    EvidenceAssertion,
    ObservationScope,
    ProvenanceKind,
    ProvenanceRecord,
    RelationshipType,
    ResolutionStatus,
)
from src.research.providers import ComplexMembership, FamilyMembership, RegulatoryRelation
from src.research.relationships import build_target_response_relationships
from src.research.semantic_adapter import (
    FunctionalContextKind,
    FunctionalContextRow,
    ResponseDirection,
    ResultRole,
    ResultSource,
    SemanticRecord,
    SemanticSimulationResult,
)


def _provenance(
    source: str,
    kind: ProvenanceKind = ProvenanceKind.LOCAL_ANNOTATION,
) -> ProvenanceRecord:
    return ProvenanceRecord(
        source=source,
        kind=kind,
        version="fixture-v1",
        locator=f"fixture://{source}",
        snapshot_id="snapshot-relationships",
    )


def _resolution(entity: Entity) -> EntityResolution:
    return EntityResolution(
        query=entity.canonical_id,
        taxon_id=entity.taxon_id,
        entity_type=entity.entity_type,
        status=ResolutionStatus.EXACT,
        entity=entity,
        source="fixture",
    )


def _context(
    *,
    include_family: bool = True,
    include_complex: bool = True,
) -> tuple[LocalResearchContext, Entity, Entity]:
    snapshot_provenance = _provenance("sophiark_snapshot", ProvenanceKind.SOPHIARK_COMPUTED)
    annotation_provenance = _provenance("local_annotations")
    family_provenance = _provenance("hgnc", ProvenanceKind.CURATED_DATABASE)
    regulatory_provenance = _provenance("trrust", ProvenanceKind.CURATED_DATABASE)
    reactome_provenance = _provenance("reactome", ProvenanceKind.PATHWAY)

    target = Entity(9606, EntityType.PROTEIN, "ENSP_TARGET", "TG", "Target")
    response = Entity(9606, EntityType.PROTEIN, "ENSP_RESPONSE", "RG", "Response")
    target_record = SemanticRecord(
        source=ResultSource.REPORT,
        source_position=0,
        canonical_id=target.canonical_id,
        symbol=target.symbol,
        roles=(ResultRole.PERTURBATION_TARGET,),
        direction=ResponseDirection.UNAVAILABLE,
        values={"Topluluk_ID": 7},
    )
    response_record = SemanticRecord(
        source=ResultSource.REPORT,
        source_position=1,
        canonical_id=response.canonical_id,
        symbol=response.symbol,
        roles=(ResultRole.POSITIVE_REDISTRIBUTION,),
        direction=ResponseDirection.GAIN,
        values={"Community": 7},
    )

    target_annotations: list[LocalAnnotation] = [
        LocalAnnotation(
            "go_biological_process",
            {"id": "GO:0000001", "term": "fixture process"},
            "local GO",
            "GO Biyolojik Süreç",
            (annotation_provenance,),
        ),
        LocalAnnotation(
            "localization", "Nucleus", "Sophiark report", "Lokalizasyon",
            (snapshot_provenance,),
        ),
        LocalAnnotation(
            "reactome_pathway",
            {"reactome_id": "R-HSA-123", "pathway_name": "Fixture Reactome pathway"},
            "Reactome snapshot",
            provenance=(reactome_provenance,),
        ),
        LocalAnnotation(
            "directed_regulation",
            RegulatoryRelation(
                taxon_id=9606,
                source_symbol="TG",
                target_symbol="RG",
                database="TRRUST",
                effect="Activation",
                mechanism=None,
                directed=True,
                stimulation=True,
                inhibition=False,
                reference_ids=("PMID:123",),
                record_identifier="TRRUST:fixture-row",
                qualifiers={},
                provenance=regulatory_provenance,
            ),
            "TRRUST snapshot",
            provenance=(regulatory_provenance,),
        ),
    ]
    response_annotations: list[LocalAnnotation] = [
        LocalAnnotation(
            "go_biological_process",
            {"id": "GO:0000001", "term": "same ID, alternate label"},
            "local GO",
            "GO Biyolojik Süreç",
            (annotation_provenance,),
        ),
        LocalAnnotation(
            "localization", "nucleus", "Sophiark report", "Lokalizasyon",
            (snapshot_provenance,),
        ),
        LocalAnnotation(
            "reactome_context",
            {"id": "R-HSA-123", "name": "Same Reactome ID, alternate label"},
            "Reactome snapshot",
            provenance=(reactome_provenance,),
        ),
    ]

    if include_family:
        target_annotations.append(
            LocalAnnotation(
                "authoritative_family_membership",
                FamilyMembership(
                    9606, "TG", "TG", "HGNC:100", "Fixture family", "approved_symbol",
                    (), None, None, family_provenance,
                ),
                "HGNC",
                provenance=(family_provenance,),
            )
        )
        response_annotations.append(
            LocalAnnotation(
                "authoritative_family_membership",
                FamilyMembership(
                    9606, "RG", "RG", "HGNC:100", "UPDATED DISPLAY NAME",
                    "approved_symbol", (), None, None, family_provenance,
                ),
                "HGNC",
                provenance=(family_provenance,),
            )
        )

    if include_complex:
        for annotations, entity in (
            (target_annotations, target),
            (response_annotations, response),
        ):
            annotations.append(
                LocalAnnotation(
                    "complex_membership",
                    ComplexMembership(
                        9606,
                        entity.canonical_id,
                        "ensembl_protein",
                        "CPX-1",
                        "Fixture complex",
                        annotation_provenance,
                    ),
                    "Complex Portal",
                    provenance=(annotation_provenance,),
                )
            )

    kegg_row = FunctionalContextRow(
        kind=FunctionalContextKind.KEGG,
        source_position=0,
        contract=CONTRACTS["kegg_enrichment"],
        term="Fixture KEGG pathway",
        source_library="KEGG_2021_Human",
        pathway_family="KEGG",
        members=(target.canonical_id, response.canonical_id),
        values={"Genes": f"{target.canonical_id};{response.canonical_id}"},
    )
    scope = ObservationScope(
        selection_mode="top_n",
        tested_count=2,
        returned_count=1,
        threshold=None,
        top_n=1,
        tissue="Pancreas",
        species="Homo sapiens",
        taxon_id=9606,
    )
    semantic = SemanticSimulationResult(
        snapshot_id="snapshot-relationships",
        simulation_id="simulation-relationships",
        scope=scope,
        target_ids=(target.canonical_id,),
        target_records=(target_record,),
        positive_redistribution=(response_record,),
        network_losses=(),
        fdr_supported=(),
        go_terms=(),
        kegg_pathways=(kegg_row,),
        unclassified_enrichment=(),
    )
    target_context = LocalEntityContext(
        entity=target,
        resolution=_resolution(target),
        roles=(ResultRole.PERTURBATION_TARGET,),
        directions=(ResponseDirection.UNAVAILABLE,),
        metrics=(),
        annotations=tuple(target_annotations),
        source_records=(target_record,),
        provenance=(snapshot_provenance,),
    )
    response_context = LocalEntityContext(
        entity=response,
        resolution=_resolution(response),
        roles=(ResultRole.POSITIVE_REDISTRIBUTION,),
        directions=(ResponseDirection.GAIN,),
        metrics=(),
        annotations=tuple(response_annotations),
        source_records=(response_record,),
        provenance=(snapshot_provenance,),
    )
    context = LocalResearchContext(
        snapshot_id=semantic.snapshot_id,
        simulation_id=semantic.simulation_id,
        semantic_result=semantic,
        entities=(target_context, response_context),
        unresolved_entities=(),
        roles={
            target.logical_key: target_context.roles,
            response.logical_key: response_context.roles,
        },
        annotations={
            target.logical_key: target_context.annotations,
            response.logical_key: response_context.annotations,
        },
        functional_context=(kegg_row,),
        source_statuses=(),
        provenance=(snapshot_provenance,),
    )
    return context, target, response


class RelationshipBuilderTests(unittest.TestCase):
    def test_local_evidence_types_are_exact_and_reactome_is_supported(self) -> None:
        context, _, _ = _context()

        result = build_target_response_relationships(context)

        self.assertEqual(
            {item.relationship_type for item in result.assertions},
            {
                RelationshipType.DIRECTED_REGULATION,
                RelationshipType.SHARED_PATHWAY,
                RelationshipType.SHARED_GO,
                RelationshipType.SHARED_FAMILY,
                RelationshipType.SHARED_COMPARTMENT,
                RelationshipType.SHARED_COMMUNITY,
            },
        )
        pathways = result.of_type(RelationshipType.SHARED_PATHWAY)
        self.assertEqual(len(pathways), 2)
        self.assertEqual(
            {item.qualifiers["basis"] for item in pathways},
            {"exact_shared_reactome_annotation", "existing_functional_context_membership"},
        )
        family = result.of_type(RelationshipType.SHARED_FAMILY)
        self.assertEqual(len(family), 1)
        self.assertEqual(family[0].value["family_id"], "HGNC:100")
        directed = result.of_type(RelationshipType.DIRECTED_REGULATION)
        self.assertEqual((directed[0].subject.symbol, directed[0].object.symbol), ("TG", "RG"))
        self.assertTrue(all(item.provenance for item in result.assertions))
        self.assertFalse(
            {
                RelationshipType.PPI_EDGE,
                RelationshipType.NETWORK_PATH,
                RelationshipType.CO_MENTIONED_IN_PUBLICATION,
            }
            & {item.relationship_type for item in result.assertions}
        )

    def test_complex_membership_is_never_converted_to_family(self) -> None:
        context, _, _ = _context(include_family=False, include_complex=True)

        result = build_target_response_relationships(context)

        self.assertEqual(result.of_type(RelationshipType.SHARED_FAMILY), ())

    def test_ppi_go_path_and_publication_types_remain_distinct(self) -> None:
        context, target, response = _context()
        provenance = _provenance("explicit", ProvenanceKind.CURATED_DATABASE)
        literature = _provenance("publication", ProvenanceKind.LITERATURE)
        # The alternate display symbol proves logical-key matching is followed
        # by canonicalisation to the frozen LocalResearchContext entity.
        equivalent_target = Entity(
            target.taxon_id,
            target.entity_type,
            target.canonical_id,
            "ALTERNATE",
        )
        explicit = (
            EvidenceAssertion(
                "input-ppi",
                equivalent_target,
                RelationshipType.PPI_EDGE,
                response,
                "A curated PPI edge exists for this exact pair.",
                (provenance,),
            ),
            EvidenceAssertion(
                "input-path",
                target,
                RelationshipType.NETWORK_PATH,
                response,
                "An already-computed current-network path exists for this exact pair.",
                (provenance,),
            ),
            EvidenceAssertion(
                "input-publication",
                target,
                RelationshipType.CO_MENTIONED_IN_PUBLICATION,
                response,
                "Both entities are mentioned in PMID:123.",
                (literature,),
                qualifiers={"publication_id": "PMID:123"},
            ),
        )
        before = (context.entities, context.annotations, context.functional_context, explicit)

        first = build_target_response_relationships(context, explicit_evidence=explicit)
        second = build_target_response_relationships(context, explicit_evidence=explicit)

        self.assertEqual(first, second)
        self.assertEqual(
            before,
            (context.entities, context.annotations, context.functional_context, explicit),
        )
        self.assertEqual(len(first.of_type(RelationshipType.PPI_EDGE)), 1)
        self.assertEqual(len(first.of_type(RelationshipType.SHARED_GO)), 1)
        self.assertEqual(len(first.of_type(RelationshipType.NETWORK_PATH)), 1)
        self.assertEqual(len(first.of_type(RelationshipType.CO_MENTIONED_IN_PUBLICATION)), 1)
        self.assertIs(first.of_type(RelationshipType.PPI_EDGE)[0].subject, target)
        self.assertNotEqual(first.of_type(RelationshipType.PPI_EDGE)[0].assertion_id, "input-ppi")
        self.assertEqual(
            [item.assertion_id for item in first.assertions],
            [item.assertion_id for item in second.assertions],
        )
        self.assertNotIn(
            RelationshipType.PPI_EDGE,
            {
                item.relationship_type
                for item in first.assertions
                if "PMID:123" in item.statement
            },
        )
        with self.assertRaises(FrozenInstanceError):
            first.notices = ()

    def test_cross_taxon_and_out_of_pair_explicit_evidence_is_rejected(self) -> None:
        context, target, response = _context()
        provenance = _provenance("explicit", ProvenanceKind.CURATED_DATABASE)
        mouse = Entity(10090, EntityType.PROTEIN, "ENSMUSP_MOUSE", "MouseGene")
        unrelated = Entity(9606, EntityType.PROTEIN, "ENSP_UNRELATED", "Other")
        explicit = (
            EvidenceAssertion(
                "mouse-edge",
                target,
                RelationshipType.PPI_EDGE,
                mouse,
                "Cross-taxon fixture edge.",
                (provenance,),
            ),
            EvidenceAssertion(
                "unrelated-edge",
                target,
                RelationshipType.PPI_EDGE,
                unrelated,
                "Out-of-pair fixture edge.",
                (provenance,),
            ),
            EvidenceAssertion(
                "valid-reversed-edge",
                response,
                RelationshipType.PPI_EDGE,
                target,
                "Same-taxon exact response-to-target fixture edge.",
                (provenance,),
            ),
        )

        result = build_target_response_relationships(context, explicit_evidence=explicit)

        self.assertEqual(result.rejected_explicit_count, 2)
        self.assertEqual(len(result.of_type(RelationshipType.PPI_EDGE)), 1)
        self.assertEqual(result.of_type(RelationshipType.PPI_EDGE)[0].subject, response)
        self.assertTrue(any("Rejected 2" in notice for notice in result.notices))


if __name__ == "__main__":
    unittest.main()
