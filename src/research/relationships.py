"""Deterministic target/response relationship assertions.

This module consumes only an immutable :class:`LocalResearchContext` plus
explicitly supplied, already-computed relationship evidence.  It never reads a
graph, computes a path, opens a file, inspects session state, calls a provider,
or performs network I/O.
"""

from __future__ import annotations

from dataclasses import dataclass, fields, is_dataclass
from enum import Enum
import hashlib
import json
from typing import Any, Iterable, Mapping

from .local_context import LocalAnnotation, LocalEntityContext, LocalResearchContext
from .models import (
    Entity,
    EvidenceAssertion,
    ProvenanceRecord,
    RelationshipType,
    immutable_mapping,
)
from .providers import FamilyMembership, RegulatoryRelation
from .semantic_adapter import FunctionalContextKind, FunctionalContextRow, ResultRole


_RESPONSE_ROLES = frozenset(
    {
        ResultRole.POSITIVE_REDISTRIBUTION,
        ResultRole.NETWORK_LOSS,
        ResultRole.FDR_SUPPORTED,
    }
)

_RELATIONSHIP_ORDER = {
    relationship_type: position
    for position, relationship_type in enumerate(
        (
            RelationshipType.PPI_EDGE,
            RelationshipType.DIRECTED_REGULATION,
            RelationshipType.SHARED_PATHWAY,
            RelationshipType.SHARED_GO,
            RelationshipType.SHARED_FAMILY,
            RelationshipType.SHARED_COMPARTMENT,
            RelationshipType.SHARED_COMMUNITY,
            RelationshipType.NETWORK_PATH,
            RelationshipType.CO_MENTIONED_IN_PUBLICATION,
        )
    )
}

_GO_ANNOTATION_TYPES = frozenset(
    {
        "go_biological_process",
        "go_cellular_component",
        "go_molecular_function",
    }
)

_PATHWAY_ANNOTATION_TYPES = frozenset(
    {
        "reactome_pathway",
        "reactome_context",
    }
)

_COMMUNITY_FIELDS = (
    "Topluluk_ID",
    "Community",
    "community",
    "Community_Context",
)


@dataclass(frozen=True, slots=True)
class RelationshipBuildResult:
    simulation_id: str
    snapshot_id: str
    assertions: tuple[EvidenceAssertion, ...]
    notices: tuple[str, ...] = ()
    rejected_explicit_count: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "assertions", tuple(self.assertions))
        object.__setattr__(self, "notices", tuple(self.notices))
        if not self.simulation_id or not self.snapshot_id:
            raise ValueError("relationship result requires simulation and snapshot IDs")
        if self.rejected_explicit_count < 0:
            raise ValueError("rejected_explicit_count cannot be negative")

    def of_type(self, relationship_type: RelationshipType | str) -> tuple[EvidenceAssertion, ...]:
        requested = RelationshipType(relationship_type)
        return tuple(
            assertion
            for assertion in self.assertions
            if assertion.relationship_type is requested
        )


def _stable(value: Any) -> Any:
    if isinstance(value, Enum):
        return value.value
    if isinstance(value, Entity):
        return value.logical_key
    if isinstance(value, Mapping):
        return {
            str(key): _stable(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (tuple, list)):
        return [_stable(item) for item in value]
    if isinstance(value, (set, frozenset)):
        return sorted((_stable(item) for item in value), key=repr)
    if is_dataclass(value):
        return {
            field.name: _stable(getattr(value, field.name))
            for field in fields(value)
        }
    if isinstance(value, (str, int, float, bool)) or value is None:
        return value
    return str(value)


def _deterministic_id(
    *,
    simulation_id: str,
    relationship_type: RelationshipType,
    subject: Entity,
    object_entity: Entity,
    basis: Any,
) -> str:
    payload = {
        "simulation_id": simulation_id,
        "relationship_type": relationship_type.value,
        "subject": subject.logical_key,
        "object": object_entity.logical_key,
        "basis": _stable(basis),
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")
    return f"relationship:{hashlib.sha256(encoded).hexdigest()[:24]}"


def _provenance_key(record: ProvenanceRecord) -> str:
    return json.dumps(_stable(record), ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _combine_provenance(
    *groups: Iterable[ProvenanceRecord],
    fallback: Iterable[ProvenanceRecord] = (),
) -> tuple[ProvenanceRecord, ...]:
    unique: dict[str, ProvenanceRecord] = {}
    for group in groups:
        for record in group:
            if isinstance(record, ProvenanceRecord):
                unique[_provenance_key(record)] = record
    if not unique:
        for record in fallback:
            if isinstance(record, ProvenanceRecord):
                unique[_provenance_key(record)] = record
    return tuple(unique[key] for key in sorted(unique))


def _entity_identifiers(context: LocalEntityContext) -> frozenset[str]:
    values = {
        context.entity.canonical_id,
        context.entity.symbol,
        *(record.canonical_id for record in context.source_records),
        *(record.symbol for record in context.source_records),
    }
    return frozenset(
        str(value).strip().casefold()
        for value in values
        if value is not None and str(value).strip()
    )


def _annotations(context: LocalEntityContext, annotation_type: str) -> tuple[LocalAnnotation, ...]:
    return tuple(
        annotation
        for annotation in context.annotations
        if annotation.annotation_type == annotation_type
    )


def _items(value: Any) -> tuple[Any, ...]:
    if isinstance(value, (tuple, list, set, frozenset)):
        return tuple(value)
    return (value,)


def _term_identity(value: Any) -> tuple[str, str] | None:
    if isinstance(value, Mapping):
        identifier = (
            value.get("id")
            or value.get("go_id")
            or value.get("reactome_id")
            or value.get("pathway_id")
            or value.get("term_id")
        )
        label = (
            value.get("term")
            or value.get("name")
            or value.get("pathway_name")
            or identifier
        )
        if identifier is not None and str(identifier).strip():
            return f"id:{str(identifier).strip().casefold()}", str(label or identifier).strip()
        if label is not None and str(label).strip():
            return f"term:{str(label).strip().casefold()}", str(label).strip()
        return None
    text = str(value or "").strip()
    if not text or text.casefold() in {"nan", "none", "unknown", "—"}:
        return None
    return f"term:{text.casefold()}", text


def _annotation_provenance(
    annotation: LocalAnnotation,
    entity_context: LocalEntityContext,
    research_context: LocalResearchContext,
) -> tuple[ProvenanceRecord, ...]:
    return _combine_provenance(
        annotation.provenance,
        fallback=(*entity_context.provenance, *research_context.provenance),
    )


def _make_assertion(
    *,
    context: LocalResearchContext,
    relationship_type: RelationshipType,
    subject: Entity,
    object_entity: Entity,
    statement: str,
    provenance: tuple[ProvenanceRecord, ...],
    value: Any,
    qualifiers: Mapping[str, Any],
    basis: Any,
) -> EvidenceAssertion | None:
    if not provenance:
        return None
    return EvidenceAssertion(
        assertion_id=_deterministic_id(
            simulation_id=context.simulation_id,
            relationship_type=relationship_type,
            subject=subject,
            object_entity=object_entity,
            basis=basis,
        ),
        subject=subject,
        relationship_type=relationship_type,
        object=object_entity,
        statement=statement,
        provenance=provenance,
        value=value,
        qualifiers=qualifiers,
    )


def _family_assertions(
    context: LocalResearchContext,
    target: LocalEntityContext,
    response: LocalEntityContext,
) -> tuple[EvidenceAssertion, ...]:
    target_rows: dict[tuple[str, str], list[tuple[FamilyMembership, LocalAnnotation]]] = {}
    response_rows: dict[tuple[str, str], list[tuple[FamilyMembership, LocalAnnotation]]] = {}
    for entity_context, output in ((target, target_rows), (response, response_rows)):
        for annotation in _annotations(entity_context, "authoritative_family_membership"):
            for value in _items(annotation.value):
                if not isinstance(value, FamilyMembership) or value.taxon_id != entity_context.entity.taxon_id:
                    continue
                group_id = value.group_id.strip()
                group_name = value.group_name.strip()
                key = (
                    ("id", group_id.casefold())
                    if group_id
                    else ("name", group_name.casefold())
                )
                output.setdefault(key, []).append((value, annotation))

    assertions: list[EvidenceAssertion] = []
    for key in sorted(set(target_rows) & set(response_rows)):
        left, right = target_rows[key][0], response_rows[key][0]
        provenance = _combine_provenance(
            _annotation_provenance(left[1], target, context),
            _annotation_provenance(right[1], response, context),
            (left[0].provenance, right[0].provenance),
        )
        basis = {"source": "HGNC", "group_id": left[0].group_id, "group_name": left[0].group_name}
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.SHARED_FAMILY,
            subject=target.entity,
            object_entity=response.entity,
            statement="Target and returned response share an authoritative HGNC family annotation.",
            provenance=provenance,
            value=immutable_mapping(
                {"family_id": left[0].group_id, "family_name": left[0].group_name}
            ),
            qualifiers=immutable_mapping({"basis": "exact_authoritative_family_membership"}),
            basis=basis,
        )
        if assertion is not None:
            assertions.append(assertion)
    return tuple(assertions)


def _go_annotation_entries(
    context: LocalResearchContext,
    entity_context: LocalEntityContext,
) -> Mapping[tuple[str, str], tuple[str, LocalAnnotation, tuple[ProvenanceRecord, ...]]]:
    entries: dict[tuple[str, str], tuple[str, LocalAnnotation, tuple[ProvenanceRecord, ...]]] = {}
    for annotation in entity_context.annotations:
        if annotation.annotation_type not in _GO_ANNOTATION_TYPES:
            continue
        for item in _items(annotation.value):
            identity = _term_identity(item)
            if identity is None:
                continue
            key, label = identity
            entries[(annotation.annotation_type, key)] = (
                label,
                annotation,
                _annotation_provenance(annotation, entity_context, context),
            )
    return entries


def _functional_entries(
    context: LocalResearchContext,
    entity_context: LocalEntityContext,
) -> Mapping[tuple[Any, ...], tuple[FunctionalContextRow, tuple[ProvenanceRecord, ...]]]:
    entries: dict[tuple[Any, ...], tuple[FunctionalContextRow, tuple[ProvenanceRecord, ...]]] = {}
    identifiers = _entity_identifiers(entity_context)
    for annotation in entity_context.annotations:
        if annotation.annotation_type not in {"go_enrichment", "kegg_enrichment"}:
            continue
        for value in _items(annotation.value):
            if not isinstance(value, FunctionalContextRow):
                continue
            key = (
                value.kind.value,
                value.source_library or "",
                value.pathway_family or "",
                value.term or "",
            )
            entries[key] = (
                value,
                _annotation_provenance(annotation, entity_context, context),
            )

    # A valid LocalResearchContext normally attaches these rows to members.
    # This exact-membership fallback handles contexts built by older adapters;
    # it does not calculate enrichment or expand pathway membership.
    for value in context.functional_context:
        member_keys = {str(member).strip().casefold() for member in value.members if str(member).strip()}
        if not (identifiers & member_keys):
            continue
        key = (
            value.kind.value,
            value.source_library or "",
            value.pathway_family or "",
            value.term or "",
        )
        entries.setdefault(
            key,
            (
                value,
                _combine_provenance(
                    entity_context.provenance,
                    fallback=context.provenance,
                ),
            ),
        )
    return entries


def _functional_assertions(
    context: LocalResearchContext,
    target: LocalEntityContext,
    response: LocalEntityContext,
) -> tuple[EvidenceAssertion, ...]:
    assertions: list[EvidenceAssertion] = []
    target_go = _go_annotation_entries(context, target)
    response_go = _go_annotation_entries(context, response)
    for key in sorted(set(target_go) & set(response_go)):
        left, right = target_go[key], response_go[key]
        provenance = _combine_provenance(left[2], right[2])
        basis = {"annotation_type": key[0], "annotation_key": key[1], "label": left[0]}
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.SHARED_GO,
            subject=target.entity,
            object_entity=response.entity,
            statement="Target and returned response share an exact existing GO annotation.",
            provenance=provenance,
            value=immutable_mapping({"annotation_type": key[0], "term": left[0]}),
            qualifiers=immutable_mapping({"basis": "exact_shared_go_annotation"}),
            basis=basis,
        )
        if assertion is not None:
            assertions.append(assertion)

    target_pathways = _pathway_annotation_entries(context, target)
    response_pathways = _pathway_annotation_entries(context, response)
    for key in sorted(set(target_pathways) & set(response_pathways)):
        left, right = target_pathways[key], response_pathways[key]
        provenance = _combine_provenance(left[2], right[2])
        basis = {"pathway_key": key, "label": left[0]}
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.SHARED_PATHWAY,
            subject=target.entity,
            object_entity=response.entity,
            statement=(
                "Target and returned response share an exact existing Reactome "
                "pathway annotation."
            ),
            provenance=provenance,
            value=immutable_mapping({"pathway": left[0]}),
            qualifiers=immutable_mapping({"basis": "exact_shared_reactome_annotation"}),
            basis=basis,
        )
        if assertion is not None:
            assertions.append(assertion)

    target_functional = _functional_entries(context, target)
    response_functional = _functional_entries(context, response)
    for key in sorted(set(target_functional) & set(response_functional)):
        left, right = target_functional[key], response_functional[key]
        row = left[0]
        relationship_type = {
            FunctionalContextKind.GO: RelationshipType.SHARED_GO,
            FunctionalContextKind.KEGG: RelationshipType.SHARED_PATHWAY,
        }.get(row.kind)
        if relationship_type is None:
            continue
        provenance = _combine_provenance(left[1], right[1])
        basis = {
            "kind": row.kind.value,
            "source_library": row.source_library,
            "pathway_family": row.pathway_family,
            "term": row.term,
        }
        assertion = _make_assertion(
            context=context,
            relationship_type=relationship_type,
            subject=target.entity,
            object_entity=response.entity,
            statement=(
                "Target and returned response occur in the same existing GO result row."
                if relationship_type is RelationshipType.SHARED_GO
                else "Target and returned response occur in the same existing pathway result row."
            ),
            provenance=provenance,
            value=immutable_mapping(
                {
                    "term": row.term,
                    "source_library": row.source_library,
                    "pathway_family": row.pathway_family,
                }
            ),
            qualifiers=immutable_mapping({"basis": "existing_functional_context_membership"}),
            basis=basis,
        )
        if assertion is not None:
            assertions.append(assertion)
    return tuple(assertions)


def _pathway_annotation_entries(
    context: LocalResearchContext,
    entity_context: LocalEntityContext,
) -> Mapping[str, tuple[str, LocalAnnotation, tuple[ProvenanceRecord, ...]]]:
    """Return exact, already-attached Reactome pathway annotations."""

    entries: dict[str, tuple[str, LocalAnnotation, tuple[ProvenanceRecord, ...]]] = {}
    for annotation in entity_context.annotations:
        if annotation.annotation_type not in _PATHWAY_ANNOTATION_TYPES:
            continue
        for item in _items(annotation.value):
            identity = _term_identity(item)
            if identity is None:
                continue
            key, label = identity
            entries[key] = (
                label,
                annotation,
                _annotation_provenance(annotation, entity_context, context),
            )
    return entries


def _localization_entries(
    context: LocalResearchContext,
    entity_context: LocalEntityContext,
) -> Mapping[str, tuple[str, tuple[ProvenanceRecord, ...]]]:
    entries: dict[str, tuple[str, tuple[ProvenanceRecord, ...]]] = {}
    for annotation in _annotations(entity_context, "localization"):
        for value in _items(annotation.value):
            text = str(value or "").strip()
            if text and text.casefold() not in {"nan", "none", "unknown", "—"}:
                entries[text.casefold()] = (
                    text,
                    _annotation_provenance(annotation, entity_context, context),
                )
    return entries


def _compartment_assertions(
    context: LocalResearchContext,
    target: LocalEntityContext,
    response: LocalEntityContext,
) -> tuple[EvidenceAssertion, ...]:
    target_values = _localization_entries(context, target)
    response_values = _localization_entries(context, response)
    assertions: list[EvidenceAssertion] = []
    for key in sorted(set(target_values) & set(response_values)):
        provenance = _combine_provenance(target_values[key][1], response_values[key][1])
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.SHARED_COMPARTMENT,
            subject=target.entity,
            object_entity=response.entity,
            statement="Target and returned response share an exact existing localization annotation.",
            provenance=provenance,
            value=immutable_mapping({"localization": target_values[key][0]}),
            qualifiers=immutable_mapping({"basis": "exact_shared_localization"}),
            basis={"localization": key},
        )
        if assertion is not None:
            assertions.append(assertion)
    return tuple(assertions)


def _community_entries(entity_context: LocalEntityContext) -> Mapping[str, Any]:
    entries: dict[str, Any] = {}
    for record in entity_context.source_records:
        for field_name in _COMMUNITY_FIELDS:
            if field_name not in record.values:
                continue
            value = record.values[field_name]
            text = str(value or "").strip()
            if text.casefold() in {"", "nan", "none", "unknown", "n/a", "-1", "—"}:
                continue
            entries[text.casefold()] = value
    return entries


def _community_assertions(
    context: LocalResearchContext,
    target: LocalEntityContext,
    response: LocalEntityContext,
) -> tuple[EvidenceAssertion, ...]:
    target_values = _community_entries(target)
    response_values = _community_entries(response)
    assertions: list[EvidenceAssertion] = []
    for key in sorted(set(target_values) & set(response_values)):
        provenance = _combine_provenance(
            target.provenance,
            response.provenance,
            fallback=context.provenance,
        )
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.SHARED_COMMUNITY,
            subject=target.entity,
            object_entity=response.entity,
            statement="Target and returned response have the same existing Sophiark community identifier.",
            provenance=provenance,
            value=immutable_mapping({"community": target_values[key]}),
            qualifiers=immutable_mapping({"basis": "existing_community_identifier"}),
            basis={"community": key},
        )
        if assertion is not None:
            assertions.append(assertion)
    return tuple(assertions)


def _regulatory_records(
    entity_context: LocalEntityContext,
) -> tuple[tuple[RegulatoryRelation, LocalAnnotation], ...]:
    records: list[tuple[RegulatoryRelation, LocalAnnotation]] = []
    for annotation in _annotations(entity_context, "directed_regulation"):
        for value in _items(annotation.value):
            if isinstance(value, RegulatoryRelation):
                records.append((value, annotation))
    return tuple(records)


def _regulatory_assertions(
    context: LocalResearchContext,
    target: LocalEntityContext,
    response: LocalEntityContext,
) -> tuple[EvidenceAssertion, ...]:
    target_ids = _entity_identifiers(target)
    response_ids = _entity_identifiers(response)
    seen: dict[tuple[Any, ...], tuple[RegulatoryRelation, LocalAnnotation, LocalEntityContext]] = {}
    for owner in (target, response):
        for relation, annotation in _regulatory_records(owner):
            key = (
                relation.taxon_id,
                relation.database,
                relation.source_symbol.casefold(),
                relation.target_symbol.casefold(),
                relation.effect,
                relation.mechanism,
                relation.reference_ids,
                relation.record_identifier,
            )
            seen[key] = (relation, annotation, owner)

    assertions: list[EvidenceAssertion] = []
    for key in sorted(seen, key=repr):
        relation, annotation, owner = seen[key]
        if relation.taxon_id != context.semantic_result.scope.taxon_id:
            continue
        source_key = relation.source_symbol.casefold()
        target_key = relation.target_symbol.casefold()
        if source_key in target_ids and target_key in response_ids:
            subject, object_entity = target.entity, response.entity
        elif source_key in response_ids and target_key in target_ids:
            subject, object_entity = response.entity, target.entity
        else:
            continue
        provenance = _combine_provenance(
            (relation.provenance,),
            _annotation_provenance(annotation, owner, context),
        )
        basis = {
            "database": relation.database,
            "source_symbol": relation.source_symbol,
            "target_symbol": relation.target_symbol,
            "record_identifier": relation.record_identifier,
            "reference_ids": relation.reference_ids,
        }
        assertion = _make_assertion(
            context=context,
            relationship_type=RelationshipType.DIRECTED_REGULATION,
            subject=subject,
            object_entity=object_entity,
            statement="An existing directed-regulatory database record connects the target and returned response.",
            provenance=provenance,
            value=immutable_mapping(
                {
                    "database": relation.database,
                    "effect": relation.effect,
                    "mechanism": relation.mechanism,
                    "reference_ids": relation.reference_ids,
                    "record_identifier": relation.record_identifier,
                }
            ),
            qualifiers=immutable_mapping({"basis": "existing_directed_database_record"}),
            basis=basis,
        )
        if assertion is not None:
            assertions.append(assertion)
    return tuple(assertions)


def _normalise_explicit(
    context: LocalResearchContext,
    assertion: EvidenceAssertion,
    *,
    subject: Entity,
    object_entity: Entity,
) -> EvidenceAssertion:
    assert assertion.object is not None
    basis = {
        "statement": assertion.statement,
        "value": assertion.value,
        "qualifiers": assertion.qualifiers,
        "provenance": assertion.provenance,
    }
    return EvidenceAssertion(
        assertion_id=_deterministic_id(
            simulation_id=context.simulation_id,
            relationship_type=assertion.relationship_type,
            subject=subject,
            object_entity=object_entity,
            basis=basis,
        ),
        subject=subject,
        relationship_type=assertion.relationship_type,
        object=object_entity,
        statement=assertion.statement,
        provenance=assertion.provenance,
        value=assertion.value,
        qualifiers=assertion.qualifiers,
    )


def build_target_response_relationships(
    context: LocalResearchContext,
    *,
    explicit_evidence: Iterable[EvidenceAssertion] = (),
) -> RelationshipBuildResult:
    """Build evidence-only relationships for exact target/response pairs.

    PPI edges, network paths and publication co-mentions are admitted only from
    ``explicit_evidence``.  No relationship type is converted to another type.
    """

    if not isinstance(context, LocalResearchContext):
        raise TypeError("context must be a LocalResearchContext")
    taxon_id = context.semantic_result.scope.taxon_id
    cross_taxon_context = tuple(
        entity_context.entity.logical_key
        for entity_context in context.entities
        if entity_context.entity.taxon_id != taxon_id
    )
    if cross_taxon_context:
        raise ValueError("LocalResearchContext contains entities outside its scientific taxon scope")

    targets = tuple(
        entity_context
        for entity_context in context.entities
        if ResultRole.PERTURBATION_TARGET in entity_context.roles
    )
    target_keys = {item.entity.logical_key for item in targets}
    responses = tuple(
        entity_context
        for entity_context in context.entities
        if item_has_response_role(entity_context)
        and entity_context.entity.logical_key not in target_keys
    )
    response_keys = {item.entity.logical_key for item in responses}
    scoped_entities = {
        item.entity.logical_key: item.entity
        for item in (*targets, *responses)
    }

    notices: list[str] = []
    assertions: list[EvidenceAssertion] = []
    if not targets:
        notices.append("No resolved perturbation target is available for relationship assembly.")
    if not responses:
        notices.append("No resolved returned-response entity is available for relationship assembly.")

    for target in targets:
        for response in responses:
            generated = (
                *_regulatory_assertions(context, target, response),
                *_functional_assertions(context, target, response),
                *_family_assertions(context, target, response),
                *_compartment_assertions(context, target, response),
                *_community_assertions(context, target, response),
            )
            assertions.extend(generated)

    rejected_explicit = 0
    for candidate in tuple(explicit_evidence):
        if not isinstance(candidate, EvidenceAssertion) or candidate.object is None:
            rejected_explicit += 1
            continue
        if (
            candidate.subject.taxon_id != taxon_id
            or candidate.object.taxon_id != taxon_id
        ):
            rejected_explicit += 1
            continue
        subject_key = candidate.subject.logical_key
        object_key = candidate.object.logical_key
        is_target_response = (
            (subject_key in target_keys and object_key in response_keys)
            or (object_key in target_keys and subject_key in response_keys)
        )
        if not is_target_response:
            rejected_explicit += 1
            continue
        assertions.append(
            _normalise_explicit(
                context,
                candidate,
                subject=scoped_entities[subject_key],
                object_entity=scoped_entities[object_key],
            )
        )

    if rejected_explicit:
        notices.append(
            f"Rejected {rejected_explicit} explicit relationship item(s) outside the exact target/response taxon scope."
        )

    unique = {assertion.assertion_id: assertion for assertion in assertions}
    target_position = {item.entity.logical_key: index for index, item in enumerate(targets)}
    response_position = {item.entity.logical_key: index for index, item in enumerate(responses)}

    def sort_key(assertion: EvidenceAssertion) -> tuple[int, int, int, str]:
        object_entity = assertion.object
        assert object_entity is not None
        subject_key = assertion.subject.logical_key
        object_key = object_entity.logical_key
        if subject_key in target_position:
            target_key, response_key = subject_key, object_key
        else:
            target_key, response_key = object_key, subject_key
        return (
            target_position.get(target_key, len(target_position)),
            response_position.get(response_key, len(response_position)),
            _RELATIONSHIP_ORDER[assertion.relationship_type],
            assertion.assertion_id,
        )

    ordered = tuple(sorted(unique.values(), key=sort_key))
    if not ordered and targets and responses:
        notices.append("No exact existing relationship evidence matched the target/response pairs.")
    return RelationshipBuildResult(
        simulation_id=context.simulation_id,
        snapshot_id=context.snapshot_id,
        assertions=ordered,
        notices=tuple(dict.fromkeys(notices)),
        rejected_explicit_count=rejected_explicit,
    )


def item_has_response_role(entity_context: LocalEntityContext) -> bool:
    """Return whether an entity is in an existing returned-response role."""

    return bool(set(entity_context.roles) & _RESPONSE_ROLES)


__all__ = [
    "RelationshipBuildResult",
    "build_target_response_relationships",
    "item_has_response_role",
]
