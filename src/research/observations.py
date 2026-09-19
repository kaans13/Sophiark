"""Deterministic, local-only Observation Detector for Research Context.

The detector consumes an immutable :class:`LocalResearchContext` and optional
explicit relationship assertions.  It never reads a graph, session state,
files, caches, providers, or the network, and it never computes a scientific
rank, score, p-value, enrichment, or significance claim.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass, field
from datetime import datetime, timezone
import hashlib
import json
import math
from typing import Any

from .config import DEFAULT_RESEARCH_CONFIG, ResearchConfig
from .local_context import LocalAnnotation, LocalEntityContext, LocalResearchContext
from .models import (
    Entity,
    EvidenceAssertion,
    Observation,
    ObservationType,
    RelationshipType,
    immutable_mapping,
)
from .providers import FamilyMembership
from .semantic_adapter import FunctionalContextKind, FunctionalContextRow, ResultRole


_RESPONSE_ROLES = frozenset({
    ResultRole.POSITIVE_REDISTRIBUTION,
    ResultRole.NETWORK_LOSS,
    ResultRole.FDR_SUPPORTED,
})
_RELATIONSHIP_ORDER = {
    relationship_type: position
    for position, relationship_type in enumerate(RelationshipType)
}
_FUNCTIONAL_TYPES = {
    "go_biological_process": ("GO", "GO biological process"),
    "go_cellular_component": ("GO", "GO cellular component"),
    "go_molecular_function": ("GO", "GO molecular function"),
    "go_annotation": ("GO", "GO annotation"),
    "kegg": ("KEGG", "KEGG annotation"),
    "kegg_pathway": ("KEGG", "KEGG pathway"),
    "kegg_annotation": ("KEGG", "KEGG annotation"),
    "reactome": ("Reactome", "Reactome annotation"),
    "reactome_pathway": ("Reactome", "Reactome pathway"),
    "reactome_annotation": ("Reactome", "Reactome annotation"),
}


@dataclass(frozen=True, slots=True)
class DetectedObservation:
    """A base observation plus typed, immutable pattern details."""

    observation: Observation
    details: Mapping[str, Any] = field(default_factory=immutable_mapping)
    assertions: tuple[EvidenceAssertion, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "details", immutable_mapping(self.details))
        object.__setattr__(self, "assertions", tuple(self.assertions))

    @property
    def observation_id(self) -> str:
        return self.observation.observation_id

    @property
    def type(self) -> ObservationType:
        return self.observation.type

    @property
    def members(self) -> tuple[Entity, ...]:
        return self.observation.members

    @property
    def basis(self) -> str:
        return self.observation.basis

    @property
    def source_fields(self) -> tuple[str, ...]:
        return self.observation.source_fields

    @property
    def created_at(self) -> str:
        return self.observation.created_at

    @property
    def limitations(self) -> str:
        return self.observation.limitations


@dataclass(frozen=True, slots=True)
class ObservationDetectionResult:
    simulation_id: str
    snapshot_id: str
    observations: tuple[DetectedObservation, ...]
    notices: tuple[str, ...] = ()
    rejected_relationship_count: int = 0

    def __post_init__(self) -> None:
        object.__setattr__(self, "observations", tuple(self.observations))
        object.__setattr__(self, "notices", tuple(map(str, self.notices)))
        if self.rejected_relationship_count < 0:
            raise ValueError("rejected_relationship_count cannot be negative")

    def of_type(self, observation_type: ObservationType | str) -> tuple[DetectedObservation, ...]:
        resolved = ObservationType(observation_type)
        return tuple(item for item in self.observations if item.type is resolved)

    @property
    def base_observations(self) -> tuple[Observation, ...]:
        return tuple(item.observation for item in self.observations)


def _unique(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value) for value in values if str(value)))


def _missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, str):
        return value.strip().casefold() in {"", "nan", "none", "<na>", "—"}
    if isinstance(value, (tuple, list, set, frozenset, Mapping)):
        return len(value) == 0
    return False


def _strict_true(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value == 1
    return isinstance(value, str) and value.strip().casefold() in {"true", "1", "yes", "evet", "✓"}


def _normalized_created_at(
    context: LocalResearchContext, value: datetime | str | None,
) -> str:
    candidate: datetime | str | None = value
    if candidate is None:
        candidate = next(
            (record.retrieved_at for record in context.provenance if record.retrieved_at),
            None,
        )
    if candidate is None:
        raise ValueError(
            "created_at must be explicit when local context provenance has no snapshot timestamp"
        )
    if isinstance(candidate, datetime):
        timestamp = candidate
    else:
        text = str(candidate).strip()
        if text.endswith("Z"):
            text = f"{text[:-1]}+00:00"
        timestamp = datetime.fromisoformat(text)
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    return timestamp.astimezone(timezone.utc).isoformat()


def _observation_id(
    *, simulation_id: str, observation_type: ObservationType,
    discriminator: str, members: Iterable[Entity], details: Mapping[str, Any],
    assertion_ids: Iterable[str] = (),
) -> str:
    def stable(value: Any) -> Any:
        if isinstance(value, Mapping):
            return {
                str(key): stable(item)
                for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
            }
        if isinstance(value, (tuple, list)):
            return [stable(item) for item in value]
        if isinstance(value, (set, frozenset)):
            return sorted((stable(item) for item in value), key=repr)
        if isinstance(value, (ObservationType, RelationshipType, ResultRole)):
            return value.value
        if isinstance(value, Entity):
            return value.logical_key
        if isinstance(value, (str, int, float, bool)) or value is None:
            return value
        return str(value)

    payload = {
        "simulation_id": simulation_id,
        "type": observation_type.value,
        "discriminator": discriminator,
        "members": sorted(entity.logical_key for entity in members),
        "assertions": sorted(map(str, assertion_ids)),
        "details": stable(details),
    }
    encoded = json.dumps(
        payload, ensure_ascii=False, allow_nan=False,
        sort_keys=True, separators=(",", ":"),
    ).encode("utf-8")
    digest = hashlib.sha256(encoded).hexdigest()[:24]
    return f"observation:{observation_type.value.casefold()}:{digest}"


def _detected(
    context: LocalResearchContext,
    *,
    observation_type: ObservationType,
    discriminator: str,
    members: Iterable[Entity],
    basis: str,
    source_fields: Iterable[str],
    created_at: str,
    limitations: str,
    details: Mapping[str, Any],
    assertions: Iterable[EvidenceAssertion] = (),
) -> DetectedObservation:
    member_tuple = tuple(members)
    assertion_tuple = tuple(assertions)
    observation = Observation(
        observation_id=_observation_id(
            simulation_id=context.simulation_id,
            observation_type=observation_type,
            discriminator=discriminator,
            members=member_tuple,
            details=details,
            assertion_ids=(item.assertion_id for item in assertion_tuple),
        ),
        simulation_id=context.simulation_id,
        type=observation_type,
        members=member_tuple,
        basis=basis,
        scope=context.semantic_result.scope,
        source_fields=_unique(source_fields),
        created_at=created_at,
        limitations=limitations,
    )
    return DetectedObservation(observation, details, assertion_tuple)


def _responses(context: LocalResearchContext) -> tuple[LocalEntityContext, ...]:
    return tuple(
        item for item in context.entities
        if any(role in _RESPONSE_ROLES for role in item.roles)
    )


def _annotation_items(value: Any) -> tuple[Any, ...]:
    if _missing(value):
        return ()
    if isinstance(value, (tuple, list, set, frozenset)):
        return tuple(value)
    return (value,)


def _annotation_identity(value: Any) -> tuple[str, str] | None:
    """Return a stable identity and display label without biological inference."""

    if _missing(value):
        return None
    if isinstance(value, Mapping):
        identifier = next(
            (str(value[key]).strip() for key in ("id", "accession", "term_id", "pathway_id")
             if key in value and not _missing(value[key])),
            None,
        )
        label = next(
            (str(value[key]).strip() for key in ("term", "name", "label", "description")
             if key in value and not _missing(value[key])),
            identifier,
        )
        if identifier or label:
            return identifier or label or "", label or identifier or ""
        return None
    text = str(value).strip()
    return (text, text) if text and not _missing(text) else None


def _family_observations(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    config: ResearchConfig,
    created_at: str,
) -> tuple[DetectedObservation, ...]:
    groups: dict[tuple[int, str], dict[str, Any]] = {}
    for entity_context in responses:
        for annotation in entity_context.annotations:
            membership = annotation.value
            if (
                annotation.annotation_type != "authoritative_family_membership"
                or not isinstance(membership, FamilyMembership)
                or membership.evidence_type != "AUTHORITATIVE_FAMILY_MEMBERSHIP"
                or membership.taxon_id != context.semantic_result.scope.taxon_id
            ):
                continue
            group_id = str(membership.group_id).strip()
            if not group_id:
                continue
            key = (membership.taxon_id, group_id)
            group = groups.setdefault(key, {
                "membership": membership,
                "member_keys": set(),
                "source_fields": [],
            })
            group["member_keys"].add(entity_context.entity.logical_key)
            group["source_fields"].extend((
                annotation.annotation_type,
                annotation.source_field or "group_id/group_name",
            ))

    minimum = config.presentation_limits.family_min_members
    observations: list[DetectedObservation] = []
    for (_, group_id), group in groups.items():
        member_keys = group["member_keys"]
        members = tuple(
            item.entity for item in responses
            if item.entity.logical_key in member_keys
        )
        if len(members) < minimum:
            continue
        membership: FamilyMembership = group["membership"]
        observations.append(_detected(
            context,
            observation_type=ObservationType.FAMILY_COOCCURRENCE,
            discriminator=f"family:{membership.taxon_id}:{group_id}",
            members=members,
            basis=(
                f"{len(members)} returned response entities share the same "
                f"authoritative family annotation ({membership.group_name})."
            ),
            source_fields=group["source_fields"],
            created_at=created_at,
            limitations=(
                "This presentation-only co-occurrence does not demonstrate shared activation, "
                "a shared mechanism, causality, enrichment, significance, or experimental validation."
            ),
            details={
                "family_id": group_id,
                "family_name": membership.group_name,
                "member_count": len(members),
                "presentation_min_members": minimum,
                "presentation_only": config.presentation_limits.family_threshold_presentation_only,
            },
        ))
    return tuple(observations)


def _functional_descriptors(
    annotation: LocalAnnotation,
) -> tuple[tuple[tuple[str, str, str], Mapping[str, Any]], ...]:
    value = annotation.value
    if annotation.annotation_type in {"go_enrichment", "kegg_enrichment"}:
        if not isinstance(value, FunctionalContextRow):
            return ()
        kind = {
            FunctionalContextKind.GO: "GO",
            FunctionalContextKind.KEGG: "KEGG",
        }.get(value.kind)
        if kind is None or not value.term:
            return ()
        identity = value.term
        return (((kind, annotation.annotation_type, identity), {
            "context_kind": kind,
            "term": value.term,
            "annotation_type": annotation.annotation_type,
            "source_library": value.source_library,
            "source_fields": (
                annotation.annotation_type,
                annotation.source_field or "Genes",
                "Term", "Pathway_Family", "Kaynak",
            ),
        }),)

    descriptor = _FUNCTIONAL_TYPES.get(annotation.annotation_type)
    if descriptor is None:
        return ()
    kind, descriptor_name = descriptor
    results: list[tuple[tuple[str, str, str], Mapping[str, Any]]] = []
    for item in _annotation_items(value):
        identity = _annotation_identity(item)
        if identity is None:
            continue
        stable_identity, label = identity
        results.append(((kind, annotation.annotation_type, stable_identity), {
            "context_kind": kind,
            "term": label,
            "annotation_type": descriptor_name,
            "source_library": annotation.source,
            "source_fields": (
                annotation.annotation_type,
                annotation.source_field or annotation.annotation_type,
            ),
        }))
    return tuple(results)


def _functional_observations(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    created_at: str,
) -> tuple[DetectedObservation, ...]:
    groups: dict[tuple[str, str, str], dict[str, Any]] = {}
    for entity_context in responses:
        seen_for_entity: set[tuple[str, str, str]] = set()
        for annotation in entity_context.annotations:
            for key, descriptor in _functional_descriptors(annotation):
                if key in seen_for_entity:
                    continue
                seen_for_entity.add(key)
                group = groups.setdefault(key, {
                    "descriptor": descriptor,
                    "member_keys": set(),
                    "source_fields": [],
                })
                group["member_keys"].add(entity_context.entity.logical_key)
                group["source_fields"].extend(descriptor["source_fields"])

    observations: list[DetectedObservation] = []
    for key, group in groups.items():
        members = tuple(
            item.entity for item in responses
            if item.entity.logical_key in group["member_keys"]
        )
        if len(members) < 2:
            continue
        descriptor = group["descriptor"]
        observations.append(_detected(
            context,
            observation_type=ObservationType.FUNCTIONAL_COOCCURRENCE,
            discriminator="functional:" + ":".join(key),
            members=members,
            basis=(
                f"{len(members)} returned response entities share the existing "
                f"{descriptor['context_kind']} annotation '{descriptor['term']}'."
            ),
            source_fields=group["source_fields"],
            created_at=created_at,
            limitations=(
                "This is shared annotation observed within the returned set. It does not "
                "establish activation, suppression, causality, a dominant mechanism, or new enrichment."
            ),
            details={
                "context_kind": descriptor["context_kind"],
                "term": descriptor["term"],
                "annotation_type": descriptor["annotation_type"],
                "source_library": descriptor["source_library"],
                "member_count": len(members),
            },
        ))
    return tuple(observations)


def _localization_labels(annotation: LocalAnnotation) -> tuple[str, ...]:
    if annotation.annotation_type not in {"localization", "go_cellular_component"}:
        return ()
    labels: list[str] = []
    for item in _annotation_items(annotation.value):
        identity = _annotation_identity(item)
        if identity is not None:
            labels.append(identity[1])
    return _unique(labels)


def _localization_observation(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    created_at: str,
) -> DetectedObservation | None:
    distribution: dict[str, int] = {}
    annotated_keys: set[str] = set()
    source_fields: list[str] = []
    for entity_context in responses:
        labels_for_entity: set[str] = set()
        for annotation in entity_context.annotations:
            labels = _localization_labels(annotation)
            if not labels:
                continue
            labels_for_entity.update(labels)
            source_fields.extend((
                annotation.annotation_type,
                annotation.source_field or annotation.annotation_type,
            ))
        if labels_for_entity:
            annotated_keys.add(entity_context.entity.logical_key)
        for label in sorted(labels_for_entity):
            distribution[label] = distribution.get(label, 0) + 1
    if not distribution:
        return None
    members = tuple(
        item.entity for item in responses
        if item.entity.logical_key in annotated_keys
    )
    annotated_count = len(members)
    fractions = {
        label: count / annotated_count
        for label, count in distribution.items()
    }
    return _detected(
        context,
        observation_type=ObservationType.LOCALIZATION_PATTERN,
        discriminator="localization-distribution",
        members=members,
        basis=(
            f"Existing localization annotations are available for {annotated_count} of "
            f"{len(responses)} returned response entities."
        ),
        source_fields=source_fields,
        created_at=created_at,
        limitations=(
            "This descriptive distribution does not show that perturbation occurred in a "
            "compartment and does not establish enrichment, significance, activity, or causality."
        ),
        details={
            "distribution": distribution,
            "fractions_of_annotated_members": fractions,
            "total_count": annotated_count,
            "returned_response_count": len(responses),
        },
    )


def _contrast_descriptors(
    annotation: LocalAnnotation,
) -> tuple[tuple[tuple[str, str], Mapping[str, str]], ...]:
    if (
        annotation.annotation_type == "authoritative_family_membership"
        and isinstance(annotation.value, FamilyMembership)
        and annotation.value.evidence_type == "AUTHORITATIVE_FAMILY_MEMBERSHIP"
    ):
        membership = annotation.value
        label = f"{membership.group_name} ({membership.group_id})"
        return ((("authoritative_family_membership", membership.group_id), {
            "annotation_type": "authoritative_family_membership",
            "annotation": label,
            "source_field": annotation.source_field or "group_id/group_name",
        }),)

    functional = _functional_descriptors(annotation)
    if functional:
        return tuple(
            ((f"functional:{key[0]}:{key[1]}", key[2]), {
                "annotation_type": descriptor["annotation_type"],
                "annotation": str(descriptor["term"]),
                "source_field": annotation.source_field or annotation.annotation_type,
            })
            for key, descriptor in functional
        )

    labels = _localization_labels(annotation)
    return tuple(
        (("localization", label), {
            "annotation_type": "localization",
            "annotation": label,
            "source_field": annotation.source_field or annotation.annotation_type,
        })
        for label in labels
    )


def _positive_negative_observation(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    created_at: str,
) -> DetectedObservation | None:
    positive = tuple(item for item in responses if ResultRole.POSITIVE_REDISTRIBUTION in item.roles)
    negative = tuple(item for item in responses if ResultRole.NETWORK_LOSS in item.roles)
    if not positive or not negative:
        return None

    rows: dict[tuple[str, str], dict[str, Any]] = {}
    for role_name, entity_contexts in (("positive", positive), ("negative", negative)):
        for entity_context in entity_contexts:
            seen_for_entity: set[tuple[str, str]] = set()
            for annotation in entity_context.annotations:
                for key, descriptor in _contrast_descriptors(annotation):
                    if key in seen_for_entity:
                        continue
                    seen_for_entity.add(key)
                    row = rows.setdefault(key, {
                        **descriptor,
                        "positive_keys": set(),
                        "negative_keys": set(),
                    })
                    row[f"{role_name}_keys"].add(entity_context.entity.logical_key)
    if not rows:
        return None

    detail_rows: list[Mapping[str, Any]] = []
    source_fields: list[str] = ["semantic.roles"]
    for row in rows.values():
        positive_count = len(row["positive_keys"])
        negative_count = len(row["negative_keys"])
        source_fields.append(row["source_field"])
        detail_rows.append({
            "annotation_type": row["annotation_type"],
            "annotation": row["annotation"],
            "positive_count": positive_count,
            "negative_count": negative_count,
            "positive_fraction": positive_count / len(positive),
            "negative_fraction": negative_count / len(negative),
        })
    member_keys = {
        item.entity.logical_key for item in (*positive, *negative)
    }
    members = tuple(
        item.entity for item in responses if item.entity.logical_key in member_keys
    )
    return _detected(
        context,
        observation_type=ObservationType.POSITIVE_NEGATIVE_CONTRAST,
        discriminator="positive-negative-annotation-distribution",
        members=members,
        basis=(
            f"Existing annotations were compared descriptively across {len(positive)} positive "
            f"and {len(negative)} network-loss response entities."
        ),
        source_fields=source_fields,
        created_at=created_at,
        limitations=(
            "Counts and fractions are descriptive only. No p-value, enrichment, significance, "
            "ranking, activation, suppression, or causal conclusion was computed."
        ),
        details={
            "annotations": tuple(detail_rows),
            "positive_total": len(positive),
            "negative_total": len(negative),
        },
    )


def _gateway_observation(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    created_at: str,
) -> DetectedObservation | None:
    gateway_keys: set[str] = set()
    source_fields: list[str] = []
    for entity_context in responses:
        gateway_annotations = tuple(
            annotation for annotation in entity_context.annotations
            if annotation.annotation_type == "gateway" and _strict_true(annotation.value)
        )
        if not gateway_annotations:
            continue
        gateway_keys.add(entity_context.entity.logical_key)
        source_fields.extend(
            annotation.source_field or "Gümrük_Kapisi"
            for annotation in gateway_annotations
        )
    if not gateway_keys:
        return None

    role_map: dict[str, tuple[str, ...]] = {}
    role_specs = (
        ResultRole.POSITIVE_REDISTRIBUTION,
        ResultRole.NETWORK_LOSS,
        ResultRole.FDR_SUPPORTED,
    )
    for role in role_specs:
        if role is ResultRole.FDR_SUPPORTED and context.semantic_result.scope.selection_mode.casefold() != "null_fdr":
            continue
        keys = tuple(
            item.entity.logical_key for item in responses
            if item.entity.logical_key in gateway_keys and role in item.roles
        )
        if keys:
            role_map[role.value] = keys
    if not role_map:
        return None
    members = tuple(
        item.entity for item in responses
        if item.entity.logical_key in gateway_keys
        and any(item.entity.logical_key in keys for keys in role_map.values())
    )
    return _detected(
        context,
        observation_type=ObservationType.GATEWAY_RESPONSE_OVERLAP,
        discriminator="gateway-response-role-overlap",
        members=members,
        basis=(
            f"{len(members)} returned response entities match the existing Compartment Bottleneck rule and "
            "one or more preserved response roles."
        ),
        source_fields=(*source_fields, "semantic.roles"),
        created_at=created_at,
        limitations=(
            "Compartment Bottleneck is an existing topology rule. This overlap does not establish biological "
            "necessity, therapeutic suitability, significance, causality, or experimental validation."
        ),
        details={
            "overlap_roles": role_map,
            "member_count": len(members),
        },
    )


def _fdr_observation(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    *,
    created_at: str,
) -> DetectedObservation | None:
    selection_mode = context.semantic_result.scope.selection_mode
    if selection_mode.casefold() != "null_fdr":
        return None
    members = tuple(
        item.entity for item in responses
        if ResultRole.FDR_SUPPORTED in item.roles
    )
    if not members:
        return None
    return _detected(
        context,
        observation_type=ObservationType.FDR_SUPPORTED_SUBSET,
        discriminator="stored-fdr-supported-subset",
        members=members,
        basis=(
            f"{len(members)} returned response entities retain the stored FDR-supported role "
            "from the current null_fdr result."
        ),
        source_fields=("significant_redistribution", "semantic.roles"),
        created_at=created_at,
        limitations=(
            "The detector used only the stored semantic FDR role. It did not recompute p-values, "
            "q-values, multiple-testing correction, significance, ranking, or causal interpretation."
        ),
        details={
            "member_count": len(members),
            "selection_mode": selection_mode,
            "stored_role": ResultRole.FDR_SUPPORTED.value,
        },
    )


def _relationship_input(
    value: Iterable[EvidenceAssertion] | Any,
) -> tuple[tuple[EvidenceAssertion, ...], tuple[str, ...], int]:
    wrapper_notices: tuple[str, ...] = ()
    prior_rejected = 0
    if value is None:
        candidates: Iterable[Any] = ()
    elif hasattr(value, "assertions") and not isinstance(value, EvidenceAssertion):
        candidates = getattr(value, "assertions")
        wrapper_notices = tuple(map(str, getattr(value, "notices", ())))
        prior_rejected = int(getattr(value, "rejected_explicit_count", 0))
    else:
        candidates = value
    assertions = tuple(candidates)
    if any(not isinstance(item, EvidenceAssertion) for item in assertions):
        raise TypeError("relationship_assertions must contain only EvidenceAssertion objects")
    return assertions, wrapper_notices, prior_rejected


def _relationship_observations(
    context: LocalResearchContext,
    responses: tuple[LocalEntityContext, ...],
    assertions: tuple[EvidenceAssertion, ...],
    *,
    created_at: str,
) -> tuple[tuple[DetectedObservation, ...], int]:
    targets = tuple(
        item for item in context.entities
        if ResultRole.PERTURBATION_TARGET in item.roles
    )
    response_by_key = {item.entity.logical_key: item for item in responses}
    target_by_key = {item.entity.logical_key: item for item in targets}
    pair_assertions: dict[tuple[str, str], list[EvidenceAssertion]] = {}
    rejected = 0
    seen_ids: set[str] = set()
    for assertion in assertions:
        if assertion.assertion_id in seen_ids:
            continue
        seen_ids.add(assertion.assertion_id)
        if assertion.object is None:
            rejected += 1
            continue
        subject_key = assertion.subject.logical_key
        object_key = assertion.object.logical_key
        if subject_key in target_by_key and object_key in response_by_key:
            pair = (subject_key, object_key)
        elif object_key in target_by_key and subject_key in response_by_key:
            pair = (object_key, subject_key)
        else:
            rejected += 1
            continue
        if pair[0] == pair[1]:
            rejected += 1
            continue
        pair_assertions.setdefault(pair, []).append(assertion)

    observations: list[DetectedObservation] = []
    for target in targets:
        for response in responses:
            pair = (target.entity.logical_key, response.entity.logical_key)
            evidence = pair_assertions.get(pair)
            if not evidence:
                continue
            ordered = tuple(sorted(
                evidence,
                key=lambda item: (
                    _RELATIONSHIP_ORDER[item.relationship_type], item.assertion_id,
                ),
            ))
            relationship_types = _unique(
                item.relationship_type.value for item in ordered
            )
            observations.append(_detected(
                context,
                observation_type=ObservationType.TARGET_RESPONSE_RELATIONSHIP,
                discriminator=f"relationship:{pair[0]}:{pair[1]}",
                members=(target.entity, response.entity),
                basis=(
                    f"{len(ordered)} explicit relationship assertion(s) connect the current "
                    "perturbation target and returned response entity; relationship types remain distinct."
                ),
                source_fields=("relationship_assertions", *relationship_types),
                created_at=created_at,
                limitations=(
                    "Relationship types are not interchangeable: co-mention is not interaction, "
                    "shared annotation is not a network edge, and none alone establishes causality."
                ),
                details={
                    "relationship_types": relationship_types,
                    "assertion_ids": tuple(item.assertion_id for item in ordered),
                    "target": target.entity.logical_key,
                    "response": response.entity.logical_key,
                },
                assertions=ordered,
            ))
    return tuple(observations), rejected


def detect_observations(
    context: LocalResearchContext,
    relationship_assertions: Iterable[EvidenceAssertion] | Any = (),
    *,
    config: ResearchConfig = DEFAULT_RESEARCH_CONFIG,
    created_at: datetime | str | None = None,
) -> ObservationDetectionResult:
    """Detect exactly the seven V1 pattern types from immutable local inputs."""

    if not isinstance(context, LocalResearchContext):
        raise TypeError("context must be a LocalResearchContext")
    if not isinstance(config, ResearchConfig):
        raise TypeError("config must be a ResearchConfig")
    resolved_created_at = _normalized_created_at(context, created_at)
    responses = _responses(context)
    explicit, relationship_notices, prior_rejected = _relationship_input(relationship_assertions)

    observations: list[DetectedObservation] = []
    observations.extend(_family_observations(
        context, responses, config=config, created_at=resolved_created_at,
    ))
    observations.extend(_functional_observations(
        context, responses, created_at=resolved_created_at,
    ))
    localization = _localization_observation(
        context, responses, created_at=resolved_created_at,
    )
    if localization is not None:
        observations.append(localization)
    contrast = _positive_negative_observation(
        context, responses, created_at=resolved_created_at,
    )
    if contrast is not None:
        observations.append(contrast)
    gateway = _gateway_observation(
        context, responses, created_at=resolved_created_at,
    )
    if gateway is not None:
        observations.append(gateway)
    fdr = _fdr_observation(
        context, responses, created_at=resolved_created_at,
    )
    if fdr is not None:
        observations.append(fdr)
    relationship_observations, rejected = _relationship_observations(
        context, responses, explicit, created_at=resolved_created_at,
    )
    observations.extend(relationship_observations)

    notices = list(relationship_notices)
    total_rejected = prior_rejected + rejected
    if total_rejected:
        notices.append(
            f"{total_rejected} relationship assertion(s) were not exact target-response pairs and were ignored."
        )
    unresolved_response_count = sum(
        1 for item in context.unresolved_entities
        if any(role in _RESPONSE_ROLES for role in item.roles)
    )
    if unresolved_response_count:
        notices.append(
            f"{unresolved_response_count} ambiguity-safe unresolved response entity record(s) "
            "were excluded from entity-based observations without guessing an identity."
        )
    if not responses:
        notices.append("No resolved returned response entities are available for local observation detection.")
    return ObservationDetectionResult(
        simulation_id=context.simulation_id,
        snapshot_id=context.snapshot_id,
        observations=tuple(observations),
        notices=tuple(notices),
        rejected_relationship_count=total_rejected,
    )


# Concise architecture alias.
detect_local_observations = detect_observations


__all__ = [
    "DetectedObservation",
    "ObservationDetectionResult",
    "detect_local_observations",
    "detect_observations",
]
