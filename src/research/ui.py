"""Pure offline projections and an injected Research Explorer renderer.

The module imports no UI framework and performs no file, graph, provider,
session, or network access.  It projects only an already-built offline bundle;
the caller supplies a Streamlit-like object when rendering is requested.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from enum import Enum
import csv
import hashlib
import io
from typing import TYPE_CHECKING, Any
import zipfile

import pandas as pd

from .config import DEFAULT_RESEARCH_CONFIG, ResearchConfig
from .contracts import format_metric
from .local_context import LocalEntityContext, LocalResearchContext
from .models import EvidenceAssertion, ObservationType, ProvenanceRecord, RelationshipType
from .providers import FamilyMembership
from .semantic_adapter import ResultRole

if TYPE_CHECKING:
    from .facade import OfflineResearchBundle


class ExplorerSection(str, Enum):
    OVERVIEW = "Overview"
    FAMILIES = "Families"
    FUNCTIONS = "Functions"
    RELATIONSHIPS = "Relationships"


class ResearchExportFormat(str, Enum):
    XLSX = "Excel (.xlsx)"
    CSV_ZIP = "CSV paketi (.zip)"


_SECTIONS = tuple(section.value for section in ExplorerSection)
_SECTION_LABELS = {
    ExplorerSection.OVERVIEW: "Genel görünüm",
    ExplorerSection.FAMILIES: "Aileler",
    ExplorerSection.FUNCTIONS: "Fonksiyonlar",
    ExplorerSection.RELATIONSHIPS: "İlişkiler",
}
_COUNT_LABELS = {
    "Detected observations": "Saptanan gözlem",
    "Family patterns": "Aile örüntüsü",
    "Functional patterns": "Fonksiyon örüntüsü",
    "Compartment Bottleneck overlaps": "Bölme darboğazı örtüşmesi",
    "Relationships": "İlişki",
}
_GO_TYPES = frozenset({
    "go_biological_process", "go_cellular_component", "go_molecular_function",
})
_REACTOME_TYPES = frozenset({
    "reactome", "reactome_pathway", "reactome_annotation", "reactome_context",
})
_NETWORK_RELATIONSHIPS = frozenset({
    RelationshipType.PPI_EDGE,
    RelationshipType.SHARED_COMMUNITY,
    RelationshipType.NETWORK_PATH,
})
_LITERATURE_RELATIONSHIPS = frozenset({
    RelationshipType.CO_MENTIONED_IN_PUBLICATION,
})
_RELATIONSHIP_SUMMARIES = {
    RelationshipType.PPI_EDGE: "An existing PPI edge record covers this exact target/response pair.",
    RelationshipType.DIRECTED_REGULATION: (
        "An existing directed-regulation record covers this exact pair and orientation."
    ),
    RelationshipType.SHARED_PATHWAY: "The pair shares an exact existing pathway annotation.",
    RelationshipType.SHARED_GO: "The pair shares an exact existing GO annotation.",
    RelationshipType.SHARED_FAMILY: "The pair shares an authoritative family annotation.",
    RelationshipType.SHARED_COMPARTMENT: "The pair shares an exact localization annotation.",
    RelationshipType.SHARED_COMMUNITY: "The pair has the same recorded Sophiark community identifier.",
    RelationshipType.NETWORK_PATH: "An already-computed current Sophiark network path is recorded for the pair.",
    RelationshipType.CO_MENTIONED_IN_PUBLICATION: (
        "The pair is co-mentioned in an existing publication record; "
        "this does not establish interaction."
    ),
}


@dataclass(frozen=True, slots=True)
class ExplorerMember:
    identifier: str
    label: str
    roles: tuple[str, ...] = ()
    directions: tuple[str, ...] = ()
    fdr_supported: str = "—"
    localizations: tuple[str, ...] = ()
    mygene_name: str = "—"
    go_biological_processes: tuple[str, ...] = ()
    go_molecular_functions: tuple[str, ...] = ()
    go_cellular_components: tuple[str, ...] = ()
    shared_functions: tuple[str, ...] = ()
    gateway: str = "—"

    def __post_init__(self) -> None:
        for field_name in (
            "roles", "directions", "localizations", "go_biological_processes",
            "go_molecular_functions", "go_cellular_components", "shared_functions",
        ):
            object.__setattr__(self, field_name, tuple(getattr(self, field_name)))


@dataclass(frozen=True, slots=True)
class ExplorerItem:
    item_id: str
    title: str
    kind: str
    summary: str
    fields: tuple[tuple[str, str], ...] = ()
    members: tuple[ExplorerMember, ...] = ()
    sources: tuple[str, ...] = ()
    provenance: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "fields", tuple(tuple(pair) for pair in self.fields))
        object.__setattr__(self, "members", tuple(self.members))
        object.__setattr__(self, "sources", tuple(self.sources))
        object.__setattr__(self, "provenance", tuple(self.provenance))


@dataclass(frozen=True, slots=True)
class ResearchExplorerView:
    selected_section: ExplorerSection
    counts: tuple[tuple[str, int], ...] = ()
    overview_members: tuple[ExplorerMember, ...] = ()
    items: tuple[ExplorerItem, ...] = ()
    selected_item: ExplorerItem | None = None
    status_badges: tuple[str, ...] = ()
    notices: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    unavailable_message: str | None = None
    truncated: bool = False

    def __post_init__(self) -> None:
        object.__setattr__(self, "selected_section", ExplorerSection(self.selected_section))
        object.__setattr__(self, "counts", tuple(tuple(pair) for pair in self.counts))
        for field_name in ("overview_members", "items", "status_badges", "notices", "errors"):
            object.__setattr__(self, field_name, tuple(getattr(self, field_name)))


def _unique(values: Iterable[Any]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(str(value).strip() for value in values if str(value).strip()))


def _items(value: Any) -> tuple[Any, ...]:
    return tuple(value) if isinstance(value, (tuple, list, set, frozenset)) else (value,)


def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value == 1
    return str(value).strip().casefold() in {"true", "1", "yes", "evet", "✓"}


def _term(value: Any) -> tuple[str, str] | None:
    if isinstance(value, Mapping):
        identifier = next((value.get(key) for key in (
            "id", "go_id", "reactome_id", "pathway_id", "term_id",
        ) if value.get(key)), None)
        label = next((value.get(key) for key in (
            "term", "name", "pathway_name",
        ) if value.get(key)), identifier)
        if identifier:
            return str(identifier).strip(), str(label or identifier).strip()
        if label:
            text = str(label).strip()
            return text.casefold(), text
        return None
    text = str(value or "").strip()
    if not text or text.casefold() in {"nan", "none", "unknown", "—"}:
        return None
    return text.casefold(), text


def _provenance_labels(records: Iterable[ProvenanceRecord]) -> tuple[str, ...]:
    return _unique(
        " · ".join(part for part in (record.source, record.version, record.locator) if part)
        for record in records
    )


def _status_badges(bundle: Any, context: LocalResearchContext | None) -> tuple[str, ...]:
    badges: list[str] = []
    if context is not None:
        for status in context.source_statuses:
            state = getattr(status.status, "value", str(status.status))
            badges.append(f"{status.source_id}: {state}")
        recorded = {status.source_id.casefold() for status in context.source_statuses}
        for record in context.provenance:
            if record.source.casefold() not in recorded:
                badges.append(f"{record.source}: recorded")
    if not badges and getattr(bundle, "errors", ()):
        badges.append("offline research context: unavailable")
    signaling = getattr(bundle, "signaling_result", None)
    if signaling is not None:
        state = getattr(getattr(signaling, "status", None), "value", getattr(signaling, "status", "UNAVAILABLE"))
        badges.append(f"omnipath_signaling: {state}")
    return _unique(badges)


def _signaling_export_rows(bundle: Any) -> tuple[dict[str, Any], ...]:
    result = getattr(bundle, "signaling_result", None)
    rows = tuple(getattr(result, "rows", ()) or ())
    return tuple({
        "perturbation_target": item.perturbation_target,
        "response_entity": item.response_entity,
        "response_type": item.response_type,
        "network_response_value": item.network_response_value,
        "directed_relation": item.directed_relation,
        "directed_distance": item.directed_distance,
        "sign": item.sign,
        "example_path": " → ".join(item.example_path),
        "path_count_within_depth": item.path_count_within_depth,
        "positive_path_count": item.positive_path_count,
        "negative_path_count": item.negative_path_count,
        "unsigned_path_count": item.unsigned_path_count,
        "uncertain_path_count": item.uncertain_path_count,
        "alternative_route_candidate": item.alternative_route_candidate,
        "tissue_support": item.tissue_support,
        "evidence_resources": " | ".join(item.evidence_resources),
        "references": " | ".join(item.references),
        "dataset_version": item.dataset_version,
    } for item in rows)


def _convergence_export_rows(bundle: Any) -> tuple[dict[str, Any], ...]:
    result = getattr(bundle, "signaling_result", None)
    records = tuple(getattr(result, "convergence", ()) or ())
    return tuple({
        "convergence_node": item.convergence_node,
        "converging_targets": " | ".join(item.converging_targets),
        "distance_per_target": json.dumps(dict(item.distances), sort_keys=True),
        "signed_relation_per_target": json.dumps(dict(item.signed_relations), sort_keys=True),
        "evidence_resources": " | ".join(item.resources),
        "references": " | ".join(item.references),
    } for item in records)


def _annotation_values(entity: LocalEntityContext, annotation_type: str) -> tuple[Any, ...]:
    return tuple(
        value
        for annotation in entity.annotations
        if annotation.annotation_type == annotation_type
        for value in _items(annotation.value)
    )


def _annotation_text(entity: LocalEntityContext, annotation_type: str, *, limit: int = 4) -> tuple[str, ...]:
    """Return recorded local annotation text without deriving a new claim."""

    values: list[str] = []
    for value in _annotation_values(entity, annotation_type):
        identity = _term(value)
        if identity is not None:
            values.append(identity[1])
            continue
        text = str(value or "").strip()
        if text and text.casefold() not in {"nan", "none", "unknown", "—"}:
            values.append(text)
    return _unique(values)[:limit]


def _member(
    entity: LocalEntityContext,
    *,
    shared_functions: Iterable[str] = (),
) -> ExplorerMember:
    gateway = any(_truthy(value) for value in _annotation_values(entity, "gateway"))
    localizations = _unique(_annotation_values(entity, "localization"))
    return ExplorerMember(
        identifier=entity.entity.canonical_id,
        label=entity.entity.symbol or entity.entity.canonical_id,
        roles=tuple(role.value for role in entity.roles),
        directions=tuple(direction.value for direction in entity.directions),
        fdr_supported=format_metric(
            "significant_redistribution",
            ResultRole.FDR_SUPPORTED in entity.roles,
        ),
        localizations=localizations,
        mygene_name=", ".join(_annotation_text(entity, "gene_name", limit=1)) or "—",
        go_biological_processes=_annotation_text(entity, "go_biological_process"),
        go_molecular_functions=_annotation_text(entity, "go_molecular_function"),
        go_cellular_components=_annotation_text(entity, "go_cellular_component"),
        shared_functions=_unique(shared_functions),
        gateway=format_metric("Gümrük_Kapisi", gateway),
    )


def _entity_identifiers(entity: LocalEntityContext) -> frozenset[str]:
    return frozenset(
        str(value).strip().casefold()
        for value in (
            entity.entity.canonical_id,
            entity.entity.symbol,
            *(record.canonical_id for record in entity.source_records),
            *(record.symbol for record in entity.source_records),
        )
        if value is not None and str(value).strip()
    )


def _shared_functions(context: LocalResearchContext, entity: LocalEntityContext, cap: int) -> tuple[str, ...]:
    identifiers = _entity_identifiers(entity)
    labels: list[str] = []
    for row in context.functional_context:
        members = {str(member).strip().casefold() for member in row.members}
        if identifiers & members:
            labels.append(str(row.term or row.source_library or row.pathway_family or "").strip())
    for annotation in entity.annotations:
        if annotation.annotation_type not in (_GO_TYPES | _REACTOME_TYPES):
            continue
        for value in _items(annotation.value):
            identity = _term(value)
            if identity:
                labels.append(identity[1])
    return _unique(labels)[:cap]


def _observation_item(item: Any) -> ExplorerItem:
    observation = getattr(item, "observation", item)
    observation_type = getattr(observation.type, "value", str(observation.type))
    members = tuple(
        ExplorerMember(
            identifier=member.canonical_id,
            label=member.symbol or member.canonical_id,
        )
        for member in observation.members
    )
    return ExplorerItem(
        item_id=observation.observation_id,
        title=observation_type.replace("_", " ").title(),
        kind=observation_type,
        summary=observation.basis,
        fields=(("Scope", observation.scope.selection_mode), ("Limitations", observation.limitations)),
        members=members,
        sources=tuple(observation.source_fields),
    )


def _overview(
    bundle: Any,
    context: LocalResearchContext | None,
    limit: int,
) -> tuple[tuple[tuple[str, int], ...], tuple[ExplorerMember, ...], tuple[ExplorerItem, ...], bool]:
    observation_result = getattr(bundle, "observation_result", None)
    observations = tuple(getattr(observation_result, "observations", ()) or ())
    relationships = tuple(getattr(getattr(bundle, "relationships", None), "assertions", ()) or ())
    counts_by_type = {
        observation_type: sum(
            getattr(item, "type", getattr(getattr(item, "observation", None), "type", None))
            is observation_type
            for item in observations
        )
        for observation_type in ObservationType
    }
    counts = (
        ("Detected observations", len(observations)),
        ("Family patterns", counts_by_type[ObservationType.FAMILY_COOCCURRENCE]),
        ("Functional patterns", counts_by_type[ObservationType.FUNCTIONAL_COOCCURRENCE]),
        ("Compartment Bottleneck overlaps", counts_by_type[ObservationType.GATEWAY_RESPONSE_OVERLAP]),
        ("Relationships", len(relationships)),
    )
    overview_members = () if context is None else tuple(
        _member(
            entity,
            shared_functions=_shared_functions(context, entity, 10),
        )
        for entity in context.entities[:50]
    )
    return (
        counts,
        overview_members,
        tuple(_observation_item(item) for item in observations[:limit]),
        len(observations) > limit,
    )


def _family_items(
    context: LocalResearchContext,
    *,
    selected_id: str | None,
    config: ResearchConfig,
    include_all: bool = False,
) -> tuple[tuple[ExplorerItem, ...], ExplorerItem | None, bool]:
    grouped: dict[str, dict[str, Any]] = {}
    for entity in context.entities:
        for annotation in entity.annotations:
            if annotation.annotation_type != "authoritative_family_membership":
                continue
            for value in _items(annotation.value):
                if not isinstance(value, FamilyMembership) or value.taxon_id != entity.entity.taxon_id:
                    continue
                identity = value.group_id.strip() or value.group_name.strip().casefold()
                item_id = f"family:{entity.entity.taxon_id}:{identity}"
                group = grouped.setdefault(item_id, {
                    "id": value.group_id, "name": value.group_name,
                    "members": [], "provenance": [], "sources": [],
                })
                if entity not in group["members"]:
                    group["members"].append(entity)
                group["provenance"].append(value.provenance)
                group["sources"].append(annotation.source)

    minimum = config.presentation_limits.family_min_members
    eligible = [pair for pair in grouped.items() if len(pair[1]["members"]) >= minimum]
    limit = len(eligible) if include_all else config.presentation_limits.max_family_rows
    visible = eligible[:limit]
    summaries = tuple(
        ExplorerItem(
            item_id=item_id,
            title=group["name"] or group["id"],
            kind="AUTHORITATIVE_FAMILY",
            summary=f"{len(group['members'])} current member(s) share this authoritative family annotation.",
            fields=(("Family ID", group["id"] or "—"), ("Member count", str(len(group["members"])))),
            sources=_unique(group["sources"]),
            provenance=_provenance_labels(group["provenance"]),
        )
        for item_id, group in visible
    )
    detail: ExplorerItem | None = None
    if selected_id in dict(visible):
        group = dict(visible)[selected_id]
        member_cap = config.presentation_limits.max_family_rows
        members = tuple(
            _member(
                entity,
                shared_functions=_shared_functions(
                    context, entity, min(config.presentation_limits.max_function_rows, 10),
                ),
            )
            for entity in group["members"][:member_cap]
        )
        detail = ExplorerItem(
            item_id=selected_id,
            title=group["name"] or group["id"],
            kind="AUTHORITATIVE_FAMILY",
            summary="Exact authoritative family membership from the recorded source.",
            fields=(
                ("Family ID", group["id"] or "—"),
                ("Displayed members", str(len(members))),
            ),
            members=members,
            sources=_unique(group["sources"]),
            provenance=_provenance_labels(group["provenance"]),
        )
    return summaries, detail, len(eligible) > limit


def _function_rows(context: LocalResearchContext) -> list[dict[str, Any]]:
    rows: dict[str, dict[str, Any]] = {}
    entity_lookup: dict[str, LocalEntityContext] = {}
    for entity in context.entities:
        for identifier in _entity_identifiers(entity):
            entity_lookup.setdefault(identifier, entity)

    for row in context.functional_context:
        label = str(row.term or row.source_library or row.pathway_family or "").strip()
        if not label:
            continue
        identity = "|".join((row.kind.value, row.source_library or "", row.pathway_family or "", label))
        item_id = f"function:{hashlib.sha256(identity.encode('utf-8')).hexdigest()[:16]}"
        members: list[LocalEntityContext] = []
        raw_members: list[str] = []
        for identifier in row.members:
            entity = entity_lookup.get(str(identifier).strip().casefold())
            if entity is not None and entity not in members:
                members.append(entity)
            elif entity is None:
                raw_members.append(str(identifier))
        rows[item_id] = {
            "kind": row.kind.value, "term": label, "source": row.source_library or "Sophiark",
            "family": row.pathway_family, "members": members, "raw_members": raw_members,
            "provenance": list(context.provenance),
        }

    for entity in context.entities:
        for annotation in entity.annotations:
            if annotation.annotation_type in _GO_TYPES:
                kind = "GO"
            elif annotation.annotation_type in _REACTOME_TYPES:
                kind = "Reactome"
            else:
                continue
            for value in _items(annotation.value):
                term = _term(value)
                if term is None:
                    continue
                identity = f"{kind}|{annotation.source}|{term[0]}"
                item_id = f"function:{hashlib.sha256(identity.encode('utf-8')).hexdigest()[:16]}"
                row = rows.setdefault(item_id, {
                    "kind": kind, "term": term[1], "source": annotation.source,
                    "family": annotation.annotation_type, "members": [], "raw_members": [],
                    "provenance": [],
                })
                if entity not in row["members"]:
                    row["members"].append(entity)
                row["provenance"].extend(annotation.provenance)
    return list(rows.values())


def _function_items(
    context: LocalResearchContext,
    *,
    selected_id: str | None,
    config: ResearchConfig,
    include_all: bool = False,
    _prepared_rows: list[dict[str, Any]] | None = None,
    _prepared_summaries: tuple[ExplorerItem, ...] | None = None,
) -> tuple[tuple[ExplorerItem, ...], ExplorerItem | None, bool]:
    rows = _function_rows(context) if _prepared_rows is None else _prepared_rows
    limit = len(rows) if include_all else min(config.presentation_limits.max_function_rows, 50)
    visible = rows[:limit]

    def item_id(row: dict[str, Any]) -> str:
        identity = f"{row['kind']}|{row['source']}|{row['family']}|{row['term']}"
        return f"function:{hashlib.sha256(identity.encode('utf-8')).hexdigest()[:16]}"

    summaries = _prepared_summaries if _prepared_summaries is not None else tuple(
        ExplorerItem(
            item_id=item_id(row),
            title=row["term"],
            kind=row["kind"],
            summary=f"{len(row['members']) + len(row['raw_members'])} current Sophiark member(s).",
            fields=(("Source", row["source"]), ("Context", str(row["family"] or "—"))),
            sources=(row["source"],),
            provenance=_provenance_labels(row["provenance"]),
        )
        for row in visible
    )
    detail: ExplorerItem | None = None
    selected_row = next((row for row in visible if item_id(row) == selected_id), None)
    if selected_row is not None:
        members = tuple(_member(entity) for entity in selected_row["members"][:limit])
        member_roles = lambda role: ", ".join(
            member.label for member in members if role in member.roles
        ) or "—"
        target_members = member_roles(ResultRole.PERTURBATION_TARGET.value)
        detail = ExplorerItem(
            item_id=selected_id or "",
            title=selected_row["term"],
            kind=selected_row["kind"],
            summary="Existing exact functional context; no new enrichment was calculated.",
            fields=(
                ("Source", selected_row["source"]),
                ("Current members", ", ".join(member.label for member in members) or "—"),
                ("Positive", member_roles(ResultRole.POSITIVE_REDISTRIBUTION.value)),
                ("Loss", member_roles(ResultRole.NETWORK_LOSS.value)),
                ("Compartment Bottleneck", ", ".join(member.label for member in members if member.gateway == "✓") or "—"),
                ("FDR supported", member_roles(ResultRole.FDR_SUPPORTED.value)),
                ("Target", target_members),
            ),
            members=members,
            sources=(selected_row["source"],),
            provenance=_provenance_labels(selected_row["provenance"]),
        )
    return summaries, detail, len(rows) > limit


def _relationship_category(relationship_type: RelationshipType) -> str:
    if relationship_type in _NETWORK_RELATIONSHIPS:
        return "Sophiark network"
    if relationship_type in _LITERATURE_RELATIONSHIPS:
        return "Literature context"
    return "Biological annotations"


def _relationship_items(
    assertions: Iterable[EvidenceAssertion],
    *,
    selected_id: str | None,
    config: ResearchConfig,
    include_all: bool = False,
    _prepared_summaries: tuple[ExplorerItem, ...] | None = None,
) -> tuple[tuple[ExplorerItem, ...], ExplorerItem | None, bool]:
    all_assertions = tuple(assertions)
    limit = len(all_assertions) if include_all else min(config.presentation_limits.max_relationship_rows, 50)
    visible = all_assertions[:limit]

    def projection(assertion: EvidenceAssertion, *, detail: bool) -> ExplorerItem:
        relationship_type = assertion.relationship_type
        category = _relationship_category(relationship_type)
        subject = assertion.subject.symbol or assertion.subject.canonical_id
        object_label = (
            assertion.object.symbol or assertion.object.canonical_id
            if assertion.object is not None else "—"
        )
        summary = _RELATIONSHIP_SUMMARIES[relationship_type]
        provenance = _provenance_labels(assertion.provenance)
        return ExplorerItem(
            item_id=assertion.assertion_id,
            title=f"{relationship_type.value}: {subject} ↔ {object_label}",
            kind=relationship_type.value,
            summary=summary,
            fields=(
                ("Category", category),
                ("Relationship type", relationship_type.value),
                ("Subject", subject),
                ("Object", object_label),
            ) if detail else (("Category", category), ("Relationship type", relationship_type.value)),
            members=(
                ExplorerMember(assertion.subject.canonical_id, subject),
                *(() if assertion.object is None else (
                    ExplorerMember(assertion.object.canonical_id, object_label),
                )),
            ) if detail else (),
            sources=_unique(record.source for record in assertion.provenance),
            provenance=provenance,
        )

    summaries = _prepared_summaries if _prepared_summaries is not None else tuple(projection(assertion, detail=False) for assertion in visible)
    selected = next((item for item in visible if item.assertion_id == selected_id), None)
    detail = projection(selected, detail=True) if selected is not None else None
    return summaries, detail, len(all_assertions) > limit


def build_research_explorer_view(
    bundle: "OfflineResearchBundle",
    *,
    selected_section: ExplorerSection | str = ExplorerSection.OVERVIEW,
    selected_family_id: str | None = None,
    selected_function_id: str | None = None,
    selected_relationship_id: str | None = None,
    config: ResearchConfig | None = None,
) -> ResearchExplorerView:
    """Build one immutable, selected-section-only offline projection."""

    section = ExplorerSection(selected_section)
    resolved_config = config or getattr(bundle, "config", None) or DEFAULT_RESEARCH_CONFIG
    if not isinstance(resolved_config, ResearchConfig):
        raise TypeError("config must be a ResearchConfig")
    context = getattr(bundle, "local_context", None)
    if context is not None and not isinstance(context, LocalResearchContext):
        raise TypeError("bundle.local_context must be a LocalResearchContext")
    notices = _unique((
        *tuple(getattr(bundle, "notices", ()) or ()),
        *tuple(getattr(getattr(bundle, "observation_result", None), "notices", ()) or ()),
        *tuple(getattr(getattr(bundle, "relationships", None), "notices", ()) or ()),
    ))
    errors = _unique(getattr(bundle, "errors", ()) or ())
    common = {
        "selected_section": section,
        "status_badges": _status_badges(bundle, context),
        "notices": notices,
        "errors": errors,
    }

    if section is ExplorerSection.OVERVIEW:
        limit = resolved_config.presentation_limits.max_overview_observations
        counts, overview_members, items, truncated = _overview(bundle, context, limit)
        has_recorded_result = bool(items or any(value for _, value in counts))
        unavailable = (
            "Offline research context is unavailable."
            if errors and not has_recorded_result
            else None
        )
        overview_notices = notices
        if not has_recorded_result and not errors:
            overview_notices = _unique((
                *notices,
                "No deterministic observations were detected in the current returned set.",
            ))
        return ResearchExplorerView(
            **{**common, "notices": overview_notices}, counts=counts,
            overview_members=overview_members, items=items,
            unavailable_message=unavailable, truncated=truncated,
        )
    if context is None:
        message = {
            ExplorerSection.FAMILIES: "Authoritative family context is unavailable.",
            ExplorerSection.FUNCTIONS: "GO/KEGG/Reactome functional context is unavailable.",
            ExplorerSection.RELATIONSHIPS: "Relationship evidence is unavailable.",
        }[section]
        return ResearchExplorerView(**common, unavailable_message=message)
    if section is ExplorerSection.FAMILIES:
        items, detail, truncated = _family_items(
            context, selected_id=selected_family_id, config=resolved_config,
        )
        return ResearchExplorerView(
            **common, items=items, selected_item=detail, truncated=truncated,
            unavailable_message=None if items else "Authoritative family context is unavailable.",
        )
    if section is ExplorerSection.FUNCTIONS:
        items, detail, truncated = _function_items(
            context, selected_id=selected_function_id, config=resolved_config,
        )
        return ResearchExplorerView(
            **common, items=items, selected_item=detail, truncated=truncated,
            unavailable_message=None if items else "GO/KEGG/Reactome functional context is unavailable.",
        )
    assertions = tuple(getattr(getattr(bundle, "relationships", None), "assertions", ()) or ())
    items, detail, truncated = _relationship_items(
        assertions, selected_id=selected_relationship_id, config=resolved_config,
    )
    return ResearchExplorerView(
        **common, items=items, selected_item=detail, truncated=truncated,
        unavailable_message=None if items else "Relationship evidence is unavailable.",
    )


def _export_item_row(item: ExplorerItem, *, context_label: str | None = None) -> dict[str, str]:
    row = {
        "ID": item.item_id,
        "Başlık": item.title,
        "Tür": item.kind,
        "Özet": item.summary,
        "Kaynaklar": " | ".join(item.sources) or "—",
        "Provenance": " | ".join(item.provenance) or "—",
    }
    if context_label is not None:
        row = {"Bağlam": context_label, **row}
    for label, value in item.fields:
        row.setdefault(label, value)
    return row


def _export_member_rows(
    items: Iterable[ExplorerItem],
    detail_for,
    *,
    context_label: str,
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in items:
        detail = detail_for(item.item_id)
        if detail is None:
            continue
        for member in detail.members:
            # _member_row intentionally keeps Compartment Bottleneck as the final column.
            rows.append({context_label: item.title, **_member_row(member)})
    return rows


def build_research_export_tables(
    bundle: "OfflineResearchBundle",
    *,
    config: ResearchConfig | None = None,
) -> dict[str, tuple[dict[str, str], ...]]:
    """Build all existing Research Explorer projections for read-only export.

    The export intentionally contains only Research Context projections. It
    neither reads the raw scientific tables nor triggers providers, refreshes,
    files, graph work, or any scientific recalculation.
    """

    resolved_config = config or getattr(bundle, "config", None) or DEFAULT_RESEARCH_CONFIG
    if not isinstance(resolved_config, ResearchConfig):
        raise TypeError("config must be a ResearchConfig")
    context = getattr(bundle, "local_context", None)
    if context is not None and not isinstance(context, LocalResearchContext):
        raise TypeError("bundle.local_context must be a LocalResearchContext")

    observations = tuple(getattr(getattr(bundle, "observation_result", None), "observations", ()) or ())
    counts, overview_members, observation_items, _ = _overview(bundle, context, len(observations))
    tables: dict[str, tuple[dict[str, str], ...]] = {
        "Overview": tuple({"Metrik": label, "Değer": str(value)} for label, value in counts),
        "General View": tuple(_member_row(member) for member in overview_members),
        "Observations": tuple(_export_item_row(item) for item in observation_items),
        "Directed Signaling Context": _signaling_export_rows(bundle),
        "Directed Convergence": _convergence_export_rows(bundle),
    }

    if context is not None:
        family_items, _, _ = _family_items(
            context,
            selected_id=None,
            config=resolved_config,
            include_all=True,
        )
        tables["Families"] = tuple(_export_item_row(item) for item in family_items)
        tables["Family Members"] = tuple(_export_member_rows(
            family_items,
            lambda item_id: _family_items(
                context, selected_id=item_id, config=resolved_config, include_all=True,
            )[1],
            context_label="Aile",
        ))

        function_rows = _function_rows(context)
        function_items, _, _ = _function_items(
            context,
            selected_id=None,
            config=resolved_config,
            include_all=True,
            _prepared_rows=function_rows,
        )
        tables["Functions"] = tuple(_export_item_row(item) for item in function_items)
        tables["Function Members"] = tuple(_export_member_rows(
            function_items,
            lambda item_id: _function_items(
                context, selected_id=item_id, config=resolved_config, include_all=True,
                _prepared_rows=function_rows, _prepared_summaries=function_items,
            )[1],
            context_label="Fonksiyon",
        ))
    else:
        tables["Families"] = ()
        tables["Family Members"] = ()
        tables["Functions"] = ()
        tables["Function Members"] = ()

    assertions = tuple(getattr(getattr(bundle, "relationships", None), "assertions", ()) or ())
    relationship_items, _, _ = _relationship_items(
        assertions,
        selected_id=None,
        config=resolved_config,
        include_all=True,
    )
    tables["Relationships"] = tuple(_export_item_row(item) for item in relationship_items)
    tables["Relationship Members"] = tuple(_export_member_rows(
        relationship_items,
        lambda item_id: _relationship_items(
            assertions, selected_id=item_id, config=resolved_config, include_all=True,
            _prepared_summaries=relationship_items,
        )[1],
        context_label="İlişki",
    ))
    tables["Sources and Notices"] = tuple(
        {"Tür": "Kaynak durumu", "Mesaj": value}
        for value in _status_badges(bundle, context)
    ) + tuple(
        {"Tür": "Bildirim", "Mesaj": value}
        for value in _unique((
            *tuple(getattr(bundle, "notices", ()) or ()),
            *tuple(getattr(getattr(bundle, "observation_result", None), "notices", ()) or ()),
            *tuple(getattr(getattr(bundle, "relationships", None), "notices", ()) or ()),
            *tuple(getattr(bundle, "errors", ()) or ()),
        ))
    ) + tuple(
        {"Tür": "OmniPath signaling metadata", "Mesaj": f"{key}: {value}"}
        for key, value in dict(getattr(getattr(bundle, "signaling_result", None), "metadata", {}) or {}).items()
    )
    return tables


def _table_columns(rows: Iterable[Mapping[str, str]]) -> tuple[str, ...]:
    columns: list[str] = []
    for row in rows:
        for key in row:
            if key not in columns:
                columns.append(key)
    return tuple(columns) or ("Durum",)


def _csv_zip_payload(tables: Mapping[str, tuple[dict[str, str], ...]]) -> bytes:
    buffer = io.BytesIO()
    with zipfile.ZipFile(buffer, mode="w", compression=zipfile.ZIP_DEFLATED) as archive:
        for table_name, rows in tables.items():
            columns = _table_columns(rows)
            text = io.StringIO(newline="")
            writer = csv.DictWriter(text, fieldnames=columns, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
            archive.writestr(f"{table_name.casefold().replace(' ', '_')}.csv", text.getvalue().encode("utf-8-sig"))
    return buffer.getvalue()


def _xlsx_payload(tables: Mapping[str, tuple[dict[str, str], ...]]) -> bytes:
    failures: list[str] = []
    for engine in ("openpyxl", "xlsxwriter"):
        try:
            buffer = io.BytesIO()
            with pd.ExcelWriter(buffer, engine=engine) as writer:
                for table_name, rows in tables.items():
                    columns = _table_columns(rows)
                    pd.DataFrame(rows, columns=columns).to_excel(
                        writer,
                        index=False,
                        sheet_name=table_name[:31],
                    )
            return buffer.getvalue()
        except Exception as error:
            failures.append(f"{engine}: {type(error).__name__}")
    raise RuntimeError("Research Explorer Excel export could not be created: " + "; ".join(failures))


def research_explorer_download_payload(
    bundle: "OfflineResearchBundle",
    format: ResearchExportFormat | str = ResearchExportFormat.XLSX,
    *,
    config: ResearchConfig | None = None,
) -> tuple[str, str, bytes]:
    """Return one multi-table download payload without performing I/O."""

    selected_format = ResearchExportFormat(format)
    tables = build_research_export_tables(bundle, config=config)
    snapshot_id = str(getattr(getattr(bundle, "snapshot", None), "snapshot_id", "research"))
    digest = hashlib.sha256(snapshot_id.encode("utf-8")).hexdigest()[:12]
    stem = f"research-explorer-{digest}"
    if selected_format is ResearchExportFormat.XLSX:
        return (
            f"{stem}.xlsx",
            "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            _xlsx_payload(tables),
        )
    return f"{stem}.zip", "application/zip", _csv_zip_payload(tables)


def _emit(ui: Any, method: str, *args: Any, **kwargs: Any) -> Any:
    function = getattr(ui, method, None)
    if not callable(function):
        return None
    try:
        return function(*args, **kwargs)
    except TypeError:
        return function(*args)


def _navigation(ui: Any) -> ExplorerSection:
    selected = None
    labels = tuple(_SECTION_LABELS[section] for section in ExplorerSection)
    if callable(getattr(ui, "segmented_control", None)):
        selected = _emit(
            ui, "segmented_control", "Araştırma bölümü", labels,
            default=_SECTION_LABELS[ExplorerSection.OVERVIEW], key="research_explorer_section",
        )
    elif callable(getattr(ui, "radio", None)):
        selected = _emit(
            ui, "radio", "Araştırma bölümü", labels,
            index=0, horizontal=True, key="research_explorer_section",
        )
    reverse = {label: section for section, label in _SECTION_LABELS.items()}
    if selected in _SECTIONS:
        return ExplorerSection(selected)
    return reverse.get(selected, ExplorerSection.OVERVIEW)


def _member_row(member: ExplorerMember) -> dict[str, str]:
    """One stable, human-readable row; Gateway is intentionally final."""

    return {
        "Gen": member.label,
        "ENSP": member.identifier,
        "MyGene Adı": member.mygene_name,
        "Roller": ", ".join(member.roles) or "—",
        "Yön": ", ".join(member.directions) or "—",
        "FDR": member.fdr_supported,
        "Lokalizasyon": ", ".join(member.localizations) or "—",
        "GO Biyolojik Süreç": "; ".join(member.go_biological_processes) or "—",
        "GO Moleküler İşlev": "; ".join(member.go_molecular_functions) or "—",
        "GO Hücresel Bileşen": "; ".join(member.go_cellular_components) or "—",
        "Paylaşılan Fonksiyonlar": "; ".join(member.shared_functions) or "—",
        "Compartment Bottleneck": member.gateway,
    }


def _render_member_table(ui: Any, members: Iterable[ExplorerMember], *, label: str) -> None:
    rows = [_member_row(member) for member in members]
    if not rows:
        return
    _emit(ui, "caption", label)
    _emit(
        ui,
        "dataframe",
        rows,
        hide_index=True,
        width="stretch",
    )


def render_research_explorer_download(
    bundle: "OfflineResearchBundle",
    *,
    ui: Any,
    config: ResearchConfig | None = None,
    format: ResearchExportFormat | str | None = None,
) -> None:
    """Render one explicit all-output download control through the injected UI."""

    selected = format
    if selected is None:
        selected = _emit(
            ui,
            "radio",
            "Tüm Research Explorer çıktıları için format",
            tuple(item.value for item in ResearchExportFormat),
            horizontal=True,
            key="research_explorer_download_format",
        ) or ResearchExportFormat.XLSX.value
    selected_format = ResearchExportFormat(selected)
    try:
        filename, mime, data = research_explorer_download_payload(
            bundle,
            selected_format,
            config=config,
        )
    except RuntimeError as error:
        if selected_format is not ResearchExportFormat.XLSX:
            raise
        _emit(ui, "warning", "Excel export is unavailable; CSV package is offered instead.")
        filename, mime, data = research_explorer_download_payload(
            bundle,
            ResearchExportFormat.CSV_ZIP,
            config=config,
        )
    _emit(
        ui,
        "download_button",
        "Tüm Research Explorer çıktılarını indir",
        data=data,
        file_name=filename,
        mime=mime,
        key="research_explorer_download_all",
        width="stretch",
    )


def _render_item(ui: Any, item: ExplorerItem, *, show_members: bool = True) -> None:
    _emit(ui, "markdown", f"#### {item.title}")
    _emit(ui, "write", item.summary)
    for label, value in item.fields:
        _emit(ui, "write", f"{label}: {value}")
    if show_members:
        _render_member_table(ui, item.members, label="Üyeler · yerel açıklamalar")


def render_research_explorer(
    bundle: "OfflineResearchBundle",
    *,
    ui: Any,
    config: ResearchConfig | None = None,
    selected_section: ExplorerSection | str | None = None,
    include_download: bool = True,
) -> ResearchExplorerView:
    """Render only the selected Explorer section through an injected UI object."""

    section = ExplorerSection(selected_section) if selected_section is not None else _navigation(ui)
    view = build_research_explorer_view(
        bundle, selected_section=section, config=config,
    )
    _emit(ui, "subheader", "Araştırma Gezgini")
    for error in view.errors[:6]:
        _emit(ui, "error", error)
    if view.unavailable_message:
        _emit(ui, "info", view.unavailable_message)
        if include_download:
            render_research_explorer_download(bundle, ui=ui, config=config)
        return view

    if section is ExplorerSection.OVERVIEW:
        columns = _emit(ui, "columns", len(view.counts))
        if isinstance(columns, (tuple, list)) and len(columns) >= len(view.counts):
            for column, (label, value) in zip(columns, view.counts):
                _emit(column, "metric", _COUNT_LABELS.get(label, label), value)
        else:
            for label, value in view.counts:
                _emit(ui, "metric", _COUNT_LABELS.get(label, label), value)
        _render_member_table(
            ui,
            view.overview_members,
            label="Genel görünüm · mevcut sonuç, MyGene ve GO bağlamı",
        )
        for item in view.items:
            _render_item(ui, item, show_members=False)
        signaling_rows = _signaling_export_rows(bundle)
        if signaling_rows:
            _emit(ui, "markdown", "#### Yönlü sinyalleme bağlamı")
            _emit(
                ui, "caption",
                "OmniPath direction/sign context is separate from Sophiark network redistribution; "
                "it does not establish activity, causality, compensation or tissue-specific signaling.",
            )
            visible = pd.DataFrame(signaling_rows).rename(columns={
                "response_entity": "Entity",
                "network_response_value": "Network Response",
                "directed_relation": "Relation",
                "directed_distance": "Distance",
                "sign": "Sign",
                "example_path": "Example Path",
                "evidence_resources": "Evidence",
            })
            columns = [
                "Entity", "Network Response", "Relation", "Distance", "Sign",
                "Example Path", "Evidence",
            ]
            _emit(ui, "dataframe", visible[columns], hide_index=True, width="stretch")
    else:
        if section is ExplorerSection.RELATIONSHIPS:
            for category in ("Sophiark network", "Biological annotations", "Literature context"):
                count = sum(dict(item.fields).get("Category") == category for item in view.items)
                _emit(ui, "write", f"{category}: {count}")
        option_ids = tuple(item.item_id for item in view.items)
        labels = {item.item_id: item.title for item in view.items}
        selected_id = _emit(
            ui, "selectbox", f"{_SECTION_LABELS[section]} içinde seçim", option_ids,
            format_func=lambda value: labels[value], key=f"research_{section.value.casefold()}_selection",
        ) or option_ids[0]
        selection_kwargs = {
            "selected_family_id": selected_id if section is ExplorerSection.FAMILIES else None,
            "selected_function_id": selected_id if section is ExplorerSection.FUNCTIONS else None,
            "selected_relationship_id": selected_id if section is ExplorerSection.RELATIONSHIPS else None,
        }
        view = build_research_explorer_view(
            bundle, selected_section=section, config=config, **selection_kwargs,
        )
        if view.selected_item is not None:
            _render_item(ui, view.selected_item)
    # Provider status, provenance and implementation notices remain available
    # in the downloadable audit tables. They are intentionally not shown in
    # the user-facing explorer because they do not explain a biological result.
    if include_download:
        render_research_explorer_download(bundle, ui=ui, config=config)
    return view


__all__ = [
    "ExplorerItem",
    "ExplorerMember",
    "ExplorerSection",
    "ResearchExportFormat",
    "ResearchExplorerView",
    "build_research_explorer_view",
    "build_research_export_tables",
    "research_explorer_download_payload",
    "render_research_explorer_download",
    "render_research_explorer",
]
