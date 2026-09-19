"""Local-only orchestration for immutable Research Context views.

This module consumes an already-frozen simulation snapshot, its semantic view,
the species-aware entity resolver, and explicitly injected local providers.  It
does not import Streamlit, inspect session state, read a graph, perform a
scientific calculation, discover project files, or make a network request.
"""

from __future__ import annotations

from dataclasses import dataclass, field
import math
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from .entity_resolution import EntityResolver
from .models import (
    Entity,
    EntityResolution,
    ProvenanceRecord,
    ResolutionStatus,
    SimulationResultSnapshot,
    immutable_mapping,
)
from .providers import (
    ComplexMembership,
    FamilyMembership,
    GeneAnnotation,
    LocalEvidenceAdapters,
    LocalProviderResult,
    LocalProviderStatus,
    PartialUniProtContext,
    RegulatoryRelation,
)
from .semantic_adapter import (
    FunctionalContextKind,
    FunctionalContextRow,
    MetricValue,
    ResponseDirection,
    ResultRole,
    ResultSource,
    SemanticRecord,
    SemanticSimulationResult,
    adapt_snapshot,
)


_SEMANTIC_ANNOTATION_FIELDS: Mapping[str, str] = MappingProxyType({
    "Lokalizasyon": "localization",
    "Gümrük_Kapisi": "gateway",
    "significant_redistribution": "fdr_supported",
    "empirical_p": "empirical_p",
    "q_value": "q_value",
    "Selection_Reason": "fdr_selection_reason",
    "Düzenleyici_TFler": "tf_regulators",
    "Hedef Gen Etkisi": "target_regulatory_effect",
    "GO_CC_Terimleri": "go_cellular_component",
    "GO_MF_Terimleri": "go_molecular_function",
    "GO Biyolojik Süreç": "go_biological_process",
    "GO Moleküler İşlev": "go_molecular_function",
    "GO Hücresel Bileşen": "go_cellular_component",
})


def _freeze(value: Any) -> Any:
    if isinstance(value, Mapping):
        return immutable_mapping(value)
    if isinstance(value, (list, tuple)):
        return tuple(_freeze(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze(item) for item in value)
    return value


def _missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, str):
        return value.strip().casefold() in {"", "nan", "none", "<na>"}
    if isinstance(value, (tuple, list, set, frozenset, Mapping)):
        return len(value) == 0
    return False


@dataclass(frozen=True, slots=True)
class LocalAnnotation:
    """One unchanged local/source value with explicit provenance."""

    annotation_type: str
    value: Any
    source: str
    source_field: str | None = None
    provenance: tuple[ProvenanceRecord, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "annotation_type", str(self.annotation_type).strip())
        object.__setattr__(self, "source", str(self.source).strip())
        object.__setattr__(self, "value", _freeze(self.value))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        if not self.annotation_type or not self.source:
            raise ValueError("annotation_type and source are required")


@dataclass(frozen=True, slots=True)
class LocalSourceStatus:
    source_id: str
    status: LocalProviderStatus
    requested_count: int = 0
    matched_count: int = 0
    missing_sources: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_id", str(self.source_id).strip())
        object.__setattr__(self, "status", LocalProviderStatus(self.status))
        object.__setattr__(self, "missing_sources", tuple(self.missing_sources))
        object.__setattr__(self, "errors", tuple(self.errors))
        if not self.source_id or min(self.requested_count, self.matched_count) < 0:
            raise ValueError("source status requires an ID and non-negative counts")

    @property
    def provider_id(self) -> str:
        return self.source_id


@dataclass(frozen=True, slots=True)
class LocalEntityContext:
    entity: Entity
    resolution: EntityResolution
    roles: tuple[ResultRole, ...]
    directions: tuple[ResponseDirection, ...]
    metrics: tuple[MetricValue, ...]
    annotations: tuple[LocalAnnotation, ...]
    source_records: tuple[SemanticRecord, ...]
    provenance: tuple[ProvenanceRecord, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "roles", tuple(ResultRole(role) for role in self.roles))
        object.__setattr__(self, "directions", tuple(ResponseDirection(item) for item in self.directions))
        object.__setattr__(self, "metrics", tuple(self.metrics))
        object.__setattr__(self, "annotations", tuple(self.annotations))
        object.__setattr__(self, "source_records", tuple(self.source_records))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        if self.resolution.entity != self.entity:
            raise ValueError("entity context must retain its exact resolution entity")

    def annotations_of_type(self, annotation_type: str) -> tuple[LocalAnnotation, ...]:
        return tuple(item for item in self.annotations if item.annotation_type == annotation_type)


@dataclass(frozen=True, slots=True)
class UnresolvedLocalEntity:
    query: str
    resolution: EntityResolution
    roles: tuple[ResultRole, ...]
    source_records: tuple[SemanticRecord, ...]
    annotations: tuple[LocalAnnotation, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "roles", tuple(ResultRole(role) for role in self.roles))
        object.__setattr__(self, "source_records", tuple(self.source_records))
        object.__setattr__(self, "annotations", tuple(self.annotations))
        if self.resolution.status not in {ResolutionStatus.AMBIGUOUS, ResolutionStatus.UNRESOLVED}:
            raise ValueError("unresolved context requires AMBIGUOUS or UNRESOLVED status")


@dataclass(frozen=True, slots=True)
class LocalResearchContext:
    snapshot_id: str
    simulation_id: str
    semantic_result: SemanticSimulationResult
    entities: tuple[LocalEntityContext, ...]
    unresolved_entities: tuple[UnresolvedLocalEntity, ...]
    roles: Mapping[str, tuple[ResultRole, ...]]
    annotations: Mapping[str, tuple[LocalAnnotation, ...]]
    functional_context: tuple[FunctionalContextRow, ...]
    source_statuses: tuple[LocalSourceStatus, ...]
    provenance: tuple[ProvenanceRecord, ...]
    notices: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "entities", tuple(self.entities))
        object.__setattr__(self, "unresolved_entities", tuple(self.unresolved_entities))
        object.__setattr__(self, "roles", immutable_mapping(self.roles))
        object.__setattr__(self, "annotations", immutable_mapping(self.annotations))
        object.__setattr__(self, "functional_context", tuple(self.functional_context))
        object.__setattr__(self, "source_statuses", tuple(self.source_statuses))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        object.__setattr__(self, "notices", tuple(self.notices))
        object.__setattr__(self, "errors", tuple(self.errors))
        if self.snapshot_id != self.semantic_result.snapshot_id:
            raise ValueError("local context and semantic result snapshot IDs must match")
        if self.simulation_id != self.semantic_result.simulation_id:
            raise ValueError("local context and semantic result simulation IDs must match")

    def entity_for(self, logical_key: str) -> LocalEntityContext | None:
        return next((item for item in self.entities if item.entity.logical_key == logical_key), None)

    def status_for(self, source_id: str) -> LocalSourceStatus | None:
        return next((item for item in self.source_statuses if item.source_id == source_id), None)


@dataclass(slots=True)
class _Seed:
    canonical_id: str | None
    symbol: str | None
    roles: list[ResultRole] = field(default_factory=list)
    records: list[SemanticRecord] = field(default_factory=list)


@dataclass(slots=True)
class _Accumulator:
    entity: Entity
    resolution: EntityResolution
    roles: list[ResultRole] = field(default_factory=list)
    directions: list[ResponseDirection] = field(default_factory=list)
    metrics: list[MetricValue] = field(default_factory=list)
    annotations: list[LocalAnnotation] = field(default_factory=list)
    records: list[SemanticRecord] = field(default_factory=list)
    provenance: list[ProvenanceRecord] = field(default_factory=list)


def _append_unique(target: list[Any], values: Iterable[Any], *, key) -> None:
    seen = {key(item) for item in target}
    for value in values:
        identity = key(value)
        if identity not in seen:
            target.append(value)
            seen.add(identity)


def _provenance_key(item: ProvenanceRecord) -> tuple[Any, ...]:
    return (
        item.source,
        item.kind.value,
        item.version,
        item.retrieved_at,
        item.locator,
        item.snapshot_id,
        item.attribution,
        repr(dict(item.details)),
    )


def _annotation_key(item: LocalAnnotation) -> tuple[Any, ...]:
    return item.annotation_type, item.source, item.source_field, repr(item.value)


def _record_key(item: SemanticRecord) -> tuple[Any, ...]:
    return item.source.value, item.source_position, item.canonical_id, item.symbol


def _metric_key(item: MetricValue) -> tuple[Any, ...]:
    return item.metric_id, item.source_field, repr(item.raw_value)


def _snapshot_provenance_for(
    snapshot: SimulationResultSnapshot, source: ResultSource
) -> tuple[ProvenanceRecord, ...]:
    table_name = {
        ResultSource.REPORT: "report",
        ResultSource.SIGNED_REDISTRIBUTION: "signed_redistribution",
        ResultSource.ENRICHMENT: "enrichment",
    }.get(source)
    if table_name is None:
        return tuple(snapshot.provenance)
    matched = tuple(
        record
        for record in snapshot.provenance
        if record.details.get("table") == table_name
    )
    return matched or tuple(snapshot.provenance)


def _annotations_from_records(
    records: Iterable[SemanticRecord], snapshot: SimulationResultSnapshot
) -> tuple[LocalAnnotation, ...]:
    annotations: list[LocalAnnotation] = []
    for record in records:
        provenance = _snapshot_provenance_for(snapshot, record.source)
        for source_field, annotation_type in _SEMANTIC_ANNOTATION_FIELDS.items():
            if source_field not in record.values:
                continue
            if (
                source_field == "significant_redistribution"
                and ResultRole.FDR_SUPPORTED not in record.roles
            ):
                # The stored flag is authoritative only in the semantic
                # adapter's null-FDR role.  A stale top-n value must never be
                # relabelled as FDR-supported evidence here.
                continue
            value = record.values[source_field]
            if _missing(value):
                continue
            annotations.append(LocalAnnotation(
                annotation_type=annotation_type,
                value=value,
                source=f"sophiark_{record.source.value}",
                source_field=source_field,
                provenance=provenance,
            ))
    unique: list[LocalAnnotation] = []
    _append_unique(unique, annotations, key=_annotation_key)
    return tuple(unique)


def _semantic_seeds(view: SemanticSimulationResult) -> tuple[_Seed, ...]:
    seeds: dict[tuple[str | None, str | None], _Seed] = {}
    records = (
        *view.target_records,
        *view.positive_redistribution,
        *view.network_losses,
        *view.fdr_supported,
    )
    for record in records:
        key = (record.canonical_id, record.symbol)
        seed = seeds.setdefault(key, _Seed(record.canonical_id, record.symbol))
        _append_unique(seed.roles, record.roles, key=lambda role: role.value)
        _append_unique(seed.records, (record,), key=_record_key)

    represented_targets = {
        seed.canonical_id for seed in seeds.values()
        if ResultRole.PERTURBATION_TARGET in seed.roles
    }
    for target_id in view.target_ids:
        if target_id in represented_targets:
            continue
        seed = seeds.setdefault((target_id, None), _Seed(target_id, None))
        _append_unique(seed.roles, (ResultRole.PERTURBATION_TARGET,), key=lambda role: role.value)
    return tuple(seeds.values())


def _safe_provider_call(
    source_id: str,
    method: Any,
    values: tuple[str, ...],
    *,
    taxon_id: int,
) -> LocalProviderResult[Any]:
    if not values:
        return LocalProviderResult(
            provider_id=source_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.NO_MATCH,
            requested_count=0,
            message="No resolved entities required this local lookup.",
        )
    try:
        result = method(values, taxon_id=taxon_id)
    except Exception as exc:
        return LocalProviderResult(
            provider_id=source_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.PROVIDER_ERROR,
            requested_count=len(values),
            errors=(f"{type(exc).__name__}: {exc}",),
            message="Local provider failed; semantic Sophiark results were retained.",
        )
    if not isinstance(result, LocalProviderResult):
        return LocalProviderResult(
            provider_id=source_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.PROVIDER_ERROR,
            requested_count=len(values),
            errors=("Invalid local provider result contract.",),
            message="Local provider returned an invalid result; semantic Sophiark results were retained.",
        )
    if result.taxon_id != taxon_id:
        return LocalProviderResult(
            provider_id=source_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.PROVIDER_ERROR,
            requested_count=len(values),
            errors=("Local provider result crossed the requested taxon boundary.",),
            message="Cross-species local evidence was rejected.",
        )
    cross_taxon_records = tuple(
        record
        for record in result.records
        if getattr(record, "taxon_id", taxon_id) != taxon_id
    )
    if cross_taxon_records:
        return LocalProviderResult(
            provider_id=result.provider_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.PROVIDER_ERROR,
            provenance=result.provenance,
            requested_count=result.requested_count,
            missing_sources=result.missing_sources,
            errors=(*result.errors, "Local provider records crossed the requested taxon boundary."),
            message="Cross-species local evidence records were rejected.",
        )
    return result


def _status(result: LocalProviderResult[Any]) -> LocalSourceStatus:
    return LocalSourceStatus(
        source_id=result.provider_id,
        status=result.status,
        requested_count=result.requested_count,
        matched_count=result.matched_count,
        missing_sources=result.missing_sources,
        errors=result.errors,
        message=result.message,
    )


def _semantic_source_statuses(snapshot: SimulationResultSnapshot) -> tuple[LocalSourceStatus, ...]:
    return (
        LocalSourceStatus(
            "sophiark_report", LocalProviderStatus.AVAILABLE,
            requested_count=snapshot.report.row_count,
            matched_count=snapshot.report.row_count,
        ),
        LocalSourceStatus(
            "sophiark_signed_redistribution",
            LocalProviderStatus.AVAILABLE if snapshot.signed_redistribution is not None else LocalProviderStatus.UNAVAILABLE,
            requested_count=snapshot.signed_redistribution.row_count if snapshot.signed_redistribution else 0,
            matched_count=snapshot.signed_redistribution.row_count if snapshot.signed_redistribution else 0,
            message=None if snapshot.signed_redistribution is not None else "Signed redistribution source is unavailable.",
        ),
        LocalSourceStatus(
            "sophiark_enrichment",
            LocalProviderStatus.AVAILABLE if snapshot.enrichment is not None else LocalProviderStatus.UNAVAILABLE,
            requested_count=snapshot.enrichment.row_count if snapshot.enrichment else 0,
            matched_count=snapshot.enrichment.row_count if snapshot.enrichment else 0,
            message=None if snapshot.enrichment is not None else "GO/KEGG enrichment source is unavailable.",
        ),
    )


def build_local_research_context(
    snapshot: SimulationResultSnapshot,
    *,
    resolver: EntityResolver,
    providers: LocalEvidenceAdapters,
    semantic_result: SemanticSimulationResult | None = None,
) -> LocalResearchContext:
    """Build local context without external fallback or scientific mutation."""

    if not isinstance(snapshot, SimulationResultSnapshot):
        raise TypeError("snapshot must be a SimulationResultSnapshot")
    view = semantic_result if semantic_result is not None else adapt_snapshot(snapshot)
    if not isinstance(view, SemanticSimulationResult):
        raise TypeError("semantic_result must be a SemanticSimulationResult")
    if view.snapshot_id != snapshot.snapshot_id or view.simulation_id != snapshot.simulation_id:
        raise ValueError("semantic_result does not belong to the supplied snapshot")

    accumulators: dict[str, _Accumulator] = {}
    unresolved: list[UnresolvedLocalEntity] = []
    resolution_count = 0

    for seed in _semantic_seeds(view):
        query = seed.canonical_id or seed.symbol or "<missing-identifier>"
        try:
            resolution = resolver.resolve(
                query,
                taxon_id=view.scope.taxon_id,
                allow_external=False,
            )
        except Exception as exc:
            resolution = EntityResolution(
                query=query,
                taxon_id=view.scope.taxon_id,
                entity_type="protein",
                status=ResolutionStatus.UNRESOLVED,
                source="entity_resolver",
                message=f"Entity resolution failed safely: {type(exc).__name__}.",
            )
        semantic_annotations = _annotations_from_records(seed.records, snapshot)
        if resolution.entity is None:
            unresolved.append(UnresolvedLocalEntity(
                query=query,
                resolution=resolution,
                roles=tuple(seed.roles),
                source_records=tuple(seed.records),
                annotations=semantic_annotations,
            ))
            continue

        resolution_count += 1
        key = resolution.entity.logical_key
        accumulator = accumulators.get(key)
        if accumulator is None:
            accumulator = _Accumulator(resolution.entity, resolution)
            accumulators[key] = accumulator
        _append_unique(accumulator.roles, seed.roles, key=lambda role: role.value)
        _append_unique(accumulator.records, seed.records, key=_record_key)
        _append_unique(
            accumulator.directions,
            (record.direction for record in seed.records),
            key=lambda direction: direction.value,
        )
        _append_unique(
            accumulator.metrics,
            (metric for record in seed.records for metric in record.metrics),
            key=_metric_key,
        )
        _append_unique(accumulator.annotations, semantic_annotations, key=_annotation_key)
        for record in seed.records:
            _append_unique(
                accumulator.provenance,
                _snapshot_provenance_for(snapshot, record.source),
                key=_provenance_key,
            )
        if not seed.records:
            _append_unique(accumulator.provenance, snapshot.provenance, key=_provenance_key)

    canonical_lookup: dict[str, list[_Accumulator]] = {}
    symbol_lookup: dict[str, list[_Accumulator]] = {}
    for accumulator in accumulators.values():
        canonical_lookup.setdefault(accumulator.entity.canonical_id.casefold(), []).append(accumulator)
        source_symbols = (
            accumulator.entity.symbol,
            *(record.symbol for record in accumulator.records),
        )
        for symbol in source_symbols:
            if symbol:
                symbol_lookup.setdefault(symbol.casefold(), []).append(accumulator)

    canonical_ids = tuple(accumulator.entity.canonical_id for accumulator in accumulators.values())
    symbols = tuple(
        dict.fromkeys(
            symbol
            for accumulator in accumulators.values()
            for symbol in (
                accumulator.entity.symbol,
                *(record.symbol for record in accumulator.records),
            )
            if symbol
        )
    )
    canonical_and_symbols = tuple(dict.fromkeys((*canonical_ids, *symbols)))

    provider_specs = (
        ("local_mygene", getattr(getattr(providers, "mygene", None), "get_annotations", None), canonical_ids),
        ("local_hgnc_family", getattr(getattr(providers, "families", None), "get_families", None), symbols),
        ("local_complex_membership", getattr(getattr(providers, "complexes", None), "get_memberships", None), canonical_and_symbols),
        ("local_uniprot_partial", getattr(getattr(providers, "uniprot", None), "get_context", None), canonical_ids),
        ("local_regulatory", getattr(getattr(providers, "regulatory", None), "get_relations", None), symbols),
    )
    provider_results: list[LocalProviderResult[Any]] = []
    for source_id, method, values in provider_specs:
        if method is None:
            result = LocalProviderResult(
                provider_id=source_id,
                taxon_id=view.scope.taxon_id,
                status=LocalProviderStatus.PROVIDER_ERROR,
                requested_count=len(values),
                errors=("Local provider method is unavailable.",),
                message="Local provider facade is incomplete; Sophiark results were retained.",
            )
        else:
            result = _safe_provider_call(
                source_id, method, values, taxon_id=view.scope.taxon_id,
            )
        provider_results.append(result)

    def matching_accumulators(*identifiers: object) -> tuple[_Accumulator, ...]:
        matched: dict[str, _Accumulator] = {}
        for identifier in identifiers:
            text = str(identifier or "").strip().casefold()
            if not text:
                continue
            for accumulator in (*canonical_lookup.get(text, ()), *symbol_lookup.get(text, ())):
                matched[accumulator.entity.logical_key] = accumulator
        return tuple(matched.values())

    def attach(
        targets: Iterable[_Accumulator], annotation_type: str, value: Any,
        source: str, source_field: str | None, provenance: tuple[ProvenanceRecord, ...],
    ) -> None:
        annotation = LocalAnnotation(
            annotation_type=annotation_type,
            value=value,
            source=source,
            source_field=source_field,
            provenance=provenance,
        )
        for accumulator in targets:
            _append_unique(accumulator.annotations, (annotation,), key=_annotation_key)
            _append_unique(accumulator.provenance, provenance, key=_provenance_key)

    for result in provider_results:
        for record in result.records:
            if isinstance(record, GeneAnnotation):
                targets = matching_accumulators(record.canonical_id, record.symbol)
                for annotation_type, value, field_name in (
                    ("gene_name", record.name, "name"),
                    ("go_biological_process", record.go_bp, "go_bp"),
                    ("go_cellular_component", record.go_cc, "go_cc"),
                    ("go_molecular_function", record.go_mf, "go_mf"),
                ):
                    if not _missing(value):
                        attach(targets, annotation_type, value, result.provider_id, field_name, (record.provenance,))
            elif isinstance(record, FamilyMembership):
                attach(
                    matching_accumulators(record.query_identifier, record.approved_symbol),
                    "authoritative_family_membership", record, result.provider_id,
                    "group_id/group_name", (record.provenance,),
                )
            elif isinstance(record, ComplexMembership):
                attach(
                    matching_accumulators(record.member_identifier),
                    "complex_membership", record, result.provider_id,
                    "complex_id", (record.provenance,),
                )
            elif isinstance(record, PartialUniProtContext):
                attach(
                    matching_accumulators(record.ensembl_protein_id),
                    "partial_uniprot_context", record, result.provider_id,
                    "uniprot_accession", (record.provenance,),
                )
            elif isinstance(record, RegulatoryRelation):
                attach(
                    matching_accumulators(record.source_symbol, record.target_symbol),
                    "directed_regulation", record, result.provider_id,
                    "source_symbol/target_symbol", (record.provenance,),
                )

    enrichment_provenance = _snapshot_provenance_for(snapshot, ResultSource.ENRICHMENT)
    for functional in view.functional_rows:
        annotation_type = {
            FunctionalContextKind.GO: "go_enrichment",
            FunctionalContextKind.KEGG: "kegg_enrichment",
            FunctionalContextKind.UNCLASSIFIED: "unclassified_functional_context",
        }[functional.kind]
        matched: dict[str, _Accumulator] = {}
        for member in functional.members:
            for accumulator in matching_accumulators(member):
                matched[accumulator.entity.logical_key] = accumulator
        attach(
            matched.values(), annotation_type, functional, "sophiark_enrichment",
            "Genes", enrichment_provenance,
        )

    entity_contexts = tuple(
        LocalEntityContext(
            entity=accumulator.entity,
            resolution=accumulator.resolution,
            roles=tuple(accumulator.roles),
            directions=tuple(accumulator.directions),
            metrics=tuple(accumulator.metrics),
            annotations=tuple(accumulator.annotations),
            source_records=tuple(accumulator.records),
            provenance=tuple(accumulator.provenance),
        )
        for accumulator in accumulators.values()
    )
    role_map = {
        item.entity.logical_key: item.roles
        for item in entity_contexts
    }
    annotation_map = {
        item.entity.logical_key: item.annotations
        for item in entity_contexts
    }

    all_provenance: list[ProvenanceRecord] = []
    _append_unique(all_provenance, snapshot.provenance, key=_provenance_key)
    for result in provider_results:
        _append_unique(all_provenance, result.provenance, key=_provenance_key)
        for record in result.records:
            record_provenance = getattr(record, "provenance", None)
            if isinstance(record_provenance, ProvenanceRecord):
                _append_unique(all_provenance, (record_provenance,), key=_provenance_key)

    resolved_status = (
        LocalProviderStatus.AVAILABLE
        if resolution_count and not unresolved
        else LocalProviderStatus.PARTIAL
        if resolution_count
        else LocalProviderStatus.NO_MATCH
    )
    statuses = [*_semantic_source_statuses(snapshot), LocalSourceStatus(
        "entity_resolution",
        resolved_status,
        requested_count=resolution_count + len(unresolved),
        matched_count=resolution_count,
        message=(
            f"{len(unresolved)} entity records remain ambiguity-safe and unresolved."
            if unresolved else None
        ),
    )]
    statuses.extend(_status(result) for result in provider_results)

    errors = tuple(
        f"{result.provider_id}: {error}"
        for result in provider_results
        for error in result.errors
    )
    provider_notices = tuple(
        f"{result.provider_id}: {result.message}"
        for result in provider_results
        if result.message
    )
    return LocalResearchContext(
        snapshot_id=snapshot.snapshot_id,
        simulation_id=snapshot.simulation_id,
        semantic_result=view,
        entities=entity_contexts,
        unresolved_entities=tuple(unresolved),
        roles=role_map,
        annotations=annotation_map,
        functional_context=view.functional_rows,
        source_statuses=tuple(statuses),
        provenance=tuple(all_provenance),
        notices=tuple((*view.notices, *provider_notices)),
        errors=errors,
    )


# Concise orchestration alias.
build_local_context = build_local_research_context


__all__ = [
    "LocalAnnotation",
    "LocalEntityContext",
    "LocalResearchContext",
    "LocalSourceStatus",
    "UnresolvedLocalEntity",
    "build_local_context",
    "build_local_research_context",
]
