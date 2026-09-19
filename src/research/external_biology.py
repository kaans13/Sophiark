"""Local-first planning and cache-first loading for external biology context.

This module is deliberately downstream of :mod:`src.research.query_planner`.
It never mutates a simulation snapshot and cannot run unless the caller marks
the load as an explicit Research Context action and enables each provider.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
import re
from types import MappingProxyType
from typing import Iterable, Mapping

from .config import ProviderTimeouts, ResearchConfig
from .models import Entity, EntityType
from .provider_cache import CacheEvidenceClass, CacheFirstProviderClient
from .provider_runtime import ProviderRequest, ProviderResult, ProviderStatus
from .providers.external_biology import (
    INTERPRO_SCHEMA_VERSION,
    QUICKGO_SCHEMA_VERSION,
    REACTOME_SCHEMA_VERSION,
    UNIPROT_SCHEMA_VERSION,
    InterProProvider,
    QuickGOProvider,
    ReactomeProvider,
    UniProtProvider,
)
from .providers.local import PartialUniProtContext
from .query_planner import LookupIntent, LookupOperation, QueryPlan


EXTERNAL_BIOLOGY_ROUTER_VERSION = "external-biology-router-v1"
_ACCESSION = re.compile(r"^[A-Z0-9][A-Z0-9]{5,9}(?:-[1-9][0-9]*)?$")
_GO_ID = re.compile(r"GO:[0-9]{7}", re.IGNORECASE)
_UNIPROT_FIELDS = frozenset({
    "protein_name",
    "function",
    "subcellular_location",
    "cross_references",
})


@dataclass(frozen=True, slots=True)
class ExternalBiologyPlan:
    """A bounded set of exact provider requests derived from one QueryPlan."""

    query_plan_id: str
    snapshot_id: str
    taxon_id: int
    explicit_request: bool
    max_provider_requests: int
    requests: tuple[ProviderRequest, ...] = ()
    exact_accessions: tuple[str, ...] = ()
    local_records_used: tuple[str, ...] = ()
    skipped: Mapping[str, int] = field(default_factory=lambda: MappingProxyType({}))
    router_version: str = EXTERNAL_BIOLOGY_ROUTER_VERSION

    def __post_init__(self) -> None:
        object.__setattr__(self, "requests", tuple(self.requests))
        object.__setattr__(self, "exact_accessions", tuple(self.exact_accessions))
        object.__setattr__(self, "local_records_used", tuple(self.local_records_used))
        object.__setattr__(
            self,
            "skipped",
            MappingProxyType({str(key): int(value) for key, value in self.skipped.items()}),
        )
        if len(self.requests) > self.max_provider_requests:
            raise ValueError("external biology plan exceeds the query plan provider budget")
        if any(request.taxon_id != self.taxon_id for request in self.requests):
            raise ValueError("external biology requests must preserve the query plan taxon")


@dataclass(frozen=True, slots=True)
class ExternalBiologyLoad:
    """Results retain the exact plan order and never replace scientific data."""

    plan: ExternalBiologyPlan
    results: tuple[ProviderResult, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "results", tuple(self.results))
        if len(self.results) != len(self.plan.requests):
            raise ValueError("one provider result is required for every planned request")
        for request, result in zip(self.plan.requests, self.results):
            if request.cache_key != result.request.cache_key:
                raise ValueError("external biology result identity does not match its plan")


def _local_index(
    records: Iterable[PartialUniProtContext],
    *,
    taxon_id: int,
) -> dict[str, PartialUniProtContext]:
    index: dict[str, PartialUniProtContext] = {}
    for record in records:
        if not isinstance(record, PartialUniProtContext):
            raise TypeError("local_uniprot_records must contain PartialUniProtContext values")
        if record.taxon_id != taxon_id:
            continue
        index.setdefault(record.ensembl_protein_id.casefold(), record)
        index.setdefault(record.uniprot_accession.casefold(), record)
    return index


def _exact_accession(
    entity: Entity,
    local_index: Mapping[str, PartialUniProtContext],
) -> tuple[str | None, PartialUniProtContext | None]:
    local = local_index.get(entity.canonical_id.casefold())
    if local is not None:
        return local.uniprot_accession.upper(), local
    canonical = entity.canonical_id.upper()
    if entity.entity_type is EntityType.PROTEIN and _ACCESSION.fullmatch(canonical):
        return canonical, None
    return None, None


def _request(
    *,
    provider: str,
    operation: str,
    taxon_id: int,
    query: str,
    schema_version: str,
    query_version: str,
    params: Mapping[str, object] | None = None,
) -> ProviderRequest:
    return ProviderRequest(
        provider=provider,
        operation=operation,
        taxon_id=taxon_id,
        canonical_query=query,
        provider_schema_version=schema_version,
        query_version=query_version,
        params=params or {},
    )


def _enabled(config: ResearchConfig, provider: str) -> bool:
    return bool(config.research_context_enabled and config.provider_enabled(provider))


def _entity_requests(
    intent: LookupIntent,
    entity: Entity,
    *,
    local_index: Mapping[str, PartialUniProtContext],
    config: ResearchConfig,
    skipped: Counter[str],
) -> tuple[ProviderRequest, ...]:
    accession, local = _exact_accession(entity, local_index)
    if accession is None:
        skipped["no_exact_uniprot_accession"] += 1
        return ()
    requests: list[ProviderRequest] = []
    local_fields = frozenset(local.available_fields) if local is not None else frozenset()
    if _enabled(config, "uniprot"):
        if _UNIPROT_FIELDS.issubset(local_fields):
            skipped["uniprot_satisfied_locally"] += 1
        else:
            requests.append(_request(
                provider="uniprot",
                operation=UniProtProvider.operation,
                taxon_id=intent.taxon_id,
                query=accession,
                schema_version=UNIPROT_SCHEMA_VERSION,
                query_version=intent.query_version,
                params={"fields": sorted(_UNIPROT_FIELDS - local_fields)},
            ))
    if _enabled(config, "interpro"):
        requests.append(_request(
            provider="interpro",
            operation=InterProProvider.operation,
            taxon_id=intent.taxon_id,
            query=accession,
            schema_version=INTERPRO_SCHEMA_VERSION,
            query_version=intent.query_version,
            params={"page_size": 20},
        ))
    if _enabled(config, "quickgo"):
        requests.append(_request(
            provider="quickgo",
            operation=QuickGOProvider.annotation_operation,
            taxon_id=intent.taxon_id,
            query=accession,
            schema_version=QUICKGO_SCHEMA_VERSION,
            query_version=intent.query_version,
            params={"limit": 25},
        ))
    if _enabled(config, "reactome"):
        requests.append(_request(
            provider="reactome",
            operation=ReactomeProvider.operation,
            taxon_id=intent.taxon_id,
            query=accession,
            schema_version=REACTOME_SCHEMA_VERSION,
            query_version=intent.query_version,
            params={"limit": 50},
        ))
    return tuple(requests)


def plan_external_biology(
    query_plan: QueryPlan,
    *,
    local_uniprot_records: Iterable[PartialUniProtContext] = (),
    config: ResearchConfig,
    explicitly_requested: bool,
) -> ExternalBiologyPlan:
    """Translate provider-agnostic intents without symbol or species inference."""

    if not isinstance(query_plan, QueryPlan):
        raise TypeError("query_plan must be QueryPlan")
    if not isinstance(config, ResearchConfig):
        raise TypeError("config must be ResearchConfig")
    if not explicitly_requested:
        return ExternalBiologyPlan(
            query_plan_id=query_plan.plan_id,
            snapshot_id=query_plan.snapshot_id,
            taxon_id=query_plan.taxon_id,
            explicit_request=False,
            max_provider_requests=query_plan.budget.max_provider_requests,
            skipped={"not_explicitly_requested": 1},
        )
    if not config.research_context_enabled or not config.external_context_enabled:
        return ExternalBiologyPlan(
            query_plan_id=query_plan.plan_id,
            snapshot_id=query_plan.snapshot_id,
            taxon_id=query_plan.taxon_id,
            explicit_request=True,
            max_provider_requests=query_plan.budget.max_provider_requests,
            skipped={"external_context_disabled": 1},
        )

    local_index = _local_index(local_uniprot_records, taxon_id=query_plan.taxon_id)
    skipped: Counter[str] = Counter()
    candidates: list[ProviderRequest] = []
    accessions: list[str] = []
    used_local: list[str] = []

    for intent in query_plan.intents:
        if intent.operation is LookupOperation.FUNCTION_CONTEXT and _enabled(config, "quickgo"):
            go_match = _GO_ID.search(intent.canonical_query)
            if go_match:
                candidates.append(_request(
                    provider="quickgo",
                    operation=QuickGOProvider.term_operation,
                    taxon_id=intent.taxon_id,
                    query=go_match.group(0).upper(),
                    schema_version=QUICKGO_SCHEMA_VERSION,
                    query_version=intent.query_version,
                ))
            else:
                skipped["function_not_exact_go_id"] += 1

        if intent.operation not in {LookupOperation.ENTITY_CONTEXT, LookupOperation.FAMILY_CONTEXT}:
            continue
        for entity in intent.entities:
            accession, local = _exact_accession(entity, local_index)
            if accession is not None and accession not in accessions:
                accessions.append(accession)
            if local is not None and local.ensembl_protein_id not in used_local:
                used_local.append(local.ensembl_protein_id)
            if intent.operation is LookupOperation.FAMILY_CONTEXT:
                if accession is None:
                    skipped["family_no_exact_uniprot_accession"] += 1
                elif _enabled(config, "interpro"):
                    candidates.append(_request(
                        provider="interpro",
                        operation=InterProProvider.operation,
                        taxon_id=intent.taxon_id,
                        query=accession,
                        schema_version=INTERPRO_SCHEMA_VERSION,
                        query_version=intent.query_version,
                        params={"page_size": 20},
                    ))
                continue
            candidates.extend(_entity_requests(
                intent,
                entity,
                local_index=local_index,
                config=config,
                skipped=skipped,
            ))

    unique: dict[str, ProviderRequest] = {}
    for request in candidates:
        if request.cache_key in unique:
            skipped["duplicate_request"] += 1
            continue
        unique[request.cache_key] = request
    requests = tuple(unique.values())
    maximum = query_plan.budget.max_provider_requests
    reserved_literature = 0
    if _enabled(config, "europepmc"):
        reserved_literature = min(
            query_plan.budget.max_literature_queries,
            sum(
                intent.operation is LookupOperation.LITERATURE_SEARCH
                for intent in query_plan.intents
            ),
            maximum,
        )
        if reserved_literature:
            skipped["reserved_for_literature"] = reserved_literature
    available_maximum = max(0, maximum - reserved_literature)
    if len(requests) > available_maximum:
        skipped["provider_budget"] += len(requests) - available_maximum
        requests = requests[:available_maximum]
    return ExternalBiologyPlan(
        query_plan_id=query_plan.plan_id,
        snapshot_id=query_plan.snapshot_id,
        taxon_id=query_plan.taxon_id,
        explicit_request=True,
        max_provider_requests=maximum,
        requests=requests,
        exact_accessions=tuple(accessions),
        local_records_used=tuple(used_local),
        skipped=skipped,
    )


class ExternalBiologyLoader:
    """Execute only through injected cache-first clients."""

    def __init__(self, clients: Mapping[str, CacheFirstProviderClient]) -> None:
        normalized: dict[str, CacheFirstProviderClient] = {}
        for name, client in clients.items():
            if not isinstance(client, CacheFirstProviderClient):
                raise TypeError("external biology clients must be CacheFirstProviderClient values")
            normalized[str(name).strip().casefold()] = client
        self._clients = MappingProxyType(normalized)

    def load(
        self,
        plan: ExternalBiologyPlan,
        *,
        timeouts: ProviderTimeouts,
    ) -> ExternalBiologyLoad:
        if not isinstance(plan, ExternalBiologyPlan):
            raise TypeError("plan must be ExternalBiologyPlan")
        if not isinstance(timeouts, ProviderTimeouts):
            raise TypeError("timeouts must be ProviderTimeouts")
        results: list[ProviderResult] = []
        for request in plan.requests:
            client = self._clients.get(request.provider)
            if client is None:
                results.append(ProviderResult(
                    request,
                    ProviderStatus.DISABLED,
                    attempts=0,
                    message="No cache-first client was injected for this provider",
                    metadata={"reason": "client_not_configured"},
                ))
                continue
            timeout = getattr(timeouts, request.provider)
            evidence_class = (
                CacheEvidenceClass.PATHWAY_ANNOTATION
                if request.provider == "reactome"
                else CacheEvidenceClass.BIOLOGICAL_ANNOTATION
            )
            results.append(client.execute(request, timeout, evidence_class=evidence_class))
        return ExternalBiologyLoad(plan=plan, results=tuple(results))


__all__ = [
    "EXTERNAL_BIOLOGY_ROUTER_VERSION",
    "ExternalBiologyLoad",
    "ExternalBiologyLoader",
    "ExternalBiologyPlan",
    "plan_external_biology",
]
