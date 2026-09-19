"""Phase 7 local-first routing and cache-first integration tests."""

from __future__ import annotations

from dataclasses import replace

from src.research.config import ProviderTimeouts, ResearchConfig, RetrievalBudget
from src.research.external_biology import ExternalBiologyLoader, plan_external_biology
from src.research.http_transport import JsonHttpResponse
from src.research.models import Entity, EntityType, ProvenanceKind, ProvenanceRecord
from src.research.provider_cache import CacheFirstProviderClient, ProviderCache
from src.research.provider_runtime import ProviderExecutionPolicy, ProviderExecutor, ProviderStatus
from src.research.providers.external_biology import UniProtProvider
from src.research.providers.local import PartialUniProtContext
from src.research.query_planner import (
    CURRENT_RESULT_ORDER_NOTE,
    LookupIntent,
    LookupOperation,
    QueryPlan,
    RetrievalLevel,
)


class CountingTransport:
    def __init__(self, data) -> None:
        self.data = data
        self.calls = 0

    def get_json(self, url, *, params, timeout, headers=None):
        self.calls += 1
        return JsonHttpResponse(200, self.data, {"x-uniprot-release": "fixture"}, url)


def _budget(max_provider_requests: int = 6) -> RetrievalBudget:
    return RetrievalBudget(20, 10, 20, 40, max_provider_requests, 3, 10)


def _entity() -> Entity:
    return Entity(9606, EntityType.PROTEIN, "ENSP00000362353", symbol="GLP1R")


def _intent(operation: LookupOperation, query: str, entities=()) -> LookupIntent:
    return LookupIntent(
        intent_id=f"intent:{operation.value}:{query}",
        operation=operation,
        taxon_id=9606,
        canonical_query=query,
        query_version="research-context-query-v1",
        snapshot_id="snapshot-1",
        entities=tuple(entities),
    )


def _plan(*intents: LookupIntent, max_provider_requests: int = 6) -> QueryPlan:
    return QueryPlan(
        plan_id="plan-1",
        snapshot_id="snapshot-1",
        simulation_id="simulation-1",
        taxon_id=9606,
        query_version="research-context-query-v1",
        level=RetrievalLevel.QUICK_CONTEXT,
        budget=_budget(max_provider_requests),
        intents=tuple(intents),
        selected_relationship_entities=(_entity(),),
        relationship_entity_reasons={_entity().logical_key: ("current target",)},
        selection_note=CURRENT_RESULT_ORDER_NOTE,
        raw_intent_count=len(intents),
        coalesced_intent_count=len(intents),
        duplicate_intent_count=0,
        dropped_intent_count=0,
    )


def _local(*, complete: bool = False) -> PartialUniProtContext:
    fields = (
        ("protein_name", "function", "subcellular_location", "cross_references")
        if complete
        else ("uniprot_accession", "structure.mean_plddt")
    )
    return PartialUniProtContext(
        taxon_id=9606,
        ensembl_protein_id="ENSP00000362353",
        uniprot_accession="P43220",
        structure={"mean_plddt": 91.0},
        ligand_summary={},
        available_fields=fields,
        unavailable_fields=() if complete else ("protein_name", "function", "subcellular_location", "cross_references"),
        provenance=ProvenanceRecord("local fixture", ProvenanceKind.LOCAL_ANNOTATION),
    )


def _enabled_config(**changes) -> ResearchConfig:
    values = dict(
        research_context_enabled=True,
        external_context_enabled=True,
        provider_uniprot_enabled=True,
        provider_interpro_enabled=True,
        provider_quickgo_enabled=True,
        provider_reactome_enabled=True,
    )
    values.update(changes)
    return ResearchConfig(**values)


def test_external_plan_requires_explicit_action_and_disabled_default_never_routes() -> None:
    query_plan = _plan(_intent(LookupOperation.ENTITY_CONTEXT, _entity().logical_key, (_entity(),)))
    implicit = plan_external_biology(query_plan, local_uniprot_records=(_local(),), config=_enabled_config(), explicitly_requested=False)
    disabled = plan_external_biology(query_plan, local_uniprot_records=(_local(),), config=ResearchConfig(), explicitly_requested=True)
    assert implicit.requests == ()
    assert implicit.skipped == {"not_explicitly_requested": 1}
    assert disabled.requests == ()
    assert disabled.skipped == {"external_context_disabled": 1}


def test_local_exact_mapping_drives_bounded_provider_requests_without_symbol_guessing() -> None:
    entity_intent = _intent(LookupOperation.ENTITY_CONTEXT, _entity().logical_key, (_entity(),))
    function_intent = _intent(LookupOperation.FUNCTION_CONTEXT, "function:go:GO:0005886", (_entity(),))
    plan = _plan(function_intent, entity_intent, max_provider_requests=3)

    external = plan_external_biology(plan, local_uniprot_records=(_local(),), config=_enabled_config(), explicitly_requested=True)

    assert len(external.requests) == 3
    assert [request.provider for request in external.requests] == ["quickgo", "uniprot", "interpro"]
    assert external.requests[0].canonical_query == "GO:0005886"
    assert all(request.taxon_id == 9606 for request in external.requests)
    assert external.exact_accessions == ("P43220",)
    assert external.local_records_used == ("ENSP00000362353",)
    assert external.skipped["provider_budget"] == 2
    assert all("GLP1R" not in request.canonical_query for request in external.requests)


def test_complete_local_fields_suppress_uniprot_but_not_other_exact_context() -> None:
    query_plan = _plan(_intent(LookupOperation.ENTITY_CONTEXT, _entity().logical_key, (_entity(),)))
    external = plan_external_biology(query_plan, local_uniprot_records=(_local(complete=True),), config=_enabled_config(), explicitly_requested=True)
    assert "uniprot" not in [request.provider for request in external.requests]
    assert [request.provider for request in external.requests] == ["interpro", "quickgo", "reactome"]
    assert external.skipped["uniprot_satisfied_locally"] == 1


def test_symbol_only_entity_is_not_converted_to_an_accession() -> None:
    symbol_entity = Entity(9606, EntityType.PROTEIN, "GLP1R", symbol="GLP1R")
    query_plan = _plan(_intent(LookupOperation.ENTITY_CONTEXT, symbol_entity.logical_key, (symbol_entity,)))
    external = plan_external_biology(query_plan, config=_enabled_config(), explicitly_requested=True)
    assert external.requests == ()
    assert external.skipped["no_exact_uniprot_accession"] == 1


def test_loader_is_cache_first_and_second_load_uses_zero_transport_calls() -> None:
    query_plan = _plan(_intent(LookupOperation.ENTITY_CONTEXT, _entity().logical_key, (_entity(),)), max_provider_requests=1)
    config = _enabled_config(
        provider_interpro_enabled=False,
        provider_quickgo_enabled=False,
        provider_reactome_enabled=False,
    )
    external = plan_external_biology(query_plan, local_uniprot_records=(_local(),), config=config, explicitly_requested=True)
    payload = {
        "primaryAccession": "P43220",
        "organism": {"taxonId": 9606, "scientificName": "Homo sapiens"},
        "proteinDescription": {"recommendedName": {"fullName": {"value": "GLP1R"}}},
    }
    transport = CountingTransport(payload)
    provider = UniProtProvider(transport=transport, enabled=True)
    executor = ProviderExecutor(provider, policy=ProviderExecutionPolicy(max_retries=0))
    client = CacheFirstProviderClient(cache=ProviderCache(), executor=executor)
    loader = ExternalBiologyLoader({"uniprot": client})
    try:
        first = loader.load(external, timeouts=ProviderTimeouts())
        second = loader.load(external, timeouts=ProviderTimeouts())
    finally:
        executor.close()
    assert first.results[0].status is ProviderStatus.AVAILABLE
    assert second.results[0].status is ProviderStatus.AVAILABLE
    assert transport.calls == 1


def test_loader_missing_client_is_fail_soft_and_preserves_request_identity() -> None:
    query_plan = _plan(_intent(LookupOperation.ENTITY_CONTEXT, _entity().logical_key, (_entity(),)), max_provider_requests=1)
    config = _enabled_config(provider_interpro_enabled=False, provider_quickgo_enabled=False, provider_reactome_enabled=False)
    external = plan_external_biology(query_plan, local_uniprot_records=(_local(),), config=config, explicitly_requested=True)
    result = ExternalBiologyLoader({}).load(external, timeouts=ProviderTimeouts()).results[0]
    assert result.status is ProviderStatus.DISABLED
    assert result.attempts == 0
    assert result.request.cache_key == external.requests[0].cache_key

