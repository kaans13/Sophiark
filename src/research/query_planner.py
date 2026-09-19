"""Deterministic, provider-agnostic lookup planning for Research Context.

The planner consumes one already-built :class:`OfflineResearchBundle`.  It
does not call providers, caches, files, UI state, graphs, or scientific code.
It only translates local observations and preserved result order into a small,
auditable set of lookup intents constrained by centralized retrieval budgets.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass, field, replace
from enum import Enum
import hashlib
import json
from types import MappingProxyType
from typing import Any

from .config import ResearchConfig, RetrievalBudget
from .models import Entity, ObservationType, immutable_mapping
from .observations import DetectedObservation
from .semantic_adapter import ResultRole, SemanticRecord
from .service import OfflineResearchBundle


QUERY_PLANNER_VERSION = "research-context-query-v1"
CURRENT_RESULT_ORDER_NOTE = "Selected from current result order."


class RetrievalLevel(str, Enum):
    QUICK_CONTEXT = "QUICK_CONTEXT"
    EXPAND_RESEARCH = "EXPAND_RESEARCH"


class LookupOperation(str, Enum):
    ENTITY_CONTEXT = "entity_context"
    FAMILY_CONTEXT = "family_context"
    FUNCTION_CONTEXT = "function_context"
    RELATIONSHIP_CONTEXT = "relationship_context"
    LITERATURE_SEARCH = "literature_search"


def _unique_text(values: Iterable[Any]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(
        str(value).strip()
        for value in values
        if value is not None and str(value).strip()
    ))


def _unique_entities(values: Iterable[Entity]) -> tuple[Entity, ...]:
    unique: dict[str, Entity] = {}
    for entity in values:
        if not isinstance(entity, Entity):
            raise TypeError("lookup intent entities must be Entity values")
        unique.setdefault(entity.logical_key, entity)
    return tuple(unique.values())


def _query_text(value: Any) -> str:
    text = " ".join(str(value).strip().split())
    if not text:
        raise ValueError("canonical_query cannot be empty")
    return text


def _stable_json(value: Any) -> str:
    def plain(item: Any) -> Any:
        if isinstance(item, Mapping):
            return {
                str(key): plain(child)
                for key, child in sorted(item.items(), key=lambda pair: str(pair[0]))
            }
        if isinstance(item, (tuple, list)):
            return [plain(child) for child in item]
        if isinstance(item, Enum):
            return item.value
        if isinstance(item, Entity):
            return item.logical_key
        if isinstance(item, (str, int, float, bool)) or item is None:
            return item
        return str(item)

    return json.dumps(
        plain(value), ensure_ascii=False, sort_keys=True,
        separators=(",", ":"), allow_nan=False,
    )


@dataclass(frozen=True, slots=True)
class LookupIntent:
    """One semantic lookup request, not yet assigned to a provider."""

    intent_id: str
    operation: LookupOperation
    taxon_id: int
    canonical_query: str
    query_version: str
    snapshot_id: str
    entities: tuple[Entity, ...] = ()
    observation_ids: tuple[str, ...] = ()
    reasons: tuple[str, ...] = ()
    selection_note: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "operation", LookupOperation(self.operation))
        object.__setattr__(self, "canonical_query", _query_text(self.canonical_query))
        object.__setattr__(self, "query_version", _query_text(self.query_version))
        object.__setattr__(self, "snapshot_id", _query_text(self.snapshot_id))
        object.__setattr__(self, "entities", _unique_entities(self.entities))
        object.__setattr__(self, "observation_ids", _unique_text(self.observation_ids))
        object.__setattr__(self, "reasons", _unique_text(self.reasons))
        object.__setattr__(self, "metadata", immutable_mapping(self.metadata))
        if not str(self.intent_id).strip():
            raise ValueError("intent_id cannot be empty")
        if isinstance(self.taxon_id, bool) or not isinstance(self.taxon_id, int):
            raise TypeError("taxon_id must be an integer")
        if self.taxon_id <= 0:
            raise ValueError("taxon_id must be positive")
        if any(entity.taxon_id != self.taxon_id for entity in self.entities):
            raise ValueError("lookup intent entities must match the intent taxon")
        if self.selection_note is not None:
            object.__setattr__(self, "selection_note", str(self.selection_note).strip() or None)

    @property
    def coalescing_key(self) -> tuple[str, int, str]:
        return (
            self.operation.value,
            self.taxon_id,
            self.canonical_query.casefold(),
        )


@dataclass(frozen=True, slots=True)
class QueryPlan:
    """One immutable, snapshot-bound and budget-bounded retrieval plan."""

    plan_id: str
    snapshot_id: str
    simulation_id: str
    taxon_id: int
    query_version: str
    level: RetrievalLevel
    budget: RetrievalBudget
    intents: tuple[LookupIntent, ...]
    selected_relationship_entities: tuple[Entity, ...]
    relationship_entity_reasons: Mapping[str, tuple[str, ...]]
    selection_note: str
    raw_intent_count: int
    coalesced_intent_count: int
    duplicate_intent_count: int
    dropped_intent_count: int
    dropped_by_budget: Mapping[str, int] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "level", RetrievalLevel(self.level))
        object.__setattr__(self, "intents", tuple(self.intents))
        object.__setattr__(
            self,
            "selected_relationship_entities",
            _unique_entities(self.selected_relationship_entities),
        )
        object.__setattr__(
            self,
            "relationship_entity_reasons",
            immutable_mapping(self.relationship_entity_reasons),
        )
        object.__setattr__(self, "dropped_by_budget", immutable_mapping(self.dropped_by_budget))
        if not self.plan_id or not self.snapshot_id or not self.simulation_id or not self.query_version:
            raise ValueError("query plan identity cannot be empty")
        if not isinstance(self.budget, RetrievalBudget):
            raise TypeError("budget must be a RetrievalBudget")
        if self.selection_note != CURRENT_RESULT_ORDER_NOTE:
            raise ValueError("relationship selection must disclose current result order")
        counts = (
            self.raw_intent_count,
            self.coalesced_intent_count,
            self.duplicate_intent_count,
            self.dropped_intent_count,
        )
        if any(isinstance(value, bool) or not isinstance(value, int) or value < 0 for value in counts):
            raise ValueError("query plan counts must be non-negative integers")
        if self.raw_intent_count - self.coalesced_intent_count != self.duplicate_intent_count:
            raise ValueError("duplicate intent count is inconsistent")
        if self.coalesced_intent_count - len(self.intents) != self.dropped_intent_count:
            raise ValueError("dropped intent count is inconsistent")
        if len(self.intents) > self.budget.max_provider_requests:
            raise ValueError("query plan exceeds max_provider_requests")
        if len(self.selected_relationship_entities) > self.budget.max_entities:
            raise ValueError("relationship selection exceeds max_entities")
        if any(entity.taxon_id != self.taxon_id for entity in self.selected_relationship_entities):
            raise ValueError("selected relationship entities must match plan taxon")
        for intent in self.intents:
            if (
                intent.taxon_id != self.taxon_id
                or intent.snapshot_id != self.snapshot_id
                or intent.query_version != self.query_version
            ):
                raise ValueError("query plan contains a stale or cross-taxon intent")

    @property
    def provider_request_count(self) -> int:
        return len(self.intents)

    def intents_for(self, operation: LookupOperation | str) -> tuple[LookupIntent, ...]:
        selected = LookupOperation(operation)
        return tuple(intent for intent in self.intents if intent.operation is selected)


def _intent(
    *,
    operation: LookupOperation,
    taxon_id: int,
    canonical_query: str,
    query_version: str,
    snapshot_id: str,
    entities: Iterable[Entity] = (),
    observations: Iterable[DetectedObservation] = (),
    reasons: Iterable[str] = (),
    selection_note: str | None = None,
    metadata: Mapping[str, Any] | None = None,
) -> LookupIntent:
    entity_tuple = _unique_entities(entities)
    observation_tuple = tuple(observations)
    observation_ids = _unique_text(item.observation_id for item in observation_tuple)
    reason_tuple = _unique_text((*reasons, *(item.basis for item in observation_tuple)))
    normalized_query = _query_text(canonical_query)
    identity = {
        "operation": operation.value,
        "taxon_id": taxon_id,
        "canonical_query": normalized_query.casefold(),
        "query_version": query_version,
        "snapshot_id": snapshot_id,
        "entities": tuple(item.logical_key for item in entity_tuple),
        "observation_ids": observation_ids,
        "reasons": reason_tuple,
    }
    digest = hashlib.sha256(_stable_json(identity).encode("utf-8")).hexdigest()[:24]
    return LookupIntent(
        intent_id=f"lookup:{operation.value}:{digest}",
        operation=operation,
        taxon_id=taxon_id,
        canonical_query=normalized_query,
        query_version=query_version,
        snapshot_id=snapshot_id,
        entities=entity_tuple,
        observation_ids=observation_ids,
        reasons=reason_tuple,
        selection_note=selection_note,
        metadata=metadata or {},
    )


def _entity_indexes(bundle: OfflineResearchBundle) -> tuple[dict[str, Entity], dict[str, tuple[ResultRole, ...]]]:
    entities: dict[str, Entity] = {}
    roles: dict[str, tuple[ResultRole, ...]] = {}
    for context in bundle.local_context.entities:
        values: list[str | None] = [context.entity.canonical_id, context.entity.symbol]
        for record in context.source_records:
            values.extend((record.canonical_id, record.symbol))
        for value in values:
            if value is not None and str(value).strip():
                entities.setdefault(str(value).strip().casefold(), context.entity)
        roles[context.entity.logical_key] = context.roles
    return entities, roles


def _entities_for_records(
    records: Sequence[SemanticRecord], index: Mapping[str, Entity],
) -> tuple[Entity, ...]:
    selected: list[Entity] = []
    for record in records:
        entity = None
        for value in (record.canonical_id, record.symbol):
            if value is not None:
                entity = index.get(str(value).strip().casefold())
            if entity is not None:
                break
        if entity is not None:
            selected.append(entity)
    return _unique_entities(selected)


def _relationship_selection(
    bundle: OfflineResearchBundle,
    budget: RetrievalBudget,
) -> tuple[tuple[Entity, ...], Mapping[str, tuple[str, ...]]]:
    if budget.max_entities == 0:
        return (), MappingProxyType({})

    entity_index, _ = _entity_indexes(bundle)
    semantic = bundle.semantic_result
    per_bucket = max(1, budget.max_entities // 4)
    target_entities = _unique_entities((
        *_entities_for_records(semantic.target_records, entity_index),
        *(
            entity_index[target_id.strip().casefold()]
            for target_id in semantic.target_ids
            if target_id.strip().casefold() in entity_index
        ),
    ))
    buckets: tuple[tuple[str, tuple[Entity, ...]], ...] = (
        ("current target", target_entities),
        (
            "positive response from current result order",
            _entities_for_records(semantic.positive_redistribution, entity_index)[:per_bucket],
        ),
        (
            "network loss from current result order",
            _entities_for_records(semantic.network_losses, entity_index)[:per_bucket],
        ),
        (
            "gateway overlap",
            _unique_entities(
                member
                for observation in bundle.observation_result.of_type(
                    ObservationType.GATEWAY_RESPONSE_OVERLAP
                )
                for member in observation.members
            )[:per_bucket],
        ),
        (
            "stored FDR-supported selection",
            _entities_for_records(semantic.fdr_supported, entity_index)[:per_bucket],
        ),
    )

    selected: dict[str, Entity] = {}
    reasons: dict[str, list[str]] = {}
    for reason, entities in buckets:
        for entity in entities:
            reasons.setdefault(entity.logical_key, []).append(reason)
            if entity.logical_key in selected:
                continue
            if len(selected) >= budget.max_entities:
                continue
            selected[entity.logical_key] = entity
    return (
        tuple(selected.values()),
        MappingProxyType({key: _unique_text(value) for key, value in reasons.items() if key in selected}),
    )


def _pair_entities(
    observation: DetectedObservation,
    roles: Mapping[str, tuple[ResultRole, ...]],
) -> tuple[Entity, Entity] | None:
    if len(observation.members) != 2:
        return None
    target = next(
        (
            member for member in observation.members
            if ResultRole.PERTURBATION_TARGET in roles.get(member.logical_key, ())
        ),
        None,
    )
    if target is None:
        return None
    response = next(member for member in observation.members if member != target)
    return target, response


def _raw_intents(
    bundle: OfflineResearchBundle,
    *,
    selected_entities: tuple[Entity, ...],
    selection_reasons: Mapping[str, tuple[str, ...]],
    query_version: str,
    budget: RetrievalBudget,
) -> tuple[LookupIntent, ...]:
    taxon_id = bundle.snapshot.taxon_id
    snapshot_id = bundle.snapshot_id
    tissue = _query_text(bundle.snapshot.tissue)
    target_scope = ",".join(bundle.snapshot.targets)
    target_clause = f"|targets:{target_scope}" if target_scope else ""
    observations = bundle.observation_result.observations
    _, role_index = _entity_indexes(bundle)
    intents: list[LookupIntent] = []

    for observation in observations:
        if observation.type is ObservationType.FAMILY_COOCCURRENCE:
            family_id = observation.details.get("family_id") or observation.details.get("group_id")
            family_name = observation.details.get("family_name") or family_id
            if family_id:
                family_query = f"family:{_query_text(family_id)}"
                intents.append(_intent(
                    operation=LookupOperation.FAMILY_CONTEXT,
                    taxon_id=taxon_id,
                    canonical_query=family_query,
                    query_version=query_version,
                    snapshot_id=snapshot_id,
                    entities=observation.members,
                    observations=(observation,),
                ))
                intents.append(_intent(
                    operation=LookupOperation.LITERATURE_SEARCH,
                    taxon_id=taxon_id,
                    canonical_query=f"{family_query}{target_clause}|tissue:{tissue}",
                    query_version=query_version,
                    snapshot_id=snapshot_id,
                    entities=observation.members,
                    observations=(observation,),
                    metadata={
                        "family_name": _query_text(family_name),
                        "max_publications": budget.max_publications_per_query,
                    },
                ))

        elif observation.type is ObservationType.FUNCTIONAL_COOCCURRENCE:
            kind = observation.details.get("context_kind")
            term = observation.details.get("term")
            if kind and term:
                function_query = f"function:{_query_text(kind)}:{_query_text(term)}"
                intents.append(_intent(
                    operation=LookupOperation.FUNCTION_CONTEXT,
                    taxon_id=taxon_id,
                    canonical_query=function_query,
                    query_version=query_version,
                    snapshot_id=snapshot_id,
                    entities=observation.members,
                    observations=(observation,),
                ))
                intents.append(_intent(
                    operation=LookupOperation.LITERATURE_SEARCH,
                    taxon_id=taxon_id,
                    canonical_query=f"{function_query}{target_clause}|tissue:{tissue}",
                    query_version=query_version,
                    snapshot_id=snapshot_id,
                    entities=observation.members,
                    observations=(observation,),
                    metadata={"max_publications": budget.max_publications_per_query},
                ))

        elif observation.type is ObservationType.TARGET_RESPONSE_RELATIONSHIP:
            pair = _pair_entities(observation, role_index)
            if pair is not None:
                target, response = pair
                intents.append(_intent(
                    operation=LookupOperation.RELATIONSHIP_CONTEXT,
                    taxon_id=taxon_id,
                    canonical_query=f"{target.logical_key}->{response.logical_key}",
                    query_version=query_version,
                    snapshot_id=snapshot_id,
                    entities=pair,
                    observations=(observation,),
                    selection_note=CURRENT_RESULT_ORDER_NOTE,
                ))

    target_entities = tuple(
        entity for entity in selected_entities
        if ResultRole.PERTURBATION_TARGET in role_index.get(entity.logical_key, ())
    )
    response_entities = tuple(entity for entity in selected_entities if entity not in target_entities)
    observations_by_member: dict[str, list[DetectedObservation]] = {}
    for observation in observations:
        for member in observation.members:
            observations_by_member.setdefault(member.logical_key, []).append(observation)

    # After direct observation intents, retain a small amount of entity context
    # before expanding all selected target/response pairs.  This preserves local
    # observation priority while preventing relationship fan-out from consuming
    # the entire QUICK_CONTEXT provider budget.
    for entity in selected_entities:
        intents.append(_intent(
            operation=LookupOperation.ENTITY_CONTEXT,
            taxon_id=taxon_id,
            canonical_query=entity.logical_key,
            query_version=query_version,
            snapshot_id=snapshot_id,
            entities=(entity,),
            observations=observations_by_member.get(entity.logical_key, ()),
            reasons=selection_reasons.get(entity.logical_key, ()),
            selection_note=CURRENT_RESULT_ORDER_NOTE,
        ))

    for target in target_entities:
        for response in response_entities:
            related = _unique_observations((
                *observations_by_member.get(target.logical_key, ()),
                *observations_by_member.get(response.logical_key, ()),
            ))
            intents.append(_intent(
                operation=LookupOperation.RELATIONSHIP_CONTEXT,
                taxon_id=taxon_id,
                canonical_query=f"{target.logical_key}->{response.logical_key}",
                query_version=query_version,
                snapshot_id=snapshot_id,
                entities=(target, response),
                observations=related,
                reasons=selection_reasons.get(response.logical_key, ()),
                selection_note=CURRENT_RESULT_ORDER_NOTE,
            ))
    return tuple(intents)


def _unique_observations(values: Iterable[DetectedObservation]) -> tuple[DetectedObservation, ...]:
    unique: dict[str, DetectedObservation] = {}
    for observation in values:
        unique.setdefault(observation.observation_id, observation)
    return tuple(unique.values())


def _coalesce(intents: Iterable[LookupIntent]) -> tuple[LookupIntent, ...]:
    merged: dict[tuple[str, int, str], LookupIntent] = {}
    for intent in intents:
        key = intent.coalescing_key
        existing = merged.get(key)
        if existing is None:
            merged[key] = intent
            continue
        combined_entities = _unique_entities((*existing.entities, *intent.entities))
        combined_observations = _unique_text((*existing.observation_ids, *intent.observation_ids))
        combined_reasons = _unique_text((*existing.reasons, *intent.reasons))
        metadata = dict(existing.metadata)
        metadata.update({key: value for key, value in intent.metadata.items() if key not in metadata})
        combined = replace(
            existing,
            entities=combined_entities,
            observation_ids=combined_observations,
            reasons=combined_reasons,
            selection_note=existing.selection_note or intent.selection_note,
            metadata=metadata,
        )
        identity = {
            "operation": combined.operation.value,
            "taxon_id": combined.taxon_id,
            "canonical_query": combined.canonical_query.casefold(),
            "query_version": combined.query_version,
            "snapshot_id": combined.snapshot_id,
            "entities": tuple(item.logical_key for item in combined.entities),
            "observation_ids": combined.observation_ids,
            "reasons": combined.reasons,
        }
        digest = hashlib.sha256(_stable_json(identity).encode("utf-8")).hexdigest()[:24]
        merged[key] = replace(combined, intent_id=f"lookup:{combined.operation.value}:{digest}")
    return tuple(merged.values())


def _apply_budget(
    intents: tuple[LookupIntent, ...], budget: RetrievalBudget,
) -> tuple[tuple[LookupIntent, ...], Mapping[str, int]]:
    limits = {
        LookupOperation.ENTITY_CONTEXT: budget.max_entities,
        LookupOperation.FAMILY_CONTEXT: budget.max_families,
        LookupOperation.FUNCTION_CONTEXT: budget.max_functional_terms,
        LookupOperation.RELATIONSHIP_CONTEXT: budget.max_relationships,
        LookupOperation.LITERATURE_SEARCH: budget.max_literature_queries,
    }
    seen: Counter[LookupOperation] = Counter()
    dropped: Counter[str] = Counter()
    category_bounded: list[LookupIntent] = []
    for intent in intents:
        if seen[intent.operation] >= limits[intent.operation]:
            dropped[intent.operation.value] += 1
            continue
        seen[intent.operation] += 1
        category_bounded.append(intent)

    allowed = budget.max_provider_requests
    if len(category_bounded) > allowed:
        dropped["max_provider_requests"] += len(category_bounded) - allowed
    return tuple(category_bounded[:allowed]), MappingProxyType(dict(dropped))


class QueryPlanner:
    """Translate local observations to bounded lookup intents without I/O."""

    def __init__(
        self,
        config: ResearchConfig | None = None,
        *,
        query_version: str = QUERY_PLANNER_VERSION,
    ) -> None:
        if config is not None and not isinstance(config, ResearchConfig):
            raise TypeError("config must be a ResearchConfig or None")
        self._config = config
        self._query_version = _query_text(query_version)

    @property
    def query_version(self) -> str:
        return self._query_version

    def plan(
        self,
        bundle: OfflineResearchBundle,
        *,
        expected_snapshot_id: str,
        expected_simulation_id: str | None = None,
        level: RetrievalLevel | str = RetrievalLevel.QUICK_CONTEXT,
    ) -> QueryPlan:
        if not isinstance(bundle, OfflineResearchBundle):
            raise TypeError("bundle must be an OfflineResearchBundle")
        if str(expected_snapshot_id).strip() != bundle.snapshot_id:
            raise ValueError("stale Research Context snapshot; query planning rejected")
        if (
            expected_simulation_id is not None
            and str(expected_simulation_id).strip() != bundle.simulation_id
        ):
            raise ValueError("query plan simulation does not match the current result")
        if bundle.snapshot.taxon_id != bundle.semantic_result.scope.taxon_id:
            raise ValueError("query plan taxon does not match the semantic result")

        selected_level = RetrievalLevel(level)
        config = self._config or bundle.config
        budget = (
            config.quick_context_budget
            if selected_level is RetrievalLevel.QUICK_CONTEXT
            else config.deep_context_budget
        )
        selected_entities, selection_reasons = _relationship_selection(bundle, budget)
        raw = _raw_intents(
            bundle,
            selected_entities=selected_entities,
            selection_reasons=selection_reasons,
            query_version=self.query_version,
            budget=budget,
        )
        coalesced = _coalesce(raw)
        bounded, dropped = _apply_budget(coalesced, budget)
        plan_identity = {
            "snapshot_id": bundle.snapshot_id,
            "simulation_id": bundle.simulation_id,
            "taxon_id": bundle.snapshot.taxon_id,
            "query_version": self.query_version,
            "level": selected_level.value,
            "intent_ids": tuple(intent.intent_id for intent in bounded),
            "selected_entities": tuple(entity.logical_key for entity in selected_entities),
        }
        digest = hashlib.sha256(_stable_json(plan_identity).encode("utf-8")).hexdigest()[:24]
        return QueryPlan(
            plan_id=f"query-plan:{digest}",
            snapshot_id=bundle.snapshot_id,
            simulation_id=bundle.simulation_id,
            taxon_id=bundle.snapshot.taxon_id,
            query_version=self.query_version,
            level=selected_level,
            budget=budget,
            intents=bounded,
            selected_relationship_entities=selected_entities,
            relationship_entity_reasons=selection_reasons,
            selection_note=CURRENT_RESULT_ORDER_NOTE,
            raw_intent_count=len(raw),
            coalesced_intent_count=len(coalesced),
            duplicate_intent_count=len(raw) - len(coalesced),
            dropped_intent_count=len(coalesced) - len(bounded),
            dropped_by_budget=dropped,
        )


def plan_research_queries(
    bundle: OfflineResearchBundle,
    *,
    expected_snapshot_id: str,
    expected_simulation_id: str | None = None,
    level: RetrievalLevel | str = RetrievalLevel.QUICK_CONTEXT,
    config: ResearchConfig | None = None,
    query_version: str = QUERY_PLANNER_VERSION,
) -> QueryPlan:
    """Convenience wrapper; QUICK_CONTEXT remains the explicit default."""

    return QueryPlanner(config, query_version=query_version).plan(
        bundle,
        expected_snapshot_id=expected_snapshot_id,
        expected_simulation_id=expected_simulation_id,
        level=level,
    )


__all__ = [
    "CURRENT_RESULT_ORDER_NOTE",
    "LookupIntent",
    "LookupOperation",
    "QUERY_PLANNER_VERSION",
    "QueryPlan",
    "QueryPlanner",
    "RetrievalLevel",
    "plan_research_queries",
]
