"""Phase 6 acceptance tests for the provider-agnostic Query Planner."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
import time

import pandas as pd
import pytest

from src.research.config import ResearchConfig, RetrievalBudget
from src.research.facade import ResearchContextService, build_research_snapshot
from src.research.models import ObservationType
from src.research.observations import DetectedObservation
from src.research.query_planner import (
    CURRENT_RESULT_ORDER_NOTE,
    LookupOperation,
    QueryPlanner,
    RetrievalLevel,
    plan_research_queries,
)
from tests.test_research_service import _adapters, _resolver, _service, _snapshot


def _bundle(*, config: ResearchConfig | None = None):
    bundle = _service().build_offline_bundle(_snapshot())
    return bundle if config is None else replace(bundle, config=config)


def _wide_budget(max_provider_requests: int = 100) -> RetrievalBudget:
    return RetrievalBudget(
        max_entities=100,
        max_families=50,
        max_functional_terms=100,
        max_relationships=250,
        max_provider_requests=max_provider_requests,
        max_literature_queries=20,
        max_publications_per_query=50,
    )


def _bundle_with_duplicate_observations():
    config = ResearchConfig(
        research_context_enabled=True,
        deep_context_budget=_wide_budget(),
    )
    bundle = _bundle(config=config)
    functional = bundle.observation_result.of_type(
        ObservationType.FUNCTIONAL_COOCCURRENCE
    )[0]
    members = bundle.local_context.entities[:2]

    family_base = replace(
        functional.observation,
        observation_id="observation:family:first",
        type=ObservationType.FAMILY_COOCCURRENCE,
        members=tuple(item.entity for item in members),
        basis="First local family observation.",
        source_fields=("family_id",),
    )
    family_first = DetectedObservation(
        family_base,
        {"family_id": "RAMP", "family_name": "RAMP family"},
    )
    family_second = DetectedObservation(
        replace(
            family_base,
            observation_id="observation:family:second",
            basis="Second local family observation.",
            members=tuple(item.entity for item in bundle.local_context.entities[1:3]),
        ),
        {"family_id": "RAMP", "family_name": "RAMP family"},
    )
    functional_second = DetectedObservation(
        replace(
            functional.observation,
            observation_id="observation:function:second",
            basis="Second local functional observation.",
        ),
        functional.details,
        functional.assertions,
    )
    observations = (
        family_first,
        functional,
        family_second,
        functional_second,
        *bundle.observation_result.observations[1:],
    )
    result = replace(bundle.observation_result, observations=observations)
    return replace(bundle, observation_result=result)


def test_default_is_quick_and_expand_research_requires_explicit_level() -> None:
    bundle = _bundle()
    quick = plan_research_queries(
        bundle,
        expected_snapshot_id=bundle.snapshot_id,
    )
    deep = plan_research_queries(
        bundle,
        expected_snapshot_id=bundle.snapshot_id,
        level=RetrievalLevel.EXPAND_RESEARCH,
    )

    assert quick.level is RetrievalLevel.QUICK_CONTEXT
    assert quick.budget is bundle.config.quick_context_budget
    assert deep.level is RetrievalLevel.EXPAND_RESEARCH
    assert deep.budget is bundle.config.deep_context_budget
    assert quick.provider_request_count <= quick.budget.max_provider_requests == 6
    assert deep.provider_request_count >= quick.provider_request_count
    assert {
        LookupOperation.FUNCTION_CONTEXT,
        LookupOperation.LITERATURE_SEARCH,
        LookupOperation.RELATIONSHIP_CONTEXT,
        LookupOperation.ENTITY_CONTEXT,
    } <= {intent.operation for intent in quick.intents}


def test_relationship_entities_preserve_each_current_result_bucket_order() -> None:
    config = ResearchConfig(
        research_context_enabled=True,
        quick_context_budget=_wide_budget(),
    )
    bundle = _bundle(config=config)
    plan = QueryPlanner(config).plan(
        bundle,
        expected_snapshot_id=bundle.snapshot_id,
    )

    assert [entity.canonical_id for entity in plan.selected_relationship_entities] == [
        "ENSP_TARGET",
        "ENSP_GAIN_B",  # Existing presentation order, not magnitude order.
        "ENSP_GAIN_A",
        "ENSP_LOSS_SHALLOW",  # Existing loss order, not magnitude order.
        "ENSP_LOSS_DEEP",
        "ENSP_LOSS_REPORT",  # Added by gateway/FDR overlap without reranking.
    ]
    assert plan.selection_note == CURRENT_RESULT_ORDER_NOTE
    assert all(
        intent.selection_note == CURRENT_RESULT_ORDER_NOTE
        for intent in plan.intents
        if intent.operation in {
            LookupOperation.ENTITY_CONTEXT,
            LookupOperation.RELATIONSHIP_CONTEXT,
        }
    )


def test_duplicate_family_function_and_literature_queries_are_coalesced() -> None:
    bundle = _bundle_with_duplicate_observations()
    plan = plan_research_queries(
        bundle,
        expected_snapshot_id=bundle.snapshot_id,
        level=RetrievalLevel.EXPAND_RESEARCH,
    )

    family = tuple(
        intent for intent in plan.intents_for(LookupOperation.FAMILY_CONTEXT)
        if intent.canonical_query == "family:RAMP"
    )
    function = tuple(
        intent for intent in plan.intents_for(LookupOperation.FUNCTION_CONTEXT)
        if "Fixture pathway" in intent.canonical_query
    )
    family_literature = tuple(
        intent for intent in plan.intents_for(LookupOperation.LITERATURE_SEARCH)
        if intent.canonical_query.startswith("family:RAMP|")
    )
    function_literature = tuple(
        intent for intent in plan.intents_for(LookupOperation.LITERATURE_SEARCH)
        if intent.canonical_query.startswith("function:KEGG:Fixture pathway|")
    )
    assert len(family) == len(function) == 1
    assert len(family_literature) == len(function_literature) == 1
    assert family[0].observation_ids == (
        "observation:family:first",
        "observation:family:second",
    )
    assert len(family[0].entities) == 3
    assert "First local family observation." in family[0].reasons
    assert "Second local family observation." in family[0].reasons
    assert plan.duplicate_intent_count >= 4


def test_category_and_global_budgets_drop_and_report_without_scientific_thresholds() -> None:
    tiny = RetrievalBudget(
        max_entities=2,
        max_families=0,
        max_functional_terms=1,
        max_relationships=1,
        max_provider_requests=2,
        max_literature_queries=0,
        max_publications_per_query=3,
    )
    config = ResearchConfig(research_context_enabled=True, quick_context_budget=tiny)
    bundle = _bundle_with_duplicate_observations()
    bundle = replace(bundle, config=config)
    plan = QueryPlanner(config).plan(bundle, expected_snapshot_id=bundle.snapshot_id)

    assert plan.provider_request_count <= 2
    assert len(plan.selected_relationship_entities) <= 2
    assert plan.dropped_intent_count == plan.coalesced_intent_count - len(plan.intents)
    assert plan.dropped_by_budget[LookupOperation.FAMILY_CONTEXT.value] >= 1
    assert plan.dropped_by_budget[LookupOperation.LITERATURE_SEARCH.value] >= 1
    assert "threshold" not in " ".join(plan.dropped_by_budget).casefold()


def test_snapshot_simulation_taxon_and_query_version_are_strictly_bound() -> None:
    bundle = _bundle()
    planner = QueryPlanner(query_version="fixture-query-v7")
    with pytest.raises(ValueError, match="stale"):
        planner.plan(bundle, expected_snapshot_id="old-snapshot")
    with pytest.raises(ValueError, match="simulation"):
        planner.plan(
            bundle,
            expected_snapshot_id=bundle.snapshot_id,
            expected_simulation_id="different-simulation",
        )

    plan = planner.plan(
        bundle,
        expected_snapshot_id=bundle.snapshot_id,
        expected_simulation_id=bundle.simulation_id,
    )
    assert plan.snapshot_id == bundle.snapshot_id
    assert plan.simulation_id == bundle.simulation_id
    assert plan.taxon_id == 9606
    assert plan.query_version == "fixture-query-v7"
    assert all(intent.snapshot_id == bundle.snapshot_id for intent in plan.intents)
    assert all(intent.taxon_id == 9606 for intent in plan.intents)
    assert all(intent.query_version == "fixture-query-v7" for intent in plan.intents)


def test_snapshot_target_without_report_row_is_still_selected_by_exact_id() -> None:
    report = _snapshot().report.to_frame()
    report = report.loc[report["gene"] != "ENSP_TARGET"].reset_index(drop=True)
    snapshot = _snapshot(report=report)
    bundle = _service().build_offline_bundle(snapshot)
    plan = QueryPlanner().plan(bundle, expected_snapshot_id=bundle.snapshot_id)

    assert bundle.semantic_result.target_records == ()
    assert plan.selected_relationship_entities[0].canonical_id == "ENSP_TARGET"
    assert plan.relationship_entity_reasons[
        plan.selected_relationship_entities[0].logical_key
    ] == ("current target",)


def test_models_are_immutable_and_same_snapshot_plan_is_deterministic() -> None:
    bundle = _bundle()
    planner = QueryPlanner()
    first = planner.plan(bundle, expected_snapshot_id=bundle.snapshot_id)
    second = planner.plan(bundle, expected_snapshot_id=bundle.snapshot_id)

    assert first == second
    assert first.plan_id == second.plan_id
    with pytest.raises(FrozenInstanceError):
        first.plan_id = "changed"  # type: ignore[misc]
    with pytest.raises(TypeError):
        first.dropped_by_budget["changed"] = 1  # type: ignore[index]
    if first.intents:
        with pytest.raises(FrozenInstanceError):
            first.intents[0].canonical_query = "changed"  # type: ignore[misc]


def test_two_hundred_entities_do_not_expand_to_entities_times_providers() -> None:
    count = 200
    identifiers = [(f"ENSP_{index:04d}", f"G{index:04d}") for index in range(count)]
    report = pd.DataFrame({
        "gene": [item[0] for item in identifiers],
        "Symbol": [item[1] for item in identifiers],
        "Hasar_Tipi": ["Birincil_Hasar_Hedef", *("Üçüncül_Hasar_Stresli",) * (count - 1)],
        "Response_Direction": [None, *("Influence Gain",) * (count - 1)],
        "Delta_PageRank_Pct": [None, *(float(index) for index in range(1, count))],
    })
    snapshot = build_research_snapshot(
        report,
        None,
        None,
        species="Homo sapiens",
        taxon_id=9606,
        tissue="Pancreas",
        targets=(identifiers[0][0],),
        attenuation=0.01,
        selection_mode="top_n",
        test_limit=None,
        tested_count=200,
        returned_count=199,
        top_n=200,
        created_at="2026-01-01T00:00:00Z",
    )
    service = ResearchContextService(
        _resolver(identifiers),
        _adapters(),
        ResearchConfig(research_context_enabled=True),
    )
    bundle = service.build_offline_bundle(snapshot)
    started = time.perf_counter()
    plan = QueryPlanner().plan(bundle, expected_snapshot_id=bundle.snapshot_id)
    elapsed = time.perf_counter() - started

    assert elapsed < 1.0
    assert plan.provider_request_count <= 6
    assert len(plan.selected_relationship_entities) <= 20
    assert [item.canonical_id for item in plan.selected_relationship_entities[:6]] == [
        "ENSP_0000", "ENSP_0001", "ENSP_0002", "ENSP_0003", "ENSP_0004", "ENSP_0005",
    ]
    assert all(intent.operation in LookupOperation for intent in plan.intents)


def test_planner_source_has_no_provider_io_ui_session_or_scientific_dependency() -> None:
    source_path = Path(__file__).parents[1] / "src" / "research" / "query_planner.py"
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)
    imports = {
        alias.name
        for node in ast.walk(tree)
        if isinstance(node, ast.Import)
        for alias in node.names
    }
    imports.update(
        str(node.module or "")
        for node in ast.walk(tree)
        if isinstance(node, ast.ImportFrom)
    )
    forbidden = (
        "requests", "httpx", "urllib", "sqlite3", "streamlit",
        "src.state", "src.scientific", "src.graph_engine",
        "src.research.provider_runtime", "src.research.provider_cache",
    )
    assert not any(name.startswith(forbidden) for name in imports)
    for text in ("session_state", "urlopen(", "read_csv(", "pagerank(", "provider.execute("):
        assert text not in source
