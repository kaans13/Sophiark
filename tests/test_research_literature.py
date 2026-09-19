"""Phase 8 query, cache/search, passage, and loader integration tests."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory

from src.research.config import ResearchConfig
from src.research.evidence_store import EvidenceStore
from src.research.external_biology import plan_external_biology
from src.research.http_transport import JsonHttpResponse
from src.research.literature import (
    LITERATURE_QUERY_VERSION,
    NO_MATCHING_PUBLICATIONS_MESSAGE,
    LiteratureLoader,
    LiteratureQueryBuilder,
    LiteratureRepository,
    Publication,
    extract_relevant_passages,
    plan_literature_requests,
)
from src.research.models import ProvenanceKind, ProvenanceRecord
from src.research.provider_cache import CacheFirstProviderClient, ProviderCache
from src.research.provider_runtime import ProviderExecutionPolicy, ProviderExecutor, ProviderStatus
from src.research.providers.literature import EuropePMCProvider
from src.research.providers.local import PartialUniProtContext
from src.research.query_planner import LookupOperation, QueryPlanner
from tests.test_research_query_planner import _bundle


def _query_plan():
    bundle = _bundle()
    return QueryPlanner().plan(bundle, expected_snapshot_id=bundle.snapshot_id)


def _config(**changes) -> ResearchConfig:
    values = dict(
        research_context_enabled=True,
        external_context_enabled=True,
        provider_europepmc_enabled=True,
    )
    values.update(changes)
    return ResearchConfig(**values)


def _publication(identifier: str = "MED:1", *, title: str = "GLP1R and RAMP1 in pancreas") -> Publication:
    return Publication(
        publication_id=identifier,
        pmid=identifier.split(":")[-1],
        title=title,
        year=2024,
        journal="Fixture Journal",
        authors_summary="A. Author et al.",
        abstract=(
            "GLP1R and RAMP1 were measured in pancreatic cells. "
            "An unrelated sentence discusses sample preparation. "
            "RAMP1 localization was also reported in the pancreas."
        ),
        source="Europe PMC",
        retrieved_at="2026-08-26T00:00:00+00:00",
        provenance={"source": "Europe PMC", "source_identifier": identifier},
    )


def test_query_builder_is_deterministic_contextual_bounded_and_preserves_reason() -> None:
    plan = _query_plan()
    first, dropped_first = LiteratureQueryBuilder().build(plan)
    second, dropped_second = LiteratureQueryBuilder().build(plan)
    assert first == second
    assert dropped_first == dropped_second
    assert len(first) <= plan.budget.max_literature_queries
    assert first
    query = first[0]
    assert query.query_version == LITERATURE_QUERY_VERSION
    assert "Pancreas" in query.query_text
    assert " AND " in query.query_text
    assert query.matched_entities
    assert query.reason_shown.startswith("Shown because the current result")
    assert query.max_publications <= plan.budget.max_publications_per_query
    assert not any(word in query.reason_shown.casefold() for word in ("validated", "proven", "causal"))
    # Two-character fixture symbols are deliberately excluded as ambiguous;
    # the pathway+tissue anchors keep the query contextual and non-broad.
    assert '"GA"' not in query.query_text and '"GB"' not in query.query_text


def test_query_builder_drops_single_anchor_broad_query() -> None:
    plan = _query_plan()
    literature_intent = next(
        intent for intent in plan.intents
        if intent.operation is LookupOperation.LITERATURE_SEARCH
    )
    broad = replace(
        literature_intent,
        canonical_query="family:RAMP",
        entities=(),
        metadata={},
    )
    narrowed_plan = replace(
        plan,
        intents=(broad,),
        raw_intent_count=1,
        coalesced_intent_count=1,
        duplicate_intent_count=0,
        dropped_intent_count=0,
    )
    queries, dropped = LiteratureQueryBuilder().build(narrowed_plan)
    assert queries == ()
    assert dropped == {"too_broad": 1}


def test_literature_requires_explicit_enablement_and_respects_combined_provider_cap() -> None:
    plan = _query_plan()
    implicit = plan_literature_requests(plan, config=_config(), explicitly_requested=False)
    disabled = plan_literature_requests(plan, config=ResearchConfig(), explicitly_requested=True)
    assert implicit.requests == ()
    assert implicit.dropped == {"not_explicitly_requested": 1}
    assert disabled.requests == ()
    assert disabled.dropped == {"literature_provider_disabled": 1}

    enabled = plan_literature_requests(
        plan,
        config=_config(),
        explicitly_requested=True,
        prior_provider_request_count=plan.budget.max_provider_requests,
    )
    assert enabled.requests == ()
    assert enabled.dropped["combined_provider_budget"] >= 1


def test_external_biology_reserves_literature_budget_and_combined_count_never_exceeds_cap() -> None:
    plan = _query_plan()
    provenance = ProvenanceRecord("fixture", ProvenanceKind.LOCAL_ANNOTATION)
    records = tuple(
        PartialUniProtContext(
            9606,
            entity.canonical_id,
            f"P{index:05d}",
            {},
            {},
            ("uniprot_accession",),
            ("protein_name", "function", "subcellular_location", "cross_references"),
            provenance,
        )
        for index, entity in enumerate(plan.selected_relationship_entities)
    )
    config = _config(
        provider_uniprot_enabled=True,
        provider_interpro_enabled=True,
        provider_quickgo_enabled=True,
        provider_reactome_enabled=True,
    )
    biology = plan_external_biology(
        plan,
        local_uniprot_records=records,
        config=config,
        explicitly_requested=True,
    )
    literature = plan_literature_requests(
        plan,
        config=config,
        explicitly_requested=True,
        prior_provider_request_count=len(biology.requests),
    )
    assert biology.skipped["reserved_for_literature"] == 1
    assert len(biology.requests) + len(literature.requests) <= plan.budget.max_provider_requests


def test_relevant_passages_are_verbatim_ranked_sentences_with_explained_overlap() -> None:
    publication = _publication()
    passages = extract_relevant_passages(
        publication,
        query_text='"GLP1R" AND "RAMP1" AND "pancreas"',
        entity_terms=("GLP1R", "RAMP1"),
        limit=2,
    )
    assert len(passages) == 2
    assert passages[0].rank == 1
    assert passages[0].text in (publication.title, *publication.abstract.split(". ")) or passages[0].text.rstrip(".") in publication.abstract
    assert passages[0].matched_entities
    assert "overlap" in passages[0].basis
    assert all(len(item.text) <= 600 for item in passages)


def test_publication_repository_fts_and_tfidf_fallback_are_idempotent_and_rank_only() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        repository = LiteratureRepository(store)
        assert repository.initialize_search() in {"fts5", "tfidf_fallback"}
        first = repository.store_publications((_publication(), _publication("MED:2", title="Unrelated baseline")))
        second = repository.store_publications((_publication(),))
        assert first.stored_count == 2
        assert second.stored_count == 1
        result = repository.search("GLP1R pancreas", limit=10)
        assert result.method in {"fts5_bm25", "tfidf_fallback"}
        assert result.hits
        assert result.hits[0].publication.publication_id == "MED:1"
        assert result.hits[0].rank == 1
        assert not hasattr(result.hits[0], "score")
        with store.read_connection() as connection:
            assert connection.execute("SELECT COUNT(*) FROM publications").fetchone()[0] == 2

    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        repository = LiteratureRepository(store)
        assert repository.store_publications((_publication(),)).search_backend == "tfidf_fallback"
        result = repository.search("GLP1R pancreas")
        assert result.method == "tfidf_fallback"
        assert result.hits[0].publication.publication_id == "MED:1"


def test_repository_empty_search_uses_required_non_novelty_message() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        result = LiteratureRepository(store).search("no-such-term")
        assert result.hits == ()
        assert result.message == NO_MATCHING_PUBLICATIONS_MESSAGE
        assert "novel" not in result.message.casefold()


def test_uninitialized_repository_is_fail_soft_and_does_not_create_database() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "not-created.sqlite"
        repository = LiteratureRepository(EvidenceStore(path))
        stored = repository.store_publications((_publication(),))
        searched = repository.search("GLP1R")
        assert stored.search_backend == "unavailable"
        assert stored.errors
        assert searched.method == "unavailable"
        assert path.exists() is False


class CountingTransport:
    def __init__(self) -> None:
        self.calls = 0

    def get_json(self, url, *, params, timeout, headers=None):
        self.calls += 1
        return JsonHttpResponse(
            200,
            {
                "version": "fixture",
                "hitCount": 1,
                "resultList": {"result": [{
                    "id": "123", "source": "MED", "pmid": "123",
                    "title": "GB and GA in a Fixture pathway",
                    "pubYear": "2024",
                    "abstractText": "GB and GA were co-mentioned in pancreatic tissue.",
                }]},
            },
            {},
            url,
        )


def test_literature_loader_is_cache_first_stores_normalized_records_and_never_summarizes() -> None:
    plan = _query_plan()
    literature_plan = plan_literature_requests(plan, config=_config(), explicitly_requested=True)
    transport = CountingTransport()
    provider = EuropePMCProvider(transport=transport, enabled=True)
    executor = ProviderExecutor(provider, policy=ProviderExecutionPolicy(max_retries=0))
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        repository = LiteratureRepository(store)
        client = CacheFirstProviderClient(cache=ProviderCache(store=store), executor=executor)
        loader = LiteratureLoader(client, repository=repository)
        try:
            first = loader.load(literature_plan, timeout=_config().provider_timeouts.europepmc)
            second = loader.load(literature_plan, timeout=_config().provider_timeouts.europepmc)
        finally:
            executor.close()
        assert first.results[0].status is ProviderStatus.AVAILABLE
        assert first.matches[0].match_label == "Co-mentioned"
        assert second.matches[0].publication.publication_id == "MED:123"
        assert transport.calls == 1
        assert first.store_result.stored_count == 1
        assert "summary" not in first.matches[0].__dataclass_fields__
        with store.read_connection() as connection:
            assert connection.execute("SELECT COUNT(*) FROM publications").fetchone()[0] == 1
