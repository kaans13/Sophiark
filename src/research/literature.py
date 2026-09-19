"""Deterministic literature queries, normalized cache, search, and passages.

The layer is explanatory only.  Query matches are labelled as co-mentions and
never as validation, proof, mechanistic support, contradiction, or novelty.
No LLM, embedding model, vector database, PDF, or full-text download is used.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime, timezone
import json
import math
import re
import sqlite3
from types import MappingProxyType
from typing import Any, Iterable, Mapping, Sequence

from .config import ProviderTimeout, ResearchConfig
from .evidence_store import EvidenceStore, EvidenceStoreError
from .models import Entity, immutable_mapping
from .provider_cache import CacheEvidenceClass, CacheFirstProviderClient
from .provider_runtime import ProviderRequest, ProviderResult, ProviderStatus
from .providers.literature import EUROPEPMC_SCHEMA_VERSION, EuropePMCProvider
from .query_planner import LookupIntent, LookupOperation, QueryPlan


LITERATURE_QUERY_VERSION = "literature-query-v1"
NO_MATCHING_PUBLICATIONS_MESSAGE = (
    "No matching publications were retrieved under the current queries."
)

_WORD = re.compile(r"[A-Za-z0-9]+(?:[-_:][A-Za-z0-9]+)*")
_SENTENCE = re.compile(r"(?<=[.!?])\s+(?=[A-Z0-9])")
_GO = re.compile(r"GO:[0-9]{7}", re.IGNORECASE)
_STOPWORDS = frozenset({
    "and", "or", "the", "a", "an", "of", "in", "to", "for", "with",
    "on", "by", "from", "is", "are", "was", "were", "current", "result",
})


def _text(value: Any, field_name: str) -> str:
    text = " ".join(str(value).strip().split())
    if not text:
        raise ValueError(f"{field_name} cannot be empty")
    return text


def _optional_text(value: Any) -> str | None:
    if value is None:
        return None
    text = " ".join(str(value).strip().split())
    return text or None


def _unique(values: Iterable[Any]) -> tuple[str, ...]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = _optional_text(value)
        if text is None or text.casefold() in seen:
            continue
        seen.add(text.casefold())
        result.append(text)
    return tuple(result)


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _quote(term: str) -> str:
    return f'"{term.replace(chr(34), " ")}"'


@dataclass(frozen=True, slots=True)
class LiteratureQuery:
    query_text: str
    query_version: str
    matched_entities: tuple[str, ...]
    reason_shown: str
    source_intent_id: str
    max_publications: int

    def __post_init__(self) -> None:
        object.__setattr__(self, "query_text", _text(self.query_text, "query_text"))
        object.__setattr__(self, "query_version", _text(self.query_version, "query_version"))
        object.__setattr__(self, "matched_entities", _unique(self.matched_entities))
        object.__setattr__(self, "reason_shown", _text(self.reason_shown, "reason_shown"))
        object.__setattr__(self, "source_intent_id", _text(self.source_intent_id, "source_intent_id"))
        if isinstance(self.max_publications, bool) or not isinstance(self.max_publications, int):
            raise TypeError("max_publications must be an integer")
        if not 1 <= self.max_publications <= 100:
            raise ValueError("max_publications must be between 1 and 100")

    def to_dict(self) -> dict[str, Any]:
        return {
            "query_text": self.query_text,
            "query_version": self.query_version,
            "matched_entities": list(self.matched_entities),
            "reason_shown": self.reason_shown,
            "source_intent_id": self.source_intent_id,
            "max_publications": self.max_publications,
        }


@dataclass(frozen=True, slots=True)
class LiteraturePlan:
    query_plan_id: str
    snapshot_id: str
    taxon_id: int
    explicit_request: bool
    prior_provider_request_count: int
    queries: tuple[LiteratureQuery, ...] = ()
    requests: tuple[ProviderRequest, ...] = ()
    dropped: Mapping[str, int] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "queries", tuple(self.queries))
        object.__setattr__(self, "requests", tuple(self.requests))
        object.__setattr__(self, "dropped", immutable_mapping(self.dropped))
        if (
            isinstance(self.prior_provider_request_count, bool)
            or not isinstance(self.prior_provider_request_count, int)
            or self.prior_provider_request_count < 0
        ):
            raise ValueError("prior_provider_request_count must be a non-negative integer")
        if self.taxon_id <= 0:
            raise ValueError("taxon_id must be positive")
        if len(self.queries) != len(self.requests):
            raise ValueError("every literature query requires exactly one provider request")
        for query, request in zip(self.queries, self.requests):
            if request.canonical_query != query.query_text or request.taxon_id != self.taxon_id:
                raise ValueError("literature request identity does not match its query")


def _segments(canonical_query: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for segment in canonical_query.split("|"):
        if ":" not in segment:
            continue
        key, value = segment.split(":", 1)
        if key.strip() and value.strip():
            parsed[key.strip().casefold()] = value.strip()
    return parsed


def _entity_terms(entities: Sequence[Entity]) -> tuple[str, ...]:
    selected: list[str] = []
    for entity in entities:
        symbol = _optional_text(entity.symbol)
        display_name = _optional_text(entity.display_name)
        canonical_id = entity.canonical_id
        if symbol is not None and len(symbol) >= 3:
            selected.append(symbol)
        elif display_name is not None and len(display_name) >= 3:
            selected.append(display_name)
        elif re.fullmatch(r"ENS[A-Z]*[GPT][0-9]{6,}", canonical_id, re.IGNORECASE):
            selected.append(canonical_id)
        elif re.fullmatch(r"[A-Z0-9][A-Z0-9]{5,9}(?:-[1-9][0-9]*)?", canonical_id, re.IGNORECASE):
            selected.append(canonical_id)
    return _unique(selected)[:6]


def _query_from_intent(intent: LookupIntent, maximum: int) -> LiteratureQuery | None:
    if intent.operation is not LookupOperation.LITERATURE_SEARCH:
        return None
    segments = _segments(intent.canonical_query)
    entities = _entity_terms(intent.entities)
    context: list[str] = []
    family = segments.get("family")
    function = segments.get("function")
    tissue = segments.get("tissue")
    family_name = _optional_text(intent.metadata.get("family_name"))
    if family_name:
        context.append(family_name)
    elif family:
        context.append(family)
    if function:
        go_match = _GO.search(function)
        context.append(go_match.group(0).upper() if go_match else function)
    if tissue:
        context.append(tissue)
    context = list(_unique(context))

    clauses: list[str] = []
    if entities:
        quoted = [_quote(term) for term in entities[:4]]
        clauses.append(quoted[0] if len(quoted) == 1 else f"({' OR '.join(quoted)})")
    clauses.extend(_quote(term) for term in context[:2])
    # One unqualified term is too broad for automatic retrieval.
    if len(clauses) < 2:
        return None
    reason = (
        "Shown because the current result contains a locally detected family pattern."
        if family
        else "Shown because the current result contains a locally detected functional pattern."
    )
    requested_maximum = intent.metadata.get("max_publications", maximum)
    if isinstance(requested_maximum, bool) or not isinstance(requested_maximum, int):
        requested_maximum = maximum
    return LiteratureQuery(
        query_text=" AND ".join(clauses),
        query_version=LITERATURE_QUERY_VERSION,
        matched_entities=tuple(entity.logical_key for entity in intent.entities),
        reason_shown=reason,
        source_intent_id=intent.intent_id,
        max_publications=min(max(requested_maximum, 1), maximum, 100),
    )


class LiteratureQueryBuilder:
    """Convert bounded planner intents to deduplicated, non-broad searches."""

    def build(self, query_plan: QueryPlan) -> tuple[tuple[LiteratureQuery, ...], Mapping[str, int]]:
        if not isinstance(query_plan, QueryPlan):
            raise TypeError("query_plan must be QueryPlan")
        dropped: Counter[str] = Counter()
        unique: dict[str, LiteratureQuery] = {}
        maximum = min(query_plan.budget.max_publications_per_query, 100)
        if maximum <= 0:
            return (), MappingProxyType({"publication_budget": 1})
        for intent in query_plan.intents:
            if intent.operation is not LookupOperation.LITERATURE_SEARCH:
                continue
            query = _query_from_intent(intent, maximum)
            if query is None:
                dropped["too_broad"] += 1
                continue
            key = query.query_text.casefold()
            if key in unique:
                dropped["duplicate_query"] += 1
                continue
            unique[key] = query
        queries = tuple(unique.values())
        limit = query_plan.budget.max_literature_queries
        if len(queries) > limit:
            dropped["literature_query_budget"] += len(queries) - limit
            queries = queries[:limit]
        return queries, MappingProxyType(dict(dropped))


def plan_literature_requests(
    query_plan: QueryPlan,
    *,
    config: ResearchConfig,
    explicitly_requested: bool,
    prior_provider_request_count: int = 0,
    builder: LiteratureQueryBuilder | None = None,
) -> LiteraturePlan:
    """Build Europe PMC requests while preserving the combined provider cap."""

    if not isinstance(query_plan, QueryPlan):
        raise TypeError("query_plan must be QueryPlan")
    if not isinstance(config, ResearchConfig):
        raise TypeError("config must be ResearchConfig")
    if isinstance(prior_provider_request_count, bool) or not isinstance(prior_provider_request_count, int):
        raise TypeError("prior_provider_request_count must be an integer")
    if prior_provider_request_count < 0:
        raise ValueError("prior_provider_request_count cannot be negative")
    base = dict(
        query_plan_id=query_plan.plan_id,
        snapshot_id=query_plan.snapshot_id,
        taxon_id=query_plan.taxon_id,
        prior_provider_request_count=prior_provider_request_count,
    )
    if not explicitly_requested:
        return LiteraturePlan(**base, explicit_request=False, dropped={"not_explicitly_requested": 1})
    if not (
        config.research_context_enabled
        and config.external_context_enabled
        and config.provider_enabled("europepmc")
    ):
        return LiteraturePlan(**base, explicit_request=True, dropped={"literature_provider_disabled": 1})
    queries, dropped_mapping = (builder or LiteratureQueryBuilder()).build(query_plan)
    dropped = Counter(dropped_mapping)
    remaining = max(0, query_plan.budget.max_provider_requests - prior_provider_request_count)
    if len(queries) > remaining:
        dropped["combined_provider_budget"] += len(queries) - remaining
        queries = queries[:remaining]
    requests = tuple(
        ProviderRequest(
            provider="europepmc",
            operation=EuropePMCProvider.operation,
            taxon_id=query_plan.taxon_id,
            canonical_query=query.query_text,
            provider_schema_version=EUROPEPMC_SCHEMA_VERSION,
            query_version=query.query_version,
            params={
                "max_publications": query.max_publications,
                "matched_entities": list(query.matched_entities),
                "reason_shown": query.reason_shown,
                "source_intent_id": query.source_intent_id,
            },
        )
        for query in queries
    )
    return LiteraturePlan(
        **base,
        explicit_request=True,
        queries=queries,
        requests=requests,
        dropped=dropped,
    )


@dataclass(frozen=True, slots=True)
class Publication:
    publication_id: str
    title: str
    source: str
    retrieved_at: str
    pmid: str | None = None
    pmcid: str | None = None
    doi: str | None = None
    year: int | None = None
    journal: str | None = None
    authors_summary: str | None = None
    abstract: str | None = None
    provenance: Mapping[str, Any] = field(default_factory=immutable_mapping)
    metadata: Mapping[str, Any] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        for field_name in ("publication_id", "title", "source", "retrieved_at"):
            object.__setattr__(self, field_name, _text(getattr(self, field_name), field_name))
        for field_name in ("pmid", "pmcid", "doi", "journal", "authors_summary", "abstract"):
            object.__setattr__(self, field_name, _optional_text(getattr(self, field_name)))
        if self.year is not None:
            if isinstance(self.year, bool) or not isinstance(self.year, int):
                raise TypeError("year must be an integer or None")
            if not 1500 <= self.year <= 3000:
                raise ValueError("year is outside the accepted publication range")
        object.__setattr__(self, "provenance", immutable_mapping(self.provenance))
        object.__setattr__(self, "metadata", immutable_mapping(self.metadata))

    @classmethod
    def from_provider_record(cls, value: Mapping[str, Any]) -> "Publication":
        return cls(
            publication_id=value["publication_id"],
            pmid=value.get("pmid"),
            pmcid=value.get("pmcid"),
            doi=value.get("doi"),
            title=value["title"],
            year=value.get("year"),
            journal=value.get("journal"),
            authors_summary=value.get("authors_summary"),
            abstract=value.get("abstract"),
            source=value.get("source", "Europe PMC"),
            retrieved_at=value["retrieved_at"],
            provenance=value.get("provenance", {}),
            metadata={
                "is_open_access": value.get("is_open_access"),
                "underlying_source": value.get("source"),
            },
        )


def publications_from_result(result: ProviderResult) -> tuple[Publication, ...]:
    if result.status is not ProviderStatus.AVAILABLE or not isinstance(result.data, Mapping):
        return ()
    records = result.data.get("records")
    if not isinstance(records, (tuple, list)):
        return ()
    publications: list[Publication] = []
    for record in records:
        if not isinstance(record, Mapping):
            continue
        try:
            publications.append(Publication.from_provider_record(record))
        except (KeyError, TypeError, ValueError):
            continue
    return tuple(publications)


@dataclass(frozen=True, slots=True)
class RelevantPassage:
    publication_id: str
    text: str
    source_section: str
    matched_entities: tuple[str, ...]
    matched_query_terms: tuple[str, ...]
    rank: int
    basis: str


def _tokens(text: str | None) -> tuple[str, ...]:
    return tuple(token.casefold() for token in _WORD.findall(text or ""))


def extract_relevant_passages(
    publication: Publication,
    *,
    query_text: str,
    entity_terms: Iterable[str],
    limit: int = 3,
) -> tuple[RelevantPassage, ...]:
    """Rank source sentences by transparent lexical overlap; no summarization."""

    if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
        raise ValueError("limit must be a positive integer")
    entity_map = {term.casefold(): term for term in _unique(entity_terms)}
    query_terms = {
        token for token in _tokens(query_text)
        if len(token) >= 2 and token not in _STOPWORDS and token not in {"and", "or"}
    }
    candidates: list[tuple[tuple[int, int, int, int], str, str, tuple[str, ...], tuple[str, ...]]] = []
    sections = (("title", publication.title), ("abstract", publication.abstract or ""))
    for section, source_text in sections:
        sentences = (source_text,) if section == "title" else tuple(_SENTENCE.split(source_text))
        for order, sentence in enumerate(sentences):
            sentence = " ".join(sentence.split())[:600]
            if not sentence:
                continue
            token_set = set(_tokens(sentence))
            entity_hits = tuple(
                original for folded, original in entity_map.items()
                if folded in token_set
            )
            query_hits = tuple(sorted(query_terms & token_set))
            if not entity_hits and not query_hits:
                continue
            sort_key = (
                len(entity_hits),
                len(query_hits),
                1 if section == "title" else 0,
                -order,
            )
            candidates.append((sort_key, section, sentence, entity_hits, query_hits))
    candidates.sort(key=lambda item: item[0], reverse=True)
    passages = []
    seen: set[str] = set()
    for _, section, sentence, entity_hits, query_hits in candidates:
        if sentence.casefold() in seen:
            continue
        seen.add(sentence.casefold())
        basis_parts = []
        if entity_hits:
            basis_parts.append("entity-term overlap: " + ", ".join(entity_hits))
        if query_hits:
            basis_parts.append("query-term overlap: " + ", ".join(query_hits))
        passages.append(RelevantPassage(
            publication_id=publication.publication_id,
            text=sentence,
            source_section=section,
            matched_entities=entity_hits,
            matched_query_terms=query_hits,
            rank=len(passages) + 1,
            basis="; ".join(basis_parts),
        ))
        if len(passages) >= limit:
            break
    return tuple(passages)


@dataclass(frozen=True, slots=True)
class LiteratureSearchHit:
    publication: Publication
    rank: int
    matched_terms: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "matched_terms", _unique(self.matched_terms))
        if isinstance(self.rank, bool) or not isinstance(self.rank, int) or self.rank <= 0:
            raise ValueError("rank must be a positive integer")


@dataclass(frozen=True, slots=True)
class LiteratureSearchResult:
    query_text: str
    method: str
    hits: tuple[LiteratureSearchHit, ...]
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "query_text", _text(self.query_text, "query_text"))
        object.__setattr__(self, "method", _text(self.method, "method"))
        object.__setattr__(self, "hits", tuple(self.hits))


@dataclass(frozen=True, slots=True)
class PublicationStoreResult:
    stored_count: int
    search_backend: str
    errors: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "errors", tuple(self.errors))


@dataclass(frozen=True, slots=True)
class PublicationMatch:
    publication: Publication
    query: LiteratureQuery
    relevant_passages: tuple[RelevantPassage, ...]
    match_label: str = "Co-mentioned"

    def __post_init__(self) -> None:
        object.__setattr__(self, "relevant_passages", tuple(self.relevant_passages))
        if self.match_label != "Co-mentioned":
            raise ValueError("literature matches must use the non-causal Co-mentioned label")


@dataclass(frozen=True, slots=True)
class LiteratureLoad:
    plan: LiteraturePlan
    results: tuple[ProviderResult, ...]
    matches: tuple[PublicationMatch, ...]
    store_result: PublicationStoreResult | None = None
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "results", tuple(self.results))
        object.__setattr__(self, "matches", tuple(self.matches))
        if len(self.results) != len(self.plan.requests):
            raise ValueError("one result is required for every literature request")


class LiteratureLoader:
    """Cache-first provider execution followed by optional normalized storage."""

    def __init__(
        self,
        client: CacheFirstProviderClient,
        *,
        repository: "LiteratureRepository | None" = None,
    ) -> None:
        if not isinstance(client, CacheFirstProviderClient):
            raise TypeError("client must be CacheFirstProviderClient")
        if repository is not None and not isinstance(repository, LiteratureRepository):
            raise TypeError("repository must be LiteratureRepository or None")
        self._client = client
        self._repository = repository

    def load(self, plan: LiteraturePlan, *, timeout: ProviderTimeout) -> LiteratureLoad:
        if not isinstance(plan, LiteraturePlan):
            raise TypeError("plan must be LiteraturePlan")
        if not isinstance(timeout, ProviderTimeout):
            raise TypeError("timeout must be ProviderTimeout")
        results: list[ProviderResult] = []
        matches: list[PublicationMatch] = []
        unique_publications: dict[str, Publication] = {}
        for query, request in zip(plan.queries, plan.requests):
            result = self._client.execute(
                request,
                timeout,
                evidence_class=CacheEvidenceClass.LITERATURE,
            )
            results.append(result)
            entity_terms = tuple(
                logical_key.rsplit(":", 1)[-1]
                for logical_key in query.matched_entities
            )
            for publication in publications_from_result(result):
                unique_publications.setdefault(publication.publication_id, publication)
                matches.append(PublicationMatch(
                    publication=publication,
                    query=query,
                    relevant_passages=extract_relevant_passages(
                        publication,
                        query_text=query.query_text,
                        entity_terms=entity_terms,
                    ),
                ))
        store_result = None
        if self._repository is not None and unique_publications:
            store_result = self._repository.store_publications(
                unique_publications.values(),
                legal_metadata=EuropePMCProvider.legal_metadata,
            )
        return LiteratureLoad(
            plan=plan,
            results=tuple(results),
            matches=tuple(matches),
            store_result=store_result,
            message=None if matches else NO_MATCHING_PUBLICATIONS_MESSAGE,
        )


class LiteratureRepository:
    """Normalized publication cache with optional FTS5 and TF-IDF fallback."""

    def __init__(self, store: EvidenceStore) -> None:
        if not isinstance(store, EvidenceStore):
            raise TypeError("store must be EvidenceStore")
        self._store = store

    def initialize_search(self) -> str:
        """Explicitly create the optional research-only FTS index."""

        try:
            with self._store.transaction() as connection:
                connection.execute(
                    "CREATE VIRTUAL TABLE IF NOT EXISTS publication_fts "
                    "USING fts5(publication_id UNINDEXED, title, abstract)"
                )
                connection.execute("DELETE FROM publication_fts")
                connection.execute(
                    "INSERT INTO publication_fts(rowid, publication_id, title, abstract) "
                    "SELECT id, publication_id, title, COALESCE(abstract, '') FROM publications"
                )
            return "fts5"
        except EvidenceStoreError:
            return "tfidf_fallback"

    @staticmethod
    def _fts_exists(connection: sqlite3.Connection) -> bool:
        return connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='publication_fts'"
        ).fetchone() is not None

    def store_publications(
        self,
        publications: Iterable[Publication],
        *,
        legal_metadata: Mapping[str, Any] | None = None,
    ) -> PublicationStoreResult:
        frozen = tuple(publications)
        if not all(isinstance(item, Publication) for item in frozen):
            raise TypeError("publications must contain Publication values")
        try:
            source_id = self._store.upsert_source(
                source_key="external:europepmc",
                source_name="Europe PMC",
                provider="europepmc",
                attribution="Europe PMC",
                license_note=(legal_metadata or {}).get("license_note"),
                metadata=dict(legal_metadata or {}),
            )
            stored = 0
            with self._store.transaction() as connection:
                fts = self._fts_exists(connection)
                for publication in frozen:
                    row = connection.execute(
                        """
                        SELECT id, publication_id FROM publications
                        WHERE publication_id = ?
                           OR (? IS NOT NULL AND pmid = ?)
                           OR (? IS NOT NULL AND pmcid = ?)
                           OR (? IS NOT NULL AND doi = ?)
                        ORDER BY id LIMIT 1
                        """,
                        (
                            publication.publication_id,
                            publication.pmid, publication.pmid,
                            publication.pmcid, publication.pmcid,
                            publication.doi, publication.doi,
                        ),
                    ).fetchone()
                    metadata_json = json.dumps(
                        {"provenance": dict(publication.provenance), **dict(publication.metadata)},
                        ensure_ascii=False,
                        sort_keys=True,
                        separators=(",", ":"),
                        allow_nan=False,
                    )
                    values = (
                        publication.publication_id, publication.pmid, publication.pmcid,
                        publication.doi, publication.title, publication.year,
                        publication.journal, publication.authors_summary, publication.abstract,
                        source_id, publication.retrieved_at, metadata_json,
                    )
                    if row is None:
                        cursor = connection.execute(
                            """
                            INSERT INTO publications(
                                publication_id, pmid, pmcid, doi, title, year, journal,
                                authors_summary, abstract, source_id, retrieved_at, metadata_json
                            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                            """,
                            values,
                        )
                        row_id = int(cursor.lastrowid)
                        stored_publication_id = publication.publication_id
                    else:
                        row_id = int(row[0])
                        stored_publication_id = str(row["publication_id"])
                        connection.execute(
                            """
                            UPDATE publications SET
                                pmid=?, pmcid=?, doi=?, title=?, year=?, journal=?,
                                authors_summary=?, abstract=?, source_id=?, retrieved_at=?, metadata_json=?
                            WHERE id=?
                            """,
                            (*values[1:], row_id),
                        )
                    if fts:
                        connection.execute("DELETE FROM publication_fts WHERE rowid = ?", (row_id,))
                        connection.execute(
                            "INSERT INTO publication_fts(rowid, publication_id, title, abstract) VALUES (?, ?, ?, ?)",
                            (row_id, stored_publication_id, publication.title, publication.abstract or ""),
                        )
                    stored += 1
            return PublicationStoreResult(stored, "fts5" if fts else "tfidf_fallback")
        except (EvidenceStoreError, TypeError, ValueError) as exc:
            return PublicationStoreResult(0, "unavailable", (f"{type(exc).__name__}: {exc}",))

    @staticmethod
    def _publication(row: sqlite3.Row) -> Publication:
        metadata = json.loads(row["metadata_json"])
        provenance = metadata.pop("provenance", {})
        underlying_source = metadata.pop("underlying_source", None)
        return Publication(
            publication_id=row["publication_id"], pmid=row["pmid"], pmcid=row["pmcid"],
            doi=row["doi"], title=row["title"], year=row["year"], journal=row["journal"],
            authors_summary=row["authors_summary"], abstract=row["abstract"],
            source=underlying_source or row["source_name"], retrieved_at=row["retrieved_at"],
            provenance=provenance, metadata=metadata,
        )

    @staticmethod
    def _select_rows(connection: sqlite3.Connection, ids: Sequence[int]) -> dict[int, sqlite3.Row]:
        if not ids:
            return {}
        placeholders = ",".join("?" for _ in ids)
        rows = connection.execute(
            f"""
            SELECT p.*, s.source_name FROM publications p
            JOIN sources s ON s.id = p.source_id
            WHERE p.id IN ({placeholders})
            """,
            tuple(ids),
        ).fetchall()
        return {int(row["id"]): row for row in rows}

    def search(self, query_text: str, *, limit: int = 20) -> LiteratureSearchResult:
        query = _text(query_text, "query_text")
        if isinstance(limit, bool) or not isinstance(limit, int) or limit <= 0:
            raise ValueError("limit must be a positive integer")
        limit = min(limit, 100)
        terms = _unique(
            token for token in _tokens(query)
            if len(token) >= 2 and token not in _STOPWORDS
        )
        if not terms:
            return LiteratureSearchResult(query, "none", (), NO_MATCHING_PUBLICATIONS_MESSAGE)
        try:
            with self._store.read_connection() as connection:
                if self._fts_exists(connection):
                    try:
                        expression = " OR ".join(_quote(term) for term in terms)
                        ranked = connection.execute(
                            "SELECT rowid FROM publication_fts WHERE publication_fts MATCH ? "
                            "ORDER BY bm25(publication_fts) LIMIT ?",
                            (expression, limit),
                        ).fetchall()
                        ids = [int(row[0]) for row in ranked]
                        rows = self._select_rows(connection, ids)
                        hits = tuple(
                            LiteratureSearchHit(
                                self._publication(rows[row_id]),
                                rank=index + 1,
                                matched_terms=tuple(
                                    term for term in terms
                                    if term.casefold() in set(_tokens(
                                        f"{rows[row_id]['title']} {rows[row_id]['abstract'] or ''}"
                                    ))
                                ),
                            )
                            for index, row_id in enumerate(ids)
                            if row_id in rows
                        )
                        return LiteratureSearchResult(
                            query, "fts5_bm25", hits,
                            None if hits else NO_MATCHING_PUBLICATIONS_MESSAGE,
                        )
                    except sqlite3.Error:
                        pass
                rows = connection.execute(
                    """
                    SELECT p.*, s.source_name FROM publications p
                    JOIN sources s ON s.id = p.source_id
                    ORDER BY p.retrieved_at DESC, p.id DESC LIMIT 500
                    """
                ).fetchall()
        except (EvidenceStoreError, sqlite3.Error) as exc:
            return LiteratureSearchResult(query, "unavailable", (), str(exc))

        documents = [set(_tokens(f"{row['title']} {row['abstract'] or ''}")) for row in rows]
        document_frequency = Counter(
            term for document in documents for term in set(terms) & document
        )
        ranked_rows = []
        for index, (row, document) in enumerate(zip(rows, documents)):
            matched = tuple(term for term in terms if term in document)
            if not matched:
                continue
            tfidf = sum(math.log((len(rows) + 1) / (document_frequency[term] + 1)) + 1 for term in matched)
            ranked_rows.append((tfidf, -index, row, matched))
        ranked_rows.sort(key=lambda item: (item[0], item[1]), reverse=True)
        hits = tuple(
            LiteratureSearchHit(self._publication(row), rank=index + 1, matched_terms=matched)
            for index, (_, __, row, matched) in enumerate(ranked_rows[:limit])
        )
        return LiteratureSearchResult(
            query, "tfidf_fallback", hits,
            None if hits else NO_MATCHING_PUBLICATIONS_MESSAGE,
        )


__all__ = [
    "LITERATURE_QUERY_VERSION",
    "NO_MATCHING_PUBLICATIONS_MESSAGE",
    "LiteraturePlan",
    "LiteratureLoad",
    "LiteratureLoader",
    "LiteratureQuery",
    "LiteratureQueryBuilder",
    "LiteratureRepository",
    "LiteratureSearchHit",
    "LiteratureSearchResult",
    "Publication",
    "PublicationMatch",
    "PublicationStoreResult",
    "RelevantPassage",
    "extract_relevant_passages",
    "plan_literature_requests",
    "publications_from_result",
]
