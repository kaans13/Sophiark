from __future__ import annotations

import ast
from concurrent.futures import ThreadPoolExecutor
from dataclasses import replace
from datetime import datetime, timedelta, timezone
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

import pytest

from src.research.config import CacheTTLs, ProviderTimeout
from src.research.evidence_store import EvidenceStore
from src.research.provider_cache import (
    CacheEvidenceClass,
    CacheFirstProviderClient,
    ProviderCache,
    ProviderCacheKey,
    ProviderTTLPolicy,
)
from src.research.provider_runtime import ProviderRequest, ProviderResult, ProviderStatus


class MutableClock:
    def __init__(self) -> None:
        self.value = datetime(2026, 8, 25, 12, 0, tzinfo=timezone.utc)

    def __call__(self) -> datetime:
        return self.value

    def advance(self, seconds: int) -> None:
        self.value += timedelta(seconds=seconds)


class CountingExecutor:
    def __init__(self, result_factory):
        self.calls = 0
        self.provider = SimpleNamespace(provider_name="test-provider")
        self._result_factory = result_factory

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        self.calls += 1
        assert timeout.total_seconds == 8.0
        return self._result_factory(request)


def _request(**changes) -> ProviderRequest:
    values = {
        "provider": "uniprot",
        "operation": "protein_annotation",
        "taxon_id": 9606,
        "canonical_query": "UniProtKB:P43220",
        "provider_schema_version": "2",
        "query_version": "research-query-v4",
    }
    values.update(changes)
    return ProviderRequest(**values)


def _result(
    request: ProviderRequest,
    *,
    status: ProviderStatus = ProviderStatus.AVAILABLE,
    data=None,
) -> ProviderResult:
    if data is None and status is ProviderStatus.AVAILABLE:
        data = {"accession": "P43220", "names": ["GLP1R"]}
    return ProviderResult(
        request=request,
        status=status,
        data=data,
        attempts=1,
        message=None if status is ProviderStatus.AVAILABLE else status.value,
        provider_version="2026_03",
        metadata={"source": "fixture"},
    )


def _short_policy() -> ProviderTTLPolicy:
    return ProviderTTLPolicy(
        CacheTTLs(
            entity_resolution_seconds=40,
            biological_annotation_seconds=30,
            pathway_annotation_seconds=20,
            literature_seconds=10,
            negative_result_seconds=5,
            transient_failure_seconds=2,
        )
    )


def test_cache_key_is_stable_and_every_required_identity_field_separates() -> None:
    base = _request()
    base_key = ProviderCacheKey.from_request(base)
    assert base_key.digest == ProviderCacheKey.from_request(_request()).digest
    assert base_key.cache_key == base_key.digest
    assert dict(base_key.identity) == {
        "provider": "uniprot",
        "operation": "protein_annotation",
        "taxon_id": 9606,
        "canonical_query": "UniProtKB:P43220",
        "provider_schema_version": "2",
        "query_version": "research-query-v4",
        "params": {},
    }
    assert base_key.digest == base.cache_key

    variants = (
        replace(base, provider="interpro"),
        replace(base, operation="domains"),
        replace(base, taxon_id=10090),
        replace(base, canonical_query="UniProtKB:P11111"),
        replace(base, provider_schema_version="3"),
        replace(base, query_version="research-query-v5"),
        replace(base, params={"fields": ["domains"]}),
    )
    assert all(ProviderCacheKey.from_request(item).digest != base_key.digest for item in variants)
    assert len({ProviderCacheKey.from_request(item).digest for item in (*variants, base)}) == 8


def test_l1_hit_returns_exact_status_and_immutable_json_payload() -> None:
    request = _request()
    cache = ProviderCache(max_memory_entries=4)
    original = _result(request)

    assert cache.put(original) is True
    cached = cache.get(request)

    assert cached is original
    assert cached.request == request
    assert cached.status is ProviderStatus.AVAILABLE
    assert cached.to_dict()["data"] == {
        "accession": "P43220",
        "names": ["GLP1R"],
    }
    assert cache.stats.l1_hits == 1
    assert cache.stats.l2_hits == 0


def test_l2_persists_across_instances_and_promotes_to_l1() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()  # Cache itself must never initialize this database.
        request = _request()
        first = ProviderCache(store=store)
        assert first.put(_result(request)) is True
        first.clear_memory()

        second = ProviderCache(store=store)
        l2_result = second.get(request)
        assert l2_result is not None
        assert l2_result.request == request
        assert l2_result.status is ProviderStatus.AVAILABLE
        assert l2_result.data["accession"] == "P43220"
        assert dict(l2_result.metadata) == {"source": "fixture"}
        assert second.stats.l2_hits == 1

        assert second.get(request) is l2_result
        assert second.stats.l1_hits == 1


def test_same_six_fields_with_different_params_are_independent_in_l1_and_l2() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        first_request = _request(params={"fields": ["names"]})
        second_request = _request(params={"fields": ["domains"]})
        cache = ProviderCache(store=store)
        cache.put(_result(first_request, data={"kind": "names"}))
        cache.put(_result(second_request, data={"kind": "domains"}))

        assert cache.get(first_request).data["kind"] == "names"
        assert cache.get(second_request).data["kind"] == "domains"
        cache.clear_memory()

        reader = ProviderCache(store=store)
        assert reader.get(first_request).data["kind"] == "names"
        assert reader.get(second_request).data["kind"] == "domains"
        with store.read_connection() as connection:
            rows = connection.execute(
                "SELECT request_params_json, request_params_hash FROM provider_cache"
            ).fetchall()
        assert len(rows) == 2
        assert len({row["request_params_hash"] for row in rows}) == 2


def test_expired_l1_and_l2_entries_are_misses() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        clock = MutableClock()
        cache = ProviderCache(store=store, ttl_policy=_short_policy(), clock=clock)
        request = _request()
        cache.put(_result(request))

        clock.advance(30)
        assert cache.get(request) is None
        assert cache.stats.misses == 1
        assert cache.stats.l1_hits == 0
        assert cache.stats.l2_hits == 0


def test_no_match_uses_short_negative_cache_without_repeat_provider_call() -> None:
    clock = MutableClock()
    cache = ProviderCache(ttl_policy=_short_policy(), clock=clock)
    executor = CountingExecutor(
        lambda request: _result(request, status=ProviderStatus.NO_MATCH, data=None)
    )
    client = CacheFirstProviderClient(cache=cache, executor=executor)
    request = _request(canonical_query="UniProtKB:UNKNOWN")

    first = client.execute(request, ProviderTimeout())
    second = client.execute(request, ProviderTimeout())
    assert first.status is ProviderStatus.NO_MATCH
    assert second.status is ProviderStatus.NO_MATCH
    assert executor.calls == 1

    clock.advance(5)
    third = client.execute(request, ProviderTimeout())
    assert third.status is ProviderStatus.NO_MATCH
    assert executor.calls == 2


def test_no_match_persists_with_utc_negative_expiry_in_l2() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        clock = MutableClock()
        request = _request(canonical_query="UniProtKB:UNKNOWN")
        writer = ProviderCache(store=store, ttl_policy=_short_policy(), clock=clock)
        writer.put(_result(request, status=ProviderStatus.NO_MATCH, data=None))
        writer.clear_memory()

        reader = ProviderCache(store=store, ttl_policy=_short_policy(), clock=clock)
        restored = reader.get(request)
        assert restored is not None
        assert restored.status is ProviderStatus.NO_MATCH
        with store.read_connection() as connection:
            row = connection.execute(
                "SELECT checked_at, expires_at, status FROM provider_cache"
            ).fetchone()
        checked_at = datetime.fromisoformat(row["checked_at"])
        expires_at = datetime.fromisoformat(row["expires_at"])
        assert checked_at.utcoffset() == timedelta(0)
        assert expires_at.utcoffset() == timedelta(0)
        assert (expires_at - checked_at).total_seconds() == 5
        assert row["status"] == ProviderStatus.NO_MATCH.value


@pytest.mark.parametrize(
    "status",
    [
        ProviderStatus.OFFLINE,
        ProviderStatus.TIMEOUT,
        ProviderStatus.RATE_LIMIT,
        ProviderStatus.PROVIDER_ERROR,
    ],
)
def test_transient_failures_use_central_short_ttl(status: ProviderStatus) -> None:
    clock = MutableClock()
    cache = ProviderCache(ttl_policy=_short_policy(), clock=clock)
    executor = CountingExecutor(lambda request: _result(request, status=status, data=None))
    client = CacheFirstProviderClient(cache=cache, executor=executor)
    request = _request()

    assert client.execute(request, ProviderTimeout()).status is status
    assert client.execute(request, ProviderTimeout()).status is status
    assert executor.calls == 1
    clock.advance(2)
    assert client.execute(request, ProviderTimeout()).status is status
    assert executor.calls == 2


def test_ttl_policy_has_distinct_configurable_positive_classes() -> None:
    policy = _short_policy()
    request = _request()
    result = _result(request)
    assert policy.ttl_seconds(result) == 30
    assert policy.ttl_seconds(result, evidence_class=CacheEvidenceClass.ENTITY_RESOLUTION) == 40
    assert policy.ttl_seconds(result, evidence_class=CacheEvidenceClass.PATHWAY_ANNOTATION) == 20
    assert policy.ttl_seconds(result, evidence_class=CacheEvidenceClass.LITERATURE) == 10
    assert policy.ttl_seconds(_result(request, status=ProviderStatus.NO_MATCH)) == 5
    assert policy.ttl_seconds(_result(request, status=ProviderStatus.TIMEOUT)) == 2


def test_provider_specific_ttls_cover_positive_literature_negative_and_transient() -> None:
    uniprot_ttls = CacheTTLs(
        entity_resolution_seconds=11,
        biological_annotation_seconds=12,
        pathway_annotation_seconds=13,
        literature_seconds=14,
        negative_result_seconds=15,
        transient_failure_seconds=16,
    )
    europepmc_ttls = CacheTTLs(
        entity_resolution_seconds=21,
        biological_annotation_seconds=22,
        pathway_annotation_seconds=23,
        literature_seconds=24,
        negative_result_seconds=25,
        transient_failure_seconds=26,
    )
    policy = ProviderTTLPolicy(
        ttls=_short_policy().ttls,
        provider_ttls={
            " UniProt ": uniprot_ttls,
            "EUROPEPMC": europepmc_ttls,
        },
    )
    uniprot = _request(provider="UNIPROT")
    uniprot_literature = _request(
        provider="uniprot",
        operation="literature_search",
    )
    europepmc = _request(
        provider="EuropePMC",
        operation="literature_search",
    )

    assert policy.ttl_seconds(_result(uniprot)) == 12
    assert policy.ttl_seconds(_result(uniprot_literature)) == 14
    assert policy.ttl_seconds(_result(europepmc)) == 24
    assert policy.ttl_seconds(
        _result(uniprot, status=ProviderStatus.NO_MATCH)
    ) == 15
    assert policy.ttl_seconds(
        _result(europepmc, status=ProviderStatus.NO_MATCH)
    ) == 25
    assert policy.ttl_seconds(
        _result(uniprot, status=ProviderStatus.TIMEOUT)
    ) == 16
    assert policy.ttl_seconds(
        _result(europepmc, status=ProviderStatus.RATE_LIMIT)
    ) == 26
    # Providers without an override retain the existing global fallback.
    assert policy.ttl_seconds(_result(_request(provider="quickgo"))) == 30


def test_provider_ttl_mapping_is_detached_immutable_normalized_and_validated() -> None:
    override = CacheTTLs(
        entity_resolution_seconds=71,
        biological_annotation_seconds=72,
        pathway_annotation_seconds=73,
        literature_seconds=74,
        negative_result_seconds=75,
        transient_failure_seconds=76,
    )
    source = {" UniProt ": override}
    policy = ProviderTTLPolicy(provider_ttls=source)
    source.clear()

    assert policy.ttls_for("UNIPROT") is override
    assert tuple(policy.provider_ttls) == ("uniprot",)
    with pytest.raises(TypeError):
        policy.provider_ttls["quickgo"] = override
    with pytest.raises(TypeError, match="must be a mapping"):
        ProviderTTLPolicy(provider_ttls=[])
    with pytest.raises(ValueError, match="non-empty"):
        ProviderTTLPolicy(provider_ttls={" ": override})
    with pytest.raises(TypeError, match="CacheTTLs"):
        ProviderTTLPolicy(provider_ttls={"uniprot": "not-ttls"})
    with pytest.raises(ValueError, match="Duplicate normalized"):
        ProviderTTLPolicy(
            provider_ttls={"UniProt": override, " uniprot ": CacheTTLs()}
        )


def test_provider_ttl_addition_preserves_positional_operation_class_api() -> None:
    policy = ProviderTTLPolicy(
        _short_policy().ttls,
        {"custom_search": CacheEvidenceClass.LITERATURE},
    )
    request = _request(operation="custom_search")
    assert policy.ttl_seconds(_result(request)) == 10


def test_disabled_result_is_never_stored_in_memory_or_sqlite() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        cache = ProviderCache(store=store)
        request = _request()

        assert cache.put(_result(request, status=ProviderStatus.DISABLED, data=None)) is False
        assert cache.memory_size == 0
        with store.read_connection() as connection:
            count = connection.execute("SELECT COUNT(*) FROM provider_cache").fetchone()[0]
        assert count == 0


def test_missing_l2_is_fail_soft_never_initialized_and_allows_l3() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "not-initialized.sqlite"
        cache = ProviderCache(store=EvidenceStore(path))
        executor = CountingExecutor(_result)
        client = CacheFirstProviderClient(cache=cache, executor=executor)
        request = _request()

        assert client.execute(request, ProviderTimeout()).status is ProviderStatus.AVAILABLE
        assert executor.calls == 1
        assert path.exists() is False
        assert cache.stats.l2_failures >= 1
        assert client.execute(request, ProviderTimeout()).status is ProviderStatus.AVAILABLE
        assert executor.calls == 1


def test_corrupt_l2_is_fail_soft_and_does_not_block_l3_or_l1() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "corrupt.sqlite"
        original = b"not a sqlite database"
        path.write_bytes(original)
        cache = ProviderCache(store=EvidenceStore(path))
        executor = CountingExecutor(_result)
        client = CacheFirstProviderClient(cache=cache, executor=executor)
        request = _request()

        assert client.execute(request, ProviderTimeout()).status is ProviderStatus.AVAILABLE
        assert client.execute(request, ProviderTimeout()).status is ProviderStatus.AVAILABLE
        assert executor.calls == 1
        assert path.read_bytes() == original
        assert cache.stats.l2_failures >= 1


def test_malformed_l2_row_is_a_miss_not_an_exception() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        request = _request()
        writer = ProviderCache(store=store)
        writer.put(_result(request))
        writer.clear_memory()
        with store.transaction() as connection:
            connection.execute(
                "UPDATE provider_cache SET expires_at = 'invalid timestamp'"
            )

        reader = ProviderCache(store=store)
        assert reader.get(request) is None
        assert reader.stats.l2_failures == 1
        assert reader.stats.misses == 1


@pytest.mark.parametrize("bad_data", [{"nan": float("nan")}, {"object": object()}])
def test_non_json_provider_payload_is_rejected_before_cache_mutation(bad_data) -> None:
    cache = ProviderCache()
    request = _request()
    with pytest.raises((TypeError, ValueError)):
        ProviderResult(
            request=request,
            status=ProviderStatus.AVAILABLE,
            data=bad_data,
        )
    assert cache.memory_size == 0


def test_nested_json_payload_and_metadata_round_trip_safely() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research.sqlite")
        store.initialize()
        request = _request(params={"include": ["names", "domains"]})
        result = ProviderResult(
            request=request,
            status=ProviderStatus.AVAILABLE,
            data={"items": [{"id": "P43220", "flags": [True, None]}]},
            metadata={"page": {"number": 1, "cursor": None}},
        )
        ProviderCache(store=store).put(result)

        restored = ProviderCache(store=store).get(request)
        assert restored is not None
        assert restored.to_dict()["data"] == {
            "items": [{"id": "P43220", "flags": [True, None]}]
        }
        assert restored.to_dict()["metadata"] == {
            "page": {"number": 1, "cursor": None}
        }


def test_l1_is_bounded_and_evicts_least_recently_used_entry() -> None:
    cache = ProviderCache(max_memory_entries=2)
    requests = [_request(canonical_query=f"UniProtKB:P{index}") for index in range(3)]
    for request in requests[:2]:
        cache.put(_result(request))

    # Refresh the first entry, making the second the least recently used.
    assert cache.get(requests[0]) is not None
    cache.put(_result(requests[2]))

    assert cache.memory_size == 2
    assert cache.get(requests[1]) is None
    assert cache.get(requests[0]) is not None
    assert cache.get(requests[2]) is not None
    assert cache.stats.evictions == 1


def test_bounded_l1_operations_are_thread_safe() -> None:
    cache = ProviderCache(max_memory_entries=16)

    def cache_one(index: int) -> None:
        request = _request(canonical_query=f"UniProtKB:THREAD-{index}")
        cache.put(_result(request, data={"index": index}))
        cache.get(request)

    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(cache_one, range(100)))

    assert cache.memory_size <= 16
    assert cache.stats.evictions >= 84


def test_client_calls_l3_once_then_uses_cache_without_status_rewrite() -> None:
    cache = ProviderCache()
    request = _request()
    executor = CountingExecutor(_result)
    client = CacheFirstProviderClient(cache=cache, executor=executor)

    first = client.execute(request, ProviderTimeout())
    second = client.execute(request, ProviderTimeout())
    assert executor.calls == 1
    assert second is first
    assert second.request == request
    assert second.status is ProviderStatus.AVAILABLE


def test_provider_cache_imports_no_ui_network_or_scientific_modules() -> None:
    source_path = Path(__file__).parents[1] / "src" / "research" / "provider_cache.py"
    tree = ast.parse(source_path.read_text(encoding="utf-8"))
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)

    forbidden = ("streamlit", "requests", "httpx", "urllib", "scientific")
    assert not any(
        name == blocked or name.startswith(blocked + ".")
        for name in imported
        for blocked in forbidden
    )
