"""Acceptance tests for the transport-agnostic Phase 6 provider runtime."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError
import json
from pathlib import Path
from threading import Lock
import time
from typing import Any, Callable

import pytest

from src.research.config import ProviderTimeout
from src.research.provider_runtime import (
    EvidenceProvider,
    ProviderCapability,
    ProviderExecutionPolicy,
    ProviderExecutor,
    ProviderRequest,
    ProviderResult,
    ProviderStatus,
    SessionCircuitBreaker,
)


_TIMEOUT = ProviderTimeout(
    connect_seconds=0.01,
    read_seconds=0.01,
    total_seconds=0.03,
)


class ScriptedProvider(EvidenceProvider):
    def __init__(
        self,
        actions: list[Any] | None = None,
        *,
        enabled: bool = True,
        capabilities: tuple[ProviderCapability, ...] | None = None,
        dynamic: Callable[[ProviderRequest], Any] | None = None,
    ) -> None:
        super().__init__(
            "fixture",
            capabilities or (ProviderCapability("lookup", (9606,)),),
            provider_version="fixture-2026.1",
            enabled=enabled,
        )
        self.actions = list(actions or [])
        self.dynamic = dynamic
        self.calls = 0
        self.active = 0
        self.max_active = 0
        self._lock = Lock()

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        assert timeout is _TIMEOUT
        with self._lock:
            self.calls += 1
            self.active += 1
            self.max_active = max(self.max_active, self.active)
            action = self.dynamic(request) if self.dynamic is not None else self.actions.pop(0)
        try:
            if isinstance(action, tuple) and action[0] == "sleep":
                time.sleep(float(action[1]))
                action = action[2]
            if isinstance(action, BaseException):
                raise action
            if callable(action):
                action = action(request)
            if isinstance(action, ProviderResult):
                return action
            status = ProviderStatus(action)
            return ProviderResult(
                request=request,
                status=status,
                data={"query": request.canonical_query} if status is ProviderStatus.AVAILABLE else None,
                provider_version=self.provider_version,
            )
        finally:
            with self._lock:
                self.active -= 1


def _request(
    query: str = "ENSP00000407554",
    *,
    provider: str = "fixture",
    operation: str = "lookup",
    taxon_id: int = 9606,
    params: dict[str, Any] | None = None,
) -> ProviderRequest:
    return ProviderRequest(
        provider=provider,
        operation=operation,
        taxon_id=taxon_id,
        canonical_query=query,
        provider_schema_version="provider-schema-2",
        query_version="query-3",
        params=params or {"fields": ["family", "function"]},
    )


def _executor(
    provider: ScriptedProvider,
    *,
    retries: int = 0,
    threshold: int = 3,
    concurrency: int = 2,
    breaker: SessionCircuitBreaker | None = None,
) -> ProviderExecutor:
    return ProviderExecutor(
        provider,
        policy=ProviderExecutionPolicy(
            max_retries=retries,
            failure_threshold=threshold,
            max_concurrency=concurrency,
        ),
        circuit_breaker=breaker,
    )


def test_provider_status_contract_is_exact_and_distinct() -> None:
    assert tuple(status.value for status in ProviderStatus) == (
        "AVAILABLE",
        "NO_MATCH",
        "OFFLINE",
        "TIMEOUT",
        "RATE_LIMIT",
        "PROVIDER_ERROR",
        "DISABLED",
    )
    assert len(set(ProviderStatus)) == 7
    assert ProviderStatus.NO_MATCH is not ProviderStatus.PROVIDER_ERROR


@pytest.mark.parametrize("status", tuple(ProviderStatus))
def test_every_status_round_trips_without_collapsing(status: ProviderStatus) -> None:
    request = _request()
    result = ProviderResult(
        request=request,
        status=status,
        data={"nested": [1, {"ok": True}]},
        attempts=0,
        message=f"fixture {status.value}",
        provider_version="v1",
        metadata={"source": {"release": 3}},
    )

    encoded = json.dumps(result.to_dict(), allow_nan=False, sort_keys=True)
    restored = ProviderResult.from_dict(json.loads(encoded))
    assert restored.status is status
    assert restored.request.cache_key == request.cache_key
    assert restored.to_dict() == result.to_dict()


def test_request_identity_and_result_payload_are_detached_and_immutable() -> None:
    caller_params = {"fields": ["family"], "nested": {"limit": 5}}
    caller_data = {"records": [{"id": "P43220"}]}
    request = _request(params=caller_params)
    result = ProviderResult(request, ProviderStatus.AVAILABLE, caller_data, metadata={"page": 1})

    caller_params["fields"].append("mutated")
    caller_params["nested"]["limit"] = 999
    caller_data["records"][0]["id"] = "MUTATED"
    assert request.to_dict()["params"] == {"fields": ["family"], "nested": {"limit": 5}}
    assert result.to_dict()["data"] == {"records": [{"id": "P43220"}]}
    assert request.identity["provider_schema_version"] == "provider-schema-2"
    assert request.identity["query_version"] == "query-3"
    with pytest.raises(TypeError):
        request.params["new"] = True  # type: ignore[index]
    with pytest.raises(TypeError):
        result.data["new"] = True  # type: ignore[index]
    with pytest.raises(FrozenInstanceError):
        result.status = ProviderStatus.NO_MATCH  # type: ignore[misc]
    with pytest.raises((TypeError, ValueError)):
        _request(params={"bad": float("nan")})


def test_cache_key_uses_every_required_identity_field_and_params() -> None:
    base = _request()
    variants = (
        _request(provider="other"),
        _request(operation="other"),
        _request(taxon_id=10090),
        _request(query="OTHER"),
        ProviderRequest(**{**base.to_dict(), "provider_schema_version": "different"}),
        ProviderRequest(**{**base.to_dict(), "query_version": "different"}),
        _request(params={"fields": ["other"]}),
    )
    assert len({base.cache_key, *(variant.cache_key for variant in variants)}) == len(variants) + 1


def test_capabilities_are_explicit_taxon_safe_and_healthcheck_is_local() -> None:
    capabilities = (
        ProviderCapability("lookup", (9606, 10090), "entity annotation"),
        ProviderCapability("family", (), "taxon-independent family description"),
    )
    provider = ScriptedProvider([ProviderStatus.AVAILABLE], capabilities=capabilities)

    assert provider.capabilities == capabilities
    assert provider.supports(_request())
    assert provider.supports(_request(taxon_id=10090))
    assert not provider.supports(_request(taxon_id=10116))
    assert provider.supports(_request(operation="family", taxon_id=10116))
    assert provider.healthcheck() is ProviderStatus.AVAILABLE
    assert ScriptedProvider(enabled=False).healthcheck() is ProviderStatus.DISABLED


@pytest.mark.parametrize("invalid", (True, 9606.5, "9606"))
def test_taxon_identity_rejects_lossy_or_implicit_coercion(invalid: Any) -> None:
    with pytest.raises(TypeError, match="integer"):
        _request(taxon_id=invalid)  # type: ignore[arg-type]
    with pytest.raises(TypeError, match="integer"):
        ProviderCapability("lookup", (invalid,))  # type: ignore[arg-type]


@pytest.mark.parametrize(
    ("provider_request", "enabled", "reason"),
    (
        (_request(provider="other"), True, "provider_mismatch"),
        (_request(operation="unsupported"), True, "unsupported_capability"),
        (_request(taxon_id=10090), True, "unsupported_capability"),
        (_request(), False, "provider_disabled"),
    ),
)
def test_disabled_mismatched_and_unsupported_requests_never_execute(
    provider_request: ProviderRequest,
    enabled: bool,
    reason: str,
) -> None:
    provider = ScriptedProvider([ProviderStatus.AVAILABLE], enabled=enabled)
    executor = _executor(provider)
    try:
        result = executor.execute(provider_request, _TIMEOUT)
    finally:
        executor.close()

    assert result.status is ProviderStatus.DISABLED
    assert result.attempts == 0
    assert result.request.cache_key == provider_request.cache_key
    assert result.metadata["reason"] == reason
    assert provider.calls == 0


def test_deadline_timeout_returns_promptly_and_does_not_escape() -> None:
    provider = ScriptedProvider([("sleep", 0.12, ProviderStatus.AVAILABLE)])
    executor = _executor(provider)
    started = time.perf_counter()
    result = executor.execute(_request(), _TIMEOUT)
    elapsed = time.perf_counter() - started

    assert result.status is ProviderStatus.TIMEOUT
    assert result.attempts == 1
    assert result.metadata["reason"] == "deadline_exceeded"
    assert elapsed < 0.1
    executor.close()


def test_context_cleanup_does_not_wait_for_an_overdue_provider_thread() -> None:
    provider = ScriptedProvider([("sleep", 0.20, ProviderStatus.AVAILABLE)])
    started = time.perf_counter()
    with _executor(provider) as executor:
        result = executor.execute(_request(), _TIMEOUT)
    elapsed = time.perf_counter() - started

    assert result.status is ProviderStatus.TIMEOUT
    assert elapsed < 0.12


@pytest.mark.parametrize(
    ("exception", "expected"),
    (
        (TimeoutError("socket read timed out"), ProviderStatus.TIMEOUT),
        (ConnectionError("network unavailable"), ProviderStatus.OFFLINE),
        (RuntimeError("malformed provider response"), ProviderStatus.PROVIDER_ERROR),
    ),
)
def test_provider_exceptions_map_to_distinct_fail_soft_statuses(
    exception: Exception,
    expected: ProviderStatus,
) -> None:
    provider = ScriptedProvider([exception])
    executor = _executor(provider)
    try:
        result = executor.execute(_request(), _TIMEOUT)
    finally:
        executor.close()
    assert result.status is expected
    assert result.attempts == 1
    assert result.metadata["exception_type"] == type(exception).__name__


def test_default_policy_allows_exactly_one_retry_and_success_resets_failures() -> None:
    provider = ScriptedProvider([
        ProviderStatus.PROVIDER_ERROR,
        ProviderStatus.AVAILABLE,
    ])
    breaker = SessionCircuitBreaker(failure_threshold=3)
    executor = _executor(provider, retries=1, breaker=breaker)
    try:
        result = executor.execute(_request(), _TIMEOUT)
    finally:
        executor.close()

    assert result.status is ProviderStatus.AVAILABLE
    assert result.attempts == 2
    assert provider.calls == 2
    assert breaker.failure_count("fixture") == 0
    assert not breaker.is_open("fixture")


def test_retry_never_exceeds_two_total_attempts() -> None:
    provider = ScriptedProvider([
        ProviderStatus.PROVIDER_ERROR,
        ProviderStatus.PROVIDER_ERROR,
        ProviderStatus.AVAILABLE,
    ])
    executor = _executor(provider, retries=1, threshold=10)
    try:
        result = executor.execute(_request(), _TIMEOUT)
    finally:
        executor.close()
    assert result.status is ProviderStatus.PROVIDER_ERROR
    assert result.attempts == 2
    assert provider.calls == 2
    with pytest.raises(ValueError, match="zero or one"):
        ProviderExecutionPolicy(max_retries=2)


def test_circuit_opens_on_consecutive_timeout_or_provider_errors_and_can_reset() -> None:
    provider = ScriptedProvider([
        ProviderStatus.TIMEOUT,
        ProviderStatus.PROVIDER_ERROR,
        ProviderStatus.TIMEOUT,
        ProviderStatus.AVAILABLE,
    ])
    breaker = SessionCircuitBreaker(failure_threshold=3)
    executor = _executor(provider, breaker=breaker)
    try:
        assert executor.execute(_request("one"), _TIMEOUT).status is ProviderStatus.TIMEOUT
        assert executor.execute(_request("two"), _TIMEOUT).status is ProviderStatus.PROVIDER_ERROR
        assert executor.execute(_request("three"), _TIMEOUT).status is ProviderStatus.TIMEOUT
        assert breaker.is_open("fixture")

        blocked = executor.execute(_request("blocked"), _TIMEOUT)
        assert blocked.status is ProviderStatus.OFFLINE
        assert blocked.attempts == 0
        assert blocked.message == "Provider temporarily unavailable."
        assert blocked.metadata["temporarily_unavailable"] is True
        assert provider.calls == 3

        breaker.reset("fixture")
        recovered = executor.execute(_request("recovered"), _TIMEOUT)
        assert recovered.status is ProviderStatus.AVAILABLE
        assert provider.calls == 4
        assert not breaker.is_open("fixture")
    finally:
        executor.close()


def test_no_match_and_every_other_nonfailure_break_consecutive_failure_sequence() -> None:
    breaker = SessionCircuitBreaker(failure_threshold=2)
    nonfailures = (
        ProviderStatus.AVAILABLE,
        ProviderStatus.NO_MATCH,
        ProviderStatus.OFFLINE,
        ProviderStatus.RATE_LIMIT,
        ProviderStatus.DISABLED,
    )
    for status in nonfailures:
        breaker.record("fixture", ProviderStatus.PROVIDER_ERROR)
        assert breaker.failure_count("fixture") == 1
        breaker.record("fixture", status)
        assert breaker.failure_count("fixture") == 0
        assert not breaker.is_open("fixture")

    provider = ScriptedProvider([
        ProviderStatus.PROVIDER_ERROR,
        ProviderStatus.NO_MATCH,
        ProviderStatus.PROVIDER_ERROR,
    ])
    executor = _executor(provider, threshold=2)
    try:
        assert executor.execute(_request("first"), _TIMEOUT).status is ProviderStatus.PROVIDER_ERROR
        assert executor.execute(_request("valid absence"), _TIMEOUT).status is ProviderStatus.NO_MATCH
        assert executor.execute(_request("third"), _TIMEOUT).status is ProviderStatus.PROVIDER_ERROR
        assert not executor.circuit_breaker.is_open("fixture")
    finally:
        executor.close()


def test_execute_many_is_bounded_independent_and_preserves_input_order() -> None:
    delays = {
        "slow-a": 0.02,
        "fast-b": 0.001,
        "middle-c": 0.01,
        "fast-d": 0.001,
        "slow-e": 0.02,
    }

    def dynamic(request: ProviderRequest) -> tuple[str, float, ProviderStatus]:
        return ("sleep", delays[request.canonical_query], ProviderStatus.AVAILABLE)

    provider = ScriptedProvider(dynamic=dynamic)
    executor = _executor(provider, concurrency=2)
    requests = tuple(_request(query) for query in delays)
    try:
        results = executor.execute_many(requests, _TIMEOUT)
    finally:
        executor.close()

    assert tuple(result.request.canonical_query for result in results) == tuple(delays)
    assert tuple(result.data["query"] for result in results) == tuple(delays)
    assert all(result.status is ProviderStatus.AVAILABLE for result in results)
    assert provider.calls == len(requests)
    assert 1 < provider.max_active <= 2


def test_invalid_provider_results_are_isolated_and_request_identity_is_preserved() -> None:
    other_request = _request("different")
    mismatched = ProviderResult(other_request, ProviderStatus.AVAILABLE, {"bad": True})
    provider = ScriptedProvider([mismatched])
    executor = _executor(provider)
    original = _request("original")
    try:
        result = executor.execute(original, _TIMEOUT)
    finally:
        executor.close()

    assert result.status is ProviderStatus.PROVIDER_ERROR
    assert result.request.cache_key == original.cache_key
    assert result.metadata["reason"] == "request_identity_mismatch"


def test_policy_is_single_immutable_source_for_retry_breaker_and_concurrency_limits() -> None:
    policy = ProviderExecutionPolicy()
    assert policy.max_retries == 1
    assert policy.failure_threshold == 3
    assert policy.max_concurrency == 4
    with pytest.raises(FrozenInstanceError):
        policy.max_retries = 0  # type: ignore[misc]
    with pytest.raises(ValueError, match="between 1 and 16"):
        ProviderExecutionPolicy(max_concurrency=100)
    with pytest.raises(ValueError, match="threshold must match"):
        ProviderExecutor(
            ScriptedProvider([ProviderStatus.AVAILABLE]),
            policy=ProviderExecutionPolicy(failure_threshold=3),
            circuit_breaker=SessionCircuitBreaker(failure_threshold=4),
        )


def test_runtime_has_no_forbidden_dependencies_or_import_time_io() -> None:
    source_path = Path(__file__).parents[1] / "src" / "research" / "provider_runtime.py"
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
    forbidden_prefixes = (
        "requests",
        "httpx",
        "urllib",
        "streamlit",
        "sqlite3",
        "src.state",
        "src.scientific",
        "src.graph_engine",
    )
    assert not any(name.startswith(forbidden_prefixes) for name in imports)
    for forbidden_text in (
        "session_state",
        "read_csv(",
        "read_pickle(",
        "urlopen(",
        "requests.get(",
    ):
        assert forbidden_text not in source

    # Importing/reloading only defines contracts; it creates no worker threads
    # because ThreadPoolExecutor starts workers lazily on first submission.
    import src.research.provider_runtime as runtime_module

    assert runtime_module.__all__ == [
        "EvidenceProvider",
        "ProviderCapability",
        "ProviderExecutionPolicy",
        "ProviderExecutor",
        "ProviderRequest",
        "ProviderResult",
        "ProviderStatus",
        "SessionCircuitBreaker",
    ]
