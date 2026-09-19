"""Fail-soft execution primitives for optional external evidence providers.

This module is deliberately transport agnostic.  It performs no HTTP, file,
session-state, UI, or scientific work; concrete providers supply the I/O in a
later layer.  The runtime only enforces immutable request identity, explicit
status semantics, finite timeouts, conservative retries, a session-local
circuit breaker, and bounded independent concurrency.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from concurrent.futures import ThreadPoolExecutor, TimeoutError as FutureTimeoutError
from dataclasses import dataclass, field, replace
from enum import Enum
import hashlib
import json
import math
from threading import Lock
from types import MappingProxyType
from typing import Any, Iterable, Mapping

from .config import ProviderTimeout


class ProviderStatus(str, Enum):
    """Mutually distinct outcomes for one provider request."""

    AVAILABLE = "AVAILABLE"
    NO_MATCH = "NO_MATCH"
    OFFLINE = "OFFLINE"
    TIMEOUT = "TIMEOUT"
    RATE_LIMIT = "RATE_LIMIT"
    PROVIDER_ERROR = "PROVIDER_ERROR"
    DISABLED = "DISABLED"


_CIRCUIT_FAILURE_STATUSES = frozenset({
    ProviderStatus.TIMEOUT,
    ProviderStatus.PROVIDER_ERROR,
})


def _required_text(value: object, field_name: str, *, folded: bool = False) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{field_name} cannot be empty")
    return text.casefold() if folded else text


def _positive_taxon_id(value: object, field_name: str = "taxon_id") -> int:
    """Reject lossy coercion so requests cannot cross species accidentally."""

    if isinstance(value, bool) or not isinstance(value, int):
        raise TypeError(f"{field_name} must be an integer")
    if value <= 0:
        raise ValueError(f"{field_name} must be a positive NCBI taxonomy identifier")
    return value


def _freeze_json(value: Any, *, path: str = "value") -> Any:
    """Detach and recursively freeze one strictly JSON-compatible value."""

    if value is None or isinstance(value, (str, bool, int)):
        return value
    if isinstance(value, float):
        if not math.isfinite(value):
            raise ValueError(f"{path} cannot contain non-finite numbers")
        return value
    if isinstance(value, Mapping):
        frozen: dict[str, Any] = {}
        for key, item in value.items():
            if not isinstance(key, str):
                raise TypeError(f"{path} keys must be strings")
            frozen[key] = _freeze_json(item, path=f"{path}.{key}")
        return MappingProxyType(frozen)
    if isinstance(value, (list, tuple)):
        return tuple(
            _freeze_json(item, path=f"{path}[{index}]")
            for index, item in enumerate(value)
        )
    raise TypeError(f"{path} contains non-JSON value {type(value).__name__}")


def _thaw_json(value: Any) -> Any:
    """Return detached built-in containers suitable for ``json.dumps``."""

    if isinstance(value, Mapping):
        return {str(key): _thaw_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_thaw_json(item) for item in value]
    return value


def _immutable_json_mapping(
    value: Mapping[str, Any] | None,
    *,
    field_name: str,
) -> Mapping[str, Any]:
    frozen = _freeze_json(value or {}, path=field_name)
    if not isinstance(frozen, Mapping):  # defensive: the input annotation is runtime unchecked
        raise TypeError(f"{field_name} must be a mapping")
    return frozen


@dataclass(frozen=True, slots=True)
class ProviderCapability:
    """One explicit provider operation and its optional taxonomic scope.

    An empty ``supported_taxa`` tuple means that the operation itself is not
    taxon-restricted.  Requests still require a positive explicit taxon ID.
    """

    operation: str
    supported_taxa: tuple[int, ...] = ()
    description: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "operation", _required_text(self.operation, "operation", folded=True))
        taxa = tuple(sorted({
            _positive_taxon_id(taxon_id, "supported_taxa item")
            for taxon_id in self.supported_taxa
        }))
        object.__setattr__(self, "supported_taxa", taxa)
        if self.description is not None:
            object.__setattr__(self, "description", str(self.description).strip() or None)

    def supports(self, request: "ProviderRequest") -> bool:
        return (
            request.operation == self.operation
            and (not self.supported_taxa or request.taxon_id in self.supported_taxa)
        )


@dataclass(frozen=True, slots=True)
class ProviderRequest:
    """Canonical, immutable identity for one provider lookup."""

    provider: str
    operation: str
    taxon_id: int
    canonical_query: str
    provider_schema_version: str = "1"
    query_version: str = "1"
    params: Mapping[str, Any] = field(default_factory=lambda: MappingProxyType({}))

    def __post_init__(self) -> None:
        object.__setattr__(self, "provider", _required_text(self.provider, "provider", folded=True))
        object.__setattr__(self, "operation", _required_text(self.operation, "operation", folded=True))
        object.__setattr__(self, "taxon_id", _positive_taxon_id(self.taxon_id))
        object.__setattr__(self, "canonical_query", _required_text(self.canonical_query, "canonical_query"))
        object.__setattr__(
            self,
            "provider_schema_version",
            _required_text(self.provider_schema_version, "provider_schema_version"),
        )
        object.__setattr__(self, "query_version", _required_text(self.query_version, "query_version"))
        object.__setattr__(
            self,
            "params",
            _immutable_json_mapping(self.params, field_name="params"),
        )

    @property
    def identity(self) -> Mapping[str, Any]:
        """The exact cache/provenance identity, including immutable params."""

        return MappingProxyType({
            "provider": self.provider,
            "operation": self.operation,
            "taxon_id": self.taxon_id,
            "canonical_query": self.canonical_query,
            "provider_schema_version": self.provider_schema_version,
            "query_version": self.query_version,
            "params": self.params,
        })

    @property
    def cache_key(self) -> str:
        encoded = json.dumps(
            self.to_dict(),
            ensure_ascii=False,
            allow_nan=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
        return hashlib.sha256(encoded).hexdigest()

    def to_dict(self) -> dict[str, Any]:
        return {
            "provider": self.provider,
            "operation": self.operation,
            "taxon_id": self.taxon_id,
            "canonical_query": self.canonical_query,
            "provider_schema_version": self.provider_schema_version,
            "query_version": self.query_version,
            "params": _thaw_json(self.params),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ProviderRequest":
        return cls(
            provider=value["provider"],
            operation=value["operation"],
            taxon_id=value["taxon_id"],
            canonical_query=value["canonical_query"],
            provider_schema_version=value.get("provider_schema_version", "1"),
            query_version=value.get("query_version", "1"),
            params=value.get("params", {}),
        )


@dataclass(frozen=True, slots=True)
class ProviderResult:
    """Immutable provider outcome bound to the request's exact identity."""

    request: ProviderRequest
    status: ProviderStatus
    data: Any = None
    attempts: int = 1
    message: str | None = None
    provider_version: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=lambda: MappingProxyType({}))

    def __post_init__(self) -> None:
        if not isinstance(self.request, ProviderRequest):
            raise TypeError("request must be a ProviderRequest")
        object.__setattr__(self, "status", ProviderStatus(self.status))
        object.__setattr__(self, "data", _freeze_json(self.data, path="data"))
        object.__setattr__(self, "attempts", int(self.attempts))
        object.__setattr__(
            self,
            "metadata",
            _immutable_json_mapping(self.metadata, field_name="metadata"),
        )
        if self.attempts < 0:
            raise ValueError("attempts cannot be negative")
        if self.message is not None:
            object.__setattr__(self, "message", str(self.message).strip() or None)
        if self.provider_version is not None:
            object.__setattr__(
                self,
                "provider_version",
                str(self.provider_version).strip() or None,
            )

    @property
    def provider(self) -> str:
        return self.request.provider

    @property
    def operation(self) -> str:
        return self.request.operation

    @property
    def request_identity(self) -> Mapping[str, Any]:
        return self.request.identity

    @property
    def available(self) -> bool:
        return self.status is ProviderStatus.AVAILABLE

    def to_dict(self) -> dict[str, Any]:
        """Return a detached representation accepted by strict JSON encoders."""

        return {
            "request": self.request.to_dict(),
            "request_identity": self.request.to_dict(),
            "status": self.status.value,
            "data": _thaw_json(self.data),
            "attempts": self.attempts,
            "message": self.message,
            "provider_version": self.provider_version,
            "metadata": _thaw_json(self.metadata),
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ProviderResult":
        request_value = value.get("request", value.get("request_identity"))
        if not isinstance(request_value, Mapping):
            raise TypeError("serialized ProviderResult requires request identity")
        identity_value = value.get("request_identity")
        if identity_value is not None:
            if not isinstance(identity_value, Mapping):
                raise TypeError("request_identity must be a mapping")
            request_key = ProviderRequest.from_dict(request_value).cache_key
            identity_key = ProviderRequest.from_dict(identity_value).cache_key
            if request_key != identity_key:
                raise ValueError("serialized request and request_identity do not match")
        return cls(
            request=ProviderRequest.from_dict(request_value),
            status=value["status"],
            data=value.get("data"),
            attempts=value.get("attempts", 1),
            message=value.get("message"),
            provider_version=value.get("provider_version"),
            metadata=value.get("metadata", {}),
        )


class EvidenceProvider(ABC):
    """Base contract for one optional external evidence source.

    Concrete providers implement transport-specific ``execute``.  Callers use
    :class:`ProviderExecutor`, which checks provider identity, enablement, and
    capability before the transport method can run.
    """

    def __init__(
        self,
        provider_name: str,
        capabilities: Iterable[ProviderCapability],
        *,
        provider_version: str | None = None,
        enabled: bool = True,
    ) -> None:
        self._provider_name = _required_text(provider_name, "provider_name", folded=True)
        frozen_capabilities = tuple(capabilities)
        if not frozen_capabilities:
            raise ValueError("at least one provider capability is required")
        if not all(isinstance(item, ProviderCapability) for item in frozen_capabilities):
            raise TypeError("capabilities must contain ProviderCapability values")
        capability_keys = {
            (item.operation, item.supported_taxa)
            for item in frozen_capabilities
        }
        if len(capability_keys) != len(frozen_capabilities):
            raise ValueError("duplicate provider capabilities are not allowed")
        self._capabilities = frozen_capabilities
        self._provider_version = (
            str(provider_version).strip() or None
            if provider_version is not None
            else None
        )
        self._enabled = bool(enabled)

    @property
    def provider_name(self) -> str:
        return self._provider_name

    @property
    def provider_version(self) -> str | None:
        return self._provider_version

    @property
    def capabilities(self) -> tuple[ProviderCapability, ...]:
        return self._capabilities

    @property
    def enabled(self) -> bool:
        return self._enabled

    def supports(self, request: ProviderRequest) -> bool:
        return (
            self.enabled
            and request.provider == self.provider_name
            and any(capability.supports(request) for capability in self.capabilities)
        )

    def healthcheck(self) -> ProviderStatus:
        """Return local availability without requiring a network probe."""

        return ProviderStatus.AVAILABLE if self.enabled else ProviderStatus.DISABLED

    @abstractmethod
    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        """Perform exactly one transport attempt for a supported request."""

        raise NotImplementedError


@dataclass(frozen=True, slots=True)
class ProviderExecutionPolicy:
    """Conservative, process-local limits for optional provider execution."""

    max_retries: int = 1
    failure_threshold: int = 3
    max_concurrency: int = 4

    def __post_init__(self) -> None:
        for field_name in ("max_retries", "failure_threshold", "max_concurrency"):
            value = getattr(self, field_name)
            if isinstance(value, bool) or not isinstance(value, int):
                raise TypeError(f"{field_name} must be an integer")
        if not 0 <= self.max_retries <= 1:
            raise ValueError("max_retries must be zero or one")
        if self.failure_threshold < 1:
            raise ValueError("failure_threshold must be positive")
        if not 1 <= self.max_concurrency <= 16:
            raise ValueError("max_concurrency must be between 1 and 16")


class SessionCircuitBreaker:
    """Thread-safe, session-lifetime breaker keyed by provider name."""

    def __init__(self, failure_threshold: int = 3) -> None:
        if isinstance(failure_threshold, bool) or not isinstance(failure_threshold, int):
            raise TypeError("failure_threshold must be an integer")
        if failure_threshold < 1:
            raise ValueError("failure_threshold must be positive")
        self._failure_threshold = failure_threshold
        self._consecutive_failures: dict[str, int] = {}
        self._open: set[str] = set()
        self._lock = Lock()

    @property
    def failure_threshold(self) -> int:
        return self._failure_threshold

    def is_open(self, provider_name: str) -> bool:
        key = _required_text(provider_name, "provider_name", folded=True)
        with self._lock:
            return key in self._open

    def failure_count(self, provider_name: str) -> int:
        key = _required_text(provider_name, "provider_name", folded=True)
        with self._lock:
            return self._consecutive_failures.get(key, 0)

    def record(self, provider_name: str, status: ProviderStatus) -> None:
        """Record one completed attempt; only two statuses count as failures."""

        key = _required_text(provider_name, "provider_name", folded=True)
        normalized_status = ProviderStatus(status)
        with self._lock:
            if normalized_status in _CIRCUIT_FAILURE_STATUSES:
                count = self._consecutive_failures.get(key, 0) + 1
                self._consecutive_failures[key] = count
                if count >= self.failure_threshold:
                    self._open.add(key)
                return
            # NO_MATCH is a valid provider answer, and every other non-failure
            # status breaks a sequence of consecutive transport/server errors.
            self._consecutive_failures.pop(key, None)
            self._open.discard(key)

    def reset(self, provider_name: str | None = None) -> None:
        """Close one provider circuit, or all circuits for a new session."""

        with self._lock:
            if provider_name is None:
                self._consecutive_failures.clear()
                self._open.clear()
                return
            key = _required_text(provider_name, "provider_name", folded=True)
            self._consecutive_failures.pop(key, None)
            self._open.discard(key)

    def snapshot(self) -> Mapping[str, Mapping[str, Any]]:
        with self._lock:
            names = set(self._consecutive_failures) | set(self._open)
            return MappingProxyType({
                name: MappingProxyType({
                    "consecutive_failures": self._consecutive_failures.get(name, 0),
                    "open": name in self._open,
                })
                for name in sorted(names)
            })


class ProviderExecutor:
    """Apply timeout, retry, circuit-breaker, and concurrency policy safely."""

    def __init__(
        self,
        provider: EvidenceProvider,
        *,
        policy: ProviderExecutionPolicy | None = None,
        circuit_breaker: SessionCircuitBreaker | None = None,
    ) -> None:
        if not isinstance(provider, EvidenceProvider):
            raise TypeError("provider must implement EvidenceProvider")
        self.provider = provider
        self.policy = policy or ProviderExecutionPolicy()
        if (
            circuit_breaker is not None
            and circuit_breaker.failure_threshold != self.policy.failure_threshold
        ):
            raise ValueError(
                "circuit_breaker threshold must match ProviderExecutionPolicy"
            )
        self.circuit_breaker = circuit_breaker or SessionCircuitBreaker(
            self.policy.failure_threshold
        )
        # ThreadPoolExecutor creates threads lazily.  No work or I/O occurs at
        # construction/import time, while active provider calls remain bounded.
        self._attempt_pool = ThreadPoolExecutor(
            max_workers=self.policy.max_concurrency,
            thread_name_prefix=f"research-{self.provider.provider_name}",
        )
        self._closed = False
        self._state_lock = Lock()

    def _result(
        self,
        request: ProviderRequest,
        status: ProviderStatus,
        *,
        attempts: int,
        message: str,
        metadata: Mapping[str, Any] | None = None,
    ) -> ProviderResult:
        return ProviderResult(
            request=request,
            status=status,
            attempts=attempts,
            message=message,
            provider_version=self.provider.provider_version,
            metadata=metadata or {},
        )

    def _preflight(self, request: ProviderRequest) -> ProviderResult | None:
        if request.provider != self.provider.provider_name:
            return self._result(
                request,
                ProviderStatus.DISABLED,
                attempts=0,
                message="Request is assigned to a different provider.",
                metadata={"reason": "provider_mismatch"},
            )
        if not self.provider.enabled:
            return self._result(
                request,
                ProviderStatus.DISABLED,
                attempts=0,
                message="Provider is disabled.",
                metadata={"reason": "provider_disabled"},
            )
        if not self.provider.supports(request):
            return self._result(
                request,
                ProviderStatus.DISABLED,
                attempts=0,
                message="Provider does not support this operation or taxon.",
                metadata={"reason": "unsupported_capability"},
            )
        if self.circuit_breaker.is_open(self.provider.provider_name):
            return self._result(
                request,
                ProviderStatus.OFFLINE,
                attempts=0,
                message="Provider temporarily unavailable.",
                metadata={"reason": "circuit_open", "temporarily_unavailable": True},
            )
        return None

    def _one_attempt(
        self,
        request: ProviderRequest,
        timeout: ProviderTimeout,
    ) -> ProviderResult:
        try:
            future = self._attempt_pool.submit(self.provider.execute, request, timeout)
        except RuntimeError as exc:
            with self._state_lock:
                closed = self._closed
            return self._result(
                request,
                ProviderStatus.DISABLED if closed else ProviderStatus.PROVIDER_ERROR,
                attempts=0,
                message=str(exc) or "Provider executor could not schedule the request.",
                metadata={
                    "reason": "executor_closed" if closed else "executor_submit_failed",
                    "exception_type": type(exc).__name__,
                },
            )
        try:
            result = future.result(timeout=timeout.total_seconds)
        except FutureTimeoutError as exc:
            # ``concurrent.futures.TimeoutError`` aliases built-in
            # ``TimeoutError``.  A completed future therefore means the
            # provider itself raised a timeout; an incomplete one exceeded our
            # centrally enforced deadline.  Both remain the distinct TIMEOUT
            # status, while provenance explains which boundary detected it.
            if future.done():
                return self._result(
                    request,
                    ProviderStatus.TIMEOUT,
                    attempts=1,
                    message=str(exc) or "Provider timed out.",
                    metadata={"exception_type": type(exc).__name__},
                )
            future.cancel()
            return self._result(
                request,
                ProviderStatus.TIMEOUT,
                attempts=1,
                message=f"Provider exceeded the {timeout.total_seconds:g}s total timeout.",
                metadata={"reason": "deadline_exceeded"},
            )
        except ConnectionError as exc:
            return self._result(
                request,
                ProviderStatus.OFFLINE,
                attempts=1,
                message=str(exc) or "Provider is offline.",
                metadata={"exception_type": type(exc).__name__},
            )
        except Exception as exc:  # provider isolation is intentionally fail-soft
            return self._result(
                request,
                ProviderStatus.PROVIDER_ERROR,
                attempts=1,
                message=str(exc) or "Provider execution failed.",
                metadata={"exception_type": type(exc).__name__},
            )

        if not isinstance(result, ProviderResult):
            return self._result(
                request,
                ProviderStatus.PROVIDER_ERROR,
                attempts=1,
                message="Provider returned an invalid result type.",
                metadata={"reason": "invalid_result_type"},
            )
        if result.request.cache_key != request.cache_key:
            return self._result(
                request,
                ProviderStatus.PROVIDER_ERROR,
                attempts=1,
                message="Provider result request identity does not match the request.",
                metadata={"reason": "request_identity_mismatch"},
            )
        return replace(
            result,
            attempts=1,
            provider_version=result.provider_version or self.provider.provider_version,
        )

    def execute(
        self,
        request: ProviderRequest,
        timeout: ProviderTimeout,
    ) -> ProviderResult:
        """Execute one request with at most one retry and no propagated failure."""

        if not isinstance(request, ProviderRequest):
            raise TypeError("request must be a ProviderRequest")
        if not isinstance(timeout, ProviderTimeout):
            raise TypeError("timeout must be a centrally configured ProviderTimeout")
        with self._state_lock:
            closed = self._closed
        if closed:
            return self._result(
                request,
                ProviderStatus.DISABLED,
                attempts=0,
                message="Provider executor is closed.",
                metadata={"reason": "executor_closed"},
            )

        preflight = self._preflight(request)
        if preflight is not None:
            return preflight

        attempts = 0
        last_result: ProviderResult | None = None
        for _ in range(self.policy.max_retries + 1):
            if self.circuit_breaker.is_open(self.provider.provider_name):
                if last_result is not None:
                    return replace(last_result, attempts=attempts)
                return self._preflight(request) or self._result(
                    request,
                    ProviderStatus.OFFLINE,
                    attempts=attempts,
                    message="Provider temporarily unavailable.",
                )

            attempts += 1
            result = self._one_attempt(request, timeout)
            result = replace(result, attempts=attempts)
            self.circuit_breaker.record(self.provider.provider_name, result.status)
            last_result = result
            if result.status not in _CIRCUIT_FAILURE_STATUSES:
                return result
            if attempts > self.policy.max_retries:
                return result

        # The range is finite and always executes at least once; this branch is
        # defensive and keeps the method total if policy implementation changes.
        return last_result or self._result(
            request,
            ProviderStatus.PROVIDER_ERROR,
            attempts=attempts,
            message="Provider execution produced no result.",
        )

    def execute_many(
        self,
        requests: Iterable[ProviderRequest],
        timeout: ProviderTimeout,
    ) -> tuple[ProviderResult, ...]:
        """Execute independent requests with bounded concurrency and stable order."""

        if not isinstance(timeout, ProviderTimeout):
            raise TypeError("timeout must be a centrally configured ProviderTimeout")
        frozen_requests = tuple(requests)
        if not all(isinstance(request, ProviderRequest) for request in frozen_requests):
            raise TypeError("requests must contain only ProviderRequest values")
        if not frozen_requests:
            return ()

        with ThreadPoolExecutor(
            max_workers=self.policy.max_concurrency,
            thread_name_prefix=f"research-batch-{self.provider.provider_name}",
        ) as batch_pool:
            futures = [
                batch_pool.submit(self.execute, request, timeout)
                for request in frozen_requests
            ]
            # Reading futures in source order preserves the planner's ordering;
            # work still proceeds independently in the bounded worker pool.
            return tuple(future.result() for future in futures)

    def close(self, *, wait: bool = False) -> None:
        """Release executor workers without reintroducing an unbounded wait.

        Python cannot forcibly stop a provider call that is already running in
        a thread.  The default therefore cancels queued work and returns
        immediately; callers performing controlled process shutdown may opt
        into ``wait=True`` explicitly.
        """

        with self._state_lock:
            if self._closed:
                return
            self._closed = True
        self._attempt_pool.shutdown(wait=wait, cancel_futures=True)

    def __enter__(self) -> "ProviderExecutor":
        return self

    def __exit__(self, exc_type, exc, traceback) -> None:
        del exc_type, exc, traceback
        self.close(wait=False)


__all__ = [
    "EvidenceProvider",
    "ProviderCapability",
    "ProviderExecutionPolicy",
    "ProviderExecutor",
    "ProviderRequest",
    "ProviderResult",
    "ProviderStatus",
    "SessionCircuitBreaker",
]
