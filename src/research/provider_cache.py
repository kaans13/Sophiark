"""Layered, cache-first storage for optional Research Context providers.

The cache is deliberately isolated from Sophiark's scientific path.  It never
initializes an :class:`~src.research.evidence_store.EvidenceStore`, never
performs provider I/O on its own, and treats an unavailable L2 SQLite store as
a cache miss.  External execution is only possible through the explicitly
injected :class:`~src.research.provider_runtime.ProviderExecutor` used by
``CacheFirstProviderClient``.
"""

from __future__ import annotations

from collections import OrderedDict
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from enum import Enum
import hashlib
import json
import sqlite3
from threading import RLock
from types import MappingProxyType
from typing import Any, Callable, Mapping

from .config import CacheTTLs, ProviderTimeout
from .evidence_store import EvidenceStore, EvidenceStoreError
from .provider_runtime import (
    EvidenceProvider,
    ProviderExecutor,
    ProviderRequest,
    ProviderResult,
    ProviderStatus,
)


UTCClock = Callable[[], datetime]


def _utc_now() -> datetime:
    return datetime.now(timezone.utc)


def _as_utc(value: datetime, *, field_name: str) -> datetime:
    if not isinstance(value, datetime):
        raise TypeError(f"{field_name} must be a datetime")
    if value.tzinfo is None or value.utcoffset() is None:
        raise ValueError(f"{field_name} must be timezone-aware")
    return value.astimezone(timezone.utc)


def _parse_utc(value: str, *, field_name: str) -> datetime:
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f"{field_name} must be a non-empty UTC timestamp")
    parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    return _as_utc(parsed, field_name=field_name)


def _json_text(value: Any, *, field_name: str) -> str:
    """Return deterministic strict JSON or reject unsafe cache content."""

    try:
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field_name} must contain JSON-safe data") from exc


def _thaw_json(value: Any) -> Any:
    """Convert immutable runtime JSON containers to detached built-ins."""

    if isinstance(value, Mapping):
        return {str(key): _thaw_json(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_thaw_json(item) for item in value]
    return value


def _status_value(status: ProviderStatus) -> str:
    if not isinstance(status, ProviderStatus):
        raise TypeError("ProviderResult.status must be a ProviderStatus")
    return str(status.value)


@dataclass(frozen=True, slots=True)
class ProviderCacheKey:
    """Collision-resistant identity for one provider operation.

    No request parameter is silently substituted for the six required cache
    identity fields.  In particular, a missing taxon remains distinct from any
    concrete taxon and query/provider schema versions are first-class inputs.
    """

    provider: str
    operation: str
    taxon_id: int | None
    canonical_query: str
    provider_schema_version: str
    query_version: str
    request_params_json: str = "{}"

    def __post_init__(self) -> None:
        for name in (
            "provider",
            "operation",
            "canonical_query",
            "provider_schema_version",
            "query_version",
        ):
            value = getattr(self, name)
            if not isinstance(value, str) or not value:
                raise ValueError(f"{name} must be a non-empty string")
        if self.taxon_id is not None and (
            isinstance(self.taxon_id, bool) or not isinstance(self.taxon_id, int)
        ):
            raise TypeError("taxon_id must be an integer or None")
        try:
            params = json.loads(self.request_params_json)
        except (TypeError, ValueError, json.JSONDecodeError) as exc:
            raise ValueError("request_params_json must be valid JSON") from exc
        if not isinstance(params, dict):
            raise ValueError("request_params_json must encode a JSON object")
        canonical_params = _json_text(params, field_name="request params")
        object.__setattr__(self, "request_params_json", canonical_params)

    @classmethod
    def from_request(cls, request: ProviderRequest) -> "ProviderCacheKey":
        return cls(
            provider=request.provider,
            operation=request.operation,
            taxon_id=request.taxon_id,
            canonical_query=request.canonical_query,
            provider_schema_version=request.provider_schema_version,
            query_version=request.query_version,
            request_params_json=_json_text(
                _thaw_json(request.params), field_name="request params"
            ),
        )

    @property
    def identity(self) -> Mapping[str, Any]:
        return MappingProxyType(
            {
                "provider": self.provider,
                "operation": self.operation,
                "taxon_id": self.taxon_id,
                "canonical_query": self.canonical_query,
                "provider_schema_version": self.provider_schema_version,
                "query_version": self.query_version,
                "params": json.loads(self.request_params_json),
            }
        )

    @property
    def request_params_hash(self) -> str:
        return hashlib.sha256(self.request_params_json.encode("utf-8")).hexdigest()

    @property
    def digest(self) -> str:
        canonical = _json_text(dict(self.identity), field_name="cache key")
        return hashlib.sha256(canonical.encode("utf-8")).hexdigest()

    @property
    def cache_key(self) -> str:
        """Alias matching the SQLite schema's ``cache_key`` terminology."""

        return self.digest


class CacheEvidenceClass(str, Enum):
    """Central TTL categories; these are not evidence or scientific scores."""

    ENTITY_RESOLUTION = "entity_resolution"
    BIOLOGICAL_ANNOTATION = "biological_annotation"
    PATHWAY_ANNOTATION = "pathway_annotation"
    LITERATURE = "literature"


_DEFAULT_OPERATION_CLASSES: Mapping[str, CacheEvidenceClass] = MappingProxyType(
    {
        "resolve_entity": CacheEvidenceClass.ENTITY_RESOLUTION,
        "entity_resolution": CacheEvidenceClass.ENTITY_RESOLUTION,
        "pathway": CacheEvidenceClass.PATHWAY_ANNOTATION,
        "pathways": CacheEvidenceClass.PATHWAY_ANNOTATION,
        "pathway_annotation": CacheEvidenceClass.PATHWAY_ANNOTATION,
        "reactome_pathways": CacheEvidenceClass.PATHWAY_ANNOTATION,
        "literature": CacheEvidenceClass.LITERATURE,
        "literature_search": CacheEvidenceClass.LITERATURE,
        "publication_search": CacheEvidenceClass.LITERATURE,
        "search_publications": CacheEvidenceClass.LITERATURE,
    }
)


_TRANSIENT_STATUSES = frozenset(
    {
        ProviderStatus.OFFLINE,
        ProviderStatus.TIMEOUT,
        ProviderStatus.RATE_LIMIT,
        ProviderStatus.PROVIDER_ERROR,
    }
)


@dataclass(frozen=True, slots=True)
class ProviderTTLPolicy:
    """Single configurable policy for positive, negative, and failure TTLs."""

    ttls: CacheTTLs = field(default_factory=CacheTTLs)
    operation_classes: Mapping[str, CacheEvidenceClass | str] = field(
        default_factory=lambda: _DEFAULT_OPERATION_CLASSES
    )
    provider_ttls: Mapping[str, CacheTTLs] = field(
        default_factory=lambda: MappingProxyType({})
    )

    def __post_init__(self) -> None:
        if not isinstance(self.ttls, CacheTTLs):
            raise TypeError("ttls must be a CacheTTLs instance")
        if not isinstance(self.provider_ttls, Mapping):
            raise TypeError("provider_ttls must be a mapping")
        normalized_provider_ttls: dict[str, CacheTTLs] = {}
        for provider_name, provider_ttls in dict(self.provider_ttls).items():
            if not isinstance(provider_name, str) or not provider_name.strip():
                raise ValueError("TTL provider names must be non-empty strings")
            normalized_name = provider_name.strip().casefold()
            if normalized_name in normalized_provider_ttls:
                raise ValueError(
                    f"Duplicate normalized TTL provider name: {normalized_name}"
                )
            if not isinstance(provider_ttls, CacheTTLs):
                raise TypeError("provider_ttls values must be CacheTTLs instances")
            normalized_provider_ttls[normalized_name] = provider_ttls

        normalized: dict[str, CacheEvidenceClass] = {}
        for operation, evidence_class in dict(self.operation_classes).items():
            if not isinstance(operation, str) or not operation.strip():
                raise ValueError("TTL operation names must be non-empty strings")
            normalized[operation.strip().casefold()] = CacheEvidenceClass(evidence_class)
        object.__setattr__(
            self,
            "provider_ttls",
            MappingProxyType(normalized_provider_ttls),
        )
        object.__setattr__(self, "operation_classes", MappingProxyType(normalized))

    def ttls_for(self, provider_name: str) -> CacheTTLs:
        """Return one normalized provider override or the global fallback."""

        if not isinstance(provider_name, str) or not provider_name.strip():
            raise ValueError("provider_name must be a non-empty string")
        return self.provider_ttls.get(provider_name.strip().casefold(), self.ttls)

    def evidence_class_for(self, request: ProviderRequest) -> CacheEvidenceClass:
        return self.operation_classes.get(
            request.operation.strip().casefold(),
            CacheEvidenceClass.BIOLOGICAL_ANNOTATION,
        )

    def ttl_seconds(
        self,
        result: ProviderResult,
        *,
        evidence_class: CacheEvidenceClass | str | None = None,
    ) -> int | None:
        """Return ``None`` only for statuses that must not be cached."""

        selected_ttls = self.ttls_for(result.request.provider)
        if result.status is ProviderStatus.DISABLED:
            return None
        if result.status is ProviderStatus.NO_MATCH:
            return selected_ttls.negative_result_seconds
        if result.status in _TRANSIENT_STATUSES:
            return selected_ttls.transient_failure_seconds
        if result.status is not ProviderStatus.AVAILABLE:
            return None

        selected = (
            self.evidence_class_for(result.request)
            if evidence_class is None
            else CacheEvidenceClass(evidence_class)
        )
        values = {
            CacheEvidenceClass.ENTITY_RESOLUTION: selected_ttls.entity_resolution_seconds,
            CacheEvidenceClass.BIOLOGICAL_ANNOTATION: selected_ttls.biological_annotation_seconds,
            CacheEvidenceClass.PATHWAY_ANNOTATION: selected_ttls.pathway_annotation_seconds,
            CacheEvidenceClass.LITERATURE: selected_ttls.literature_seconds,
        }
        return values[selected]


@dataclass(frozen=True, slots=True)
class ProviderCacheEntry:
    """One immutable result together with its absolute validity interval."""

    key: ProviderCacheKey
    result: ProviderResult
    checked_at: datetime
    expires_at: datetime

    def __post_init__(self) -> None:
        checked_at = _as_utc(self.checked_at, field_name="checked_at")
        expires_at = _as_utc(self.expires_at, field_name="expires_at")
        if expires_at <= checked_at:
            raise ValueError("expires_at must be later than checked_at")
        if ProviderCacheKey.from_request(self.result.request) != self.key:
            raise ValueError("ProviderResult request does not match the cache key")
        object.__setattr__(self, "checked_at", checked_at)
        object.__setattr__(self, "expires_at", expires_at)

    def is_expired(self, now: datetime) -> bool:
        return self.expires_at <= _as_utc(now, field_name="now")


@dataclass(frozen=True, slots=True)
class ProviderCacheStats:
    l1_hits: int
    l2_hits: int
    misses: int
    l2_failures: int
    writes: int
    evictions: int


class ProviderCache:
    """Bounded L1 memory cache backed by an optional research-only L2 store."""

    def __init__(
        self,
        *,
        store: EvidenceStore | None = None,
        max_memory_entries: int = 512,
        ttl_policy: ProviderTTLPolicy | None = None,
        clock: UTCClock = _utc_now,
    ) -> None:
        if isinstance(max_memory_entries, bool) or max_memory_entries <= 0:
            raise ValueError("max_memory_entries must be a positive integer")
        self._store = store
        self._max_memory_entries = int(max_memory_entries)
        self._ttl_policy = ttl_policy or ProviderTTLPolicy()
        self._clock = clock
        self._memory: OrderedDict[str, ProviderCacheEntry] = OrderedDict()
        self._lock = RLock()
        self._l1_hits = 0
        self._l2_hits = 0
        self._misses = 0
        self._l2_failures = 0
        self._writes = 0
        self._evictions = 0

    @property
    def max_memory_entries(self) -> int:
        return self._max_memory_entries

    @property
    def memory_size(self) -> int:
        with self._lock:
            return len(self._memory)

    @property
    def stats(self) -> ProviderCacheStats:
        with self._lock:
            return ProviderCacheStats(
                l1_hits=self._l1_hits,
                l2_hits=self._l2_hits,
                misses=self._misses,
                l2_failures=self._l2_failures,
                writes=self._writes,
                evictions=self._evictions,
            )

    def clear_memory(self) -> None:
        """Clear only process memory; persistent evidence is never deleted."""

        with self._lock:
            self._memory.clear()

    def _remember(self, entry: ProviderCacheEntry) -> None:
        digest = entry.key.digest
        with self._lock:
            self._memory[digest] = entry
            self._memory.move_to_end(digest)
            while len(self._memory) > self._max_memory_entries:
                self._memory.popitem(last=False)
                self._evictions += 1

    def _memory_lookup(
        self, key: ProviderCacheKey, *, now: datetime
    ) -> ProviderCacheEntry | None:
        with self._lock:
            entry = self._memory.get(key.digest)
            if entry is None:
                return None
            if entry.key != key or entry.is_expired(now):
                self._memory.pop(key.digest, None)
                return None
            self._memory.move_to_end(key.digest)
            self._l1_hits += 1
            return entry

    @staticmethod
    def _result_metadata(result: ProviderResult) -> Mapping[str, Any]:
        return {
            "attempts": result.attempts,
            "message": result.message,
            "provider_version": result.provider_version,
            "result_metadata": _thaw_json(result.metadata),
        }

    @staticmethod
    def _validated_entry(entry: ProviderCacheEntry) -> tuple[str, str]:
        payload_json = _json_text(
            _thaw_json(entry.result.data), field_name="provider payload"
        )
        metadata_json = _json_text(
            ProviderCache._result_metadata(entry.result),
            field_name="provider result metadata",
        )
        return payload_json, metadata_json

    def _l2_lookup(
        self,
        request: ProviderRequest,
        key: ProviderCacheKey,
        *,
        now: datetime,
    ) -> ProviderCacheEntry | None:
        if self._store is None:
            return None
        try:
            with self._store.read_connection() as connection:
                row = connection.execute(
                    """
                    SELECT cache_key, provider, operation, taxon_id,
                           canonical_query, provider_schema_version, query_version,
                           request_params_json, request_params_hash,
                           status, payload_json, checked_at, expires_at, metadata_json
                    FROM provider_cache WHERE cache_key = ?
                    """,
                    (key.digest,),
                ).fetchone()
            if row is None:
                return None

            stored_key = ProviderCacheKey(
                provider=row["provider"],
                operation=row["operation"],
                taxon_id=row["taxon_id"],
                canonical_query=row["canonical_query"],
                provider_schema_version=row["provider_schema_version"],
                query_version=row["query_version"],
                request_params_json=row["request_params_json"],
            )
            if (
                stored_key != key
                or row["cache_key"] != key.digest
                or row["request_params_hash"] != key.request_params_hash
            ):
                return None
            checked_at = _parse_utc(row["checked_at"], field_name="checked_at")
            if row["expires_at"] is None:
                return None
            expires_at = _parse_utc(row["expires_at"], field_name="expires_at")
            if expires_at <= now:
                return None

            metadata_value = json.loads(row["metadata_json"])
            payload_value = (
                None if row["payload_json"] is None else json.loads(row["payload_json"])
            )
            if not isinstance(metadata_value, dict):
                return None
            result_metadata = metadata_value.get("result_metadata", {})
            if not isinstance(result_metadata, dict):
                return None
            attempts = metadata_value.get("attempts", 1)
            if isinstance(attempts, bool) or not isinstance(attempts, int) or attempts < 0:
                return None

            result = ProviderResult(
                request=request,
                status=ProviderStatus(row["status"]),
                data=payload_value,
                attempts=attempts,
                message=metadata_value.get("message"),
                provider_version=metadata_value.get("provider_version"),
                metadata=result_metadata,
            )
            return ProviderCacheEntry(
                key=key,
                result=result,
                checked_at=checked_at,
                expires_at=expires_at,
            )
        except (
            EvidenceStoreError,
            OSError,
            sqlite3.Error,
            TypeError,
            ValueError,
            KeyError,
            json.JSONDecodeError,
        ):
            with self._lock:
                self._l2_failures += 1
            return None

    def get(self, request: ProviderRequest) -> ProviderResult | None:
        """Resolve strictly through L1, then L2, otherwise return a miss."""

        key = ProviderCacheKey.from_request(request)
        now = _as_utc(self._clock(), field_name="clock result")
        entry = self._memory_lookup(key, now=now)
        if entry is not None:
            return entry.result

        entry = self._l2_lookup(request, key, now=now)
        if entry is not None:
            self._remember(entry)
            with self._lock:
                self._l2_hits += 1
            return entry.result

        with self._lock:
            self._misses += 1
        return None

    def put(
        self,
        result: ProviderResult,
        *,
        evidence_class: CacheEvidenceClass | str | None = None,
    ) -> bool:
        """Cache a policy-eligible result without leaking L2 failures.

        Invalid JSON is rejected before either cache level is changed.  A valid
        result remains useful in L1 even if the optional SQLite store is down.
        """

        if not isinstance(result, ProviderResult):
            raise TypeError("result must be a ProviderResult")
        ttl_seconds = self._ttl_policy.ttl_seconds(
            result, evidence_class=evidence_class
        )
        if ttl_seconds is None:
            return False

        key = ProviderCacheKey.from_request(result.request)
        checked_at = _as_utc(self._clock(), field_name="clock result")
        entry = ProviderCacheEntry(
            key=key,
            result=result,
            checked_at=checked_at,
            expires_at=checked_at + timedelta(seconds=ttl_seconds),
        )
        payload_json, metadata_json = self._validated_entry(entry)
        self._remember(entry)
        with self._lock:
            self._writes += 1

        if self._store is None:
            return True
        try:
            with self._store.transaction() as connection:
                connection.execute(
                    """
                    INSERT INTO provider_cache(
                        cache_key, provider, operation, taxon_id,
                        canonical_query, provider_schema_version, query_version,
                        request_params_json, request_params_hash,
                        status, payload_json, checked_at, expires_at, metadata_json
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(cache_key) DO UPDATE SET
                        provider = excluded.provider,
                        operation = excluded.operation,
                        taxon_id = excluded.taxon_id,
                        canonical_query = excluded.canonical_query,
                        provider_schema_version = excluded.provider_schema_version,
                        query_version = excluded.query_version,
                        request_params_json = excluded.request_params_json,
                        request_params_hash = excluded.request_params_hash,
                        status = excluded.status,
                        payload_json = excluded.payload_json,
                        checked_at = excluded.checked_at,
                        expires_at = excluded.expires_at,
                        metadata_json = excluded.metadata_json
                    """,
                    (
                        key.digest,
                        key.provider,
                        key.operation,
                        key.taxon_id,
                        key.canonical_query,
                        key.provider_schema_version,
                        key.query_version,
                        key.request_params_json,
                        key.request_params_hash,
                        _status_value(result.status),
                        payload_json,
                        entry.checked_at.isoformat(),
                        entry.expires_at.isoformat(),
                        metadata_json,
                    ),
                )
        except (EvidenceStoreError, OSError, sqlite3.Error):
            with self._lock:
                self._l2_failures += 1
        return True


class CacheFirstProviderClient:
    """L1 -> L2 -> explicitly injected L3 execution orchestration."""

    def __init__(
        self,
        *,
        cache: ProviderCache,
        executor: ProviderExecutor,
    ) -> None:
        if not isinstance(cache, ProviderCache):
            raise TypeError("cache must be a ProviderCache")
        self._cache = cache
        self._executor = executor

    @property
    def provider(self) -> EvidenceProvider:
        """Expose the injected provider without constructing or calling one."""

        return self._executor.provider

    def execute(
        self,
        request: ProviderRequest,
        timeout: ProviderTimeout,
        *,
        evidence_class: CacheEvidenceClass | str | None = None,
    ) -> ProviderResult:
        cached = self._cache.get(request)
        if cached is not None:
            return cached

        result = self._executor.execute(request, timeout)
        if ProviderCacheKey.from_request(result.request) != ProviderCacheKey.from_request(
            request
        ):
            raise ValueError("Provider executor returned a result for a different request")
        try:
            self._cache.put(result, evidence_class=evidence_class)
        except (TypeError, ValueError):
            # Provider data can remain usable even when it is not persistable.
            pass
        return result


__all__ = [
    "CacheEvidenceClass",
    "CacheFirstProviderClient",
    "ProviderCache",
    "ProviderCacheEntry",
    "ProviderCacheKey",
    "ProviderCacheStats",
    "ProviderTTLPolicy",
]
