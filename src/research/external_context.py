"""Immutable external-context snapshots with explicit refresh and offline restore.

Snapshots capture provider outputs downstream of one immutable simulation
result.  Saving a refresh always inserts a new row; existing snapshots are
never rewritten and restoring performs no provider or network operation.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
import sqlite3
from types import MappingProxyType
from typing import Any, Iterable, Mapping
from uuid import uuid4

from .evidence_store import EvidenceStore, EvidenceStoreError
from .external_biology import ExternalBiologyLoad
from .literature import LiteratureLoad
from .models import SimulationResultSnapshot, immutable_mapping
from .provider_runtime import ProviderResult, ProviderStatus
from .query_planner import QueryPlan


EXTERNAL_CONTEXT_SCHEMA_VERSION = "external-context-snapshot-v1"
MAX_EXTERNAL_SNAPSHOT_BYTES = 16 * 1024 * 1024


class ExternalSnapshotStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    PARTIAL = "PARTIAL"
    UNAVAILABLE = "UNAVAILABLE"
    EMPTY = "EMPTY"


class SnapshotRestoreStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    NOT_FOUND = "NOT_FOUND"
    UNAVAILABLE = "UNAVAILABLE"
    INVALID = "INVALID"


def _plain(value: Any) -> Any:
    if isinstance(value, Mapping):
        return {str(key): _plain(item) for key, item in value.items()}
    if isinstance(value, (tuple, list)):
        return [_plain(item) for item in value]
    if isinstance(value, Enum):
        return value.value
    return value


def _json(value: Any) -> str:
    return json.dumps(
        _plain(value),
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _digest(value: Any) -> str:
    return hashlib.sha256(_json(value).encode("utf-8")).hexdigest()


def _utc(value: str | datetime | None) -> str:
    if value is None:
        timestamp = datetime.now(timezone.utc)
    elif isinstance(value, datetime):
        timestamp = value
    else:
        text = str(value).strip()
        if text.endswith("Z"):
            text = text[:-1] + "+00:00"
        timestamp = datetime.fromisoformat(text)
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    return timestamp.astimezone(timezone.utc).isoformat()


def _required(value: Any, field_name: str) -> str:
    text = str(value).strip()
    if not text:
        raise ValueError(f"{field_name} cannot be empty")
    return text


def _status(results: tuple[ProviderResult, ...]) -> ExternalSnapshotStatus:
    if not results:
        return ExternalSnapshotStatus.EMPTY
    usable = {
        ProviderStatus.AVAILABLE,
        ProviderStatus.NO_MATCH,
    }
    usable_count = sum(result.status in usable for result in results)
    if usable_count == len(results):
        return ExternalSnapshotStatus.AVAILABLE
    if usable_count:
        return ExternalSnapshotStatus.PARTIAL
    return ExternalSnapshotStatus.UNAVAILABLE


def _publication(value: Any) -> dict[str, Any]:
    return {
        "publication_id": value.publication_id,
        "pmid": value.pmid,
        "pmcid": value.pmcid,
        "doi": value.doi,
        "title": value.title,
        "year": value.year,
        "journal": value.journal,
        "authors_summary": value.authors_summary,
        "abstract": value.abstract,
        "source": value.source,
        "retrieved_at": value.retrieved_at,
        "provenance": dict(value.provenance),
        "metadata": dict(value.metadata),
    }


def _literature_match(value: Any) -> dict[str, Any]:
    return {
        "publication": _publication(value.publication),
        "query": value.query.to_dict(),
        "relevant_passages": [
            {
                "publication_id": passage.publication_id,
                "text": passage.text,
                "source_section": passage.source_section,
                "matched_entities": list(passage.matched_entities),
                "matched_query_terms": list(passage.matched_query_terms),
                "rank": passage.rank,
                "basis": passage.basis,
            }
            for passage in value.relevant_passages
        ],
        "match_label": value.match_label,
    }


@dataclass(frozen=True, slots=True)
class ExternalContextSnapshot:
    context_snapshot_id: str
    simulation_id: str
    simulation_snapshot_id: str
    created_at: str
    query_version: str
    providers: tuple[str, ...]
    provider_versions: Mapping[str, str | None]
    provider_statuses: Mapping[str, tuple[str, ...]]
    literature_cutoff: str | None
    status: ExternalSnapshotStatus
    refresh_token: str
    refresh_of: str | None
    payload: Mapping[str, Any]
    schema_version: str = EXTERNAL_CONTEXT_SCHEMA_VERSION

    def __post_init__(self) -> None:
        for name in (
            "context_snapshot_id", "simulation_id", "simulation_snapshot_id",
            "query_version", "refresh_token", "schema_version",
        ):
            object.__setattr__(self, name, _required(getattr(self, name), name))
        object.__setattr__(self, "created_at", _utc(self.created_at))
        object.__setattr__(self, "providers", tuple(dict.fromkeys(
            str(item).strip().casefold() for item in self.providers if str(item).strip()
        )))
        object.__setattr__(
            self,
            "provider_versions",
            MappingProxyType({
                str(key).casefold(): (str(value).strip() or None if value is not None else None)
                for key, value in self.provider_versions.items()
            }),
        )
        object.__setattr__(
            self,
            "provider_statuses",
            MappingProxyType({
                str(key).casefold(): tuple(str(item) for item in values)
                for key, values in self.provider_statuses.items()
            }),
        )
        object.__setattr__(self, "status", ExternalSnapshotStatus(self.status))
        object.__setattr__(self, "payload", immutable_mapping(self.payload))
        if self.refresh_of is not None:
            object.__setattr__(self, "refresh_of", _required(self.refresh_of, "refresh_of"))
            if self.refresh_of == self.context_snapshot_id:
                raise ValueError("a refresh cannot reference itself")
        if self.literature_cutoff is not None:
            object.__setattr__(self, "literature_cutoff", _utc(self.literature_cutoff))
        expected = _snapshot_id(
            simulation_id=self.simulation_id,
            simulation_snapshot_id=self.simulation_snapshot_id,
            created_at=self.created_at,
            query_version=self.query_version,
            refresh_token=self.refresh_token,
            refresh_of=self.refresh_of,
            providers=self.providers,
            provider_versions=self.provider_versions,
            provider_statuses=self.provider_statuses,
            literature_cutoff=self.literature_cutoff,
            status=self.status,
            payload=self.payload,
            schema_version=self.schema_version,
        )
        if self.context_snapshot_id != expected:
            raise ValueError("external context snapshot integrity check failed")

    @property
    def payload_digest(self) -> str:
        return _digest(self.payload)

    def provider_results(self) -> tuple[ProviderResult, ...]:
        values = self.payload.get("provider_results", ())
        return tuple(
            ProviderResult.from_dict(value)
            for value in values
            if isinstance(value, Mapping)
        )

    def to_dict(self) -> dict[str, Any]:
        return {
            "context_snapshot_id": self.context_snapshot_id,
            "simulation_id": self.simulation_id,
            "simulation_snapshot_id": self.simulation_snapshot_id,
            "created_at": self.created_at,
            "query_version": self.query_version,
            "providers": list(self.providers),
            "provider_versions": dict(self.provider_versions),
            "provider_statuses": _plain(self.provider_statuses),
            "literature_cutoff": self.literature_cutoff,
            "status": self.status.value,
            "refresh_token": self.refresh_token,
            "refresh_of": self.refresh_of,
            "payload": _plain(self.payload),
            "schema_version": self.schema_version,
        }

    @classmethod
    def from_dict(cls, value: Mapping[str, Any]) -> "ExternalContextSnapshot":
        return cls(
            context_snapshot_id=value["context_snapshot_id"],
            simulation_id=value["simulation_id"],
            simulation_snapshot_id=value["simulation_snapshot_id"],
            created_at=value["created_at"],
            query_version=value["query_version"],
            providers=tuple(value.get("providers", ())),
            provider_versions=value.get("provider_versions", {}),
            provider_statuses=value.get("provider_statuses", {}),
            literature_cutoff=value.get("literature_cutoff"),
            status=value["status"],
            refresh_token=value["refresh_token"],
            refresh_of=value.get("refresh_of"),
            payload=value.get("payload", {}),
            schema_version=value.get("schema_version", EXTERNAL_CONTEXT_SCHEMA_VERSION),
        )


def _snapshot_id(**identity: Any) -> str:
    return "external-context-" + _digest(identity)[:32]


def create_external_context_snapshot(
    simulation_snapshot: SimulationResultSnapshot,
    query_plan: QueryPlan,
    *,
    biology_load: ExternalBiologyLoad | None = None,
    literature_load: LiteratureLoad | None = None,
    expected_simulation_snapshot_id: str,
    created_at: str | datetime | None = None,
    refresh_of: ExternalContextSnapshot | None = None,
    refresh_token: str | None = None,
) -> ExternalContextSnapshot:
    """Capture already-retrieved results; this function performs no I/O."""

    if not isinstance(simulation_snapshot, SimulationResultSnapshot):
        raise TypeError("simulation_snapshot must be SimulationResultSnapshot")
    if not isinstance(query_plan, QueryPlan):
        raise TypeError("query_plan must be QueryPlan")
    if str(expected_simulation_snapshot_id).strip() != simulation_snapshot.snapshot_id:
        raise ValueError("stale simulation snapshot; external capture rejected")
    if (
        query_plan.snapshot_id != simulation_snapshot.snapshot_id
        or query_plan.simulation_id != simulation_snapshot.simulation_id
        or query_plan.taxon_id != simulation_snapshot.taxon_id
    ):
        raise ValueError("query plan does not belong to the simulation snapshot")
    if biology_load is not None:
        if biology_load.plan.query_plan_id != query_plan.plan_id:
            raise ValueError("external biology load belongs to a different query plan")
        if biology_load.plan.snapshot_id != simulation_snapshot.snapshot_id:
            raise ValueError("external biology load is stale")
    if literature_load is not None:
        if literature_load.plan.query_plan_id != query_plan.plan_id:
            raise ValueError("literature load belongs to a different query plan")
        if literature_load.plan.snapshot_id != simulation_snapshot.snapshot_id:
            raise ValueError("literature load is stale")
    if refresh_of is not None:
        if (
            refresh_of.simulation_id != simulation_snapshot.simulation_id
            or refresh_of.simulation_snapshot_id != simulation_snapshot.snapshot_id
        ):
            raise ValueError("refresh source belongs to a different simulation snapshot")

    biology_results = biology_load.results if biology_load is not None else ()
    literature_results = literature_load.results if literature_load is not None else ()
    results = tuple((*biology_results, *literature_results))
    providers = tuple(dict.fromkeys(result.provider for result in results))
    versions: dict[str, str | None] = {}
    statuses: dict[str, list[str]] = {}
    for result in results:
        versions[result.provider] = result.provider_version or versions.get(result.provider)
        statuses.setdefault(result.provider, []).append(result.status.value)
    matches = (
        tuple(_literature_match(match) for match in literature_load.matches)
        if literature_load is not None
        else ()
    )
    cutoff_values = [
        str(match["publication"]["retrieved_at"])
        for match in matches
        if match["publication"].get("retrieved_at")
    ]
    literature_cutoff = max(cutoff_values) if cutoff_values else None
    payload = {
        "simulation_scientific_fingerprint": simulation_snapshot.scientific_fingerprint,
        "simulation_snapshot_fingerprint": simulation_snapshot.snapshot_fingerprint,
        "query_plan_id": query_plan.plan_id,
        "query_level": query_plan.level.value,
        "query_budget": {
            "max_provider_requests": query_plan.budget.max_provider_requests,
            "max_literature_queries": query_plan.budget.max_literature_queries,
            "max_publications_per_query": query_plan.budget.max_publications_per_query,
        },
        "provider_results": [result.to_dict() for result in results],
        "literature_matches": list(matches),
        "selection_note": query_plan.selection_note,
    }
    payload_size = len(_json(payload).encode("utf-8"))
    if payload_size > MAX_EXTERNAL_SNAPSHOT_BYTES:
        raise ValueError(
            f"external context snapshot exceeds {MAX_EXTERNAL_SNAPSHOT_BYTES} byte safety limit"
        )
    timestamp = _utc(created_at)
    token = _required(refresh_token or uuid4().hex, "refresh_token")
    fields = {
        "simulation_id": simulation_snapshot.simulation_id,
        "simulation_snapshot_id": simulation_snapshot.snapshot_id,
        "created_at": timestamp,
        "query_version": query_plan.query_version,
        "refresh_token": token,
        "refresh_of": refresh_of.context_snapshot_id if refresh_of is not None else None,
        "providers": providers,
        "provider_versions": versions,
        "provider_statuses": {key: tuple(value) for key, value in statuses.items()},
        "literature_cutoff": literature_cutoff,
        "status": _status(results),
        "payload": payload,
        "schema_version": EXTERNAL_CONTEXT_SCHEMA_VERSION,
    }
    return ExternalContextSnapshot(
        context_snapshot_id=_snapshot_id(**fields),
        **fields,
    )


def refresh_external_context_snapshot(
    previous: ExternalContextSnapshot,
    simulation_snapshot: SimulationResultSnapshot,
    query_plan: QueryPlan,
    *,
    biology_load: ExternalBiologyLoad | None = None,
    literature_load: LiteratureLoad | None = None,
    expected_simulation_snapshot_id: str,
    created_at: str | datetime | None = None,
    refresh_token: str | None = None,
) -> ExternalContextSnapshot:
    """Explicit refresh API; never modifies or overwrites ``previous``."""

    if not isinstance(previous, ExternalContextSnapshot):
        raise TypeError("previous must be ExternalContextSnapshot")
    return create_external_context_snapshot(
        simulation_snapshot,
        query_plan,
        biology_load=biology_load,
        literature_load=literature_load,
        expected_simulation_snapshot_id=expected_simulation_snapshot_id,
        created_at=created_at,
        refresh_of=previous,
        refresh_token=refresh_token,
    )


@dataclass(frozen=True, slots=True)
class SnapshotSaveResult:
    saved: bool
    context_snapshot_id: str
    message: str | None = None


@dataclass(frozen=True, slots=True)
class SnapshotRestoreResult:
    status: SnapshotRestoreStatus
    snapshot: ExternalContextSnapshot | None = None
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "status", SnapshotRestoreStatus(self.status))
        if self.status is SnapshotRestoreStatus.AVAILABLE and self.snapshot is None:
            raise ValueError("AVAILABLE restore requires a snapshot")


class ExternalContextSnapshotStore:
    """Append-only persistence and provider-free offline restore."""

    def __init__(self, store: EvidenceStore) -> None:
        if not isinstance(store, EvidenceStore):
            raise TypeError("store must be EvidenceStore")
        self._store = store

    @staticmethod
    def _simulation_metadata(snapshot: SimulationResultSnapshot) -> str:
        return _json({
            "snapshot_id": snapshot.snapshot_id,
            "scientific_fingerprint": snapshot.scientific_fingerprint,
            "snapshot_fingerprint": snapshot.snapshot_fingerprint,
            "tested_count": snapshot.tested_count,
            "returned_count": snapshot.returned_count,
            "top_n": snapshot.top_n,
        })

    @staticmethod
    def _persist_literature_matches(
        connection: sqlite3.Connection,
        snapshot: ExternalContextSnapshot,
    ) -> None:
        matches = snapshot.payload.get("literature_matches", ())
        if not matches:
            return
        connection.execute(
            """
            INSERT INTO sources(
                source_key, source_name, provider, source_identifier,
                source_version, retrieved_at, attribution, license_note,
                metadata_json, created_at
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            ON CONFLICT(source_key) DO UPDATE SET
                source_version=excluded.source_version,
                retrieved_at=excluded.retrieved_at,
                attribution=excluded.attribution,
                license_note=excluded.license_note,
                metadata_json=excluded.metadata_json
            """,
            (
                "external:europepmc", "Europe PMC", "europepmc", None,
                snapshot.provider_versions.get("europepmc"),
                snapshot.literature_cutoff, "Europe PMC",
                "Publication content remains subject to source copyright and license terms.",
                _json({"snapshot_source": True}), snapshot.created_at,
            ),
        )
        source_id = int(connection.execute(
            "SELECT id FROM sources WHERE source_key='external:europepmc'"
        ).fetchone()[0])
        fts_exists = connection.execute(
            "SELECT 1 FROM sqlite_master WHERE type='table' AND name='publication_fts'"
        ).fetchone() is not None
        for match in matches:
            if not isinstance(match, Mapping):
                continue
            publication = match.get("publication")
            query = match.get("query")
            if not isinstance(publication, Mapping) or not isinstance(query, Mapping):
                continue
            publication_id = _required(publication.get("publication_id"), "publication_id")
            title = _required(publication.get("title"), "publication title")
            pmid = publication.get("pmid")
            pmcid = publication.get("pmcid")
            doi = publication.get("doi")
            row = connection.execute(
                """
                SELECT id, publication_id FROM publications
                WHERE publication_id=?
                   OR (? IS NOT NULL AND pmid=?)
                   OR (? IS NOT NULL AND pmcid=?)
                   OR (? IS NOT NULL AND doi=?)
                ORDER BY id LIMIT 1
                """,
                (publication_id, pmid, pmid, pmcid, pmcid, doi, doi),
            ).fetchone()
            publication_metadata = _json({
                "provenance": publication.get("provenance", {}),
                **dict(publication.get("metadata", {})),
                "underlying_source": publication.get("source"),
            })
            if row is None:
                cursor = connection.execute(
                    """
                    INSERT INTO publications(
                        publication_id, pmid, pmcid, doi, title, year, journal,
                        authors_summary, abstract, source_id, retrieved_at, metadata_json
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        publication_id, pmid, pmcid, doi, title,
                        publication.get("year"), publication.get("journal"),
                        publication.get("authors_summary"), publication.get("abstract"),
                        source_id, publication.get("retrieved_at") or snapshot.created_at,
                        publication_metadata,
                    ),
                )
                publication_row_id = int(cursor.lastrowid)
                stored_publication_id = publication_id
            else:
                publication_row_id = int(row["id"])
                stored_publication_id = str(row["publication_id"])
                connection.execute(
                    """
                    UPDATE publications SET
                        pmid=?, pmcid=?, doi=?, title=?, year=?, journal=?,
                        authors_summary=?, abstract=?, source_id=?, retrieved_at=?, metadata_json=?
                    WHERE id=?
                    """,
                    (
                        pmid, pmcid, doi, title, publication.get("year"),
                        publication.get("journal"), publication.get("authors_summary"),
                        publication.get("abstract"), source_id,
                        publication.get("retrieved_at") or snapshot.created_at,
                        publication_metadata, publication_row_id,
                    ),
                )
            if fts_exists:
                connection.execute("DELETE FROM publication_fts WHERE rowid=?", (publication_row_id,))
                connection.execute(
                    "INSERT INTO publication_fts(rowid, publication_id, title, abstract) VALUES (?, ?, ?, ?)",
                    (publication_row_id, stored_publication_id, title, publication.get("abstract") or ""),
                )
            match_metadata = {
                "context_snapshot_id": snapshot.context_snapshot_id,
                "source_intent_id": query.get("source_intent_id"),
                "match_label": match.get("match_label", "Co-mentioned"),
            }
            existing_matches = connection.execute(
                """
                SELECT id, metadata_json FROM publication_matches
                WHERE publication_id=? AND simulation_id=?
                  AND query_text=? AND query_version=?
                """,
                (
                    publication_row_id, snapshot.simulation_id,
                    query.get("query_text"), query.get("query_version"),
                ),
            ).fetchall()
            existing_match_id = next(
                (
                    int(item["id"])
                    for item in existing_matches
                    if json.loads(item["metadata_json"]).get("context_snapshot_id")
                    == snapshot.context_snapshot_id
                ),
                None,
            )
            values = (
                query.get("reason_shown") or "Shown for the current explicit literature query.",
                _json(query.get("matched_entities", ())),
                _json(match.get("relevant_passages", ())),
                snapshot.created_at,
                _json(match_metadata),
            )
            if existing_match_id is None:
                connection.execute(
                    """
                    INSERT INTO publication_matches(
                        publication_id, simulation_id, observation_id, entity_id,
                        query_text, query_version, reason_shown, match_type,
                        matched_entities_json, relevant_passages_json, created_at,
                        metadata_json
                    ) VALUES (?, ?, NULL, NULL, ?, ?, ?, 'CO_MENTIONED_IN_PUBLICATION', ?, ?, ?, ?)
                    """,
                    (
                        publication_row_id, snapshot.simulation_id,
                        query.get("query_text"), query.get("query_version"),
                        *values,
                    ),
                )
            else:
                connection.execute(
                    """
                    UPDATE publication_matches SET
                        reason_shown=?, matched_entities_json=?, relevant_passages_json=?,
                        created_at=?, metadata_json=? WHERE id=?
                    """,
                    (*values, existing_match_id),
                )

    def save(
        self,
        simulation_snapshot: SimulationResultSnapshot,
        snapshot: ExternalContextSnapshot,
    ) -> SnapshotSaveResult:
        if (
            snapshot.simulation_id != simulation_snapshot.simulation_id
            or snapshot.simulation_snapshot_id != simulation_snapshot.snapshot_id
        ):
            raise ValueError("external snapshot does not belong to the simulation snapshot")
        try:
            with self._store.transaction() as connection:
                existing_simulation = connection.execute(
                    "SELECT metadata_json FROM simulation_snapshots WHERE simulation_id = ?",
                    (simulation_snapshot.simulation_id,),
                ).fetchone()
                simulation_metadata = self._simulation_metadata(simulation_snapshot)
                if existing_simulation is not None:
                    existing_metadata = json.loads(existing_simulation["metadata_json"])
                    if existing_metadata.get("snapshot_id") != simulation_snapshot.snapshot_id:
                        raise ValueError("simulation ID is already bound to a different immutable snapshot")
                else:
                    connection.execute(
                        """
                        INSERT INTO simulation_snapshots(
                            simulation_id, created_at, species, taxon_id, tissue,
                            targets_json, attenuation, selection_mode, test_limit,
                            threshold_parameters_json, result_schema_version,
                            source_reference, metadata_json
                        ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                        """,
                        (
                            simulation_snapshot.simulation_id,
                            simulation_snapshot.created_at,
                            simulation_snapshot.species,
                            simulation_snapshot.taxon_id,
                            simulation_snapshot.tissue,
                            _json(simulation_snapshot.targets),
                            simulation_snapshot.attenuation,
                            simulation_snapshot.selection_mode,
                            simulation_snapshot.test_limit,
                            _json(simulation_snapshot.threshold_parameters),
                            str(simulation_snapshot.result_schema_version),
                            simulation_snapshot.snapshot_id,
                            simulation_metadata,
                        ),
                    )
                encoded = _json({
                    "snapshot": snapshot.to_dict(),
                    "payload_digest": snapshot.payload_digest,
                })
                existing = connection.execute(
                    "SELECT metadata_json FROM external_context_snapshots WHERE context_snapshot_id = ?",
                    (snapshot.context_snapshot_id,),
                ).fetchone()
                if existing is not None:
                    if str(existing["metadata_json"]) != encoded:
                        raise ValueError("context snapshot ID collision with different content")
                    return SnapshotSaveResult(False, snapshot.context_snapshot_id, "Snapshot already stored unchanged.")
                connection.execute(
                    """
                    INSERT INTO external_context_snapshots(
                        context_snapshot_id, simulation_id, created_at, query_version,
                        providers_json, provider_versions_json, literature_cutoff,
                        status, metadata_json
                    ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        snapshot.context_snapshot_id,
                        snapshot.simulation_id,
                        snapshot.created_at,
                        snapshot.query_version,
                        _json(snapshot.providers),
                        _json(snapshot.provider_versions),
                        snapshot.literature_cutoff,
                        snapshot.status.value,
                        encoded,
                    ),
                )
                self._persist_literature_matches(connection, snapshot)
            return SnapshotSaveResult(True, snapshot.context_snapshot_id)
        except ValueError:
            raise
        except (EvidenceStoreError, sqlite3.Error, TypeError) as exc:
            return SnapshotSaveResult(False, snapshot.context_snapshot_id, f"{type(exc).__name__}: {exc}")

    @staticmethod
    def _decode(row: sqlite3.Row) -> ExternalContextSnapshot:
        metadata = json.loads(row["metadata_json"])
        snapshot = ExternalContextSnapshot.from_dict(metadata["snapshot"])
        if metadata.get("payload_digest") != snapshot.payload_digest:
            raise ValueError("stored external snapshot payload digest is invalid")
        if (
            snapshot.context_snapshot_id != row["context_snapshot_id"]
            or snapshot.simulation_id != row["simulation_id"]
            or snapshot.created_at != _utc(row["created_at"])
            or snapshot.query_version != row["query_version"]
            or snapshot.status.value != row["status"]
            or list(snapshot.providers) != json.loads(row["providers_json"])
            or dict(snapshot.provider_versions) != json.loads(row["provider_versions_json"])
            or snapshot.literature_cutoff != row["literature_cutoff"]
        ):
            raise ValueError("stored external snapshot columns do not match its payload")
        return snapshot

    def restore(self, context_snapshot_id: str) -> SnapshotRestoreResult:
        identifier = _required(context_snapshot_id, "context_snapshot_id")
        try:
            with self._store.read_connection() as connection:
                row = connection.execute(
                    "SELECT * FROM external_context_snapshots WHERE context_snapshot_id = ?",
                    (identifier,),
                ).fetchone()
            if row is None:
                return SnapshotRestoreResult(SnapshotRestoreStatus.NOT_FOUND)
            return SnapshotRestoreResult(SnapshotRestoreStatus.AVAILABLE, self._decode(row))
        except (EvidenceStoreError, sqlite3.Error) as exc:
            return SnapshotRestoreResult(SnapshotRestoreStatus.UNAVAILABLE, message=str(exc))
        except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
            return SnapshotRestoreResult(SnapshotRestoreStatus.INVALID, message=str(exc))

    def restore_latest(
        self,
        simulation_id: str,
        *,
        expected_simulation_snapshot_id: str,
    ) -> SnapshotRestoreResult:
        simulation = _required(simulation_id, "simulation_id")
        expected = _required(expected_simulation_snapshot_id, "expected_simulation_snapshot_id")
        try:
            with self._store.read_connection() as connection:
                rows = connection.execute(
                    """
                    SELECT * FROM external_context_snapshots
                    WHERE simulation_id = ? ORDER BY created_at DESC, rowid DESC
                    """,
                    (simulation,),
                ).fetchall()
            if not rows:
                return SnapshotRestoreResult(SnapshotRestoreStatus.NOT_FOUND)
            try:
                snapshot = self._decode(rows[0])
            except (KeyError, TypeError, ValueError, json.JSONDecodeError) as exc:
                return SnapshotRestoreResult(SnapshotRestoreStatus.INVALID, message=str(exc))
            if snapshot.simulation_snapshot_id == expected:
                return SnapshotRestoreResult(SnapshotRestoreStatus.AVAILABLE, snapshot)
            return SnapshotRestoreResult(
                SnapshotRestoreStatus.NOT_FOUND,
                message="No cached external context matches the current simulation snapshot.",
            )
        except (EvidenceStoreError, sqlite3.Error) as exc:
            return SnapshotRestoreResult(SnapshotRestoreStatus.UNAVAILABLE, message=str(exc))


__all__ = [
    "EXTERNAL_CONTEXT_SCHEMA_VERSION",
    "MAX_EXTERNAL_SNAPSHOT_BYTES",
    "ExternalContextSnapshot",
    "ExternalContextSnapshotStore",
    "ExternalSnapshotStatus",
    "SnapshotRestoreResult",
    "SnapshotRestoreStatus",
    "SnapshotSaveResult",
    "create_external_context_snapshot",
    "refresh_external_context_snapshot",
]
