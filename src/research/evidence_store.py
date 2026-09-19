"""Isolated SQLite storage for the read-only Research Context layer.

The module deliberately performs no filesystem or database work at import time.
Creating the database is an explicit :meth:`EvidenceStore.initialize` action,
and existing Sophiark scientific databases are rejected as store targets.
"""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3
from typing import Any, Iterator, Mapping


RESEARCH_SCHEMA_VERSION = 3
SCHEMA_VERSION_KEY = "research_schema_version"

_CORE_DATABASE_NAMES = frozenset(
    {
        "hinterland_core.db",
        "hinterland_core_physical.db",
        "mouse_hinterland_core.db",
        "structural_binding.db",
    }
)

_REQUIRED_TABLES = frozenset(
    {
        "schema_meta",
        "entities",
        "entity_aliases",
        "simulation_snapshots",
        "observations",
        "observation_members",
        "assertions",
        "sources",
        "external_context_snapshots",
        "provider_cache",
        "publications",
        "publication_matches",
        "diseases",
        "disease_aliases",
        "disease_gene_associations",
    }
)

_REQUIRED_PROVIDER_CACHE_COLUMNS = frozenset(
    {
        "cache_key",
        "provider",
        "operation",
        "taxon_id",
        "canonical_query",
        "provider_schema_version",
        "query_version",
        "request_params_json",
        "request_params_hash",
        "status",
        "payload_json",
        "checked_at",
        "expires_at",
        "metadata_json",
    }
)


class EvidenceStoreError(RuntimeError):
    """Base error raised at the research-storage boundary."""


class EvidenceStoreUnavailable(EvidenceStoreError):
    """The research store is missing, inaccessible, or corrupt."""


class EvidenceStoreMigrationError(EvidenceStoreError):
    """The database cannot be migrated without risking unrelated data."""


class UnsupportedSchemaVersion(EvidenceStoreMigrationError):
    """The database was created by a newer, unsupported schema version."""


class UnsafeEvidenceStorePath(EvidenceStoreError):
    """A scientific-core database was supplied as the research store."""


@dataclass(frozen=True)
class StoreHealth:
    """Non-throwing health result suitable for provider/debug UI."""

    available: bool
    status: str
    schema_version: int | None = None
    detail: str | None = None


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _json_text(value: Any, *, field: str) -> str:
    try:
        return json.dumps(
            value,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
            allow_nan=False,
        )
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{field} must be valid JSON data") from exc


_MIGRATION_1 = (
    """
    CREATE TABLE schema_meta (
        key TEXT PRIMARY KEY,
        value TEXT NOT NULL,
        updated_at TEXT NOT NULL
    )
    """,
    """
    CREATE TABLE entities (
        id INTEGER PRIMARY KEY,
        taxon_id INTEGER,
        entity_type TEXT NOT NULL,
        canonical_id TEXT NOT NULL,
        symbol TEXT,
        display_name TEXT,
        created_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE UNIQUE INDEX uq_entities_canonical
    ON entities(COALESCE(taxon_id, -1), entity_type, canonical_id)
    """,
    """
    CREATE INDEX idx_entities_symbol
    ON entities(taxon_id, symbol)
    """,
    """
    CREATE TABLE entity_aliases (
        id INTEGER PRIMARY KEY,
        entity_id INTEGER NOT NULL REFERENCES entities(id) ON DELETE CASCADE,
        alias_type TEXT NOT NULL,
        alias_value TEXT NOT NULL,
        source TEXT NOT NULL,
        created_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json)),
        UNIQUE(entity_id, alias_type, alias_value, source)
    )
    """,
    """
    CREATE INDEX idx_entity_alias_lookup
    ON entity_aliases(alias_type, alias_value)
    """,
    """
    CREATE TABLE simulation_snapshots (
        simulation_id TEXT PRIMARY KEY,
        created_at TEXT NOT NULL,
        species TEXT NOT NULL,
        taxon_id INTEGER NOT NULL,
        tissue TEXT,
        targets_json TEXT NOT NULL CHECK (json_valid(targets_json)),
        attenuation REAL,
        selection_mode TEXT,
        test_limit INTEGER,
        threshold_parameters_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(threshold_parameters_json)),
        result_schema_version TEXT NOT NULL,
        source_reference TEXT,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE INDEX idx_simulation_snapshots_context
    ON simulation_snapshots(taxon_id, tissue, created_at)
    """,
    """
    CREATE TABLE observations (
        observation_id TEXT PRIMARY KEY,
        simulation_id TEXT NOT NULL
            REFERENCES simulation_snapshots(simulation_id) ON DELETE CASCADE,
        type TEXT NOT NULL,
        basis TEXT NOT NULL,
        scope_json TEXT NOT NULL CHECK (json_valid(scope_json)),
        source_fields_json TEXT NOT NULL CHECK (json_valid(source_fields_json)),
        created_at TEXT NOT NULL,
        limitations TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE INDEX idx_observations_simulation_type
    ON observations(simulation_id, type)
    """,
    """
    CREATE TABLE observation_members (
        observation_id TEXT NOT NULL
            REFERENCES observations(observation_id) ON DELETE CASCADE,
        entity_id INTEGER NOT NULL REFERENCES entities(id) ON DELETE RESTRICT,
        role TEXT NOT NULL DEFAULT 'member',
        ordinal INTEGER,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json)),
        PRIMARY KEY(observation_id, entity_id, role)
    )
    """,
    """
    CREATE INDEX idx_observation_members_entity
    ON observation_members(entity_id)
    """,
    """
    CREATE TABLE sources (
        id INTEGER PRIMARY KEY,
        source_key TEXT NOT NULL UNIQUE,
        source_name TEXT NOT NULL,
        provider TEXT NOT NULL,
        source_identifier TEXT,
        source_version TEXT,
        retrieved_at TEXT,
        attribution TEXT,
        license_note TEXT,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json)),
        created_at TEXT NOT NULL
    )
    """,
    """
    CREATE INDEX idx_sources_provider
    ON sources(provider)
    """,
    """
    CREATE TABLE assertions (
        id INTEGER PRIMARY KEY,
        subject_entity_id INTEGER NOT NULL
            REFERENCES entities(id) ON DELETE CASCADE,
        predicate TEXT NOT NULL,
        object_entity_id INTEGER REFERENCES entities(id) ON DELETE CASCADE,
        literal_value TEXT,
        source_id INTEGER NOT NULL REFERENCES sources(id) ON DELETE RESTRICT,
        source_identifier TEXT,
        taxon_id INTEGER,
        evidence_type TEXT NOT NULL,
        created_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json)),
        CHECK (
            (object_entity_id IS NOT NULL AND literal_value IS NULL)
            OR (object_entity_id IS NULL AND literal_value IS NOT NULL)
        )
    )
    """,
    """
    CREATE INDEX idx_assertions_subject_predicate
    ON assertions(subject_entity_id, predicate)
    """,
    """
    CREATE INDEX idx_assertions_object_predicate
    ON assertions(object_entity_id, predicate)
    """,
    """
    CREATE INDEX idx_assertions_source
    ON assertions(source_id)
    """,
    """
    CREATE INDEX idx_assertions_taxon_evidence
    ON assertions(taxon_id, evidence_type)
    """,
    """
    CREATE TABLE external_context_snapshots (
        context_snapshot_id TEXT PRIMARY KEY,
        simulation_id TEXT NOT NULL
            REFERENCES simulation_snapshots(simulation_id) ON DELETE CASCADE,
        created_at TEXT NOT NULL,
        query_version TEXT NOT NULL,
        providers_json TEXT NOT NULL CHECK (json_valid(providers_json)),
        provider_versions_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(provider_versions_json)),
        literature_cutoff TEXT,
        status TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE INDEX idx_external_snapshots_simulation
    ON external_context_snapshots(simulation_id, created_at)
    """,
    """
    CREATE TABLE provider_cache (
        cache_key TEXT PRIMARY KEY,
        provider TEXT NOT NULL,
        operation TEXT NOT NULL,
        taxon_id INTEGER,
        canonical_query TEXT NOT NULL,
        provider_schema_version TEXT NOT NULL,
        query_version TEXT NOT NULL,
        status TEXT NOT NULL,
        payload_json TEXT CHECK (payload_json IS NULL OR json_valid(payload_json)),
        checked_at TEXT NOT NULL,
        expires_at TEXT,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE UNIQUE INDEX uq_provider_cache_identity
    ON provider_cache(
        provider,
        operation,
        COALESCE(taxon_id, -1),
        canonical_query,
        provider_schema_version,
        query_version
    )
    """,
    """
    CREATE INDEX idx_provider_cache_lookup
    ON provider_cache(
        provider,
        operation,
        taxon_id,
        canonical_query,
        provider_schema_version,
        query_version
    )
    """,
    """
    CREATE INDEX idx_provider_cache_expiry
    ON provider_cache(expires_at, status)
    """,
    """
    CREATE TABLE publications (
        id INTEGER PRIMARY KEY,
        publication_id TEXT NOT NULL UNIQUE,
        pmid TEXT,
        pmcid TEXT,
        doi TEXT,
        title TEXT NOT NULL,
        year INTEGER,
        journal TEXT,
        authors_summary TEXT,
        abstract TEXT,
        source_id INTEGER NOT NULL REFERENCES sources(id) ON DELETE RESTRICT,
        retrieved_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE UNIQUE INDEX uq_publications_pmid
    ON publications(pmid) WHERE pmid IS NOT NULL
    """,
    """
    CREATE UNIQUE INDEX uq_publications_pmcid
    ON publications(pmcid) WHERE pmcid IS NOT NULL
    """,
    """
    CREATE UNIQUE INDEX uq_publications_doi
    ON publications(doi) WHERE doi IS NOT NULL
    """,
    """
    CREATE INDEX idx_publications_source_year
    ON publications(source_id, year)
    """,
    """
    CREATE TABLE publication_matches (
        id INTEGER PRIMARY KEY,
        publication_id INTEGER NOT NULL
            REFERENCES publications(id) ON DELETE CASCADE,
        simulation_id TEXT NOT NULL
            REFERENCES simulation_snapshots(simulation_id) ON DELETE CASCADE,
        observation_id TEXT REFERENCES observations(observation_id) ON DELETE SET NULL,
        entity_id INTEGER REFERENCES entities(id) ON DELETE SET NULL,
        query_text TEXT NOT NULL,
        query_version TEXT NOT NULL,
        reason_shown TEXT NOT NULL,
        match_type TEXT NOT NULL DEFAULT 'CO_MENTIONED_IN_PUBLICATION',
        matched_entities_json TEXT NOT NULL DEFAULT '[]'
            CHECK (json_valid(matched_entities_json)),
        relevant_passages_json TEXT NOT NULL DEFAULT '[]'
            CHECK (json_valid(relevant_passages_json)),
        created_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}'
            CHECK (json_valid(metadata_json))
    )
    """,
    """
    CREATE INDEX idx_publication_matches_simulation
    ON publication_matches(simulation_id)
    """,
    """
    CREATE INDEX idx_publication_matches_publication
    ON publication_matches(publication_id)
    """,
    """
    CREATE INDEX idx_publication_matches_observation
    ON publication_matches(observation_id)
    """,
    """
    CREATE INDEX idx_publication_matches_entity
    ON publication_matches(entity_id)
    """,
)

_EMPTY_PARAMS_SHA256 = "44136fa355b3678a1146ad16f7e8649e94fb4fc21fe77e8310c060f61caaff8a"

_MIGRATION_2: tuple[str, ...] = (
    """
    ALTER TABLE provider_cache
    ADD COLUMN request_params_json TEXT NOT NULL DEFAULT '{}'
        CHECK (json_valid(request_params_json))
    """,
    f"""
    ALTER TABLE provider_cache
    ADD COLUMN request_params_hash TEXT NOT NULL DEFAULT '{_EMPTY_PARAMS_SHA256}'
    """,
    "DROP INDEX uq_provider_cache_identity",
    """
    CREATE UNIQUE INDEX uq_provider_cache_identity
    ON provider_cache(
        provider,
        operation,
        COALESCE(taxon_id, -1),
        canonical_query,
        provider_schema_version,
        query_version,
        request_params_hash
    )
    """,
    "DROP INDEX idx_provider_cache_lookup",
    """
    CREATE INDEX idx_provider_cache_lookup
    ON provider_cache(
        provider,
        operation,
        taxon_id,
        canonical_query,
        provider_schema_version,
        query_version,
        request_params_hash
    )
    """,
)

# Disease Bank remains inside the research-only database.  These tables carry
# external reference evidence only; they deliberately have no foreign key into
# a simulation or a scientific-core table.
_MIGRATION_3: tuple[str, ...] = (
    """
    CREATE TABLE diseases (
        disease_id TEXT PRIMARY KEY,
        disease_name TEXT NOT NULL,
        ontology_source TEXT NOT NULL,
        mondo_id TEXT,
        synonyms_json TEXT NOT NULL DEFAULT '[]' CHECK (json_valid(synonyms_json)),
        parents_json TEXT NOT NULL DEFAULT '[]' CHECK (json_valid(parents_json)),
        source_key TEXT NOT NULL REFERENCES sources(source_key) ON DELETE RESTRICT,
        source_version TEXT,
        retrieved_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}' CHECK (json_valid(metadata_json))
    )
    """,
    "CREATE INDEX idx_diseases_name ON diseases(disease_name)",
    "CREATE INDEX idx_diseases_mondo ON diseases(mondo_id)",
    """
    CREATE TABLE disease_aliases (
        id INTEGER PRIMARY KEY,
        disease_id TEXT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
        alias_value TEXT NOT NULL,
        alias_type TEXT NOT NULL DEFAULT 'synonym',
        source_key TEXT NOT NULL REFERENCES sources(source_key) ON DELETE RESTRICT,
        created_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}' CHECK (json_valid(metadata_json)),
        UNIQUE(disease_id, alias_value, alias_type, source_key)
    )
    """,
    "CREATE INDEX idx_disease_aliases_lookup ON disease_aliases(alias_value)",
    """
    CREATE TABLE disease_gene_associations (
        id INTEGER PRIMARY KEY,
        disease_id TEXT NOT NULL REFERENCES diseases(disease_id) ON DELETE CASCADE,
        gene_symbol TEXT NOT NULL DEFAULT '',
        ensembl_gene_id TEXT NOT NULL DEFAULT '',
        target_name TEXT,
        association_score REAL,
        evidence_count INTEGER,
        evidence_json TEXT NOT NULL DEFAULT '{}' CHECK (json_valid(evidence_json)),
        source_key TEXT NOT NULL REFERENCES sources(source_key) ON DELETE RESTRICT,
        source_version TEXT,
        retrieved_at TEXT NOT NULL,
        metadata_json TEXT NOT NULL DEFAULT '{}' CHECK (json_valid(metadata_json)),
        UNIQUE(disease_id, ensembl_gene_id, gene_symbol, source_key)
    )
    """,
    "CREATE INDEX idx_disease_associations_disease ON disease_gene_associations(disease_id)",
    "CREATE INDEX idx_disease_associations_ensembl ON disease_gene_associations(ensembl_gene_id)",
    "CREATE INDEX idx_disease_associations_symbol ON disease_gene_associations(gene_symbol)",
)

_MIGRATIONS: Mapping[int, tuple[str, ...]] = {
    1: _MIGRATION_1,
    2: _MIGRATION_2,
    3: _MIGRATION_3,
}


class EvidenceStore:
    """Explicitly initialized, research-only SQLite evidence store."""

    def __init__(self, path: str | Path, *, timeout_seconds: float = 5.0):
        self.path = Path(path)
        self.timeout_seconds = float(timeout_seconds)

    def _validate_target(self) -> None:
        if self.path.name.casefold() in _CORE_DATABASE_NAMES:
            raise UnsafeEvidenceStorePath(
                f"Research EvidenceStore refuses scientific-core database: {self.path.name}"
            )

    def _connect(self, *, read_only: bool) -> sqlite3.Connection:
        self._validate_target()
        try:
            if read_only:
                if not self.path.is_file():
                    raise EvidenceStoreUnavailable(
                        f"Research evidence store is not initialized: {self.path}"
                    )
                connection = sqlite3.connect(
                    self.path.resolve().as_uri() + "?mode=ro",
                    uri=True,
                    timeout=self.timeout_seconds,
                    isolation_level=None,
                )
            else:
                connection = sqlite3.connect(
                    str(self.path),
                    timeout=self.timeout_seconds,
                    isolation_level=None,
                )
            connection.row_factory = sqlite3.Row
            connection.execute("PRAGMA foreign_keys = ON")
            connection.execute(f"PRAGMA busy_timeout = {max(0, int(self.timeout_seconds * 1000))}")
            if read_only:
                connection.execute("PRAGMA query_only = ON")
            return connection
        except EvidenceStoreError:
            raise
        except (OSError, sqlite3.Error) as exc:
            raise EvidenceStoreUnavailable(
                f"Research evidence store is unavailable: {self.path}"
            ) from exc

    @staticmethod
    def _current_version(connection: sqlite3.Connection) -> int:
        tables = {
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        }
        if "schema_meta" not in tables:
            if tables:
                raise EvidenceStoreMigrationError(
                    "Refusing to initialize a non-empty database without research schema metadata"
                )
            return 0
        row = connection.execute(
            "SELECT value FROM schema_meta WHERE key = ?", (SCHEMA_VERSION_KEY,)
        ).fetchone()
        if row is None:
            raise EvidenceStoreMigrationError(
                "Research schema metadata exists but has no schema version"
            )
        try:
            return int(row[0])
        except (TypeError, ValueError) as exc:
            raise EvidenceStoreMigrationError(
                f"Invalid research schema version: {row[0]!r}"
            ) from exc

    @staticmethod
    def _verify_schema(connection: sqlite3.Connection) -> None:
        tables = {
            str(row[0])
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
            )
        }
        missing = sorted(_REQUIRED_TABLES - tables)
        if missing:
            raise EvidenceStoreMigrationError(
                "Research schema is incomplete; missing tables: " + ", ".join(missing)
            )
        provider_cache_columns = {
            str(row[1])
            for row in connection.execute("PRAGMA table_info(provider_cache)")
        }
        missing_provider_cache_columns = sorted(
            _REQUIRED_PROVIDER_CACHE_COLUMNS - provider_cache_columns
        )
        if missing_provider_cache_columns:
            raise EvidenceStoreMigrationError(
                "Research provider cache schema is incomplete; missing columns: "
                + ", ".join(missing_provider_cache_columns)
            )
        violations = connection.execute("PRAGMA foreign_key_check").fetchall()
        if violations:
            raise EvidenceStoreMigrationError(
                f"Research schema contains {len(violations)} foreign-key violation(s)"
            )

    def initialize(self) -> StoreHealth:
        """Create or migrate the research DB in one explicit transaction."""

        self._validate_target()
        try:
            self.path.parent.mkdir(parents=True, exist_ok=True)
        except OSError as exc:
            raise EvidenceStoreUnavailable(
                f"Research evidence directory cannot be created: {self.path.parent}"
            ) from exc

        connection: sqlite3.Connection | None = None
        try:
            connection = self._connect(read_only=False)
            connection.execute("BEGIN IMMEDIATE")
            current = self._current_version(connection)
            if current > RESEARCH_SCHEMA_VERSION:
                raise UnsupportedSchemaVersion(
                    f"Research schema {current} is newer than supported version "
                    f"{RESEARCH_SCHEMA_VERSION}"
                )
            for version in range(current + 1, RESEARCH_SCHEMA_VERSION + 1):
                statements = _MIGRATIONS.get(version)
                if statements is None:
                    raise EvidenceStoreMigrationError(
                        f"No migration registered for research schema version {version}"
                    )
                for statement in statements:
                    connection.execute(statement)
                connection.execute(
                    "INSERT OR REPLACE INTO schema_meta(key, value, updated_at) VALUES (?, ?, ?)",
                    (SCHEMA_VERSION_KEY, str(version), _utc_now()),
                )
            self._verify_schema(connection)
            connection.commit()
            return StoreHealth(True, "available", RESEARCH_SCHEMA_VERSION)
        except EvidenceStoreError:
            if connection is not None:
                connection.rollback()
            raise
        except (OSError, sqlite3.Error) as exc:
            if connection is not None:
                connection.rollback()
            raise EvidenceStoreUnavailable(
                f"Research evidence store initialization failed: {self.path}"
            ) from exc
        finally:
            if connection is not None:
                connection.close()

    def schema_version(self) -> int:
        """Read the current schema version without creating the database."""

        try:
            with self.read_connection() as connection:
                return self._current_version(connection)
        except EvidenceStoreError:
            raise
        except sqlite3.Error as exc:
            raise EvidenceStoreUnavailable(
                f"Research evidence store schema cannot be read: {self.path}"
            ) from exc

    def healthcheck(self) -> StoreHealth:
        """Return availability without leaking storage failures into Sophiark core."""

        if not self.path.is_file():
            return StoreHealth(False, "not_initialized")
        try:
            with self.read_connection() as connection:
                version = self._current_version(connection)
                if version != RESEARCH_SCHEMA_VERSION:
                    return StoreHealth(
                        False,
                        "schema_version_mismatch",
                        version,
                        f"supported={RESEARCH_SCHEMA_VERSION}",
                    )
                self._verify_schema(connection)
                quick_check = connection.execute("PRAGMA quick_check(1)").fetchone()
                if quick_check is None or str(quick_check[0]).casefold() != "ok":
                    return StoreHealth(False, "corrupt_or_unavailable", version)
                return StoreHealth(True, "available", version)
        except (EvidenceStoreError, OSError, sqlite3.Error) as exc:
            return StoreHealth(False, "corrupt_or_unavailable", detail=str(exc))

    @contextmanager
    def read_connection(self) -> Iterator[sqlite3.Connection]:
        """Yield a query-only connection; never creates a missing DB."""

        connection: sqlite3.Connection | None = None
        try:
            connection = self._connect(read_only=True)
            yield connection
        except EvidenceStoreError:
            raise
        except sqlite3.Error as exc:
            raise EvidenceStoreUnavailable(
                f"Research evidence read failed: {self.path}"
            ) from exc
        finally:
            if connection is not None:
                connection.close()

    @contextmanager
    def transaction(self) -> Iterator[sqlite3.Connection]:
        """Yield a foreign-key-safe write transaction with automatic rollback."""

        if not self.path.is_file():
            raise EvidenceStoreUnavailable(
                "EvidenceStore.initialize() must be called before writing"
            )
        connection: sqlite3.Connection | None = None
        try:
            connection = self._connect(read_only=False)
            if self._current_version(connection) != RESEARCH_SCHEMA_VERSION:
                raise EvidenceStoreMigrationError("Research schema must be migrated before writing")
            connection.execute("BEGIN IMMEDIATE")
            yield connection
            connection.commit()
        except EvidenceStoreError:
            if connection is not None:
                connection.rollback()
            raise
        except sqlite3.Error as exc:
            if connection is not None:
                connection.rollback()
            raise EvidenceStoreError("Research evidence transaction failed") from exc
        except BaseException:
            if connection is not None:
                connection.rollback()
            raise
        finally:
            if connection is not None:
                connection.close()

    def upsert_entity(
        self,
        *,
        taxon_id: int | None,
        entity_type: str,
        canonical_id: str,
        symbol: str | None = None,
        display_name: str | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> int:
        """Insert or update a canonical entity and return its stable row ID."""

        metadata_json = _json_text(dict(metadata or {}), field="entity metadata")
        with self.transaction() as connection:
            row = connection.execute(
                """
                SELECT id FROM entities
                WHERE COALESCE(taxon_id, -1) = COALESCE(?, -1)
                  AND entity_type = ? AND canonical_id = ?
                """,
                (taxon_id, entity_type, canonical_id),
            ).fetchone()
            if row is None:
                cursor = connection.execute(
                    """
                    INSERT INTO entities(
                        taxon_id, entity_type, canonical_id, symbol,
                        display_name, created_at, metadata_json
                    ) VALUES (?, ?, ?, ?, ?, ?, ?)
                    """,
                    (
                        taxon_id,
                        entity_type,
                        canonical_id,
                        symbol,
                        display_name,
                        _utc_now(),
                        metadata_json,
                    ),
                )
                return int(cursor.lastrowid)
            entity_id = int(row[0])
            connection.execute(
                """
                UPDATE entities
                SET symbol = ?, display_name = ?, metadata_json = ?
                WHERE id = ?
                """,
                (symbol, display_name, metadata_json, entity_id),
            )
            return entity_id

    def add_entity_alias(
        self,
        *,
        entity_id: int,
        alias_type: str,
        alias_value: str,
        source: str,
        metadata: Mapping[str, Any] | None = None,
    ) -> int:
        """Add an exact, source-attributed alias without symbol heuristics."""

        with self.transaction() as connection:
            connection.execute(
                """
                INSERT OR IGNORE INTO entity_aliases(
                    entity_id, alias_type, alias_value, source,
                    created_at, metadata_json
                ) VALUES (?, ?, ?, ?, ?, ?)
                """,
                (
                    entity_id,
                    alias_type,
                    alias_value,
                    source,
                    _utc_now(),
                    _json_text(dict(metadata or {}), field="alias metadata"),
                ),
            )
            row = connection.execute(
                """
                SELECT id FROM entity_aliases
                WHERE entity_id = ? AND alias_type = ?
                  AND alias_value = ? AND source = ?
                """,
                (entity_id, alias_type, alias_value, source),
            ).fetchone()
            if row is None:  # pragma: no cover - guarded by the INSERT/constraints
                raise EvidenceStoreError("Entity alias was not persisted")
            return int(row[0])

    def upsert_source(
        self,
        *,
        source_key: str,
        source_name: str,
        provider: str,
        source_identifier: str | None = None,
        source_version: str | None = None,
        retrieved_at: str | None = None,
        attribution: str | None = None,
        license_note: str | None = None,
        metadata: Mapping[str, Any] | None = None,
    ) -> int:
        """Persist provenance independently from scientific result data."""

        with self.transaction() as connection:
            connection.execute(
                """
                INSERT INTO sources(
                    source_key, source_name, provider, source_identifier,
                    source_version, retrieved_at, attribution, license_note,
                    metadata_json, created_at
                ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
                ON CONFLICT(source_key) DO UPDATE SET
                    source_name = excluded.source_name,
                    provider = excluded.provider,
                    source_identifier = excluded.source_identifier,
                    source_version = excluded.source_version,
                    retrieved_at = excluded.retrieved_at,
                    attribution = excluded.attribution,
                    license_note = excluded.license_note,
                    metadata_json = excluded.metadata_json
                """,
                (
                    source_key,
                    source_name,
                    provider,
                    source_identifier,
                    source_version,
                    retrieved_at,
                    attribution,
                    license_note,
                    _json_text(dict(metadata or {}), field="source metadata"),
                    _utc_now(),
                ),
            )
            row = connection.execute(
                "SELECT id FROM sources WHERE source_key = ?", (source_key,)
            ).fetchone()
            if row is None:  # pragma: no cover - guarded by the UPSERT
                raise EvidenceStoreError("Evidence source was not persisted")
            return int(row[0])


__all__ = [
    "EvidenceStore",
    "EvidenceStoreError",
    "EvidenceStoreMigrationError",
    "EvidenceStoreUnavailable",
    "RESEARCH_SCHEMA_VERSION",
    "SCHEMA_VERSION_KEY",
    "StoreHealth",
    "UnsafeEvidenceStorePath",
    "UnsupportedSchemaVersion",
]
