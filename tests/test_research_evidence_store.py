from __future__ import annotations

import json
from pathlib import Path
import sqlite3
from tempfile import TemporaryDirectory

import pytest

from src.research import evidence_store as evidence_store_module
from src.research.evidence_store import (
    EvidenceStore,
    EvidenceStoreError,
    EvidenceStoreMigrationError,
    RESEARCH_SCHEMA_VERSION,
    SCHEMA_VERSION_KEY,
    UnsafeEvidenceStorePath,
    UnsupportedSchemaVersion,
)


REQUIRED_TABLES = {
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


def test_store_has_no_constructor_side_effect_and_migrates_idempotently() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "nested" / "research_context.sqlite"
        store = EvidenceStore(path)

        assert not path.exists()
        assert store.healthcheck().status == "not_initialized"

        health = store.initialize()
        assert path.is_file()
        assert health.available is True
        assert health.schema_version == RESEARCH_SCHEMA_VERSION
        assert store.schema_version() == RESEARCH_SCHEMA_VERSION

        # Reopening is an idempotent no-op migration, not a destructive rebuild.
        assert store.initialize().available is True
        with store.read_connection() as connection:
            tables = {
                row[0]
                for row in connection.execute(
                    "SELECT name FROM sqlite_master "
                    "WHERE type='table' AND name NOT LIKE 'sqlite_%'"
                )
            }
            version = connection.execute(
                "SELECT value FROM schema_meta WHERE key = ?",
                (SCHEMA_VERSION_KEY,),
            ).fetchone()[0]
            assert REQUIRED_TABLES <= tables
            assert version == str(RESEARCH_SCHEMA_VERSION)
            assert connection.execute("PRAGMA foreign_keys").fetchone()[0] == 1
            assert connection.execute("PRAGMA query_only").fetchone()[0] == 1
            provider_cache_columns = {
                row[1]
                for row in connection.execute("PRAGMA table_info(provider_cache)")
            }
            assert {"request_params_json", "request_params_hash"} <= provider_cache_columns


def test_v1_provider_cache_migrates_without_losing_rows_or_parameter_identity() -> None:
    """Schema v2 fixes collisions between requests with different parameters."""

    with TemporaryDirectory() as directory:
        path = Path(directory) / "research_context.sqlite"
        connection = sqlite3.connect(path)
        for statement in evidence_store_module._MIGRATIONS[1]:
            connection.execute(statement)
        connection.execute(
            "INSERT INTO schema_meta(key, value, updated_at) VALUES (?, ?, ?)",
            (SCHEMA_VERSION_KEY, "1", "2026-01-01T00:00:00+00:00"),
        )
        connection.execute(
            """
            INSERT INTO provider_cache(
                cache_key, provider, operation, taxon_id, canonical_query,
                provider_schema_version, query_version, status, payload_json,
                checked_at, expires_at, metadata_json
            ) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)
            """,
            (
                "legacy-key", "fixture", "entity_context", 9606, "ENSP1",
                "provider-v1", "query-v1", "AVAILABLE", '{"value":1}',
                "2026-01-01T00:00:00+00:00", "2026-01-02T00:00:00+00:00", "{}",
            ),
        )
        connection.commit()
        connection.close()

        store = EvidenceStore(path)
        health = store.initialize()
        assert health.available is True
        assert health.schema_version == RESEARCH_SCHEMA_VERSION

        with store.read_connection() as migrated:
            row = migrated.execute(
                """
                SELECT cache_key, request_params_json, request_params_hash
                FROM provider_cache WHERE cache_key = 'legacy-key'
                """
            ).fetchone()
            assert tuple(row) == (
                "legacy-key",
                "{}",
                evidence_store_module._EMPTY_PARAMS_SHA256,
            )
            index_sql = migrated.execute(
                "SELECT sql FROM sqlite_master WHERE name = 'uq_provider_cache_identity'"
            ).fetchone()[0]
            assert "request_params_hash" in index_sql


def test_entity_alias_and_source_preserve_json_and_species_identity() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research_context.sqlite")
        store.initialize()

        human = store.upsert_entity(
            taxon_id=9606,
            entity_type="protein",
            canonical_id="ENSP00000407554",
            symbol="GLP1R",
            display_name="Glucagon-like peptide 1 receptor",
            metadata={"origin": "existing_project_mapping"},
        )
        mouse = store.upsert_entity(
            taxon_id=10090,
            entity_type="protein",
            canonical_id="ENSMUSP00000033715",
            symbol="Glp1r",
            metadata={"origin": "existing_project_mapping"},
        )
        assert human != mouse

        alias_id = store.add_entity_alias(
            entity_id=human,
            alias_type="uniprot",
            alias_value="P43220",
            source="structural_binding.db",
            metadata={"provenance": "local"},
        )
        # Exact duplicate aliases coalesce and do not multiply evidence.
        assert store.add_entity_alias(
            entity_id=human,
            alias_type="uniprot",
            alias_value="P43220",
            source="structural_binding.db",
            metadata={"provenance": "local"},
        ) == alias_id

        source_id = store.upsert_source(
            source_key="local:hgnc-gene-groups",
            source_name="HGNC Gene Groups",
            provider="local_hgnc",
            source_version="2026-07-07",
            attribution="HGNC",
            metadata={"sha256": "abc123"},
        )
        assert source_id > 0

        with store.read_connection() as connection:
            entity_row = connection.execute(
                "SELECT taxon_id, metadata_json FROM entities WHERE id = ?", (human,)
            ).fetchone()
            alias_row = connection.execute(
                "SELECT metadata_json FROM entity_aliases WHERE id = ?", (alias_id,)
            ).fetchone()
            source_row = connection.execute(
                "SELECT metadata_json FROM sources WHERE id = ?", (source_id,)
            ).fetchone()

        assert entity_row["taxon_id"] == 9606
        assert json.loads(entity_row["metadata_json"])["origin"] == "existing_project_mapping"
        assert json.loads(alias_row["metadata_json"]) == {"provenance": "local"}
        assert json.loads(source_row["metadata_json"]) == {"sha256": "abc123"}


def test_transactions_roll_back_and_foreign_keys_are_enforced() -> None:
    with TemporaryDirectory() as directory:
        store = EvidenceStore(Path(directory) / "research_context.sqlite")
        store.initialize()

        with pytest.raises(RuntimeError, match="abort test transaction"):
            with store.transaction() as connection:
                connection.execute(
                    """
                    INSERT INTO entities(
                        taxon_id, entity_type, canonical_id, created_at, metadata_json
                    ) VALUES (9606, 'protein', 'ROLLBACK-ME', 'now', '{}')
                    """
                )
                raise RuntimeError("abort test transaction")

        with store.read_connection() as connection:
            assert connection.execute(
                "SELECT COUNT(*) FROM entities WHERE canonical_id = 'ROLLBACK-ME'"
            ).fetchone()[0] == 0

        with pytest.raises(EvidenceStoreError):
            store.add_entity_alias(
                entity_id=999999,
                alias_type="gene_symbol",
                alias_value="MISSING",
                source="test",
            )

        # JSON metadata is checked at both the Python and SQLite boundaries.
        with pytest.raises(ValueError):
            store.upsert_entity(
                taxon_id=9606,
                entity_type="protein",
                canonical_id="BAD-JSON-PYTHON",
                metadata={"not_finite": float("nan")},
            )
        with pytest.raises(EvidenceStoreError):
            with store.transaction() as connection:
                connection.execute(
                    """
                    INSERT INTO entities(
                        taxon_id, entity_type, canonical_id, created_at, metadata_json
                    ) VALUES (9606, 'protein', 'BAD-JSON-SQL', 'now', '{bad json')
                    """
                )


def test_nonempty_unrelated_database_is_never_adopted_or_modified() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "unrelated.sqlite"
        connection = sqlite3.connect(path)
        connection.execute("CREATE TABLE scientific_results(value REAL NOT NULL)")
        connection.execute("INSERT INTO scientific_results VALUES (1.25)")
        connection.commit()
        connection.close()

        store = EvidenceStore(path)
        with pytest.raises(EvidenceStoreMigrationError, match="non-empty database"):
            store.initialize()

        connection = sqlite3.connect(path)
        tables = {
            row[0]
            for row in connection.execute(
                "SELECT name FROM sqlite_master WHERE type='table'"
            )
        }
        value = connection.execute("SELECT value FROM scientific_results").fetchone()[0]
        connection.close()
        assert tables == {"scientific_results"}
        assert value == 1.25


def test_core_database_names_are_rejected_before_creation() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "hinterland_core.db"
        store = EvidenceStore(path)

        with pytest.raises(UnsafeEvidenceStorePath):
            store.initialize()
        assert not path.exists()


def test_corrupt_database_failure_is_non_throwing_for_healthcheck() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "research_context.sqlite"
        original = b"this is not a sqlite database"
        path.write_bytes(original)
        store = EvidenceStore(path)

        health = store.healthcheck()
        assert health.available is False
        assert health.status == "corrupt_or_unavailable"
        assert path.read_bytes() == original


def test_newer_schema_is_not_downgraded() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "research_context.sqlite"
        connection = sqlite3.connect(path)
        connection.execute(
            "CREATE TABLE schema_meta(key TEXT PRIMARY KEY, value TEXT, updated_at TEXT)"
        )
        connection.execute(
            "INSERT INTO schema_meta VALUES (?, ?, 'now')",
            (SCHEMA_VERSION_KEY, str(RESEARCH_SCHEMA_VERSION + 1)),
        )
        connection.commit()
        connection.close()

        with pytest.raises(UnsupportedSchemaVersion):
            EvidenceStore(path).initialize()

        connection = sqlite3.connect(path)
        version = connection.execute(
            "SELECT value FROM schema_meta WHERE key = ?", (SCHEMA_VERSION_KEY,)
        ).fetchone()[0]
        connection.close()
        assert version == str(RESEARCH_SCHEMA_VERSION + 1)
