"""Persistence for externally retrieved Disease Bank records only."""

from __future__ import annotations

import json
from datetime import datetime, timezone
from typing import Mapping

from ..evidence_store import EvidenceStore
from .models import Disease, DiseaseGeneAssociation
from .source_registry import DISEASE_SOURCE_REGISTRY, SourceRegistryEntry


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _plain_json(value):
    """Thaw immutable provider payloads before SQLite JSON serialization."""

    if isinstance(value, Mapping):
        return {str(key): _plain_json(item) for key, item in value.items()}
    if isinstance(value, tuple):
        return [_plain_json(item) for item in value]
    if isinstance(value, list):
        return [_plain_json(item) for item in value]
    return value


class DiseaseRepository:
    def __init__(self, store: EvidenceStore) -> None:
        self._store = store

    def initialize(self) -> None:
        """Explicitly initialize/migrate only the research evidence store."""

        self._store.initialize()
        for source in DISEASE_SOURCE_REGISTRY.values():
            self.register_source(source)

    def register_source(self, source: SourceRegistryEntry, *, retrieved_at: str | None = None) -> int:
        return self._store.upsert_source(
            source_key=source.source_key, source_name=source.source_name, provider=source.source_name,
            source_version=source.source_version, retrieved_at=retrieved_at or _utc_now(), attribution=source.attribution,
            license_note=source.license, metadata={
                "license": source.license, "commercial_use": source.commercial_use,
                "attribution_required": source.attribution_required, "redistribution_allowed": source.redistribution_allowed,
                "source_url": source.source_url,
            },
        )

    def save(self, disease: Disease, associations: tuple[DiseaseGeneAssociation, ...], *, retrieved_at: str) -> None:
        source = DISEASE_SOURCE_REGISTRY["open_targets_platform"]
        self.register_source(source, retrieved_at=retrieved_at)
        with self._store.transaction() as connection:
            connection.execute(
                """INSERT INTO diseases(disease_id,disease_name,ontology_source,mondo_id,synonyms_json,parents_json,source_key,source_version,retrieved_at,metadata_json)
                   VALUES(?,?,?,?,?,?,?,?,?,?)
                   ON CONFLICT(disease_id) DO UPDATE SET disease_name=excluded.disease_name,ontology_source=excluded.ontology_source,
                   mondo_id=excluded.mondo_id,synonyms_json=excluded.synonyms_json,parents_json=excluded.parents_json,source_key=excluded.source_key,
                   source_version=excluded.source_version,retrieved_at=excluded.retrieved_at,metadata_json=excluded.metadata_json""",
                (disease.disease_id, disease.disease_name, disease.ontology_source, disease.mondo_id,
                 json.dumps(disease.synonyms), json.dumps(disease.parents), source.source_key, source.source_version, retrieved_at, "{}"),
            )
            connection.execute("DELETE FROM disease_aliases WHERE disease_id = ? AND source_key = ?", (disease.disease_id, source.source_key))
            for alias in disease.synonyms:
                connection.execute(
                    "INSERT OR IGNORE INTO disease_aliases(disease_id,alias_value,alias_type,source_key,created_at,metadata_json) VALUES(?,?,?,?,?,?)",
                    (disease.disease_id, alias, "synonym", source.source_key, retrieved_at, "{}"),
                )
            connection.execute("DELETE FROM disease_gene_associations WHERE disease_id = ? AND source_key = ?", (disease.disease_id, source.source_key))
            for item in associations:
                connection.execute(
                    """INSERT INTO disease_gene_associations(disease_id,gene_symbol,ensembl_gene_id,target_name,association_score,evidence_count,evidence_json,source_key,source_version,retrieved_at,metadata_json)
                       VALUES(?,?,?,?,?,?,?,?,?,?,?)""",
                    (item.disease_id, item.gene_symbol, item.ensembl_gene_id, item.target_name, item.association_score,
                     item.evidence_count, json.dumps(_plain_json(item.evidence), sort_keys=True), item.source_key, item.source_version,
                     item.retrieved_at, "{}"),
                )

    def get(self, disease_id: str) -> tuple[Disease, tuple[DiseaseGeneAssociation, ...]] | None:
        with self._store.read_connection() as connection:
            row = connection.execute("SELECT * FROM diseases WHERE disease_id = ?", (str(disease_id).strip(),)).fetchone()
            if row is None:
                return None
            associations = connection.execute(
                "SELECT * FROM disease_gene_associations WHERE disease_id = ? ORDER BY association_score DESC, gene_symbol, ensembl_gene_id",
                (row["disease_id"],),
            ).fetchall()
        disease = Disease(
            disease_id=row["disease_id"], disease_name=row["disease_name"], ontology_source=row["ontology_source"],
            mondo_id=row["mondo_id"], synonyms=tuple(json.loads(row["synonyms_json"])), parents=tuple(json.loads(row["parents_json"])),
        )
        return disease, tuple(DiseaseGeneAssociation(
            disease_id=item["disease_id"], gene_symbol=item["gene_symbol"], ensembl_gene_id=item["ensembl_gene_id"],
            target_name=item["target_name"] or "", association_score=item["association_score"], evidence_count=item["evidence_count"],
            evidence=json.loads(item["evidence_json"]), source_key=item["source_key"], source_version=item["source_version"] or "unknown",
            retrieved_at=item["retrieved_at"],
        ) for item in associations)


__all__ = ["DiseaseRepository"]
