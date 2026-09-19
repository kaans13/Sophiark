"""SQLite persistence for the normalized local signaling dataset."""

from __future__ import annotations

import hashlib
import json
import os
import sqlite3
import tempfile
from dataclasses import asdict
from pathlib import Path

from .models import (
    DatasetQA, LayerMembership, MappingReport, SignalingDataset, SignalingEdge,
    SignalingNode, SignStatus,
)


SCHEMA_VERSION = 2


def _json(value) -> str:
    return json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))


def _optional_bool(value):
    return None if value is None else bool(value)


def write_signaling_cache(
    dataset: SignalingDataset,
    destination: str | Path,
    *,
    source_snapshot: str | Path,
) -> Path:
    """Atomically write a normalized cache; callers must invoke explicitly."""

    destination = Path(destination)
    destination.parent.mkdir(parents=True, exist_ok=True)
    source_snapshot = Path(source_snapshot)
    source_checksum = hashlib.sha256(source_snapshot.read_bytes()).hexdigest()
    handle, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.stem}-", suffix=".sqlite.tmp", dir=destination.parent,
    )
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        connection = sqlite3.connect(temporary)
        try:
            connection.executescript("""
                PRAGMA journal_mode=OFF;
                PRAGMA synchronous=FULL;
                CREATE TABLE metadata (key TEXT PRIMARY KEY, value_json TEXT NOT NULL);
                CREATE TABLE nodes (
                    node_id TEXT PRIMARY KEY,
                    source_identifier TEXT NOT NULL,
                    symbol TEXT NOT NULL,
                    layer_membership TEXT NOT NULL,
                    layer1_id TEXT,
                    mapping_status TEXT NOT NULL
                );
                CREATE TABLE omnipath_signaling_edges (
                    source TEXT NOT NULL,
                    target TEXT NOT NULL,
                    source_identifier TEXT NOT NULL,
                    source_symbol TEXT NOT NULL,
                    target_identifier TEXT NOT NULL,
                    target_symbol TEXT NOT NULL,
                    directed INTEGER NOT NULL,
                    stimulation INTEGER NOT NULL,
                    inhibition INTEGER NOT NULL,
                    canonical_sign INTEGER NOT NULL,
                    sign_status TEXT NOT NULL,
                    interaction TEXT NOT NULL,
                    resources_json TEXT NOT NULL,
                    references_json TEXT NOT NULL,
                    evidence_metadata_json TEXT NOT NULL,
                    consensus_direction INTEGER,
                    consensus_stimulation INTEGER,
                    consensus_inhibition INTEGER,
                    dataset_version TEXT NOT NULL,
                    PRIMARY KEY (source, target)
                );
                CREATE INDEX idx_omnipath_signaling_source ON omnipath_signaling_edges(source);
                CREATE INDEX idx_omnipath_signaling_target ON omnipath_signaling_edges(target);
                CREATE INDEX idx_omnipath_signaling_pair ON omnipath_signaling_edges(source, target);
            """)
            meta = {
                "schema_version": SCHEMA_VERSION,
                "source_snapshot": str(source_snapshot),
                "source_checksum_sha256": source_checksum,
                "fingerprint": dataset.fingerprint,
                "metadata": dict(dataset.metadata),
                "qa": asdict(dataset.qa),
                "mapping": asdict(dataset.mapping),
            }
            connection.executemany(
                "INSERT INTO metadata(key, value_json) VALUES (?, ?)",
                ((key, _json(value)) for key, value in meta.items()),
            )
            connection.executemany(
                "INSERT INTO nodes VALUES (?, ?, ?, ?, ?, ?)",
                ((node.node_id, node.source_identifier, node.symbol, node.layer_membership.value,
                  node.layer1_id, node.mapping_status) for node in dataset.nodes),
            )
            connection.executemany(
                "INSERT INTO omnipath_signaling_edges VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
                ((
                    edge.source, edge.target, edge.source_identifier, edge.source_symbol,
                    edge.target_identifier, edge.target_symbol, int(edge.directed),
                    int(edge.stimulation), int(edge.inhibition), edge.canonical_sign,
                    edge.sign_status.value, edge.interaction, _json(edge.resources),
                    _json(edge.references), _json(dict(edge.evidence_metadata)),
                    None if edge.consensus_direction is None else int(edge.consensus_direction),
                    None if edge.consensus_stimulation is None else int(edge.consensus_stimulation),
                    None if edge.consensus_inhibition is None else int(edge.consensus_inhibition),
                    edge.dataset_version,
                ) for edge in dataset.edges),
            )
            connection.commit()
        finally:
            connection.close()
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()
    return destination


def load_signaling_cache(
    path: str | Path,
    *,
    source_snapshot: str | Path | None = None,
) -> SignalingDataset:
    """Read a cache and reject schema/source-version mismatches."""

    path = Path(path)
    if not path.is_file():
        raise FileNotFoundError(path)
    connection = sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True)
    try:
        meta = {
            key: json.loads(value)
            for key, value in connection.execute("SELECT key, value_json FROM metadata")
        }
        if int(meta.get("schema_version", 0)) != SCHEMA_VERSION:
            raise ValueError("Unsupported OmniPath signaling cache schema")
        if source_snapshot is not None:
            checksum = hashlib.sha256(Path(source_snapshot).read_bytes()).hexdigest()
            if checksum != meta.get("source_checksum_sha256"):
                raise ValueError("Stale OmniPath signaling cache")
        nodes = tuple(SignalingNode(
            node_id=row[0], source_identifier=row[1], symbol=row[2],
            layer_membership=LayerMembership(row[3]), layer1_id=row[4], mapping_status=row[5],
        ) for row in connection.execute(
            "SELECT node_id, source_identifier, symbol, layer_membership, layer1_id, mapping_status FROM nodes ORDER BY node_id"
        ))
        edges = tuple(SignalingEdge(
            source=row[0], target=row[1], source_identifier=row[2], source_symbol=row[3],
            target_identifier=row[4], target_symbol=row[5], directed=bool(row[6]),
            stimulation=bool(row[7]), inhibition=bool(row[8]), canonical_sign=int(row[9]),
            sign_status=SignStatus(row[10]), interaction=row[11],
            resources=tuple(json.loads(row[12])), references=tuple(json.loads(row[13])),
            evidence_metadata=json.loads(row[14]), consensus_direction=_optional_bool(row[15]),
            consensus_stimulation=_optional_bool(row[16]), consensus_inhibition=_optional_bool(row[17]),
            dataset_version=row[18],
        ) for row in connection.execute("""
            SELECT source, target, source_identifier, source_symbol, target_identifier, target_symbol,
                   directed, stimulation, inhibition, canonical_sign, sign_status, interaction,
                   resources_json, references_json, evidence_metadata_json, consensus_direction,
                   consensus_stimulation, consensus_inhibition, dataset_version
            FROM omnipath_signaling_edges ORDER BY source, target
        """))
    finally:
        connection.close()
    return SignalingDataset(
        nodes=nodes, edges=edges, qa=DatasetQA(**meta["qa"]),
        mapping=MappingReport(**meta["mapping"]), metadata=meta["metadata"],
        fingerprint=str(meta["fingerprint"]),
    )
