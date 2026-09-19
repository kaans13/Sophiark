"""Persistent, species-aware MyGene annotation cache.

Cache hitleri hem pozitif hem negatif cevaplarÄ± saklar. Bu modÃ¼l sadece
annotation/mapping iÃ§indir; aÄŸ skorlarÄ±nÄ± veya perturbation semantics'ini deÄŸiÅŸtirmez.
"""
from __future__ import annotations

import json
import sqlite3
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Callable, Iterable


@dataclass
class CacheStats:
    memory_hits: int = 0
    persistent_hits: int = 0
    misses: int = 0
    external_batches: int = 0
    external_queries: int = 0
    negative_hits: int = 0

    def as_dict(self) -> dict[str, int]:
        return self.__dict__.copy()


class PersistentMyGeneCache:
    VERSION = "mygene-sqlite-v1"

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._memory: dict[tuple[str, str, str], dict] = {}
        self.stats = CacheStats()
        con = self._connect()
        try:
            con.execute("""CREATE TABLE IF NOT EXISTS mygene_cache (
                organism TEXT NOT NULL, query_id TEXT NOT NULL, fields_key TEXT NOT NULL,
                status TEXT NOT NULL, payload TEXT, checked_at TEXT NOT NULL,
                PRIMARY KEY (organism, query_id, fields_key)
            )""")
            con.commit()
        finally:
            con.close()

    def _connect(self):
        return sqlite3.connect(str(self.path))

    @staticmethod
    def _key(organism: str, query_id: str, fields: str) -> tuple[str, str, str]:
        return organism.casefold(), str(query_id).strip(), ",".join(sorted(part.strip() for part in fields.split(",")))

    def get(self, organism: str, query_id: str, fields: str) -> dict | None:
        key = self._key(organism, query_id, fields)
        if key in self._memory:
            self.stats.memory_hits += 1
            record = self._memory[key]
        else:
            con = self._connect()
            try:
                row = con.execute("SELECT status, payload, checked_at FROM mygene_cache WHERE organism=? AND query_id=? AND fields_key=?", key).fetchone()
            finally:
                con.close()
            if row is None:
                self.stats.misses += 1
                return None
            self.stats.persistent_hits += 1
            record = {"status": row[0], "payload": json.loads(row[1]) if row[1] else None, "checked_at": row[2]}
            self._memory[key] = record
        if record["status"] != "found":
            self.stats.negative_hits += 1
        return record

    def put(self, organism: str, query_id: str, fields: str, *, status: str, payload: dict | None) -> None:
        key = self._key(organism, query_id, fields)
        record = {"status": status, "payload": payload, "checked_at": datetime.now(timezone.utc).isoformat()}
        self._memory[key] = record
        con = self._connect()
        try:
            con.execute("INSERT OR REPLACE INTO mygene_cache VALUES (?, ?, ?, ?, ?, ?)", (*key, status, json.dumps(payload, ensure_ascii=False) if payload else None, record["checked_at"]))
            con.commit()
        finally:
            con.close()

    def resolve_batch(
        self, organism: str, ids: Iterable[str], fields: str,
        remote_batch: Callable[[list[str]], list[dict]] | None,
    ) -> dict[str, dict | None]:
        requested = [str(identifier) for identifier in ids]
        resolved: dict[str, dict | None] = {}
        unresolved: list[str] = []
        for query_id in requested:
            record = self.get(organism, query_id, fields)
            if record is None:
                unresolved.append(query_id)
            else:
                resolved[query_id] = record["payload"] if record["status"] == "found" else None
        if not unresolved or remote_batch is None:
            return resolved
        self.stats.external_batches += 1
        self.stats.external_queries += len(unresolved)
        try:
            hits = remote_batch(unresolved)
            found = {str(hit.get("query")): hit for hit in hits if isinstance(hit, dict) and not hit.get("notfound") and hit.get("query")}
            for query_id in unresolved:
                hit = found.get(query_id)
                self.put(organism, query_id, fields, status="found" if hit else "not_found", payload=hit)
                resolved[query_id] = hit
        except Exception:
            # API failure cache'lenmez; sonraki uygun Ã§alÄ±ÅŸmada tekrar denenebilir.
            for query_id in unresolved:
                resolved[query_id] = None
        return resolved

    def snapshot_identifier(self) -> str:
        stat = self.path.stat() if self.path.exists() else None
        return f"{self.VERSION}:{stat.st_size if stat else 0}:{stat.st_mtime_ns if stat else 0}"
