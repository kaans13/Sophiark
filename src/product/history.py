"""Rebuildable history index backed by durable `.sophiark` analysis bundles."""

from __future__ import annotations

from dataclasses import dataclass
from contextlib import closing
from datetime import datetime
from enum import Enum
import hashlib
import logging
import os
from pathlib import Path
import sqlite3
import uuid
from typing import Mapping

from src.unified.service import UnifiedResearchReport, UnifiedStatus
from .analysis_bundle import AnalysisBundleReader, AnalysisBundleWriter, BundleError, SavedUnifiedResult
from .paths import ProjectPaths


log = logging.getLogger(__name__)


class ImportStatus(str, Enum):
    IMPORTED = "IMPORTED"
    ALREADY_IMPORTED = "ALREADY_IMPORTED"
    ID_CONFLICT = "ID_CONFLICT"


@dataclass(frozen=True, slots=True)
class HistoryEntry:
    analysis_id: str
    target: str
    tissue: str
    created_at: str
    analysis_status: str
    sophiark_version: str
    bundle_schema_version: int
    bundle_filename: str
    bundle_checksum: str
    content_fingerprint: str
    classic_status: str
    directed_status: str
    evidence_available: bool
    complex_context_available: bool


@dataclass(frozen=True, slots=True)
class RebuildReport:
    valid_bundles: int
    indexed: int
    corrupt: int
    duplicate_ids: int
    corrupt_files: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ImportResult:
    status: ImportStatus
    entry: HistoryEntry | None
    message: str = ""


def _checksum(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _safe_stem(target: str, tissue: str) -> str:
    def clean(value: str) -> str:
        filtered = "".join(character if character.isalnum() or character in "-_" else "_" for character in value)
        return filtered.strip("_")[:48] or "analysis"
    return f"{clean(target)}_{clean(tissue)}"


class HistoryService:
    def __init__(
        self,
        history_dir: str | Path | None = None,
        *,
        paths: ProjectPaths | None = None,
    ) -> None:
        active_paths = paths or ProjectPaths.discover()
        self.history_dir = Path(history_dir or active_paths.history).resolve()
        self.index_path = self.history_dir / "history.sqlite3"
        self.reader = AnalysisBundleReader()
        self.writer = AnalysisBundleWriter()

    def _connect(self, path: Path | None = None) -> sqlite3.Connection:
        self.history_dir.mkdir(parents=True, exist_ok=True)
        connection = sqlite3.connect(path or self.index_path)
        connection.row_factory = sqlite3.Row
        connection.execute("""
            CREATE TABLE IF NOT EXISTS analyses (
                analysis_id TEXT PRIMARY KEY,
                target TEXT NOT NULL,
                tissue TEXT NOT NULL,
                created_at TEXT NOT NULL,
                analysis_status TEXT NOT NULL,
                sophiark_version TEXT NOT NULL,
                bundle_schema_version INTEGER NOT NULL,
                bundle_filename TEXT NOT NULL,
                bundle_checksum TEXT NOT NULL,
                content_fingerprint TEXT NOT NULL,
                classic_status TEXT NOT NULL,
                directed_status TEXT NOT NULL,
                evidence_available INTEGER NOT NULL,
                complex_context_available INTEGER NOT NULL
            )
        """)
        connection.commit()
        return connection

    @staticmethod
    def _entry(saved: SavedUnifiedResult, path: Path) -> HistoryEntry:
        manifest = saved.manifest
        context = dict(manifest.get("context_availability") or {})
        evidence = dict(context.get("evidence_provenance") or {}).get("status") == "AVAILABLE"
        complex_context = dict(context.get("corum_complex_context") or {}).get("status") == "AVAILABLE"
        statuses = dict(manifest.get("engine_statuses") or {})
        return HistoryEntry(
            analysis_id=str(manifest["analysis_id"]), target=str(manifest["target"]),
            tissue=str(manifest["tissue"]), created_at=str(manifest["created_at"]),
            analysis_status=str(manifest["analysis_status"]),
            sophiark_version=str(manifest.get("sophiark_version") or "UNKNOWN"),
            bundle_schema_version=int(manifest["bundle_schema_version"]),
            bundle_filename=path.name, bundle_checksum=_checksum(path),
            content_fingerprint=saved.content_fingerprint,
            classic_status=str(statuses.get("classic") or saved.report.classic_status),
            directed_status=str(statuses.get("directed") or saved.report.directed_status),
            evidence_available=evidence, complex_context_available=complex_context,
        )

    @staticmethod
    def _upsert(connection: sqlite3.Connection, entry: HistoryEntry) -> None:
        connection.execute(
            """INSERT OR REPLACE INTO analyses VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                entry.analysis_id, entry.target, entry.tissue, entry.created_at,
                entry.analysis_status, entry.sophiark_version, entry.bundle_schema_version,
                entry.bundle_filename, entry.bundle_checksum, entry.content_fingerprint,
                entry.classic_status, entry.directed_status, int(entry.evidence_available),
                int(entry.complex_context_available),
            ),
        )

    def save(self, report: UnifiedResearchReport) -> HistoryEntry:
        if report.status is UnifiedStatus.FAILED:
            raise ValueError("FAILED Unified analysis is not persisted as completed history")
        now = datetime.now().astimezone()
        identifier = str(uuid.uuid4())
        filename = f"{_safe_stem(report.target, report.tissue)}_{now:%Y%m%d}_{identifier[:8]}.sophiark"
        path = self.history_dir / filename
        self.writer.write(report, path, analysis_id=identifier)
        saved = self.reader.load(path)
        entry = self._entry(saved, path)
        try:
            with closing(self._connect()) as connection:
                self._upsert(connection, entry)
                connection.commit()
        except Exception:
            log.exception("bundle persisted but history indexing failed: id=%s", identifier)
            raise
        log.info("history bundle indexed: id=%s", identifier)
        return entry

    def open(self, bundle: str | Path, *, current_dataset_fingerprints: Mapping[str, str] | None = None) -> SavedUnifiedResult:
        return self.reader.load(bundle, current_dataset_fingerprints=current_dataset_fingerprints)

    def open_entry(self, analysis_id: str, *, current_dataset_fingerprints: Mapping[str, str] | None = None) -> SavedUnifiedResult:
        entries = {entry.analysis_id: entry for entry in self.list_entries()}
        if analysis_id not in entries:
            raise KeyError(f"analysis is not indexed: {analysis_id}")
        return self.open(self.history_dir / entries[analysis_id].bundle_filename, current_dataset_fingerprints=current_dataset_fingerprints)

    def list_entries(self) -> list[HistoryEntry]:
        if not self.index_path.exists():
            self.rebuild_history_index()
        with closing(self._connect()) as connection:
            rows = connection.execute("SELECT * FROM analyses ORDER BY created_at DESC,analysis_id").fetchall()
        return [HistoryEntry(
            analysis_id=row["analysis_id"], target=row["target"], tissue=row["tissue"],
            created_at=row["created_at"], analysis_status=row["analysis_status"],
            sophiark_version=row["sophiark_version"], bundle_schema_version=row["bundle_schema_version"],
            bundle_filename=row["bundle_filename"], bundle_checksum=row["bundle_checksum"],
            content_fingerprint=row["content_fingerprint"], classic_status=row["classic_status"],
            directed_status=row["directed_status"], evidence_available=bool(row["evidence_available"]),
            complex_context_available=bool(row["complex_context_available"]),
        ) for row in rows]

    def rebuild_history_index(self) -> RebuildReport:
        self.history_dir.mkdir(parents=True, exist_ok=True)
        temporary = self.history_dir / f".history.{uuid.uuid4().hex}.sqlite3.tmp"
        valid = indexed = corrupt = duplicates = 0
        corrupt_files: list[str] = []
        seen: set[str] = set()
        try:
            with closing(self._connect(temporary)) as connection:
                for path in sorted(self.history_dir.glob("*.sophiark")):
                    try:
                        saved = self.reader.load(path)
                        valid += 1
                        identifier = str(saved.manifest["analysis_id"])
                        if identifier in seen:
                            duplicates += 1
                            continue
                        seen.add(identifier)
                        self._upsert(connection, self._entry(saved, path))
                        indexed += 1
                    except (BundleError, OSError, ValueError) as exc:
                        corrupt += 1
                        corrupt_files.append(path.name)
                        log.warning("corrupt history bundle skipped: file=%s reason=%s", path.name, exc)
                connection.commit()
            os.replace(temporary, self.index_path)
        finally:
            temporary.unlink(missing_ok=True)
        log.info("history index rebuilt: valid=%d indexed=%d corrupt=%d duplicates=%d", valid, indexed, corrupt, duplicates)
        return RebuildReport(valid, indexed, corrupt, duplicates, tuple(corrupt_files))

    def import_bundle(self, source: str | Path) -> ImportResult:
        source_path = Path(source)
        saved = self.reader.load(source_path)
        return self._import_validated(saved, source_path.read_bytes())

    def import_bytes(self, data: bytes) -> ImportResult:
        saved = self.reader.load(data)
        return self._import_validated(saved, data)

    def _import_validated(self, saved: SavedUnifiedResult, data: bytes) -> ImportResult:
        identifier = str(saved.manifest["analysis_id"])
        existing = {entry.analysis_id: entry for entry in self.list_entries()}
        # Bundles, not the disposable index, are authoritative. Include any
        # valid unindexed file left behind by an earlier index failure.
        for path in sorted(self.history_dir.glob("*.sophiark")):
            try:
                candidate = self.reader.load(path)
                candidate_id = str(candidate.manifest["analysis_id"])
                if candidate_id not in existing:
                    existing[candidate_id] = self._entry(candidate, path)
            except BundleError:
                continue
        if identifier in existing:
            entry = existing[identifier]
            if entry.content_fingerprint == saved.content_fingerprint:
                return ImportResult(ImportStatus.ALREADY_IMPORTED, entry)
            return ImportResult(ImportStatus.ID_CONFLICT, entry, "same analysis_id has different content")
        filename = f"{_safe_stem(saved.report.target, saved.report.tissue)}_{identifier[:8]}.sophiark"
        final = self.history_dir / filename
        self.history_dir.mkdir(parents=True, exist_ok=True)
        temporary = final.with_name(f".{final.name}.{uuid.uuid4().hex}.tmp")
        try:
            temporary.write_bytes(data)
            self.reader.load(temporary)
            os.replace(temporary, final)
        finally:
            temporary.unlink(missing_ok=True)
        entry = self._entry(self.reader.load(final), final)
        with closing(self._connect()) as connection:
            self._upsert(connection, entry)
            connection.commit()
        log.info("external analysis bundle imported: id=%s", identifier)
        return ImportResult(ImportStatus.IMPORTED, entry)


def rebuild_history_index(history_dir: str | Path | None = None) -> RebuildReport:
    return HistoryService(history_dir).rebuild_history_index()
