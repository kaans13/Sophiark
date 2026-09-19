"""Portable, untrusted-input-safe `.sophiark` bundle v1 reader and writer."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime, timezone
from enum import Enum
import hashlib
import io
import json
import logging
import math
import os
from pathlib import Path, PurePosixPath
import re
import uuid
import zipfile
from typing import Any, Mapping

import numpy as np
import pandas as pd

from src.unified.service import UnifiedResearchReport, UnifiedStatus
from .versions import BUNDLE_SCHEMA_VERSION, SOPHIARK_VERSION, UNIFIED_SCHEMA_VERSION


log = logging.getLogger(__name__)
MAX_MEMBERS = 32
MAX_UNCOMPRESSED_BYTES = 512 * 1024 * 1024
REQUIRED_MEMBERS = frozenset({
    "manifest.json", "integrity.json", "report.json",
    "tables/unified_candidates.json", "tables/classic_results.json",
    "tables/directed_results.json", "tables/complex_context.json",
})
TABLE_MEMBERS = {
    "candidates": "tables/unified_candidates.json",
    "classic_raw": "tables/classic_results.json",
    "directed_raw": "tables/directed_results.json",
    "complex_context": "tables/complex_context.json",
}


class BundleStatus(str, Enum):
    READY = "READY"
    BUNDLE_CORRUPT = "BUNDLE_CORRUPT"
    UNSUPPORTED_BUNDLE_VERSION = "UNSUPPORTED_BUNDLE_VERSION"


class BundleError(ValueError):
    def __init__(self, status: BundleStatus, message: str):
        super().__init__(message)
        self.status = status


@dataclass(frozen=True, slots=True)
class BundleWriteResult:
    path: Path
    analysis_id: str
    content_fingerprint: str
    bundle_checksum: str
    bytes_written: int


@dataclass(slots=True)
class SavedUnifiedResult:
    report: UnifiedResearchReport
    manifest: Mapping[str, Any]
    content_fingerprint: str
    historical_dataset_mismatches: Mapping[str, Mapping[str, str | None]]
    bundle_path: Path | None = None

    @property
    def is_historical(self) -> bool:
        return bool(self.historical_dataset_mismatches)


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _json_bytes(value: Any) -> bytes:
    return json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"),
        allow_nan=False,
    ).encode("utf-8")


def _portable(value: Any, *, key: str = "") -> Any:
    """Convert numpy/pandas values and remove machine-specific absolute paths."""
    if value is pd.NA:
        return {"__sophiark_scalar__": "pd.NA"}
    if value is pd.NaT:
        return {"__sophiark_scalar__": "pd.NaT"}
    if value is None:
        return None
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float) and not math.isfinite(value):
        return {"__sophiark_scalar__": "nan" if math.isnan(value) else "inf" if value > 0 else "-inf"}
    if isinstance(value, (datetime, pd.Timestamp)):
        return value.isoformat()
    if isinstance(value, Path):
        return value.name
    if isinstance(value, Mapping):
        return {str(name): _portable(item, key=str(name)) for name, item in value.items()}
    if isinstance(value, (list, tuple, set)):
        return [_portable(item) for item in value]
    if isinstance(value, str):
        candidate = value.strip()
        # Provenance producers use several names (path, source, file, index_path).
        # Never publish a machine-specific absolute filesystem location,
        # regardless of the surrounding metadata key.
        if Path(candidate).is_absolute() or re.match(r"^[A-Za-z]:[\\/]", candidate):
            return re.split(r"[\\/]", candidate)[-1]
        return value
    if isinstance(value, (bool, int, float)):
        return value
    return str(value)


def _decode_scalar(value: Any) -> Any:
    if isinstance(value, Mapping) and value.get("__sophiark_scalar__") in {
        "nan", "inf", "-inf", "pd.NA", "pd.NaT",
    }:
        token = value["__sophiark_scalar__"]
        return {
            "nan": float("nan"), "inf": float("inf"), "-inf": float("-inf"),
            "pd.NA": pd.NA, "pd.NaT": pd.NaT,
        }[token]
    return value


def _table_bytes(frame: pd.DataFrame) -> bytes:
    payload = {
        "schema_version": 1,
        "columns": [str(column) for column in frame.columns],
        "dtypes": {str(column): str(frame[column].dtype) for column in frame.columns},
        "data": {
            str(column): [_portable(value, key=str(column)) for value in frame[column].tolist()]
            for column in frame.columns
        },
    }
    return _json_bytes(payload)


def _table_from_bytes(data: bytes) -> pd.DataFrame:
    payload = json.loads(data)
    if payload.get("schema_version") != 1 or not isinstance(payload.get("columns"), list):
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "invalid table schema")
    columns = payload["columns"]
    values = payload.get("data")
    if not isinstance(values, dict) or set(values) != set(columns):
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "table columns and data disagree")
    lengths = {len(values[column]) for column in columns}
    if len(lengths) > 1:
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "table columns have different lengths")
    decoded = {column: [_decode_scalar(value) for value in values[column]] for column in columns}
    frame = pd.DataFrame(decoded)
    for column, dtype in payload.get("dtypes", {}).items():
        try:
            if dtype.startswith("datetime64"):
                frame[column] = pd.to_datetime(frame[column])
            elif dtype == "category":
                frame[column] = frame[column].astype("category")
            elif dtype == "object":
                # Construction can infer a narrower dtype when an object
                # column happens to contain only numbers in this analysis.
                # Rebuild from decoded values as well, since inference can
                # turn an object-column None into float NaN.
                frame[column] = pd.Series(decoded[column], dtype="object")
            else:
                frame[column] = frame[column].astype(dtype)
        except (TypeError, ValueError):
            raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"table dtype is incompatible: {column}={dtype}")
    return frame


def _safe_member_names(infos: list[zipfile.ZipInfo]) -> None:
    if len(infos) > MAX_MEMBERS:
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "bundle has too many archive members")
    names = [info.filename for info in infos]
    if len(names) != len(set(names)):
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "bundle contains duplicate archive names")
    if sum(info.file_size for info in infos) > MAX_UNCOMPRESSED_BYTES:
        raise BundleError(BundleStatus.BUNDLE_CORRUPT, "bundle exceeds the uncompressed size limit")
    for name in names:
        path = PurePosixPath(name)
        if not name or name.startswith(("/", "\\")) or "\\" in name or ".." in path.parts or ":" in path.parts[0]:
            raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"unsafe archive member: {name}")


def _manifest(report: UnifiedResearchReport, analysis_id: str, created_at: str) -> dict[str, Any]:
    provenance = _portable(dict(report.provenance))
    run_manifest = dict(provenance.get("run_manifest") or {})
    datasets = dict(run_manifest.get("datasets") or {})
    return {
        "bundle_schema_version": BUNDLE_SCHEMA_VERSION,
        "analysis_id": analysis_id,
        "created_at": created_at,
        "target": report.target,
        "tissue": report.tissue,
        "analysis_status": report.status.value,
        "sophiark_version": run_manifest.get("sophiark_version", SOPHIARK_VERSION),
        "unified_schema_version": UNIFIED_SCHEMA_VERSION,
        "engine_versions": run_manifest.get("engine_versions", {
            "classic": provenance.get("classic_engine_version"),
            "directed": provenance.get("directed_engine_version"),
            "unified": provenance.get("unified_engine_version"),
        }),
        "engine_config_fingerprints": run_manifest.get("engine_config_fingerprints", {
            "classic": provenance.get("classic_config_fingerprint"),
            "directed": provenance.get("directed_config_fingerprint"),
        }),
        "engine_statuses": {"classic": report.classic_status, "directed": report.directed_status},
        "datasets": datasets,
        "preprocessing_versions": run_manifest.get("preprocessing_versions", {}),
        "capabilities": run_manifest.get("availability", {}),
        "export_schema_version": run_manifest.get("export_schema_version"),
        "result_table_schema_version": 1,
        "evidence_weighting_status": provenance.get("evidence_weighting_status", "BIAS_WARNING"),
        "context_availability": _portable(dict(report.context_availability)),
        "legacy_data_provenance": {
            "classic_hpa": {
                "classification": ["LEGACY_PARTIALLY_REPRODUCIBLE", "SCIENTIFIC_MIGRATION_REQUIRED"],
                "dataset": datasets.get("hpa_human_tissue"),
            },
            "classic_corum": {
                "classification": "LEGACY_REPRODUCIBLE",
                "dataset": datasets.get("corum_human_legacy_lookup"),
            },
        },
    }


class AnalysisBundleWriter:
    def write(
        self,
        report: UnifiedResearchReport,
        destination: str | Path,
        *,
        analysis_id: str | None = None,
        created_at: datetime | None = None,
    ) -> BundleWriteResult:
        if report.status is UnifiedStatus.FAILED:
            raise ValueError("FAILED analyses are not saved as completed history bundles")
        identifier = str(uuid.UUID(analysis_id)) if analysis_id else str(uuid.uuid4())
        timestamp = (created_at or datetime.now(timezone.utc)).astimezone(timezone.utc).isoformat()
        manifest = _manifest(report, identifier, timestamp)
        report_payload = {
            "target": report.target, "tissue": report.tissue, "status": report.status.value,
            "classic_status": report.classic_status, "directed_status": report.directed_status,
            "agreement": _portable(dict(report.agreement)),
            "provenance": _portable(dict(report.provenance)),
            "context_availability": _portable(dict(report.context_availability)),
            "errors": _portable(dict(report.errors)),
        }
        members: dict[str, bytes] = {
            "report.json": _json_bytes(report_payload),
            **{name: _table_bytes(getattr(report, attribute)) for attribute, name in TABLE_MEMBERS.items()},
        }
        identity_payload = {"manifest": manifest, "members": {name: _sha256(data) for name, data in members.items()}}
        content_fingerprint = _sha256(_json_bytes(identity_payload))
        manifest["content_fingerprint"] = content_fingerprint
        members["manifest.json"] = _json_bytes(manifest)
        integrity = {
            "algorithm": "SHA-256",
            "files": {name: _sha256(data) for name, data in sorted(members.items())},
        }
        members["integrity.json"] = _json_bytes(integrity)
        final = Path(destination)
        if final.suffix.lower() != ".sophiark":
            raise ValueError("analysis bundle destination must use .sophiark extension")
        final.parent.mkdir(parents=True, exist_ok=True)
        temporary = final.with_name(f".{final.name}.{uuid.uuid4().hex}.tmp")
        try:
            with zipfile.ZipFile(temporary, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
                for name, data in sorted(members.items()):
                    archive.writestr(name, data)
            AnalysisBundleReader().load(temporary)
            os.replace(temporary, final)
        finally:
            temporary.unlink(missing_ok=True)
        bundle_data = final.read_bytes()
        log.info("analysis bundle saved: id=%s bytes=%d", identifier, len(bundle_data))
        return BundleWriteResult(
            final, identifier, content_fingerprint, _sha256(bundle_data), len(bundle_data),
        )


class AnalysisBundleReader:
    def load(
        self,
        source: str | Path | bytes,
        *,
        current_dataset_fingerprints: Mapping[str, str] | None = None,
    ) -> SavedUnifiedResult:
        path = None if isinstance(source, bytes) else Path(source)
        stream: str | Path | io.BytesIO = io.BytesIO(source) if isinstance(source, bytes) else path
        try:
            with zipfile.ZipFile(stream, "r") as archive:
                infos = archive.infolist()
                _safe_member_names(infos)
                names = {info.filename for info in infos}
                missing = REQUIRED_MEMBERS - names
                if missing:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"missing required members: {sorted(missing)}")
                manifest_bytes = archive.read("manifest.json")
                try:
                    manifest = json.loads(manifest_bytes)
                except (UnicodeDecodeError, json.JSONDecodeError) as exc:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"invalid manifest: {exc}") from exc
                version = manifest.get("bundle_schema_version")
                if version != BUNDLE_SCHEMA_VERSION:
                    raise BundleError(
                        BundleStatus.UNSUPPORTED_BUNDLE_VERSION,
                        f"unsupported bundle schema version: {version}; supported: {BUNDLE_SCHEMA_VERSION}",
                    )
                try:
                    uuid.UUID(str(manifest["analysis_id"]))
                except (KeyError, ValueError, TypeError) as exc:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, "invalid analysis_id") from exc
                try:
                    integrity = json.loads(archive.read("integrity.json"))
                    expected = integrity["files"]
                except (UnicodeDecodeError, json.JSONDecodeError, KeyError, TypeError) as exc:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"invalid integrity manifest: {exc}") from exc
                listed = names - {"integrity.json"}
                if set(expected) != listed:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, "integrity member list does not match archive")
                for name, digest in expected.items():
                    if _sha256(archive.read(name)) != digest:
                        raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"integrity mismatch: {name}")
                try:
                    payload = json.loads(archive.read("report.json"))
                    tables = {
                        attribute: _table_from_bytes(archive.read(name))
                        for attribute, name in TABLE_MEMBERS.items()
                    }
                    status = UnifiedStatus(payload["status"])
                except BundleError:
                    raise
                except (UnicodeDecodeError, json.JSONDecodeError, KeyError, ValueError, TypeError) as exc:
                    raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"invalid saved result: {exc}") from exc
        except BundleError:
            raise
        except (OSError, zipfile.BadZipFile, RuntimeError) as exc:
            raise BundleError(BundleStatus.BUNDLE_CORRUPT, f"invalid analysis archive: {exc}") from exc
        report = UnifiedResearchReport(
            target=str(payload["target"]), tissue=str(payload["tissue"]), status=status,
            classic_status=str(payload["classic_status"]), directed_status=str(payload["directed_status"]),
            candidates=tables["candidates"], classic_raw=tables["classic_raw"],
            directed_raw=tables["directed_raw"], complex_context=tables["complex_context"],
            agreement=dict(payload.get("agreement") or {}), provenance=dict(payload.get("provenance") or {}),
            context_availability=dict(payload.get("context_availability") or {}),
            errors=dict(payload.get("errors") or {}),
        )
        saved = {
            name: details.get("fingerprint")
            for name, details in dict(manifest.get("datasets") or {}).items()
            if isinstance(details, Mapping)
        }
        current = dict(current_dataset_fingerprints or {})
        mismatches = {
            name: {"saved": fingerprint, "current": current.get(name)}
            for name, fingerprint in saved.items()
            if name in current and current[name] != fingerprint
        }
        log.info("analysis bundle loaded: id=%s historical=%s", manifest["analysis_id"], bool(mismatches))
        return SavedUnifiedResult(
            report=report, manifest=manifest,
            content_fingerprint=str(manifest.get("content_fingerprint") or ""),
            historical_dataset_mismatches=mismatches, bundle_path=path,
        )


def migrate_to_current_schema(manifest: Mapping[str, Any]) -> Mapping[str, Any]:
    """Explicit migration boundary; v1 is current and needs no migration."""
    if manifest.get("bundle_schema_version") != BUNDLE_SCHEMA_VERSION:
        raise BundleError(BundleStatus.UNSUPPORTED_BUNDLE_VERSION, "no migration is available")
    return manifest
