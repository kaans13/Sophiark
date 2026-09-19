"""Small shared lifecycle helpers for deterministic shadow builds."""

from __future__ import annotations

from dataclasses import asdict, dataclass, field
import json
from pathlib import Path
from typing import Any, Mapping

from src.product.fingerprints import sha256_file, stable_fingerprint
from src.product.paths import ProjectPaths
from src.product.status import ProductStatus


class SourceBuildError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class SourceBuildReport:
    dataset_id: str
    source_version: str
    processing_version: str
    schema_version: str
    source_records: int
    valid_records: int
    mapped_records: int
    unmapped_records: int
    usable_records: int
    unique_source_entities: int
    mapped_entities: int
    unmapped_entities: int
    mapping_coverage: float
    raw_fingerprint: str
    processed_fingerprint: str
    output_artifacts: tuple[str, ...]
    status: ProductStatus
    warnings: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    details: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        counts = (
            self.source_records, self.valid_records, self.mapped_records,
            self.unmapped_records, self.usable_records, self.unique_source_entities,
            self.mapped_entities, self.unmapped_entities,
        )
        if any(value < 0 for value in counts):
            raise ValueError("build report counts cannot be negative")
        if self.valid_records > self.source_records:
            raise ValueError("valid_records cannot exceed source_records")
        if self.mapped_records + self.unmapped_records != self.valid_records:
            raise ValueError("mapped + unmapped records must equal valid records")
        if self.usable_records > self.valid_records:
            raise ValueError("usable_records cannot exceed valid_records")
        if self.mapped_entities + self.unmapped_entities != self.unique_source_entities:
            raise ValueError("mapped + unmapped entities must equal unique source entities")
        expected = self.mapped_entities / self.unique_source_entities if self.unique_source_entities else 0.0
        if abs(self.mapping_coverage - expected) > 1e-12:
            raise ValueError("mapping coverage is inconsistent with entity counts")
        if not 0.0 <= self.mapping_coverage <= 1.0:
            raise ValueError("mapping coverage must be in [0, 1]")

    def as_dict(self) -> dict[str, Any]:
        payload = asdict(self)
        payload["status"] = self.status.value
        return payload


def raw_identity(path: str | Path) -> str:
    item = Path(path)
    if not item.is_file():
        raise SourceBuildError(f"raw source is missing: {item}")
    return sha256_file(item)


def shadow_directory(
    *, dataset_id: str, raw_fingerprint: str, processing_version: str,
    output_root: str | Path | None = None, parameters: Mapping[str, Any] | None = None,
) -> Path:
    root = Path(output_root) if output_root is not None else ProjectPaths.discover().derived / "staging"
    build_id = stable_fingerprint({
        "dataset_id": dataset_id,
        "raw_fingerprint": raw_fingerprint,
        "processing_version": processing_version,
        "parameters": dict(parameters or {}),
    })[:16]
    destination = root / dataset_id / build_id
    destination.mkdir(parents=True, exist_ok=True)
    return destination


def processed_identity(
    *, normalized_content_fingerprint: str, source_version: str,
    processing_version: str, schema_version: str,
    parameters: Mapping[str, Any] | None = None,
) -> str:
    return stable_fingerprint({
        "normalized_content": normalized_content_fingerprint,
        "source_version": source_version,
        "processing_version": processing_version,
        "schema_version": schema_version,
        "parameters": dict(parameters or {}),
    })


def write_report(report: SourceBuildReport, destination: str | Path | None = None) -> Path:
    path = Path(destination) if destination is not None else Path(report.output_artifacts[0]).parent / "build-report.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        json.dumps(report.as_dict(), ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )
    return path
