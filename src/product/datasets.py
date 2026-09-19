"""Central dataset catalog, schema checks and health reporting."""

from __future__ import annotations

from dataclasses import dataclass, field
import csv
import json
from pathlib import Path
import sqlite3
from typing import Iterable, Mapping

from .fingerprints import sha256_file
from .paths import ProjectPaths
from .status import ProductStatus


@dataclass(frozen=True, slots=True)
class DatasetSpec:
    name: str
    source: str
    version: str
    organism: str
    roles: tuple[str, ...]
    required: bool
    raw_location: str | None = None
    processed_location: str | None = None
    schema_version: str = "1"
    processing_version: str = "not-applicable"
    expected_columns: tuple[str, ...] = ()
    sqlite_tables: tuple[str, ...] = ()
    fingerprint_policy: str = "sha256-content"
    capabilities: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for field_name in ("name", "source", "version", "organism"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"dataset {field_name} must be non-empty text")
        if not isinstance(self.required, bool):
            raise TypeError("dataset required must be a boolean")
        if not self.roles or not all(isinstance(role, str) and role for role in self.roles):
            raise ValueError("dataset roles must contain non-empty role names")
        if self.location is None:
            raise ValueError("dataset must declare a raw or processed location")

    @property
    def location(self) -> str | None:
        return self.processed_location or self.raw_location


@dataclass(frozen=True, slots=True)
class DatasetHealth:
    name: str
    status: ProductStatus
    required: bool
    version: str
    path: str | None
    records: int | None = None
    mapped: int | None = None
    unmapped: int | None = None
    usable: int | None = None
    fingerprint: str | None = None
    warnings: tuple[str, ...] = ()

    def as_dict(self) -> dict[str, object]:
        return {
            "name": self.name, "status": self.status.value, "required": self.required,
            "version": self.version, "path": self.path, "records": self.records,
            "mapped": self.mapped, "unmapped": self.unmapped, "usable": self.usable,
            "fingerprint": self.fingerprint, "warnings": list(self.warnings),
        }


class DataRegistry:
    def __init__(self, root: str | Path, datasets: Iterable[DatasetSpec]):
        self.root = Path(root).resolve()
        items = tuple(datasets)
        names = [item.name for item in items]
        if len(names) != len(set(names)):
            raise ValueError("dataset names must be unique")
        self._datasets = {item.name: item for item in items}

    @classmethod
    def from_json(cls, path: str | Path, *, root: str | Path) -> "DataRegistry":
        payload = json.loads(Path(path).read_text(encoding="utf-8"))
        specs = []
        for raw in payload.get("datasets", []):
            item = dict(raw)
            item["roles"] = tuple(item.get("roles", ()))
            item["expected_columns"] = tuple(item.get("expected_columns", ()))
            item["sqlite_tables"] = tuple(item.get("sqlite_tables", ()))
            item["capabilities"] = tuple(item.get("capabilities", ()))
            specs.append(DatasetSpec(**item))
        return cls(root, specs)

    def names(self) -> tuple[str, ...]:
        return tuple(self._datasets)

    def get(self, name: str) -> DatasetSpec:
        return self._datasets[name]

    def resolve(self, spec: DatasetSpec) -> Path | None:
        if not spec.location:
            return None
        path = Path(spec.location)
        return path if path.is_absolute() else self.root / path

    def _schema_warnings(self, spec: DatasetSpec, path: Path) -> tuple[str, ...]:
        warnings: list[str] = []
        suffix = path.suffix.casefold()
        if spec.sqlite_tables and suffix in {".db", ".sqlite", ".sqlite3"}:
            try:
                with sqlite3.connect(f"file:{path.as_posix()}?mode=ro", uri=True) as connection:
                    tables = {row[0] for row in connection.execute(
                        "SELECT name FROM sqlite_master WHERE type='table'"
                    )}
                missing = sorted(set(spec.sqlite_tables) - tables)
                if missing:
                    warnings.append("missing SQLite tables: " + ", ".join(missing))
            except sqlite3.Error as exc:
                warnings.append(f"SQLite validation failed: {exc}")
        elif spec.expected_columns and suffix in {".csv", ".tsv", ".txt"}:
            delimiter = "\t" if suffix in {".tsv", ".txt"} else ","
            try:
                with path.open("r", encoding="utf-8-sig", newline="") as handle:
                    columns = next(csv.reader(handle, delimiter=delimiter))
                missing = sorted(set(spec.expected_columns) - set(columns))
                if missing:
                    warnings.append("missing columns: " + ", ".join(missing))
            except (OSError, UnicodeError, StopIteration) as exc:
                warnings.append(f"schema inspection failed: {exc}")
        return tuple(warnings)

    def inspect(self, name: str, *, fingerprint: bool = False) -> DatasetHealth:
        spec = self.get(name)
        path = self.resolve(spec)
        if path is None or not path.is_file():
            return DatasetHealth(
                name, ProductStatus.DATA_UNAVAILABLE, spec.required, spec.version,
                str(path) if path else None, warnings=("dataset file is missing",),
            )
        warnings = self._schema_warnings(spec, path)
        status = ProductStatus.INCOMPATIBLE_SCHEMA if warnings else ProductStatus.READY
        return DatasetHealth(
            name, status, spec.required, spec.version, str(path),
            fingerprint=sha256_file(path) if fingerprint else None, warnings=warnings,
        )

    def health_report(self, *, fingerprint: bool = False) -> dict[str, DatasetHealth]:
        return {name: self.inspect(name, fingerprint=fingerprint) for name in self.names()}

    def overall_status(self) -> ProductStatus:
        health = self.health_report()
        if any(item.required and item.status is not ProductStatus.READY for item in health.values()):
            return ProductStatus.DATA_UNAVAILABLE
        if any(item.status is not ProductStatus.READY for item in health.values()):
            return ProductStatus.PARTIAL
        return ProductStatus.READY


def default_registry(paths: ProjectPaths | None = None) -> DataRegistry:
    active_paths = paths or ProjectPaths.discover()
    return DataRegistry.from_json(active_paths.config / "datasets.json", root=active_paths.root)


def inspect_build_report(path: str | Path) -> DatasetHealth:
    """Translate a shadow builder report into the shared health vocabulary."""
    report_path = Path(path)
    if not report_path.is_file():
        return DatasetHealth(
            report_path.parent.name or "unknown", ProductStatus.BUILD_REQUIRED,
            False, "UNKNOWN", str(report_path), warnings=("build report is missing",),
        )
    try:
        payload = json.loads(report_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return DatasetHealth(
            report_path.parent.name or "unknown", ProductStatus.BUILD_FAILED,
            False, "UNKNOWN", str(report_path), warnings=(f"invalid build report: {exc}",),
        )
    required = {
        "dataset_id", "source_version", "processing_version", "schema_version",
        "source_records", "mapped_records", "unmapped_records", "usable_records",
        "processed_fingerprint", "output_artifacts", "status",
    }
    missing = sorted(required - set(payload))
    if missing:
        return DatasetHealth(
            str(payload.get("dataset_id") or report_path.parent.name),
            ProductStatus.INCOMPATIBLE_SCHEMA, False,
            str(payload.get("source_version") or "UNKNOWN"), str(report_path),
            warnings=("build report missing fields: " + ", ".join(missing),),
        )
    artifacts = [Path(item) for item in payload["output_artifacts"]]
    absent = [str(item) for item in artifacts if not item.is_file()]
    errors = tuple(str(item) for item in payload.get("errors", ()))
    declared = str(payload["status"])
    if declared != ProductStatus.READY.value or errors or absent:
        warnings = (*tuple(str(item) for item in payload.get("warnings", ())), *errors)
        if absent:
            warnings = (*warnings, "missing build artifacts: " + ", ".join(absent))
        status = ProductStatus.BUILD_FAILED
    else:
        warnings = tuple(str(item) for item in payload.get("warnings", ()))
        status = ProductStatus.READY
    return DatasetHealth(
        name=str(payload["dataset_id"]), status=status, required=False,
        version=str(payload["source_version"]), path=str(report_path),
        records=int(payload["source_records"]), mapped=int(payload["mapped_records"]),
        unmapped=int(payload["unmapped_records"]), usable=int(payload["usable_records"]),
        fingerprint=str(payload["processed_fingerprint"]), warnings=warnings,
    )
