"""Deterministic, source-specific data build contracts.

This module produces processed artifacts only. It never imports or executes a
calculation engine.
"""

from __future__ import annotations

import csv
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import tempfile
from typing import Mapping

from .entities import normalize_string_protein_id
from .fingerprints import sha256_file
from .status import ProductStatus
from .versions import PREPROCESSING_VERSIONS


class DataBuildError(ValueError):
    pass


@dataclass(frozen=True, slots=True)
class BuildReport:
    dataset: str
    detected_version: str
    status: ProductStatus
    source_records: int
    valid_records: int
    mapped_records: int
    unmapped_records: int
    usable_records: int
    dropped_records: int
    mapping_coverage: float
    schema_warnings: tuple[str, ...]
    previous_fingerprint: str | None
    new_fingerprint: str
    derived_artifacts_rebuilt: tuple[str, ...]
    affected_caches: tuple[str, ...]
    smoke_tests: tuple[str, ...]
    processing_version: str
    built_at_utc: str

    def as_dict(self) -> dict[str, object]:
        result = asdict(self)
        result["status"] = self.status.value
        return result


def _delimiter(path: Path) -> str:
    if path.suffix.casefold() in {".tsv", ".txt"}:
        return "\t"
    return ","


def _atomic_write_csv(path: Path, fieldnames: list[str], rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent, prefix=f".{path.name}.", suffix=".tmp", text=True
    )
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        os.replace(temporary_name, path)
    except Exception:
        try:
            os.unlink(temporary_name)
        except FileNotFoundError:
            pass
        raise


def write_build_report(report: BuildReport, path: str | Path) -> None:
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(
        json.dumps(report.as_dict(), ensure_ascii=False, indent=2, sort_keys=True),
        encoding="utf-8",
    )


def build_string_dataset(
    source: str | Path,
    destination: str | Path,
    *,
    dataset_version: str,
    taxon_id: int = 9606,
    canonical_mapping: Mapping[str, str] | None = None,
    built_at: datetime | None = None,
) -> BuildReport:
    """Validate and build a canonical STRING edge table.

    `combined_score` remains on its source 0..1000 scale. No graph, weight or
    ranking semantics are introduced here.
    """
    source_path = Path(source)
    output_path = Path(destination)
    if not source_path.is_file():
        raise DataBuildError(f"STRING source is missing: {source_path}")
    previous = sha256_file(output_path) if output_path.is_file() else None
    with source_path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter=_delimiter(source_path))
        fieldnames = list(reader.fieldnames or ())
        required = {"protein1", "protein2", "combined_score"}
        missing = sorted(required - set(fieldnames))
        if missing:
            raise DataBuildError("incompatible STRING schema; missing columns: " + ", ".join(missing))
        source_rows = list(reader)

    normalized: list[dict[str, str]] = []
    mapped = 0
    unmapped = 0
    seen: set[tuple[str, str]] = set()
    for line_number, row in enumerate(source_rows, start=2):
        try:
            first = normalize_string_protein_id(row["protein1"], taxon_id=taxon_id)
            second = normalize_string_protein_id(row["protein2"], taxon_id=taxon_id)
            score = int(row["combined_score"])
        except (TypeError, ValueError) as exc:
            raise DataBuildError(f"malformed STRING record at line {line_number}: {exc}") from exc
        if not 0 <= score <= 1000:
            raise DataBuildError(f"combined_score out of range at line {line_number}: {score}")
        if first == second:
            raise DataBuildError(f"self interaction at line {line_number}: {first}")
        key = tuple(sorted((first, second)))
        if key in seen:
            raise DataBuildError(f"duplicate STRING interaction at line {line_number}: {key}")
        seen.add(key)
        if canonical_mapping is not None:
            mapped_first = canonical_mapping.get(first)
            mapped_second = canonical_mapping.get(second)
            if mapped_first is None or mapped_second is None:
                unmapped += 1
                continue
            first, second = mapped_first, mapped_second
        mapped += 1
        output = dict(row)
        output["protein1"] = first
        output["protein2"] = second
        output["combined_score"] = str(score)
        normalized.append(output)

    normalized.sort(key=lambda item: (item["protein1"], item["protein2"], item["combined_score"]))
    _atomic_write_csv(output_path, fieldnames, normalized)
    current = sha256_file(output_path)
    count = len(source_rows)
    now = built_at or datetime.now(timezone.utc)
    return BuildReport(
        dataset="STRING",
        detected_version=dataset_version,
        status=ProductStatus.READY,
        source_records=count,
        valid_records=count,
        mapped_records=mapped,
        unmapped_records=unmapped,
        usable_records=mapped,
        dropped_records=unmapped,
        mapping_coverage=(mapped / count if count else 0.0),
        schema_warnings=(),
        previous_fingerprint=previous,
        new_fingerprint=current,
        derived_artifacts_rebuilt=(str(output_path),),
        affected_caches=("content-addressed cache keys change with the new fingerprint",),
        smoke_tests=("schema", "identifier_format", "taxon", "score_range", "duplicate_key", "determinism"),
        processing_version=PREPROCESSING_VERSIONS["string"],
        built_at_utc=now.astimezone(timezone.utc).isoformat(),
    )
