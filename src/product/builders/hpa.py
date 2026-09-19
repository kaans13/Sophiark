"""HPA TPM shadow builder preserving the runtime's max-per-tissue contract."""

from __future__ import annotations

import os
from pathlib import Path
import sqlite3
import tempfile

import pandas as pd

from src.product.fingerprints import sha256_file
from src.product.status import ProductStatus
from src.product.versions import PREPROCESSING_VERSIONS

from .common import (
    SourceBuildError, SourceBuildReport, processed_identity, raw_identity,
    shadow_directory, write_report,
)


DATASET_ID = "hpa_human_tissue"
SCHEMA_VERSION = "hpa-tissue-expression-1"
REQUIRED_COLUMNS = {"Gene", "Tissue", "TPM"}


def _mapping(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str, usecols=["ENSP", "ENSG"])
    if frame["ENSG"].duplicated().any():
        raise SourceBuildError("ENSG mapping is ambiguous; duplicate ENSG entries found")
    return dict(zip(frame["ENSG"].str.strip(), frame["ENSP"].str.strip()))


def _write_sqlite(frame: pd.DataFrame, destination: Path) -> None:
    handle, temporary_name = tempfile.mkstemp(
        prefix=f".{destination.stem}-", suffix=".sqlite.tmp", dir=destination.parent,
    )
    os.close(handle)
    temporary = Path(temporary_name)
    try:
        connection = sqlite3.connect(temporary)
        try:
            frame.to_sql("tissue_expression", connection, if_exists="replace", index=False)
            connection.execute("CREATE INDEX idx_hpa_protein ON tissue_expression(protein_id)")
            connection.execute("CREATE INDEX idx_hpa_tissue ON tissue_expression(tissue)")
            connection.commit()
        finally:
            connection.close()
        os.replace(temporary, destination)
    finally:
        if temporary.exists():
            temporary.unlink()


def build_hpa_shadow(
    source: str | Path,
    *,
    mapping_path: str | Path,
    output_root: str | Path | None = None,
    source_version: str = "UNKNOWN",
) -> SourceBuildReport:
    source_path = Path(source)
    raw_fingerprint = raw_identity(source_path)
    try:
        raw = pd.read_csv(source_path, sep="\t", low_memory=False)
    except Exception as exc:
        raise SourceBuildError(f"HPA source could not be read: {exc}") from exc
    missing = sorted(REQUIRED_COLUMNS - set(raw.columns))
    if missing:
        raise SourceBuildError("incompatible HPA schema; missing columns: " + ", ".join(missing))
    source_records = len(raw)
    work = raw[["Gene", "Tissue", "TPM"]].copy()
    work["Gene"] = work["Gene"].astype("string").str.strip()
    work["Tissue"] = work["Tissue"].astype("string").str.strip()
    work["TPM"] = pd.to_numeric(work["TPM"], errors="coerce")
    valid_mask = work["Gene"].str.fullmatch(r"ENSG\d+") & work["Tissue"].ne("") & work["TPM"].notna() & work["TPM"].ge(0)
    valid = work.loc[valid_mask].copy()
    if valid.empty:
        raise SourceBuildError("HPA source contains no valid tissue-expression records")
    valid = valid.groupby(["Gene", "Tissue"], as_index=False, sort=True)["TPM"].max()
    entity_map = _mapping(Path(mapping_path))
    unique_entities = set(valid["Gene"].astype(str))
    mapped_entity_set = unique_entities & set(entity_map)
    valid["protein_id"] = valid["Gene"].map(entity_map)
    mapped = valid[valid["protein_id"].notna()].copy()
    mapped["tissue"] = mapped["Tissue"].str.title()
    mapped["expression_level"] = mapped["TPM"].astype(float)
    normalized = mapped[["protein_id", "tissue", "expression_level"]].sort_values(
        ["protein_id", "tissue"], kind="mergesort"
    ).reset_index(drop=True)
    directory = shadow_directory(
        dataset_id=DATASET_ID, raw_fingerprint=raw_fingerprint,
        processing_version=PREPROCESSING_VERSIONS["hpa"], output_root=output_root,
        parameters={"taxon_id": 9606, "value_field": "TPM", "duplicate_policy": "max"},
    )
    csv_path = directory / "tissue_expression.csv"
    sqlite_path = directory / "tissue_expression.sqlite"
    normalized.to_csv(csv_path, index=False, encoding="utf-8", lineterminator="\n")
    _write_sqlite(normalized, sqlite_path)
    content_fingerprint = sha256_file(csv_path)
    processed_fingerprint = processed_identity(
        normalized_content_fingerprint=content_fingerprint,
        source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["hpa"],
        schema_version=SCHEMA_VERSION,
        parameters={"taxon_id": 9606, "value_field": "TPM", "duplicate_policy": "max"},
    )
    report = SourceBuildReport(
        dataset_id=DATASET_ID, source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["hpa"], schema_version=SCHEMA_VERSION,
        source_records=source_records, valid_records=len(valid), mapped_records=len(mapped),
        unmapped_records=len(valid) - len(mapped), usable_records=len(normalized),
        unique_source_entities=len(unique_entities), mapped_entities=len(mapped_entity_set),
        unmapped_entities=len(unique_entities - mapped_entity_set),
        mapping_coverage=len(mapped_entity_set) / len(unique_entities) if unique_entities else 0.0,
        raw_fingerprint=raw_fingerprint, processed_fingerprint=processed_fingerprint,
        output_artifacts=(str(csv_path), str(sqlite_path)), status=ProductStatus.READY,
        warnings=(), details={
            "organism": "Homo sapiens", "taxon_id": 9606,
            "source_tissues": int(valid["Tissue"].nunique()),
            "normalized_tissues": int(normalized["tissue"].nunique()),
            "value_semantics": "HPA TPM; max duplicate per ENSG/tissue; no threshold or multiplier change",
        },
    )
    write_report(report)
    return report


def audit_hpa_shadow(production_db: str | Path, shadow_db: str | Path) -> dict[str, int | float | None]:
    connection = sqlite3.connect(f"file:{Path(production_db).resolve().as_posix()}?mode=ro", uri=True)
    try:
        connection.execute("ATTACH DATABASE ? AS shadow", (str(Path(shadow_db).resolve()),))
        row = connection.execute("""
            WITH production AS (
                SELECT protein_id, tissue, MAX(expression_level) AS expression_level
                FROM tissue_expression GROUP BY protein_id, tissue
            )
            SELECT COUNT(*) AS shadow_rows,
                   SUM(CASE WHEN p.protein_id IS NOT NULL THEN 1 ELSE 0 END) AS matched_rows,
                   SUM(CASE WHEN p.expression_level = s.expression_level THEN 1 ELSE 0 END) AS exact_rows,
                   MAX(ABS(p.expression_level - s.expression_level)) AS max_delta
            FROM shadow.tissue_expression s
            LEFT JOIN production p USING(protein_id, tissue)
        """).fetchone()
    finally:
        connection.close()
    return {"shadow_rows": row[0], "matched_rows": row[1], "exact_rows": row[2], "max_delta": row[3]}
