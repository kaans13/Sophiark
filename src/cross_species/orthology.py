from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import pandas as pd


class OrthologyClass(str, Enum):
    ONE_TO_ONE = "ONE_TO_ONE"
    ONE_TO_MANY = "ONE_TO_MANY"
    MANY_TO_MANY = "MANY_TO_MANY"
    NO_ORTHOLOG = "NO_ORTHOLOG"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class OrthologRecord:
    human_gene: str
    mouse_gene: str | None
    human_ensembl_gene_id: str | None
    mouse_ensembl_gene_id: str | None
    orthology_class: OrthologyClass
    source: str
    source_snapshot: str
    confidence: int | None = None


_EXPECTED_COLUMNS = {
    "Gene stable ID": "human_id",
    "Gene name": "human_gene",
    "Mouse gene stable ID": "mouse_id",
    "Mouse gene name": "mouse_gene",
    "Mouse homology type": "homology_type",
    "Mouse orthology confidence [0 low, 1 high]": "confidence",
}


class OrthologyIndex:
    """Immutable, provenance-checked Human→Mouse Ensembl Compara index."""

    def __init__(self, records: dict[str, tuple[OrthologRecord, ...]], metadata: dict):
        self._records = records
        self.metadata = dict(metadata)

    @classmethod
    def from_snapshot(cls, tsv_path: Path, metadata_path: Path) -> "OrthologyIndex":
        content = tsv_path.read_bytes()
        metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
        actual_sha = hashlib.sha256(content).hexdigest()
        if actual_sha != metadata.get("sha256"):
            raise ValueError("ORTHOLOGY_SNAPSHOT_CHECKSUM_MISMATCH")
        if int(metadata.get("human_taxon_id", 0)) != 9606 or int(metadata.get("mouse_taxon_id", 0)) != 10090:
            raise ValueError("ORTHOLOGY_SNAPSHOT_SPECIES_MISMATCH")

        frame = pd.read_csv(tsv_path, sep="\t", dtype=str).rename(columns=_EXPECTED_COLUMNS)
        missing = set(_EXPECTED_COLUMNS.values()) - set(frame.columns)
        if missing:
            raise ValueError(f"ORTHOLOGY_SNAPSHOT_COLUMNS_MISSING: {sorted(missing)}")
        source = str(metadata.get("source", "Ensembl Compara"))
        snapshot = actual_sha
        grouped: dict[str, list[OrthologRecord]] = {}
        def clean(value) -> str:
            if pd.isna(value):
                return ""
            return str(value).strip()

        for row in frame.itertuples(index=False):
            human_gene = clean(getattr(row, "human_gene", ""))
            if not human_gene:
                continue
            mouse_gene = clean(getattr(row, "mouse_gene", "")) or None
            raw_type = clean(getattr(row, "homology_type", ""))
            if not mouse_gene:
                cls_value = OrthologyClass.NO_ORTHOLOG
            elif raw_type == "ortholog_one2one":
                cls_value = OrthologyClass.ONE_TO_ONE
            elif raw_type == "ortholog_one2many":
                cls_value = OrthologyClass.ONE_TO_MANY
            elif raw_type == "ortholog_many2many":
                cls_value = OrthologyClass.MANY_TO_MANY
            else:
                cls_value = OrthologyClass.UNRESOLVED
            confidence_raw = clean(getattr(row, "confidence", ""))
            record = OrthologRecord(
                human_gene=human_gene,
                mouse_gene=mouse_gene,
                human_ensembl_gene_id=clean(getattr(row, "human_id", "")) or None,
                mouse_ensembl_gene_id=clean(getattr(row, "mouse_id", "")) or None,
                orthology_class=cls_value,
                source=source,
                source_snapshot=snapshot,
                confidence=int(float(confidence_raw)) if confidence_raw else None,
            )
            grouped.setdefault(human_gene.casefold(), []).append(record)

        records: dict[str, tuple[OrthologRecord, ...]] = {}
        for key, values in grouped.items():
            mapped = [item for item in values if item.mouse_gene]
            if len({item.mouse_gene for item in mapped}) > 1:
                adjusted = []
                for item in values:
                    inferred = (
                        OrthologyClass.MANY_TO_MANY
                        if item.orthology_class is OrthologyClass.MANY_TO_MANY
                        else OrthologyClass.ONE_TO_MANY
                    )
                    adjusted.append(OrthologRecord(**{**item.__dict__, "orthology_class": inferred}))
                values = adjusted
            records[key] = tuple(values)
        return cls(records, metadata)

    def lookup(self, human_gene: str) -> tuple[OrthologRecord, ...]:
        query = str(human_gene or "").strip().casefold()
        return self._records.get(query, ())

    def one_to_one(self, human_gene: str) -> OrthologRecord | None:
        records = [r for r in self.lookup(human_gene) if r.orthology_class is OrthologyClass.ONE_TO_ONE]
        return records[0] if len(records) == 1 else None

    def statistics(self) -> dict[str, int]:
        counts = {item.value: 0 for item in OrthologyClass}
        for records in self._records.values():
            classes = {record.orthology_class for record in records}
            for item in classes:
                counts[item.value] += 1
        counts["human_symbol_count"] = len(self._records)
        counts["relationship_count"] = sum(len(items) for items in self._records.values())
        return counts
