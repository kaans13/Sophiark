"""Safe construction of immutable snapshots from explicit result objects."""

from __future__ import annotations

from dataclasses import asdict, is_dataclass
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from .models import (
    FrozenResultTable,
    ProvenanceKind,
    ProvenanceRecord,
    SimulationResultSnapshot,
)


RESULT_SCHEMA_VERSION = 1
_KNOWN_TAXA = {
    "human": 9606,
    "homo sapiens": 9606,
    "insan": 9606,
    "insan (homo sapiens)": 9606,
    "mouse": 10090,
    "mus musculus": 10090,
    "fare": 10090,
    "fare (mus musculus)": 10090,
}


def _canonical_json(value: Any) -> str:
    return json.dumps(
        value,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    )


def _fingerprint(value: Any) -> str:
    return hashlib.sha256(_canonical_json(value).encode("utf-8")).hexdigest()


def _json_value(value: Any) -> Any:
    """Convert table cells to deterministic JSON without retaining references."""

    if value is None or value is pd.NA:
        return None
    if isinstance(value, np.generic):
        return _json_value(value.item())
    if isinstance(value, float):
        if math.isnan(value):
            return None
        if math.isinf(value):
            return {"__sophiark_float__": "Infinity" if value > 0 else "-Infinity"}
        return value
    if isinstance(value, (str, int, bool)):
        return value
    if isinstance(value, (datetime, pd.Timestamp)):
        timestamp = value.to_pydatetime() if isinstance(value, pd.Timestamp) else value
        if timestamp.tzinfo is None:
            timestamp = timestamp.replace(tzinfo=timezone.utc)
        return {"__sophiark_datetime__": timestamp.astimezone(timezone.utc).isoformat()}
    if isinstance(value, Path):
        return str(value)
    if is_dataclass(value):
        return _json_value(asdict(value))
    if isinstance(value, Mapping):
        return {
            str(key): _json_value(item)
            for key, item in sorted(value.items(), key=lambda pair: str(pair[0]))
        }
    if isinstance(value, (list, tuple)):
        return [_json_value(item) for item in value]
    if isinstance(value, (set, frozenset)):
        items = [_json_value(item) for item in value]
        return sorted(items, key=_canonical_json)
    try:
        if bool(pd.isna(value)):
            return None
    except (TypeError, ValueError):
        pass
    return str(value)


def _restore_json_value(value: Any) -> Any:
    if isinstance(value, list):
        return [_restore_json_value(item) for item in value]
    if isinstance(value, dict):
        if value.keys() == {"__sophiark_float__"}:
            return float(value["__sophiark_float__"])
        if value.keys() == {"__sophiark_datetime__"}:
            return value["__sophiark_datetime__"]
        return {key: _restore_json_value(item) for key, item in value.items()}
    return value


def freeze_result_table(frame: pd.DataFrame, *, name: str) -> FrozenResultTable:
    """Detach and serialize a DataFrame without mutating or retaining it."""

    if not isinstance(frame, pd.DataFrame):
        raise TypeError(f"{name} must be a pandas DataFrame")
    detached = frame.copy(deep=True)
    columns = tuple(map(str, detached.columns))
    if len(set(columns)) != len(columns):
        raise ValueError(f"{name} contains duplicate column names")
    detached.columns = list(columns)
    dtypes = tuple(str(detached[column].dtype) for column in columns)
    records = [
        {column: _json_value(row[position]) for position, column in enumerate(columns)}
        for row in detached.itertuples(index=False, name=None)
    ]
    rows_json = _canonical_json(records)
    table_fingerprint = _fingerprint({
        "columns": columns,
        "dtypes": dtypes,
        "records": records,
    })
    return FrozenResultTable(
        name=name,
        columns=columns,
        dtypes=dtypes,
        rows_json=rows_json,
        row_count=len(records),
        fingerprint=table_fingerprint,
    )


def thaw_result_table(table: FrozenResultTable) -> pd.DataFrame:
    """Materialize a fresh DataFrame, restoring tagged non-finite values."""

    records = [
        {key: _restore_json_value(value) for key, value in record.items()}
        for record in table.records()
    ]
    return pd.DataFrame(records, columns=list(table.columns))


def _utc_iso(value: datetime | str | None) -> str:
    if value is None:
        timestamp = datetime.now(timezone.utc)
    elif isinstance(value, datetime):
        timestamp = value
    else:
        text = str(value).strip()
        if text.endswith("Z"):
            text = f"{text[:-1]}+00:00"
        timestamp = datetime.fromisoformat(text)
    if timestamp.tzinfo is None:
        timestamp = timestamp.replace(tzinfo=timezone.utc)
    return timestamp.astimezone(timezone.utc).isoformat()


def _validate_species_taxon(species: str, taxon_id: int) -> None:
    normalized = " ".join(str(species).strip().casefold().split())
    if not normalized:
        raise ValueError("species must be explicit")
    expected = _KNOWN_TAXA.get(normalized)
    if expected is not None and expected != int(taxon_id):
        raise ValueError(
            f"species/taxon mismatch: {species!r} requires taxon_id={expected}, got {taxon_id}"
        )
    if int(taxon_id) <= 0:
        raise ValueError("taxon_id must be positive")


def _provenance_identity(record: ProvenanceRecord) -> dict[str, Any]:
    """Stable provenance identity; retrieval time is telemetry, not content."""

    return {
        "source": record.source,
        "kind": record.kind.value,
        "version": record.version,
        "locator": record.locator,
        "snapshot_id": record.snapshot_id,
        "attribution": record.attribution,
        "details": _json_value(record.details),
    }


def _default_provenance(
    report: FrozenResultTable,
    signed: FrozenResultTable | None,
    enrichment: FrozenResultTable | None,
    *,
    created_at: str,
) -> tuple[ProvenanceRecord, ...]:
    records = [
        ProvenanceRecord(
            source="Sophiark forward simulation report",
            kind=ProvenanceKind.SOPHIARK_COMPUTED,
            retrieved_at=created_at,
            snapshot_id=report.fingerprint,
            details={"table": report.name, "rows": report.row_count},
        )
    ]
    if signed is not None:
        records.append(ProvenanceRecord(
            source="Sophiark signed redistribution result",
            kind=ProvenanceKind.SOPHIARK_COMPUTED,
            retrieved_at=created_at,
            snapshot_id=signed.fingerprint,
            details={"table": signed.name, "rows": signed.row_count},
        ))
    if enrichment is not None:
        records.append(ProvenanceRecord(
            source="Sophiark post-simulation enrichment result",
            kind=ProvenanceKind.SOPHIARK_COMPUTED,
            retrieved_at=created_at,
            snapshot_id=enrichment.fingerprint,
            details={"table": enrichment.name, "rows": enrichment.row_count},
        ))
    return tuple(records)


def build_simulation_result_snapshot(
    report: pd.DataFrame,
    signed_redistribution: pd.DataFrame | None,
    enrichment: pd.DataFrame | None,
    *,
    species: str,
    taxon_id: int,
    tissue: str,
    targets: Sequence[str],
    attenuation: float | None,
    selection_mode: str,
    test_limit: int | None,
    tested_count: int | None,
    returned_count: int,
    top_n: int | None = None,
    threshold_parameters: Mapping[str, Any] | None = None,
    result_schema_version: int = RESULT_SCHEMA_VERSION,
    provenance: Iterable[ProvenanceRecord] | None = None,
    simulation_id: str | None = None,
    created_at: datetime | str | None = None,
) -> SimulationResultSnapshot:
    """Build a content-addressed read-only snapshot from explicit results.

    The function never reads Streamlit session state, output files, the graph,
    or external providers.  Input frames are defensively copied and serialized.
    Creation/retrieval timestamps are intentionally excluded from stable content
    fingerprints.
    """

    _validate_species_taxon(species, int(taxon_id))
    if not str(selection_mode).strip():
        raise ValueError("selection_mode must be explicit")
    if returned_count < 0 or (tested_count is not None and tested_count < 0):
        raise ValueError("tested_count and returned_count cannot be negative")
    if test_limit is not None and test_limit < 0:
        raise ValueError("test_limit cannot be negative")

    created_at_iso = _utc_iso(created_at)
    frozen_report = freeze_result_table(report, name="report")
    frozen_signed = (
        freeze_result_table(signed_redistribution, name="signed_redistribution")
        if isinstance(signed_redistribution, pd.DataFrame)
        else None
    )
    if signed_redistribution is not None and frozen_signed is None:
        raise TypeError("signed_redistribution must be a pandas DataFrame or None")
    frozen_enrichment = (
        freeze_result_table(enrichment, name="enrichment")
        if isinstance(enrichment, pd.DataFrame)
        else None
    )
    if enrichment is not None and frozen_enrichment is None:
        raise TypeError("enrichment must be a pandas DataFrame or None")

    target_ids = tuple(str(target).strip() for target in targets if str(target).strip())
    thresholds = {
        str(key): _json_value(value)
        for key, value in sorted((threshold_parameters or {}).items(), key=lambda pair: str(pair[0]))
    }
    scientific_payload = {
        "species": str(species).strip(),
        "taxon_id": int(taxon_id),
        "tissue": str(tissue).strip() or "None",
        "targets": target_ids,
        "attenuation": _json_value(attenuation),
        "selection_mode": str(selection_mode).strip(),
        "test_limit": test_limit,
        "tested_count": tested_count,
        "returned_count": returned_count,
        "top_n": top_n,
        "threshold_parameters": thresholds,
        "result_schema_version": int(result_schema_version),
        "tables": {
            "report": frozen_report.fingerprint,
            "signed_redistribution": frozen_signed.fingerprint if frozen_signed else None,
            "enrichment": frozen_enrichment.fingerprint if frozen_enrichment else None,
        },
    }
    scientific_fingerprint = _fingerprint(scientific_payload)
    resolved_simulation_id = (
        str(simulation_id).strip()
        if simulation_id is not None and str(simulation_id).strip()
        else f"sim-{scientific_fingerprint[:24]}"
    )

    provenance_records = tuple(provenance or _default_provenance(
        frozen_report,
        frozen_signed,
        frozen_enrichment,
        created_at=created_at_iso,
    ))
    snapshot_fingerprint = _fingerprint({
        "scientific_fingerprint": scientific_fingerprint,
        "provenance": [_provenance_identity(record) for record in provenance_records],
    })
    snapshot_id = f"snapshot-{_fingerprint({'simulation_id': resolved_simulation_id, 'snapshot': snapshot_fingerprint})[:24]}"

    return SimulationResultSnapshot(
        simulation_id=resolved_simulation_id,
        snapshot_id=snapshot_id,
        created_at=created_at_iso,
        species=str(species).strip(),
        taxon_id=int(taxon_id),
        tissue=str(tissue).strip() or "None",
        targets=target_ids,
        attenuation=float(attenuation) if attenuation is not None else None,
        selection_mode=str(selection_mode).strip(),
        test_limit=test_limit,
        threshold_parameters=thresholds,
        tested_count=tested_count,
        returned_count=returned_count,
        top_n=top_n,
        result_schema_version=int(result_schema_version),
        report=frozen_report,
        signed_redistribution=frozen_signed,
        enrichment=frozen_enrichment,
        provenance=provenance_records,
        scientific_fingerprint=scientific_fingerprint,
        snapshot_fingerprint=snapshot_fingerprint,
    )


# Concise facade alias for callers that do not need the longer type name.
build_snapshot = build_simulation_result_snapshot
