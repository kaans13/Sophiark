"""Structural PPI'dan ayrÄ±, versioned directed-regulatory evidence katmanÄ±."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import pandas as pd


PARSER_VERSION = "omnipath-tsv-v1"
OVERLAY_COLUMNS = [
    "source", "target", "direction", "stimulation", "inhibition",
    "consensus_direction", "consensus_stimulation", "consensus_inhibition",
    "provenance", "references", "organism", "snapshot_version",
    "evidence_status", "evidence_layer",
]


@dataclass(frozen=True)
class DirectedEvidence:
    source: str
    target: str
    direction: str | None
    stimulation: bool | None
    inhibition: bool | None
    consensus_direction: bool | None
    consensus_stimulation: bool | None
    consensus_inhibition: bool | None
    provenance: tuple[str, ...]
    evidence_layer: str = "directed_regulatory"


def _bool(value):
    if pd.isna(value):
        return None
    if isinstance(value, bool):
        return value
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def _text_list(value) -> str:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    if isinstance(value, (tuple, list, set)):
        return ";".join(sorted(map(str, value)))
    return str(value)


def _status(stimulation, inhibition, directed) -> str:
    if stimulation and inhibition:
        return "conflicting_evidence"
    if stimulation:
        return "consensus_stimulation"
    if inhibition:
        return "consensus_inhibition"
    if directed:
        return "direction_only"
    return "unsupported_direction"


def trrust_records(target_to_regulators: dict) -> pd.DataFrame:
    records: list[dict] = []
    for target, regulators in (target_to_regulators or {}).items():
        for item in regulators:
            source = str(item[0])
            relation = str(item[1]) if len(item) > 1 else "Unknown"
            lowered = relation.casefold()
            stimulation = True if ("activ" in lowered or "uyar" in lowered) else None
            inhibition = True if ("repress" in lowered or "inhib" in lowered or "bask" in lowered) else None
            row = asdict(DirectedEvidence(
                source=source, target=str(target), direction="source_to_target",
                stimulation=stimulation, inhibition=inhibition, consensus_direction=True,
                consensus_stimulation=stimulation, consensus_inhibition=inhibition,
                provenance=("TRRUST",),
            ))
            row.update({
                "references": "", "organism": "Homo sapiens", "snapshot_version": "local_processed",
                "evidence_status": _status(stimulation, inhibition, True),
            })
            records.append(row)
    return pd.DataFrame(records, columns=OVERLAY_COLUMNS)


def _read_metadata(snapshot: Path) -> dict:
    metadata_path = snapshot.with_name(f"{snapshot.stem}_metadata.json")
    if not metadata_path.exists():
        # Backward compatibility with the first single-snapshot layout.
        metadata_path = snapshot.with_name("omnipath_metadata.json")
    if not metadata_path.exists():
        return {"available": True, "metadata_available": False}
    try:
        return json.loads(metadata_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        return {"available": True, "metadata_available": False, "metadata_error": str(exc)}


def omnipath_snapshot_metadata(path: str | Path | None) -> dict[str, object]:
    """BÃ¼yÃ¼k TSV'yi okumadan analysis metadata iÃ§in snapshot durumunu verir."""
    if not path or not Path(path).exists():
        return {"available": False, "version": "unavailable"}
    snapshot = Path(path)
    metadata = _read_metadata(snapshot)
    return {
        "available": True,
        "version": metadata.get("release") or metadata.get("snapshot_id") or "local_snapshot",
        "checksum_sha256": metadata.get("checksum_sha256") or hashlib.sha256(snapshot.read_bytes()).hexdigest(),
        "path": str(snapshot),
    }


def _normalise_snapshot(frame: pd.DataFrame, metadata: dict) -> pd.DataFrame:
    # OmniPath'in standart TSV'si hem UniProt (`source`) hem gene symbol
    # (`source_genesymbol`) iÃ§erir. Sophiark TRRUST ile aynÄ± namespace'te
    # aÃ§Ä±klama yapabilsin diye gene symbol'u tercih eder.
    if {"source_genesymbol", "target_genesymbol"}.issubset(frame.columns):
        frame = frame.drop(columns=[column for column in ("source", "target") if column in frame.columns]).rename(
            columns={"source_genesymbol": "source", "target_genesymbol": "target"}
        )
    aliases = {
        "sources": "provenance", "is_directed": "consensus_direction",
        "is_stimulation": "stimulation", "is_inhibition": "inhibition",
        "consensus_stimulation": "consensus_stimulation",
        "consensus_inhibition": "consensus_inhibition",
    }
    frame = frame.rename(columns={
        old: new for old, new in aliases.items()
        if old in frame.columns and new not in frame.columns
    })
    for field in ("stimulation", "inhibition", "consensus_direction", "consensus_stimulation", "consensus_inhibition"):
        frame[field] = frame[field].map(_bool) if field in frame else None
    frame["direction"] = frame["consensus_direction"].map(lambda x: "source_to_target" if x else None)
    frame["provenance"] = frame.get("provenance", "").map(_text_list) if "provenance" in frame else ""
    frame["references"] = frame.get("references", "").map(_text_list) if "references" in frame else ""
    frame["organism"] = frame.get("organism", metadata.get("organism", "Homo sapiens"))
    frame["snapshot_version"] = metadata.get("release") or metadata.get("snapshot_id") or "local_snapshot"
    frame["evidence_status"] = [
        _status(s, i, d) for s, i, d in zip(frame["consensus_stimulation"], frame["consensus_inhibition"], frame["consensus_direction"])
    ]
    frame["evidence_layer"] = "directed_regulatory"
    return frame.reindex(columns=OVERLAY_COLUMNS)


def load_omnipath_overlay(path: str | Path | None) -> tuple[pd.DataFrame, dict[str, object]]:
    """Local snapshot yÃ¼kler; bulunmaması negatif evidence anlamÄ±na gelmez."""
    if not path or not Path(path).exists():
        return pd.DataFrame(columns=OVERLAY_COLUMNS), {
            "available": False,
            "status": "OmniPath unavailable",
            "parser_version": PARSER_VERSION,
            # ASCII escape keeps this compatibility text stable on legacy Windows encodings.
            "reason": "Local OmniPath snapshot unavailable; negatif evidence de" + chr(0x011F) + "ildir.",
        }
    snapshot = Path(path)
    metadata = _read_metadata(snapshot)
    frame = _normalise_snapshot(pd.read_csv(snapshot, sep="\t"), metadata)
    checksum = hashlib.sha256(snapshot.read_bytes()).hexdigest()
    conflicts = int((frame["evidence_status"] == "conflicting_evidence").sum())
    return frame, {
        **metadata,
        "available": True,
        "status": "loaded",
        "records": len(frame),
        "conflicting_records": conflicts,
        "path": str(snapshot),
        "checksum_sha256": checksum,
        "parser_version": PARSER_VERSION,
    }
