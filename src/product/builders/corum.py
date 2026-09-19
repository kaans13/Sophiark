"""CORUM shadow builder; membership remains context, never a new score."""

from __future__ import annotations

from collections import defaultdict
import pickle
from pathlib import Path

import pandas as pd

from src.product.fingerprints import sha256_file
from src.product.status import ProductStatus
from src.product.versions import PREPROCESSING_VERSIONS

from .common import SourceBuildError, SourceBuildReport, processed_identity, raw_identity, shadow_directory, write_report


DATASET_ID = "corum_human"
SCHEMA_VERSION = "corum-membership-1"
REQUIRED_COLUMNS = {"ComplexID", "ComplexName", "Organism", "subunits(Gene name)"}


def _unique_symbol_map(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str, usecols=["gene", "Symbol"]).dropna()
    counts = frame.groupby("Symbol")["gene"].nunique()
    unique = set(counts[counts.eq(1)].index)
    selected = frame[frame["Symbol"].isin(unique)]
    return dict(zip(selected["Symbol"], selected["gene"]))


def build_corum_shadow(
    source: str | Path,
    *, mapping_path: str | Path, output_root: str | Path | None = None,
    source_version: str = "UNKNOWN",
) -> SourceBuildReport:
    source_path = Path(source)
    raw_fingerprint = raw_identity(source_path)
    try:
        raw = pd.read_csv(source_path, sep="\t", dtype=str, keep_default_na=False)
    except Exception as exc:
        raise SourceBuildError(f"CORUM source could not be read: {exc}") from exc
    missing = sorted(REQUIRED_COLUMNS - set(raw.columns))
    if missing:
        raise SourceBuildError("incompatible CORUM schema; missing columns: " + ", ".join(missing))
    human = raw[raw["Organism"].str.strip().eq("Human")].copy()
    if human.empty:
        raise SourceBuildError("CORUM source contains no Human records")
    symbol_map = _unique_symbol_map(Path(mapping_path))
    rows: list[dict[str, object]] = []
    legacy: dict[str, list[str]] = defaultdict(list)
    seen: set[tuple[str, str]] = set()
    for record in human.to_dict("records"):
        complex_id = str(record["ComplexID"]).strip()
        complex_name = str(record["ComplexName"]).strip()
        if not complex_id or not complex_name:
            raise SourceBuildError("CORUM contains an empty complex ID or name")
        members = [
            item.strip() for item in str(record["subunits(Gene name)"]).split(";")
            if item.strip() and item.strip() != "None"
        ]
        for symbol in members:
            key = (complex_id, symbol)
            if key in seen:
                continue
            seen.add(key)
            legacy[symbol].append(complex_name)
            rows.append({
                "complex_id": complex_id, "complex_name": complex_name,
                "member_symbol": symbol, "entity_id": symbol_map.get(symbol),
                "organism": "Homo sapiens", "taxon_id": 9606,
            })
    if not rows:
        raise SourceBuildError("CORUM source contains no usable Human complex members")
    memberships = pd.DataFrame(rows).sort_values(
        ["complex_id", "member_symbol"], kind="mergesort"
    ).reset_index(drop=True)
    mapped_mask = memberships["entity_id"].notna()
    symbols = set(memberships["member_symbol"])
    mapped_symbols = symbols & set(symbol_map)
    directory = shadow_directory(
        dataset_id=DATASET_ID, raw_fingerprint=raw_fingerprint,
        processing_version=PREPROCESSING_VERSIONS["corum"], output_root=output_root,
        parameters={"taxon_id": 9606, "organism_label": "Human", "member_field": "subunits(Gene name)"},
    )
    csv_path = directory / "complex_memberships.csv"
    pickle_path = directory / "complex_members.pkl"
    memberships.to_csv(csv_path, index=False, encoding="utf-8", lineterminator="\n")
    with pickle_path.open("wb") as handle:
        pickle.dump(dict(legacy), handle, protocol=4)
    fingerprint = processed_identity(
        normalized_content_fingerprint=sha256_file(csv_path), source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["corum"], schema_version=SCHEMA_VERSION,
        parameters={"taxon_id": 9606, "member_field": "subunits(Gene name)"},
    )
    report = SourceBuildReport(
        dataset_id=DATASET_ID, source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["corum"], schema_version=SCHEMA_VERSION,
        source_records=len(memberships), valid_records=len(memberships),
        mapped_records=int(mapped_mask.sum()), unmapped_records=int((~mapped_mask).sum()),
        usable_records=len(memberships), unique_source_entities=len(symbols),
        mapped_entities=len(mapped_symbols), unmapped_entities=len(symbols - mapped_symbols),
        mapping_coverage=len(mapped_symbols) / len(symbols) if symbols else 0.0,
        raw_fingerprint=raw_fingerprint, processed_fingerprint=fingerprint,
        output_artifacts=(str(csv_path), str(pickle_path)), status=ProductStatus.READY,
        warnings=("Classic production complex_members.pkl has legacy/uncertain build provenance; shadow output is not promoted.",),
        details={
            "source_organisms": sorted(set(raw["Organism"].astype(str))),
            "selected_organism": "Human", "excluded_records": len(raw) - len(human),
            "source_file_records": len(raw),
            "human_complexes": int(human["ComplexID"].nunique()),
            "membership_semantics": "context-only membership; no pairwise score generated",
        },
    )
    write_report(report)
    return report


def audit_corum_shadow(*, production_pickle: str | Path, shadow_pickle: str | Path) -> dict[str, object]:
    with Path(production_pickle).open("rb") as handle:
        current = pickle.load(handle)
    with Path(shadow_pickle).open("rb") as handle:
        shadow = pickle.load(handle)
    current_keys, shadow_keys = set(current), set(shadow)
    common = current_keys & shadow_keys
    exact = sum(current[key] == shadow[key] for key in common)
    set_equal = sum(set(current[key]) == set(shadow[key]) for key in common)
    return {
        "production_entities": len(current), "shadow_entities": len(shadow),
        "common_entities": len(common), "production_only": len(current_keys - shadow_keys),
        "shadow_only": len(shadow_keys - current_keys), "exact_ordered_memberships": exact,
        "set_equal_memberships": set_equal,
        "authoritative_source_status": "AUTHORITATIVE_SOURCE_UNCERTAIN_FOR_CLASSIC_LEGACY_PICKLE",
    }
