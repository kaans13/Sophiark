"""TRRUST shadow builder reproducing the current legacy lookup semantics."""

from __future__ import annotations

from collections import defaultdict
import pickle
from pathlib import Path

import pandas as pd

from src.product.fingerprints import sha256_file
from src.product.status import ProductStatus
from src.product.versions import PREPROCESSING_VERSIONS

from .common import SourceBuildError, SourceBuildReport, processed_identity, raw_identity, shadow_directory, write_report


DATASET_ID = "trrust_human"
SCHEMA_VERSION = "trrust-lookup-1"
RELATION_LABELS = {"Activation": "Aktivasyon", "Repression": "Baskılama", "Unknown": "Bilinmiyor"}


def _unique_symbol_map(path: Path) -> dict[str, str]:
    frame = pd.read_csv(path, dtype=str, usecols=["gene", "Symbol"]).dropna()
    counts = frame.groupby("Symbol")["gene"].nunique()
    unique = set(counts[counts.eq(1)].index)
    return dict(zip(frame.loc[frame["Symbol"].isin(unique), "Symbol"], frame.loc[frame["Symbol"].isin(unique), "gene"]))


def build_trrust_shadow(
    source: str | Path,
    *, mapping_path: str | Path, output_root: str | Path | None = None,
    source_version: str = "UNKNOWN",
) -> SourceBuildReport:
    source_path = Path(source)
    raw_fingerprint = raw_identity(source_path)
    try:
        raw = pd.read_csv(source_path, sep="\t", header=None, names=["TF", "Target", "Relation", "PMID"], dtype=str)
    except Exception as exc:
        raise SourceBuildError(f"TRRUST source could not be read: {exc}") from exc
    if raw.shape[1] != 4 or raw.empty:
        raise SourceBuildError("incompatible TRRUST schema; expected four non-header TSV columns")
    raw = raw.fillna("")
    valid_mask = raw["TF"].str.strip().ne("") & raw["Target"].str.strip().ne("") & raw["Relation"].isin(RELATION_LABELS)
    valid = raw.loc[valid_mask].copy()
    if len(valid) != len(raw):
        raise SourceBuildError("TRRUST contains malformed or unsupported relation records")
    tf_targets: dict[str, set[str]] = defaultdict(set)
    target_regulators: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for tf, target, relation in valid[["TF", "Target", "Relation"]].itertuples(index=False, name=None):
        tf, target = str(tf).strip(), str(target).strip()
        tf_targets[tf].add(target)
        target_regulators[target].append((tf, RELATION_LABELS[relation]))
    tf_payload = {key: sorted(values) for key, values in tf_targets.items()}
    regulator_payload = dict(target_regulators)
    symbols = set(valid["TF"]) | set(valid["Target"])
    symbol_map = _unique_symbol_map(Path(mapping_path))
    mapped_symbols = symbols & set(symbol_map)
    both_mapped = valid["TF"].isin(symbol_map) & valid["Target"].isin(symbol_map)
    canonical = valid.loc[both_mapped, ["TF", "Target", "Relation", "PMID"]].copy()
    canonical.insert(0, "regulator_id", canonical["TF"].map(symbol_map))
    canonical.insert(1, "target_id", canonical["Target"].map(symbol_map))
    directory = shadow_directory(
        dataset_id=DATASET_ID, raw_fingerprint=raw_fingerprint,
        processing_version=PREPROCESSING_VERSIONS["trrust"], output_root=output_root,
        parameters={"taxon_id": 9606, "unknown_policy": "preserve"},
    )
    tf_path = directory / "tf_targets.pkl"
    regulators_path = directory / "target_to_regulators.pkl"
    canonical_path = directory / "canonical_relations.csv"
    with tf_path.open("wb") as handle:
        pickle.dump(tf_payload, handle, protocol=4)
    with regulators_path.open("wb") as handle:
        pickle.dump(regulator_payload, handle, protocol=4)
    canonical.sort_values(["regulator_id", "target_id", "Relation", "PMID"], kind="mergesort").to_csv(
        canonical_path, index=False, encoding="utf-8", lineterminator="\n"
    )
    content = sha256_file(canonical_path)
    fingerprint = processed_identity(
        normalized_content_fingerprint=content, source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["trrust"], schema_version=SCHEMA_VERSION,
        parameters={"taxon_id": 9606, "unknown_policy": "preserve"},
    )
    report = SourceBuildReport(
        dataset_id=DATASET_ID, source_version=source_version,
        processing_version=PREPROCESSING_VERSIONS["trrust"], schema_version=SCHEMA_VERSION,
        source_records=len(raw), valid_records=len(valid), mapped_records=int(both_mapped.sum()),
        unmapped_records=int((~both_mapped).sum()), usable_records=len(valid),
        unique_source_entities=len(symbols), mapped_entities=len(mapped_symbols),
        unmapped_entities=len(symbols - mapped_symbols),
        mapping_coverage=len(mapped_symbols) / len(symbols) if symbols else 0.0,
        raw_fingerprint=raw_fingerprint, processed_fingerprint=fingerprint,
        output_artifacts=(str(tf_path), str(regulators_path), str(canonical_path)),
        status=ProductStatus.READY,
        details={"organism": "Homo sapiens", "taxon_id": 9606, "relations_preserved": len(valid)},
    )
    write_report(report)
    return report


def audit_trrust_shadow(
    *, production_tf: str | Path, production_regulators: str | Path,
    shadow_tf: str | Path, shadow_regulators: str | Path,
) -> dict[str, object]:
    def load(path: str | Path):
        with Path(path).open("rb") as handle:
            return pickle.load(handle)
    current_tf, current_reg = load(production_tf), load(production_regulators)
    built_tf, built_reg = load(shadow_tf), load(shadow_regulators)
    return {
        "tf_targets_exact": current_tf == built_tf,
        "target_regulators_exact": current_reg == built_reg,
        "production_tf_count": len(current_tf), "shadow_tf_count": len(built_tf),
        "production_target_count": len(current_reg), "shadow_target_count": len(built_reg),
    }
