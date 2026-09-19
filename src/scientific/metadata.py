"""Analizlerin yeniden üretilebilirlik metadata sözleşmesi."""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
from typing import Any

from src.product.versions import SOPHIARK_VERSION


def file_snapshot(path: str | Path | None) -> str:
    if not path:
        return "not_configured"
    item = Path(path)
    if not item.exists():
        return "missing"
    stat = item.stat()
    payload = f"{item.resolve()}|{stat.st_size}|{stat.st_mtime_ns}".encode()
    return hashlib.sha256(payload).hexdigest()[:16]


def build_reproducibility_metadata(
    *,
    organism: str,
    tissue: str | None,
    tissue_normalization_mode: str,
    attenuation_fraction: float,
    attenuation_severity: str,
    pagerank_damping: float,
    bc_mode: str,
    bc_sample_size: int | None,
    null_iterations: int,
    random_seed: int,
    db_path: str | Path | None,
    efficiency_metadata: dict[str, Any] | None = None,
    string_version: str = "database snapshot; upstream version not encoded",
    expression_source_version: str = "database snapshot; upstream version not encoded",
    trrust_version: str = "local processed snapshot; release not encoded",
    omnipath_version: str = "not_loaded",
    edge_evidence_weighting_mode: str = "structural_only",
    cache_snapshot_identifiers: dict[str, str] | None = None,
) -> dict[str, Any]:
    return {
        "sophiark_version": SOPHIARK_VERSION,
        "analysis_timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "organism": organism,
        "tissue": tissue or "None",
        "tissue_normalization_mode": tissue_normalization_mode,
        "string_version": string_version,
        "expression_source_version": expression_source_version,
        "trrust_version": trrust_version,
        "omnipath_version": omnipath_version,
        "edge_evidence_weighting_mode": edge_evidence_weighting_mode,
        "enrichment_library_version": "reported by enrichment result",
        "attenuation_fraction": float(attenuation_fraction),
        "attenuation_severity": attenuation_severity,
        "pagerank_damping": float(pagerank_damping),
        "bc_mode": bc_mode,
        "bc_sample_size": bc_sample_size,
        "null_model_iterations": int(null_iterations),
        "random_seed": int(random_seed),
        "data_snapshot_identifier": file_snapshot(db_path),
        "cache_snapshot_identifiers": cache_snapshot_identifiers or {},
        "efficiency": efficiency_metadata or {},
        "deprecated_fields": {
            "Efficiency_Kayip_Pct": "Systemic_Network_Shift_Pct aliasıdır; global efficiency değildir",
            "Local_Efficiency_Kayip_Pct": "Top_Positive_PageRank_Mean_Pct aliasıdır; local efficiency değildir",
        },
    }


def write_metadata(metadata: dict[str, Any], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(metadata, ensure_ascii=False, indent=2), encoding="utf-8")


def attach_metadata_columns(frame, metadata: dict[str, Any]):
    """CSV/Excel exportunun sidecar olmadan da yeniden üretilebilir kalması."""
    mapping = {
        "Sophiark_Version": "sophiark_version",
        "Analysis_Timestamp_UTC": "analysis_timestamp_utc",
        "Organism": "organism",
        "Tissue": "tissue",
        "Tissue_Normalization_Mode": "tissue_normalization_mode",
        "STRING_Version": "string_version",
        "Expression_Source_Version": "expression_source_version",
        "TRRUST_Version": "trrust_version",
        "OmniPath_Version": "omnipath_version",
        "Enrichment_Library_Version": "enrichment_library_version",
        "PageRank_Damping": "pagerank_damping",
        "BC_Mode": "bc_mode",
        "BC_Sample_Size": "bc_sample_size",
        "Null_Model_Iterations": "null_model_iterations",
        "Random_Seed": "random_seed",
        "Data_Snapshot_ID": "data_snapshot_identifier",
        "Edge_Evidence_Weighting_Mode": "edge_evidence_weighting_mode",
        "Cache_Snapshot_IDs": "cache_snapshot_identifiers",
    }
    result = frame.copy()
    for column, key in mapping.items():
        value = metadata.get(key)
        if isinstance(value, (dict, list, tuple)):
            value = json.dumps(value, ensure_ascii=False, sort_keys=True)
        result[column] = value
    return result
