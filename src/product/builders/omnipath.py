"""Thin shadow-build adapter around the validated OmniPath ingestion path."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from src.product.status import ProductStatus
from src.product.versions import PREPROCESSING_VERSIONS
from src.research.entity_resolution import EntityResolver
from src.signaling.ingestion import PARSER_VERSION, load_omnipath_signaling_dataset
from src.signaling.storage import load_signaling_cache, write_signaling_cache

from .common import SourceBuildReport, raw_identity, shadow_directory, write_report


DATASET_ID = "omnipath_human"
SCHEMA_VERSION = "omnipath-signaling-cache-2"


def _metadata(source: Path) -> dict[str, Any]:
    path = source.with_name(f"{source.stem}_metadata.json")
    if not path.is_file():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def build_omnipath_shadow(
    source: str | Path,
    *, project_root: str | Path, output_root: str | Path | None = None,
    taxon_id: int = 9606, minimum_rows: int = 100,
) -> SourceBuildReport:
    source_path = Path(source)
    raw_fingerprint = raw_identity(source_path)
    metadata = _metadata(source_path)
    declared_taxon = metadata.get("organism_taxid")
    if declared_taxon is not None and int(declared_taxon) != int(taxon_id):
        raise ValueError(f"OmniPath metadata taxon {declared_taxon} does not match selected {taxon_id}")
    resolver = EntityResolver.from_project_data(
        project_root, include_human=taxon_id == 9606, include_mouse=taxon_id == 10090,
    )
    dataset = load_omnipath_signaling_dataset(
        source_path, resolver=resolver, taxon_id=taxon_id, minimum_rows=minimum_rows,
    )
    directory = shadow_directory(
        dataset_id=DATASET_ID if taxon_id == 9606 else f"omnipath_{taxon_id}",
        raw_fingerprint=raw_fingerprint,
        processing_version=PREPROCESSING_VERSIONS["omnipath"], output_root=output_root,
        parameters={"taxon_id": taxon_id, "parser": PARSER_VERSION},
    )
    cache_path = directory / f"omnipath_signaling_{taxon_id}.sqlite"
    write_signaling_cache(dataset, cache_path, source_snapshot=source_path)
    mapping = dataset.mapping
    report = SourceBuildReport(
        dataset_id=DATASET_ID if taxon_id == 9606 else f"omnipath_{taxon_id}",
        source_version=str(dataset.metadata.get("dataset_version", "UNKNOWN")),
        processing_version=PREPROCESSING_VERSIONS["omnipath"], schema_version=SCHEMA_VERSION,
        source_records=dataset.qa.raw_interactions, valid_records=dataset.qa.raw_interactions,
        mapped_records=mapping.fully_layer1_mapped_edges,
        unmapped_records=dataset.qa.raw_interactions - mapping.fully_layer1_mapped_edges,
        usable_records=dataset.qa.normalized_interactions,
        unique_source_entities=mapping.unique_omnipath_entities,
        mapped_entities=mapping.mapped_layer1_entities,
        unmapped_entities=mapping.unique_omnipath_entities - mapping.mapped_layer1_entities,
        mapping_coverage=mapping.entity_mapping_rate,
        raw_fingerprint=raw_fingerprint, processed_fingerprint=dataset.fingerprint,
        output_artifacts=(str(cache_path),), status=ProductStatus.READY,
        warnings=(() if mapping.acceptance == "PREFERRED" else (f"mapping acceptance: {mapping.acceptance}",)),
        details={
            "taxon_id": taxon_id, "source_organism": metadata.get("organism", "UNKNOWN"),
            "selected_organism": taxon_id, "directed_pairs": dataset.qa.unique_directed_edges,
            "activating_edges": dataset.qa.activating_edges,
            "inhibitory_edges": dataset.qa.inhibitory_edges,
            "unsigned_edges": dataset.qa.unsigned_edges,
            "conflicting_edges": dataset.qa.conflicting_edges,
            "direction_policy": "preserved by existing signaling ingestion; no numerical weight created",
        },
    )
    write_report(report)
    return report


def audit_omnipath_shadow(*, production_cache: str | Path, shadow_cache: str | Path, source: str | Path) -> dict[str, object]:
    current = load_signaling_cache(production_cache, source_snapshot=source)
    shadow = load_signaling_cache(shadow_cache, source_snapshot=source)
    return {
        "fingerprint_equal": current.fingerprint == shadow.fingerprint,
        "nodes_equal": current.nodes == shadow.nodes,
        "edges_equal": current.edges == shadow.edges,
        "qa_equal": current.qa == shadow.qa,
        "mapping_equal": current.mapping == shadow.mapping,
        "production_nodes": len(current.nodes), "shadow_nodes": len(shadow.nodes),
        "production_edges": len(current.edges), "shadow_edges": len(shadow.edges),
    }
