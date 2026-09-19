"""Local-first OmniPath ingestion and existing-resolver mapping."""

from __future__ import annotations

import hashlib
import json
from collections import defaultdict
from pathlib import Path
from typing import Any

import pandas as pd

from src.research.entity_resolution import EntityResolver
from src.research.models import EntityType, ResolutionStatus

from .models import (
    DatasetQA, LayerMembership, MappingReport, SignalingDataset, SignalingEdge,
    SignalingNode, SignStatus,
)


PARSER_VERSION = "omnipath-signaling-v2"
REQUIRED_COLUMNS = {
    "source", "target", "source_genesymbol", "target_genesymbol",
    "is_directed", "is_stimulation", "is_inhibition", "sources", "references",
}


def _boolean(value: Any) -> bool | None:
    if value is None or pd.isna(value):
        return None
    if isinstance(value, bool):
        return value
    return str(value).strip().casefold() in {"1", "true", "yes", "y"}


def _items(value: Any) -> tuple[str, ...]:
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ()
    return tuple(sorted({part.strip() for part in str(value).replace("|", ";").split(";") if part.strip()}))


def _metadata_path(snapshot: Path) -> Path:
    candidate = snapshot.with_name(f"{snapshot.stem}_metadata.json")
    return candidate if candidate.exists() else snapshot.with_name("omnipath_metadata.json")


def _read_metadata(snapshot: Path) -> dict[str, Any]:
    path = _metadata_path(snapshot)
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}


def _sign(stimulation: bool, inhibition: bool) -> tuple[SignStatus, int]:
    if stimulation and inhibition:
        return SignStatus.CONFLICTING, 0
    if stimulation:
        return SignStatus.ACTIVATION, 1
    if inhibition:
        return SignStatus.INHIBITION, -1
    return SignStatus.UNSIGNED, 0


def _mapping_acceptance(entity_rate: float, edge_rate: float) -> str:
    if edge_rate < .70:
        return "FAILURE"
    if entity_rate < .90 or edge_rate < .85:
        return "WARNING"
    return "PREFERRED"


def load_omnipath_signaling_dataset(
    path: str | Path,
    *,
    resolver: EntityResolver,
    taxon_id: int,
    minimum_rows: int = 100,
) -> SignalingDataset:
    """Load, normalize and map a local snapshot without network access."""

    snapshot = Path(path)
    if not snapshot.is_file():
        raise FileNotFoundError(f"Local OmniPath snapshot unavailable: {snapshot}")
    frame = pd.read_csv(snapshot, sep="\t", low_memory=False)
    missing = REQUIRED_COLUMNS - set(frame.columns)
    if missing:
        raise ValueError(f"OmniPath snapshot missing required columns: {sorted(missing)}")
    if len(frame) < minimum_rows:
        raise ValueError(f"OmniPath snapshot is abnormally small: {len(frame)} rows")

    metadata = _read_metadata(snapshot)
    checksum = hashlib.sha256(snapshot.read_bytes()).hexdigest()
    expected_checksum = metadata.get("checksum_sha256")
    if expected_checksum and expected_checksum != checksum:
        raise ValueError("OmniPath snapshot checksum does not match metadata")
    version = str(metadata.get("release") or checksum[:16])

    endpoints: dict[tuple[str, str], dict[str, Any]] = {}
    for identifier_column, symbol_column in (
        ("source", "source_genesymbol"), ("target", "target_genesymbol")
    ):
        for identifier, symbol in frame[[identifier_column, symbol_column]].itertuples(index=False, name=None):
            stable_id, label = str(identifier or "").strip(), str(symbol or "").strip()
            if stable_id and label:
                endpoints[(stable_id, label)] = {}

    resolved: dict[tuple[str, str], Any] = {}
    nodes: dict[str, SignalingNode] = {}
    ambiguous = unresolved = mapped = 0
    atomic_entities = mapped_atomic_entities = 0
    for endpoint in sorted(endpoints):
        stable_id, symbol = endpoint
        result = resolver.resolve(symbol, taxon_id=taxon_id, entity_type=EntityType.PROTEIN)
        resolved[endpoint] = result
        is_atomic = not stable_id.startswith("COMPLEX:")
        if is_atomic:
            atomic_entities += 1
        if result.status is ResolutionStatus.AMBIGUOUS:
            ambiguous += 1
        if result.entity is not None:
            mapped += 1
            if is_atomic:
                mapped_atomic_entities += 1
            node_id = result.entity.canonical_id
            membership = LayerMembership.LAYER_1_AND_2
            layer1_id = node_id
        else:
            if result.status is ResolutionStatus.UNRESOLVED:
                unresolved += 1
            node_id = f"omnipath:{stable_id}"
            membership = LayerMembership.SIGNALING_ONLY
            layer1_id = None
        existing = nodes.get(node_id)
        candidate = SignalingNode(
            node_id=node_id, source_identifier=stable_id, symbol=symbol,
            layer_membership=membership, layer1_id=layer1_id,
            mapping_status=result.status.value,
        )
        if existing is None or (existing.layer1_id is None and candidate.layer1_id is not None):
            nodes[node_id] = candidate

    aggregate: dict[tuple[str, str], dict[str, Any]] = defaultdict(lambda: {
        "stim": False, "inhib": False, "resources": set(), "references": set(),
        "raw_rows": 0, "consensus_direction": False,
        "consensus_stimulation": False, "consensus_inhibition": False,
        "curation_effort": 0, "source_ids": set(), "source_symbols": set(),
        "target_ids": set(), "target_symbols": set(),
        "raw_row_numbers": [],
    })
    unusable = 0
    fully_mapped = partially_mapped = signaling_only_edges = 0
    atomic_edges = mapped_atomic_edges = 0
    for raw_row_number, row in enumerate(frame.to_dict(orient="records"), start=2):
        if not _boolean(row.get("is_directed")):
            continue
        source_key = (str(row.get("source") or "").strip(), str(row.get("source_genesymbol") or "").strip())
        target_key = (str(row.get("target") or "").strip(), str(row.get("target_genesymbol") or "").strip())
        if source_key not in resolved or target_key not in resolved:
            unusable += 1
            continue
        source_result, target_result = resolved[source_key], resolved[target_key]
        is_atomic_edge = not source_key[0].startswith("COMPLEX:") and not target_key[0].startswith("COMPLEX:")
        if is_atomic_edge:
            atomic_edges += 1
        source_node = source_result.entity.canonical_id if source_result.entity else f"omnipath:{source_key[0]}"
        target_node = target_result.entity.canonical_id if target_result.entity else f"omnipath:{target_key[0]}"
        if source_result.entity and target_result.entity:
            fully_mapped += 1
            if is_atomic_edge:
                mapped_atomic_edges += 1
        elif source_result.entity or target_result.entity:
            partially_mapped += 1
        else:
            signaling_only_edges += 1
        item = aggregate[(source_node, target_node)]
        stimulation = _boolean(row.get("consensus_stimulation"))
        inhibition = _boolean(row.get("consensus_inhibition"))
        item["stim"] = item["stim"] or bool(stimulation if stimulation is not None else _boolean(row.get("is_stimulation")))
        item["inhib"] = item["inhib"] or bool(inhibition if inhibition is not None else _boolean(row.get("is_inhibition")))
        item["consensus_direction"] = item["consensus_direction"] or bool(_boolean(row.get("consensus_direction")))
        item["consensus_stimulation"] = item["consensus_stimulation"] or bool(stimulation)
        item["consensus_inhibition"] = item["consensus_inhibition"] or bool(inhibition)
        item["resources"].update(_items(row.get("sources")))
        item["references"].update(_items(row.get("references")))
        item["source_ids"].add(source_key[0])
        item["source_symbols"].add(source_key[1])
        item["target_ids"].add(target_key[0])
        item["target_symbols"].add(target_key[1])
        item["raw_rows"] += 1
        item["raw_row_numbers"].append(raw_row_number)
        try:
            item["curation_effort"] += int(row.get("curation_effort") or 0)
        except (TypeError, ValueError):
            pass

    edges: list[SignalingEdge] = []
    for (source_node, target_node), item in sorted(aggregate.items()):
        source, target = nodes[source_node], nodes[target_node]
        status, canonical = _sign(item["stim"], item["inhib"])
        edges.append(SignalingEdge(
            source=source_node, target=target_node,
            source_identifier=";".join(sorted(item["source_ids"])),
            source_symbol=";".join(sorted(item["source_symbols"])),
            target_identifier=";".join(sorted(item["target_ids"])),
            target_symbol=";".join(sorted(item["target_symbols"])),
            directed=True, stimulation=item["stim"], inhibition=item["inhib"],
            canonical_sign=canonical, sign_status=status, interaction="signaling/regulation",
            resources=tuple(sorted(item["resources"])), references=tuple(sorted(item["references"])),
            evidence_metadata={
                "raw_rows": item["raw_rows"], "raw_row_numbers": tuple(item["raw_row_numbers"]),
                "curation_effort": item["curation_effort"],
            },
            consensus_direction=item["consensus_direction"],
            consensus_stimulation=item["consensus_stimulation"],
            consensus_inhibition=item["consensus_inhibition"], dataset_version=version,
        ))

    raw_count = len(frame)
    entity_rate = mapped / len(endpoints) if endpoints else 0.0
    edge_rate = fully_mapped / raw_count if raw_count else 0.0
    atomic_entity_rate = mapped_atomic_entities / atomic_entities if atomic_entities else 0.0
    atomic_edge_rate = mapped_atomic_edges / atomic_edges if atomic_edges else 0.0
    mapping = MappingReport(
        unique_omnipath_entities=len(endpoints), mapped_layer1_entities=mapped,
        signaling_only_entities=len(endpoints) - mapped, ambiguous_entities=ambiguous,
        unresolved_entities=unresolved, total_source_edges=raw_count,
        fully_layer1_mapped_edges=fully_mapped, partially_layer1_mapped_edges=partially_mapped,
        signaling_only_edges=signaling_only_edges, unusable_edges=unusable,
        entity_mapping_rate=entity_rate, edge_mapping_rate=edge_rate,
        eligible_atomic_entities=atomic_entities, mapped_atomic_entities=mapped_atomic_entities,
        eligible_atomic_edges=atomic_edges, mapped_atomic_edges=mapped_atomic_edges,
        atomic_entity_mapping_rate=atomic_entity_rate, atomic_edge_mapping_rate=atomic_edge_rate,
        acceptance=_mapping_acceptance(atomic_entity_rate, atomic_edge_rate),
    )
    counts = {status: 0 for status in SignStatus}
    for edge in edges:
        counts[edge.sign_status] += 1
    qa = DatasetQA(
        raw_interactions=raw_count, normalized_interactions=len(edges), unique_entities=len(nodes),
        unique_directed_edges=len(edges), activating_edges=counts[SignStatus.ACTIVATION],
        inhibitory_edges=counts[SignStatus.INHIBITION], unsigned_edges=counts[SignStatus.UNSIGNED],
        conflicting_edges=counts[SignStatus.CONFLICTING], duplicates_removed=raw_count - len(edges),
        mapped_edges=fully_mapped, unmapped_entities=len(endpoints) - mapped,
    )
    dataset_metadata = {
        **metadata, "path": str(snapshot), "checksum_sha256": checksum,
        "parser_version": PARSER_VERSION, "dataset_version": version,
        "mapping_acceptance": mapping.acceptance,
    }
    fingerprint = hashlib.sha256(json.dumps({
        "checksum": checksum, "parser": PARSER_VERSION, "taxon_id": taxon_id,
        "edges": len(edges), "mapped": mapped,
    }, sort_keys=True).encode("utf-8")).hexdigest()
    return SignalingDataset(
        nodes=tuple(sorted(nodes.values(), key=lambda item: item.node_id)),
        edges=tuple(edges), qa=qa, mapping=mapping, metadata=dataset_metadata,
        fingerprint=fingerprint,
    )
