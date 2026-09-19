"""Immutable contracts for the optional OmniPath signaling layer."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping


class SignStatus(str, Enum):
    ACTIVATION = "activation"
    INHIBITION = "inhibition"
    UNSIGNED = "unsigned"
    CONFLICTING = "conflicting"


class LayerMembership(str, Enum):
    LAYER_1_AND_2 = "layer_1_and_2"
    SIGNALING_ONLY = "signaling_only"


class SignalingStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    PARTIAL = "PARTIAL"
    UNAVAILABLE = "UNAVAILABLE"


def frozen_mapping(value: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    return MappingProxyType(dict(value or {}))


@dataclass(frozen=True, slots=True)
class SignalingNode:
    node_id: str
    source_identifier: str
    symbol: str
    layer_membership: LayerMembership
    layer1_id: str | None = None
    mapping_status: str = "UNRESOLVED"


@dataclass(frozen=True, slots=True)
class SignalingEdge:
    source: str
    target: str
    source_identifier: str
    source_symbol: str
    target_identifier: str
    target_symbol: str
    directed: bool
    stimulation: bool
    inhibition: bool
    canonical_sign: int
    sign_status: SignStatus
    interaction: str
    resources: tuple[str, ...]
    references: tuple[str, ...]
    evidence_metadata: Mapping[str, Any] = field(default_factory=frozen_mapping)
    consensus_direction: bool | None = None
    consensus_stimulation: bool | None = None
    consensus_inhibition: bool | None = None
    dataset_version: str = "local_snapshot"

    def __post_init__(self) -> None:
        object.__setattr__(self, "sign_status", SignStatus(self.sign_status))
        object.__setattr__(self, "resources", tuple(self.resources))
        object.__setattr__(self, "references", tuple(self.references))
        object.__setattr__(self, "evidence_metadata", frozen_mapping(self.evidence_metadata))
        if self.canonical_sign not in {-1, 0, 1}:
            raise ValueError("canonical_sign must be -1, 0, or 1")
        if self.sign_status in {SignStatus.UNSIGNED, SignStatus.CONFLICTING} and self.canonical_sign != 0:
            raise ValueError("unsigned/conflicting edges must have canonical_sign=0")


@dataclass(frozen=True, slots=True)
class MappingReport:
    unique_omnipath_entities: int
    mapped_layer1_entities: int
    signaling_only_entities: int
    ambiguous_entities: int
    unresolved_entities: int
    total_source_edges: int
    fully_layer1_mapped_edges: int
    partially_layer1_mapped_edges: int
    signaling_only_edges: int
    unusable_edges: int
    entity_mapping_rate: float
    edge_mapping_rate: float
    eligible_atomic_entities: int
    mapped_atomic_entities: int
    eligible_atomic_edges: int
    mapped_atomic_edges: int
    atomic_entity_mapping_rate: float
    atomic_edge_mapping_rate: float
    acceptance: str


@dataclass(frozen=True, slots=True)
class DatasetQA:
    raw_interactions: int
    normalized_interactions: int
    unique_entities: int
    unique_directed_edges: int
    activating_edges: int
    inhibitory_edges: int
    unsigned_edges: int
    conflicting_edges: int
    duplicates_removed: int
    mapped_edges: int
    unmapped_entities: int


@dataclass(frozen=True, slots=True)
class SignalingDataset:
    nodes: tuple[SignalingNode, ...]
    edges: tuple[SignalingEdge, ...]
    qa: DatasetQA
    mapping: MappingReport
    metadata: Mapping[str, Any]
    fingerprint: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "nodes", tuple(self.nodes))
        object.__setattr__(self, "edges", tuple(self.edges))
        object.__setattr__(self, "metadata", frozen_mapping(self.metadata))


@dataclass(frozen=True, slots=True)
class PathSummary:
    shortest_distance: int | None
    example_path: tuple[str, ...]
    path_count: int
    positive_paths: int
    negative_paths: int
    unsigned_paths: int
    uncertain_paths: int
    net_sign: str


@dataclass(frozen=True, slots=True)
class SignalingContextRow:
    perturbation_target: str
    response_entity: str
    response_type: str
    network_response_value: float | None
    directed_relation: str
    directed_distance: int | None
    sign: str
    example_path: tuple[str, ...]
    path_count_within_depth: int
    positive_path_count: int
    negative_path_count: int
    unsigned_path_count: int
    uncertain_path_count: int
    alternative_route_candidate: bool
    tissue_support: str
    evidence_resources: tuple[str, ...]
    references: tuple[str, ...]
    dataset_version: str


@dataclass(frozen=True, slots=True)
class ConvergenceRecord:
    convergence_node: str
    converging_targets: tuple[str, ...]
    distances: Mapping[str, int]
    signed_relations: Mapping[str, str]
    resources: tuple[str, ...]
    references: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "distances", frozen_mapping(self.distances))
        object.__setattr__(self, "signed_relations", frozen_mapping(self.signed_relations))


@dataclass(frozen=True, slots=True)
class SignalingAnalysisResult:
    status: SignalingStatus
    rows: tuple[SignalingContextRow, ...] = ()
    convergence: tuple[ConvergenceRecord, ...] = ()
    metadata: Mapping[str, Any] = field(default_factory=frozen_mapping)
    notices: tuple[str, ...] = ()
    error: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "status", SignalingStatus(self.status))
        object.__setattr__(self, "rows", tuple(self.rows))
        object.__setattr__(self, "convergence", tuple(self.convergence))
        object.__setattr__(self, "metadata", frozen_mapping(self.metadata))
        object.__setattr__(self, "notices", tuple(self.notices))
