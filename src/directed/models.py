"""Contracts for the optional hybrid directed calculation engine."""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from types import MappingProxyType
from typing import Any, Mapping

import igraph as ig
import pandas as pd


class DirectedEngineStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    PARTIAL = "PARTIAL"
    UNAVAILABLE = "UNAVAILABLE"


class PerturbationStrategy(str, Enum):
    INCIDENT_EDGE_ATTENUATION = "incident_edge_attenuation"


def _frozen(value: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    return MappingProxyType(dict(value or {}))


@dataclass(frozen=True, slots=True)
class HybridGraphProvenance:
    string_edges: int
    direction_resolved_edges: int
    single_direction_edges: int
    bidirectional_evidence_edges: int
    fallback_string_edges: int
    directed_arcs: int
    omnipath_layer1_edges_not_in_string: int
    omnipath_mapping_misses: int
    affected_unique_nodes: int
    dropped_string_edges: int
    direction_coverage: float
    tissue_context: str
    string_dataset_version: str
    omnipath_fingerprint: str
    direction_policy: str
    fallback_policy: str
    cache_key: str


@dataclass(slots=True)
class HybridDirectedGraph:
    graph: ig.Graph
    provenance: HybridGraphProvenance


@dataclass(slots=True)
class DirectedCalculationResult:
    status: DirectedEngineStatus
    report: pd.DataFrame = field(default_factory=pd.DataFrame)
    directed_graph: HybridDirectedGraph | None = None
    metadata: Mapping[str, Any] = field(default_factory=_frozen)
    timings: Mapping[str, float] = field(default_factory=_frozen)
    error: str | None = None

    def __post_init__(self) -> None:
        self.status = DirectedEngineStatus(self.status)
        self.metadata = _frozen(self.metadata)
        self.timings = _frozen(self.timings)


@dataclass(slots=True)
class EngineComparisonResult:
    table: pd.DataFrame
    validation: Mapping[str, Any]

    def __post_init__(self) -> None:
        self.validation = _frozen(self.validation)
