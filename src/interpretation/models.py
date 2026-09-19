"""Read-only data models for deterministic biological interpretation."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Any


@dataclass(frozen=True)
class InterpretationCluster:
    """A transparent gene grouping; it never changes scientific ranking."""

    cluster_id: str
    cluster_type: str
    name: str
    members: tuple[str, ...]
    member_count: int
    confidence_score: float
    confidence_label: str
    evidence: str
    source: str
    heuristic_used: bool
    possible_topology_amplification: bool
    warning: str | None
    metrics: dict[str, float | int | None]
    topology_amplification_risk: str

    def to_dict(self) -> dict[str, Any]:
        """Expose a UI/export-safe representation without hidden state."""
        value = asdict(self)
        value["members"] = list(self.members)
        return value
