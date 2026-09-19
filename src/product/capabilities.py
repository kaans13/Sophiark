"""Capability detection derived exclusively from the central dataset registry."""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import Mapping

from .datasets import DataRegistry, DatasetHealth
from .status import ProductStatus


class Capability(str, Enum):
    HAS_PPI = "HAS_PPI"
    HAS_TISSUE_CONTEXT = "HAS_TISSUE_CONTEXT"
    HAS_DIRECTION = "HAS_DIRECTION"
    HAS_COMPLEX_CONTEXT = "HAS_COMPLEX_CONTEXT"
    HAS_EVIDENCE_CHANNELS = "HAS_EVIDENCE_CHANNELS"
    HAS_PHYSICAL_SUPPORT = "HAS_PHYSICAL_SUPPORT"
    HAS_CANONICAL_MAPPING = "HAS_CANONICAL_MAPPING"


ROLE_CAPABILITIES = {
    "PPI": Capability.HAS_PPI,
    "tissue_context": Capability.HAS_TISSUE_CONTEXT,
    "direction": Capability.HAS_DIRECTION,
    "complex_context": Capability.HAS_COMPLEX_CONTEXT,
    "evidence_channels": Capability.HAS_EVIDENCE_CHANNELS,
    "physical_support": Capability.HAS_PHYSICAL_SUPPORT,
    "canonical_mapping": Capability.HAS_CANONICAL_MAPPING,
}


@dataclass(frozen=True, slots=True)
class CapabilityReport:
    statuses: Mapping[Capability, ProductStatus]
    datasets: Mapping[str, DatasetHealth]

    def available(self, capability: Capability) -> bool:
        return self.statuses.get(capability) is ProductStatus.READY

    def as_dict(self) -> dict[str, object]:
        return {
            "capabilities": {key.value: value.value for key, value in self.statuses.items()},
            "datasets": {key: value.as_dict() for key, value in self.datasets.items()},
        }


def detect_capabilities(registry: DataRegistry) -> CapabilityReport:
    health = registry.health_report()
    by_capability: dict[Capability, list[ProductStatus]] = {
        capability: [] for capability in Capability
    }
    for name, item in health.items():
        for role in registry.get(name).roles:
            capability = ROLE_CAPABILITIES.get(role)
            if capability is not None:
                by_capability[capability].append(item.status)
    statuses = {
        capability: (
            ProductStatus.READY
            if any(status is ProductStatus.READY for status in values)
            else ProductStatus.DATA_UNAVAILABLE
        )
        for capability, values in by_capability.items()
    }
    return CapabilityReport(statuses, health)
