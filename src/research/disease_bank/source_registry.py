"""Central, auditable source metadata for Disease Bank providers."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from types import MappingProxyType
from typing import Mapping


@dataclass(frozen=True, slots=True)
class SourceRegistryEntry:
    source_key: str
    source_name: str
    source_version: str
    license: str
    commercial_use: bool
    attribution_required: bool
    redistribution_allowed: bool
    source_url: str
    attribution: str

    def metadata(self) -> dict[str, object]:
        """Return all stable licensing fields for ``EvidenceStore.sources``."""

        return asdict(self)


OPEN_TARGETS_SOURCE = SourceRegistryEntry(
    source_key="open_targets_platform",
    source_name="Open Targets Platform",
    source_version="GraphQL API v4",
    license="CC0 1.0 (Platform data; upstream source terms still apply)",
    commercial_use=True,
    attribution_required=True,
    redistribution_allowed=True,
    source_url="https://platform.opentargets.org/",
    attribution="Open Targets Platform; cite the Platform and respect original data-source terms.",
)

MONDO_SOURCE = SourceRegistryEntry(
    source_key="mondo",
    source_name="Mondo Disease Ontology",
    source_version="Cross-reference recorded by Open Targets",
    license="CC BY 4.0",
    commercial_use=True,
    attribution_required=True,
    redistribution_allowed=True,
    source_url="https://mondo.monarchinitiative.org/",
    attribution="Mondo Disease Ontology (CC BY 4.0).",
)

DISEASE_SOURCE_REGISTRY: Mapping[str, SourceRegistryEntry] = MappingProxyType({
    OPEN_TARGETS_SOURCE.source_key: OPEN_TARGETS_SOURCE,
    MONDO_SOURCE.source_key: MONDO_SOURCE,
})

__all__ = [
    "DISEASE_SOURCE_REGISTRY",
    "MONDO_SOURCE",
    "OPEN_TARGETS_SOURCE",
    "SourceRegistryEntry",
]
