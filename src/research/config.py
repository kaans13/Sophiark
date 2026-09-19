"""Configuration for the optional, read-only Research Context subsystem.

Importing this module performs no I/O, creates no cache, and starts no network
request.  External providers are deliberately disabled by default so the
scientific result path remains independent from research enrichment.
"""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass(frozen=True, slots=True)
class ProviderTimeout:
    """Finite network limits for one optional external provider request."""

    connect_seconds: float = 2.0
    read_seconds: float = 5.0
    total_seconds: float = 8.0

    def __post_init__(self) -> None:
        if min(self.connect_seconds, self.read_seconds, self.total_seconds) <= 0:
            raise ValueError("Provider timeouts must be positive.")
        if self.total_seconds < max(self.connect_seconds, self.read_seconds):
            raise ValueError("total_seconds cannot be shorter than a component timeout.")


@dataclass(frozen=True, slots=True)
class ProviderTimeouts:
    """Provider-specific limits; providers may fail without affecting core."""

    uniprot: ProviderTimeout = field(default_factory=ProviderTimeout)
    interpro: ProviderTimeout = field(default_factory=ProviderTimeout)
    quickgo: ProviderTimeout = field(default_factory=ProviderTimeout)
    reactome: ProviderTimeout = field(default_factory=ProviderTimeout)
    europepmc: ProviderTimeout = field(
        default_factory=lambda: ProviderTimeout(connect_seconds=2.0, read_seconds=8.0, total_seconds=10.0)
    )


@dataclass(frozen=True, slots=True)
class CacheTTLs:
    """Persistent-cache lifetimes, in seconds, by evidence class."""

    entity_resolution_seconds: int = 30 * 24 * 60 * 60
    biological_annotation_seconds: int = 30 * 24 * 60 * 60
    pathway_annotation_seconds: int = 14 * 24 * 60 * 60
    literature_seconds: int = 7 * 24 * 60 * 60
    negative_result_seconds: int = 6 * 60 * 60
    transient_failure_seconds: int = 60

    def __post_init__(self) -> None:
        if min(
            self.entity_resolution_seconds,
            self.biological_annotation_seconds,
            self.pathway_annotation_seconds,
            self.literature_seconds,
            self.negative_result_seconds,
            self.transient_failure_seconds,
        ) <= 0:
            raise ValueError("Cache TTL values must be positive.")


@dataclass(frozen=True, slots=True)
class RetrievalBudget:
    """Hard request/entity caps for an explicitly requested context load."""

    max_entities: int
    max_families: int
    max_functional_terms: int
    max_relationships: int
    max_provider_requests: int
    max_literature_queries: int
    max_publications_per_query: int

    def __post_init__(self) -> None:
        if min(
            self.max_entities,
            self.max_families,
            self.max_functional_terms,
            self.max_relationships,
            self.max_provider_requests,
            self.max_literature_queries,
            self.max_publications_per_query,
        ) < 0:
            raise ValueError("Retrieval budget values cannot be negative.")


@dataclass(frozen=True, slots=True)
class PresentationLimits:
    """UI-only limits; none of these values is a significance threshold."""

    family_min_members: int = 2
    max_returned_candidates: int = 200
    # First-paint disclosure cap, not a scientific or significance threshold.
    max_overview_observations: int = 6
    max_family_rows: int = 50
    max_function_rows: int = 100
    max_relationship_rows: int = 200
    max_literature_rows: int = 100

    def __post_init__(self) -> None:
        if self.family_min_members < 2:
            raise ValueError("family_min_members must suppress singleton presentation rows.")
        if min(
            self.max_returned_candidates,
            self.max_overview_observations,
            self.max_family_rows,
            self.max_function_rows,
            self.max_relationship_rows,
            self.max_literature_rows,
        ) <= 0:
            raise ValueError("Presentation limits must be positive.")

    @property
    def family_threshold_presentation_only(self) -> bool:
        """Make the non-scientific nature of the family threshold explicit."""

        return True


def _quick_budget() -> RetrievalBudget:
    return RetrievalBudget(
        max_entities=20,
        max_families=10,
        max_functional_terms=20,
        max_relationships=40,
        max_provider_requests=6,
        max_literature_queries=3,
        max_publications_per_query=10,
    )


def _deep_budget() -> RetrievalBudget:
    return RetrievalBudget(
        max_entities=100,
        max_families=50,
        max_functional_terms=100,
        max_relationships=250,
        max_provider_requests=30,
        max_literature_queries=20,
        max_publications_per_query=50,
    )


@dataclass(frozen=True, slots=True)
class ResearchConfig:
    """Single configuration boundary for Research Context behavior."""

    research_context_enabled: bool = False
    external_context_enabled: bool = False
    directed_signaling_enabled: bool = False
    directed_signaling_max_depth: int = 4

    provider_uniprot_enabled: bool = False
    provider_interpro_enabled: bool = False
    provider_quickgo_enabled: bool = False
    provider_reactome_enabled: bool = False
    provider_europepmc_enabled: bool = False

    provider_timeouts: ProviderTimeouts = field(default_factory=ProviderTimeouts)
    cache_ttls: CacheTTLs = field(default_factory=CacheTTLs)
    quick_context_budget: RetrievalBudget = field(default_factory=_quick_budget)
    deep_context_budget: RetrievalBudget = field(default_factory=_deep_budget)
    presentation_limits: PresentationLimits = field(default_factory=PresentationLimits)

    def __post_init__(self) -> None:
        if not 1 <= int(self.directed_signaling_max_depth) <= 6:
            raise ValueError("directed_signaling_max_depth must be between 1 and 6")

    def provider_enabled(self, provider_name: str) -> bool:
        """Return an explicit flag; unknown providers are always disabled."""

        key = str(provider_name).strip().casefold().replace("-", "").replace("_", "")
        flags = {
            "uniprot": self.provider_uniprot_enabled,
            "interpro": self.provider_interpro_enabled,
            "quickgo": self.provider_quickgo_enabled,
            "reactome": self.provider_reactome_enabled,
            "europepmc": self.provider_europepmc_enabled,
        }
        return bool(self.external_context_enabled and flags.get(key, False))


DEFAULT_RESEARCH_CONFIG = ResearchConfig()


__all__ = [
    "CacheTTLs",
    "DEFAULT_RESEARCH_CONFIG",
    "PresentationLimits",
    "ProviderTimeout",
    "ProviderTimeouts",
    "ResearchConfig",
    "RetrievalBudget",
]
