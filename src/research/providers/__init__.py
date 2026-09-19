"""Local and optional evidence-provider contracts.

The local provider exports below are definition-only imports: they perform no
dataset access, database connection, cache creation, or network request.
"""

from __future__ import annotations

from importlib import import_module
from typing import Any

from .local import (
    ComplexMembership,
    FamilyMembership,
    GeneAnnotation,
    HUMAN_TAXON_ID,
    LocalComplexProvider,
    LocalEvidenceAdapters,
    LocalHGNCFamilyProvider,
    LocalMyGeneProvider,
    LocalProviderResult,
    LocalProviderStatus,
    LocalRegulatoryProvider,
    LocalUniProtProvider,
    MOUSE_TAXON_ID,
    PartialUniProtContext,
    RegulatoryRelation,
)

_LAZY_EXTERNAL_EXPORTS = {
    "InterProProvider": ("src.research.providers.external_biology", "InterProProvider"),
    "QuickGOProvider": ("src.research.providers.external_biology", "QuickGOProvider"),
    "ReactomeProvider": ("src.research.providers.external_biology", "ReactomeProvider"),
    "UniProtProvider": ("src.research.providers.external_biology", "UniProtProvider"),
    "EuropePMCProvider": ("src.research.providers.literature", "EuropePMCProvider"),
}

__all__ = [
    "ComplexMembership",
    "EuropePMCProvider",
    "FamilyMembership",
    "GeneAnnotation",
    "HUMAN_TAXON_ID",
    "InterProProvider",
    "LocalComplexProvider",
    "LocalEvidenceAdapters",
    "LocalHGNCFamilyProvider",
    "LocalMyGeneProvider",
    "LocalProviderResult",
    "LocalProviderStatus",
    "LocalRegulatoryProvider",
    "LocalUniProtProvider",
    "MOUSE_TAXON_ID",
    "PartialUniProtContext",
    "QuickGOProvider",
    "ReactomeProvider",
    "RegulatoryRelation",
    "UniProtProvider",
]


def __getattr__(name: str) -> Any:
    target = _LAZY_EXTERNAL_EXPORTS.get(name)
    if target is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    value = getattr(import_module(target[0]), target[1])
    globals()[name] = value
    return value


def __dir__() -> list[str]:
    return sorted(set(globals()) | set(__all__))
