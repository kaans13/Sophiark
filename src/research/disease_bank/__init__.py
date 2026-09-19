"""Optional, cache-first Disease Bank for post-analysis reference comparison.

Nothing in this package imports a scientific engine.  It can only resolve
external disease evidence and compare it with already materialized results.
"""

from .benchmark import DiseaseBenchmark, benchmark_disease_reference
from .models import Disease, DiseaseGeneAssociation, DiseaseReferenceSet, MappingSummary
from .service import DiseaseBankService, DiseaseLoad
from .source_registry import MONDO_SOURCE, OPEN_TARGETS_SOURCE, SourceRegistryEntry

__all__ = [
    "Disease",
    "DiseaseBankService",
    "DiseaseBenchmark",
    "DiseaseGeneAssociation",
    "DiseaseLoad",
    "DiseaseReferenceSet",
    "MONDO_SOURCE",
    "MappingSummary",
    "OPEN_TARGETS_SOURCE",
    "SourceRegistryEntry",
    "benchmark_disease_reference",
]
