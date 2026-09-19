"""Immutable Disease Bank records, independent from scientific result models."""

from __future__ import annotations

from dataclasses import dataclass, field
from types import MappingProxyType
from typing import Mapping

from ..models import EntityResolution


def _text(value: object, field_name: str, *, required: bool = False) -> str:
    text = str(value or "").strip()
    if required and not text:
        raise ValueError(f"{field_name} is required")
    return text


@dataclass(frozen=True, slots=True)
class Disease:
    disease_id: str
    disease_name: str
    ontology_source: str = "Open Targets"
    synonyms: tuple[str, ...] = ()
    parents: tuple[str, ...] = ()
    mondo_id: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "disease_id", _text(self.disease_id, "disease_id", required=True))
        object.__setattr__(self, "disease_name", _text(self.disease_name, "disease_name", required=True))
        object.__setattr__(self, "ontology_source", _text(self.ontology_source, "ontology_source", required=True))
        object.__setattr__(self, "synonyms", tuple(dict.fromkeys(_text(item, "synonym") for item in self.synonyms if _text(item, "synonym"))))
        object.__setattr__(self, "parents", tuple(dict.fromkeys(_text(item, "parent") for item in self.parents if _text(item, "parent"))))
        mondo = _text(self.mondo_id, "mondo_id")
        object.__setattr__(self, "mondo_id", mondo or None)


@dataclass(frozen=True, slots=True)
class DiseaseGeneAssociation:
    disease_id: str
    gene_symbol: str = ""
    ensembl_gene_id: str = ""
    target_name: str = ""
    association_score: float | None = None
    evidence_count: int | None = None
    evidence: Mapping[str, object] = field(default_factory=dict)
    source_key: str = "open_targets_platform"
    source_version: str = "GraphQL API v4"
    retrieved_at: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "disease_id", _text(self.disease_id, "disease_id", required=True))
        object.__setattr__(self, "gene_symbol", _text(self.gene_symbol, "gene_symbol"))
        object.__setattr__(self, "ensembl_gene_id", _text(self.ensembl_gene_id, "ensembl_gene_id"))
        object.__setattr__(self, "target_name", _text(self.target_name, "target_name"))
        if not (self.gene_symbol or self.ensembl_gene_id):
            raise ValueError("an association requires an Ensembl gene ID or a gene symbol")
        if self.association_score is not None:
            object.__setattr__(self, "association_score", float(self.association_score))
        if self.evidence_count is not None:
            if int(self.evidence_count) < 0:
                raise ValueError("evidence_count cannot be negative")
            object.__setattr__(self, "evidence_count", int(self.evidence_count))
        object.__setattr__(self, "evidence", MappingProxyType(dict(self.evidence)))
        object.__setattr__(self, "source_key", _text(self.source_key, "source_key", required=True))
        object.__setattr__(self, "source_version", _text(self.source_version, "source_version", required=True))
        object.__setattr__(self, "retrieved_at", _text(self.retrieved_at, "retrieved_at", required=True))


@dataclass(frozen=True, slots=True)
class ResolvedDiseaseGene:
    association: DiseaseGeneAssociation
    resolution: EntityResolution
    query_used: str


@dataclass(frozen=True, slots=True)
class MappingSummary:
    total: int
    mapped: int
    unmapped: int
    ambiguous: int

    def __post_init__(self) -> None:
        if min(self.total, self.mapped, self.unmapped, self.ambiguous) < 0:
            raise ValueError("mapping counts cannot be negative")
        if self.mapped + self.unmapped + self.ambiguous != self.total:
            raise ValueError("mapping counts must partition total disease genes")


@dataclass(frozen=True, slots=True)
class DiseaseReferenceSet:
    disease: Disease
    associations: tuple[DiseaseGeneAssociation, ...]
    resolved: tuple[ResolvedDiseaseGene, ...]
    mapping_summary: MappingSummary
    source_status: str
    source_message: str | None = None

    @property
    def canonical_proteins(self) -> frozenset[str]:
        return frozenset(
            item.resolution.entity.canonical_id
            for item in self.resolved
            if item.resolution.entity is not None
        )


__all__ = [
    "Disease",
    "DiseaseGeneAssociation",
    "DiseaseReferenceSet",
    "MappingSummary",
    "ResolvedDiseaseGene",
]
