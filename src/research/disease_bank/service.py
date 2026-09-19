"""Application service for explicit Disease Bank retrieval and safe mapping."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

from ..config import ProviderTimeout
from ..entity_resolution import EntityResolver
from ..evidence_store import EvidenceStore, EvidenceStoreError
from ..models import ResolutionStatus
from ..provider_cache import CacheEvidenceClass, CacheFirstProviderClient, ProviderCache
from ..provider_runtime import ProviderExecutionPolicy, ProviderExecutor, ProviderRequest, ProviderStatus
from .models import Disease, DiseaseReferenceSet, MappingSummary, ResolvedDiseaseGene
from .opentargets import OpenTargetsDiseaseProvider, associations_from_payload, disease_from_payload
from .repository import DiseaseRepository


@dataclass(frozen=True, slots=True)
class DiseaseLoad:
    status: str
    disease: Disease | None = None
    reference: DiseaseReferenceSet | None = None
    message: str | None = None


class DiseaseBankService:
    """Optional service: network execution never constructs or calls this class."""

    def __init__(self, *, repository: DiseaseRepository, resolver: EntityResolver, client: CacheFirstProviderClient, timeout: ProviderTimeout | None = None) -> None:
        self._repository = repository
        self._resolver = resolver
        self._client = client
        self._timeout = timeout or ProviderTimeout()

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "DiseaseBankService":
        root = Path(project_root)
        store = EvidenceStore(root / "data" / "cache" / "research_evidence.sqlite3")
        provider = OpenTargetsDiseaseProvider(enabled=True)
        return cls(
            repository=DiseaseRepository(store), resolver=EntityResolver.from_project_data(root, include_mouse=False),
            client=CacheFirstProviderClient(
                cache=ProviderCache(store=store),
                executor=ProviderExecutor(provider, policy=ProviderExecutionPolicy(max_retries=0, max_concurrency=1)),
            ),
        )

    def _initialize(self) -> str | None:
        try:
            self._repository.initialize()
            return None
        except EvidenceStoreError as error:
            return str(error)

    def search(self, term: str) -> DiseaseLoad | tuple[Disease, ...]:
        text = str(term).strip()
        if not text:
            return DiseaseLoad("NO_MATCH", message="A disease search term is required.")
        error = self._initialize()
        if error:
            return DiseaseLoad("UNAVAILABLE", message=error)
        result = self._client.execute(
            ProviderRequest(provider="open_targets_disease", operation="disease_search", taxon_id=9606, canonical_query=text, provider_schema_version="v4", query_version="disease-bank-v1"),
            self._timeout, evidence_class=CacheEvidenceClass.BIOLOGICAL_ANNOTATION,
        )
        if result.status is not ProviderStatus.AVAILABLE:
            return DiseaseLoad(result.status.value, message=result.message)
        data = result.data if isinstance(result.data, Mapping) else {}
        diseases = tuple(Disease(disease_id=item["disease_id"], disease_name=item["disease_name"]) for item in data.get("diseases", ()) if isinstance(item, Mapping))
        return diseases if diseases else DiseaseLoad("NO_MATCH", message="No matching disease was returned.")

    def load_reference(self, disease_id: str) -> DiseaseLoad:
        identifier = str(disease_id).strip()
        if not identifier:
            return DiseaseLoad("NO_MATCH", message="A disease identifier is required.")
        error = self._initialize()
        if error:
            return DiseaseLoad("UNAVAILABLE", message=error)
        result = self._client.execute(
            ProviderRequest(provider="open_targets_disease", operation="disease_associations", taxon_id=9606, canonical_query=identifier, provider_schema_version="v4", query_version="disease-bank-v1", params={"limit": 500}),
            self._timeout, evidence_class=CacheEvidenceClass.BIOLOGICAL_ANNOTATION,
        )
        if result.status is not ProviderStatus.AVAILABLE:
            cached = self._repository.get(identifier)
            if cached is None:
                return DiseaseLoad(result.status.value, message=result.message)
            disease, associations = cached
            reference = self._reference(disease, associations, source_status="AVAILABLE_CACHED", source_message=result.message)
            return DiseaseLoad("AVAILABLE_CACHED", disease=disease, reference=reference, message=result.message)
        data = result.data if isinstance(result.data, Mapping) else {}
        disease_value = data.get("disease")
        if not isinstance(disease_value, Mapping):
            return DiseaseLoad("UNAVAILABLE", message="Open Targets returned no disease record.")
        retrieved_at = str(result.metadata.get("retrieved_at", "")).strip() or "1970-01-01T00:00:00+00:00"
        disease = disease_from_payload(disease_value)
        associations = associations_from_payload(data.get("associations"), retrieved_at=retrieved_at)
        try:
            self._repository.save(disease, associations, retrieved_at=retrieved_at)
        except EvidenceStoreError:
            # The loaded response remains usable this session; failed cache write is not core failure.
            pass
        reference = self._reference(disease, associations, source_status="AVAILABLE")
        return DiseaseLoad("AVAILABLE", disease=disease, reference=reference)

    def _reference(self, disease: Disease, associations, *, source_status: str, source_message: str | None = None) -> DiseaseReferenceSet:
        resolved: list[ResolvedDiseaseGene] = []
        mapped = unresolved = ambiguous = 0
        for association in associations:
            query = association.ensembl_gene_id or association.gene_symbol
            resolution = self._resolver.resolve(query, taxon_id=9606)
            # ENSG has mandatory precedence.  Symbol is considered only if that
            # exact project mapping did not resolve, never to override ambiguity.
            if resolution.status is ResolutionStatus.UNRESOLVED and association.gene_symbol and association.ensembl_gene_id:
                query = association.gene_symbol
                resolution = self._resolver.resolve(query, taxon_id=9606)
            resolved.append(ResolvedDiseaseGene(association, resolution, query))
            if resolution.entity is not None:
                mapped += 1
            elif resolution.status is ResolutionStatus.AMBIGUOUS:
                ambiguous += 1
            else:
                unresolved += 1
        return DiseaseReferenceSet(
            disease=disease, associations=tuple(associations), resolved=tuple(resolved),
            mapping_summary=MappingSummary(total=len(associations), mapped=mapped, unmapped=unresolved, ambiguous=ambiguous),
            source_status=source_status, source_message=source_message,
        )

    def resolve_perturbed_proteins(self, targets) -> tuple[str, ...]:
        """Map displayed perturbation targets through the same project resolver."""

        canonical: list[str] = []
        for target in targets:
            resolution = self._resolver.resolve(target, taxon_id=9606)
            if resolution.entity is not None:
                canonical.append(resolution.entity.canonical_id)
        return tuple(dict.fromkeys(canonical))


__all__ = ["DiseaseBankService", "DiseaseLoad"]
