from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from src.research.disease_bank.benchmark import benchmark_disease_reference
from src.research.disease_bank.models import Disease, DiseaseGeneAssociation, DiseaseReferenceSet, MappingSummary, ResolvedDiseaseGene
from src.research.disease_bank.opentargets import OpenTargetsDiseaseProvider, associations_from_payload, disease_from_payload
from src.research.disease_bank.repository import DiseaseRepository
from src.research.disease_bank.service import DiseaseBankService
from src.research.entity_resolution import EntityResolver, load_human_project_index
from src.research.evidence_store import EvidenceStore
from src.research.provider_cache import CacheFirstProviderClient, ProviderCache
from src.research.provider_runtime import ProviderExecutor, ProviderRequest, ProviderResult, ProviderStatus


class _Response:
    def __init__(self, payload): self.payload = payload
    def raise_for_status(self): return None
    def json(self): return self.payload


def _post(*_args, **kwargs):
    query = kwargs["json"]["query"]
    if "search" in query:
        return _Response({"data": {"search": {"hits": [{"id": "EFO_0002508", "name": "Parkinson disease"}]}}})
    return _Response({"data": {"disease": {
        "id": "EFO_0002508", "name": "Parkinson disease", "synonyms": ["PD"], "dbXRefs": ["MONDO:0005180"],
        "parents": [{"id": "EFO_PARENT", "name": "Neurodegenerative disease"}],
        "associatedTargets": {"count": 2, "rows": [
            {"score": 0.9, "datatypeScores": [{"id": "genetic_association", "score": 0.9}], "target": {"id": "ENSG_PINK", "approvedSymbol": "PINK1", "approvedName": "PTEN induced kinase 1"}},
            {"score": 0.7, "datatypeScores": [], "target": {"id": "ENSG_LRRK", "approvedSymbol": "LRRK2", "approvedName": "leucine rich repeat kinase 2"}},
        ]},
    }}})


def _resolver(root: Path) -> EntityResolver:
    (root / "symbols.csv").write_text("gene,Symbol\nENSP_PINK,PINK1\nENSP_LRRK,LRRK2\nENSP_OTHER,OTHER\n", encoding="utf-8")
    (root / "ensp.csv").write_text("ENSP,ENSG\nENSP_PINK,ENSG_PINK\nENSP_LRRK,ENSG_LRRK\n", encoding="utf-8")
    (root / "hgnc.tsv").write_text(
        "HGNC ID\tApproved symbol\tApproved name\tStatus\tLocus type\tPrevious symbols\tAlias symbols\tChromosome\tNCBI Gene ID\tEnsembl gene ID\n"
        "HGNC:1\tPINK1\tPINK1\tApproved\tgene\t\t\t1\t1\tENSG_PINK\n"
        "HGNC:2\tLRRK2\tLRRK2\tApproved\tgene\t\t\t2\t2\tENSG_LRRK\n", encoding="utf-8")
    return EntityResolver((load_human_project_index(root, symbol_path=root / "symbols.csv", ensp_to_ensg_path=root / "ensp.csv", hgnc_path=root / "hgnc.tsv"),))


def _reference(resolver: EntityResolver, *, include_unmapped: bool = False) -> DiseaseReferenceSet:
    disease = Disease("EFO_0002508", "Parkinson disease", mondo_id="MONDO:0005180")
    associations = [
        DiseaseGeneAssociation("EFO_0002508", "PINK1", "ENSG_PINK", retrieved_at="2026-01-01T00:00:00+00:00"),
        DiseaseGeneAssociation("EFO_0002508", "LRRK2", "ENSG_LRRK", retrieved_at="2026-01-01T00:00:00+00:00"),
    ]
    if include_unmapped:
        associations.append(DiseaseGeneAssociation("EFO_0002508", "MISS", "ENSG_MISS", retrieved_at="2026-01-01T00:00:00+00:00"))
    resolved = tuple(ResolvedDiseaseGene(item, resolver.resolve(item.ensembl_gene_id, taxon_id=9606), item.ensembl_gene_id) for item in associations)
    mapped = sum(item.resolution.entity is not None for item in resolved)
    return DiseaseReferenceSet(disease, tuple(associations), resolved, MappingSummary(len(associations), mapped, len(associations) - mapped, 0), "AVAILABLE")


def test_open_targets_search_and_mondo_normalization() -> None:
    provider = OpenTargetsDiseaseProvider(post=_post)
    search = provider.execute(ProviderRequest("open_targets_disease", "disease_search", 9606, "parkinson", "v4", "v1"), timeout=__import__("src.research.config", fromlist=["ProviderTimeout"]).ProviderTimeout())
    assert search.status is ProviderStatus.AVAILABLE
    association = provider.execute(ProviderRequest("open_targets_disease", "disease_associations", 9606, "EFO_0002508", "v4", "v1", {"limit": 20}), timeout=__import__("src.research.config", fromlist=["ProviderTimeout"]).ProviderTimeout())
    disease = disease_from_payload(association.data["disease"])
    assert disease.mondo_id == "MONDO:0005180"
    assert len(associations_from_payload(association.data["associations"], retrieved_at="2026-01-01T00:00:00+00:00")) == 2


def test_gene_mapping_and_mapping_failure_are_explicit() -> None:
    with TemporaryDirectory() as directory:
        reference = _reference(_resolver(Path(directory)), include_unmapped=True)
    assert reference.mapping_summary.total == 3
    assert reference.mapping_summary.mapped == 2
    assert reference.mapping_summary.unmapped == 1


def test_benchmark_excludes_single_and_multiple_perturbation_targets_and_classifies_rows() -> None:
    with TemporaryDirectory() as directory:
        reference = _reference(_resolver(Path(directory)))
    report = pd.DataFrame({"gene": ["ENSP_PINK", "ENSP_LRRK", "ENSP_OTHER"], "Network_Loss_Pct": [0.0, -1.0, -2.0]})
    redistribution = pd.DataFrame({"gene": ["ENSP_LRRK", "ENSP_OTHER", "ENSP_PINK"]})
    benchmark = benchmark_disease_reference(reference, report=report, redistribution=redistribution, perturbed_proteins=("ENSP_PINK", "ENSP_LRRK"), population_size=100)
    assert benchmark.perturbed_disease_targets == ("ENSP_LRRK", "ENSP_PINK")
    assert benchmark.known_response == ()
    assert benchmark.known_network_loss == ()
    assert benchmark.network_candidates_outside_reference == ("ENSP_OTHER",)
    assert benchmark.recalls[0].recall is None


def test_recall_expected_overlap_and_non_reference_candidates_are_deterministic() -> None:
    with TemporaryDirectory() as directory:
        reference = _reference(_resolver(Path(directory)))
    report = pd.DataFrame({"gene": ["ENSP_LRRK", "ENSP_OTHER", "ENSP_PINK"], "Network_Loss_Pct": [-1.0, 0.0, -2.0]})
    benchmark = benchmark_disease_reference(reference, report=report, perturbed_proteins=("ENSP_PINK",), population_size=100)
    metric = benchmark.recalls[0]
    assert (metric.overlap_count, metric.evaluation_reference_count, metric.recall) == (1, 1, 1.0)
    assert metric.expected_random_overlap == 0.03
    assert benchmark.known_network_loss == ("ENSP_LRRK",)
    assert benchmark.network_candidates_outside_reference == ("ENSP_OTHER",)


def test_cache_and_unavailable_provider_fallback_never_mutate_result_frames() -> None:
    with TemporaryDirectory() as directory:
        root = Path(directory)
        store = EvidenceStore(root / "research.sqlite")
        repository = DiseaseRepository(store)
        repository.initialize()
        resolver = _resolver(root)
        provider = OpenTargetsDiseaseProvider(post=_post)
        client = CacheFirstProviderClient(cache=ProviderCache(store=store), executor=ProviderExecutor(provider))
        service = DiseaseBankService(repository=repository, resolver=resolver, client=client)
        first = service.load_reference("EFO_0002508")
        second = service.load_reference("EFO_0002508")
        assert first.status == second.status == "AVAILABLE"
        class OfflineClient:
            def execute(self, request, *_args, **_kwargs):
                return ProviderResult(request=request, status=ProviderStatus.OFFLINE, data=None, message="offline")
        fallback = DiseaseBankService(repository=repository, resolver=resolver, client=OfflineClient()).load_reference("EFO_0002508")
        assert fallback.status == "AVAILABLE_CACHED"
        report = pd.DataFrame({"gene": ["ENSP_PINK", "ENSP_LRRK"], "Network_Loss_Pct": [0.0, -1.0]})
        before = report.copy(deep=True)
        benchmark_disease_reference(first.reference, report=report, perturbed_proteins=("ENSP_PINK",))
        pd.testing.assert_frame_equal(report, before)
