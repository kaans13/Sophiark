"""Focused tests for local-only Research Context orchestration."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
import inspect
from pathlib import Path
from types import SimpleNamespace
from tempfile import TemporaryDirectory
import time
import unittest

import pandas as pd
from pandas.testing import assert_frame_equal

from src.research.entity_resolution import EntityResolver, load_human_project_index
from src.research.local_context import build_local_research_context
from src.research.models import ProvenanceKind, ProvenanceRecord, ResolutionStatus
from src.research.providers import (
    ComplexMembership,
    FamilyMembership,
    GeneAnnotation,
    LocalEvidenceAdapters,
    LocalProviderResult,
    LocalProviderStatus,
    PartialUniProtContext,
    RegulatoryRelation,
)
from src.research.semantic_adapter import FunctionalContextKind, ResultRole, adapt_snapshot
from src.research.snapshot import build_snapshot


def _provenance(source: str) -> ProvenanceRecord:
    return ProvenanceRecord(
        source=source,
        kind=ProvenanceKind.LOCAL_ANNOTATION,
        version="fixture-v1",
        locator=f"fixture://{source}",
    )


def _result(
    provider_id: str,
    records=(),
    *,
    requested_count: int = 1,
    status: LocalProviderStatus | None = None,
    taxon_id: int = 9606,
) -> LocalProviderResult:
    records = tuple(records)
    return LocalProviderResult(
        provider_id=provider_id,
        taxon_id=taxon_id,
        status=status or (LocalProviderStatus.AVAILABLE if records else LocalProviderStatus.NO_MATCH),
        records=records,
        provenance=tuple(record.provenance for record in records),
        requested_count=requested_count,
        matched_count=len(records),
    )


class _NoMatchProvider:
    def __init__(self, provider_id: str, method_name: str) -> None:
        self.provider_id = provider_id
        setattr(self, method_name, self._lookup)

    def _lookup(self, values, *, taxon_id: int):
        return _result(
            self.provider_id,
            requested_count=len(tuple(values)),
            status=LocalProviderStatus.NO_MATCH,
            taxon_id=taxon_id,
        )


def _no_match_adapters() -> LocalEvidenceAdapters:
    return LocalEvidenceAdapters(
        mygene=_NoMatchProvider("local_mygene", "get_annotations"),
        families=_NoMatchProvider("local_hgnc_family", "get_families"),
        complexes=_NoMatchProvider("local_complex_membership", "get_memberships"),
        uniprot=_NoMatchProvider("local_uniprot_partial", "get_context"),
        regulatory=_NoMatchProvider("local_regulatory", "get_relations"),
    )


class LocalContextTests(unittest.TestCase):
    def setUp(self) -> None:
        self.temp_dir = TemporaryDirectory()
        self.root = Path(self.temp_dir.name)

    def tearDown(self) -> None:
        self.temp_dir.cleanup()

    def _resolver(self, identifiers: list[tuple[str, str]], *, external_fallback=None) -> EntityResolver:
        symbols_path = self.root / "data" / "processed" / "ensp_with_symbols.csv"
        symbols_path.parent.mkdir(parents=True, exist_ok=True)
        symbols_path.write_text(
            "gene,Symbol\n" + "".join(f"{canonical},{symbol}\n" for canonical, symbol in identifiers),
            encoding="utf-8",
        )
        return EntityResolver(
            (load_human_project_index(self.root),),
            external_fallback=external_fallback,
        )

    @staticmethod
    def _snapshot():
        report = pd.DataFrame({
            "gene": ["ENSP_TARGET", "ENSP_GAIN", "ENSP_FDR_LOSS"],
            "Symbol": ["T", "G", "FL"],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", "Üçüncül_Hasar_Stresli", "Üçüncül_Hasar_Stresli"],
            "Response_Direction": [None, "Influence Gain", "Influence Loss"],
            "Delta_PageRank_Pct": [None, 8.5, -2.5],
            "Lokalizasyon": ["Nucleus", "Cell_Membrane", "Cytosol"],
            "Gümrük_Kapisi": [True, False, True],
            "significant_redistribution": [False, True, True],
            "empirical_p": [None, 0.003, 0.004],
            "q_value": [None, 0.01, 0.02],
            "Düzenleyici_TFler": ["TF1 (Activation)", "TF2 (Repression)", "—"],
            "Hedef Gen Etkisi": ["Activation", "Repression", "—"],
            "GO_CC_Terimleri": ["nucleus", "membrane", "cytosol"],
            "GO_MF_Terimleri": ["binding", "kinase activity", "transferase activity"],
        })
        signed = pd.DataFrame({
            "gene": ["ENSP_LOSS"],
            "Symbol": ["L"],
            "Response_Direction": ["Influence Loss"],
            "Delta_PageRank_Pct": [-12.0],
            "significant_redistribution": [True],
            "q_value": [0.005],
        })
        enrichment = pd.DataFrame({
            "Pathway_Family": ["GO", "KEGG"],
            "Kaynak": ["GO_Biological_Process_2021", "KEGG_2021_Human"],
            "Term": ["membrane organization", "Signaling pathway"],
            "Adjusted P-value": [0.01, 0.02],
            "Genes": ["T;G", "G;L"],
        })
        snapshot = build_snapshot(
            report,
            signed,
            enrichment,
            species="Homo sapiens",
            taxon_id=9606,
            tissue="Pancreas",
            targets=("ENSP_TARGET",),
            attenuation=0.01,
            selection_mode="null_fdr",
            test_limit=200,
            tested_count=200,
            returned_count=3,
            top_n=None,
            threshold_parameters={"fdr_alpha": 0.05},
            created_at="2026-01-01T00:00:00Z",
        )
        return snapshot, report, signed, enrichment

    @staticmethod
    def _evidence_adapters() -> LocalEvidenceAdapters:
        mygene_prov = _provenance("local_mygene")
        family_prov = _provenance("local_hgnc_family")
        complex_prov = _provenance("local_complex_membership")
        uniprot_prov = _provenance("local_uniprot_partial")
        regulatory_prov = _provenance("local_regulatory")

        def mygene(values, *, taxon_id: int):
            record = GeneAnnotation(
                taxon_id, "ENSP_GAIN", "G", "Gain protein",
                ("response to stimulus",), ("membrane",), ("kinase activity",), mygene_prov,
            )
            return _result("local_mygene", (record,), requested_count=len(tuple(values)))

        def families(values, *, taxon_id: int):
            record = FamilyMembership(
                taxon_id, "G", "G", "HGNC:123", "Test family", "approved_symbol",
                (), None, None, family_prov,
            )
            return _result("local_hgnc_family", (record,), requested_count=len(tuple(values)))

        def complexes(values, *, taxon_id: int):
            record = ComplexMembership(
                taxon_id, "ENSP_TARGET", "ensembl_protein", "CPX-1", "Test complex", complex_prov,
            )
            return _result("local_complex_membership", (record,), requested_count=len(tuple(values)))

        def uniprot(values, *, taxon_id: int):
            record = PartialUniProtContext(
                taxon_id, "ENSP_LOSS", "P00001", {"pdb": ("1ABC",)}, {},
                ("structure",), ("ligands",), uniprot_prov,
            )
            return _result("local_uniprot_partial", (record,), requested_count=len(tuple(values)))

        def regulatory(values, *, taxon_id: int):
            record = RegulatoryRelation(
                taxon_id, "T", "G", "TRRUST", "Activation", "transcriptional",
                True, True, False, ("12345",), "TRRUST-1", {}, regulatory_prov,
            )
            return _result("local_regulatory", (record,), requested_count=len(tuple(values)))

        return LocalEvidenceAdapters(
            mygene=SimpleNamespace(get_annotations=mygene),
            families=SimpleNamespace(get_families=families),
            complexes=SimpleNamespace(get_memberships=complexes),
            uniprot=SimpleNamespace(get_context=uniprot),
            regulatory=SimpleNamespace(get_relations=regulatory),
        )

    def test_combines_semantic_roles_annotations_relationships_and_provenance(self) -> None:
        snapshot, report, signed, enrichment = self._snapshot()
        originals = tuple(frame.copy(deep=True) for frame in (report, signed, enrichment))
        fallback_calls: list[str] = []

        def forbidden_external(**kwargs):
            fallback_calls.append(kwargs["query"])
            raise AssertionError("external fallback must not run")

        resolver = self._resolver(
            [("ENSP_TARGET", "T"), ("ENSP_GAIN", "G"), ("ENSP_FDR_LOSS", "FL"), ("ENSP_LOSS", "L")],
            external_fallback=forbidden_external,
        )
        semantic = adapt_snapshot(snapshot)
        context = build_local_research_context(
            snapshot,
            resolver=resolver,
            providers=self._evidence_adapters(),
            semantic_result=semantic,
        )

        self.assertIs(context.semantic_result, semantic)
        self.assertEqual(fallback_calls, [])
        self.assertEqual(len(context.entities), 4)
        self.assertEqual(context.unresolved_entities, ())
        by_id = {item.entity.canonical_id: item for item in context.entities}
        self.assertIn(ResultRole.PERTURBATION_TARGET, by_id["ENSP_TARGET"].roles)
        self.assertIn(ResultRole.POSITIVE_REDISTRIBUTION, by_id["ENSP_GAIN"].roles)
        self.assertIn(ResultRole.FDR_SUPPORTED, by_id["ENSP_GAIN"].roles)
        self.assertIn(ResultRole.NETWORK_LOSS, by_id["ENSP_LOSS"].roles)
        self.assertIn(ResultRole.NETWORK_LOSS, by_id["ENSP_FDR_LOSS"].roles)
        self.assertIn(ResultRole.FDR_SUPPORTED, by_id["ENSP_FDR_LOSS"].roles)

        target_types = {item.annotation_type for item in by_id["ENSP_TARGET"].annotations}
        gain_types = {item.annotation_type for item in by_id["ENSP_GAIN"].annotations}
        loss_types = {item.annotation_type for item in by_id["ENSP_LOSS"].annotations}
        self.assertTrue({"localization", "gateway", "tf_regulators", "target_regulatory_effect"} <= target_types)
        self.assertTrue({
            "fdr_supported", "empirical_p", "q_value", "go_cellular_component",
            "go_molecular_function", "go_biological_process", "authoritative_family_membership",
            "directed_regulation", "go_enrichment", "kegg_enrichment",
        } <= gain_types)
        self.assertIn("partial_uniprot_context", loss_types)
        self.assertIn("complex_membership", target_types)
        self.assertEqual(
            [row.kind for row in context.functional_context],
            [FunctionalContextKind.GO, FunctionalContextKind.KEGG],
        )
        self.assertTrue(context.provenance)
        self.assertTrue(all(annotation.provenance for annotation in by_id["ENSP_GAIN"].annotations))
        self.assertIs(context.status_for("entity_resolution").status, LocalProviderStatus.AVAILABLE)
        self.assertIs(context.status_for("local_mygene").status, LocalProviderStatus.AVAILABLE)

        for original, current in zip(originals, (report, signed, enrichment)):
            assert_frame_equal(original, current)
        with self.assertRaises(TypeError):
            context.roles[by_id["ENSP_GAIN"].entity.logical_key] = ()
        with self.assertRaises(FrozenInstanceError):
            context.entities[0].roles = ()

    def test_ambiguous_symbol_remains_unresolved_and_never_calls_external(self) -> None:
        report = pd.DataFrame({
            "gene": ["ENSP_TARGET", ""],
            "Symbol": ["T", "AMB"],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", "Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [None, 1.0],
        })
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=("ENSP_TARGET",), attenuation=None, selection_mode="top_n",
            test_limit=None, tested_count=2, returned_count=2,
        )
        external_calls: list[str] = []
        resolver = self._resolver(
            [("ENSP_TARGET", "T"), ("ENSP_A", "AMB"), ("ENSP_B", "AMB")],
            external_fallback=lambda **kwargs: external_calls.append(kwargs["query"]),
        )
        context = build_local_research_context(snapshot, resolver=resolver, providers=_no_match_adapters())

        self.assertEqual(external_calls, [])
        self.assertEqual(len(context.entities), 1)
        self.assertEqual(len(context.unresolved_entities), 1)
        unresolved = context.unresolved_entities[0]
        self.assertIs(unresolved.resolution.status, ResolutionStatus.AMBIGUOUS)
        self.assertEqual({candidate.canonical_id for candidate in unresolved.resolution.candidates}, {"ENSP_A", "ENSP_B"})
        self.assertIs(context.status_for("entity_resolution").status, LocalProviderStatus.PARTIAL)
        self.assertIs(context.status_for("sophiark_signed_redistribution").status, LocalProviderStatus.UNAVAILABLE)
        self.assertIs(context.status_for("sophiark_enrichment").status, LocalProviderStatus.UNAVAILABLE)

    def test_top_n_stale_significance_flag_is_not_exposed_as_fdr_evidence(self) -> None:
        report = pd.DataFrame({
            "gene": ["ENSP_GAIN"],
            "Symbol": ["G"],
            "Hasar_Tipi": ["Üçüncül_Hasar_Stresli"],
            "Delta_PageRank_Pct": [2.0],
            "significant_redistribution": [True],
            "q_value": [0.001],
        })
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=(), attenuation=None, selection_mode="top_n",
            test_limit=None, tested_count=1, returned_count=1, top_n=1,
        )
        resolver = self._resolver([("ENSP_GAIN", "G")])
        context = build_local_research_context(snapshot, resolver=resolver, providers=_no_match_adapters())

        entity = context.entities[0]
        self.assertNotIn(ResultRole.FDR_SUPPORTED, entity.roles)
        self.assertNotIn("fdr_supported", {item.annotation_type for item in entity.annotations})
        self.assertIn("q_value", {item.annotation_type for item in entity.annotations})

    def test_provider_failures_and_cross_species_records_are_isolated(self) -> None:
        snapshot, *_ = self._snapshot()
        resolver = self._resolver(
            [("ENSP_TARGET", "T"), ("ENSP_GAIN", "G"), ("ENSP_FDR_LOSS", "FL"), ("ENSP_LOSS", "L")],
        )

        def raises(*args, **kwargs):
            raise RuntimeError("fixture provider unavailable")

        wrong_taxon = GeneAnnotation(
            10090, "ENSP_GAIN", "G", "wrong species", (), (), (), _provenance("wrong_taxon"),
        )
        providers = LocalEvidenceAdapters(
            mygene=SimpleNamespace(get_annotations=lambda values, *, taxon_id: _result(
                "local_mygene", (wrong_taxon,), requested_count=len(tuple(values)), taxon_id=taxon_id,
            )),
            families=SimpleNamespace(get_families=raises),
            complexes=SimpleNamespace(get_memberships=raises),
            uniprot=SimpleNamespace(get_context=raises),
            regulatory=SimpleNamespace(get_relations=raises),
        )
        context = build_local_research_context(snapshot, resolver=resolver, providers=providers)

        self.assertEqual(len(context.entities), 4)
        self.assertTrue(context.errors)
        self.assertIs(context.status_for("local_mygene").status, LocalProviderStatus.PROVIDER_ERROR)
        self.assertTrue(all(
            context.status_for(source).status is LocalProviderStatus.PROVIDER_ERROR
            for source in (
                "local_hgnc_family", "local_complex_membership",
                "local_uniprot_partial", "local_regulatory",
            )
        ))
        gain = next(item for item in context.entities if item.entity.canonical_id == "ENSP_GAIN")
        self.assertNotIn("gene_name", {annotation.annotation_type for annotation in gain.annotations})
        self.assertIn("localization", {annotation.annotation_type for annotation in gain.annotations})

    def test_mismatched_semantic_result_is_rejected(self) -> None:
        snapshot, *_ = self._snapshot()
        semantic = replace(adapt_snapshot(snapshot), snapshot_id="snapshot-from-another-run")
        resolver = self._resolver([("ENSP_TARGET", "T")])
        with self.assertRaisesRegex(ValueError, "does not belong"):
            build_local_research_context(
                snapshot,
                resolver=resolver,
                providers=_no_match_adapters(),
                semantic_result=semantic,
            )

    def test_module_has_no_session_ui_or_network_dependency(self) -> None:
        import src.research.local_context as module

        tree = ast.parse(inspect.getsource(module))
        imported_roots = {
            alias.name.split(".", 1)[0]
            for node in ast.walk(tree)
            if isinstance(node, (ast.Import, ast.ImportFrom))
            for alias in node.names
        }
        self.assertTrue(imported_roots.isdisjoint({
            "streamlit", "requests", "httpx", "urllib", "socket",
        }))

    def test_200_entity_local_build_budget(self) -> None:
        size = 200
        identifiers = [(f"ENSP_{index:03d}", f"S{index:03d}") for index in range(size)]
        report = pd.DataFrame({
            "gene": [canonical for canonical, _ in identifiers],
            "Symbol": [symbol for _, symbol in identifiers],
            "Hasar_Tipi": ["Üçüncül_Hasar_Stresli"] * size,
            "Response_Direction": ["Influence Gain"] * size,
            "Delta_PageRank_Pct": [1.0] * size,
            "Lokalizasyon": ["Nucleus"] * size,
            "Gümrük_Kapisi": [False] * size,
        })
        snapshot = build_snapshot(
            report, None, None, species="human", taxon_id=9606, tissue="None",
            targets=(), attenuation=None, selection_mode="top_n", test_limit=None,
            tested_count=size, returned_count=size,
        )
        resolver = self._resolver(identifiers)

        started = time.perf_counter()
        context = build_local_research_context(snapshot, resolver=resolver, providers=_no_match_adapters())
        elapsed = time.perf_counter() - started

        self.assertEqual(len(context.entities), size)
        self.assertLess(elapsed, 2.0, f"200-entity local context build took {elapsed:.3f}s")


if __name__ == "__main__":
    unittest.main()
