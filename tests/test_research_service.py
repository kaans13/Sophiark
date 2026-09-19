"""Focused acceptance tests for the offline Research Context service boundary."""

from __future__ import annotations

import ast
from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from tempfile import TemporaryDirectory
import time
import unittest
from unittest.mock import patch

import pandas as pd
from pandas.testing import assert_frame_equal

from src.research.config import ResearchConfig
from src.research.entity_resolution import (
    EntityResolver,
    ExternalCandidateSet,
    external_cache_key,
)
from src.research.facade import (
    OfflineResearchBundle,
    ResearchContextService,
    ServiceInitializationMode,
    build_offline_research_bundle,
    build_research_snapshot,
)
from src.research.models import Entity, EntityType
from src.research.providers import (
    LocalEvidenceAdapters,
    LocalProviderResult,
    LocalProviderStatus,
)


class _NoMatchProvider:
    def __init__(self, provider_id: str, method_name: str) -> None:
        self.provider_id = provider_id
        setattr(self, method_name, self._lookup)

    def _lookup(self, values, *, taxon_id: int) -> LocalProviderResult:
        requested = len(tuple(values))
        return LocalProviderResult(
            provider_id=self.provider_id,
            taxon_id=taxon_id,
            status=LocalProviderStatus.NO_MATCH,
            requested_count=requested,
        )


class _FailingProvider:
    def __init__(self, method_name: str) -> None:
        setattr(self, method_name, self._lookup)

    @staticmethod
    def _lookup(values, *, taxon_id: int):
        del values, taxon_id
        raise RuntimeError("fixture provider failure")


def _adapters(*, failing: bool = False) -> LocalEvidenceAdapters:
    provider = _FailingProvider if failing else _NoMatchProvider
    if failing:
        return LocalEvidenceAdapters(
            mygene=provider("get_annotations"),
            families=provider("get_families"),
            complexes=provider("get_memberships"),
            uniprot=provider("get_context"),
            regulatory=provider("get_relations"),
        )
    return LocalEvidenceAdapters(
        mygene=provider("local_mygene", "get_annotations"),
        families=provider("local_hgnc_family", "get_families"),
        complexes=provider("local_complex_membership", "get_memberships"),
        uniprot=provider("local_uniprot_partial", "get_context"),
        regulatory=provider("local_regulatory", "get_relations"),
    )


def _resolver(identifiers: list[tuple[str, str]], *, fallback=None) -> EntityResolver:
    cache = {}
    for canonical_id, symbol in identifiers:
        entity = Entity(
            taxon_id=9606,
            entity_type=EntityType.PROTEIN,
            canonical_id=canonical_id,
            symbol=symbol,
            display_name=f"Fixture {symbol}",
        )
        cache[external_cache_key(canonical_id, 9606)] = ExternalCandidateSet(
            candidates=(entity,),
            source="test_local_cache",
        )
    return EntityResolver(local_external_cache=cache, external_fallback=fallback)


def _frames() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    report = pd.DataFrame({
        "gene": ["ENSP_TARGET", "ENSP_GAIN_B", "ENSP_LOSS_REPORT", "ENSP_GAIN_A"],
        "Symbol": ["T", "GB", "LR", "GA"],
        "Hasar_Tipi": [
            "Birincil_Hasar_Hedef",
            "Üçüncül_Hasar_Stresli",
            "Üçüncül_Hasar_Stresli",
            "Üçüncül_Hasar_Stresli",
        ],
        "Response_Direction": [None, "Influence Gain", "Influence Loss", "Influence Gain"],
        "Delta_PageRank_Pct": [None, 1.0, -4.0, 9.0],
        "PageRank_Baseline": [0.2, 0.1, 0.3, 0.4],
        "Hinterland_Skoru": [50.0, 20.0, 30.0, 40.0],
        "BC_Skoru": [5.0, 2.0, 3.0, 4.0],
        "Gümrük_Kapisi": [True, False, True, False],
        "Lokalizasyon": ["Nucleus", "Membrane", "Cytosol", "Membrane"],
        "q_value": [None, 0.001, 0.002, 0.9],
        "significant_redistribution": [False, True, True, False],
    })
    # Deliberately not ordered by magnitude: the service must retain this
    # already-bounded presentation order exactly.
    signed = pd.DataFrame({
        "gene": ["ENSP_LOSS_SHALLOW", "ENSP_GAIN_SIGNED", "ENSP_LOSS_DEEP"],
        "Symbol": ["LS", "GS", "LD"],
        "Response_Direction": ["Influence Loss", "Influence Gain", "Influence Loss"],
        "Delta_PageRank_Pct": [-1.0, 7.0, -20.0],
        "q_value": [0.9, 0.01, 0.001],
        "significant_redistribution": [False, True, True],
    })
    enrichment = pd.DataFrame({
        "Pathway_Family": ["KEGG", "GO:BP"],
        "Kaynak": ["KEGG_2021_Human", "GO_Biological_Process_2021"],
        "Term": ["Fixture pathway", "Fixture process"],
        "Adjusted P-value": [0.02, 0.01],
        "Genes": ["GB;GA", "T;GB"],
    })
    return report, signed, enrichment


def _snapshot(
    report: pd.DataFrame | None = None,
    signed: pd.DataFrame | None = None,
    enrichment: pd.DataFrame | None = None,
):
    defaults = _frames()
    report = defaults[0] if report is None else report
    signed = defaults[1] if signed is None else signed
    enrichment = defaults[2] if enrichment is None else enrichment
    return build_research_snapshot(
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
        returned_count=len(report),
        threshold_parameters={"fdr_alpha": 0.05},
        created_at="2026-01-01T00:00:00Z",
    )


def _service(*, failing_providers: bool = False, fallback=None) -> ResearchContextService:
    identifiers = [
        ("ENSP_TARGET", "T"),
        ("ENSP_GAIN_B", "GB"),
        ("ENSP_LOSS_REPORT", "LR"),
        ("ENSP_GAIN_A", "GA"),
        ("ENSP_LOSS_SHALLOW", "LS"),
        ("ENSP_GAIN_SIGNED", "GS"),
        ("ENSP_LOSS_DEEP", "LD"),
    ]
    return ResearchContextService(
        resolver=_resolver(identifiers, fallback=fallback),
        providers=_adapters(failing=failing_providers),
        config=ResearchConfig(research_context_enabled=True),
    )


class ResearchContextServiceTests(unittest.TestCase):
    def test_bundle_is_complete_immutable_and_bound_to_one_snapshot(self) -> None:
        snapshot = _snapshot()
        bundle = _service().build_offline_bundle(snapshot, snapshot_build_ms=1.25)

        self.assertIsInstance(bundle, OfflineResearchBundle)
        self.assertEqual(bundle.snapshot_id, snapshot.snapshot_id)
        self.assertEqual(bundle.simulation_id, snapshot.simulation_id)
        self.assertEqual(bundle.semantic_result.snapshot_id, snapshot.snapshot_id)
        self.assertEqual(bundle.local_context.snapshot_id, snapshot.snapshot_id)
        self.assertEqual(bundle.relationships.snapshot_id, snapshot.snapshot_id)
        self.assertEqual(bundle.observation_result.snapshot_id, snapshot.snapshot_id)
        self.assertTrue(bundle.bundle_id.startswith(f"offline:{snapshot.snapshot_id[:16]}:"))
        self.assertGreaterEqual(bundle.timings.snapshot_build_ms, 1.25)
        self.assertGreaterEqual(
            bundle.timings.total_ms,
            sum((
                bundle.timings.snapshot_build_ms,
                bundle.timings.semantic_adapter_ms,
                bundle.timings.local_context_ms,
                bundle.timings.relationship_ms,
                bundle.timings.observation_ms,
            )),
        )
        self.assertTrue(bundle.healthy)
        with self.assertRaises(FrozenInstanceError):
            bundle.bundle_id = "stale"  # type: ignore[misc]
        with self.assertRaises(TypeError):
            bundle.overview()["snapshot_id"] = "stale"  # type: ignore[index]
        with self.assertRaises(ValueError):
            replace(bundle, snapshot=_snapshot(report=_frames()[0].assign(Delta_PageRank_Pct=[None, 2, -4, 9])))

    def test_pipeline_is_offline_and_service_module_has_no_forbidden_dependencies(self) -> None:
        fallback_calls: list[str] = []

        def forbidden_external(**kwargs):
            fallback_calls.append(str(kwargs.get("query")))
            raise AssertionError("external entity fallback must remain outside the offline path")

        service = _service(fallback=forbidden_external)
        with patch("socket.create_connection", side_effect=AssertionError("network forbidden")), patch(
            "urllib.request.urlopen", side_effect=AssertionError("network forbidden"),
        ):
            bundle = service.build_offline_bundle(_snapshot())
        self.assertTrue(bundle.healthy)
        self.assertEqual(fallback_calls, [])

        source_path = Path(__file__).parents[1] / "src" / "research" / "service.py"
        tree = ast.parse(source_path.read_text(encoding="utf-8"))
        imports = {
            alias.name
            for node in ast.walk(tree)
            if isinstance(node, ast.Import)
            for alias in node.names
        }
        imports.update(
            str(node.module or "")
            for node in ast.walk(tree)
            if isinstance(node, ast.ImportFrom)
        )
        self.assertFalse(any(name.startswith(("streamlit", "requests", "urllib", "httpx")) for name in imports))
        source_text = source_path.read_text(encoding="utf-8")
        for forbidden in ("session_state", "read_csv(", "read_pickle(", "pagerank(", "betweenness_centrality("):
            self.assertNotIn(forbidden, source_text)

    def test_same_snapshot_is_deterministic_and_changed_snapshot_is_distinct(self) -> None:
        service = _service()
        snapshot = _snapshot()
        first = service.build_offline_bundle(snapshot)
        second = service.build_offline_bundle(snapshot)
        self.assertEqual(first.bundle_id, second.bundle_id)
        self.assertEqual(first.semantic_result, second.semantic_result)
        self.assertEqual(
            tuple(item.assertion_id for item in first.relationships.assertions),
            tuple(item.assertion_id for item in second.relationships.assertions),
        )
        self.assertEqual(
            tuple(item.observation_id for item in first.observation_result.observations),
            tuple(item.observation_id for item in second.observation_result.observations),
        )

        changed_report = _frames()[0].copy(deep=True)
        changed_report.loc[1, "Delta_PageRank_Pct"] = 1.25
        changed_snapshot = _snapshot(report=changed_report)
        changed = service.build_offline_bundle(changed_snapshot)
        self.assertNotEqual(changed_snapshot.snapshot_id, snapshot.snapshot_id)
        self.assertNotEqual(changed.bundle_id, first.bundle_id)

    def test_inputs_are_defensively_frozen_and_signed_order_is_not_reranked(self) -> None:
        report, signed, enrichment = _frames()
        originals = tuple(frame.copy(deep=True) for frame in (report, signed, enrichment))
        snapshot = _snapshot(report, signed, enrichment)
        bundle = _service().build_offline_bundle(snapshot)

        for actual, expected in zip((report, signed, enrichment), originals):
            assert_frame_equal(actual, expected)
        self.assertEqual(
            [row.canonical_id for row in bundle.semantic_result.network_losses],
            ["ENSP_LOSS_SHALLOW", "ENSP_LOSS_DEEP"],
        )
        self.assertEqual(
            [row.source_position for row in bundle.semantic_result.network_losses],
            [0, 2],
        )
        report.loc[1, "gene"] = "MUTATED_AFTER_SNAPSHOT"
        self.assertEqual(snapshot.report.records()[1]["gene"], "ENSP_GAIN_B")

    def test_stage_and_local_provider_failures_are_fail_soft(self) -> None:
        snapshot = _snapshot()
        provider_bundle = _service(failing_providers=True).build_offline_bundle(snapshot)
        self.assertEqual(provider_bundle.snapshot_id, snapshot.snapshot_id)
        self.assertFalse(provider_bundle.healthy)
        self.assertTrue(provider_bundle.semantic_result.positive_redistribution)
        self.assertTrue(any("fixture provider failure" in item for item in provider_bundle.errors))

        service = _service()
        with patch("src.research.service.adapt_snapshot", side_effect=RuntimeError("semantic boom")):
            fallback_bundle = service.build_offline_bundle(snapshot)
        self.assertEqual(fallback_bundle.snapshot_id, snapshot.snapshot_id)
        self.assertFalse(fallback_bundle.healthy)
        self.assertEqual(fallback_bundle.semantic_result.positive_redistribution, ())
        self.assertTrue(any("semantic_adapter: RuntimeError" in item for item in fallback_bundle.errors))

    def test_config_requires_explicit_opt_in_and_all_external_flags_off(self) -> None:
        resolver = _resolver([])
        adapters = _adapters()
        with self.assertRaisesRegex(ValueError, "research_context_enabled"):
            ResearchContextService(resolver, adapters, ResearchConfig())
        with self.assertRaisesRegex(ValueError, "external flags"):
            ResearchContextService(
                resolver,
                adapters,
                ResearchConfig(research_context_enabled=True, external_context_enabled=True),
            )
        with self.assertRaisesRegex(ValueError, "provider_uniprot_enabled"):
            ResearchContextService(
                resolver,
                adapters,
                ResearchConfig(research_context_enabled=True, provider_uniprot_enabled=True),
            )

    def test_explicit_project_root_constructor_is_the_only_local_io_boundary(self) -> None:
        resolver = _resolver([])
        adapters = _adapters()
        config = ResearchConfig(research_context_enabled=True)
        with TemporaryDirectory() as directory, patch.object(
            EntityResolver, "from_project_data", return_value=resolver,
        ) as resolver_factory, patch.object(
            LocalEvidenceAdapters, "from_project_root", return_value=adapters,
        ) as adapter_factory:
            service = ResearchContextService.from_project_root(directory, config=config)
            expected_root = Path(directory).resolve()
            resolver_factory.assert_called_once_with(
                expected_root,
                include_human=True,
                include_mouse=True,
                external_fallback=None,
            )
            adapter_factory.assert_called_once_with(expected_root)
        self.assertIs(service.initialization_mode, ServiceInitializationMode.PROJECT_ROOT_READ_ONLY_IO)
        self.assertEqual(service.project_root, str(expected_root))

    def test_facade_helpers_preserve_exact_scope_and_delegate_without_stale_cache(self) -> None:
        report, signed, enrichment = _frames()
        service = _service()
        first = service.build_from_results(
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
            returned_count=4,
            threshold_parameters={"fdr_alpha": 0.05},
            created_at="2026-01-01T00:00:00Z",
        )
        second = build_offline_research_bundle(service, first.snapshot)
        self.assertEqual(first.snapshot.scope, second.snapshot.scope)
        self.assertEqual(first.bundle_id, second.bundle_id)
        self.assertGreaterEqual(first.timings.snapshot_build_ms, 0.0)

    def test_two_hundred_entities_remain_bounded_and_preserve_source_count(self) -> None:
        count = 200
        identifiers = [(f"ENSP_{index:04d}", f"G{index:04d}") for index in range(count)]
        report = pd.DataFrame({
            "gene": [item[0] for item in identifiers],
            "Symbol": [item[1] for item in identifiers],
            "Hasar_Tipi": ["Birincil_Hasar_Hedef", *("\u00dcçüncül_Hasar_Stresli",) * (count - 1)],
            "Response_Direction": [None, *("Influence Gain",) * (count - 1)],
            "Delta_PageRank_Pct": [None, *(float(index) / 10.0 for index in range(1, count))],
        })
        snapshot = build_research_snapshot(
            report,
            None,
            None,
            species="Homo sapiens",
            taxon_id=9606,
            tissue="Pancreas",
            targets=(identifiers[0][0],),
            attenuation=0.01,
            selection_mode="top_n",
            test_limit=200,
            tested_count=200,
            returned_count=200,
            top_n=200,
            created_at="2026-01-01T00:00:00Z",
        )
        service = ResearchContextService(
            _resolver(identifiers),
            _adapters(),
            ResearchConfig(research_context_enabled=True),
        )
        started = time.perf_counter()
        bundle = service.build_offline_bundle(snapshot)
        elapsed = time.perf_counter() - started
        self.assertLess(elapsed, 3.0)
        self.assertEqual(len(bundle.local_context.entities), count)
        self.assertEqual(len(bundle.semantic_result.positive_redistribution), count - 1)
        self.assertEqual(bundle.snapshot.scope.returned_count, count)
        self.assertTrue(bundle.healthy)


if __name__ == "__main__":
    unittest.main()
