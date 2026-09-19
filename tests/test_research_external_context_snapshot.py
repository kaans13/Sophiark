"""Phase 9 acceptance tests for immutable external snapshots and restore."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pytest

from src.research import external_context as external_context_module
from src.research.config import ResearchConfig
from src.research.evidence_store import EvidenceStore
from src.research.external_biology import ExternalBiologyLoad, ExternalBiologyPlan
from src.research.external_context import (
    ExternalContextSnapshot,
    ExternalContextSnapshotStore,
    ExternalSnapshotStatus,
    SnapshotRestoreStatus,
    create_external_context_snapshot,
    refresh_external_context_snapshot,
)
from src.research.literature import (
    LiteratureLoad,
    PublicationMatch,
    extract_relevant_passages,
    plan_literature_requests,
    publications_from_result,
)
from src.research.provider_runtime import ProviderRequest, ProviderResult, ProviderStatus
from src.research.query_planner import QueryPlanner
from tests.test_research_query_planner import _bundle


def _fixture_loads(*, biology_status: ProviderStatus = ProviderStatus.AVAILABLE):
    bundle = _bundle()
    query_plan = QueryPlanner().plan(bundle, expected_snapshot_id=bundle.snapshot_id)
    biology_request = ProviderRequest(
        "uniprot", "protein_record", 9606, "P43220",
        "uniprot-json-v1", query_plan.query_version,
    )
    biology_plan = ExternalBiologyPlan(
        query_plan_id=query_plan.plan_id,
        snapshot_id=query_plan.snapshot_id,
        taxon_id=9606,
        explicit_request=True,
        max_provider_requests=query_plan.budget.max_provider_requests,
        requests=(biology_request,),
        exact_accessions=("P43220",),
    )
    biology_result = ProviderResult(
        biology_request,
        biology_status,
        data={
            "records": [{
                "accession": "P43220",
                "protein_name": "GLP1R",
                "provenance": {
                    "source": "UniProtKB",
                    "source_identifier": "P43220",
                    "retrieved_at": "2026-08-26T00:00:00+00:00",
                    "provider": "uniprot",
                    "source_version": "2026_03",
                },
            }]
        } if biology_status is ProviderStatus.AVAILABLE else None,
        provider_version="2026_03",
    )
    biology_load = ExternalBiologyLoad(biology_plan, (biology_result,))

    config = ResearchConfig(
        research_context_enabled=True,
        external_context_enabled=True,
        provider_europepmc_enabled=True,
    )
    literature_plan = plan_literature_requests(
        query_plan,
        config=config,
        explicitly_requested=True,
        prior_provider_request_count=1,
    )
    literature_request = literature_plan.requests[0]
    publication_record = {
        "publication_id": "MED:123",
        "pmid": "123",
        "pmcid": "PMC123",
        "doi": "10.1000/fixture",
        "title": "Fixture pathway in pancreatic cells",
        "year": 2024,
        "journal": "Fixture Journal",
        "authors_summary": "A. Author et al.",
        "abstract": "Fixture pathway members were co-mentioned in pancreatic cells.",
        "source": "MED",
        "retrieved_at": "2026-08-26T01:00:00+00:00",
        "is_open_access": False,
        "provenance": {
            "source": "Europe PMC",
            "source_identifier": "MED:123",
            "retrieved_at": "2026-08-26T01:00:00+00:00",
            "provider": "europepmc",
            "source_version": "6.9",
        },
    }
    literature_result = ProviderResult(
        literature_request,
        ProviderStatus.AVAILABLE,
        data={
            "query": literature_plan.queries[0].to_dict(),
            "records": [publication_record],
            "matched_publication_count": 1,
        },
        provider_version="6.9",
    )
    publication = publications_from_result(literature_result)[0]
    match = PublicationMatch(
        publication,
        literature_plan.queries[0],
        extract_relevant_passages(
            publication,
            query_text=literature_plan.queries[0].query_text,
            entity_terms=("Fixture", "pancreatic"),
        ),
    )
    literature_load = LiteratureLoad(
        literature_plan,
        (literature_result,),
        (match,),
    )
    return bundle, query_plan, biology_load, literature_load


def _snapshot(*, biology_status=ProviderStatus.AVAILABLE, created_at="2026-08-26T02:00:00Z", refresh_of=None, token="capture-1"):
    bundle, query_plan, biology, literature = _fixture_loads(biology_status=biology_status)
    snapshot = create_external_context_snapshot(
        bundle.snapshot,
        query_plan,
        biology_load=biology,
        literature_load=literature,
        expected_simulation_snapshot_id=bundle.snapshot_id,
        created_at=created_at,
        refresh_of=refresh_of,
        refresh_token=token,
    )
    return bundle, query_plan, snapshot


def test_snapshot_has_required_identity_provenance_status_and_immutable_payload() -> None:
    bundle, _, snapshot = _snapshot()
    assert snapshot.context_snapshot_id.startswith("external-context-")
    assert snapshot.simulation_id == bundle.simulation_id
    assert snapshot.simulation_snapshot_id == bundle.snapshot_id
    assert snapshot.query_version
    assert snapshot.providers == ("uniprot", "europepmc")
    assert snapshot.provider_versions == {"uniprot": "2026_03", "europepmc": "6.9"}
    assert snapshot.provider_statuses == {
        "uniprot": ("AVAILABLE",),
        "europepmc": ("AVAILABLE",),
    }
    assert snapshot.literature_cutoff == "2026-08-26T01:00:00+00:00"
    assert snapshot.status is ExternalSnapshotStatus.AVAILABLE
    assert len(snapshot.provider_results()) == 2
    assert snapshot.payload["simulation_scientific_fingerprint"] == bundle.snapshot.scientific_fingerprint
    with pytest.raises(TypeError):
        snapshot.payload["changed"] = True
    with pytest.raises(FrozenInstanceError):
        snapshot.status = ExternalSnapshotStatus.PARTIAL
    assert ExternalContextSnapshot.from_dict(snapshot.to_dict()) == snapshot


def test_partial_and_empty_statuses_are_explicit() -> None:
    _, _, partial = _snapshot(biology_status=ProviderStatus.TIMEOUT)
    assert partial.status is ExternalSnapshotStatus.PARTIAL
    bundle = _bundle()
    query_plan = QueryPlanner().plan(bundle, expected_snapshot_id=bundle.snapshot_id)
    empty = create_external_context_snapshot(
        bundle.snapshot,
        query_plan,
        expected_simulation_snapshot_id=bundle.snapshot_id,
        created_at="2026-08-26T02:00:00Z",
        refresh_token="empty",
    )
    assert empty.status is ExternalSnapshotStatus.EMPTY
    assert empty.providers == ()


def test_capture_rejects_stale_cross_simulation_or_cross_plan_inputs() -> None:
    bundle, query_plan, biology, literature = _fixture_loads()
    with pytest.raises(ValueError, match="stale simulation snapshot"):
        create_external_context_snapshot(
            bundle.snapshot,
            query_plan,
            biology_load=biology,
            literature_load=literature,
            expected_simulation_snapshot_id="wrong",
        )
    stale_plan = replace(
        query_plan,
        snapshot_id="other-snapshot",
        intents=tuple(
            replace(intent, snapshot_id="other-snapshot")
            for intent in query_plan.intents
        ),
    )
    with pytest.raises(ValueError, match="query plan"):
        create_external_context_snapshot(
            bundle.snapshot,
            stale_plan,
            expected_simulation_snapshot_id=bundle.snapshot_id,
        )


def test_save_restore_is_append_only_idempotent_and_offline() -> None:
    bundle, _, snapshot = _snapshot()
    with TemporaryDirectory() as directory:
        evidence = EvidenceStore(Path(directory) / "research.sqlite")
        evidence.initialize()
        store = ExternalContextSnapshotStore(evidence)
        assert store.save(bundle.snapshot, snapshot).saved is True
        assert store.save(bundle.snapshot, snapshot).saved is False

        with patch("requests.Session.get", side_effect=AssertionError("network forbidden during restore")):
            restored = store.restore(snapshot.context_snapshot_id)
            latest = store.restore_latest(
                bundle.simulation_id,
                expected_simulation_snapshot_id=bundle.snapshot_id,
            )
        assert restored.status is SnapshotRestoreStatus.AVAILABLE
        assert restored.snapshot == snapshot
        assert latest.snapshot == snapshot
        with evidence.read_connection() as connection:
            assert connection.execute("SELECT COUNT(*) FROM simulation_snapshots").fetchone()[0] == 1
            assert connection.execute("SELECT COUNT(*) FROM external_context_snapshots").fetchone()[0] == 1
            assert connection.execute("SELECT COUNT(*) FROM publications").fetchone()[0] == 1
            match = connection.execute(
                "SELECT match_type, query_text, query_version, reason_shown, metadata_json "
                "FROM publication_matches"
            ).fetchone()
        assert match["match_type"] == "CO_MENTIONED_IN_PUBLICATION"
        assert match["query_text"]
        assert match["query_version"]
        assert match["reason_shown"].startswith("Shown because")
        assert snapshot.context_snapshot_id in match["metadata_json"]


def test_refresh_creates_new_snapshot_without_changing_old_snapshot() -> None:
    bundle, query_plan, first = _snapshot(created_at="2026-08-26T02:00:00Z", token="first")
    _, _, biology, literature = _fixture_loads()
    second = refresh_external_context_snapshot(
        first,
        bundle.snapshot,
        query_plan,
        biology_load=biology,
        literature_load=literature,
        expected_simulation_snapshot_id=bundle.snapshot_id,
        created_at="2026-08-27T02:00:00Z",
        refresh_token="second",
    )
    assert second.context_snapshot_id != first.context_snapshot_id
    assert second.refresh_of == first.context_snapshot_id
    with TemporaryDirectory() as directory:
        evidence = EvidenceStore(Path(directory) / "research.sqlite")
        evidence.initialize()
        store = ExternalContextSnapshotStore(evidence)
        assert store.save(bundle.snapshot, first).saved
        assert store.save(bundle.snapshot, second).saved
        assert store.restore(first.context_snapshot_id).snapshot == first
        assert store.restore_latest(
            bundle.simulation_id,
            expected_simulation_snapshot_id=bundle.snapshot_id,
        ).snapshot == second
        with evidence.read_connection() as connection:
            assert connection.execute("SELECT COUNT(*) FROM external_context_snapshots").fetchone()[0] == 2
            assert connection.execute("SELECT COUNT(*) FROM publication_matches").fetchone()[0] == 2


def test_restore_missing_store_is_fail_soft_and_never_creates_it() -> None:
    with TemporaryDirectory() as directory:
        path = Path(directory) / "missing.sqlite"
        result = ExternalContextSnapshotStore(EvidenceStore(path)).restore_latest(
            "sim-1",
            expected_simulation_snapshot_id="snapshot-1",
        )
        assert result.status is SnapshotRestoreStatus.UNAVAILABLE
        assert path.exists() is False


def test_restore_detects_tampered_snapshot_and_does_not_fall_back_to_live_data() -> None:
    bundle, _, snapshot = _snapshot()
    with TemporaryDirectory() as directory:
        evidence = EvidenceStore(Path(directory) / "research.sqlite")
        evidence.initialize()
        store = ExternalContextSnapshotStore(evidence)
        store.save(bundle.snapshot, snapshot)
        with evidence.transaction() as connection:
            connection.execute(
                "UPDATE external_context_snapshots SET metadata_json='{}' "
                "WHERE context_snapshot_id=?",
                (snapshot.context_snapshot_id,),
            )
        result = store.restore(snapshot.context_snapshot_id)
        assert result.status is SnapshotRestoreStatus.INVALID
        assert result.snapshot is None


def test_snapshot_payload_has_a_hard_size_limit(monkeypatch) -> None:
    bundle, query_plan, biology, literature = _fixture_loads()
    monkeypatch.setattr(external_context_module, "MAX_EXTERNAL_SNAPSHOT_BYTES", 1)
    with pytest.raises(ValueError, match="safety limit"):
        create_external_context_snapshot(
            bundle.snapshot,
            query_plan,
            biology_load=biology,
            literature_load=literature,
            expected_simulation_snapshot_id=bundle.snapshot_id,
        )
