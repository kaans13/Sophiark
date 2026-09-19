from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import uuid
import zipfile

import pandas as pd
import pytest
from pandas.testing import assert_frame_equal

from src.product.analysis_bundle import (
    AnalysisBundleReader, AnalysisBundleWriter, BundleError, BundleStatus,
)
from src.product.history import HistoryService, ImportStatus
from src.unified.service import UnifiedResearchReport, UnifiedStatus


def _report(target: str = "ENSP_CFTR", value: float = 1.2345678901234567) -> UnifiedResearchReport:
    candidates = pd.DataFrame({
        "entity_id": pd.Series(["ENSP_A", "ENSP_B"], dtype="string"),
        "finding_class": pd.Series(["ROBUST_CROSS_MODEL", "DIRECTION_SENSITIVE"], dtype="string"),
        "rank": pd.Series([1, pd.NA], dtype="Int64"),
        "response": pd.Series([value, -0.0], dtype="float64"),
        "available": pd.Series([True, False], dtype="boolean"),
        "evidence_status": pd.Series(["AVAILABLE", "UNAVAILABLE"], dtype="string"),
    })
    classic = pd.DataFrame({
        "gene": ["ENSP_A", "ENSP_B"], "PageRank_Baseline": [1e-12, 0.2],
        "PageRank_Perturbed": [2e-12, 0.1], "Delta_PageRank_Pct": [100.0, -50.0],
        # Mirrors live evidence columns that intentionally retain object dtype
        # even when every populated value in a particular run is numeric.
        "Evidence_combined_score": pd.Series([0.75, pd.NA], dtype="object"),
    })
    directed = pd.DataFrame({
        "entity": ["ENSP_A", "ENSP_B"], "Directed_PageRank_Baseline": [1e-12, 0.2],
        "Directed_PageRank_Perturbed": [2e-12, 0.1], "Directed_Redistribution_Pct": [100.0, -50.0],
        "Directed_BC": [0.0, 42.25],
    })
    context = pd.DataFrame({"complex_id": ["C1"], "member_count": pd.Series([2], dtype="Int64")})
    return UnifiedResearchReport(
        target=target, tissue="Lung", status=UnifiedStatus.COMPLETE,
        classic_status="COMPLETE", directed_status="AVAILABLE",
        candidates=candidates, classic_raw=classic, directed_raw=directed,
        complex_context=context,
        agreement={"status": "AVAILABLE", "signed_spearman": 0.9999999999999999},
        provenance={
            "evidence_weighting_status": "BIAS_WARNING",
            "corum_context_source": {"path": "C:/Users/person/private/coreComplexes.txt", "bytes": 12},
            "dataset_source": "C:/Users/person/private/source.tsv",
            "run_manifest": {
                "sophiark_version": "5.3.0", "engine_versions": {"classic": "v1"},
                "engine_config_fingerprints": {"classic": "cfg"},
                "datasets": {
                    "hpa_human_tissue": {"version": "UNKNOWN", "fingerprint": "hpa-old"},
                    "corum_human_legacy_lookup": {"version": "UNKNOWN", "fingerprint": "corum"},
                },
                "preprocessing_versions": {"hpa": "v1"}, "availability": {"HAS_PPI": "READY"},
                "export_schema_version": "sophiark-export-v1",
            },
        },
        context_availability={
            "evidence_provenance": {"status": "AVAILABLE"},
            "corum_complex_context": {"status": "AVAILABLE"},
        },
        errors={},
    )


def _bundle(tmp_path: Path, report: UnifiedResearchReport | None = None, *, analysis_id: str | None = None) -> Path:
    path = tmp_path / f"{uuid.uuid4().hex}.sophiark"
    AnalysisBundleWriter().write(
        report or _report(), path, analysis_id=analysis_id,
        created_at=datetime(2026, 8, 29, tzinfo=timezone.utc),
    )
    return path


def _rewrite(path: Path, changes: dict[str, bytes], additions: dict[str, bytes] | None = None) -> None:
    with zipfile.ZipFile(path, "r") as source:
        members = {name: source.read(name) for name in source.namelist()}
    members.update(changes)
    members.update(additions or {})
    with zipfile.ZipFile(path, "w", compression=zipfile.ZIP_DEFLATED) as target:
        for name, data in members.items():
            target.writestr(name, data)


def test_v1_roundtrip_preserves_tables_metadata_missing_false_and_zero(tmp_path: Path):
    original = _report()
    path = _bundle(tmp_path, original)
    saved = AnalysisBundleReader().load(path)
    assert saved.report.target == original.target and saved.report.tissue == original.tissue
    assert saved.report.status is UnifiedStatus.COMPLETE
    assert_frame_equal(saved.report.candidates, original.candidates)
    assert_frame_equal(saved.report.classic_raw, original.classic_raw)
    assert_frame_equal(saved.report.directed_raw, original.directed_raw)
    assert_frame_equal(saved.report.complex_context, original.complex_context)
    assert saved.report.agreement == original.agreement
    assert saved.report.candidates.loc[1, "available"] == False  # noqa: E712
    assert saved.report.candidates.loc[1, "response"] == 0.0
    assert pd.isna(saved.report.candidates.loc[1, "rank"])


def test_bundle_load_never_calls_engines_and_marks_dataset_mismatch(monkeypatch, tmp_path: Path):
    path = _bundle(tmp_path)
    import src.unified.service as service
    monkeypatch.setattr(service, "run_unified_analysis", lambda **_: (_ for _ in ()).throw(AssertionError("engine called")))
    saved = AnalysisBundleReader().load(path, current_dataset_fingerprints={"hpa_human_tissue": "hpa-new"})
    assert saved.report.status is UnifiedStatus.COMPLETE
    assert saved.is_historical
    assert saved.historical_dataset_mismatches["hpa_human_tissue"] == {"saved": "hpa-old", "current": "hpa-new"}


def test_corrupt_result_file_is_rejected(tmp_path: Path):
    path = _bundle(tmp_path)
    _rewrite(path, {"tables/classic_results.json": b"{}"})
    with pytest.raises(BundleError) as error:
        AnalysisBundleReader().load(path)
    assert error.value.status is BundleStatus.BUNDLE_CORRUPT


def test_unsupported_future_schema_is_rejected_before_guessing(tmp_path: Path):
    path = _bundle(tmp_path)
    with zipfile.ZipFile(path) as archive:
        manifest = json.loads(archive.read("manifest.json"))
    manifest["bundle_schema_version"] = 999
    _rewrite(path, {"manifest.json": json.dumps(manifest).encode()})
    with pytest.raises(BundleError) as error:
        AnalysisBundleReader().load(path)
    assert error.value.status is BundleStatus.UNSUPPORTED_BUNDLE_VERSION


def test_zip_path_traversal_is_rejected_without_extraction(tmp_path: Path):
    path = _bundle(tmp_path)
    outside = tmp_path.parent / "outside.txt"
    _rewrite(path, {}, {"../../outside.txt": b"unsafe"})
    with pytest.raises(BundleError) as error:
        AnalysisBundleReader().load(path)
    assert error.value.status is BundleStatus.BUNDLE_CORRUPT
    assert not outside.exists()


def test_manifest_excludes_personal_absolute_paths_and_pickle(tmp_path: Path):
    path = _bundle(tmp_path)
    with zipfile.ZipFile(path) as archive:
        assert not any(name.endswith((".pkl", ".pickle", ".joblib")) for name in archive.namelist())
        text = archive.read("report.json").decode()
        assert "C:/Users/person" not in text and "coreComplexes.txt" in text


def test_atomic_failure_leaves_no_final_or_temporary_bundle(monkeypatch, tmp_path: Path):
    final = tmp_path / "atomic.sophiark"
    monkeypatch.setattr(AnalysisBundleReader, "load", lambda *_, **__: (_ for _ in ()).throw(BundleError(BundleStatus.BUNDLE_CORRUPT, "forced")))
    with pytest.raises(BundleError):
        AnalysisBundleWriter().write(_report(), final)
    assert not final.exists()
    assert list(tmp_path.glob(".*.tmp")) == []


def test_history_index_is_disposable_and_rebuilds_semantically(tmp_path: Path):
    service = HistoryService(tmp_path / "history")
    entries = [service.save(_report(target=f"ENSP_{index}")) for index in range(3)]
    before = {(entry.analysis_id, entry.target, entry.bundle_filename) for entry in service.list_entries()}
    service.index_path.unlink()
    rebuilt = service.rebuild_history_index()
    after = {(entry.analysis_id, entry.target, entry.bundle_filename) for entry in service.list_entries()}
    assert rebuilt.valid_bundles == rebuilt.indexed == 3
    assert rebuilt.corrupt == rebuilt.duplicate_ids == 0
    assert before == after == {(entry.analysis_id, entry.target, entry.bundle_filename) for entry in entries}


def test_external_import_deduplicates_same_content_and_rejects_id_conflict(tmp_path: Path):
    identifier = str(uuid.uuid4())
    external_a = _bundle(tmp_path, _report(value=1.0), analysis_id=identifier)
    external_b = _bundle(tmp_path, _report(value=2.0), analysis_id=identifier)
    service = HistoryService(tmp_path / "history")
    first = service.import_bundle(external_a)
    same = service.import_bundle(external_a)
    conflict = service.import_bundle(external_b)
    assert first.status is ImportStatus.IMPORTED
    assert same.status is ImportStatus.ALREADY_IMPORTED
    assert conflict.status is ImportStatus.ID_CONFLICT


def test_rebuild_reports_corrupt_bundle_without_deleting_it(tmp_path: Path):
    service = HistoryService(tmp_path / "history")
    valid = _bundle(service.history_dir)
    corrupt = service.history_dir / "corrupt.sophiark"
    corrupt.write_bytes(b"not a zip")
    report = service.rebuild_history_index()
    assert report.valid_bundles == report.indexed == 1
    assert report.corrupt == 1 and report.corrupt_files == ("corrupt.sophiark",)
    assert valid.exists() and corrupt.exists()
