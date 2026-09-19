from __future__ import annotations

from datetime import datetime, timezone
import json
from pathlib import Path
import sqlite3

import pytest

from src.product.build import DataBuildError, build_string_dataset, write_build_report
from src.product.capabilities import Capability, detect_capabilities
from src.product.datasets import DataRegistry, DatasetSpec
from src.product.engine_freeze import EngineFreezeViolation, verify_engine_freeze
from src.product.entities import map_entities, normalize_string_protein_id
from src.product.fingerprints import cache_key, sha256_file
from src.product.paths import ProjectPaths
from src.product.preflight import AnalysisRequest, preflight
from src.product.provenance import build_run_manifest
from src.product.status import ProductStatus


def test_string_identifier_normalization_is_strict_and_idempotent():
    canonical = "ENSP00000354587"
    assert normalize_string_protein_id(f"9606.{canonical}") == canonical
    assert normalize_string_protein_id(canonical) == canonical
    assert normalize_string_protein_id(normalize_string_protein_id(f"9606.{canonical}")) == canonical
    with pytest.raises(ValueError, match="taxon"):
        normalize_string_protein_id(f"10090.{canonical}")
    with pytest.raises(ValueError, match="unsupported"):
        normalize_string_protein_id("TP53")


def test_mapping_report_never_hides_unmapped_records():
    values = ["9606.ENSP00000000001", "ENSP00000000002", "bad"]
    resolved, report = map_entities(values, {"ENSP00000000001": "TP53"})
    assert resolved == ["TP53", None, None]
    assert report.source_records == 3
    assert report.mapped_records == 1
    assert report.unmapped_records == 2
    assert report.usable_records == 2
    assert report.coverage == pytest.approx(1 / 3)


def _registry(root: Path, *, direction: bool = True) -> DataRegistry:
    raw = root / "data" / "raw"
    processed = root / "data" / "processed"
    raw.mkdir(parents=True)
    processed.mkdir(parents=True)
    database = raw / "core.db"
    with sqlite3.connect(database) as connection:
        connection.execute("CREATE TABLE interactions (protein1 TEXT, protein2 TEXT, combined_score INTEGER)")
        connection.execute("CREATE TABLE tissue_expression (protein_id TEXT, tissue TEXT, expression_level REAL)")
    (processed / "mapping.csv").write_text("gene,Symbol\nENSP00000000001,TP53\n", encoding="utf-8")
    if direction:
        (processed / "direction.tsv").write_text("source\ttarget\nA\tB\n", encoding="utf-8")
    return DataRegistry(root, (
        DatasetSpec("core", "fixture", "1", "human", ("PPI", "tissue_context"), True,
                    raw_location="data/raw/core.db", sqlite_tables=("interactions", "tissue_expression")),
        DatasetSpec("mapping", "fixture", "1", "human", ("canonical_mapping",), True,
                    processed_location="data/processed/mapping.csv", expected_columns=("gene", "Symbol")),
        DatasetSpec("direction", "fixture", "1", "human", ("direction",), False,
                    processed_location="data/processed/direction.tsv", expected_columns=("source", "target")),
        DatasetSpec("evidence", "fixture", "1", "human", ("evidence_channels",), False,
                    processed_location="data/processed/missing.sqlite"),
    ))


def test_registry_and_capabilities_share_one_availability_truth(tmp_path: Path):
    registry = _registry(tmp_path)
    health = registry.health_report()
    assert health["core"].status is ProductStatus.READY
    assert health["evidence"].status is ProductStatus.DATA_UNAVAILABLE
    assert registry.overall_status() is ProductStatus.PARTIAL
    capabilities = detect_capabilities(registry)
    assert capabilities.available(Capability.HAS_PPI)
    assert capabilities.available(Capability.HAS_TISSUE_CONTEXT)
    assert capabilities.available(Capability.HAS_DIRECTION)
    assert not capabilities.available(Capability.HAS_EVIDENCE_CHANNELS)


def test_missing_optional_direction_does_not_disable_classic(tmp_path: Path):
    registry = _registry(tmp_path, direction=False)
    classic = preflight(AnalysisRequest("ENSP00000000001", "Lung", "classic"), registry)
    directed = preflight(AnalysisRequest("ENSP00000000001", "Lung", "directed"), registry)
    assert classic.status is ProductStatus.READY
    assert directed.status is ProductStatus.DATA_UNAVAILABLE
    assert directed.missing == (Capability.HAS_DIRECTION,)


def test_incompatible_required_schema_is_explicit(tmp_path: Path):
    registry = _registry(tmp_path)
    (tmp_path / "data" / "processed" / "mapping.csv").write_text("wrong,column\n1,2\n", encoding="utf-8")
    health = registry.inspect("mapping")
    assert health.status is ProductStatus.INCOMPATIBLE_SCHEMA
    assert "missing columns" in health.warnings[0]
    assert registry.overall_status() is ProductStatus.DATA_UNAVAILABLE


def test_string_build_is_deterministic_and_preserves_score_scale(tmp_path: Path):
    source = tmp_path / "string.tsv"
    destination = tmp_path / "processed" / "string.csv"
    source.write_text(
        "protein1\tprotein2\tcombined_score\texperiments\n"
        "9606.ENSP00000000002\t9606.ENSP00000000003\t875\t100\n"
        "9606.ENSP00000000001\tENSP00000000002\t500\t0\n",
        encoding="utf-8",
    )
    moment = datetime(2026, 1, 1, tzinfo=timezone.utc)
    first = build_string_dataset(source, destination, dataset_version="12.0", built_at=moment)
    first_bytes = destination.read_bytes()
    second = build_string_dataset(source, destination, dataset_version="12.0", built_at=moment)
    assert destination.read_bytes() == first_bytes
    assert first.new_fingerprint == second.new_fingerprint
    assert second.previous_fingerprint == first.new_fingerprint
    text = destination.read_text(encoding="utf-8")
    assert "9606." not in text
    assert ",875," in text
    assert first.mapping_coverage == 1.0


def test_string_build_reports_mapping_coverage_and_writes_json(tmp_path: Path):
    source = tmp_path / "string.csv"
    destination = tmp_path / "processed.csv"
    source.write_text(
        "protein1,protein2,combined_score\n"
        "ENSP00000000001,ENSP00000000002,500\n"
        "ENSP00000000003,ENSP00000000004,600\n",
        encoding="utf-8",
    )
    mapping = {"ENSP00000000001": "ENSP00000000001", "ENSP00000000002": "ENSP00000000002"}
    report = build_string_dataset(source, destination, dataset_version="fixture", canonical_mapping=mapping)
    assert report.source_records == 2
    assert report.mapped_records == 1
    assert report.unmapped_records == 1
    assert report.usable_records == 1
    assert report.dropped_records == 1
    assert report.mapping_coverage == 0.5
    report_path = tmp_path / "report.json"
    write_build_report(report, report_path)
    payload = json.loads(report_path.read_text(encoding="utf-8"))
    assert payload["status"] == "READY"
    assert payload["new_fingerprint"] == sha256_file(destination)


@pytest.mark.parametrize("content, message", [
    ("protein1,protein2\nA,B\n", "missing columns"),
    ("protein1,protein2,combined_score\nENSP00000000001,ENSP00000000002,1001\n", "out of range"),
    ("protein1,protein2,combined_score\nENSP00000000001,ENSP00000000002,500\nENSP00000000002,ENSP00000000001,500\n", "duplicate"),
])
def test_string_build_fails_fast_on_bad_source(tmp_path: Path, content: str, message: str):
    source = tmp_path / "bad.csv"
    source.write_text(content, encoding="utf-8")
    with pytest.raises(DataBuildError, match=message):
        build_string_dataset(source, tmp_path / "out.csv", dataset_version="fixture")


def test_cache_key_changes_for_every_scientific_identity_dimension():
    base = {
        "dataset_fingerprints": {"string": "a"},
        "preprocessing_versions": {"string": "1"},
        "engine_version": "classic-1",
        "engine_config": {"damping": 0.85},
        "target": "ENSP1",
        "tissue": "Lung",
        "analysis_type": "classic",
    }
    reference = cache_key(**base)
    for field, replacement in {
        "dataset_fingerprints": {"string": "b"},
        "preprocessing_versions": {"string": "2"},
        "engine_version": "classic-2",
        "engine_config": {"damping": 0.9},
        "target": "ENSP2",
        "tissue": "Liver",
        "analysis_type": "directed",
    }.items():
        assert cache_key(**{**base, field: replacement}) != reference


def test_project_paths_are_rooted_and_do_not_depend_on_cwd(tmp_path: Path, monkeypatch):
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    paths = ProjectPaths.discover(tmp_path / "project")
    assert paths.raw == (tmp_path / "project" / "data" / "raw").resolve()
    assert paths.outputs == (tmp_path / "project" / "outputs").resolve()


def test_run_manifest_keeps_engine_data_and_export_versions_separate(tmp_path: Path):
    capabilities = detect_capabilities(_registry(tmp_path))
    manifest = build_run_manifest(
        target="ENSP00000000001", tissue="Lung",
        dataset_versions={"string": "12"}, dataset_fingerprints={"string": "abc"},
        capabilities=capabilities,
        timestamp=datetime(2026, 1, 1, tzinfo=timezone.utc),
    )
    assert manifest["engine_versions"]["classic"].startswith("classic-")
    assert manifest["datasets"]["string"] == {"version": "12", "fingerprint": "abc"}
    assert manifest["export_schema_version"]
    assert manifest["timestamp_utc"] == "2026-01-01T00:00:00+00:00"


def test_current_scientific_engine_sources_match_freeze_manifest():
    frozen = verify_engine_freeze()
    assert "src/biology_logic.py" in frozen
    assert "src/directed/engine.py" in frozen
    assert "src/evidence/engine.py" in frozen


def test_engine_freeze_detects_drift_without_updating_reference(tmp_path: Path):
    (tmp_path / "src").mkdir()
    source = tmp_path / "src" / "engine.py"
    source.write_text("reference\n", encoding="utf-8")
    config = tmp_path / "config"
    config.mkdir()
    manifest = config / "engine-freeze.json"
    manifest.write_text(json.dumps({"files": [{"path": "src/engine.py", "sha256": "wrong"}]}), encoding="utf-8")
    with pytest.raises(EngineFreezeViolation, match="freeze violated"):
        verify_engine_freeze(root=tmp_path, manifest_path=manifest)


def test_engine_freeze_rejects_duplicate_and_generated_entries(tmp_path: Path):
    config = tmp_path / "config"
    config.mkdir()
    manifest = config / "engine-freeze.json"
    duplicate = {"path": "src/engine.py", "sha256": "0" * 64}
    manifest.write_text(json.dumps({"files": [duplicate, duplicate]}), encoding="utf-8")
    with pytest.raises(EngineFreezeViolation, match="duplicate"):
        verify_engine_freeze(root=tmp_path, manifest_path=manifest)
    manifest.write_text(json.dumps({
        "files": [{"path": "outputs/cache/generated.py", "sha256": "0" * 64}]
    }), encoding="utf-8")
    with pytest.raises(EngineFreezeViolation, match="unsafe/non-source"):
        verify_engine_freeze(root=tmp_path, manifest_path=manifest)
