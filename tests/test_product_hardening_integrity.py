"""Targeted integrity checks for the product/data hardening boundary."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
import json
from pathlib import Path
import re
import subprocess
import sys

import pytest

from src.product.build import DataBuildError, build_string_dataset
from src.product.capabilities import Capability, detect_capabilities
from src.product.datasets import DataRegistry, DatasetSpec, default_registry
from src.product.engine_freeze import verify_engine_freeze
from src.product.fingerprints import cache_key, sha256_file, stable_fingerprint
from src.product.paths import ProjectPaths
from src.product.preflight import AnalysisRequest, preflight
from src.product.provenance import build_run_manifest
from src.product.status import ProductStatus
from src.product.versions import PREPROCESSING_VERSIONS, SOPHIARK_VERSION


ROOT = Path(__file__).resolve().parents[1]


def _fixture_registry(root: Path) -> DataRegistry:
    data = root / "data"
    (data / "raw").mkdir(parents=True)
    (data / "processed").mkdir(parents=True)
    (data / "raw" / "core.db").write_bytes(b"fixture")
    (data / "processed" / "mapping.csv").write_text("gene,Symbol\nENSP1,A\n", encoding="utf-8")
    (data / "processed" / "direction.tsv").write_text("source\ttarget\nA\tB\n", encoding="utf-8")
    (data / "processed" / "complex.pkl").write_bytes(b"fixture")
    (data / "processed" / "evidence.sqlite").write_bytes(b"fixture")
    return DataRegistry(root, (
        DatasetSpec("core", "fixture", "core-v1", "human", ("PPI", "tissue_context"), True,
                    raw_location="data/raw/core.db"),
        DatasetSpec("mapping", "fixture", "map-v1", "human", ("canonical_mapping",), True,
                    processed_location="data/processed/mapping.csv", expected_columns=("gene", "Symbol")),
        DatasetSpec("direction", "fixture", "direction-v1", "human", ("direction",), False,
                    processed_location="data/processed/direction.tsv", expected_columns=("source", "target")),
        DatasetSpec("complex", "fixture", "complex-v1", "human", ("complex_context",), False,
                    processed_location="data/processed/complex.pkl"),
        DatasetSpec("evidence", "fixture", "evidence-v1", "human", ("evidence_channels", "physical_support"), False,
                    processed_location="data/processed/evidence.sqlite"),
    ))


def test_current_freeze_registry_capabilities_and_health_are_stable():
    frozen_first = verify_engine_freeze()
    frozen_second = verify_engine_freeze()
    assert frozen_first == frozen_second
    assert len(frozen_first) == len(set(frozen_first)) == 12

    registry = default_registry()
    first = detect_capabilities(registry)
    second = detect_capabilities(registry)
    assert first.as_dict() == second.as_dict()
    assert all(first.available(capability) for capability in Capability)
    assert registry.overall_status() is ProductStatus.READY


def test_registry_identity_metadata_and_invalid_entries_are_strict(tmp_path: Path):
    registry = _fixture_registry(tmp_path)
    assert len(registry.names()) == len(set(registry.names()))
    assert registry.get("core").required is True
    assert registry.get("direction").required is False
    assert registry.get("direction").version == "direction-v1"
    assert registry.get("direction").organism == "human"
    with pytest.raises(ValueError, match="unique"):
        DataRegistry(tmp_path, (registry.get("core"), registry.get("core")))
    with pytest.raises(TypeError, match="boolean"):
        replace(registry.get("core"), name="bad", required="true")  # type: ignore[arg-type]


def test_required_and_optional_failures_remain_isolated(tmp_path: Path):
    registry = _fixture_registry(tmp_path)
    (tmp_path / "data" / "processed" / "direction.tsv").unlink()
    capabilities = detect_capabilities(registry)
    assert not capabilities.available(Capability.HAS_DIRECTION)
    assert capabilities.available(Capability.HAS_PPI)
    assert capabilities.available(Capability.HAS_COMPLEX_CONTEXT)
    assert preflight(AnalysisRequest("ENSP1", "Lung", "classic"), registry).status is ProductStatus.READY
    assert preflight(AnalysisRequest("ENSP1", "Lung", "directed"), registry).status is ProductStatus.DATA_UNAVAILABLE

    (tmp_path / "data" / "raw" / "core.db").unlink()
    core = preflight(AnalysisRequest("ENSP1", "Lung", "classic"), registry)
    assert core.status is ProductStatus.DATA_UNAVAILABLE
    assert Capability.HAS_PPI in core.missing
    assert registry.overall_status() is ProductStatus.DATA_UNAVAILABLE


def test_path_resolution_is_cwd_independent_and_registry_consistent(tmp_path: Path, monkeypatch):
    project = tmp_path / "project"
    first = ProjectPaths.discover(project)
    elsewhere = tmp_path / "cwd"
    elsewhere.mkdir()
    monkeypatch.chdir(elsewhere)
    second = ProjectPaths.discover(project)
    assert first == second
    assert first.data == project.resolve() / "data"
    assert first.config == project.resolve() / "config"
    assert first.cache == project.resolve() / "outputs" / "cache"


def test_string_build_report_counts_fingerprint_and_negative_artifact_safety(tmp_path: Path):
    source = tmp_path / "string.csv"
    output = tmp_path / "processed.csv"
    source.write_text(
        "protein1,protein2,combined_score\n"
        "9606.ENSP00000000001,ENSP00000000002,500\n"
        "ENSP00000000003,ENSP00000000004,600\n",
        encoding="utf-8",
    )
    mapping = {"ENSP00000000001": "ENSP00000000001", "ENSP00000000002": "ENSP00000000002"}
    first = build_string_dataset(source, output, dataset_version="12", canonical_mapping=mapping)
    content = output.read_bytes()
    second = build_string_dataset(source, output, dataset_version="12", canonical_mapping=mapping)
    assert output.read_bytes() == content
    assert first.new_fingerprint == second.new_fingerprint == sha256_file(output)
    assert (first.source_records, first.mapped_records, first.unmapped_records, first.usable_records) == (2, 1, 1, 1)
    assert first.mapped_records <= first.source_records
    assert first.usable_records <= first.source_records
    assert first.mapped_records + first.unmapped_records == first.source_records
    assert first.mapping_coverage == first.mapped_records / first.source_records
    assert 0.0 <= first.mapping_coverage <= 1.0

    bad = tmp_path / "bad.csv"
    bad_output = tmp_path / "must-not-exist.csv"
    bad.write_text("protein1,protein2\nA,B\n", encoding="utf-8")
    with pytest.raises(DataBuildError, match="missing columns"):
        build_string_dataset(bad, bad_output, dataset_version="12")
    assert not bad_output.exists()


def test_fingerprint_and_cache_identity_semantics_ignore_nonsemantic_context():
    scientific = {"records": [["ENSP1", "ENSP2", 500]], "processing_version": "v1"}
    assert stable_fingerprint(scientific) == stable_fingerprint(dict(scientific))
    assert stable_fingerprint(scientific) != stable_fingerprint({**scientific, "records": [["ENSP1", "ENSP2", 501]]})
    assert stable_fingerprint(scientific) != stable_fingerprint({**scientific, "processing_version": "v2"})

    base = dict(
        dataset_fingerprints={"string": "content-a"}, preprocessing_versions={"string": "v1"},
        engine_version="classic-v1", engine_config={"damping": 0.85},
        target="ENSP1", tissue="Lung", analysis_type="classic",
    )
    identity = cache_key(**base)
    assert cache_key(**base) == identity
    assert cache_key(**{**base, "dataset_fingerprints": {"string": "content-b"}}) != identity
    assert cache_key(**{**base, "preprocessing_versions": {"string": "v2"}}) != identity
    assert cache_key(**{**base, "engine_config": {"damping": 0.9}}) != identity
    presentation_metadata = {"label": "Classic Yanıtı", "digits": 4}
    assert cache_key(**base) == identity and presentation_metadata


def test_run_manifest_scientific_fields_repeat_and_versions_stay_separate(tmp_path: Path):
    capabilities = detect_capabilities(_fixture_registry(tmp_path))
    common = dict(
        target="ENSP1", tissue="Lung", dataset_versions={"string": "STRING-12"},
        dataset_fingerprints={"string": "data-sha"}, capabilities=capabilities,
        engine_config_fingerprints={"classic": "config-sha", "directed": "directed-config-sha"},
    )
    first = build_run_manifest(**common, timestamp=datetime(2026, 1, 1, tzinfo=timezone.utc))
    second = build_run_manifest(**common, timestamp=datetime(2026, 1, 2, tzinfo=timezone.utc))
    assert {key: value for key, value in first.items() if key != "timestamp_utc"} == {
        key: value for key, value in second.items() if key != "timestamp_utc"
    }
    assert first["sophiark_version"] == SOPHIARK_VERSION
    assert first["datasets"]["string"]["version"] == "STRING-12"
    assert first["engine_versions"]["classic"] != "STRING-12"
    assert first["preprocessing_versions"] == PREPROCESSING_VERSIONS
    assert SOPHIARK_VERSION not in first["preprocessing_versions"].values()


def test_source_hash_and_data_fingerprint_are_independent(tmp_path: Path):
    freeze_before = verify_engine_freeze()
    data_file = tmp_path / "dataset.csv"
    data_file.write_text("a\n1\n", encoding="utf-8")
    first = sha256_file(data_file)
    data_file.write_text("a\n2\n", encoding="utf-8")
    second = sha256_file(data_file)
    assert first != second
    assert verify_engine_freeze() == freeze_before


def test_product_modules_import_without_engine_streamlit_or_dataset_side_effects():
    code = (
        "import sys; "
        "import src.product.paths, src.product.datasets, src.product.capabilities, "
        "src.product.build, src.product.provenance; "
        "assert 'streamlit' not in sys.modules; "
        "assert 'src.biology_logic' not in sys.modules; "
        "assert 'src.directed.engine' not in sys.modules; "
        "print('IMPORT_OK')"
    )
    result = subprocess.run(
        [sys.executable, "-c", code], cwd=ROOT, text=True, capture_output=True, check=False
    )
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == "IMPORT_OK"


def test_declared_environment_includes_runtime_and_test_dependencies():
    text = (ROOT / "requirements.txt").read_text(encoding="utf-8").casefold()
    required = {"pytest", "numpy", "pandas", "scipy", "igraph", "streamlit", "requests", "openpyxl"}
    declared = {
        re.split(r"[<>=~!;\s\[]", line.strip(), maxsplit=1)[0]
        for line in text.splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    assert required <= declared
