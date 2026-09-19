from __future__ import annotations

import hashlib
import json
from pathlib import Path
import pickle
import sqlite3

import pandas as pd
import pytest

from src.product.builders.common import SourceBuildError
from src.product.builders.corum import build_corum_shadow
from src.product.builders.hpa import build_hpa_shadow
from src.product.builders.omnipath import build_omnipath_shadow
from src.product.builders.trrust import build_trrust_shadow
from src.product.datasets import DataRegistry, inspect_build_report
from src.product.status import ProductStatus


def _symbol_mapping(root: Path) -> Path:
    path = root / "symbols.csv"
    pd.DataFrame({
        "gene": ["ENSP_A", "ENSP_B", "ENSP_C"],
        "Symbol": ["A", "B", "C"],
    }).to_csv(path, index=False)
    return path


def test_hpa_shadow_is_deterministic_idempotent_and_maps_ensp(tmp_path: Path):
    source = tmp_path / "hpa.tsv"
    pd.DataFrame([
        {"Gene": "ENSG00000000001", "Gene name": "A", "Tissue": "lung", "TPM": 2.0, "pTPM": 0, "nTPM": 0},
        {"Gene": "ENSG00000000001", "Gene name": "A", "Tissue": "lung", "TPM": 3.0, "pTPM": 0, "nTPM": 0},
        {"Gene": "ENSG00000000002", "Gene name": "B", "Tissue": "liver", "TPM": 1.5, "pTPM": 0, "nTPM": 0},
        {"Gene": "ENSG00000000003", "Gene name": "X", "Tissue": "lung", "TPM": 4.0, "pTPM": 0, "nTPM": 0},
    ]).to_csv(source, sep="\t", index=False)
    mapping = tmp_path / "mapping.csv"
    pd.DataFrame({"ENSP": ["ENSP_A", "ENSP_B"], "ENSG": ["ENSG00000000001", "ENSG00000000002"]}).to_csv(mapping, index=False)
    first = build_hpa_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage", source_version="fixture-v1")
    second = build_hpa_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage", source_version="fixture-v1")
    assert first.processed_fingerprint == second.processed_fingerprint
    assert first.output_artifacts == second.output_artifacts
    assert (first.source_records, first.valid_records, first.mapped_records, first.unmapped_records) == (4, 3, 2, 1)
    assert (first.unique_source_entities, first.mapped_entities, first.unmapped_entities) == (3, 2, 1)
    frame = pd.read_csv(first.output_artifacts[0])
    assert frame.to_dict("records") == [
        {"protein_id": "ENSP_A", "tissue": "Lung", "expression_level": 3.0},
        {"protein_id": "ENSP_B", "tissue": "Liver", "expression_level": 1.5},
    ]
    with sqlite3.connect(first.output_artifacts[1]) as connection:
        assert connection.execute("SELECT COUNT(*) FROM tissue_expression").fetchone()[0] == 2
    health = inspect_build_report(Path(first.output_artifacts[0]).parent / "build-report.json")
    assert health.status is ProductStatus.READY
    assert health.mapped == 2 and health.unmapped == 1 and health.usable == 2

    changed = pd.read_csv(source, sep="\t")
    changed.loc[0, "TPM"] = 9.0
    changed.to_csv(source, sep="\t", index=False)
    third = build_hpa_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage", source_version="fixture-v2")
    assert third.processed_fingerprint != first.processed_fingerprint


def test_hpa_schema_failure_creates_no_shadow_artifact(tmp_path: Path):
    source = tmp_path / "bad.tsv"
    source.write_text("Gene\tTissue\nENSG_A\tlung\n", encoding="utf-8")
    mapping = tmp_path / "mapping.csv"
    mapping.write_text("ENSP,ENSG\nENSP_A,ENSG_A\n", encoding="utf-8")
    with pytest.raises(SourceBuildError, match="missing columns"):
        build_hpa_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage")
    assert not (tmp_path / "stage").exists()


def test_trrust_shadow_reproduces_relation_and_unknown_semantics(tmp_path: Path):
    source = tmp_path / "trrust.tsv"
    source.write_text(
        "A\tB\tActivation\t1\nA\tC\tUnknown\t2\nB\tC\tRepression\t3\n",
        encoding="utf-8",
    )
    mapping = _symbol_mapping(tmp_path)
    first = build_trrust_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage")
    second = build_trrust_shadow(source, mapping_path=mapping, output_root=tmp_path / "stage")
    assert first.processed_fingerprint == second.processed_fingerprint
    assert first.mapped_records == first.usable_records == 3
    with Path(first.output_artifacts[0]).open("rb") as handle:
        assert pickle.load(handle) == {"A": ["B", "C"], "B": ["C"]}
    with Path(first.output_artifacts[1]).open("rb") as handle:
        assert pickle.load(handle)["C"] == [("A", "Bilinmiyor"), ("B", "Baskılama")]


def test_trrust_rejects_unknown_relation_without_guessing(tmp_path: Path):
    source = tmp_path / "trrust.tsv"
    source.write_text("A\tB\tMaybe\t1\n", encoding="utf-8")
    with pytest.raises(SourceBuildError, match="unsupported relation"):
        build_trrust_shadow(source, mapping_path=_symbol_mapping(tmp_path), output_root=tmp_path / "stage")


def test_corum_shadow_filters_human_maps_members_and_never_builds_scores(tmp_path: Path):
    source = tmp_path / "corum.tsv"
    frame = pd.DataFrame([
        {"ComplexID": "1", "ComplexName": "AB complex", "Organism": "Human", "subunits(Gene name)": "A;B"},
        {"ComplexID": "1", "ComplexName": "AB complex", "Organism": "Human", "subunits(Gene name)": "A;B"},
        {"ComplexID": "2", "ComplexName": "Mouse complex", "Organism": "Mouse", "subunits(Gene name)": "C"},
        {"ComplexID": "3", "ComplexName": "AX complex", "Organism": "Human", "subunits(Gene name)": "A;X"},
    ])
    frame.to_csv(source, sep="\t", index=False)
    first = build_corum_shadow(source, mapping_path=_symbol_mapping(tmp_path), output_root=tmp_path / "stage")
    second = build_corum_shadow(source, mapping_path=_symbol_mapping(tmp_path), output_root=tmp_path / "stage")
    assert first.processed_fingerprint == second.processed_fingerprint
    memberships = pd.read_csv(first.output_artifacts[0])
    assert len(memberships) == 4
    assert set(memberships["complex_id"].astype(str)) == {"1", "3"}
    assert "score" not in " ".join(memberships.columns).casefold()
    assert first.details["excluded_records"] == 1
    assert first.mapped_entities == 2 and first.unmapped_entities == 1


def test_corum_schema_and_organism_fail_safely(tmp_path: Path):
    source = tmp_path / "corum.tsv"
    source.write_text("ComplexID\tComplexName\tOrganism\n1\tX\tHuman\n", encoding="utf-8")
    with pytest.raises(SourceBuildError, match="missing columns"):
        build_corum_shadow(source, mapping_path=_symbol_mapping(tmp_path), output_root=tmp_path / "stage")


def _omnipath_project(root: Path) -> tuple[Path, Path]:
    processed = root / "data" / "processed"
    processed.mkdir(parents=True)
    pd.DataFrame({"gene": ["ENSP_A", "ENSP_B"], "Symbol": ["A", "B"]}).to_csv(
        processed / "ensp_with_symbols.csv", index=False
    )
    pd.DataFrame({"ENSP": ["ENSP_A", "ENSP_B"], "ENSG": ["ENSG_A", "ENSG_B"]}).to_csv(
        processed / "ensp_to_ensg_map.csv", index=False
    )
    source = root / "omnipath_9606.tsv"
    pd.DataFrame([{
        "source": "P-A", "target": "P-B", "source_genesymbol": "A", "target_genesymbol": "B",
        "is_directed": True, "is_stimulation": True, "is_inhibition": False,
        "consensus_direction": True, "consensus_stimulation": True,
        "consensus_inhibition": False, "sources": "DB", "references": "PMID:1",
    }]).to_csv(source, sep="\t", index=False)
    checksum = hashlib.sha256(source.read_bytes()).hexdigest()
    source.with_name("omnipath_9606_metadata.json").write_text(json.dumps({
        "organism": "Homo sapiens", "organism_taxid": 9606,
        "release": "fixture-v1", "checksum_sha256": checksum,
    }), encoding="utf-8")
    return source, processed


def test_omnipath_shadow_preserves_direction_sign_and_identity(tmp_path: Path):
    source, _ = _omnipath_project(tmp_path)
    first = build_omnipath_shadow(
        source, project_root=tmp_path, output_root=tmp_path / "stage", minimum_rows=1,
    )
    second = build_omnipath_shadow(
        source, project_root=tmp_path, output_root=tmp_path / "stage", minimum_rows=1,
    )
    assert first.processed_fingerprint == second.processed_fingerprint
    assert first.details["directed_pairs"] == 1
    assert first.details["activating_edges"] == 1
    assert first.details["inhibitory_edges"] == 0
    assert first.mapped_records == 1


def test_omnipath_wrong_organism_and_schema_fail_before_ready(tmp_path: Path):
    source, _ = _omnipath_project(tmp_path)
    metadata = source.with_name("omnipath_9606_metadata.json")
    payload = json.loads(metadata.read_text(encoding="utf-8"))
    payload["organism_taxid"] = 10090
    metadata.write_text(json.dumps(payload), encoding="utf-8")
    with pytest.raises(ValueError, match="does not match"):
        build_omnipath_shadow(source, project_root=tmp_path, output_root=tmp_path / "stage", minimum_rows=1)


def test_extended_production_manifest_is_portable_and_parseable():
    registry = DataRegistry.from_json(
        Path(__file__).resolve().parents[1] / "config" / "datasets.json",
        root=Path(__file__).resolve().parents[1],
    )
    assert {"hpa_human_tissue", "omnipath_human", "corum_human_context", "trrust_human_tf_targets"} <= set(registry.names())
    assert registry.get("hpa_human_tissue").required is True
    assert registry.get("omnipath_human").required is False
    assert all(not Path(registry.get(name).location or "").is_absolute() for name in registry.names())
    assert registry.inspect("hpa_human_tissue").status is ProductStatus.READY


def test_shadow_health_rejects_missing_report_and_missing_artifact(tmp_path: Path):
    assert inspect_build_report(tmp_path / "missing.json").status is ProductStatus.BUILD_REQUIRED
    report = tmp_path / "build-report.json"
    report.write_text(json.dumps({
        "dataset_id": "fixture", "source_version": "1", "processing_version": "1",
        "schema_version": "1", "source_records": 1, "mapped_records": 1,
        "unmapped_records": 0, "usable_records": 1, "processed_fingerprint": "abc",
        "output_artifacts": [str(tmp_path / "absent.csv")], "status": "READY",
    }), encoding="utf-8")
    assert inspect_build_report(report).status is ProductStatus.BUILD_FAILED
