from __future__ import annotations

import hashlib
from pathlib import Path
import pickle
import sqlite3

import pandas as pd

from src.product.legacy_provenance import audit_hpa_mixed_artifact, reconstruct_classic_corum


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def test_hpa_audit_is_read_only_and_separates_runtime_max_sources(tmp_path: Path):
    production = tmp_path / "production.sqlite"
    shadow = tmp_path / "shadow.sqlite"
    with sqlite3.connect(production) as connection:
        connection.execute("CREATE TABLE tissue_expression(protein_id TEXT,tissue TEXT,expression_level REAL)")
        connection.executemany("INSERT INTO tissue_expression VALUES (?,?,?)", [
            ("P1", "Lung", 3.0), ("P2", "Liver", 1.0),
            ("P1", "Lung", 2.0), ("P3", "Lung", 4.0),
        ])
    with sqlite3.connect(shadow) as connection:
        connection.execute("CREATE TABLE tissue_expression(protein_id TEXT,tissue TEXT,expression_level REAL)")
        connection.executemany("INSERT INTO tissue_expression VALUES (?,?,?)", [
            ("P1", "Lung", 2.0), ("P3", "Lung", 4.0),
        ])
    before = _hash(production)
    first = audit_hpa_mixed_artifact(production, shadow)
    second = audit_hpa_mixed_artifact(production, shadow)
    assert first == second
    assert first["runtime_max_winner"] == {"legacy": 1, "tpm": 0, "tie": 0}
    assert first["legacy_only_keys"] == 1
    assert first["tpm_suffix_vs_clean_shadow"]["semantic_exact"] is True
    assert _hash(production) == before


def test_corum_reconstruction_is_staging_only_deterministic_and_semantic(tmp_path: Path):
    raw = tmp_path / "coreComplexes.txt"
    pd.DataFrame([
        {"ComplexID": "1", "ComplexName": "Human C", "Organism": "Human", "subunits(Gene name)": "A;c12orf1"},
        {"ComplexID": "2", "ComplexName": "Mouse C", "Organism": "Mouse", "subunits(Gene name)": "B"},
        {"ComplexID": "3", "ComplexName": "Rat Missing", "Organism": "Rat", "subunits(Gene name)": None},
    ]).to_csv(raw, sep="\t", index=False)
    production = tmp_path / "complex_members.pkl"
    expected = {"A": ["Human C"], "C12ORF1": ["Human C"], "B": ["Mouse C"], "NAN": ["Rat Missing"]}
    with production.open("wb") as stream:
        pickle.dump(expected, stream)
    before = _hash(production)
    first = reconstruct_classic_corum(raw, production, tmp_path / "stage")
    artifact_hash = _hash(Path(first["artifact"]))
    second = reconstruct_classic_corum(raw, production, tmp_path / "stage")
    assert first == second
    assert _hash(Path(second["artifact"])) == artifact_hash
    assert first["semantic_equivalence"] is True
    assert first["classification"] == "LEGACY_REPRODUCIBLE"
    assert first["production_only_vs_human_builder"] == 3
    assert _hash(production) == before
