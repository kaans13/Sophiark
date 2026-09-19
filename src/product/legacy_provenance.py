"""Read-only provenance reconstruction for the frozen HPA and Classic CORUM artifacts."""

from __future__ import annotations

from collections import defaultdict
import hashlib
import json
from pathlib import Path
import pickle
import sqlite3
import subprocess

import pandas as pd


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _db(path: Path) -> sqlite3.Connection:
    return sqlite3.connect(f"file:{path.resolve().as_posix()}?mode=ro", uri=True)


def inspect_hpa_git_evidence(root: str | Path) -> dict[str, object]:
    """Inspect the committed v22 source and ingestion script without restoring either file."""
    repository = Path(root)
    process = subprocess.Popen(
        ["git", "show", "ed096a5:normal_tissue.tsv"], cwd=repository,
        stdout=subprocess.PIPE,
    )
    assert process.stdout is not None
    frame = pd.read_csv(process.stdout, sep="\t", dtype=str)
    if process.wait() != 0:
        raise RuntimeError("could not read committed normal_tissue.tsv")
    level_map = {"High": 3.0, "Medium": 2.0, "Low": 1.0, "Not detected": 0.1}
    frame["value"] = frame["Level"].map(level_map).fillna(0.1)
    grouped = frame.groupby(["Gene name", "Tissue"], as_index=False)["value"].max()
    script = subprocess.check_output(
        ["git", "show", "ed096a5:db_patcher.py"], cwd=repository, text=True,
        encoding="utf-8",
    )
    return {
        "git_commit": "ed096a51e260fb2c66fa19f3f001fe1872f30e29",
        "source_path": "normal_tissue.tsv",
        "declared_source_url": "https://v22.proteinatlas.org/download/normal_tissue.tsv.zip",
        "source_rows": len(frame),
        "source_genes": int(frame["Gene name"].nunique()),
        "source_tissues": int(frame["Tissue"].nunique()),
        "grouped_rows_before_mapping": len(grouped),
        "grouped_genes_before_mapping": int(grouped["Gene name"].nunique()),
        "transformation_evidence": {
            "level_map_exact": all(token in script for token in (
                "'High': 3.0", "'Medium': 2.0", "'Low': 1.0", "'Not detected': 0.1",
            )),
            "cell_type_collapse": "groupby(['Gene name', 'Tissue'])['expression_level'].max()" in script,
            "mapping": "MyGene symbol,alias -> ensembl.protein; first returned protein",
            "tissue_normalization": "str.title()",
            "output": "hinterland_core.db:tissue_expression (replace)",
        },
        "historical_mapping_snapshot_present": False,
    }


def audit_hpa_mixed_artifact(
    production_db: str | Path,
    clean_shadow_db: str | Path,
) -> dict[str, object]:
    """Quantify the two physical row blocks and runtime MAX selection."""
    production = Path(production_db)
    shadow = Path(clean_shadow_db)
    with _db(production) as prod, _db(shadow) as clean:
        total = int(prod.execute("SELECT COUNT(*) FROM tissue_expression").fetchone()[0])
        clean_rows = int(clean.execute("SELECT COUNT(*) FROM tissue_expression").fetchone()[0])
        boundary = total - clean_rows
        suffix = {
            (protein, tissue): float(value)
            for protein, tissue, value in prod.execute(
                "SELECT protein_id,tissue,expression_level FROM tissue_expression WHERE rowid>?",
                (boundary,),
            )
        }
        shadow_values = {
            (protein, tissue): float(value)
            for protein, tissue, value in clean.execute(
                "SELECT protein_id,tissue,expression_level FROM tissue_expression"
            )
        }
        block_sql = """
            SELECT COUNT(*), COUNT(DISTINCT protein_id), COUNT(DISTINCT tissue),
                   COUNT(DISTINCT protein_id||char(0)||tissue)
            FROM tissue_expression WHERE {predicate}
        """
        legacy = prod.execute(block_sql.format(predicate="rowid<=?"), (boundary,)).fetchone()
        tpm = prod.execute(block_sql.format(predicate="rowid>?"), (boundary,)).fetchone()
        base = """
            WITH legacy AS (
                SELECT protein_id,tissue,MAX(expression_level) value
                FROM tissue_expression WHERE rowid<=? GROUP BY protein_id,tissue
            ), tpm AS (
                SELECT protein_id,tissue,expression_level value
                FROM tissue_expression WHERE rowid>?
            )
        """
        overlap = prod.execute(
            base + """SELECT COUNT(*), SUM(legacy.value>tpm.value),
                       SUM(legacy.value<tpm.value), SUM(legacy.value=tpm.value)
                       FROM legacy JOIN tpm USING(protein_id,tissue)""",
            (boundary, boundary),
        ).fetchone()
        winner_values = prod.execute(
            base + """SELECT legacy.value,COUNT(*) FROM legacy JOIN tpm
                       USING(protein_id,tissue) WHERE legacy.value>tpm.value
                       GROUP BY legacy.value ORDER BY legacy.value""",
            (boundary, boundary),
        ).fetchall()
        winner_tissues = prod.execute(
            base + """SELECT legacy.tissue,COUNT(*) FROM legacy JOIN tpm
                       USING(protein_id,tissue) WHERE legacy.value>tpm.value
                       GROUP BY legacy.tissue ORDER BY COUNT(*) DESC,legacy.tissue""",
            (boundary, boundary),
        ).fetchall()
        natural_discrete = prod.execute(
            """SELECT expression_level,COUNT(*) FROM tissue_expression
               WHERE rowid>? AND expression_level IN (0.1,1,2,3)
               GROUP BY expression_level ORDER BY expression_level""",
            (boundary,),
        ).fetchall()
        legacy_only = prod.execute(
            base + """SELECT COUNT(*) FROM legacy LEFT JOIN tpm
                       USING(protein_id,tissue) WHERE tpm.protein_id IS NULL""",
            (boundary, boundary),
        ).fetchone()[0]
        tpm_only = prod.execute(
            base + """SELECT COUNT(*) FROM tpm LEFT JOIN legacy
                       USING(protein_id,tissue) WHERE legacy.protein_id IS NULL""",
            (boundary, boundary),
        ).fetchone()[0]
    return {
        "production_sha256": _sha256(production),
        "production_rows": total,
        "physical_block_boundary": boundary,
        "legacy_block": {
            "rows": legacy[0], "entities": legacy[1], "tissues": legacy[2], "keys": legacy[3],
        },
        "tpm_block": {"rows": tpm[0], "entities": tpm[1], "tissues": tpm[2], "keys": tpm[3]},
        "tpm_suffix_vs_clean_shadow": {
            "common_keys": len(suffix.keys() & shadow_values.keys()),
            "exact_values": sum(suffix[key] == shadow_values[key] for key in suffix.keys() & shadow_values.keys()),
            "suffix_only": len(suffix.keys() - shadow_values.keys()),
            "shadow_only": len(shadow_values.keys() - suffix.keys()),
            "semantic_exact": suffix == shadow_values,
        },
        "overlapping_runtime_keys": overlap[0],
        "runtime_max_winner": {"legacy": overlap[1], "tpm": overlap[2], "tie": overlap[3]},
        "legacy_only_keys": legacy_only,
        "tpm_only_keys": tpm_only,
        "legacy_wins_by_value": {str(value): count for value, count in winner_values},
        "legacy_wins_by_tissue": {tissue: count for tissue, count in winner_tissues},
        "natural_tpm_rows_equal_to_legacy_numbers": {
            str(value): count for value, count in natural_discrete
        },
        "classification": ["LEGACY_PARTIALLY_REPRODUCIBLE", "SCIENTIFIC_MIGRATION_REQUIRED"],
        "classification_reason": "v22 source and transformations are known, but the historical MyGene response snapshot is absent",
    }


def reconstruct_classic_corum(
    source_path: str | Path,
    production_pickle: str | Path,
    output_root: str | Path,
    canonical_symbols_path: str | Path | None = None,
    broad_mapping_path: str | Path | None = None,
) -> dict[str, object]:
    """Recreate the observed all-organism, uppercase, NaN-string legacy semantics."""
    source = Path(source_path)
    production_path = Path(production_pickle)
    output = Path(output_root)
    frame = pd.read_csv(source, sep="\t", dtype=str)
    memberships: defaultdict[str, list[str]] = defaultdict(list)
    for record in frame.to_dict("records"):
        # str(NaN) is intentional: it reconstructs the observed historical "NAN" key.
        for member in str(record["subunits(Gene name)"]).split(";"):
            symbol = member.strip().upper()
            if symbol:
                memberships[symbol].append(str(record["ComplexName"]))
    reconstructed = dict(memberships)
    with production_path.open("rb") as stream:
        production = pickle.load(stream)
    all_keys = set(production) | set(reconstructed)
    set_equal = sum(
        set(production.get(key, ())) == set(reconstructed.get(key, ())) for key in all_keys
    )
    human = frame[frame["Organism"].eq("Human")]
    human_primary = {
        member.strip()
        for value in human["subunits(Gene name)"].dropna()
        for member in value.split(";") if member.strip()
    }
    production_only = set(production) - human_primary
    human_upper = {symbol.upper() for symbol in human_primary}
    nonhuman_upper = {
        member.strip().upper()
        for value in frame.loc[~frame["Organism"].eq("Human"), "subunits(Gene name)"]
        for member in str(value).split(";") if member.strip()
    }
    build_id = hashlib.sha256(
        ("classic-corum-legacy-v1:" + _sha256(source)).encode("ascii")
    ).hexdigest()[:16]
    directory = output / "corum" / build_id
    directory.mkdir(parents=True, exist_ok=True)
    artifact = directory / "classic_complex_members.pkl"
    with artifact.open("wb") as stream:
        pickle.dump(reconstructed, stream, protocol=4)
    report = {
        "source_sha256": _sha256(source),
        "production_sha256": _sha256(production_path),
        "source_rows": len(frame),
        "organisms": {str(key): int(value) for key, value in frame["Organism"].value_counts().items()},
        "production_symbols": len(production),
        "reconstructed_symbols": len(reconstructed),
        "common_symbols": len(set(production) & set(reconstructed)),
        "semantic_set_equal_keys": set_equal,
        "semantic_equivalence": set_equal == len(all_keys),
        "production_only_vs_human_builder": len(production_only),
        "production_only_explained_by_human_case_normalization": len(production_only & human_upper),
        "production_only_explained_by_nonhuman_uppercase": len((production_only - human_upper) & nonhuman_upper),
        "unexplained_production_only": sorted(production_only - human_upper - nonhuman_upper),
        "nan_key": production.get("NAN"),
        "alias_expansion_required": False,
        "classification": "LEGACY_REPRODUCIBLE",
        "artifact": str(artifact.resolve()),
    }
    if canonical_symbols_path is not None:
        canonical_frame = pd.read_csv(canonical_symbols_path)
        canonical_column = "Symbol" if "Symbol" in canonical_frame else "symbol"
        canonical = {str(value).upper() for value in canonical_frame[canonical_column].dropna()}
        report["production_only_current_canonical_symbols"] = len(production_only & canonical)
    if broad_mapping_path is not None:
        broad_frame = pd.read_parquet(broad_mapping_path)
        broad_column = "symbol" if "symbol" in broad_frame else "Symbol"
        broad = {str(value).upper() for value in broad_frame[broad_column].dropna()}
        report["production_only_present_in_broad_mapping"] = len(production_only & broad)
        report["production_only_absent_from_broad_mapping"] = len(production_only - broad)
    report_path = directory / "corum-provenance-report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True), encoding="utf-8")
    return report
