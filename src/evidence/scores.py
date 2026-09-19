from __future__ import annotations

import gzip
import hashlib
import json
import sqlite3
from pathlib import Path
from typing import Iterable, Sequence

import numpy as np

from .config import FULL_LINKS_PATH, EVIDENCE_DB_PATH, PHYSICAL_DB_PATH, STRING_VERSION


EXPECTED_COLUMNS = (
    "protein1", "protein2", "neighborhood", "neighborhood_transferred",
    "fusion", "cooccurence", "homology", "coexpression",
    "coexpression_transferred", "experiments", "experiments_transferred",
    "database", "database_transferred", "textmining",
    "textmining_transferred", "combined_score",
)
COMBINED_CHANNELS = tuple(c for c in EXPECTED_COLUMNS[2:-1] if c != "homology")
NON_TEXT_CHANNELS = tuple(c for c in COMBINED_CHANNELS if not c.startswith("textmining"))


def combine_string_channels(scores: Sequence[float], prior: float = 0.041) -> float:
    """Official STRING independent-channel combiner; inputs and output are 0..1."""
    running_no_prior = 0.0
    for score in scores:
        score = float(score)
        if score <= 0:
            continue
        no_prior = max(0.0, (score - prior) / (1.0 - prior))
        running_no_prior = 1.0 - (1.0 - running_no_prior) * (1.0 - no_prior)
    return prior + (1.0 - prior) * running_no_prior


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def inspect_local_sources() -> dict[str, object]:
    return {
        "string_full_available": FULL_LINKS_PATH.exists(),
        "string_full_path": str(FULL_LINKS_PATH),
        "string_full_bytes": FULL_LINKS_PATH.stat().st_size if FULL_LINKS_PATH.exists() else None,
        "physical_available": PHYSICAL_DB_PATH.exists(),
        "physical_path": str(PHYSICAL_DB_PATH),
        "evidence_index_available": EVIDENCE_DB_PATH.exists(),
    }


def _physical_scores() -> dict[tuple[str, str], int]:
    if not PHYSICAL_DB_PATH.exists():
        return {}
    result: dict[tuple[str, str], int] = {}
    con = sqlite3.connect(PHYSICAL_DB_PATH)
    try:
        for p1, p2, score in con.execute("SELECT protein1, protein2, combined_score FROM interactions"):
            a, b = sorted((str(p1).replace("9606.", ""), str(p2).replace("9606.", "")))
            result[(a, b)] = max(int(score), result.get((a, b), 0))
    finally:
        con.close()
    return result


def build_evidence_index(*, source: Path = FULL_LINKS_PATH, destination: Path = EVIDENCE_DB_PATH,
                         threshold: int = 500, prior: float = 0.041) -> dict[str, object]:
    """Build the minimum local index: only edges that can occur in E1/E2 at Classic threshold."""
    if not source.exists():
        raise FileNotFoundError(f"Required local STRING full resource missing: {source}")
    destination.parent.mkdir(parents=True, exist_ok=True)
    physical = _physical_scores()
    tmp = destination.with_suffix(".building.sqlite")
    if tmp.exists():
        tmp.unlink()
    con = sqlite3.connect(tmp)
    columns_sql = ", ".join(f"{c} INTEGER NOT NULL" for c in EXPECTED_COLUMNS[2:])
    con.execute(f"CREATE TABLE evidence (protein1 TEXT NOT NULL, protein2 TEXT NOT NULL, {columns_sql}, reconstructed_full REAL NOT NULL, s_nontext REAL NOT NULL, text_dependency REAL NOT NULL, physical_score INTEGER, physical_status TEXT NOT NULL, PRIMARY KEY(protein1, protein2)) WITHOUT ROWID")
    insert_cols = [*EXPECTED_COLUMNS, "reconstructed_full", "s_nontext", "text_dependency", "physical_score", "physical_status"]
    placeholders = ",".join("?" for _ in insert_cols)
    sql = f"INSERT OR REPLACE INTO evidence ({','.join(insert_cols)}) VALUES ({placeholders})"
    total = kept = 0
    errors: list[float] = []
    batch = []
    with gzip.open(source, "rt", encoding="utf-8") as stream:
        header = tuple(stream.readline().split())
        if header != EXPECTED_COLUMNS:
            raise ValueError(f"Unexpected STRING full schema: {header}")
        channel_idx = [header.index(c) for c in COMBINED_CHANNELS]
        nontext_idx = [header.index(c) for c in NON_TEXT_CHANNELS]
        for line in stream:
            total += 1
            values = line.split()
            official = int(values[-1])
            if official < threshold:
                continue
            p1, p2 = values[0].replace("9606.", ""), values[1].replace("9606.", "")
            p1, p2 = sorted((p1, p2))
            numeric = [int(v) for v in values[2:]]
            normalized = [int(values[i]) / 1000.0 for i in channel_idx]
            nontext = [int(values[i]) / 1000.0 for i in nontext_idx]
            reconstructed = combine_string_channels(normalized, prior)
            s_nontext = combine_string_channels(nontext, prior)
            errors.append(abs(reconstructed - official / 1000.0))
            dependency = max(0.0, (official / 1000.0 - s_nontext) / max(official / 1000.0, 1e-12))
            ps = physical.get((p1, p2))
            batch.append((p1, p2, *numeric, reconstructed, s_nontext, dependency, ps,
                          "SUPPORTED" if ps is not None else "NOT_SUPPORTED"))
            kept += 1
            if len(batch) >= 25_000:
                con.executemany(sql, batch); con.commit(); batch.clear()
    if batch:
        con.executemany(sql, batch); con.commit()
    con.execute("CREATE INDEX evidence_s_nontext_idx ON evidence(s_nontext)")
    err = np.asarray(errors, dtype=float)
    unique_rows = int(con.execute("SELECT COUNT(*) FROM evidence").fetchone()[0])
    stats = {
        "source": str(source), "source_sha256": _sha256(source), "string_version": STRING_VERSION,
        "threshold": threshold, "source_rows": total,
        "eligible_directional_rows": kept, "indexed_unique_rows": unique_rows,
        "pearson": None, "spearman": None,
        "median_absolute_error": float(np.median(err)), "p95_error": float(np.quantile(err, .95)),
        "p99_error": float(np.quantile(err, .99)), "max_error": float(err.max()),
        "combiner": "STRING prior-corrected independent channel combination", "prior": prior,
        "excluded_from_full": ["homology"],
        "excluded_from_s_nontext": ["textmining", "textmining_transferred", "homology"],
    }
    # Exact rank/correlation values are calculated from indexed records without retaining them in RAM.
    rows = con.execute("SELECT combined_score / 1000.0, reconstructed_full FROM evidence").fetchall()
    official = np.asarray([r[0] for r in rows]); reconstructed = np.asarray([r[1] for r in rows])
    stats["pearson"] = float(np.corrcoef(official, reconstructed)[0, 1])
    try:
        from scipy.stats import spearmanr
        stats["spearman"] = float(spearmanr(official, reconstructed).statistic)
    except Exception:
        stats["spearman"] = None
    con.execute("CREATE TABLE metadata (key TEXT PRIMARY KEY, value TEXT NOT NULL)")
    con.executemany("INSERT INTO metadata VALUES (?,?)", [(k, json.dumps(v, sort_keys=True)) for k, v in stats.items()])
    con.commit(); con.close()
    tmp.replace(destination)
    return stats


def index_metadata(path: Path = EVIDENCE_DB_PATH) -> dict[str, object]:
    con = sqlite3.connect(path)
    try:
        return {k: json.loads(v) for k, v in con.execute("SELECT key,value FROM metadata")}
    finally:
        con.close()
