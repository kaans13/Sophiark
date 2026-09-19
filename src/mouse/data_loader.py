# ---------------------------------------------------------------------------
# mouse_data_loader.py
# Fare (Mus musculus) için ayrı veritabanından (mouse_hinterland_core.db)
# ham etkileşim ve doku ifade verisinin çekilmesi.
# İnsan data_loader.py'den tamamen bağımsızdır.
# ---------------------------------------------------------------------------
import sqlite3

import numpy as np
import pandas as pd

from .config import log, BASE_DIR

# Yapısal veritabanı farede kullanılmıyor; bu sabiti burada tanımlamıyoruz.
# Gerekirse ileride eklenebilir.

# ===========================================================================
# SQL Tarafında Doku Pivotu
# ===========================================================================
def _load_tissue_pivot(con: sqlite3.Connection, hedef_doku: str) -> dict:
    """Doku ifade değerlerini doğrudan SQLite mouse_tissue_expression tablosundan yükle."""
    rows = con.execute(
        """SELECT protein_id, MAX(expression_level) AS max_exp
           FROM mouse_tissue_expression
           WHERE tissue = ? AND expression_level > 0
           GROUP BY protein_id""",
        (hedef_doku,)
    ).fetchall()
    pivot = {r[0]: float(r[1]) for r in rows}
    log.info("    [Mouse Pivot] Doku=%s → %d protein (SQL)", hedef_doku, len(pivot))
    return pivot


def _attach_tissue_weights(
    edges_df: pd.DataFrame,
    tissue_pivot: dict,
) -> pd.DataFrame:
    """Doku ağırlıklarını LEFT JOIN mantığıyla ekler (insanla aynı mantık)."""
    p1_clean = edges_df["protein1"].str.replace("10090.", "", regex=False)
    p2_clean = edges_df["protein2"].str.replace("10090.", "", regex=False)

    edges_df["exp1"] = p1_clean.map(tissue_pivot).fillna(1.0).astype(np.float32)
    edges_df["exp2"] = p2_clean.map(tissue_pivot).fillna(1.0).astype(np.float32)

    known_p1 = p1_clean.isin(tissue_pivot).sum()
    known_p2 = p2_clean.isin(tissue_pivot).sum()
    missing  = (len(edges_df) - known_p1) + (len(edges_df) - known_p2)

    log.info(
        "    [Mouse LEFT JOIN] p1 eşleşme: %d/%d | p2 eşleşme: %d/%d | "
        "Nötr exp=1.0 uygulanan kenar tarafı: %d",
        known_p1, len(edges_df),
        known_p2, len(edges_df),
        missing,
    )
    return edges_df


def filter_edges_by_tissue(
    edges_df: pd.DataFrame,
    tissue_pivot: dict,
    tissue_name: str,
) -> pd.DataFrame:
    """Yalnızca seçilen dokuda pozitif ifadesi bulunan uçları tut."""
    p1_clean = edges_df["protein1"].str.replace("10090.", "", regex=False)
    p2_clean = edges_df["protein2"].str.replace("10090.", "", regex=False)
    expressed = set(tissue_pivot)
    keep_mask = p1_clean.isin(expressed) & p2_clean.isin(expressed)
    filtered = edges_df.loc[keep_mask].copy()

    log.info(
        "    [Mouse Doku filtresi] %s: %d/%d kenar tutuldu; %d kenar "
        "en az bir ucu pozitif ifade göstermediği için çıkarıldı.",
        tissue_name,
        len(filtered),
        len(edges_df),
        len(edges_df) - len(filtered),
    )
    return filtered


# ===========================================================================
# ADIM 1: Derece Özeti
# ===========================================================================
def fetch_degree_summary(db_path: str, threshold: int) -> pd.DataFrame:
    log.info("Mouse Adım 1: SQL aggregation (elite filtre >= %d)…", threshold)
    query = """
    SELECT
        REPLACE(protein, '10090.', '') AS gene,
        COUNT(DISTINCT partner) AS k_i,
        MAX(combined_score) AS max_cs
    FROM (
        SELECT protein1 AS protein, protein2 AS partner, combined_score
        FROM mouse_interactions WHERE combined_score >= :thr
        UNION ALL
        SELECT protein2 AS protein, protein1 AS partner, combined_score
        FROM mouse_interactions WHERE combined_score >= :thr
    ) AS elite_undir
    GROUP BY REPLACE(protein, '10090.', '')
    """
    con = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(query, con, params={"thr": threshold})
    finally:
        con.close()
    df["k_i"]    = df["k_i"].astype(np.int32)
    df["max_cs"] = df["max_cs"].astype(np.int32)
    log.info("  → Mouse %d benzersiz protein.", len(df))
    return df
