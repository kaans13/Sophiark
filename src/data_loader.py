# ---------------------------------------------------------------------------
# data_loader.py
# Veritabanından (SQLite) ham etkileşim ve doku ifade verisinin çekilmesi.
# ---------------------------------------------------------------------------
import pickle
import sqlite3

import numpy as np
import pandas as pd

from .config import _BASE, _ensp_to_symbol, log, BASE_DIR

# Yapısal veritabanı yolu
STRUCT_DB = BASE_DIR / "data" / "processed" / "structural_binding.db"


# ===========================================================================
# SQL Tarafında Doku Pivotu
# ===========================================================================
def _load_tissue_pivot(con: sqlite3.Connection, hedef_doku: str) -> dict:
    """Doku ifade değerlerini doğrudan SQLite tissue_expression tablosundan yükle."""
    rows = con.execute(
        """SELECT protein_id, MAX(expression_level) AS max_exp
           FROM tissue_expression
           WHERE tissue = ? AND expression_level > 0
           GROUP BY protein_id""",
        (hedef_doku,)
    ).fetchall()
    pivot = {r[0]: float(r[1]) for r in rows}
    log.info("    [Pivot] Doku=%s → %d protein (SQL)", hedef_doku, len(pivot))
    return pivot


def _attach_tissue_weights(
    edges_df: pd.DataFrame,
    tissue_pivot: dict,
) -> pd.DataFrame:
    """Doku ağırlıklarını LEFT JOIN mantığıyla ekler."""
    p1_clean = edges_df["protein1"].str.replace("9606.", "", regex=False)
    p2_clean = edges_df["protein2"].str.replace("9606.", "", regex=False)

    edges_df["exp1"] = p1_clean.map(tissue_pivot).fillna(1.0).astype(np.float32)
    edges_df["exp2"] = p2_clean.map(tissue_pivot).fillna(1.0).astype(np.float32)

    known_p1 = p1_clean.isin(tissue_pivot).sum()
    known_p2 = p2_clean.isin(tissue_pivot).sum()
    missing  = (len(edges_df) - known_p1) + (len(edges_df) - known_p2)

    log.info(
        "    [LEFT JOIN] p1 eşleşme: %d/%d | p2 eşleşme: %d/%d | "
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
    """Seçilen dokuda her iki ucu da pozitif ifade edilen kenarları tut."""
    p1_clean = edges_df["protein1"].str.replace("9606.", "", regex=False)
    p2_clean = edges_df["protein2"].str.replace("9606.", "", regex=False)
    expressed = set(tissue_pivot)
    keep_mask = p1_clean.isin(expressed) & p2_clean.isin(expressed)
    filtered = edges_df.loc[keep_mask].copy()
    log.info(
        "    [Doku filtresi] %s: %d/%d kenar tutuldu; %d kenar "
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
    log.info("Adım 1: SQL aggregation (elite filtre >= %d)…", threshold)
    query = """
    SELECT
        protein AS gene,
        COUNT(*) AS k_i,
        MAX(combined_score) AS max_cs
    FROM (
        SELECT protein1 AS protein, combined_score FROM interactions WHERE combined_score >= :thr
        UNION ALL
        SELECT protein2 AS protein, combined_score FROM interactions WHERE combined_score >= :thr
    ) AS elite_undir
    GROUP BY protein
    """
    con = sqlite3.connect(db_path)
    try:
        df = pd.read_sql_query(query, con, params={"thr": threshold})
    finally:
        con.close()
    df["k_i"]    = df["k_i"].astype(np.int32)
    df["max_cs"] = df["max_cs"].astype(np.int32)
    log.info("  → %d benzersiz protein.", len(df))
    return df


# ===========================================================================
# Yapısal Veritabanı Sorgu Fonksiyonu (UniProt eşleşmeli)
# ===========================================================================

def query_structural_db(ensp_list: list) -> dict:
    """
    ENSP listesine karşılık gelen yapısal, NDS ve bağlanma verilerini döndürür.
    ENSP → UniProt eşlemesini veritabanından yaparak AlphaFold ve BindingDB verilerini çeker.
    """
    if not STRUCT_DB.exists():
        log.warning("structural_binding.db bulunamadı.")
        return {}

    con = sqlite3.connect(str(STRUCT_DB))
    cur = con.cursor()

    # 1) ENSP → UniProt eşlemesini al
    placeholders = ','.join(['?'] * len(ensp_list))
    cur.execute(
        f"SELECT ensp, uniprot FROM ensp_uniprot_map WHERE ensp IN ({placeholders})",
        ensp_list
    )
    ensp_to_uniprot = dict(cur.fetchall())

    # 2) UniProt listesini oluştur
    uniprot_list = list(set(ensp_to_uniprot.values()))

    result = {ensp: {} for ensp in ensp_list}

    # 3) NDS skorları (ENSP ile)
    cur.execute(
        f"SELECT ensp, nds FROM nds_scores WHERE ensp IN ({placeholders})",
        ensp_list
    )
    for ensp, nds in cur.fetchall():
        result[ensp]['nds'] = nds

    # 4) AlphaFold yapısal verileri (UniProt ile)
    if uniprot_list:
        up_ph = ','.join(['?'] * len(uniprot_list))
        cur.execute(
            f"""SELECT uniprot, mean_plddt, sasa, largest_pocket_volume, hydrophobicity
               FROM structural_descriptors
               WHERE uniprot IN ({up_ph})""",
            uniprot_list
        )
        struct_by_uniprot = {row[0]: row[1:] for row in cur.fetchall()}

        for ensp, uniprot in ensp_to_uniprot.items():
            if uniprot in struct_by_uniprot:
                row = struct_by_uniprot[uniprot]
                result[ensp]['mean_plddt'] = row[0]
                result[ensp]['sasa'] = row[1]
                result[ensp]['pocket_vol'] = row[2]
                result[ensp]['hydrophobicity'] = row[3]

    # 5) BindingDB özeti (UniProt ile)
    if uniprot_list:
        up_ph = ','.join(['?'] * len(uniprot_list))
        cur.execute(
            f"""SELECT uniprot, num_ligands, avg_ki, min_ki, avg_ic50, min_ic50
               FROM bindingdb_summary
               WHERE uniprot IN ({up_ph})""",
            uniprot_list
        )
        binding_by_uniprot = {row[0]: row[1:] for row in cur.fetchall()}

        for ensp, uniprot in ensp_to_uniprot.items():
            if uniprot in binding_by_uniprot:
                row = binding_by_uniprot[uniprot]
                result[ensp]['binding_num_ligands'] = row[0]
                result[ensp]['binding_avg_ki'] = row[1]
                result[ensp]['binding_min_ki'] = row[2]
                result[ensp]['binding_avg_ic50'] = row[3]
                result[ensp]['binding_min_ic50'] = row[4]

    con.close()
    return result

def get_ligands_for_uniprot(uniprot_id: str) -> list:
    """
    Belirli bir UniProt ID'ye bağlanan tüm ligandların listesini döndürür.
    """
    if not STRUCT_DB.exists():
        return []
    
    con = sqlite3.connect(str(STRUCT_DB))
    cur = con.cursor()
    cur.execute(
        """SELECT ligand_name, smiles, ki_nm, ic50_nm, kd_nm, pmid, bindingdb_monomer_id
           FROM bindingdb_ligands
           WHERE uniprot = ?""",
        (uniprot_id,)
    )
    ligands = []
    for row in cur.fetchall():
        ligands.append({
            'ligand_name': row[0],
            'smiles': row[1],
            'ki_nm': row[2],
            'ic50_nm': row[3],
            'kd_nm': row[4],
            'pmid': row[5],
            'bindingdb_monomer_id': row[6]
        })
    con.close()
    return ligands



def get_similar_protein_ligands(ensp: str, top_n: int = 10) -> dict:
    """
    Hedef ENSP'ye yapısal olarak en benzeyen **ve BindingDB'de ligand kaydı
    bulunan** proteinleri bulur, ligandlarını listeler.
    """
    con = sqlite3.connect(str(STRUCT_DB))
    cur = con.cursor()
    
    # 1) Min-max normalizasyon için global min/max değerlerini al
    cur.execute("""SELECT MIN(mean_plddt), MAX(mean_plddt),
                          MIN(sasa), MAX(sasa),
                          MIN(hydrophobicity), MAX(hydrophobicity),
                          MIN(radius_gyration), MAX(radius_gyration)
                   FROM structural_descriptors""")
    mins_maxes = cur.fetchone()
    min_vals = np.array([mins_maxes[0], mins_maxes[2], mins_maxes[4], mins_maxes[6]], dtype=float)
    max_vals = np.array([mins_maxes[1], mins_maxes[3], mins_maxes[5], mins_maxes[7]], dtype=float)
    ranges = max_vals - min_vals
    ranges[ranges == 0] = 1.0
    
    # 2) Hedef proteinin ENSP → UniProt eşlemesini ve yapısal verilerini al
    cur.execute("SELECT uniprot FROM ensp_uniprot_map WHERE ensp = ?", (ensp,))
    row = cur.fetchone()
    if not row:
        con.close()
        return {'error': 'ENSP bulunamadı.'}
    target_uniprot = row[0]
    
    cur.execute(
        """SELECT mean_plddt, sasa, hydrophobicity, radius_gyration
           FROM structural_descriptors WHERE uniprot = ?""",
        (target_uniprot,)
    )
    struct_row = cur.fetchone()
    if not struct_row:
        con.close()
        return {'error': 'Yapısal veri bulunamadı.'}
    
    target_raw = np.array([struct_row[0], struct_row[1], struct_row[2], struct_row[3]], dtype=float)
    target_vector = (target_raw - min_vals) / ranges
    
    # 3) SADECE ligand kaydı olan proteinleri tara
    cur.execute(
        """SELECT sd.uniprot, em.ensp, sd.mean_plddt, sd.sasa, sd.hydrophobicity, sd.radius_gyration
           FROM structural_descriptors sd
           INNER JOIN ensp_uniprot_map em ON sd.uniprot = em.uniprot
           INNER JOIN bindingdb_summary bs ON sd.uniprot = bs.uniprot
           WHERE sd.mean_plddt IS NOT NULL 
             AND sd.sasa IS NOT NULL 
             AND sd.hydrophobicity IS NOT NULL 
             AND sd.radius_gyration IS NOT NULL
             AND bs.num_ligands > 0"""
    )
    
    similarities = []
    for row in cur.fetchall():
        uniprot, other_ensp, plddt, sasa, hydro, radius = row
        if other_ensp == ensp:
            continue
        other_raw = np.array([plddt, sasa, hydro, radius], dtype=float)
        other_vector = (other_raw - min_vals) / ranges
        dist = np.sqrt(np.sum((target_vector - other_vector) ** 2))
        sim = 1.0 / (1.0 + dist)
        similarities.append((other_ensp, uniprot, sim, plddt, sasa, hydro))
    
    # 4) En benzer top_n
    similarities.sort(key=lambda x: x[2], reverse=True)
    top_matches = similarities[:top_n]
    
    # 5) Sonuçları topla
    result = {
        'target': {
            'ensp': ensp,
            'uniprot': target_uniprot,
            'mean_plddt': struct_row[0],
            'sasa': struct_row[1],
            'hydrophobicity': struct_row[2],
            'radius_gyration': struct_row[3]
        },
        'similar_proteins': []
    }
    
    for other_ensp, uniprot, sim_score, plddt, sasa, hydro in top_matches:
        ligands = get_ligands_for_uniprot(uniprot)
        ligands_sorted = sorted(
            ligands,
            key=lambda x: (x.get('ki_nm') or 1e9, x.get('ic50_nm') or 1e9)
        )
        result['similar_proteins'].append({
            'ensp': other_ensp,
            'uniprot': uniprot,
            'similarity_score': round(sim_score, 4),
            'mean_plddt': plddt,
            'sasa': sasa,
            'hydrophobicity': hydro,
            'ligand_count': len(ligands),
            'top_ligands': ligands_sorted[:5]
        })
    
    con.close()
    return result
