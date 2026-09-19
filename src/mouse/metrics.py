# ---------------------------------------------------------------------------
# mouse_metrics.py
# Fare (Mus musculus) için Hinterland skoru, sınıflandırma ve Gümrük Kapısı
# tespiti. İnsan metrics.py ile aynı matematiksel prensipleri kullanır,
# ancak fare config ve fare graf motoruna bağlıdır.
# ---------------------------------------------------------------------------
import gc
import math
import random
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import igraph as ig
import numpy as np
import pandas as pd
from scipy.stats import rankdata

from .config import (
    NOISE_MAX_DEGREE,
    NOISE_MAX_PR_PERCENTILE,
    NOISE_MAX_SCORE,
    THRESHOLD_BRIDGE,
    THRESHOLD_EMPEROR,
    THRESHOLD_STRATEGIC,
    WBI_BC_PERCENTILE,
    WBI_EXT_RATIO_MIN,
    log,
)
from .graph_engine import _mask_inf_distances

# ===========================================================================
# ADIM 2c: Hinterland Skoru
# ===========================================================================

def build_hinterland_score(df: pd.DataFrame, pr: dict) -> pd.DataFrame:
    log.info("Mouse Adım 2c: Hinterland_Skoru (rank-based, PR:0.60 / deg:0.40)…")
    df = df.copy()
    df["log_degree"]      = np.log10(df["k_i"] + 1)
    df["pagerank_raw"]    = df["gene"].map(pr).fillna(0.0)
    n = len(df)
    df["pagerank_norm"]   = rankdata(df["pagerank_raw"].values, method="average") / n
    df["log_degree_norm"] = rankdata(df["log_degree"].values,   method="average") / n
    df["Hinterland_Skoru"] = (
        0.60 * df["pagerank_norm"] + 0.40 * df["log_degree_norm"]
    ) * 100.0
    log.info(
        "  → Mouse min=%.2f max=%.2f ort=%.2f",
        df["Hinterland_Skoru"].min(),
        df["Hinterland_Skoru"].max(),
        df["Hinterland_Skoru"].mean(),
    )
    return df


# ===========================================================================
# ADIM 3a: Sınıflandırma
# ===========================================================================

def classify_genes(df: pd.DataFrame) -> pd.DataFrame:
    log.info(
        "Mouse Adım 3a: Sınıflandırma (İmparator>=%d, Dağıtıcı>=%d)…",
        THRESHOLD_EMPEROR, THRESHOLD_STRATEGIC,
    )
    pr_pct_cutoff = float(
        np.percentile(df["pagerank_norm"].values, NOISE_MAX_PR_PERCENTILE)
    )

    def _label(row) -> str:
        if (
            row["k_i"] <= NOISE_MAX_DEGREE
            and row["max_cs"] < NOISE_MAX_SCORE
            and row["pagerank_norm"] <= pr_pct_cutoff
        ):
            return "Gri (İstatistiksel Gürültü)"
        s = row["Hinterland_Skoru"]
        if s >= THRESHOLD_EMPEROR:   return "Koyu Kırmızı (İmparator Hub)"
        if s >= THRESHOLD_STRATEGIC: return "Açık Kırmızı (Stratejik Dağıtıcı)"
        if s >= THRESHOLD_BRIDGE:    return "Sarı (Aktif Geçiş Yolu)"
        return "Mavi (İzole Sınır Polisi)"

    df = df.copy()
    df["Kategori"] = df.apply(_label, axis=1)
    for cat, cnt in df["Kategori"].value_counts().items():
        log.info("    Mouse %-45s : %7d  (%.2f%%)", cat, cnt, cnt / len(df) * 100)
    return df


# ===========================================================================
# GÖREV 4: Gümrük Kapısı Tespiti
# ===========================================================================

def detect_customs_gates(
    g: ig.Graph,
    df: pd.DataFrame,
    bc_sample_sources: int | None = None,
) -> pd.DataFrame:
    """Hibrit WBI Gümrük Kapısı Tespiti (fare için kompartıman-bazlı BC persentili)."""
    log.info(
        "Mouse Görev 4: Hibrit WBI Gümrük Kapısı "
        "(dış-oran>%.2f AND BC>%.0f.persentil [kompartıman-bazlı], weights=distance)…",
        WBI_EXT_RATIO_MIN, WBI_BC_PERCENTILE,
    )

    g_masked = g.copy()
    _mask_inf_distances(g_masked)

    n_vertices = g_masked.vcount()
    use_sample = (
        bc_sample_sources is not None
        and 0 < bc_sample_sources < n_vertices
    )
    if use_sample:
        source_count = min(int(bc_sample_sources), n_vertices)
        bc_random_seed = 42
        sources = np.random.default_rng(bc_random_seed).choice(
            n_vertices, size=source_count, replace=False
        ).tolist()
        log.info(
            "  Mouse BC hesaplanıyor — hızlı örneklem: %d/%d kaynak, "
            "weights='distance'…", source_count, n_vertices,
        )
        try:
            bc_values = g_masked.betweenness(
                weights="distance", directed=False, sources=sources,
            )
            bc_values = [v * n_vertices / source_count for v in bc_values]
            bc_mode = "Fast Approximate"
        except TypeError:
            log.warning(
                "  Bu igraph sürümü kaynak örneklemli BC'yi desteklemiyor; "
                "tam BC hesaplanıyor."
            )
            bc_values = g_masked.betweenness(weights="distance", directed=False)
            source_count = n_vertices
            bc_mode = "Full Validation"
    else:
        log.info("  Mouse BC hesaplanıyor — tam, weights='distance' + Inf maskeleme…")
        bc_values = g_masked.betweenness(weights="distance", directed=False)
        source_count = n_vertices
        bc_random_seed = 42
        bc_mode = "Full Validation"

    bc_arr       = np.array(bc_values)
    bc           = {g.vs[i]["name"]: bc_values[i] for i in range(g.vcount())}

    bc_threshold_global = float(np.percentile(bc_arr, WBI_BC_PERCENTILE))
    bc_95_threshold = float(np.percentile(bc_arr, 95.0)) if len(bc_arr) > 0 else float("inf")

    # Kompartıman-bazlı BC persentili hesapla
    bc_by_gene = pd.DataFrame({
        "gene": [g.vs[i]["name"] for i in range(g.vcount())],
        "BC_raw": bc_values,
        "Lokalizasyon": [g.vs[i]["lokalizasyon"] for i in range(g.vcount())]
    })
    log.info("  Mouse BC modu: %s", bc_mode)

    bc_thresholds = {}
    for loc in bc_by_gene["Lokalizasyon"].unique():
        loc_bc = bc_by_gene[bc_by_gene["Lokalizasyon"] == loc]["BC_raw"]
        if len(loc_bc) > 0:
            bc_thresholds[loc] = float(np.percentile(loc_bc, WBI_BC_PERCENTILE))
        else:
            bc_thresholds[loc] = bc_threshold_global

    log.info("  Mouse BC %.0f. persentil (global) = %.6f", WBI_BC_PERCENTILE, bc_threshold_global)
    log.info("  Mouse kompartıman-bazlı BC eşikleri: %d lokalizasyon", len(bc_thresholds))
    del g_masked

    name_to_idx  = {v["name"]: v.index for v in g.vs}
    gate_records = []

    for _, row in df.iterrows():
        gene         = row["gene"]
        self_loc_str = str(row.get("Lokalizasyon", "Unknown"))

        base_rec = {
            "gene"                 : gene,
            "Gümrük_Kapisi"        : False,
            "Dis_Komsu_Orani"      : float("nan"),
            "BC_Skoru"             : bc.get(gene, 0.0),
            "Komsu_Lokalizasyonlar": "N/A",
            "BC_Esik_Kullanilan"   : float("nan"),
        }

        if gene not in name_to_idx:
            gate_records.append(base_rec)
            continue

        vidx     = name_to_idx[gene]
        bc_score = bc.get(gene, 0.0)

        if self_loc_str == "Unknown":
            bc_threshold = bc_threshold_global
            if bc_score >= bc_95_threshold and bc_score > 0:
                base_rec["Gümrük_Kapisi"]         = True
                base_rec["BC_Skoru"]              = round(bc_score, 6)
                base_rec["BC_Esik_Kullanilan"]    = round(bc_95_threshold, 6)
                base_rec["Komsu_Lokalizasyonlar"] = "TOPOLOJİK İSTİSNA (UNKNOWN)"
            gate_records.append(base_rec)
            continue

        self_locs = self_loc_str.split(" | ")
        bc_thresholds_for_gene = [
            bc_thresholds.get(loc, bc_threshold_global)
            for loc in self_locs
        ]
        bc_threshold = min(bc_thresholds_for_gene) if bc_thresholds_for_gene else bc_threshold_global

        base_rec["BC_Esik_Kullanilan"] = round(bc_threshold, 6)

        neighbors = g.neighbors(vidx)
        if not neighbors:
            base_rec["Dis_Komsu_Orani"]      = 0.0
            base_rec["Komsu_Lokalizasyonlar"] = "no_neighbors"
            gate_records.append(base_rec)
            continue

        self_locs_set = set(self_locs)
        total_w = ext_w = 0.0
        loc_counter: dict[str, float] = {}

        for nb_idx in neighbors:
            nb_loc_str = str(g.vs[nb_idx]["lokalizasyon"])
            if nb_loc_str == "Unknown":
                continue
            nb_locs = set(nb_loc_str.split(" | "))
            eid     = g.get_eid(vidx, nb_idx, error=False)
            w       = g.es[eid]["weight"] if eid != -1 else 1.0
            total_w += w
            if not self_locs_set.intersection(nb_locs):
                ext_w += w
            for loc in nb_locs:
                loc_counter[loc] = loc_counter.get(loc, 0.0) + w

        if total_w == 0:
            base_rec["Komsu_Lokalizasyonlar"] = "all_unknown"
            gate_records.append(base_rec)
            continue

        ratio    = ext_w / total_w
        is_gate  = (ratio > WBI_EXT_RATIO_MIN) and (bc_score >= bc_threshold)
        summary  = ", ".join(
            f"{l}:{c:.2f}"
            for l, c in sorted(loc_counter.items(), key=lambda x: x[1], reverse=True)[:3]
        )
        gate_records.append({
            "gene"                 : gene,
            "Gümrük_Kapisi"        : is_gate,
            "Dis_Komsu_Orani"      : round(ratio, 4),
            "BC_Skoru"             : round(bc_score, 6),
            "BC_Esik_Kullanilan"   : round(bc_threshold, 6),
            "Komsu_Lokalizasyonlar": summary,
        })

    gate_df = pd.DataFrame(gate_records)
    gate_df["BC_Hesaplama_Modu"] = bc_mode
    gate_df["BC_Ornek_Kaynak_Sayisi"] = source_count
    gate_df["BC_Random_Seed"] = bc_random_seed
    gate_df["BC_Graf_Dugum_Sayisi"] = n_vertices
    gate_df["BC_Graf_Kenar_Sayisi"] = g.ecount()
    gate_df["BC_Weighted"] = True
    gate_df["BC_Distance_Attribute"] = "distance=1/strength"
    gate_df["BC_Timestamp_UTC"] = datetime.now(timezone.utc).isoformat()
    gate_df["BC_Implementation_Version"] = "sophiark-bc-v2"
    merged  = df.merge(gate_df, on="gene", how="left")
    log.info("  → Mouse %d Gümrük Kapısı (hibrit WBI, kompartıman-bazlı).", int(merged["Gümrük_Kapisi"].sum()))
    return merged
