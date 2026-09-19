# ---------------------------------------------------------------------------
# metrics.py
# Hinterland skoru, sınıflandırma, trafik matrisi, gümrük kapısı tespiti,
# threshold sweep, görselleştirme ve bilimsel doğrulama paketi.
# ---------------------------------------------------------------------------
import gc
import math
import random
import sqlite3
from datetime import datetime, timezone
from pathlib import Path

import igraph as ig
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import rankdata, norm as scipy_norm

from .config import (
    NOISE_MAX_DEGREE,
    NOISE_MAX_PR_PERCENTILE,
    NOISE_MAX_SCORE,
    NULL_MODEL_ITER,
    NULL_MODEL_Z_CUTOFF,
    OUT_DIR,
    PAGERANK_DAMPING,
    SWEEP_MAX,
    SWEEP_MIN,
    SWEEP_STEP,
    THRESHOLD_BRIDGE,
    THRESHOLD_EMPEROR,
    THRESHOLD_STRATEGIC,
    WBI_BC_PERCENTILE,
    WBI_EXT_RATIO_MIN,
    _DARK_BG,
    _PANEL_BG,
    log,
)
from .graph_engine import _mask_inf_distances

# ===========================================================================
# ADIM 2c: Hinterland Skoru
# ===========================================================================

def build_hinterland_score(df: pd.DataFrame, pr: dict) -> pd.DataFrame:
    log.info("Adım 2c: Hinterland_Skoru (rank-based, PR:0.60 / deg:0.40)…")
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
        "  → min=%.2f max=%.2f ort=%.2f",
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
        "Adım 3a: Sınıflandırma (İmparator>=%d, Dağıtıcı>=%d)…",
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
        log.info("    %-45s : %7d  (%.2f%%)", cat, cnt, cnt / len(df) * 100)
    return df


# ===========================================================================
# GÖREV 3: Kompartıman Trafik Matrisi
# ===========================================================================

def build_traffic_matrix(g: ig.Graph, out_dir: Path) -> pd.DataFrame:
    log.info("Görev 3: Kompartıman trafik matrisi (vektörel)…")

    loc_arr = np.array(g.vs["lokalizasyon"], dtype=object)
    elist   = np.array(g.get_edgelist(), dtype=np.int32)
    if elist.size == 0:
        log.warning("  Graf boş, trafik matrisi atlandı.")
        return pd.DataFrame()

    src_locs = loc_arr[elist[:, 0]]
    tgt_locs = loc_arr[elist[:, 1]]

    pair_a = np.where(src_locs <= tgt_locs, src_locs, tgt_locs)
    pair_b = np.where(src_locs <= tgt_locs, tgt_locs, src_locs)

    long_df = (
        pd.DataFrame({"Kompartiman_A": pair_a, "Kompartiman_B": pair_b})
        .groupby(["Kompartiman_A", "Kompartiman_B"], sort=False)
        .size()
        .reset_index(name="Kenar_Sayisi")
    )

    intra = int(
        long_df.loc[long_df["Kompartiman_A"] == long_df["Kompartiman_B"], "Kenar_Sayisi"].sum()
    )
    total = int(long_df["Kenar_Sayisi"].sum())
    inter = total - intra
    log.info(
        "  Toplam: %d | Intra: %d (%.1f%%) | Inter: %d (%.1f%%)",
        total, intra, intra / max(total, 1) * 100,
        inter, inter / max(total, 1) * 100,
    )

    comps  = sorted(
        set(long_df["Kompartiman_A"].tolist() + long_df["Kompartiman_B"].tolist())
    )
    matrix = pd.DataFrame(0, index=comps, columns=comps)
    for _, row in long_df.iterrows():
        a, b, cnt = row["Kompartiman_A"], row["Kompartiman_B"], row["Kenar_Sayisi"]
        matrix.loc[a, b] += cnt
        if a != b:
            matrix.loc[b, a] += cnt

    matrix.to_csv(out_dir / "reports" / "kompartiman_trafigi.csv", encoding="utf-8-sig")
    return long_df


# ===========================================================================
# Görselleştirme
# ===========================================================================

def _style_ax(ax, title="", xlabel="", ylabel=""):
    ax.set_facecolor(_PANEL_BG)
    ax.set_title(title, color="white", fontsize=13, pad=10)
    ax.set_xlabel(xlabel, color="white", fontsize=11)
    ax.set_ylabel(ylabel, color="white", fontsize=11)
    ax.tick_params(colors="white", labelsize=9)
    for s in ax.spines.values():
        s.set_edgecolor("#333")


def plot_localization_bar(df: pd.DataFrame, out_dir: Path) -> None:
    log.info("Görsel: Lokalizasyon dağılım grafiği…")
    dist = df.groupby(["Lokalizasyon", "Kategori"]).size().unstack(fill_value=0)
    colors = {
        "Koyu Kırmızı (İmparator Hub)"     : "#8B0000",
        "Açık Kırmızı (Stratejik Dağıtıcı)": "#FF4444",
        "Sarı (Aktif Geçiş Yolu)"          : "#FFD700",
        "Mavi (İzole Sınır Polisi)"         : "#1E90FF",
        "Gri (İstatistiksel Gürültü)"       : "#A9A9A9",
    }
    cols = [c for c in colors if c in dist.columns]
    dist = dist[cols].sort_values(cols[0] if cols else dist.columns[0], ascending=False)
    fig, ax = plt.subplots(figsize=(14, 6), facecolor=_DARK_BG)
    dist.plot(
        kind="bar", stacked=True, ax=ax,
        color=[colors.get(c, "#777") for c in dist.columns],
        edgecolor="#222", linewidth=0.5,
    )
    _style_ax(ax, "Lokalizasyon × Kategori Dağılımı", "Ana Lokalizasyon", "Protein Sayısı")
    ax.legend(loc="upper right", facecolor="#1A1A2E", labelcolor="white",
              fontsize=8, framealpha=0.9)
    ax.set_xticklabels(dist.index.tolist(), rotation=40, ha="right")
    fig.tight_layout()
    fig.savefig(
        out_dir / "reports" / "lokalizasyon_dagilim.png",
        dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor(),
    )
    plt.close(fig)


def plot_threshold_sweep(sweep_df: pd.DataFrame, out_dir: Path) -> None:
    if sweep_df.empty:
        return
    fig, axes = plt.subplots(2, 2, figsize=(14, 9), facecolor=_DARK_BG)
    fig.suptitle("Threshold Sweep — Ağ Direnç Analizi", color="white", fontsize=14, y=0.98)
    metrics = [
        ("LCC_Orani",      "LCC Oranı",           "#00BFFF"),
        ("Ort_Derece",     "Ortalama Derece",      "#FFD700"),
        ("ImparatorHub_N", "İmparator Hub Sayısı", "#FF4444"),
        ("Resilience",     "Resilience Skoru",     "#00FF88"),
    ]
    for ax, (col, label, color) in zip(axes.flat, metrics):
        ax.plot(sweep_df["Esik"], sweep_df[col],
                color=color, linewidth=2.5, marker="o", markersize=7)
        ax.fill_between(sweep_df["Esik"], sweep_df[col], alpha=0.15, color=color)
        _style_ax(ax, title=label, xlabel="Eşik", ylabel=label)
    fig.tight_layout()
    fig.savefig(
        out_dir / "reports" / "threshold_sweep_grafik.png",
        dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor(),
    )
    plt.close(fig)


def plot_dose_response(dose_df: pd.DataFrame, out_dir: Path) -> None:
    if dose_df.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 5), facecolor=_DARK_BG)
    ax.plot(
        dose_df["Survival_Fraction"], dose_df["Signal_Kayip_Pct"],
        color="#FF6B6B", linewidth=2.5, marker="o", markersize=8, label="Sinyal Kaybı %",
    )
    ic50_val = dose_df["IC50_Survival_Frac"].iloc[0] if "IC50_Survival_Frac" in dose_df.columns else None
    try:
        ic50_float = float(ic50_val)
        ax.axvline(
            ic50_float, color="#FFD700", linestyle="--", linewidth=1.5,
            label=f"IC50 (Surv.Frac) ≈ {ic50_float:.3f}",
        )
    except (TypeError, ValueError):
        pass
    ax.axhline(50, color="#888", linestyle=":", linewidth=1)
    ax.invert_xaxis()
    _style_ax(
        ax,
        "Farmakolojik Doz-Yanıt (Hill n=2)",
        "Kalan Aktivite Fraksiyonu (← artan inhibisyon)",
        "Sinyal Kaybı %",
    )
    ax.legend(facecolor="#1A1A2E", labelcolor="white", fontsize=9)
    fig.tight_layout()
    fig.savefig(
        out_dir / "reports" / "doz_yanit_egrisi.png",
        dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor(),
    )
    plt.close(fig)


def plot_congestion(congestion_df: pd.DataFrame, out_dir: Path) -> None:
    if congestion_df.empty:
        return
    top = congestion_df.nlargest(15, "Max_Congestion_Score")
    fig, ax = plt.subplots(figsize=(10, 5), facecolor=_DARK_BG)
    ax.bar(
        top["Survival_Fraction"].astype(str), top["Max_Congestion_Score"],
        color="#FF6B35", edgecolor="#222", linewidth=0.5,
    )
    _style_ax(
        ax,
        "Kalan Aktiviteye Göre Sıkışıklık (BC Artışı)",
        "Kalan Aktivite Fraksiyonu",
        "Max Congestion Score (ΔBC)",
    )
    fig.tight_layout()
    fig.savefig(
        out_dir / "reports" / "congestion_analiz.png",
        dpi=150, bbox_inches="tight", facecolor=fig.get_facecolor(),
    )
    plt.close(fig)


# ===========================================================================
# GÖREV 4: Gümrük Kapısı Tespiti
# ===========================================================================

def detect_customs_gates(
    g: ig.Graph,
    df: pd.DataFrame,
    bc_sample_sources: int | None = None,
) -> pd.DataFrame:
    """Hibrit WBI Gümrük Kapısı Tespiti (kompartıman-bazlı BC persentili)."""
    log.info(
        "Görev 4: Hibrit WBI Gümrük Kapısı "
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
        # Sabit tohum, doku aynı kaldıkça aynı aday sıralamasını üretir.
        bc_random_seed = 42
        sources = np.random.default_rng(bc_random_seed).choice(
            n_vertices, size=source_count, replace=False
        ).tolist()
        log.info(
            "  BC hesaplanıyor — hızlı örneklem: %d/%d kaynak, "
            "weights='distance'…", source_count, n_vertices,
        )
        try:
            bc_values = g_masked.betweenness(
                weights="distance", directed=False, sources=sources,
            )
            # Örneklem toplamını tam ağ ölçeğine getir. Kapı seçimindeki
            # persentiller ölçekten bağımsızdır; bu değer raporlamayı da
            # karşılaştırılabilir tutar.
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
        log.info("  BC hesaplanıyor — tam, weights='distance' + Inf maskeleme…")
        bc_values = g_masked.betweenness(weights="distance", directed=False)
        source_count = n_vertices
        bc_random_seed = 42
        bc_mode = "Full Validation"

    bc_arr       = np.array(bc_values)
    bc           = {g.vs[i]["name"]: bc_values[i] for i in range(g.vcount())}

    # Global BC persentili (yedek olarak)
    bc_threshold_global = float(np.percentile(bc_arr, WBI_BC_PERCENTILE))
    bc_95_threshold = float(np.percentile(bc_arr, 95.0)) if len(bc_arr) > 0 else float("inf")

    # Kompartıman-bazlı BC persentili hesapla
    bc_by_gene = pd.DataFrame({
        "gene": [g.vs[i]["name"] for i in range(g.vcount())],
        "BC_raw": bc_values,
        "Lokalizasyon": [g.vs[i]["lokalizasyon"] for i in range(g.vcount())]
    })
    log.info("  BC modu: %s", bc_mode)

    # Her lokalizasyon için BC persentili
    bc_thresholds = {}
    for loc in bc_by_gene["Lokalizasyon"].unique():
        loc_bc = bc_by_gene[bc_by_gene["Lokalizasyon"] == loc]["BC_raw"]
        if len(loc_bc) > 0:
            bc_thresholds[loc] = float(np.percentile(loc_bc, WBI_BC_PERCENTILE))
        else:
            bc_thresholds[loc] = bc_threshold_global

    log.info("  BC %.0f. persentil (global) = %.6f", WBI_BC_PERCENTILE, bc_threshold_global)
    log.info("  Kompartıman-bazlı BC eşikleri hesaplandı: %d lokalizasyon", len(bc_thresholds))
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

        # Kompartıman-bazlı BC eşiği seç
        if self_loc_str == "Unknown":
            bc_threshold = bc_threshold_global
            if bc_score >= bc_95_threshold and bc_score > 0:
                base_rec["Gümrük_Kapisi"]         = True
                base_rec["BC_Skoru"]              = round(bc_score, 6)
                base_rec["BC_Esik_Kullanilan"]    = round(bc_95_threshold, 6)
                base_rec["Komsu_Lokalizasyonlar"] = "TOPOLOJİK İSTİSNA (UNKNOWN)"
            gate_records.append(base_rec)
            continue

        # Birden fazla lokalizasyon varsa (örn: "Nucleus | Cytoplasm")
        # minimum eşiği kullan (en toleranslı)
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
        # Kompartıman-bazlı BC eşiği kullan
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
    log.info("  → %d Gümrük Kapısı (hibrit WBI, kompartıman-bazlı).", int(merged["Gümrük_Kapisi"].sum()))
    return merged


# ===========================================================================
# GÖREV C: Threshold Sweep
# ===========================================================================

def run_threshold_sweep(
    db_path: str,
    sweep_min: int = SWEEP_MIN,
    sweep_max: int = SWEEP_MAX,
    sweep_step: int = SWEEP_STEP,
    out_dir: Path = OUT_DIR,
    sweep_full_metrics: bool = False,
) -> pd.DataFrame:
    log.info(
        "Görev C: Threshold Sweep (%d→%d, adım=%d)…",
        sweep_min, sweep_max, sweep_step,
    )
    thresholds  = list(range(sweep_min, sweep_max + 1, sweep_step))
    records     = []
    sweep_query = "SELECT protein1, protein2, combined_score FROM interactions WHERE combined_score >= ?"

    con = sqlite3.connect(db_path)
    try:
        for thr in thresholds:
            log.info("  Eşik %d taranıyor…", thr)
            names_set: dict[str, int] = {}
            edges_tmp: list[tuple]    = []
            weights_tmp: list[float]  = []

            def _goc(nm: str) -> int:
                if nm not in names_set:
                    names_set[nm] = len(names_set)
                return names_set[nm]

            cursor = con.execute(sweep_query, (thr,))
            while True:
                rows = cursor.fetchmany(200_000)
                if not rows:
                    break
                for r in rows:
                    ui = _goc(r[0]); vi = _goc(r[1])
                    edges_tmp.append((ui, vi))
                    weights_tmp.append(float(r[2] / 1000.0) ** 2)

            n_nodes = len(names_set)
            if n_nodes == 0:
                del edges_tmp, weights_tmp, names_set
                gc.collect()
                continue

            G_temp   = ig.Graph(n=n_nodes, edges=edges_tmp, directed=False)
            G_temp.es["weight"] = weights_tmp
            G_temp   = G_temp.simplify(combine_edges={"weight": "max"})
            n_nodes  = G_temp.vcount()
            n_edges  = G_temp.ecount()
            clusters = G_temp.clusters()
            lcc_size  = max(len(c) for c in clusters) if clusters else 0
            lcc_ratio = lcc_size / n_nodes
            avg_deg   = 2 * n_edges / n_nodes
            if n_nodes > 1 and avg_deg > 0:
                log_n = math.log(n_nodes + 1)
                if log_n > 0:
                    resilience = lcc_ratio * math.sqrt(avg_deg / log_n)
                else:
                    resilience = 0.0
            else:
                resilience = 0.0

            emperor_count = 0
            if sweep_full_metrics:
                pr_vals = G_temp.pagerank(
                    damping=PAGERANK_DAMPING, weights="weight", directed=False
                )
                pr_v = np.array(pr_vals)
                n    = len(pr_v)
                pr_n = rankdata(pr_v, method="average") / n
                dg_v = np.array(G_temp.degree())
                dg_n = rankdata(np.log10(dg_v + 1), method="average") / n
                scores        = (0.60 * pr_n + 0.40 * dg_n) * 100
                emperor_count = int(np.sum(scores >= THRESHOLD_EMPEROR))
                del pr_vals, pr_v, scores

            records.append({
                "Esik"          : thr,
                "Dugum_Sayisi"  : n_nodes,
                "Kenar_Sayisi"  : n_edges,
                "LCC_Orani"     : round(lcc_ratio, 4),
                "Ort_Derece"    : round(avg_deg, 4),
                "ImparatorHub_N": emperor_count,
                "Resilience"    : round(resilience, 6),
            })
            del G_temp, edges_tmp, weights_tmp, names_set
            gc.collect()
    finally:
        con.close()

    sweep_df = pd.DataFrame(records)
    sweep_df.to_csv(out_dir / "reports" / "threshold_sweep_raporu.csv", index=False, encoding="utf-8-sig")
    log.info("  → Sweep raporu kaydedildi.")
    return sweep_df


# ===========================================================================
# BİLİMSEL DOĞRULAMA PAKETİ
# ===========================================================================

def validate_topological_importance(
    g: ig.Graph,
    df: pd.DataFrame,
    top_n: int = 10,
) -> tuple:
    log.info("─" * 60)
    log.info("[VAL-1] BC Hub Knockout Test (top_n=%d)…", top_n)

    if g.vcount() == 0:
        log.warning("[VAL-1] Boş graf, test atlandı.")
        return pd.DataFrame(), 0.0

    try:
        clusters = g.clusters()
        orig_lcc = max(len(c) for c in clusters)
    except Exception:
        log.warning("[VAL-1] LCC hesaplanamadı, test atlandı.")
        return pd.DataFrame(), 0.0

    log.info("  Başlangıç LCC: %d düğüm", orig_lcc)

    if "BC_Skoru" not in df.columns:
        log.warning("  ⚠ BC_Skoru sütunu yok — detect_customs_gates önce çalıştırılmalı.")
        return pd.DataFrame(), 0.0

    top_genes = df.nlargest(top_n, "BC_Skoru")["gene"].astype(str).tolist()
    all_names = set(g.vs["name"])
    results   = []

    for gene in top_genes:
        if gene not in all_names:
            log.warning("  ⚠ %s grafta bulunamadı, atlandı.", gene)
            continue
        temp_g = g.copy()
        try:
            idx     = temp_g.vs.find(name=gene).index
            temp_g.delete_vertices(idx)
            if temp_g.vcount() == 0:
                continue
            comps   = temp_g.clusters()
            new_lcc = max(len(c) for c in comps) if comps else 0
            drop    = (orig_lcc - new_lcc) / orig_lcc * 100
            results.append({"gene": gene, "LCC_Drop_Pct": round(drop, 6)})
            log.info("    ✂ %-30s → LCC düşüş: %.4f%%", gene, drop)
        except Exception as e:
            log.warning("  ⚠ Knockout hata (%s): %s", gene, e)
        finally:
            del temp_g
            gc.collect()

    rand_drop = 0.0
    try:
        rand_names    = random.sample(list(all_names), min(top_n, g.vcount()))
        rand_g        = g.copy()
        current_names = set(rand_g.vs["name"])
        rand_indices  = [rand_g.vs.find(name=nm).index for nm in rand_names if nm in current_names]
        rand_indices_sorted = sorted(rand_indices, reverse=True)
        rand_g.delete_vertices(rand_indices_sorted)
        comps_r   = rand_g.clusters()
        rand_lcc  = max(len(c) for c in comps_r) if comps_r else 0
        rand_drop = (orig_lcc - rand_lcc) / orig_lcc * 100
        del rand_g
        gc.collect()
    except Exception as e:
        log.warning("  ⚠ Rastgele kontrol hata: %s", e)

    hub_mean = float(np.mean([r["LCC_Drop_Pct"] for r in results])) if results else 0.0
    ratio    = hub_mean / rand_drop if rand_drop > 0 else float("inf")

    log.info("  → BC Hub knockout ort. düşüş : %.4f%%", hub_mean)
    log.info("  → Rastgele knockout düşüş    : %.4f%%", rand_drop)
    log.info("  → BC Hub / Rastgele oranı    : %.2fx  %s",
             ratio, "✅ GEÇER (>=5x)" if ratio >= 5 else "⚠ DİKKAT (<5x)")

    return pd.DataFrame(results), rand_drop


def validate_jaccard_stability(
    g: ig.Graph,
    df: pd.DataFrame,
    rewire_pct: float = 0.05,
) -> float:
    log.info("─" * 60)
    log.info("[VAL-2] Jaccard Stability Test (rewire=%.0f%%)…", rewire_pct * 100)

    try:
        orig_top = set(
            df[df["Hinterland_Skoru"] >= THRESHOLD_EMPEROR]["gene"].astype(str)
        )
        log.info("  Orijinal Emperor Hub: %d gen", len(orig_top))
        if not orig_top:
            log.warning("[VAL-2] Emperor Hub listesi boş.")
            return 0.0

        noisy_g = g.copy()
        n_swap  = min(int(g.ecount() * rewire_pct), 50_000)
        if n_swap <= 10:
            log.info("  [VAL-2] Swap sayısı çok düşük (%d), stabilite=1.0 varsayıldı.", n_swap)
            return 1.0

        log.info("  %d kenar rewire ediliyor…", n_swap)
        try:
            noisy_g.rewire(n=n_swap, mode="simple")
        except Exception as e:
            log.warning("  ⚠ Rewire kısmi başarısız: %s", e)

        # D6 düzeltmesi: orijinal formülle simetrik skor (PR:0.60 + deg:0.40)
        pr_noisy_vals = noisy_g.pagerank(
            damping=PAGERANK_DAMPING, weights="weight", directed=False
        )
        noisy_n = noisy_g.vcount()
        pr_norm = rankdata(pr_noisy_vals, method="average") / noisy_n
        deg_arr = np.array(noisy_g.degree(), dtype=float)
        dg_norm = rankdata(np.log10(deg_arr + 1), method="average") / noisy_n
        noisy_scores = (0.60 * pr_norm + 0.40 * dg_norm) * 100.0

        noisy_df = pd.DataFrame({
            "gene" : noisy_g.vs["name"],
            "score": noisy_scores,
            "degree": deg_arr,
            "pr_raw": pr_noisy_vals
        })

        # PR persentilini hesapla (0-1 arası)
        noisy_df["pr_percentile"] = noisy_df["pr_raw"].rank(pct=True)

        # Emperor Hub adaylarını al
        emperor_candidates = noisy_df[noisy_df["score"] >= THRESHOLD_EMPEROR]

        # Gürültü filtresi: Yüksek dereceli ama düşük PR persentilli genleri dışarıda bırak
        filtered_emperors = emperor_candidates[
            (emperor_candidates["degree"] < NOISE_MAX_DEGREE) |
            (emperor_candidates["pr_percentile"] > NOISE_MAX_PR_PERCENTILE / 100.0)
        ]

        noisy_top = set(filtered_emperors["gene"].astype(str))

        inter   = len(orig_top & noisy_top)
        union   = len(orig_top | noisy_top)
        j_index = inter / union if union > 0 else 0.0

        log.info("  Kesişim=%d | Birleşim=%d | Jaccard=%.4f  %s",
                 inter, union, j_index,
                 "✅ GEÇER" if j_index >= 0.70 else "⚠ KABUL" if j_index >= 0.50 else "❌ BAŞARISIZ")

        del noisy_g, noisy_df
        gc.collect()
        return float(j_index)

    except Exception as e:
        log.warning("[VAL-2] Jaccard testi hata verdi: %s", e)
        return 0.0


def validate_biological_coherence(df: pd.DataFrame, top_n: int = 50) -> pd.DataFrame:
    log.info("─" * 60)
    log.info("[VAL-3] Functional Enrichment — KEGG (top_n=%d)…", top_n)
    try:
        from .services.enrichment import run_post_simulation_enrichment
        mapping_path = BASE_DIR / "data" / "processed" / "ensp_with_symbols.csv"
        mapping = pd.read_csv(mapping_path, usecols=["gene", "Symbol"]).dropna().drop_duplicates("gene")
        symbol_by_gene = dict(zip(mapping["gene"].astype(str), mapping["Symbol"].astype(str)))
        active_genes = df["gene"].astype(str).tolist()
        top_genes = df.nlargest(top_n, "Hinterland_Skoru")["gene"].astype(str).tolist()
        symbols = [symbol_by_gene[gene] for gene in top_genes if gene in symbol_by_gene]
        universe = [symbol_by_gene[gene] for gene in active_genes if gene in symbol_by_gene]
        enrichment, notices = run_post_simulation_enrichment(
            symbols, is_mouse=False, background_symbols=universe,
        )
        for notice in notices:
            log.info("  [VAL-3] %s", notice)
        if enrichment.empty:
            log.warning("[VAL-3] Aktif ağ evrenli enrichment sonucu boş.")
            return pd.DataFrame()

        top_pathways = enrichment.nsmallest(5, "Adjusted P-value")
        sig_count = int((enrichment["Adjusted P-value"] < 0.05).sum())
        log.info("  p-adj < 0.05 olan yolak sayısı: %d", sig_count)
        for _, row in top_pathways.iterrows():
            log.info("    → %-55s (p-adj: %.2e)", row["Term"], row["Adjusted P-value"])
        log.info("  → [VAL-3] %s",
                 "✅ GEÇER" if sig_count >= 3 else "⚠ ZAYIF" if sig_count >= 1 else "❌ BAŞARISIZ")
        return top_pathways
    except Exception as e:
        log.warning("[VAL-3] Enrichment hata: %s", e)
        return pd.DataFrame()


# ===========================================================================
# validate_null_model — O5 düzeltmesi: yalnızca topoloji rewire
# ===========================================================================

def validate_null_model(
    g: ig.Graph,
    df: pd.DataFrame,
    n_iter: int = NULL_MODEL_ITER,
    z_cutoff: float = NULL_MODEL_Z_CUTOFF,
    out_dir: Path = OUT_DIR,
) -> pd.DataFrame:
    """
    Null Model Z-Score Doğrulaması.

    O5 düzeltmesi: v4.0'da hem topoloji rewire hem ağırlık shuffle yapılıyordu;
    bu çift rastgeleleştirme null dağılım varyansını şişiriyor ve Z skorlarını
    gereğinden küçük üretiyordu. v5.0'da yalnızca topoloji rewire uygulanır
    (ağırlıklar sabit kalır) — null hipotezi tutarlı: "bu BC değeri rastgele
    bir topolojide de elde edilir mi?"
    """
    log.info("─" * 60)
    log.info(
        "[VAL-NULL] Null Model Z-Score (%d iterasyon, Z-kesme=%.1f)…",
        n_iter, z_cutoff,
    )

    if "Gümrük_Kapisi" not in df.columns or "BC_Skoru" not in df.columns:
        log.warning(
            "  ⚠ 'Gümrük_Kapisi' veya 'BC_Skoru' sütunu bulunamadı. "
            "detect_customs_gates() önce çalıştırılmalı."
        )
        return pd.DataFrame()

    gate_df = df[df["Gümrük_Kapisi"] == True].copy()
    if gate_df.empty:
        log.warning("  ⚠ Hiç Gümrük Kapısı yok. Null Model atlandı.")
        return pd.DataFrame()

    gate_genes  = gate_df["gene"].astype(str).tolist()
    real_bc     = dict(zip(gate_genes, gate_df["BC_Skoru"].values))
    null_dists: dict[str, list] = {gn: [] for gn in gate_genes}

    n_swap = max(int(g.ecount() * 0.10), 10)
    log.info("  Swap/iter=%d | Gümrük Kapısı=%d", n_swap, len(gate_genes))

    for it in range(n_iter):
        if (it + 1) % 5 == 0:
            log.info("    Null iter %d/%d…", it + 1, n_iter)

        null_g = g.copy()

        try:
            null_g.rewire(n=n_swap, mode="simple")
        except Exception as e:
            log.warning("    ⚠ Rewire hatası iter=%d: %s", it + 1, e)

        # O5: ağırlıklar değiştirilmiyor — sadece topoloji rastgele
        _mask_inf_distances(null_g)

        null_bc_vals = null_g.betweenness(weights="distance", directed=False)
        null_bc      = {null_g.vs[i]["name"]: null_bc_vals[i] for i in range(null_g.vcount())}

        for gn in gate_genes:
            null_dists[gn].append(null_bc.get(gn, 0.0))

        del null_g
        gc.collect()

    results = []
    for gn in gate_genes:
        arr     = np.array(null_dists[gn])
        mu      = float(np.mean(arr))
        sigma   = float(np.std(arr))
        bc_real = real_bc.get(gn, 0.0)
        z       = (bc_real - mu) / sigma if sigma > 0 else float("inf")
        p_val   = float(1.0 - scipy_norm.cdf(z))
        is_sig  = z > z_cutoff
        results.append({
            "gene"        : gn,
            "BC_Gercek"   : round(bc_real, 6),
            "BC_Null_Ort" : round(mu, 6),
            "BC_Null_Std" : round(sigma, 6),
            "Z_Skoru"     : round(z, 4),
            "P_Deger"     : round(p_val, 6),
            "Anlamli_Z2"  : is_sig,
        })

    result_df = pd.DataFrame(results).sort_values("Z_Skoru", ascending=False)
    sig_count = int(result_df["Anlamli_Z2"].sum())
    log.info(
        "  → Toplam: %d | Z>%.1f anlamlı: %d (%.1f%%)",
        len(gate_genes), z_cutoff,
        sig_count, sig_count / max(len(gate_genes), 1) * 100,
    )

    if sig_count == 0:
        log.warning(
            "  ⚠ Hiç anlamlı Gümrük Kapısı bulunamadı (Z<%.1f). "
            "WBI_EXT_RATIO_MIN veya WBI_BC_PERCENTILE'i sorgula.",
            z_cutoff,
        )
    else:
        for _, row in result_df.head(5).iterrows():
            log.info(
                "    ★ %-30s Z=%.2f  p=%.4f  %s",
                row["gene"], row["Z_Skoru"], row["P_Deger"],
                "✅" if row["Anlamli_Z2"] else "—",
            )

    out_path = out_dir / "reports" / "null_model_z_raporu.csv"
    result_df.to_csv(out_path, index=False, encoding="utf-8-sig")
    log.info("  → Null Model raporu: %s", out_path)
    return result_df
def compute_druggable_gate_score(df: pd.DataFrame, target_type: str = "all") -> pd.DataFrame:
    """
    Proteinler için İlaçlanabilirlik Skoru (Drug_Score) hesaplar.
    
    target_type: "all" -> tüm proteinlere hesapla
                 "gates" -> sadece Gümrük Kapılarına hesapla
                 "stressed" -> sadece stresli genlere hesapla
    """
    from .data_loader import query_structural_db
    import numpy as np

    # Hangi proteinler için hesaplanacak?
    if target_type == "gates":
        targets = df[df["Gümrük_Kapisi"] == True].copy()
    elif target_type == "stressed":
        targets = df[df["Hasar_Tipi"] == "Üçüncül_Hasar_Stresli"].copy()
    else:  # "all"
        targets = df.copy()

    if targets.empty:
        df["Drug_Score"] = np.nan
        return df

    ensp_list = targets["gene"].tolist()
    struct = query_structural_db(ensp_list)

    # Normalizasyon aralıkları
    NDS_MIN, NDS_MAX = 0.20, 0.81
    PLDDT_MIN, PLDDT_MAX = 26.0, 99.0
    SASA_MIN, SASA_MAX = 1700.0, 330000.0
    COMP_MIN, COMP_MAX = 0.0, 1.0
    HYDROPHOBIC_RATIO_MIN, HYDROPHOBIC_RATIO_MAX = 0.0, 1.0
    CHARGED_RATIO_MIN, CHARGED_RATIO_MAX = 0.0, 1.0

    scores = {}
    for ensp in ensp_list:
        info = struct.get(ensp, {})

        nds = info.get("nds")
        plddt = info.get("mean_plddt")
        sasa = info.get("sasa")
        num_lig = info.get("binding_num_ligands")
        compactness = info.get("compactness")
        hydrophobic_ratio = info.get("hydrophobic_ratio")
        charged_ratio = info.get("charged_ratio")

        nds_norm = max(0, min(1, (nds - NDS_MIN) / (NDS_MAX - NDS_MIN))) if nds is not None else 0
        plddt_norm = max(0, min(1, (plddt - PLDDT_MIN) / (PLDDT_MAX - PLDDT_MIN))) if plddt is not None else 0
        sasa_norm = max(0, min(1, (sasa - SASA_MIN) / (SASA_MAX - SASA_MIN))) if sasa is not None else 0
        comp_norm = max(0, min(1, (compactness - COMP_MIN) / (COMP_MAX - COMP_MIN))) if compactness is not None else 0
        hydro_norm = max(0, min(1, (hydrophobic_ratio - HYDROPHOBIC_RATIO_MIN) / (HYDROPHOBIC_RATIO_MAX - HYDROPHOBIC_RATIO_MIN))) if hydrophobic_ratio is not None else 0
        charged_norm = max(0, min(1, (charged_ratio - CHARGED_RATIO_MIN) / (CHARGED_RATIO_MAX - CHARGED_RATIO_MIN))) if charged_ratio is not None else 0

        if num_lig and num_lig > 0:
            ligand_bonus = min(1.0, np.log10(num_lig + 1) / np.log10(1000))
        else:
            ligand_bonus = 0

        score = (
            0.25 * nds_norm +
            0.20 * plddt_norm +
            0.15 * sasa_norm +
            0.20 * ligand_bonus +
            0.10 * comp_norm +
            0.05 * hydro_norm +
            0.05 * charged_norm
        )
        scores[ensp] = score * 100

    df["Drug_Score"] = df["gene"].map(scores)
    return df


# ===========================================================================
# metrics.py içine, validate_null_model() fonksiyonunun ALTINA ekle.
# Bu fonksiyon, config.py'deki GENE_ESSENTIALITY sözlüğünü (BioGRID ORCS)
# kullanarak Hinterland_Skoru'nun gerçek biyolojik essentiality ile
# korelasyonunu test eder — validate_against_known_biology'nin (sentinel
# test) istatistiksel/ölçeklenebilir versiyonu.
# ===========================================================================

def validate_essentiality_correlation(
    df: pd.DataFrame,
    symbol_col: str = "Symbol",
    min_test_sayisi: int = 3,
    out_dir: Path = OUT_DIR,
) -> dict:
    """
    Hinterland_Skoru ile BioGRID ORCS essentiality (Hit_Orani) arasındaki
    Spearman korelasyonunu hesaplar.

    Mantık: eğer ağ topolojisi (PageRank + derece) gerçekten biyolojik
    önemi yakalıyorsa, yüksek Hinterland_Skoru'na sahip genlerin gerçek
    CRISPR ekranlarında da daha sık "essential/hit" çıkması beklenir.
    Bu, sentinel testten farklı olarak elle seçilmiş birkaç gen yerine
    TÜM veri setinde istatistiksel bir dış doğrulama sağlar.

    NOT: BioGRID ORCS MIT lisanslıdır, ticari kullanım dahil serbesttir.
    """
    from scipy.stats import spearmanr
    from .config import GENE_ESSENTIALITY, ESSENTIALITY_ENABLED

    log.info("─" * 60)
    log.info("[VAL-ESSENTIALITY] Hinterland_Skoru vs BioGRID ORCS essentiality…")

    if not ESSENTIALITY_ENABLED or not GENE_ESSENTIALITY:
        log.warning("  ⚠ biogrid_essentiality.pkl yüklenmemiş. Test atlandı.")
        return {"durum": "VERI_YOK"}

    if symbol_col not in df.columns or "Hinterland_Skoru" not in df.columns:
        log.warning(f"  ⚠ '{symbol_col}' veya 'Hinterland_Skoru' sütunu yok. Test atlandı.")
        return {"durum": "SUTUN_YOK"}

    calisma_df = df[["gene", symbol_col, "Hinterland_Skoru"]].copy()
    calisma_df["Hit_Orani"] = calisma_df[symbol_col].map(
        lambda s: GENE_ESSENTIALITY.get(s, {}).get("hit_orani")
        if GENE_ESSENTIALITY.get(s, {}).get("test_sayisi", 0) >= min_test_sayisi
        else None
    )
    gecerli = calisma_df.dropna(subset=["Hit_Orani"])

    log.info(
        "  Essentiality verisi olan gen: %d / %d (min_test_sayisi>=%d)",
        len(gecerli), len(calisma_df), min_test_sayisi,
    )

    if len(gecerli) < 20:
        log.warning("  ⚠ Yeterli örtüşen gen yok (<20). Test güvenilir değil, atlandı.")
        return {"durum": "YETERSIZ_ORTUSME", "n": len(gecerli)}

    rho, p_val = spearmanr(gecerli["Hinterland_Skoru"], gecerli["Hit_Orani"])

    if p_val < 0.05 and rho > 0.2:
        sonuc_str = "✅ GEÇER — Hinterland_Skoru gerçek essentiality ile anlamlı pozitif ilişkili"
    elif p_val < 0.05 and rho <= 0.2:
        sonuc_str = "⚠ ZAYIF — istatistiksel anlamlı ama korelasyon zayıf"
    else:
        sonuc_str = "❌ ANLAMSIZ — istatistiksel olarak anlamlı ilişki bulunamadı"

    log.info("  n=%d | Spearman rho=%.4f | p=%.4e", len(gecerli), rho, p_val)
    log.info("  → %s", sonuc_str)

    sonuc = {
        "durum": "TAMAMLANDI",
        "n": len(gecerli),
        "spearman_rho": round(float(rho), 4),
        "p_degeri": float(p_val),
        "yorum": sonuc_str,
    }

    gecerli.to_csv(
        out_dir / "reports" / "essentiality_korelasyon_raporu.csv",
        index=False, encoding="utf-8-sig",
    )
    log.info("  → Rapor: %s", out_dir / "reports" / "essentiality_korelasyon_raporu.csv")

    return sonuc
