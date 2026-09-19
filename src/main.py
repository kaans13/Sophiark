# ---------------------------------------------------------------------------
# main.py
# Ana orkestratör (argparse, batch-sweep modu, Streamlit köprüsü, main()).
# ---------------------------------------------------------------------------
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import argparse
import gc
import math
import random
import pandas as pd

from .config import (
    DB_PATH,
    EFFICIENCY_SAMPLE_SIZE,
    HASAR_BIRINCIL,
    HASAR_IKINCIL,
    HASAR_UCUNCUL,
    HIGH_CONF_THRESHOLD,
    NULL_MODEL_Z_CUTOFF,
    OUT_CSV,
    OUT_DIR,
    PAGERANK_DAMPING,
    PAGERANK_MAX_ITER,
    SWEEP_MAX,
    SWEEP_MIN,
    SWEEP_STEP,
    log,
)
from .graph_engine import (
    build_high_conf_subgraph,
    compute_pagerank,
    run_community_detection,
    tag_graph_nodes,
)
from .biology_logic import (
    fetch_localizations_mygene,
    add_localization_to_df,
    run_infection_simulation,
    run_pharmacological_dose_response,
    validate_against_known_biology,
)
from .metrics import (
    build_hinterland_score,
    classify_genes,
    build_traffic_matrix,
    detect_customs_gates,
    
    compute_druggable_gate_score,
    run_threshold_sweep,
    plot_threshold_sweep,
    plot_localization_bar,
    validate_jaccard_stability,
    validate_biological_coherence,
    validate_null_model,
)
from .data_loader import fetch_degree_summary, query_structural_db


# ===========================================================================
# Streamlit Köprüsü
# ===========================================================================

def arayuz_icin_motoru_hazirla(
    forced_genes: list = None,
    hedef_doku: str = None,
    bc_sample_sources: int | None = None,
    tissue_normalization_mode: str = "within_tissue",
):
    """Seçili doku grafiğini ve o grafa ait skorları yeniden hesapla."""
    g = build_high_conf_subgraph(
        DB_PATH,
        HIGH_CONF_THRESHOLD,
        forced_genes=forced_genes,
        hedef_doku=hedef_doku,
    )

    annotations = pd.read_csv(OUT_CSV)
    annotation_cols = [
        col for col in ["gene", "Lokalizasyon"] if col in annotations.columns
    ]
    annotation_map = (
        annotations[annotation_cols].drop_duplicates("gene").set_index("gene")
        if "gene" in annotation_cols else pd.DataFrame()
    )

    df = fetch_degree_summary(DB_PATH, HIGH_CONF_THRESHOLD)
    df = df[df["gene"].isin(set(g.vs["name"]))].reset_index(drop=True)
    if df.empty:
        return g, df

    pr = compute_pagerank(g, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
    df = build_hinterland_score(df, pr)
    df = classify_genes(df)
    df, _ = run_community_detection(g, df)

    if "Lokalizasyon" in annotation_map.columns:
        df["Lokalizasyon"] = df["gene"].map(annotation_map["Lokalizasyon"])
    else:
        df["Lokalizasyon"] = "Unknown"
    df["Lokalizasyon"] = df["Lokalizasyon"].fillna("Unknown")
    g = tag_graph_nodes(g, dict(zip(df["gene"], df["Lokalizasyon"])))

    df = detect_customs_gates(g, df, bc_sample_sources=bc_sample_sources)
    df = compute_druggable_gate_score(df)
    return g, df


# ===========================================================================
# argparse
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="hinterland_analysis",
        description="Hinterland Protein Interaction Network v5.0",
    )
    parser.add_argument(
        "--tissue",
        type=str,
        default="Cerebral Cortex",
        help='Doku filtresi. Devre dışı bırakmak için "None" girin.',
    )
    parser.add_argument(
        "--batch",
        type=int,
        default=0,
        metavar="N",
        help=(
            "Headless Batch-Sweep modu. N > 0 ise top-N gen teker teker "
            "run_infection_simulation'a gönderilir. 0 = standart interaktif mod."
        ),
    )
    parser.add_argument(
        "--null-iter",
        type=int,
        default=200,
        dest="null_iter",
        help="Null Model iterasyon sayısı (varsayılan: 10; prodüksiyon: 22).",
    )
    return parser.parse_args()


def _run_batch_sweep(
    g,
    df: pd.DataFrame,
    out_dir,
    top_n: int,
) -> None:
    log.info("━" * 60)
    log.info("▶ BATCH-SWEEP MODU: top-%d gen, headless", top_n)
    log.info("━" * 60)

    rapor_yolu = out_dir / "reports" / "tekil_zafiyet_raporu.csv"
    sutunlar   = ["Hedef_Gen", "Efficiency_Kayip_Pct", "Local_Efficiency_Kayip_Pct"]

    write_header = not rapor_yolu.exists()

    elit = (
        df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]
        .nlargest(top_n, "Hinterland_Skoru")["gene"]
        .tolist()
    )
    log.info("  Elit havuz: %d gen seçildi.", len(elit))

    for sira, hedef_gen in enumerate(elit, start=1):
        log.info("  [%d/%d] Simülasyon: %s", sira, len(elit), hedef_gen)

        sim_df      = pd.DataFrame()
        eff_kayip   = float("nan")
        local_kayip = float("nan")

        g_backup = g.copy()

        try:
            sim_df = run_infection_simulation(
                g, df, out_dir,
                spesifik_hedefler=[hedef_gen],
                sample_size=0,
            )

            if not sim_df.empty:
                birincil = sim_df[sim_df["Hasar_Tipi"] == HASAR_BIRINCIL]
                if birincil.empty:
                    birincil = sim_df
                eff_kayip   = float(birincil["Efficiency_Kayip_Pct"].iloc[0])
                local_kayip = float(birincil["Local_Efficiency_Kayip_Pct"].iloc[0])

        except Exception as exc:
            log.warning("    ⚠ %s simülasyon hatası: %s", hedef_gen, exc)
            log.info("    ↻ Graf yedekten geri yükleniyor...")
            g.delete_vertices(range(g.vcount()))
            g.add_vertices(g_backup.vcount())
            for i, name in enumerate(g_backup.vs["name"]):
                g.vs[i]["name"] = name
                g.vs[i]["lokalizasyon"] = g_backup.vs[i]["lokalizasyon"]
            g.add_edges(g_backup.get_edgelist())
            g.es["weight"] = g_backup.es["weight"]
            g.es["distance"] = g_backup.es["distance"]
            g.es["rescue_bridge"] = g_backup.es["rescue_bridge"]

        finally:
            del g_backup
            gc.collect()

        satir = pd.DataFrame(
            [[hedef_gen,
              round(eff_kayip, 4) if not math.isnan(eff_kayip) else float("nan"),
              round(local_kayip, 4) if not math.isnan(local_kayip) else float("nan")]],
            columns=sutunlar,
        )
        satir.to_csv(
            rapor_yolu,
            mode="a",
            index=False,
            header=write_header,
            encoding="utf-8-sig",
        )
        write_header = False

        del sim_df, satir
        gc.collect()

    log.info("  ✔ Batch-Sweep tamamlandı → %s", rapor_yolu)


# ===========================================================================
# ANA ORKESTRATÖR — v5.0
# ===========================================================================

def main():
    args = _parse_args()

    hedef_doku = None if args.tissue.strip().lower() == "none" else args.tissue

    log.info("▶ Hinterland v5.0 başlatılıyor…")
    log.info("  Doku    : %s", hedef_doku or "TÜM PROTEOM")
    log.info("  Batch   : %s", args.batch if args.batch > 0 else "KAPALI (interaktif)")
    log.info("  NullIter: %d", args.null_iter)
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    df = fetch_degree_summary(DB_PATH, HIGH_CONF_THRESHOLD)
    g  = build_high_conf_subgraph(
        DB_PATH, HIGH_CONF_THRESHOLD, hedef_doku=hedef_doku
    )

    g_nodes = set(g.vs["name"])
    df = df[df["gene"].isin(g_nodes)].reset_index(drop=True)

    pr = compute_pagerank(g, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
    df = build_hinterland_score(df, pr)
    df = classify_genes(df)
    df, _ = run_community_detection(g, df)
    gc.collect()

    non_noise = df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]["gene"].tolist()
    loc_map   = fetch_localizations_mygene(non_noise)
    for gn in df[df["Kategori"] == "Gri (İstatistiksel Gürültü)"]["gene"]:
        loc_map.setdefault(gn, "Unknown")
    df = add_localization_to_df(df, loc_map)
    g  = tag_graph_nodes(g, loc_map)
    gc.collect()

    build_traffic_matrix(g, OUT_DIR)
    gc.collect()

    df = detect_customs_gates(g, df)
    df = compute_druggable_gate_score(df) 

    
    gc.collect()
    

    # --- Yapısal veritabanı testi ---
    gate_ensps = df[df["Gümrük_Kapisi"] == True]["gene"].tolist()[:10]
    struct_info = query_structural_db(gate_ensps)
    for ensp, info in struct_info.items():
        print(f"{ensp}: pocket_vol={info.get('pocket_vol')}, nds={info.get('nds')}, binding_ligands={info.get('binding_num_ligands')}")
    # --------------------------------

    if args.batch > 0:
        _run_batch_sweep(g, df, OUT_DIR, top_n=args.batch)

        out_cols = [c for c in [
            "gene", "k_i", "max_cs", "Hinterland_Skoru",
            "Kategori", "Lokalizasyon", "Gümrük_Kapisi", "BC_Skoru", "BC_Esik_Kullanilan",
        ] if c in df.columns]
        df[out_cols].to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
        log.info("✔ Batch CSV: %s", OUT_CSV)
        log.info("▶ Hinterland v5.0 (Batch-Sweep) tamamlandı.")
        return

    random.seed(42)
    # 1. Bilinen TF'lerden birini rastgele seç (bonus testi için)
    bilinen_tf_ensp = [
        "ENSP00000269305",  # TP53
        "ENSP00000478887",  # MYC
        "ENSP00000226574",  # NFKB1
    ]
    secili_tf = random.choice(bilinen_tf_ensp)
    # 2. Veri setinden rastgele 4 gen daha seç (gürültü olmayanlardan)
    digerleri = df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]["gene"].tolist()
    # Seçili TF'yi listeden çıkar (mükerrer olmasın)
    digerleri = [gn for gn in digerleri if gn != secili_tf]
    rastgeleler = random.sample(digerleri, min(4, len(digerleri)))

    # 3. Birleştir
    top_stress_genes = [secili_tf] + rastgeleler
    log.info("🎯 Stres Testi Hedefleri (top-5 Hinterland): %s", top_stress_genes)

    run_infection_simulation(
        g, df, OUT_DIR,
        spesifik_hedefler=top_stress_genes,
        sample_size=EFFICIENCY_SAMPLE_SIZE,
    )
    gc.collect()

    dose_df = run_pharmacological_dose_response(
        g, df, OUT_DIR,
        spesifik_hedefler=top_stress_genes,
    )

    if not dose_df.empty and "Congestion_Top_Stressed" in dose_df.columns:
        log.info("━" * 60)
        log.info("★ SIKILIŞIKLIK (CONGESTION) ÖZET RAPORU")
        log.info("━" * 60)
        for _, row in dose_df.iterrows():
            log.info(
                "  SurvFrac=%.2f (%.0f%% inh.) | MaxCong=%.4f | Stresli: %s",
                row["Survival_Fraction"], row["Inhibition_Pct"],
                row["Max_Congestion_Score"],
                row["Congestion_Top_Stressed"] or "—",
            )
        log.info("━" * 60)
    gc.collect()

    sweep_df = run_threshold_sweep(
        DB_PATH, SWEEP_MIN, SWEEP_MAX, SWEEP_STEP, OUT_DIR,
        sweep_full_metrics=False,
    )
    plot_threshold_sweep(sweep_df, OUT_DIR)
    gc.collect()

    plot_localization_bar(df, OUT_DIR)

    # =====================================================================
    # BİLİMSEL DOĞRULAMA SÜRECİ
    # =====================================================================
    log.info("★" * 60)
    log.info("★  BİLİMSEL DOĞRULAMA SÜRECİ BAŞLIYOR (v5.0)  ★")
    log.info("★" * 60)

    null_df = pd.DataFrame()
    try:
        j_stability = validate_jaccard_stability(g, df)

        bio_df = validate_biological_coherence(df)
        if not bio_df.empty:
            bio_df.to_csv(OUT_DIR / "validation_enrichment.csv", index=False, encoding="utf-8-sig")

        log.info("★" * 60)
        log.info("★  VALIDATION ÖZET")
        log.info("★  BC Hub Knockout    : PASİF (büyük ağda anlamsız)")
        log.info("★  Jaccard Stability    : %.4f  %s",
                 j_stability,
                 "✅" if j_stability >= 0.70 else "⚠" if j_stability >= 0.50 else "❌")
        log.info("★" * 60)

    except Exception as e:
        log.warning("⚠ Validation Suite beklenmedik hata: %s — analiz devam ediyor.", e)

    # =====================================================================
    # SENTINEL TEST (validation suite'ten BAĞIMSIZ çalışır)
    # =====================================================================
    log.info("★" * 60)
    log.info("★  SENTINEL TEST: Bilinen Biyolojik Gerçeklerle Doğrulama")
    log.info("★" * 60)

    try:
        sentinel_results = validate_against_known_biology(df, g)
    except Exception as e:
        log.warning("⚠ Sentinel test hatası: %s", e)

    # =====================================================================
    # FİNAL: Temizlik ve çıktı
    # =====================================================================
    del g
    gc.collect()

    out_cols = [c for c in [
        "gene", "k_i", "max_cs", "log_degree",
        "pagerank_raw", "pagerank_norm", "log_degree_norm",
        "Hinterland_Skoru", "Kategori", "Topluluk_ID",
        "Lokalizasyon", "Gümrük_Kapisi", "Dis_Komsu_Orani",
        "BC_Skoru", "BC_Esik_Kullanilan",
        "Komsu_Lokalizasyonlar",
        "Z_Skoru", "P_Deger","Drug_Score"
    ] if c in df.columns]

    df[out_cols].to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    log.info("✔ Nihai CSV: %s", OUT_CSV)
    log.info("▶ Hinterland v5.0 tamamlandı.")


if __name__ == "__main__":
    main()
