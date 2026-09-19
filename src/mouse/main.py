# ---------------------------------------------------------------------------
# mouse_main.py
# Fare (Mus musculus) için ana orkestratör ve Streamlit köprüsü.
# İnsan main.py'den tamamen bağımsızdır; ayrı çıktı klasörü ve ayrı veritabanı kullanır.
# ---------------------------------------------------------------------------
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
    OUT_CSV,
    OUT_DIR,
    PAGERANK_DAMPING,
    PAGERANK_MAX_ITER,
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
)
from .metrics import (
    build_hinterland_score,
    classify_genes,
    detect_customs_gates,
)
from .data_loader import fetch_degree_summary
from .config import (
    MOUSE_MYGENE_CACHE,
    MOUSE_SYMBOL_MAP,
    format_essentiality,
    format_regulators,
)


def _add_mouse_evidence_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Fare çıktısını sembol, GO ve translasyonel kanıtla zenginleştir."""
    df = df.copy()
    df["Symbol"] = df["gene"].map(MOUSE_SYMBOL_MAP).fillna(df["gene"])
    cache_records = df["gene"].map(MOUSE_MYGENE_CACHE)
    df["Protein_Adi"] = cache_records.apply(
        lambda record: record.get("name", "—") if isinstance(record, dict) else "—"
    )
    df["GO_CC_Terimleri"] = cache_records.apply(
        lambda record: " | ".join(record.get("go_cc", [])[:8]) if isinstance(record, dict) else "—"
    )
    df["GO_MF_Terimleri"] = cache_records.apply(
        lambda record: " | ".join(record.get("go_mf", [])[:8]) if isinstance(record, dict) else "—"
    )
    df["Düzenleyici_TFler"] = df["Symbol"].apply(format_regulators)
    df["Essentiality"] = df["Symbol"].apply(format_essentiality)
    df["Yapisal_Kanit"] = "Fareye özgü yapısal/ligand veritabanı henüz entegre değil"
    return df


# ===========================================================================
# Streamlit Köprüsü
# ===========================================================================

def arayuz_icin_motoru_hazirla(
    forced_genes: list = None,
    hedef_doku: str = None,
    bc_sample_sources: int | None = None,
    tissue_normalization_mode: str = "within_tissue",
):
    """Fare için seçili dokuya ait grafiği ve skor tablosunu kur."""
    g = build_high_conf_subgraph(
        DB_PATH,
        HIGH_CONF_THRESHOLD,
        forced_genes=forced_genes,
        hedef_doku=hedef_doku,
        tissue_normalization_mode=tissue_normalization_mode,
    )

    df = fetch_degree_summary(DB_PATH, HIGH_CONF_THRESHOLD)
    df = df[df["gene"].isin(set(g.vs["name"]))].reset_index(drop=True)
    if df.empty:
        return g, df

    pr = compute_pagerank(g, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
    df = build_hinterland_score(df, pr)
    df = classify_genes(df)
    df, _ = run_community_detection(g, df)

    # Lokalizasyonları MyGene'den al
    non_noise = df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]["gene"].tolist()
    loc_map = fetch_localizations_mygene(non_noise)
    for gn in df[df["Kategori"] == "Gri (İstatistiksel Gürültü)"]["gene"]:
        loc_map.setdefault(gn, "Unknown")
    df = add_localization_to_df(df, loc_map)
    g = tag_graph_nodes(g, dict(zip(df["gene"], df["Lokalizasyon"])))

    df = detect_customs_gates(g, df, bc_sample_sources=bc_sample_sources)


    # Farede yapısal veritabanı yok; Drug_Score boş bırak
    df["Drug_Score"] = None
    return g, _add_mouse_evidence_columns(df)


# ===========================================================================
# argparse
# ===========================================================================

def _parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        prog="mouse_hinterland_analysis",
        description="Mouse Hinterland Protein Interaction Network v1.0",
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
        help="Headless Batch-Sweep modu. N > 0 ise top-N gen teker teker simüle edilir.",
    )
    return parser.parse_args()


def _run_batch_sweep(g, df: pd.DataFrame, out_dir, top_n: int) -> None:
    log.info("━" * 60)
    log.info("▶ MOUSE BATCH-SWEEP MODU: top-%d gen, headless", top_n)
    log.info("━" * 60)

    rapor_yolu = out_dir / "reports" / "mouse_tekil_zafiyet_raporu.csv"
    sutunlar = [
        "Hedef_Gen", "Systemic_Network_Shift_Pct", "Top_Positive_PageRank_Mean_Pct",
        "Efficiency_Kayip_Pct", "Local_Efficiency_Kayip_Pct",
    ]
    write_header = not rapor_yolu.exists()

    elit = (
        df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]
        .nlargest(top_n, "Hinterland_Skoru")["gene"]
        .tolist()
    )
    log.info("  Mouse elit havuz: %d gen seçildi.", len(elit))

    for sira, hedef_gen in enumerate(elit, start=1):
        log.info("  [%d/%d] Mouse Simülasyon: %s", sira, len(elit), hedef_gen)
        g_backup = g.copy()

        try:
            sim_df = run_infection_simulation(
                g, df, out_dir,
                spesifik_hedefler=[hedef_gen],
                sample_size=0,
                redistribution_mode="top_n", null_iterations=0,
                compute_structural_metrics=False,
            )
            if not sim_df.empty:
                birincil = sim_df[sim_df["Hasar_Tipi"] == HASAR_BIRINCIL]
                if birincil.empty:
                    birincil = sim_df
                eff_kayip   = float(birincil["Efficiency_Kayip_Pct"].iloc[0])
                local_kayip = float(birincil["Local_Efficiency_Kayip_Pct"].iloc[0])
            else:
                eff_kayip = float("nan")
                local_kayip = float("nan")
        except Exception as exc:
            log.warning("    ⚠ %s simülasyon hatası: %s", hedef_gen, exc)
            eff_kayip = float("nan")
            local_kayip = float("nan")
        finally:
            del g_backup
            gc.collect()

        eff_value = round(eff_kayip, 4) if not math.isnan(eff_kayip) else float("nan")
        local_value = round(local_kayip, 4) if not math.isnan(local_kayip) else float("nan")
        satir = pd.DataFrame(
            [[hedef_gen, eff_value, local_value, eff_value, local_value]],
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

    log.info("  ✔ Mouse Batch-Sweep tamamlandı → %s", rapor_yolu)


# ===========================================================================
# ANA ORKESTRATÖR
# ===========================================================================

def main():
    args = _parse_args()
    hedef_doku = None if args.tissue.strip().lower() == "none" else args.tissue

    log.info("▶ Mouse Hinterland v1.0 başlatılıyor…")
    log.info("  Doku: %s", hedef_doku or "TÜM PROTEOM")
    log.info("  Batch: %s", args.batch if args.batch > 0 else "KAPALI")
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
    g  = tag_graph_nodes(g, dict(zip(df["gene"], df["Lokalizasyon"])))
    gc.collect()

    df = detect_customs_gates(g, df, bc_sample_sources=1200)
    df["Drug_Score"] = None  # farede yapısal veri yok
    df = _add_mouse_evidence_columns(df)

    if args.batch > 0:
        _run_batch_sweep(g, df, OUT_DIR, top_n=args.batch)
        out_cols = [c for c in [
            "gene", "k_i", "max_cs", "Hinterland_Skoru",
            "Kategori", "Lokalizasyon", "Gümrük_Kapisi", "BC_Skoru",
        ] if c in df.columns]
        df[out_cols].to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
        log.info("✔ Mouse Batch CSV: %s", OUT_CSV)
        log.info("▶ Mouse Hinterland v1.0 (Batch) tamamlandı.")
        return

    # Interaktif varsayılan: rastgele 5 non-noise gen seç
    random.seed(42)
    adaylar = df[df["Kategori"] != "Gri (İstatistiksel Gürültü)"]["gene"].tolist()
    top_stress_genes = random.sample(adaylar, min(5, len(adaylar)))
    log.info("🎯 Mouse Stres Testi Hedefleri: %s", top_stress_genes)

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
    if not dose_df.empty:
        log.info("Mouse Doz-Yanıt tamamlandı.")

    del g
    gc.collect()

    out_cols = [c for c in [
        "gene", "k_i", "max_cs", "log_degree",
        "pagerank_raw", "pagerank_norm", "log_degree_norm",
        "Hinterland_Skoru", "Kategori", "Topluluk_ID",
        "Lokalizasyon", "Gümrük_Kapisi", "BC_Skoru", "Symbol",
        "Protein_Adi", "GO_CC_Terimleri", "GO_MF_Terimleri",
        "Düzenleyici_TFler", "Essentiality", "Yapisal_Kanit",
    ] if c in df.columns]

    df[out_cols].to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    log.info("✔ Mouse Nihai CSV: %s", OUT_CSV)
    log.info("▶ Mouse Hinterland v1.0 tamamlandı.")


if __name__ == "__main__":
    main()
