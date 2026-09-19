# ---------------------------------------------------------------------------
# biology_logic.py
# Biyolojik bonus/ceza sistemi, GO-CC/MF lokalizasyon çözümlemesi,
# farmakolojik doz-yanıt simülasyonu, enfeksiyon şok dalgası ve
# bilinen biyolojiyle sentinel doğrulama.
# ---------------------------------------------------------------------------
import gc
import math
import time
from pathlib import Path
from typing import Optional, Union

import igraph as ig
import numpy as np
import pandas as pd
import requests

from .config import (
    API_CHUNK_SIZE,
    API_RETRY_MAX,
    API_RETRY_WAIT,
    API_TIMEOUT,
    BASE_DIR,
    BLAST_SCORE_MIN,
    BLOCK_WEIGHT_FRACTION,
    DOSE_RESPONSE_DAMPING,
    DOZLAR,
    EFFICIENCY_SAMPLE_SIZE,
    FORCED_DISTANCE_SCALE,
    GO_CC_HIERARCHY,
    HASAR_BIRINCIL,
    HASAR_UCUNCUL,
    HILL_COEFFICIENT,
    MYGENE_FIELDS,
    MYGENE_SPECIES,
    MYGENE_URL,
    OMNIPATH_SNAPSHOT_PATH,
    PAGERANK_DAMPING,
    PARALOG_BOOST_FACTOR,
    BONUS_ENABLED,
    TF_TARGETS,
    COMPLEX_MEMBERS,
    _ensp_to_symbol,
    _save_symbol_cache,
    _ENSPSYMBOL_CACHE,
    log,
)
from .graph_engine import (
    _calculate_congestion,
    _get_paralogs,
    _mask_inf_distances,
    _restore_edges,
)
from .scientific.efficiency import graph_efficiency_summary, percent_change, strength_to_distance
from .scientific.metadata import attach_metadata_columns, build_reproducibility_metadata, file_snapshot, write_metadata
from .scientific.perturbation import severity_name
from .scientific.redistribution import classify_redistribution, top_positive_mean_pct
from .scientific.statistics import empirical_pvalues_degree_matched
from .scientific.regulatory import omnipath_snapshot_metadata
from .scientific.sentinel import resolve_sentinels, sentinel_summary


# ===========================================================================
# BİYOLOJİK BONUS / CEZA
# ===========================================================================

def apply_biological_bonus(
    target_genes: list,
    g: ig.Graph,
    bonus_factor: float = 0.20,
    snapshot: dict = None
) -> None:
    """
    Biyolojik bonus: TF-Hedef + Kompleks ilişkilerini kenar ağırlıklarına yansıt.
    g.es["weight"] doğrudan güncellenir.
    """
    if not BONUS_ENABLED:
        return

    bonus_applied = 0

    for tg in target_genes:
        tg_symbol = _ensp_to_symbol(tg)
        if not tg_symbol:
            continue

        # Kural 1: TF-Hedef bonusu
        if tg_symbol in TF_TARGETS:
            known_targets = TF_TARGETS[tg_symbol]
            for v in g.vs:
                v_symbol = _ensp_to_symbol(v["name"])
                if v_symbol and v_symbol in known_targets:
                    tg_idx = g.vs.find(name=tg).index if tg in g.vs["name"] else -1
                    if tg_idx != -1:
                        v_idx = g.vs.find(name=v["name"]).index
                        eid = g.get_eid(tg_idx, v_idx, error=False)
                        if eid != -1:
                            if snapshot is not None and eid not in snapshot:
                                snapshot[eid] = (g.es[eid]["weight"], g.es[eid]["distance"])
                            g.es[eid]["weight"] *= (1.0 + bonus_factor)
                            g.es[eid]["distance"] = strength_to_distance(g.es[eid]["weight"])
                            bonus_applied += 1

        # Kural 2: Kompleks bonusu
        if tg_symbol in COMPLEX_MEMBERS:
            tg_complexes = COMPLEX_MEMBERS[tg_symbol]
            for v in g.vs:
                v_symbol = _ensp_to_symbol(v["name"])
                if v_symbol and v_symbol in COMPLEX_MEMBERS:
                    v_complexes = COMPLEX_MEMBERS[v_symbol]
                    if set(tg_complexes) & set(v_complexes):

                        tg_idx = g.vs.find(name=tg).index if tg in g.vs["name"] else -1
                        if tg_idx != -1:
                            v_idx = g.vs.find(name=v["name"]).index
                            eid = g.get_eid(tg_idx, v_idx, error=False)
                            if eid != -1:
                                if snapshot is not None and eid not in snapshot:
                                    snapshot[eid] = (g.es[eid]["weight"], g.es[eid]["distance"])
                                g.es[eid]["weight"] *= (1.0 + bonus_factor)
                                g.es[eid]["distance"] = strength_to_distance(g.es[eid]["weight"])
                                bonus_applied += 1

    if bonus_applied > 0:
        log.info("  [Bonus] %d kenara biyolojik bonus uygulandı.", bonus_applied)

    # Çalışma sonunda cache'i kaydet
    _save_symbol_cache(_ENSPSYMBOL_CACHE)


def apply_functional_penalty(target_genes: list, g: ig.Graph, penalty_factor: float = 0.30, snapshot: dict = None) -> None:

    """Aynı GO-MF sınıfından komşulara ceza (kompanzasyon gürültüsü)."""
    if not hasattr(fetch_localizations_mygene, "_mf_cache"):
        return

    mf_cache = fetch_localizations_mygene._mf_cache
    penalty_applied = 0

    for tg in target_genes:
        tg_mf = mf_cache.get(tg, "Unknown")
        if tg_mf in ("Unknown", "Other"):
            continue

        tg_idx = g.vs.find(name=tg).index if tg in g.vs["name"] else -1
        if tg_idx == -1:
            continue

        for nb_idx in g.neighbors(tg_idx):
            nb_name = g.vs[nb_idx]["name"]
            nb_mf = mf_cache.get(nb_name, "Unknown")

            if nb_mf == tg_mf:
                eid = g.get_eid(tg_idx, nb_idx, error=False)
                if eid != -1:
                    if snapshot is not None and eid not in snapshot:
                        snapshot[eid] = (g.es[eid]["weight"], g.es[eid]["distance"])
                    g.es[eid]["weight"] *= (1.0 - penalty_factor)
                    g.es[eid]["distance"] = strength_to_distance(g.es[eid]["weight"])
                    penalty_applied += 1

    if penalty_applied > 0:
        log.info("  [Penalty] %d kenara fonksiyonel ceza uygulandı.", penalty_applied)


# ===========================================================================
# GÖREV 1: GO-CC API Pipeline
# ===========================================================================

def _resolve_localization(go_cc_data) -> str:
    if not go_cc_data:
        return "Unknown"
    items = [go_cc_data] if isinstance(go_cc_data, dict) else go_cc_data
    terms_lower = [
        str(item.get("term", "")).lower()
        for item in items
        if isinstance(item, dict) and "term" in item
    ]
    if not terms_lower:
        return "Unknown"
    found = {
        canonical for kw, canonical in GO_CC_HIERARCHY
        if any(kw in t for t in terms_lower)
    }
    return " | ".join(sorted(found)) if found else "Unknown"


def _resolve_molecular_function(go_mf_data) -> str:
    """GO Molecular Function verisini basitleştirilmiş kategoriye çevir."""
    if not go_mf_data:
        return "Unknown"
    items = [go_mf_data] if isinstance(go_mf_data, dict) else go_mf_data
    terms_lower = [
        str(item.get("term", "")).lower()
        for item in items
        if isinstance(item, dict) and "term" in item
    ]
    if not terms_lower:
        return "Unknown"

    for term in terms_lower:
        if "channel" in term or "transporter" in term:  return "Ion_Channel"
        if "transcription" in term or "dna binding" in term: return "Transcription_Factor"
        if "kinase" in term:    return "Kinase"
        if "receptor" in term:  return "Receptor"
        if "structural" in term or "cytoskeleton" in term: return "Structural"
        if "enzyme" in term or "catalytic" in term: return "Enzyme"
        if "ligand" in term or "hormone" in term: return "Signaling"
        if "rna" in term:       return "RNA_Binding"
        if "ubiquitin" in term: return "Ubiquitin_Related"
    return "Other"


def fetch_localizations_mygene(genes: list) -> dict:
    """
    GO-CC lokalizasyon verisi — mygene.info üzerinden.
    Offline Fallback: API başarısız olursa yerel 'lokalizasyon_yedek.csv' kullanılır.
    """
    log.info("Görev 1: MyGene.info GO-CC (%d gen)…", len(genes))
    result: dict[str, str] = {}
    total_chunks = (len(genes) + API_CHUNK_SIZE - 1) // API_CHUNK_SIZE

    for ci in range(total_chunks):
        chunk     = genes[ci * API_CHUNK_SIZE:(ci + 1) * API_CHUNK_SIZE]
        log.info("  Chunk %d/%d (%d gen)…", ci + 1, total_chunks, len(chunk))
        map_clean = {str(g_).replace("9606.", ""): str(g_) for g_ in chunk}
        payload   = {
            "q"      : ",".join(map_clean.keys()),
            "scopes" : "ensemblprotein",
            "fields" : MYGENE_FIELDS,
            "species": MYGENE_SPECIES,
            "size"   : str(len(map_clean)),
        }
        response = None
        for attempt in range(1, API_RETRY_MAX + 1):
            try:
                response = requests.post(MYGENE_URL, data=payload, timeout=API_TIMEOUT)
                if response.status_code == 200:
                    break
                time.sleep(API_RETRY_WAIT * attempt)
            except requests.exceptions.RequestException:
                time.sleep(API_RETRY_WAIT * attempt)

        if response is None or response.status_code != 200:
            for gn in chunk:
                result[gn] = "Unknown"
            continue

        try:
            hits = response.json()
            if isinstance(hits, dict):
                hits = hits.get("hits", [])
            for hit in hits:
                if not isinstance(hit, dict) or hit.get("notfound"):
                    continue
                orig = map_clean.get(hit.get("query"))
                if orig:
                    go_cc        = hit.get("go", {}).get("CC") if "go" in hit else None
                    result[orig] = _resolve_localization(go_cc)

                    # ─── GO-MF kaydı ───
                    go_mf = hit.get("go", {}).get("MF") if "go" in hit else None
                    mf_result = _resolve_molecular_function(go_mf)
                    if not hasattr(fetch_localizations_mygene, "_mf_cache"):
                        fetch_localizations_mygene._mf_cache = {}
                    fetch_localizations_mygene._mf_cache[orig] = mf_result

        except Exception as e:
            log.error("JSON Error: %s", e)
            for gn in chunk:
                result[gn] = "Unknown"

    for gn in genes:
        result.setdefault(gn, "Unknown")

    fallback_path = BASE_DIR / "data" / "raw" / "lokalizasyon_yedek.csv"
    if fallback_path.exists():
        try:
            fallback_df = pd.read_csv(fallback_path)
            if "gene" in fallback_df.columns and "Lokalizasyon" in fallback_df.columns:
                fallback_map = dict(zip(
                    fallback_df["gene"].astype(str),
                    fallback_df["Lokalizasyon"].astype(str),
                ))
                kurtarilan = 0
                for gn in genes:
                    if result[gn] == "Unknown" and gn in fallback_map:
                        result[gn] = fallback_map[gn]
                        kurtarilan += 1
                if kurtarilan > 0:
                    log.info("    [Offline Fallback] %d gen yerel dosyadan kurtarıldı.", kurtarilan)
        except Exception as e:
            log.warning("    ⚠ Offline fallback okunamadı: %s", e)

    return result


def add_localization_to_df(df: pd.DataFrame, loc_map: dict) -> pd.DataFrame:
    df = df.copy()
    df["Lokalizasyon"] = df["gene"].map(loc_map).fillna("Unknown")
    log.info(
        "  Lokalizasyon dağılımı:\n%s",
        df["Lokalizasyon"].value_counts().to_string(),
    )
    return df


# ===========================================================================
# Hill İnhibisyon
# ===========================================================================

def _hill_inhibition(
    baseline: float,
    survival_fraction: float,
    n: int = HILL_COEFFICIENT
) -> float:
    """
    Hill inhibisyon modeli (düzeltilmiş).

    survival_fraction : Kalan aktivite fraksiyonu (0–1 arası).
                        0.90 → %10 inhibisyon; 0.01 → %99 inhibisyon.
    IC50 (bu kodda)   : sinyal kaybının %50 olduğu survival_fraction değeri.
                        IC50 büyükse (ör. 0.80) sistem hassas;
                        IC50 küçükse (ör. 0.05) sistem dirençlidir.

    Hill denklemi:
        inhibition_frac = 1 - survival_fraction
        response = inhibition_frac^n / (IC50^n + inhibition_frac^n)

    survival_fraction=1.0 → response=0 → baseline korunur.
    """
    inhibition_frac = 1.0 - survival_fraction
    ic50_frac       = 0.5  # IC50 = %50 inhibisyon noktası
    hill_factor     = (inhibition_frac ** n) / (ic50_frac ** n + inhibition_frac ** n)
    return baseline * (1.0 - hill_factor)


def _get_paralog_boost_factor(
    survival_fraction: float,
    base_boost: float = PARALOG_BOOST_FACTOR
) -> float:
    """
    Sürekli paralog boost faktörü.

    Hafif inhibisyonda düşük boost, ağır inhibisyonda yüksek boost.
    sf=1.0 → boost=0 (inhibisyon yok)
    sf=0.01 → boost=base_boost (maksimum inhibisyon)
    """
    inhibition_frac = 1.0 - survival_fraction
    # Sigmoid-benzeri sürekli fonksiyon
    return base_boost * (inhibition_frac ** 2) / (0.25 + inhibition_frac ** 2)


def _estimate_ic50(dose_df: pd.DataFrame) -> Union[str, float]:
    """
    IC50 tahmini: sinyal kaybının %50 olduğu survival_fraction noktası.

    Veriler artan inhibisyon sırasıyla (azalan survival_fraction) sıralanır.
    İkili interpolasyonla %50 kayıp noktası bulunur.
    Hiçbir noktada %50'ye ulaşılamamışsa "Resistant" döner.
    """
    df = dose_df.sort_values("Survival_Fraction").reset_index(drop=True)
    losses = df["Signal_Kayip_Pct"].values
    sfs = df["Survival_Fraction"].values

    # Edge case: ilk dozda bile %50'den fazla kayıp → çok hassas sistem
    if losses[0] > 50.0:
        return float(sfs[0])

    # Edge case: son dozda bile %50'den az kayıp → dirençli sistem
    if losses[-1] < 50.0:
        return "Resistant"

    # Standart interpolasyon
    for i in range(len(df) - 1):
        if losses[i] <= 50.0 <= losses[i + 1]:
            frac = (50.0 - losses[i]) / (losses[i + 1] - losses[i] + 1e-12)
            return float(sfs[i] + frac * (sfs[i + 1] - sfs[i]))

    return "Resistant"


# ===========================================================================
# GÖREV A: Farmakolojik Doz-Yanıt Simülasyonu
# ===========================================================================

def run_pharmacological_dose_response(
    g: ig.Graph,
    df: pd.DataFrame,
    out_dir: Path,
    spesifik_hedefler: list = None,
    suffix: str = "",
    hedef_lokalizasyon: str = None,
    dozlar: list = None,
    sample_size: int = EFFICIENCY_SAMPLE_SIZE,
    hill_n: int = HILL_COEFFICIENT,
    paralog_boost: float = PARALOG_BOOST_FACTOR,
    damping: float = DOSE_RESPONSE_DAMPING,
) -> pd.DataFrame:
    """
    Farmakolojik Doz-Yanıt + Congestion Simülasyonu (düzeltilmiş).

    dozlar  : kalan aktivite fraksiyonu listesi (survival fraction).
    reset_vec her PR çağrısından önce dinamik olarak g.vcount() ile üretiliyor.
    Paralog boost: sürekli doz-bağımlı (sf düştükçe artar).
    """
    # Görselleştirme fonksiyonları metrics.py'da; döngüsel import'tan kaçınmak için
    # burada içeride import ediliyor.
    from .metrics import plot_congestion, plot_dose_response

    if dozlar is None:
        dozlar = DOZLAR

    log.info(
        "Görev A: Farmakolojik Doz-Yanıt + Congestion "
        "(Hill n=%d, kalan_aktivite_fraksiyonları=%s, sürekli paralog boost)…",
        hill_n, dozlar,
    )

    # Legacy protein IDs local symbol registry ile gÃ¼ncel graph node'una Ã§Ã¶zÃ¼lÃ¼r.
    # UNMAPPED / AMBIGUOUS durumlarÄ± FAIL olarak sayÄ±lmaz.
    sentinel_frame = resolve_sentinels(
        df, BASE_DIR / "data" / "processed" / "ensp_with_symbols.csv"
    )
    summary = sentinel_summary(sentinel_frame)
    for record in sentinel_frame.to_dict("records"):
        log.info(
            "  [SENTINEL] %s -> %s | %s | %s",
            record["original_identifier"], record["resolved_network_node_id"] or "-",
            record["status"], record["mapping_source"],
        )
    log.info(
        "  Sentinel: %d/%d PASS (deÄŸerlendirilen), coverage %d/%d",
        summary["evaluated_success_count"], summary["evaluated_count"],
        summary["evaluated_count"], summary["total_count"],
    )
    all_names = set(g.vs["name"])

    if spesifik_hedefler:
        target_genes = []
        seen_targets = set()
        for gn in spesifik_hedefler:
            resolved = None
            if gn in all_names:
                resolved = gn
            elif f"9606.{gn}" in all_names:
                resolved = f"9606.{gn}"
            elif str(gn).replace("9606.", "") in all_names:
                resolved = str(gn).replace("9606.", "")
            if resolved is not None and resolved not in seen_targets:
                target_genes.append(resolved)
                seen_targets.add(resolved)
        log.info("[HEDEF] SNIPER MODU: %s", target_genes)
    elif hedef_lokalizasyon:
        if "Gümrük_Kapisi" in df.columns:
            gates = df[df["Gümrük_Kapisi"] == True]
        else:
            gates = df
        col = "Lokalizasyon"
        if col not in df.columns:
            log.warning("  ⚠ 'Lokalizasyon' sütunu yok.")
            return pd.DataFrame()
        target_genes = gates[
            gates[col].str.contains(hedef_lokalizasyon, na=False)
        ]["gene"].tolist()
        log.info("  💊 LOKALIZASYON MODU: %s (%d hedef)", hedef_lokalizasyon, len(target_genes))
    else:
        log.warning("  ⚠ Hedef belirtilmedi!")
        return pd.DataFrame()

    if not target_genes:
        log.warning("  ⚠ Doz-Yanıt iptal: Hiçbir hedef ağda bulunamadı.")
        return pd.DataFrame()
    bonus_snapshot = {}
    apply_biological_bonus(target_genes, g, snapshot=bonus_snapshot)
    apply_functional_penalty(target_genes, g, snapshot=bonus_snapshot)

    name_to_idx = {v["name"]: v.index for v in g.vs}
    n_targets   = len(target_genes)

    def _make_reset_vec() -> list:
        vec = [0.0] * g.vcount()
        found = 0
        for tg in target_genes:
            idx = name_to_idx.get(tg, -1)
            if idx != -1 and idx < g.vcount():
                vec[idx] = 1.0 / n_targets
                found += 1
        if found < n_targets:
            log.warning(
                "  ⚠ reset_vec: %d/%d hedef ağda bulunamadı.", n_targets - found, n_targets
            )
        return vec

    reset_vec        = _make_reset_vec()
    pr_baseline_vals = g.personalized_pagerank(
        reset=reset_vec,
        damping=damping,
        weights="weight",
        directed=False,
    )
    pr_baseline     = {g.vs[i]["name"]: pr_baseline_vals[i] for i in range(g.vcount())}
    baseline_signal = sum(pr_baseline.get(tg, 0.0) for tg in target_genes)
    log.info("  Baseline sinyal (RWR hedef PR toplamı): %.6f", baseline_signal)

    original_g = g.copy()

    gate_genes_list = (
        df[df["Gümrük_Kapisi"] == True]["gene"].astype(str).tolist()
        if "Gümrük_Kapisi" in df.columns else []
    )

    records = []
    for sf in sorted(dozlar, reverse=True):
        inhibition_pct = round((1 - sf) * 100, 1)
        log.info(
            "  Kalan Aktivite=%.2f uygulanıyor (%.0f%% inhibisyon)…",
            sf, inhibition_pct,
        )

        snapshot: dict[int, tuple] = {}
        paralog_snapshot: dict[int, float] = {}

        # SÜREKLİ PARALOG BOOST - her dozda uygulanır
        paralog_factor = _get_paralog_boost_factor(sf, paralog_boost)

        for tg in target_genes:
            if tg not in name_to_idx:
                continue
            tg_idx = name_to_idx[tg]

            # Ana hedef kenarlarını inhibe et
            for nb_idx in g.neighbors(tg_idx):
                eid = g.get_eid(tg_idx, nb_idx, error=False)
                if eid == -1:
                    continue
                if eid not in snapshot:
                    snapshot[eid] = (g.es[eid]["weight"], g.es[eid]["distance"])
                ow, od = snapshot[eid]
                g.es[eid]["weight"]   = _hill_inhibition(ow, sf, hill_n)
                g.es[eid]["distance"] = od / max(sf, 1e-6)

            # SÜREKLİ: Paralog boost her dozda, dozla orantılı
            if paralog_factor > 0.001:  # Eşik değer - çok küçük boost'ları atla
                paralogs = _get_paralogs(g, tg_idx, BLAST_SCORE_MIN)
                for pl_idx in paralogs:
                    for nb2_idx in g.neighbors(pl_idx):
                        eid2 = g.get_eid(pl_idx, nb2_idx, error=False)
                        if eid2 == -1 or eid2 in snapshot:
                            continue
                        if eid2 not in paralog_snapshot:
                            paralog_snapshot[eid2] = g.es[eid2]["weight"]
                        g.es[eid2]["weight"] *= (1.0 + paralog_factor)

        reset_vec_now = _make_reset_vec()
        pr_vals = g.personalized_pagerank(
            reset=reset_vec_now,
            damping=damping,
            weights="weight",
            directed=False,
        )
        pr_dict       = {g.vs[i]["name"]: pr_vals[i] for i in range(g.vcount())}
        target_signal = sum(pr_dict.get(tg, 0.0) for tg in target_genes)
        signal_loss   = (
            (baseline_signal - target_signal) / baseline_signal * 100
            if baseline_signal > 0 else 0.0
        )

        # Congestion hesapla (doğrudan g'yi kullan - snapshot ile geri yüklenecek)
        delta_bc = _calculate_congestion(
            original_g, g,
            gate_genes=gate_genes_list
        )

        stressed_idx   = np.argsort(delta_bc)[-5:]
        stressed_genes = [
            g.vs[i]["name"] for i in stressed_idx if delta_bc[i] > 0
        ]
        max_cong = float(np.max(delta_bc)) if len(delta_bc) > 0 else 0.0

        log.info(
            "    SurvFrac=%.2f (%.0f%% inh.) | Paralogs=%.4f | RWR=%.6f | Kayıp=%.4f%% | MaxCong=%.4f",
            sf, inhibition_pct, paralog_factor, target_signal, signal_loss, max_cong,
        )

        records.append({
            "Survival_Fraction"      : sf,
            "Inhibition_Pct"         : inhibition_pct,
            "Hill_n"                 : hill_n,
            "Paralog_Boost_Factor"   : round(paralog_factor, 6),
            "Hedef_Sayisi"           : n_targets,
            "Baskılanan_Kenar"       : len(snapshot),
            "Paralog_Boost_Kenar"    : len(paralog_snapshot),
            "RWR_Baseline"           : round(baseline_signal, 8),
            "RWR_Signal"             : round(target_signal, 8),
            "Signal_Kayip_Pct"       : round(signal_loss, 4),
            "Congestion_Top_Stressed": ", ".join(stressed_genes),
            "Max_Congestion_Score"   : round(max_cong, 4),
        })

        _restore_edges(g, snapshot)
        for eid2, ow2 in paralog_snapshot.items():
            g.es[eid2]["weight"] = ow2

    _restore_edges(g, bonus_snapshot)
    del original_g
    gc.collect()

    dose_df                       = pd.DataFrame(records)
    ic50                          = _estimate_ic50(dose_df)
    dose_df["IC50_Survival_Frac"] = ic50

    if ic50 == "Resistant":
        log.info(
            "  📊 IC50: Sistem dirençli — %%50 sinyal kaybına ulaşılamadı."
        )
    else:
        log.info(
            "  📊 IC50 (kalan aktivite noktası): %.4f "
            "(bu noktada %%50 sinyal kaybı gerçekleşiyor)",
            ic50,
        )

    dosya_adi = f"farmakolojik_doz_yanit_{suffix}.csv" if suffix else "farmakolojik_doz_yanit_raporu.csv"
    out_path  = out_dir / dosya_adi
    dose_df.to_csv(out_path, index=False, encoding="utf-8-sig")
    log.info("  -> Doz-Yanıt raporu: %s", out_path)

    try:
        plot_dose_response(dose_df, out_dir)
        plot_congestion(dose_df, out_dir)
    except Exception as e:
        log.warning("  ⚠ Grafik hatası: %s", e)

    return dose_df


# ===========================================================================
# GÖREV 9: Enfeksiyon Şok Dalgası
# ===========================================================================

def run_infection_simulation(
    g: ig.Graph,
    df: pd.DataFrame,
    out_dir: Path,
    hedef_lokalizasyon: str = None,
    spesifik_hedefler: list = None,
    sample_size: int = EFFICIENCY_SAMPLE_SIZE,
    block_weight_fraction: float = BLOCK_WEIGHT_FRACTION,  # ← yeni parametre
    damping: float = PAGERANK_DAMPING,
    return_comparison_metrics: bool = False,
    comparison_genes: list | None = None,
    redistribution_mode: str = "top_n",
    exploratory_top_n: int = 100,
    redistribution_percentile: float = 99.0,
    null_iterations: int = 99,
    fdr_alpha: float = 0.05,
    random_seed: int = 42,
    compute_structural_metrics: bool = True,
    efficiency_sample_sources: int = 256,
    local_efficiency_sample_nodes: int = 128,
    tissue: str | None = None,
    tissue_normalization_mode: str = "within_tissue",
    bc_mode: str = "not_computed_in_simulation",
    bc_sample_size: int | None = None,
    edge_evidence_weighting_mode: str = "legacy_edge_modifiers",

) -> pd.DataFrame | tuple[pd.DataFrame, dict]:
    """Difüzyon Şok Dalgası — Trafik Kayması Modeli."""
    global _LATEST_SIGNED_REDISTRIBUTION
    out_dir.mkdir(parents=True, exist_ok=True)
    log.info("Görev 9: Difüzyon Şok Dalgası v5.0")

    all_names = set(g.vs["name"])

    if spesifik_hedefler:
        target_genes = []
        seen_targets = set()
        for gn in spesifik_hedefler:
            resolved = None
            if gn in all_names:
                resolved = gn
            elif f"9606.{gn}" in all_names:
                resolved = f"9606.{gn}"
            elif str(gn).replace("9606.", "") in all_names:
                resolved = str(gn).replace("9606.", "")
            if resolved is not None and resolved not in seen_targets:
                target_genes.append(resolved)
                seen_targets.add(resolved)
        log.info("[HEDEF] SNIPER MODU: %s", target_genes)
    elif hedef_lokalizasyon:
        gates        = df[df["Gümrük_Kapisi"] == True] if "Gümrük_Kapisi" in df.columns else df
        target_genes = gates[
            gates["Lokalizasyon"].str.contains(hedef_lokalizasyon, na=False)
        ]["gene"].tolist()
        log.info("💊 LOKALIZASYON MODU: %s", hedef_lokalizasyon)
    else:
        log.warning("Hedef belirtilmedi!")
        return pd.DataFrame()

    if not target_genes:
        log.warning("⚠️ Simülasyon iptal: Hedef ağda bulunamadı.")
        bos = pd.DataFrame(columns=["gene", "Hasar_Tipi", "Hinterland_Skoru", "Lokalizasyon"])
        bos.to_csv(out_dir / "enfeksiyon_sok_dalgasi_raporu.csv", index=False)
        return bos
    bonus_snapshot = {}
    if edge_evidence_weighting_mode == "legacy_edge_modifiers":
        apply_biological_bonus(target_genes, g, snapshot=bonus_snapshot)
        apply_functional_penalty(target_genes, g, snapshot=bonus_snapshot)
    elif edge_evidence_weighting_mode != "structural_only":
        raise ValueError("edge_evidence_weighting_mode structural_only veya legacy_edge_modifiers olmalıdır")

    name_to_idx = {v["name"]: v.index for v in g.vs}

    # forced_ids: spesifik_hedefler içinden ağda bulunanlar
    forced_ids = set()
    if spesifik_hedefler:
        for gn in spesifik_hedefler:
            clean = str(gn).replace("9606.", "")
            if clean in all_names:
                forced_ids.add(clean)

    pr_before_vals = g.pagerank(damping=damping, weights="weight", directed=False)
    pr_before = {g.vs[i]["name"]: pr_before_vals[i] for i in range(g.vcount())}

    baseline_efficiency = (
        graph_efficiency_summary(
            g, target_genes,
            global_sample_sources=efficiency_sample_sources,
            local_sample_nodes=local_efficiency_sample_nodes,
            random_seed=random_seed,
        ) if compute_structural_metrics else None
    )

    normalized_edges: dict[int, float] = {}
    if forced_ids:
        lcc_dists = g.es["distance"]
        lcc_avg_d = float(np.mean(lcc_dists)) if lcc_dists else 1.0
        for fn in forced_ids:
            if fn not in name_to_idx:
                continue
            fn_idx = name_to_idx[fn]
            for nb_idx in g.neighbors(fn_idx):
                eid = g.get_eid(fn_idx, nb_idx, error=False)
                if eid == -1 or eid in normalized_edges:
                    continue
                od = g.es[eid]["distance"]
                if od > lcc_avg_d * FORCED_DISTANCE_SCALE:
                    normalized_edges[eid] = od
                    g.es[eid]["distance"] = lcc_avg_d * FORCED_DISTANCE_SCALE

    # Çoklu-gen hasarında ortak bir kenar yalnızca bir kez susturulur.
    # Aksi hâlde ikinci hedef için alınan ara snapshot geri yükleme sırasında
    # grafiği kalıcı olarak zayıf bırakabiliyordu.
    bloklanan: dict[int, tuple] = {}
    for tg in target_genes:
        if tg not in name_to_idx:
            continue
        tg_idx = name_to_idx[tg]
        for nb_idx in g.neighbors(tg_idx):
            eid = g.get_eid(tg_idx, nb_idx, error=False)
            if eid == -1:
                continue
            if eid in bloklanan:
                continue
            orig_w = g.es[eid]["weight"]
            orig_d = g.es[eid]["distance"]
            new_w  = orig_w * block_weight_fraction
            new_d = strength_to_distance(new_w)
            g.es[eid]["weight"]   = new_w
            g.es[eid]["distance"] = new_d
            bloklanan[eid] = (orig_w, orig_d)

    log.info(
        "  %d hedef için %d benzersiz kenar bloke edildi.",
        len(target_genes), len(bloklanan),
    )

    pr_after_vals = g.pagerank(damping=damping, weights="weight", directed=False)
    pr_after = {g.vs[i]["name"]: pr_after_vals[i] for i in range(g.vcount())}

    perturbed_efficiency = (
        graph_efficiency_summary(
            g, target_genes,
            global_sample_sources=efficiency_sample_sources,
            local_sample_nodes=local_efficiency_sample_nodes,
            random_seed=random_seed,
        ) if compute_structural_metrics else None
    )

    comparison_metrics = {}
    if return_comparison_metrics and comparison_genes:
        for gene in comparison_genes:
            resolved = gene if gene in name_to_idx else f"9606.{gene}"
            if resolved not in name_to_idx:
                continue
            before = float(pr_before.get(resolved, 0.0))
            after = float(pr_after.get(resolved, 0.0))
            comparison_metrics[resolved] = {
                "pagerank_before": before,
                "pagerank_after": after,
                "delta_pagerank": after - before,
                "delta_pagerank_pct": ((after - before) / before * 100.0) if before else 0.0,
            }

    for eid, (ow, od) in bloklanan.items():
        g.es[eid]["weight"]   = ow
        g.es[eid]["distance"] = od

    observed_abs_delta = {
        name: abs(float(pr_after.get(name, 0.0)) - float(pr_before.get(name, 0.0)))
        for name in pr_before
    }
    empirical_p, null_metadata = empirical_pvalues_degree_matched(
        g, pr_before, observed_abs_delta, target_genes,
        attenuation_fraction=block_weight_fraction,
        damping=damping,
        iterations=null_iterations if redistribution_mode == "null_fdr" else 0,
        random_seed=random_seed,
    )
    redistribution = classify_redistribution(
        pr_before, pr_after,
        empirical_p=empirical_p,
        target_names=target_genes,
        selection_mode=redistribution_mode,
        # Eski pozitif-artış seçimini uygulayabilmek için Top-N modunda
        # önce tüm düğümlerin iki-yönlü tablosunu oluştur; aşağıda
        # pozitif eşik ve kullanıcı N limiti deterministik olarak uygulanır.
        exploratory_top_n=(len(pr_before) if redistribution_mode == "top_n" else exploratory_top_n),
        percentile=redistribution_percentile,
        fdr_alpha=fdr_alpha,
    )
    # UI/yorumlama katmanı için iki-yönlü ham sonucu koru. Aşağıdaki eski
    # uyumluluk filtresi ve disk CSV'si bilinçli olarak yalnız pozitif listeyi
    # korur; bu ek metadata bilimsel hesaplama veya seçimi değiştirmez.
    signed_redistribution = df.merge(redistribution, on="gene", how="inner")
    # Geri dönüş uyumluluğu: bugünün başındaki Forward Simulation,
    # "stresli" listeyi pozitif PageRank artışı > %0,05 olan genlerden
    # yüzde artışa göre ilk N olarak kuruyordu. Null/FDR ve percentile
    # modları ayrıca seçildiğinde yeni iki-yönlü davranış korunur.
    if redistribution_mode == "top_n" and not redistribution.empty:
        redistribution = (
            redistribution.loc[
                (redistribution["Response_Direction"] == "Influence Gain")
                & (redistribution["Delta_PageRank_Pct"] > 0.05)
            ]
            .sort_values("Delta_PageRank_Pct", ascending=False, kind="mergesort")
            .head(exploratory_top_n)
            .copy()
        )
        redistribution["Selection_Reason"] = (
            f"Eski uyumlu pozitif PageRank Top-{exploratory_top_n}"
        )
    for eid, od in normalized_edges.items():
        g.es[eid]["distance"] = od
    _restore_edges(g, bonus_snapshot)

    # YENİ: Etkilenen düğüm oranı + Total Variation Distance
    pr_diffs = [abs(pr_before.get(n, 0) - pr_after.get(n, 0)) for n in g.vs["name"]]
    tvd = sum(pr_diffs) / 2.0  # Total Variation Distance (0-1 arası)

    # Etkilenen düğüm sayısı (PR değişimi > 1e-6 olanlar)
    etkilenen_düğüm = sum(1 for diff in pr_diffs if diff > 1e-6)
    etkilenen_oran = etkilenen_düğüm / g.vcount() * 100 if g.vcount() > 0 else 0.0

    # Hibrit sistemik kayma: TVD + etkilenen oran
    sistemik_kayma = (tvd * 100 * 0.5) + (etkilenen_oran * 0.5)

    ucuncul_genler = redistribution["gene"].astype(str).tolist() if not redistribution.empty else []
    sok_siddeti = top_positive_mean_pct(redistribution, top_n=100)

    log.info(
        "  BİRİNCİL: %d gen | ÜÇÜNCÜL: %d gen | "
        "Trafik Kayması: %.4f%% | Şok Şiddeti: +%.4f%%",
        len(target_genes), len(ucuncul_genler), sistemik_kayma, sok_siddeti,
    )

    # Ağ ve derece özeti canonical ENSP kimlikleri kullanır; prefiksli ikinci
    # bir eşleme kümesi oluşturmak çoklu hedeflerde yinelenen/boş kayıt üretir.
    target_set = set(target_genes)
    hedef_df   = df[df["gene"].isin(target_set)].copy()
    hedef_df["Hasar_Tipi"] = HASAR_BIRINCIL
    rapor = hedef_df

    if ucuncul_genler:
        ucuncul_set = set(ucuncul_genler)
        uc_df       = df[df["gene"].isin(ucuncul_set)].copy()
        uc_df["Hasar_Tipi"] = HASAR_UCUNCUL
        uc_df = uc_df.merge(redistribution, on="gene", how="left")
        rapor = pd.concat([rapor, uc_df])

    rapor = rapor.sort_values(["Hasar_Tipi", "Hinterland_Skoru"], ascending=[True, False])
    rapor["Systemic_Network_Shift_Pct"] = round(sistemik_kayma, 4)
    rapor["Top_Positive_PageRank_Mean_Pct"] = round(sok_siddeti, 4)
    # Geriye dönük aliaslar: UI'da gösterilmez; yeni gerçek efficiency alanlarına
    # eşitlenmez çünkü tarihsel olarak farklı heuristicleri temsil eder.
    rapor["Efficiency_Kayip_Pct"] = rapor["Systemic_Network_Shift_Pct"]
    rapor["Local_Efficiency_Kayip_Pct"] = rapor["Top_Positive_PageRank_Mean_Pct"]
    if baseline_efficiency and perturbed_efficiency:
        rapor["Global_Efficiency_Baseline"] = baseline_efficiency.global_efficiency
        rapor["Global_Efficiency_Perturbed"] = perturbed_efficiency.global_efficiency
        rapor["Global_Efficiency_Change"] = perturbed_efficiency.global_efficiency - baseline_efficiency.global_efficiency
        rapor["Global_Efficiency_Change_Pct"] = percent_change(baseline_efficiency.global_efficiency, perturbed_efficiency.global_efficiency)
        rapor["Mean_Local_Efficiency_Baseline"] = baseline_efficiency.mean_local_efficiency
        rapor["Mean_Local_Efficiency_Perturbed"] = perturbed_efficiency.mean_local_efficiency
        rapor["Mean_Local_Efficiency_Change_Pct"] = percent_change(baseline_efficiency.mean_local_efficiency, perturbed_efficiency.mean_local_efficiency)
        rapor["Target_Local_Efficiency_Baseline"] = baseline_efficiency.target_local_efficiency
        rapor["Target_Local_Efficiency_Perturbed"] = perturbed_efficiency.target_local_efficiency
        rapor["Target_Local_Efficiency_Change_Pct"] = percent_change(baseline_efficiency.target_local_efficiency, perturbed_efficiency.target_local_efficiency)
        rapor["Global_Efficiency_Mode"] = baseline_efficiency.metadata["global"]["mode"]
        rapor["Mean_Local_Efficiency_Mode"] = baseline_efficiency.metadata["mean_local"]["mode"]
        rapor["Target_Local_Efficiency_Mode"] = baseline_efficiency.metadata["target_local"]["mode"]
    rapor["Kombine_Hasar"]               = len(target_genes) > 1
    rapor["Hasar_Hedefleri"]             = " | ".join(target_genes)
    rapor["Benzersiz_Bloklanan_Kenar"]  = len(bloklanan)
    rapor["Perturbation_Severity"] = severity_name(block_weight_fraction)
    rapor["Edge_Attenuation_Fraction"] = float(block_weight_fraction)
    rapor["Redistribution_Selection_Mode"] = redistribution_mode
    rapor["Edge_Evidence_Weighting_Mode"] = edge_evidence_weighting_mode

    efficiency_metadata = {
        "baseline": baseline_efficiency.metadata if baseline_efficiency else {"mode": "skipped"},
        "perturbed": perturbed_efficiency.metadata if perturbed_efficiency else {"mode": "skipped"},
    }
    omnipath_snapshot = omnipath_snapshot_metadata(OMNIPATH_SNAPSHOT_PATH)
    metadata = build_reproducibility_metadata(
        organism="Homo sapiens", tissue=tissue,
        tissue_normalization_mode=tissue_normalization_mode,
        attenuation_fraction=block_weight_fraction,
        attenuation_severity=severity_name(block_weight_fraction),
        pagerank_damping=damping, bc_mode=bc_mode, bc_sample_size=bc_sample_size,
        null_iterations=int(null_metadata["iterations"]), random_seed=random_seed,
        db_path=BASE_DIR / "data" / "raw" / "hinterland_core.db",
        efficiency_metadata=efficiency_metadata,
        string_version="STRING v12.0 (local 9606 snapshot)",
        expression_source_version="Human Protein Atlas RNA tissue consensus (local snapshot; release metadata absent)",
        trrust_version="TRRUST human local snapshot (release metadata absent)",
        omnipath_version=str(omnipath_snapshot["version"]),
        edge_evidence_weighting_mode=edge_evidence_weighting_mode,
        cache_snapshot_identifiers={
            "symbol_cache": file_snapshot(BASE_DIR / "outputs" / "cache" / "ensp_symbol_cache.pkl"),
            "trrust_regulators": file_snapshot(BASE_DIR / "data" / "processed" / "target_to_regulators.pkl"),
            "omnipath_snapshot": str(omnipath_snapshot.get("checksum_sha256", "missing")),
        },
    )
    rapor = attach_metadata_columns(rapor, metadata)
    # Pandas attrs içine DataFrame koymak concat/nlargest sırasında attrs
    # karşılaştırmasını bozabildiği için signed sonuç modül belleğinde tutulur.
    _LATEST_SIGNED_REDISTRIBUTION = signed_redistribution.copy(deep=True)

    csv_path = out_dir / "enfeksiyon_sok_dalgasi_raporu.csv"
    rapor.to_csv(csv_path, index=False, encoding="utf-8-sig")
    write_metadata(metadata, out_dir / "enfeksiyon_sok_dalgasi_metadata.json")
    log.info("  -> Şok dalgası raporu: %s", csv_path)
    if return_comparison_metrics:
        return rapor, {
            "observed": comparison_metrics,
            "top_pagerank_gains": [
                {"gene": name, "delta_pagerank": float(pr_after[name] - pr_before[name]),
                 "delta_pagerank_pct": float((pr_after[name] - pr_before[name]) / pr_before[name] * 100.0) if pr_before[name] else 0.0}
                for name in sorted(pr_before, key=lambda item: pr_after[item] - pr_before[item], reverse=True)[:exploratory_top_n]
                if pr_after[name] > pr_before[name]
            ],
            "system_shift_pct": float(sistemik_kayma),
            "local_stress_pct": float(sok_siddeti),
            "blocked_edges": len(bloklanan),
            "redistribution": redistribution.to_dict("records"),
            "efficiency": efficiency_metadata,
            "metadata": metadata,
        }
    return rapor


_LATEST_SIGNED_REDISTRIBUTION: pd.DataFrame | None = None


def latest_signed_redistribution() -> pd.DataFrame | None:
    """Son çalışmanın mevcut, önceden hesaplanmış iki-yönlü sonucunu döndür."""
    return _LATEST_SIGNED_REDISTRIBUTION.copy(deep=True) if isinstance(_LATEST_SIGNED_REDISTRIBUTION, pd.DataFrame) else None


# ===========================================================================
# SENTINEL TEST - Bilinen Biyolojik Gerçeklerle Doğrulama
# ===========================================================================

def validate_against_known_biology(df: pd.DataFrame, g: ig.Graph) -> dict:
    """
    Sentinel Test: Kodun çıktısını, literatürdeki bilinen gerçeklerle karşılaştırır.

    Bu test, algoritmanın biyolojik olarak anlamlı sonuçlar üretip
    üretmediğini kontrol eder.
    """
    log.info("=" * 60)
    log.info("[SENTINEL] Bilinen Biyolojik Gerçeklerle Karşılaştırma")
    log.info("=" * 60)

    sentinel_genes = {
        "ENSP00000269305": {
            "min_skor": 90.0,
            "beklenen": "İmparator Hub veya Stratejik Dağıtıcı",
            "gerekçe": "Genom bekçisi, en yüksek bağlantılı proteinlerden",
        },
        "ENSP00000478887": {
            "min_skor": 85.0,
            "beklenen": "İmparator Hub veya Stratejik Dağıtıcı",
            "gerekçe": "Ana transkripsiyon düzenleyici, kanserde kritik",
        },
        "ENSP00000418960": {
            "min_skor": 70.0,
            "beklenen": "En az Stratejik Dağıtıcı",
            "gerekçe": "DNA onarımında kritik rol (BRCA1)",
        },
        "ENSP00000275493": {
            "min_skor": 75.0,
            "beklenen": "En az Stratejik Dağıtıcı",
            "gerekçe": "Sinyal iletiminde ana reseptör",
        },
        "ENSP00000335153": {   # HSP90AA1
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Şaperon, tüm protein ağlarında en yüksek bağlantılılardan",
},
"ENSP00000344818": {   # UBC (Ubiquitin C)
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Binlerce proteinin yıkımını düzenler, ağın en merkezi düğümlerinden",
},
"ENSP00000344456": {   # CTNNB1 (β‑catenin)
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Wnt sinyali ve hücre adezyonunda kritik, çok sayıda etkileşimi var",
},
"ENSP00000382883": {   # AKT1
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "PI3K/AKT yolağının ana kinazı, aşırı bağlantılı",
},
"ENSP00000263253": {   # EP300 (p300)
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Histon asetiltransferaz, transkripsiyon düzenleyici hub",
},
"ENSP00000346847": {   # FN1 (Fibronektin)
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Hücre dışı matrisin ana bileşeni, çok sayıda integrin bağlar",
},
"ENSP00000395213": {   # CDK1
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Hücre döngüsünün ana düzenleyicisi, birçok substratı var",
},
"ENSP00000253840": {   # ESR1 (Östrojen reseptörü)
    "min_skor": 80.0,
    "beklenen": "İmparator Hub",
    "gerekçe": "Nükleer reseptör, transkripsiyon ağlarını yönetir",
},
    }

    all_names = set(g.vs["name"])
    results = {}
    gecti = 0
    kaldi = 0
    veri_yok = 0

    for gene_symbol, kriter in sentinel_genes.items():
        # Ensembl ID'sini bulmaya çalış
        matching_rows = df[df["gene"].str.contains(gene_symbol, case=False, na=False)]

        if matching_rows.empty:
            log.warning("  ⚠ %s: Veri setinde bulunamadı (test atlandı)", gene_symbol)
            results[gene_symbol] = "VERI_YOK"
            veri_yok += 1
            continue

        row = matching_rows.iloc[0]
        skor = row["Hinterland_Skoru"]
        kategori = row["Kategori"]

        if skor >= kriter["min_skor"]:
            log.info("  ✅ %-8s | Skor=%6.1f | %-40s | Beklenen: %s",
                     gene_symbol, skor, kategori, kriter["beklenen"])
            results[gene_symbol] = "GECTI"
            gecti += 1
        else:
            log.warning("  ❌ %-8s | Skor=%6.1f | %-40s | Beklenen: %s | %s",
                       gene_symbol, skor, kategori,
                       kriter["beklenen"], kriter["gerekçe"])
            results[gene_symbol] = "KALDI"
            kaldi += 1

    degerlendirilen = gecti + kaldi
    tanimli_toplam = len(sentinel_genes)
    if degerlendirilen > 0:
        basari_orani = (gecti / degerlendirilen) * 100
        kapsama_orani = (degerlendirilen / tanimli_toplam) * 100
        log.info("  ─────────────────────────────────────────────")
        log.info(
            "  Sentinel Başarı Oranı (değerlendirilebilir): %d/%d (%%%.0f)",
            gecti,
            degerlendirilen,
            basari_orani,
        )
        log.info(
            "  Sentinel Veri Kapsaması: %d/%d (%%%.0f) | veri yok: %d",
            degerlendirilen,
            tanimli_toplam,
            kapsama_orani,
            veri_yok,
        )

        if basari_orani >= 75:
            log.info("  SONUÇ: ✅ Metodoloji biyolojik olarak TUTARLI")
        elif basari_orani >= 50:
            log.warning("  SONUÇ: ⚠ Kısmen tutarlı, parametreler gözden geçirilmeli")
        else:
            log.error("  SONUÇ: ❌ Metodoloji biyolojik gerçeklerle UYUMSUZ!")

    return results
