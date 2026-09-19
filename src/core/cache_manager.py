"""
cache_manager.py — Hinterland v3.8 Cache Yönetimi
====================================================
1. Scorched Earth Protokolü  → tüm st.cache_data + session_state temizler
2. GC Kilidi                 → memory-safe temizlik kancası
3. CSV Yaş Denetçisi        → dosya mtime kontrolü ile stale detection
"""

import gc
import os
import time
from pathlib import Path
from typing import Optional

import streamlit as st


# ─────────────────────────────────────────────────────────────────────────────
# SCORCHED EARTH: temizlenecek session_state anahtarları
# Buraya yeni state ekledikçe listeyi genişlet
# ─────────────────────────────────────────────────────────────────────────────
_STATE_KEYS_TO_PURGE = [
    "rapor_df",
    "fonksiyon_cache",
    "fonksiyon_sonuc",
    "sim_tamam",
    "sim_hedef",
    "disease_bank_search_result",
    "disease_bank_reference",
    "debug_log",
    "hayalet_genler",
    "signed_redistribution_result",
    "enrichment_results",
    "enrichment_notices",
    "enrichment_engine",
    "enrichment_candidate_fingerprint",
    "directed_run_bundle",
    "directed_engine_mode",
    "research_current_simulation_id",
    "research_current_snapshot_id",
    "research_external_snapshot_id",
    "research_simulation_scope",
    "research_snapshot",
    "research_bundle_snapshot_id",
    "research_offline_bundle",
    "research_explorer_open",
    "G_canli",
    "_mygene_query_hash",   # API Guard anahtarı
    "_last_csv_mtime",      # mtime tracker
]


def scorched_earth_reset() -> int:
    """
    Tüm Streamlit cache ve session_state'i acımasızca sil.
    Döner: temizlenen key sayısı (audit için).

    Kullanım:
        n = scorched_earth_reset()
        st.success(f"✓ {n} state anahtarı silindi. RAM serbest bırakıldı.")
        st.rerun()
    """
    # 1. Yalnızca bellek içi Streamlit cache'leri temizlenir. Research
    # evidence SQLite dosyası persistent'tır ve bu reset tarafından silinmez.
    st.cache_data.clear()
    st.cache_resource.clear()

    # 2. session_state temizliği — yalnızca bilinen anahtarlar (güvenli)
    purged = 0
    for key in _STATE_KEYS_TO_PURGE:
        if key in st.session_state:
            del st.session_state[key]
            purged += 1

    # 3. Kalan tüm state'i nükleer temizle (isteğe bağlı — daha agresif)
    # for key in list(st.session_state.keys()):
    #     del st.session_state[key]

    # 4. CPython garbage collector
    collected = force_gc()

    return purged, collected


# ─────────────────────────────────────────────────────────────────────────────
# GC KİLİDİ — Bellek Zorlama
# ─────────────────────────────────────────────────────────────────────────────

def force_gc() -> int:
    """
    3 nesil GC döngüsü çalıştır.
    Ağır DataFrame atamaları sonrasında çağır.
    Döner: toplanan nesne sayısı.
    """
    collected = 0
    for generation in range(3):
        collected += gc.collect(generation)
    return collected


def gc_after_dataload(func):
    """
    Decorator: fonksiyon bitince otomatik GC tetikler.
    csv_yukle_debug ve motoru_isit gibi ağır yükleyicilere ekle.

    Örnek:
        @gc_after_dataload
        def csv_yukle_debug(dosya_yolu): ...
    """
    def wrapper(*args, **kwargs):
        result = func(*args, **kwargs)
        force_gc()
        return result
    wrapper.__name__ = func.__name__
    return wrapper


# ─────────────────────────────────────────────────────────────────────────────
# CSV YAŞ DENETÇİSİ — Stale Cache Dedektörü
# ─────────────────────────────────────────────────────────────────────────────

def csv_stale_check(dosya_yolu: Path) -> bool:
    """
    Dosyanın mtime'ı son bilinen mtime'dan farklıysa True döner.
    Bu durumda cache kırılmalı ve dosya taze okunmalıdır.

    session_state["_last_csv_mtime"] → {str(path): mtime_float}
    """
    key = "_last_csv_mtime"
    if key not in st.session_state:
        st.session_state[key] = {}

    yol_str = str(dosya_yolu)
    try:
        current_mtime = os.path.getmtime(yol_str)
    except FileNotFoundError:
        return False

    last_mtime = st.session_state[key].get(yol_str, None)
    st.session_state[key][yol_str] = current_mtime

    if last_mtime is None:
        return False  # İlk okuma, stale değil

    return abs(current_mtime - last_mtime) > 0.5  # 0.5 saniye tolerans


def auto_invalidate_on_change(dosya_yolu: Path) -> bool:
    """
    CSV değişmişse cache'i otomatik temizle ve True döner.
    app.py'nin başına koy — her rerun'da çalışır.

    Örnek:
        if auto_invalidate_on_change(SOK_RAPORU):
            st.toast("⚡ Yeni simülasyon verisi tespit edildi. Cache yenilendi.", icon="🔄")
    """
    if csv_stale_check(dosya_yolu):
        st.cache_data.clear()
        # Sadece veri state'lerini temizle, UI state'ini koru
        for key in ["rapor_df", "fonksiyon_sonuc", "sim_tamam"]:
            if key in st.session_state:
                del st.session_state[key]
        force_gc()
        return True
    return False
