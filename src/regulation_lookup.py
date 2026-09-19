# src/regulation_lookup.py — SIGNOR regülasyon bilgisi arama (v3 - kesin çözüm)

from pathlib import Path
import pandas as pd

_REGULATION_MAP = None  # (source, target) -> effect

def _load_regulation_map():
    """SIGNOR dosyasını yükler, (source, target) -> effect sözlüğü oluşturur."""
    global _REGULATION_MAP
    if _REGULATION_MAP is not None:
        return _REGULATION_MAP

    dosya = Path(__file__).resolve().parent.parent / "data" / "processed" / "signor_relations.tsv"
    if not dosya.exists():
        _REGULATION_MAP = {}
        return _REGULATION_MAP

    df = pd.read_csv(dosya, sep="\t", dtype=str)
    _REGULATION_MAP = {}
    for _, row in df.iterrows():
        key = (row["source"].strip().upper(), row["target"].strip().upper())
        _REGULATION_MAP[key] = row["effect"]  # "Aktivasyon" veya "İnhibisyon"
    return _REGULATION_MAP

def get_regulation(source_symbol: str, target_symbol: str) -> str:
    """
    İki gen sembolü arasındaki düzenleyici ilişkiyi döndürür.

    Parametreler:
        source_symbol: Hedef genin sembolü (örn. "TP53")
        target_symbol: Stresli genin sembolü (örn. "STK17A")

    Dönüş:
        - "Aktivasyon"        : source → target aktive eder
        - "İnhibisyon"        : source → target inhibe eder
        - "Aktivasyon (ters)" : target → source aktive eder
        - "İnhibisyon (ters)" : target → source inhibe eder
        - "Bilinmiyor"        : İki yönde de kayıt yok
    """
    if not source_symbol or not target_symbol:
        return "Bilinmiyor"

    reg_map = _load_regulation_map()
    if not reg_map:
        return "Bilinmiyor"

    src = source_symbol.strip().upper()
    tgt = target_symbol.strip().upper()

    # 1. Hedef → Stresli yönü
    forward = (src, tgt)
    if forward in reg_map:
        return reg_map[forward]  # "Aktivasyon" veya "İnhibisyon"

    # 2. Stresli → Hedef yönü (ters)
    reverse = (tgt, src)
    if reverse in reg_map:
        return f"{reg_map[reverse]} (ters)"

    return "Bilinmiyor" 