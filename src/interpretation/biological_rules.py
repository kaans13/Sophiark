"""Conservative wording and schema helpers for the interpretation layer.

Nothing in this module calculates, ranks, filters, or alters a scientific
result.  It only reads already-produced values and selects safe descriptions.
"""

from __future__ import annotations

from typing import Iterable

import pandas as pd


LIMITATIONS = (
    "Bu sonuçlar hesaplamalı ağ yeniden dağılımını tanımlar. Pozitif ΔPageRank; "
    "transkripsiyonel aktivasyon, artmış protein bolluğu, biyokimyasal telafi "
    "veya fonksiyonel kurtarma anlamına gelmez. Biyolojik nedensellik için "
    "deneysel doğrulama gereklidir."
)


def first_number(frame: pd.DataFrame, columns: Iterable[str]) -> float | None:
    """Return the first available finite numeric value without modifying frame."""
    if frame is None or frame.empty:
        return None
    for column in columns:
        if column not in frame.columns:
            continue
        value = pd.to_numeric(frame[column], errors="coerce").iloc[0]
        if pd.notna(value):
            return float(value)
    return None


def first_text(frame: pd.DataFrame, columns: Iterable[str]) -> str | None:
    if frame is None or frame.empty:
        return None
    for column in columns:
        if column not in frame.columns:
            continue
        value = frame[column].iloc[0]
        if pd.notna(value) and str(value).strip():
            return str(value).strip()
    return None


def fmt_percent(value: float | None, digits: int = 4, signed: bool = False) -> str:
    if value is None:
        return "veri yok"
    sign = "+" if signed and value > 0 else ""
    return f"{sign}{value:.{digits}f}%"


def directional_effect(row: pd.Series) -> str:
    """Describe only a direction explicitly present in the evidence record."""
    stimulation = bool(row.get("consensus_stimulation") or row.get("stimulation"))
    inhibition = bool(row.get("consensus_inhibition") or row.get("inhibition"))
    if stimulation and inhibition:
        return "çelişkili yönlü kanıt"
    if stimulation:
        return "uyarıcı düzenleme kaydı"
    if inhibition:
        return "baskılayıcı düzenleme kaydı"
    if bool(row.get("consensus_direction")) or str(row.get("direction") or "").strip():
        return "yönlü düzenleme kaydı (işaret belirtilmemiş)"
    return "yön bilgisi olmayan kayıt"
