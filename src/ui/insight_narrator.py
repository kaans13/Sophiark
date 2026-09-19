"""Narrates already-calculated results; it never makes new scientific claims."""

from __future__ import annotations

import pandas as pd


def en_onemli_bulguyu_ozetle(results: pd.DataFrame, *, mode: str, symbol_col: str = "Symbol") -> tuple[str, bool]:
    if results is None or results.empty:
        return "Henüz yorumlanacak bir analiz sonucu bulunmuyor.", False
    if mode == "Target Stress Search":
        top = results.iloc[0]
        directed = str(top.get("Yönlü Kanıt", "—") or "—")
        return (f"En güçlü doğrulanmış aday {top.get('Aday', '—')}: hedef üzerindeki öngörülen ağ stresi "
                f"{float(top.get('Öngörülen Hedef Stresi', 0)):.3f}. "
                f"Yönlü kanıt durumu: {directed}.", directed not in ('', '—', 'nan'))
    if mode == "Compensation Analysis":
        top = results.iloc[0]
        symbol = top.get(symbol_col, top.get("Gene", top.get("Aday", "—")))
        return (f"En yüksek ağ önemi artışı: {symbol} (+%{float(top.get('ΔPageRank (%)', 0)):.2f}).", False)
    ordered = results.copy()
    if "Drug_Score" in ordered.columns:
        ordered = ordered.sort_values("Drug_Score", ascending=False)
    elif "Hinterland_Skoru" in ordered.columns:
        ordered = ordered.sort_values("Hinterland_Skoru", ascending=False)
    top = ordered.iloc[0]
    symbol = top.get(symbol_col, top.get("gene", "—"))
    function = str(top.get("MyGene Adı", "") or "").strip()
    stressed = int((results.get("Hasar_Tipi", pd.Series(dtype=str)).astype(str).str.contains("ÜÇÜNCÜL", case=False, na=False)).sum())
    label = f"{symbol} ({function})" if function and function not in {"—", "nan", "None"} else str(symbol)
    if str(top.get("Engine", "")).casefold() == "directed":
        return (
            f"Directed yeniden-dağılım adaylarında öne çıkan değişim {label} üzerinde görüldü; "
            f"Directed PageRank değişimi %{float(top.get('Delta_PageRank_Pct', 0)):.2f}. "
            "Classic verimlilik veya Hinterland metriği bu özete taşınmadı.",
            len(results) >= 20,
        )
    if stressed == 0:
        return ("Bu perturbasyonda stresli gen sınıfına giren kayıt oluşmadı. "
                "Aşağıdaki ağ ölçümleri yine de perturbasyonun yapısal etkisini gösterir.", False)
    priority = "Drug_Score" if "Drug_Score" in ordered else "Hinterland skoru" if "Hinterland_Skoru" in ordered else "kaynak sonuç"
    return (f"Ağda {stressed} yeniden-dağılım adayı kaydedildi. Mevcut {priority} sırasındaki ilk kayıt {label}; "
            f"başlangıç ağındaki önemi {float(top.get('Hinterland_Skoru', 0)):.1f}. "
            "Bu öncelik sırası yanıt büyüklüğüyle aynı değildir; en yüksek göreli yanıtlar aşağıda ayrıca gösterilir.", False)
