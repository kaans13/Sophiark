"""Gen raporu için tekrar kullanılabilir kart görünümü."""

from __future__ import annotations

import numpy as np
import pandas as pd


def gene_card_html(row: pd.Series) -> str:
    """NaN değerlerini kullanıcıya sızdırmadan bir gen öncelik kartı oluştur."""
    gene_id = str(row.get("gene", "—"))
    symbol = row.get("Symbol", gene_id)
    if pd.isna(symbol) or str(symbol).strip().lower() in ("", "nan", "none"):
        symbol = gene_id

    def number(name: str, precision: int) -> str:
        try:
            value = float(row.get(name))
            return f"{value:.{precision}f}" if np.isfinite(value) else "—"
        except (TypeError, ValueError):
            return "—"

    drug = number("Drug_Score", 1)
    hinterland = number("Hinterland_Skoru", 1)
    bc = number("BC_Skoru", 3)
    return (
        '<div class="sk-gene-card"><div class="sk-gene-card-head">'
        f'<div class="sk-gene-card-symbol">{symbol}</div><div class="sk-gene-card-score-val">{drug}</div>'
        '</div><div class="sk-gene-card-stats">'
        f'<div>Ağdaki Önemi<br><b>{hinterland}</b></div><div>Geçiş Merkeziliği<br><b>{bc}</b></div>'
        '</div><div class="sk-gene-card-footer">'
        f'<span class="sk-gene-card-id">{gene_id}</span></div></div>'
    )
