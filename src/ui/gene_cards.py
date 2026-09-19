"""Progressive-disclosure gene cards for result surfaces."""

from __future__ import annotations

import math
import pandas as pd
import streamlit as st

from src.ui.design_tokens import render_badge
from src.ui.localization import locale, t


def _value(row: pd.Series, key: str, digits: int = 2) -> str:
    value = row.get(key)
    try:
        return f"{float(value):.{digits}f}" if math.isfinite(float(value)) else "—"
    except (TypeError, ValueError):
        return "—"


def _clean_text(value: object) -> str:
    try:
        if pd.isna(value):
            return ""
    except (TypeError, ValueError):
        pass
    text = str(value).strip()
    return "" if text.casefold() in {"", "nan", "none", "null"} else text


def render_gene_card(gene_data: pd.Series, *, key: str) -> None:
    """Visible interpretation + compact metrics; provenance stays collapsed."""
    is_en = locale() == "en"
    entity_id = _clean_text(gene_data.get("gene"))
    symbol = _clean_text(gene_data.get("Symbol"))
    if not symbol or symbol == entity_id:
        symbol = "Alias unavailable" if is_en else "Alias bulunamadı"
    score = float(gene_data.get("Drug_Score", 0) or 0)
    tier = ("Kritik", "tier_kritik") if score >= 65 else ("Orta", "tier_orta") if score >= 40 else ("Düşük", "tier_dusuk")
    mygene_name = str(gene_data.get("MyGene Adı") or "—")
    st.markdown(f"### {symbol} {render_badge(tier[0], tier[1])}", unsafe_allow_html=True)
    if mygene_name != "—":
        st.caption(f"MyGene name: {mygene_name}" if is_en else f"MyGene adı: {mygene_name}")
    cols = st.columns(2)
    is_directed = str(gene_data.get("Engine", "")).casefold() == "directed"
    if is_directed:
        metrics = ("Directed Δ PageRank (%)", "Directed BC")
        values = (_value(gene_data, "Delta_PageRank_Pct", 2), _value(gene_data, "BC_Skoru", 4))
    else:
        metrics = (t("hinterland_score"), t("gateway"))
        gate = str(gene_data.get("Gümrük_Kapisi") or "—")
        values = (_value(gene_data, "Hinterland_Skoru", 1), gate)
    for col, label, value in zip(cols, metrics, values):
        col.metric(label, value)
    with st.expander(t("evidence_details"), expanded=False):
        st.caption(f"ENSP: {entity_id or '—'}")
        st.caption(f"GO biological process: {gene_data.get('GO Biyolojik Süreç', '—')}" if is_en else f"GO biyolojik süreç: {gene_data.get('GO Biyolojik Süreç', '—')}")
        st.caption(f"GO molecular function: {gene_data.get('GO Moleküler İşlev', '—')}" if is_en else f"GO moleküler işlev: {gene_data.get('GO Moleküler İşlev', '—')}")
        st.caption(f"GO cellular component: {gene_data.get('GO Hücresel Bileşen', '—')}" if is_en else f"GO hücresel bileşen: {gene_data.get('GO Hücresel Bileşen', '—')}")
        st.caption(f"TRRUST / regulatory evidence: {gene_data.get('Düzenleyici_TFler', '—')}" if is_en else f"TRRUST / düzenleyici kanıt: {gene_data.get('Düzenleyici_TFler', '—')}")
        if is_directed:
            st.caption("Directed result: Classic Hinterland, efficiency and Compartment Bottleneck metrics are not presented as Directed metrics.")
        else:
            st.caption(f"Compartment Bottleneck: {gene_data.get('Gümrük_Kapisi', '—')} · Essentiality: {gene_data.get('Essentiality', '—')}")
        if not is_directed and locale() == "en":
            st.caption("Network Importance summarizes a gene's structural position in the analyzed baseline network without changing the underlying scientific metric.")
            st.caption("It is not a direct measure of biological essentiality, expression level, or clinical importance.")
            st.caption("Compartment Bottleneck identifies nodes that may occupy bridge-like positions between different network regions or cellular contexts. Such a position may be structurally important for redistribution of perturbation effects across the network.")
            st.caption("Compartment Bottleneck status does not by itself imply biological causality or an obligatory regulatory role.")
        elif not is_directed:
            st.caption("Ağdaki Önemi, genin başlangıç ağındaki yapısal konumunu özetler; bilimsel metrik veya hesaplama değiştirilmez.")
            st.caption("Bu skor biyolojik esansiyellik, ekspresyon düzeyi veya klinik önem ölçüsü değildir.")
            st.caption("Compartment Bottleneck, ağın farklı bölgeleri veya hücresel bağlamları arasındaki akışta köprü rolü gösterebilen düğümleri işaretler. Böyle bir konum, pertürbasyon etkisinin ağın başka bölgelerine taşınmasında yapısal olarak önemli olabilir.")
            st.caption("Bir genin Compartment Bottleneck olarak işaretlenmesi biyolojik nedensellik veya zorunlu düzenleyici rol anlamına gelmez.")
        st.json({k: v for k, v in gene_data.to_dict().items() if not isinstance(v, (list, dict))}, expanded=False)
