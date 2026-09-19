"""Presentation for cross-tissue forward simulation comparisons."""

from __future__ import annotations

import logging

import pandas as pd
import streamlit as st

from src.models import TissueDifferentialResult
from src.ui.workspace import render_insight_box
from src.ui.components import render_scientific_note, render_section_header
from src.ui.tables import render_dataframe

log = logging.getLogger(__name__)

LOCAL_EFFICIENCY_COLUMN = "Ortalama Local Efficiency Değişimi (%)"


def _local_efficiency_values(summary: pd.DataFrame) -> pd.Series | None:
    """Read the service's canonical summary metric without display-label guessing."""
    if LOCAL_EFFICIENCY_COLUMN not in summary.columns:
        return None
    values = pd.to_numeric(summary[LOCAL_EFFICIENCY_COLUMN], errors="coerce")
    return values if values.notna().any() else None


def _narrative(result: TissueDifferentialResult) -> str:
    if result.summary.empty:
        return "Karşılaştırmaya alınabilecek doku sonucu oluşmadı."
    local_efficiency = _local_efficiency_values(result.summary)
    if local_efficiency is None:
        return "Doku karşılaştırması üretildi; canonical yerel verimlilik metriği bu sonuç şemasında kullanılamıyor."
    top = result.summary.loc[local_efficiency.idxmax()]
    specific = result.critical_genes
    specific_count = int((specific.get("Yanıt Tipi", pd.Series(dtype=str)) == "Dokuya özgü kritik değişim").sum())
    return (
        f"En yüksek ölçülen ortalama local efficiency değişimi {top['Doku']} dokusunda (%{float(local_efficiency.loc[top.name]):.3f}) görüldü. "
        f"Karşılaştırmada {specific_count} dokuya özgü kritik değişim raporlandı. "
        "Bu bulgular dokuya özgü fenotip değil, doku-ağı üzerindeki hesaplamalı yanıt farklarıdır."
    )


def render_tissue_differential(result: TissueDifferentialResult, *, target_label: str, species: str) -> bool:
    st.markdown("## Dokular Arası Karşılaştırma · BETA")
    st.caption(
        f"Hedef: {target_label} · {species}. BETA karşılaştırması her doku için aynı forward simülasyonu ayrı doku ağı üzerinde çalıştırır. "
        "Bu karşılaştırma deneysel doku fenotipi veya yolak aktivasyonu göstermez."
    )
    for notice in result.notices:
        st.warning(notice)
    if result.summary.empty:
        return st.button("Simülasyon sonucuna dön", key="tissue_differential_close_empty")

    if _local_efficiency_values(result.summary) is None:
        log.warning("Tissue differential summary missing usable canonical column %r; columns=%s", LOCAL_EFFICIENCY_COLUMN, list(result.summary.columns))
        st.warning("Doku özeti gösteriliyor; ancak yerel verimlilik metriği bu sonuç şemasında bulunamadı. Hesaplama sonucu değiştirilmedi.")

    render_insight_box(_narrative(result))
    metrics = st.columns(3)
    metrics[0].metric("Karşılaştırılan doku", len(result.summary))
    metrics[1].metric("Toplam kritik kayıt", len(result.critical_genes))
    response_type = result.critical_genes.get("Yanıt Tipi", pd.Series(dtype=str))
    metrics[2].metric("Dokuya özgü kayıt", int((response_type == "Dokuya özgü kritik değişim").sum()))
    render_section_header("Tissue comparison", "Her doku için aynı forward simulation ayrı ağ üzerinde çalıştırılır.", label="CROSS-TISSUE")
    st.markdown("### Doku bazında temel yanıt")
    render_dataframe(result.summary, label="tissue_differential_summary", key="tissue_summary_table")

    st.markdown("### Kritik gen karşılaştırma matrisi")
    if not result.critical_genes.empty:
        matrix = result.critical_genes.assign(Durum="Kritik").pivot_table(
            index="Gen", columns="Doku", values="Durum", aggfunc="first", fill_value="—"
        ).reset_index()
        render_dataframe(matrix, label="tissue_differential_matrix", key="tissue_matrix_table")
        with st.expander("Kritik gen ayrıntıları", expanded=False):
            render_dataframe(result.critical_genes, label="tissue_differential_critical_genes", key="tissue_critical_table")
    if not result.similarity.empty:
        st.markdown("### Doku çiftleri arasında kritik-gen örtüşmesi")
        render_dataframe(result.similarity, label="tissue_differential_similarity", key="tissue_similarity_table")
        render_scientific_note("Jaccard yalnızca seçilen kritik-gen listelerinin örtüşmesidir; doku fenotipi veya birleşik biyolojik benzerlik skoru değildir.")
    with st.expander("Export", expanded=False):
        render_dataframe(result.summary, label="tissue_differential_summary", key="tissue_summary_export", row_limit=500)
        render_dataframe(result.critical_genes, label="tissue_differential_critical_genes", key="tissue_critical_export", row_limit=500)
        render_dataframe(result.similarity, label="tissue_differential_similarity", key="tissue_similarity_export", row_limit=500)
    return st.button("Simülasyon sonucuna dön", key="tissue_differential_close")
