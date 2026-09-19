"""Presentation-only view for the Classic + Directed Unified Research Report."""
from __future__ import annotations

import json

import pandas as pd
import streamlit as st

from src.unified import UnifiedResearchReport
from src.ui.metric_presentation import display_status, format_metric
from src.ui.tables import render_dataframe


def _subset(frame: pd.DataFrame, finding: str) -> pd.DataFrame:
    if frame.empty or "finding_class" not in frame:
        return pd.DataFrame()
    return frame.loc[frame["finding_class"].eq(finding)].copy()


def render_unified_research_report(report: UnifiedResearchReport) -> None:
    """Render synthesis/context only; it never recomputes or reorders engine data."""
    st.title("Unified Araştırma Raporu")
    st.caption(f"Hedef: {report.target} · Doku: {report.tissue} · Durum: {display_status(report.status.value)}")
    metrics = st.columns(5)
    metrics[0].metric("Classic", display_status(report.classic_status))
    metrics[1].metric("Directed", display_status(report.directed_status))
    metrics[2].metric("Eşleşen entity", report.agreement.get("aligned_entities", "—"))
    metrics[3].metric("İşaret uyumu", format_metric("sign_agreement_rate", report.agreement.get("sign_agreement_rate")) if report.agreement.get("status") == "AVAILABLE" else "—")
    metrics[4].metric("Evidence provenance", display_status(report.context_availability.get("evidence_provenance", {}).get("status", "UNAVAILABLE")))
    st.info("Classic ve Directed ayrı hesaplamalı topoloji varsayımlarıdır. Evidence yalnızca açıklayıcı provenance sağlar; mevcut Evidence weighting durumu BIAS_WARNING olarak kalır.")
    robust, sensitive, dominant, low_agreement, all_candidates, provenance, complexes, agreement, raw = st.tabs([
        "Sağlam bulgular", "Yön duyarlı", "Motor farkları", "Düşük uyum", "Tüm adaylar", "Evidence provenance", "CORUM bağlamı", "Uyum", "Ham motor sonuçları",
    ])
    with robust:
        rows = _subset(report.candidates, "ROBUST_CROSS_MODEL")
        st.caption("Sağlam bulgu, Classic ve Directed varsayımları arasında tutarlılık demektir; deneysel doğrulama değildir.")
        render_dataframe(rows, label="unified_robust_findings", key="unified_robust", row_limit=100) if not rows.empty else st.info("Seçili sunum kuralında ortak, belirgin ve aynı işaretli bulgu yok.")
    with sensitive:
        rows = _subset(report.candidates, "DIRECTION_SENSITIVE")
        st.caption("Directed modelde yön duyarlıdır; aktivasyon, inhibisyon veya nedensellik kanıtlamaz.")
        render_dataframe(rows, label="unified_direction_sensitive", key="unified_sensitive", row_limit=100) if not rows.empty else st.info("Gösterilen adaylarda yön duyarlı bulgu yok.")
    with dominant:
        rows = report.candidates.loc[report.candidates.get("finding_class", pd.Series(dtype=str)).isin(["CLASSIC_DOMINANT", "DIRECTED_EMERGENT"])].copy()
        render_dataframe(rows, label="unified_engine_emergent", key="unified_emergent", row_limit=100) if not rows.empty else st.info("Gösterilen adaylarda motor-baskın bulgu yok.")
    with low_agreement:
        rows = _subset(report.candidates, "LOW_AGREEMENT")
        st.caption("Düşük uyum, iki motorun bu aday için uyumlu bir bulgu sınıfı üretmediğini gösterir; sonuç sırası veya ham değerler yeniden hesaplanmaz.")
        render_dataframe(rows, label="unified_low_agreement", key="unified_low_agreement", row_limit=100, paginate=True) if not rows.empty else st.info("Bu raporda düşük uyumlu aday yok.")
    with all_candidates:
        st.caption("Canonical Unified aday evreni sınıflandırmaya göre filtrelenmeden, korunmuş motor sırasıyla gezilebilir ve dışa aktarılabilir.")
        render_dataframe(report.candidates, label="unified_all_candidates", key="unified_all_candidates", row_limit=100, paginate=True) if not report.candidates.empty else st.info("Bu raporda aday yok.")
    with provenance:
        columns = [column for column in ["symbol", "entity_id", "finding_class", "Evidence_combined_score", "Evidence_s_nontext", "Evidence_text_dependency", "Evidence_experiments", "Evidence_database", "Evidence_coexpression", "Evidence_Transferred_Channel_Sum", "Evidence_physical_status", "Evidence_Provenance_Status", "CORUM_Target_Complex_Context"] if column in report.candidates]
        st.caption("STRING bileşimi bağlamsal metadata’dır; yanıt, sıra, yüzdelik veya bulgu sınıfını değiştirmez.")
        render_dataframe(report.candidates[columns], label="unified_evidence_provenance", key="unified_provenance", row_limit=100) if columns else st.info("Bu aday kümesi için provenance alanı yok.")
    with complexes:
        st.caption("Kompleks ilişkili ağ yanıtı, fiziksel kompleks aktivasyonu, inhibisyonu veya yıkımı anlamına gelmez.")
        render_dataframe(report.complex_context, label="unified_complex_context", key="unified_complex", row_limit=100) if not report.complex_context.empty else st.info("Hedefi içeren CORUM kompleks bağlamı bulunamadı.")
    with agreement:
        st.json(dict(report.agreement))
    with raw:
        classic, directed = st.tabs(["Classic ham", "Directed ham"])
        with classic:
            render_dataframe(report.classic_raw, label="unified_classic_raw", key="unified_classic_raw", row_limit=500, paginate=True)
        with directed:
            render_dataframe(report.directed_raw, label="unified_directed_raw", key="unified_directed_raw", row_limit=500, paginate=True)
    with st.expander("Makine tarafından okunabilir provenance", expanded=False):
        st.code(json.dumps(dict(report.provenance), indent=2, default=str), language="json")
