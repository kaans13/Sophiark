"""Saved-run and human/mouse comparison workspace."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from src.state import invalidate_saved_unified_result, load_saved_unified_result
from src.product.analysis_bundle import AnalysisBundleReader, BundleError
from src.product.history import HistoryService, ImportStatus
from src.services.report_context import cross_species_comparison
from src.ui.metric_presentation import display_status
from src.ui.tables import render_dataframe
from src.ui.workspace import render_insight_box



def render_saved_analysis_history():
    history = HistoryService()
    entries = history.list_entries()
    st.caption(".sophiark dosyaları kalıcı bilimsel kayıttır; SQLite indeks yalnız listeleme için yeniden oluşturulabilir.")
    if st.button("Geçmiş indeksini yeniden oluştur", key="rebuild_sophiark_history"):
        rebuilt = history.rebuild_history_index()
        st.success(f"İndeks yenilendi: {rebuilt.indexed} bundle · {rebuilt.corrupt} bozuk · {rebuilt.duplicate_ids} yinelenen ID")
        entries = history.list_entries()
    if entries:
        labels = {
            f"{entry.target} · {entry.tissue} · {entry.created_at[:19]} · {display_status(entry.analysis_status)}": entry
            for entry in entries
        }
        selected = st.selectbox("Geçmiş analiz", list(labels), index=None, placeholder="Bir kayıt seçin")
        if selected:
            entry = labels[selected]
            bundle_path = history.history_dir / entry.bundle_filename
            controls = st.columns(2)
            if controls[0].button("Kayıtlı raporu aç", key=f"open_{entry.analysis_id}"):
                st.session_state["saved_unified_result"] = history.open_entry(entry.analysis_id)
            controls[1].download_button(
                "Analiz dosyasını paylaş", bundle_path.read_bytes(), entry.bundle_filename,
                "application/vnd.sophiark.analysis+zip", key=f"share_{entry.analysis_id}",
            )
    else:
        st.info("Henüz kalıcı Unified analiz bundle’ı yok.")

    uploaded = st.file_uploader("Analiz Dosyası Aç / İçe Aktar", type=["sophiark"], key="open_sophiark_file")
    if uploaded is not None:
        data = uploaded.getvalue()
        open_col, import_col = st.columns(2)
        if open_col.button("Kopyalamadan aç", key="open_external_sophiark"):
            try:
                load_saved_unified_result(lambda: AnalysisBundleReader().load(data))
            except BundleError as error:
                st.error(f"Analiz dosyası açılamadı: {error.status.value}")
        if import_col.button("Geçmişe içe aktar", key="import_external_sophiark"):
            invalidate_saved_unified_result()
            try:
                imported = history.import_bytes(data)
                if imported.status is ImportStatus.ID_CONFLICT:
                    st.error("Aynı analysis_id farklı içerikle zaten mevcut; dosya değiştirilmedi.")
                else:
                    st.success("Bundle içe aktarıldı." if imported.status is ImportStatus.IMPORTED else "Bu bundle zaten geçmişte mevcut.")
                    if imported.entry is not None:
                        st.session_state["saved_unified_result"] = history.open_entry(imported.entry.analysis_id)
            except BundleError as error:
                st.error(f"Analiz dosyası açılamadı: {error.status.value}")

    saved = st.session_state.get("saved_unified_result")
    if saved is not None:
        if saved.is_historical:
            st.warning("Bu analiz mevcut dataset/config kimliklerinden daha eski bir bilimsel kayıttır; yeniden hesaplanmadan gösteriliyor.")
        from src.ui.unified_research_report import render_unified_research_report
        render_unified_research_report(saved.report)

    snapshots = pd.DataFrame(st.session_state.get("analysis_snapshots", []))
    if not snapshots.empty:
        with st.expander("Legacy oturum özetleri", expanded=False):
            st.caption("LEGACY_HISTORY_NOT_BUNDLE_COMPATIBLE: bu satırlar tam scientific result tablolarını içermez ve silinmemiştir.")
            render_dataframe(snapshots, label="legacy_analiz_gecmisi", key="legacy_analysis_history")
