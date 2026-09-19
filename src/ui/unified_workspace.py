"""Single-page entry to the existing Unified and saved-report workflows."""
from __future__ import annotations

import logging
import streamlit as st

from src.ui.unified_research_report import render_unified_research_report

log = logging.getLogger(__name__)


def render_unified_workspace(*, targets, tissue, is_mouse, output_dir):
    st.header("Modeller arası bulgular ve kayıtlı analizler")
    view = st.radio("Rapor görünümü", ["Unified · Classic + Directed", "Kayıtlı analizler"], horizontal=True, key="unified_workspace_view")
    if view == "Unified · Classic + Directed":
        st.caption("Aynı hedef ve doku için mevcut Classic ve Directed modellerini ayrı çalıştırır. Ham değerleri ortalamaz; Evidence ve CORUM yalnız bağlam sağlar.")
        context = dict(targets=tuple(targets or ()), tissue=tissue, is_mouse=is_mouse,
            block_strength=float(st.session_state.get("block_strength", .001)),
            damping=float(st.session_state.get("damping", .85)),
            bc_sample_sources=st.session_state.get("bc_sample_sources", 1200))
        available = not is_mouse and len(targets or ()) == 1
        if not available:
            st.info("Unified için insan türünde tam olarak bir hedef seçin. Kayıtlı analizler hedef seçiminden bağımsız açılabilir.")
        if st.button("Unified analizi başlat", type="primary", disabled=not available, key="root_run_unified"):
            from src.unified import run_unified_analysis, export_unified_report
            from src.product.history import HistoryService
            st.session_state["unified_research_report"] = None
            st.session_state["unified_research_context"] = None
            try:
                with st.status("Classic ve Directed hesaplanıyor…", expanded=True) as status:
                    report = run_unified_analysis(target=targets[0], tissue=tissue,
                        block_strength=context["block_strength"], damping=context["damping"],
                        bc_sample_sources=context["bc_sample_sources"])
                    export_unified_report(report, output_dir / "unified_reports" / targets[0] / tissue.replace("/", "_"))
                    st.session_state["unified_research_report"] = report
                    st.session_state["unified_research_context"] = context
                    if report.status.value != "FAILED":
                        try:
                            entry = HistoryService().save(report)
                            st.session_state["unified_history_entry"] = entry
                            status.write("Tam analiz dosyası kayıtlı analizlere eklendi.")
                        except Exception:
                            log.exception("Unified history save failed")
                            status.write("Sonuç hazır; kalıcı geçmiş kaydı oluşturulamadı.")
                    status.update(label=f"Unified · {report.status.value}", state="complete" if report.status.value == "COMPLETE" else "error", expanded=False)
            except Exception:
                log.exception("Unified UI execution failed")
                st.error("Unified analizi tamamlanamadı. Hedef/doku ve veri kaynaklarını kontrol edin; ayrıntılar yerel günlükte.")
        report = st.session_state.get("unified_research_report")
        if report is not None and st.session_state.get("unified_research_context") == context:
            render_unified_research_report(report)
        elif report is not None:
            st.info("Seçimler kayıtlı Unified sonucuyla eşleşmiyor. Sonucu güncellemek için analizi yeniden çalıştırın.")
    else:
        from src.ui.saved_analysis_history import render_saved_analysis_history
        render_saved_analysis_history()
