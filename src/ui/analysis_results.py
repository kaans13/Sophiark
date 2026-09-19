"""Product-facing forward result layout built over existing output objects."""

from __future__ import annotations

import logging

import pandas as pd
import streamlit as st

from src.services.report_context import critical_rows, symbol_first_report, symbol_map
from src.ui.network_visualization import render_network_views
from src.ui.tables import render_dataframe
from src.ui.workspace import render_insight_box
from src.ui.components import render_context_header, render_scientific_note, render_section_header
from src.ui.localization import locale

log = logging.getLogger(__name__)


def _metric(report: pd.DataFrame, name: str) -> float:
    try:
        return float(report[name].iloc[0])
    except (KeyError, IndexError, TypeError, ValueError):
        return 0.0


def _summary(report: pd.DataFrame, critical: pd.DataFrame, mapping: dict[str, str], targets: list[str]) -> str:
    target_labels = ", ".join(mapping.get(value, value) for value in targets) or "Seçilen lokalizasyon"
    if critical.empty:
        return f"{target_labels} perturbasyonu tamamlandı; tanımlı güçlü yeniden-dağılım eşiğini geçen gen bulunmadı. Bu, biyolojik etkisizlik kanıtı değildir."
    top = critical.iloc[0]
    top_label = mapping.get(str(top.get("gene", "")), str(top.get("Symbol") or "bir gen"))
    return (
        f"{target_labels} perturbasyonu sonrasında {len(critical)} gen güçlü ağ yeniden-dağılımı gösterdi. "
        f"En belirgin değişim {top_label} üzerinde görüldü; bu sonuç hesaplamalı ağ etkisidir, deneysel nedensellik iddiası değildir."
    )


def render_forward_results(*, report: pd.DataFrame, graph, scores: pd.DataFrame, targets: list[str], symbols: pd.DataFrame, runtime, is_mouse: bool, tissue: str) -> None:
    mapping = symbol_map(symbols)
    critical = critical_rows(report)
    display = symbol_first_report(report, mapping)
    target_label = ", ".join(mapping.get(value, value) for value in targets) or "Localization selection"
    render_context_header(target=target_label, tissue=tissue, species="Fare" if is_mouse else "İnsan", attenuation=float(st.session_state.get("block_strength", .001)), candidate_limit=len(targets))
    render_insight_box(_summary(report, critical, mapping, targets))
    metrics = st.columns(4)
    metrics[0].metric("Güçlü değişim", f"{len(critical):,}")
    metrics[1].metric("Global verimlilik değişimi", f"%{_metric(report, 'Global_Efficiency_Change_Pct'):.3f}")
    metrics[2].metric("Yerel verimlilik değişimi", f"%{_metric(report, 'Mean_Local_Efficiency_Change_Pct'):.3f}")
    metrics[3].metric("Sistemik ağ kayması", f"%{_metric(report, 'Systemic_Network_Shift_Pct'):.3f}")

    overview, network, changes, biology, advanced, export = st.tabs(
        ["Özet", "Network", "Redistribution", "Biyolojik bağlam", "Advanced", "Export"]
    )
    with overview:
        render_section_header("Summary", "Perturbation sonucu ve baseline topology alanları görsel olarak ayrıdır.", label="ANALYSIS")
        columns = [c for c in ["Gen Sembolü", "Protein Kimliği", "Delta_PageRank", "q_value", "BC_Skoru", "Hinterland_Skoru", "Community"] if c in display]
        render_dataframe(symbol_first_report(critical.head(6), mapping)[columns], label="kritik_degisim_ozeti", key="critical_summary", row_limit=6) if not critical.empty else st.info("Kritik değişim listesi boş.")
        render_scientific_note("Ağ önemi değişimi pertürbasyon sonucudur. Geçiş Merkeziliği, Ağdaki Önemi ve topluluk hedef zayıflatması öncesindeki başlangıç topolojisi bağlamıdır.")
    with network:
        render_network_views(graph, report, targets, mapping)
    with changes:
        render_section_header("Ağ yeniden dağılımı", "Hedef zayıflatması sonrası göreli ağ önemi değişimi.", label="PERTÜRBASYON SONUCU")
        st.markdown(f"### {len(critical):,} redistribution candidate")
        render_dataframe(symbol_first_report(critical, mapping), label="kritik_degisimler", key="critical_changes", row_limit=500)
        from src.interpretation import select_redistribution_losses
        signed = st.session_state.get("signed_redistribution_result")
        limit = int(st.session_state.get("negative_redistribution_top_n", 100))
        threshold = float(st.session_state.get("negative_min_abs_delta_pct", .05))
        losses = select_redistribution_losses(signed, min_abs_delta_pct=threshold, limit=limit)
        render_section_header("Ağ Önemi Kaybı", "Pertürbasyon sonrasında göreli yapısal önem kaybı; ekspresyon veya aktivite kaybı değildir.", label="İKİ YÖNLÜ PERTÜRBASYON SONUCU")
        if losses.empty:
            st.info("Bu oturumda signed iki-yönlü sonuç erişilebilir değil veya seçili eşikte negatif aday yok.")
        else:
            render_dataframe(symbol_first_report(losses, mapping), label="redistribution_losses", key="redistribution_losses", row_limit=max(limit, 100))
    with biology:
        from src.core.analysis_runtime import directed_evidence_for_selection
        from src.interpretation import interpret_forward_results
        from src.ui.biological_interpretation import render_biological_interpretation
        directed, directed_metadata = directed_evidence_for_selection(is_mouse)
        st.caption(
            f"Yapısal katman: STRING yönsüz PPI · TRRUST: {directed_metadata.get('trrust_records', 0):,} kayıt · "
            f"OmniPath: {'kullanılabilir' if directed_metadata.get('available') else 'snapshot kullanılamıyor'}"
        )
        enrichment = st.session_state.get("enrichment_results")
        interpretation_candidates = critical.copy(deep=True)
        if "gene" in interpretation_candidates:
            interpretation_candidates["Symbol"] = interpretation_candidates["gene"].astype(str).map(mapping).fillna("")
        interpretation = interpret_forward_results(
            report=report, candidates=interpretation_candidates, targets=targets, gene_to_symbol=mapping,
            tissue=tissue, organism="Mus musculus" if is_mouse else "Homo sapiens",
            edge_attenuation=float(st.session_state.get("block_strength", .001)),
            network_nodes=graph.vcount(), network_edges=graph.ecount(),
            directed_evidence=directed,
            enrichment=enrichment if isinstance(enrichment, pd.DataFrame) else None,
            loss_candidates=select_redistribution_losses(
                st.session_state.get("signed_redistribution_result"),
                min_abs_delta_pct=float(st.session_state.get("negative_min_abs_delta_pct", .05)),
                limit=int(st.session_state.get("negative_redistribution_top_n", 100)),
            ),
        )
        render_biological_interpretation(interpretation)
        for notice in st.session_state.get("enrichment_notices", []):
            st.caption(notice)
    with advanced:
        dose_tab, sweep_tab, propagation_tab, sensitivity_tab = st.tabs(["Doz-Yanıt", "Threshold Sweep", "Propagation Trace", "Hassasiyet"])
        with dose_tab:
            if not targets:
                st.info("Doz-yanıt için gen hedefiyle bir simülasyon çalıştırın.")
            elif st.button("Doz-yanıt simülasyonunu başlat", key="run_dose_response"):
                with st.status("Doz-yanıt senaryoları hesaplanıyor…") as status:
                    try:
                        st.session_state["dose_df"] = runtime.motor_module.run_pharmacological_dose_response(
                            graph, scores, runtime.output_dir, spesifik_hedefler=targets,
                            hill_n=float(st.session_state.get("hill_n", 2.0)),
                            paralog_boost=float(st.session_state.get("paralog_boost", 10)) / 100.0,
                            damping=float(st.session_state.get("damping", .85)),
                        )
                        status.update(label="Doz-yanıt tamamlandı", state="complete")
                    except Exception as exc:
                        log.exception("Doz-yanıt çalıştırılamadı")
                        status.update(label="Doz-yanıt tamamlanamadı", state="error")
                        st.error("Doz-yanıt analizi bu koşulda tamamlanamadı.")
            dose = st.session_state.get("dose_df")
            if isinstance(dose, pd.DataFrame) and not dose.empty:
                if {"Survival_Fraction", "Signal_Kayip_Pct"} <= set(dose):
                    _dose_chart = dose.set_index("Survival_Fraction")[["Signal_Kayip_Pct"]].copy()
                    if locale() == "en":
                        _dose_chart = _dose_chart.rename(columns={"Signal_Kayip_Pct": "Signal loss (%)"})
                    st.line_chart(_dose_chart, color=ACCENT)
                render_dataframe(dose, label="doz_yanit", key="dose_table")
        with sweep_tab:
            if is_mouse:
                st.info("Threshold Sweep fare motorunda desteklenmiyor.")
            elif st.button("Threshold Sweep çalıştır", key="run_threshold_sweep"):
                with st.status("Ağ eşikleri taranıyor…") as status:
                    try:
                        runtime.motor_module.run_threshold_sweep(str(runtime.motor_module.DB_PATH))
                        path = runtime.output_dir / "reports" / "threshold_sweep_raporu.csv"
                        st.session_state["sweep_df"] = pd.read_csv(path) if path.exists() else pd.DataFrame()
                        status.update(label="Threshold Sweep tamamlandı", state="complete")
                    except Exception as exc:
                        log.exception("Threshold Sweep tamamlanamadı")
                        status.update(label="Threshold Sweep tamamlanamadı", state="error")
            sweep = st.session_state.get("sweep_df")
            if isinstance(sweep, pd.DataFrame) and not sweep.empty:
                chart_cols = [c for c in ["Resilience", "LCC_Orani"] if c in sweep]
                _sweep_chart = sweep.set_index("Esik")[chart_cols].copy()
                if locale() == "en":
                    _sweep_chart = _sweep_chart.rename(columns={"Resilience": "Network resilience", "LCC_Orani": "LCC ratio"})
                st.line_chart(_sweep_chart, color=[ACCENT, TIER_DUSUK][:len(chart_cols)])
                render_dataframe(sweep, label="threshold_sweep", key="sweep_table")
        with propagation_tab:
            from src.services.propagation_trace import build_propagation_trace, export_propagation_trace
            from src.ui.pages.propagation_trace import render_propagation_trace
            if not targets:
                st.info("Propagation Trace için en az bir gen hedefi gerekir.")
            elif st.button("Etkinin yayılımını incele", key="run_propagation"):
                if is_mouse:
                    from src.mouse.config import OMNIPATH_SNAPSHOT_PATH
                else:
                    from src.config import OMNIPATH_SNAPSHOT_PATH
                from src.core.analysis_runtime import regulators_for_selection
                st.session_state["propagation_trace"] = build_propagation_trace(
                    graph=graph, report=report, target_gene=targets[0], baseline_scores=scores,
                    gene_to_symbol=mapping, regulators=regulators_for_selection(is_mouse),
                    omnipath_path=OMNIPATH_SNAPSHOT_PATH, max_nodes=int(st.session_state.get("propagation_max_nodes", 40)),
                    max_routes=int(st.session_state.get("propagation_max_routes", 12)), max_hops=3,
                    trace_context={"organism":"Mus musculus" if is_mouse else "Homo sapiens","tissue":tissue},
                )
            trace = st.session_state.get("propagation_trace")
            if trace is not None:
                render_propagation_trace(trace, target_label=mapping.get(targets[0], targets[0]), tissue=tissue, species="Fare" if is_mouse else "İnsan")
                if st.button("Propagation Trace dosyalarını hazırla", key="export_trace"):
                    folder = export_propagation_trace(trace, runtime.output_dir / "propagation_trace" / targets[0])
                    st.success(f"Dosyalar hazır: {folder}")
        with sensitivity_tab:
            st.caption("Perturbasyon şiddeti ve çoklu hedef etkileşimi sonuçları burada ayrı ve açıklanabilir metriklerle doğrulanır.")
            if targets and st.button("Şiddet hassasiyetini çalıştır", key="run_severity"):
                from src.services.severity_sensitivity import run_severity_sensitivity
                st.session_state["severity_sensitivity_result"] = run_severity_sensitivity(
                    motor_module=runtime.motor_module, graph=graph, scores=scores, targets=targets,
                    output_root=runtime.output_dir, damping=float(st.session_state.get("damping", .85)),
                )
            sensitivity = st.session_state.get("severity_sensitivity_result")
            if isinstance(sensitivity, pd.DataFrame):
                render_dataframe(sensitivity, label="hassasiyet_analizi", key="sensitivity_table")
            if len(targets) > 1:
                st.markdown("#### Predicted Network Interaction")
                st.caption("A, B ve A+B forward sonuçları ayrı karşılaştırılır; bu biyolojik synergy iddiası değildir.")
                tolerance = st.slider("Additive göreli tolerans", 0.0, 0.50, 0.10, 0.01, key="network_interaction_tolerance")
                if st.button("Çoklu hedef etkileşimini doğrula", key="run_network_interaction"):
                    from src.services.network_interaction import predicted_network_interaction
                    st.session_state["network_interaction_result"] = predicted_network_interaction(
                        motor_module=runtime.motor_module, graph=graph, scores=scores, targets=targets,
                        output_root=runtime.output_dir, block_strength=float(st.session_state.get("block_strength", .001)),
                        damping=float(st.session_state.get("damping", .85)), relative_tolerance=tolerance,
                    )
                interaction = st.session_state.get("network_interaction_result")
                if isinstance(interaction, pd.DataFrame):
                    render_dataframe(interaction, label="etkilesim_analizi", key="interaction_table")
    with export:
        render_section_header("Export", "Tüm tablo satırlarını standardize edilmiş CSV veya Excel olarak indirin.", label="EXPORT")
        st.markdown("### Tam analiz raporu")
        render_dataframe(display, label="sophiark_perturbasyon_raporu", key="forward_full_export", row_limit=500, paginate=True)
        with st.expander("Metodoloji ve yeniden üretilebilirlik metadata'sı"):
            metadata_cols = [column for column in report.columns if any(token in column.lower() for token in ("version", "mode", "seed", "timestamp", "snapshot", "damping", "iteration"))]
            st.json({column: str(report[column].iloc[0]) for column in metadata_cols}, expanded=False)
