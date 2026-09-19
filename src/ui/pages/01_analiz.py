"""Forward simulation and tissue-comparison workspace.

NOT_SUPPORTED_ENTRYPOINT: launch the supported product with
``streamlit run app.py``.  This source file is an internal workspace page,
not a standalone Streamlit application contract.
"""

from __future__ import annotations

from datetime import datetime, timezone
import logging

import pandas as pd
import streamlit as st

from src.state import invalidate_live_unified_result
from src.core.analysis_runtime import prepare_active_network, run_forward_simulation, runtime_for_selection, symbol_options
from src.metrics import compute_druggable_gate_score
from src.services.enrichment import run_post_simulation_enrichment
from src.services.report_context import critical_rows, symbol_map
from src.ui.analysis_results import render_forward_results
from src.ui.feedback import user_error
from src.ui.guides import render_mode_guide

log = logging.getLogger(__name__)


def _ghost_filter(report: pd.DataFrame, graph):
    names = set(map(str, graph.vs["name"]))
    gene = report["gene"].astype(str)
    return report[gene.isin(names)].copy(), report.loc[~gene.isin(names), "gene"].astype(str).tolist()


is_mouse = st.session_state.get("ui_species", "İnsan") == "Fare"
tissue = st.session_state.get("ui_tissue", "None")
runtime = runtime_for_selection(is_mouse=is_mouse, tissue=tissue)
symbols = symbol_options(is_mouse)
mapping = symbol_map(symbols)
reverse_mapping = {str(symbol).upper(): str(gene) for gene, symbol in mapping.items()}

st.title("Analiz")
st.caption("Bir müdahaleyi simüle edin veya aynı perturbasyonu dokular arasında karşılaştırın.")
analysis_kind = st.segmented_control(
    "Analiz türü", ["Unified Research Report", "Forward Simulation", "Dokular Arası Karşılaştırma · BETA"],
    default="Unified Research Report", key="analysis_kind",
) or "Unified Research Report"

if analysis_kind == "Unified Research Report":
    st.caption("Classic + Directed modelleri ayrı çalışır; Evidence yalnızca kanıt bileşimi, CORUM yalnızca complex bağlamı sağlar.")
    options = [f"{row.Symbol} · {row.gene}" for _, row in symbols.iterrows()]
    chosen = st.selectbox("Hedef gen", options, index=None, placeholder="Bir gen seçin", key="unified_target")
    if is_mouse:
        st.info("Unified Research Report şu an insan Directed BETA altyapısını gerektirir.")
    elif st.button("Unified analizi başlat", type="primary", disabled=not chosen, key="run_unified"):
        invalidate_live_unified_result()
        from src.unified import export_unified_report, run_unified_analysis
        target = chosen.rsplit(" · ", 1)[-1]
        with st.status("Classic ve Directed ayrı graph kopyalarında çalışıyor…", expanded=True) as status:
            report = run_unified_analysis(
                target=target, tissue=tissue,
                block_strength=float(st.session_state.get("block_strength", .001)),
                damping=float(st.session_state.get("damping", .85)),
                bc_sample_sources=st.session_state.get("bc_sample_sources", 1200),
            )
            paths = export_unified_report(report, runtime.output_dir / "unified_reports" / target / tissue.replace("/", "_"))
            st.session_state["unified_research_report"] = report
            st.session_state["unified_research_context"] = {"target": target, "tissue": tissue}
            status.write(f"Export hazır: {paths['report']}")
            if report.status.value != "FAILED":
                try:
                    from src.product.history import HistoryService
                    history_entry = HistoryService().save(report)
                    st.session_state["unified_history_entry"] = history_entry
                    st.session_state["unified_history_error"] = None
                    status.write(f"Geçmiş bundle hazır: {history_entry.bundle_filename}")
                except Exception as history_error:
                    log.exception("Unified analysis completed but history save failed")
                    st.session_state["unified_history_entry"] = None
                    st.session_state["unified_history_error"] = str(history_error)
                    status.write("Analiz tamamlandı; geçmiş kaydı oluşturulamadı.")
            status.update(label=f"Unified rapor: {report.status.value}", state="complete" if report.status.value == "COMPLETE" else "error", expanded=False)
    unified = st.session_state.get("unified_research_report")
    unified_context = st.session_state.get("unified_research_context")
    if unified is not None and unified_context == {"target": (chosen.rsplit(" · ", 1)[-1] if chosen else None), "tissue": tissue}:
        from src.ui.unified_research_report import render_unified_research_report
        render_unified_research_report(unified)
    else:
        st.info("Bir target ve tissue seçerek Classic + Directed sentez raporunu başlatın.")

elif analysis_kind == "Forward Simulation":
    render_mode_guide("Perturbasyon Analizi")
    setup_mode = st.segmented_control(
        "Hedef seçim yöntemi", ["Gen", "Hastalık", "Hücresel lokalizasyon"],
        default="Gen", key="forward_setup_mode",
    ) or "Gen"
    targets: list[str] = []
    localization = "Mitochondrion"
    options = [f"{row.Symbol} · {row.gene}" for _, row in symbols.iterrows()]

    if setup_mode == "Gen":
        chosen = st.multiselect(
            "Hedef genler", options, max_selections=10, key="forward_targets",
            placeholder="Gen sembolü veya protein kimliğiyle arayın",
        )
        targets = [item.rsplit(" · ", 1)[-1] for item in chosen]
        st.caption("Tek veya çoklu hedef seçebilirsiniz. Harita ve tablolarda gen sembolü birincil gösterilir.")
    elif setup_mode == "Hastalık":
        if is_mouse:
            st.info("Open Targets hastalık hedefleri insan sembol uzayındadır. Fare için gen hedefi seçimini kullanın.")
        else:
            disease_term = st.text_input("Hastalık ara", placeholder="Örn. cystic fibrosis, Alzheimer", key="disease_query")
            if st.button("Hastalıkları bul", key="search_disease", disabled=not disease_term):
                try:
                    from src.services.opentargets import search_diseases
                    st.session_state["disease_results"] = search_diseases(disease_term)
                except Exception as exc:
                    log.warning("Open Targets araması başarısız: %s", exc)
                    st.session_state["disease_results"] = []
                    st.warning("Hastalık kaynağına şu anda erişilemiyor; gen hedefini manuel seçebilirsiniz.")
            diseases = st.session_state.get("disease_results", [])
            disease_labels = {f"{item['name']} · {item['id']}": item for item in diseases}
            selected_disease = st.selectbox("Hastalık", list(disease_labels), index=None, placeholder="Bir hastalık seçin")
            if selected_disease and st.button("Hastalık genlerini getir", key="load_disease_genes"):
                try:
                    from src.services.opentargets import disease_genes
                    found_symbols = disease_genes(disease_labels[selected_disease]["id"])
                    st.session_state["disease_targets"] = [reverse_mapping[s.upper()] for s in found_symbols if s.upper() in reverse_mapping]
                except Exception as exc:
                    log.warning("Open Targets hedefleri alınamadı: %s", exc)
                    st.warning("Hastalık hedefleri şu anda alınamadı.")
            disease_targets = st.session_state.get("disease_targets", [])
            disease_options = [f"{mapping.get(gene, 'Alias bulunamadı')} · {gene}" for gene in disease_targets]
            selected = st.multiselect("Simüle edilecek hastalık hedefleri", disease_options, default=disease_options[:10], max_selections=25)
            targets = [item.rsplit(" · ", 1)[-1] for item in selected]
    else:
        localization = st.selectbox(
            "Hücresel lokalizasyon",
            ["Mitochondrion", "Nucleus", "Cell_Membrane", "Endoplasmic_Reticulum", "Golgi_Apparatus", "Lysosome", "Cytoplasm", "Peroxisome"],
        )
        st.caption("Seçilen bölgedeki uygun proteinler mevcut motorun lokalizasyon perturbasyonu ile birlikte baskılanır.")

    with st.expander("Bu analizde kullanılan senaryo", expanded=False):
        st.write(f"Tür: **{'Fare' if is_mouse else 'İnsan'}** · Doku: **{tissue}**")
        st.write(f"Kenar zayıflatması: **{float(st.session_state.get('block_strength', .001)):.3f}** · Ağ önemi yayılımı: **{float(st.session_state.get('damping', .85)):.2f}**")
        st.caption("Tüm gelişmiş parametreler Ayarlar sayfasında açıkça görülebilir.")

    can_run = bool(targets) or setup_mode == "Hücresel lokalizasyon"
    if st.button("Simülasyonu başlat", type="primary", disabled=not can_run, key="run_forward"):
        try:
            with st.status("Aktif doku ağı hazırlanıyor…", expanded=True) as status:
                prepared = prepare_active_network(
                    runtime, forced_genes=tuple(targets) if targets else None,
                    bc_sample_sources=st.session_state.get("bc_sample_sources", 1200),
                )
                status.write(f"Ağ hazır: {prepared.graph.vcount():,} düğüm / {prepared.graph.ecount():,} kenar")
                status.write("Perturbation response hesaplanıyor…")
                report, ghosts = run_forward_simulation(
                    context=runtime, graph=prepared.graph, scores=prepared.scores, output_dir=runtime.output_dir,
                    targets=targets or None, localization="SNIPER_MODU" if targets else localization,
                    block_strength=float(st.session_state.get("block_strength", .001)),
                    damping=float(st.session_state.get("damping", .85)),
                    apply_druggability=(lambda frame: frame if is_mouse else compute_druggable_gate_score(frame, target_type="stressed")),
                    ghost_filter=_ghost_filter,
                )
                st.session_state.update(
                    workspace_report=report, workspace_graph=prepared.graph, workspace_scores=prepared.scores,
                    workspace_targets=list(targets), workspace_context={"species":"Fare" if is_mouse else "İnsan", "tissue":tissue},
                )
                species_key = "mouse" if is_mouse else "human"
                st.session_state["reports_by_species"][species_key] = report.copy()
                st.session_state["analysis_snapshots"].append({
                    "timestamp": datetime.now(timezone.utc).isoformat(), "species": species_key, "tissue": tissue,
                    "targets": list(targets), "rows": len(report), "critical": len(critical_rows(report)),
                })
                status.write("KEGG ve GO zenginleştirme analizi hazırlanıyor…")
                affected = [mapping.get(str(gene), "") for gene in critical_rows(report).get("gene", [])]
                universe = [mapping.get(str(gene), "") for gene in prepared.graph.vs["name"]]
                enrichment, notices = run_post_simulation_enrichment(affected, is_mouse=is_mouse, background_symbols=universe)
                st.session_state["enrichment_results"] = enrichment
                st.session_state["enrichment_notices"] = notices
                status.update(label="Analiz tamamlandı", state="complete", expanded=False)
        except Exception as exc:
            user_error("Analiz bu koşulda tamamlanamadı. Seçimleri kontrol edip yeniden deneyin.", exception=exc)

    report = st.session_state.get("workspace_report")
    context = st.session_state.get("workspace_context")
    if isinstance(report, pd.DataFrame) and not report.empty and context == {"species":"Fare" if is_mouse else "İnsan", "tissue":tissue}:
        render_forward_results(
            report=report, graph=st.session_state["workspace_graph"], scores=st.session_state["workspace_scores"],
            targets=st.session_state.get("workspace_targets", []), symbols=symbols, runtime=runtime,
            is_mouse=is_mouse, tissue=tissue,
        )
    else:
        st.info("Henüz bu tür ve doku bağlamında bir simülasyon çalıştırmadınız. CFTR, GLP1R veya PSEN1 ile başlayabilirsiniz.")

else:
    st.subheader("Dokular Arası Karşılaştırma · BETA")
    options = [f"{row.Symbol} · {row.gene}" for _, row in symbols.iterrows()]
    chosen = st.selectbox("Hedef gen", options, index=None, placeholder="Bir hedef seçin", key="tissue_diff_target")
    from src.core.analysis_runtime import tissue_options
    available_tissues = [value for value in tissue_options(is_mouse) if value != "None"]
    tissues = st.multiselect("Karşılaştırılacak dokular", available_tissues, max_selections=8)
    if st.button("Doku karşılaştırmasını başlat", type="primary", disabled=not chosen or len(tissues) < 2):
        from src.services.tissue_differential import run_tissue_differential
        target = chosen.rsplit(" · ", 1)[-1]
        with st.status("Her doku için forward simülasyon çalıştırılıyor…", expanded=True) as status:
            result = run_tissue_differential(
                motor_module=runtime.motor_module, tissues=tissues, targets=[target], output_root=runtime.output_dir,
                bc_sample_sources=st.session_state.get("bc_sample_sources", 1200),
                block_strength=float(st.session_state.get("block_strength", .001)), damping=float(st.session_state.get("damping", .85)),
                ghost_filter=_ghost_filter,
                apply_druggability=(lambda frame: frame if is_mouse else compute_druggable_gate_score(frame, target_type="stressed")),
            )
            st.session_state["tissue_differential_result"] = result
            st.session_state["tissue_differential_context"] = {"target": target, "species":"Fare" if is_mouse else "İnsan"}
            status.update(label="Doku karşılaştırması tamamlandı", state="complete")
    result = st.session_state.get("tissue_differential_result")
    if result is not None:
        from src.ui.pages.tissue_differential import render_tissue_differential
        context = st.session_state.get("tissue_differential_context", {})
        render_tissue_differential(result, target_label=mapping.get(context.get("target", ""), context.get("target", "")), species=context.get("species", ""))
    else:
        st.info("En az iki doku seçin. Bu analiz dokuları ortak referans normalizasyonuyla ayrı ayrı simüle eder.")
