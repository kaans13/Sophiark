"""Wire existing advanced analyses into the supported desktop entry point.

Engines retain their defaults. Results are bound to the initiating context and
every scientific call requires its own explicit action.
"""
from __future__ import annotations
import logging
from pathlib import Path
import pickle
import time
import pandas as pd
import streamlit as st
from src.ui.tables import render_dataframe
from src.ui.localization import locale

log = logging.getLogger(__name__)


@st.fragment(run_every=3)
def _dose_progress():
    job = st.session_state.get("desktop_dose_job")
    if not job:
        return
    folder = Path(job["folder"])
    if (folder / "result.pkl").exists():
        with (folder / "result.pkl").open("rb") as stream:
            result = pickle.load(stream)
        st.session_state["desktop_dose"] = (job["context"], result)
        st.session_state.pop("desktop_dose_job", None)
        st.rerun()
    elif (folder / "error.json").exists() or job["process"].poll() is not None:
        st.error("Doz–yanıt işlemi tamamlanamadı. Ayrıntılar yerel iş günlüğünde.")
        st.caption(str(folder / "worker.log"))
        st.session_state.pop("desktop_dose_job", None)
    else:
        elapsed = int(time.time() - job["started"])
        st.info(f"Doz–yanıt hesaplanıyor · {elapsed // 60} dk {elapsed % 60} sn. Ağ merkeziliği adımları uzun sürebilir; diğer sonuçları inceleyebilirsiniz.")
        if st.button("Bu doz–yanıt çalışmasını durdur", key="stop_desktop_dose"):
            job["process"].terminate()
            st.session_state.pop("desktop_dose_job", None)
            st.rerun()


def _run(key, context, action):
    st.session_state.pop(key, None)
    try:
        with st.status("Mevcut hesaplama çalıştırılıyor…", expanded=True) as status:
            result = action()
            st.session_state[key] = (dict(context), result)
            status.update(label="Analiz tamamlandı", state="complete", expanded=False)
    except Exception:
        log.exception("Advanced UI workflow failed: %s", key)
        st.error("Bu analiz tamamlanamadı. Ana sonuç korunuyor; teknik ayrıntılar yerel günlükte.")


def _result(key, context):
    stored = st.session_state.get(key)
    if stored is None:
        return None
    if stored[0] != context:
        st.info("Bu aracın parametreleri değişti; önceki sonuç gizlendi. Yeniden çalıştırın.")
        return None
    return stored[1]


def render_advanced_workflows(*, motor, graph, scores, report, targets, tissue, is_mouse, output_dir, mapping, tissues, ghost_filter, context):
    dose, sensitivity, trace_tab, tissue_tab = st.tabs(["Doz–yanıt", "Şiddet hassasiyeti", "Yayılım izi", "Dokular arası"])
    with dose:
        st.markdown("### Farklı kalan aktivite değerlerinde ağ nasıl yanıt verir?")
        st.caption("Mevcut Hill/paralog doz–yanıt modeli; deneysel doz, klinik IC50 veya ilaç etkinliği ölçümü değildir. Varsayılan doz serisi korunur.")
        dose_context = {**context, "hill_n": st.session_state.get("hill_n", 2.0), "paralog_boost": st.session_state.get("paralog_boost", 10)}
        st.caption(f"Hill n = {dose_context['hill_n']} · Paralog boost = %{dose_context['paralog_boost']}")
        if st.button("Doz-yanıt simülasyonunu başlat", key="desktop_run_dose", disabled=not targets or bool(st.session_state.get("desktop_dose_job"))):
            from src.ui.dose_job import start_dose_job
            st.session_state.pop("desktop_dose", None)
            job = start_dose_job(motor_name=motor.__name__, graph=graph.copy(), scores=scores.copy(deep=True),
                output_dir=output_dir, targets=targets, hill_n=float(dose_context["hill_n"]),
                paralog_boost=float(dose_context["paralog_boost"])/100, damping=float(context["damping"]))
            job["context"] = dose_context
            st.session_state["desktop_dose_job"] = job
        if st.session_state.get("desktop_dose_job"):
            _dose_progress()
        result = _result("desktop_dose", dose_context)
        if isinstance(result, pd.DataFrame) and not result.empty:
            if {"Survival_Fraction", "Signal_Kayip_Pct"} <= set(result):
                st.line_chart(
                    result.set_index("Survival_Fraction")[["Signal_Kayip_Pct"]],
                    x_label="Remaining activity fraction" if locale() == "en" else "Kalan aktivite fraksiyonu",
                    y_label="Signal loss (%)" if locale() == "en" else "Sinyal kaybı (%)",
                )
            render_dataframe(result, label="doz_yanit", key="desktop_dose_table", paginate=True)
        elif result is None and not st.session_state.get("desktop_dose_job"):
            st.info("Doz–yanıt henüz bu bağlamda çalıştırılmadı.")
    with sensitivity:
        st.markdown("### Bulgular pertürbasyon şiddetine ne kadar duyarlı?")
        st.caption("Mevcut şiddet senaryoları ayrı simüle edilir. Sonuçlar birleştirilmez ve yeni bir aday skoru üretilmez.")
        if st.button("Şiddet hassasiyetini çalıştır", key="desktop_run_sensitivity", disabled=not targets):
            from src.services.severity_sensitivity import run_severity_sensitivity
            _run("desktop_sensitivity", context, lambda: run_severity_sensitivity(motor_module=motor,
                graph=graph, scores=scores, targets=targets, output_root=output_dir, damping=context["damping"]))
        result = _result("desktop_sensitivity", context)
        if isinstance(result, pd.DataFrame):
            render_dataframe(result, label="hassasiyet_analizi", key="desktop_sensitivity_table", paginate=True)
        if len(targets) > 1:
            tolerance = st.slider("Additive göreli tolerans", 0.0, .50, .10, .01, key="desktop_interaction_tolerance")
            interaction_context = {**context, "relative_tolerance": tolerance}
            st.caption("A, B ve A+B ağ yanıtlarının mevcut karşılaştırması; biyolojik sinerji kanıtı değildir.")
            if st.button("Çoklu hedef etkileşimini doğrula", key="desktop_run_interaction"):
                from src.services.network_interaction import predicted_network_interaction
                _run("desktop_interaction", interaction_context, lambda: predicted_network_interaction(
                    motor_module=motor, graph=graph, scores=scores, targets=targets, output_root=output_dir,
                    block_strength=context["block_strength"], damping=context["damping"], relative_tolerance=tolerance))
            result = _result("desktop_interaction", interaction_context)
            if isinstance(result, pd.DataFrame):
                render_dataframe(result, label="etkilesim_analizi", key="desktop_interaction_table")
    with trace_tab:
        st.markdown("### Etki hangi ağ rotalarında izlenebilir?")
        st.caption("Mevcut yayılım izi, ağ topolojisi ile yönlü kanıtı ayrı gösterir; biyokimyasal sinyal iletimini doğrulamaz.")
        if st.button("Etkinin yayılımını incele", key="desktop_run_trace", disabled=not targets):
            log.info("Advanced UI workflow requested: desktop_trace")
            from src.services.propagation_trace import build_propagation_trace
            from src.core.analysis_runtime import regulators_for_selection
            if is_mouse:
                from src.mouse.config import OMNIPATH_SNAPSHOT_PATH
            else:
                from src.config import OMNIPATH_SNAPSHOT_PATH
            _run("desktop_trace", context, lambda: build_propagation_trace(graph=graph, report=report,
                target_gene=targets[0], baseline_scores=scores, gene_to_symbol=mapping,
                regulators=regulators_for_selection(is_mouse), omnipath_path=OMNIPATH_SNAPSHOT_PATH,
                max_nodes=int(st.session_state.get("propagation_max_nodes", 40)),
                max_routes=int(st.session_state.get("propagation_max_routes", 12)), max_hops=3,
                trace_context={"organism": "Mus musculus" if is_mouse else "Homo sapiens", "tissue": tissue}))
        result = _result("desktop_trace", context)
        if result is not None:
            from src.ui.pages.propagation_trace import render_propagation_trace
            if render_propagation_trace(result, target_label=mapping.get(targets[0], targets[0]), tissue=tissue, species="Fare" if is_mouse else "İnsan"):
                st.session_state.pop("desktop_trace", None)
                st.rerun()
    with tissue_tab:
        st.markdown("### Aynı hedef farklı dokularda nasıl yanıt verir?")
        chosen = st.multiselect("Karşılaştırılacak dokular", [item for item in tissues if item != "None"], max_selections=8, key="desktop_tissue_choices")
        tissue_context = {**context, "tissues": tuple(chosen)}
        st.caption("Her doku mevcut ortak referans normalizasyonuyla ayrı hesaplanır. En az iki doku ve tek hedef gerekir.")
        if st.button("Doku karşılaştırmasını başlat", key="desktop_run_tissues", disabled=len(targets) != 1 or len(chosen) < 2):
            from src.services.tissue_differential import run_tissue_differential
            from src.metrics import compute_druggable_gate_score
            _run("desktop_tissues", tissue_context, lambda: run_tissue_differential(motor_module=motor,
                tissues=chosen, targets=targets, output_root=output_dir,
                bc_sample_sources=context["bc_sample_sources"], block_strength=context["block_strength"],
                damping=context["damping"], ghost_filter=ghost_filter,
                apply_druggability=(lambda frame: frame if is_mouse else compute_druggable_gate_score(frame, target_type="stressed"))))
        result = _result("desktop_tissues", tissue_context)
        if result is not None:
            from src.ui.pages.tissue_differential import render_tissue_differential
            if render_tissue_differential(result, target_label=mapping.get(targets[0], targets[0]), species="Fare" if is_mouse else "İnsan"):
                st.session_state.pop("desktop_tissues", None)
                st.rerun()
