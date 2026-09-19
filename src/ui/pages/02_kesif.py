"""Discovery analyses using the shared prepared graph and runtime."""

from __future__ import annotations

import logging

import pandas as pd
import streamlit as st

from src.core.analysis_runtime import (
    prepare_active_network, regulators_for_selection, run_compensation_analysis,
    run_target_stress_search, runtime_for_selection, symbol_options,
)
from src.ui.feedback import user_error
from src.ui.pages.compensation import render_compensation_results
from src.ui.pages.target_stress import render_target_stress_results

log = logging.getLogger(__name__)


def _ghost_filter(frame: pd.DataFrame, graph):
    names = set(map(str, graph.vs["name"]))
    gene = frame["gene"].astype(str)
    return frame[gene.isin(names)].copy(), frame.loc[~gene.isin(names), "gene"].astype(str).tolist()


is_mouse = st.session_state.get("ui_species", "İnsan") == "Fare"
species = "Fare" if is_mouse else "İnsan"
tissue = st.session_state.get("ui_tissue", "None")
runtime = runtime_for_selection(is_mouse=is_mouse, tissue=tissue)
symbols = symbol_options(is_mouse)
mapping = dict(zip(symbols["gene"].astype(str), symbols["Symbol"].astype(str))) if not symbols.empty else {}

st.title("Keşif")
st.caption("Bir hedefi dolaylı etkileyebilecek adayları veya perturbasyon sonrası olası telafi yanıtını inceleyin.")
kind = st.segmented_control(
    "Analiz türü", ["Target Stress Search", "Compensation Analysis · BETA"],
    default="Target Stress Search", key="discovery_kind",
) or "Target Stress Search"
is_target_search = kind == "Target Stress Search"

options = [f"{row.Symbol} · {row.gene}" for _, row in symbols.iterrows()]
chosen = st.selectbox("Hedef gen", options, index=None, placeholder="Gen sembolü veya protein kimliğiyle arayın", key=f"discovery_target_{kind}")
if is_target_search:
    influence_mode = st.segmented_control(
        "Etki modu", ["Doğrudan Etki", "Ağ Aracılı Etki", "Birleşik"],
        default="Ağ Aracılı Etki", key="target_influence_mode",
    ) or "Ağ Aracılı Etki"
    candidate_limit = st.slider("Doğrulanacak aday sayısı", 3, 25, 5, key="target_candidate_limit")
    with st.expander("Gelişmiş doğrulama", expanded=False):
        deep_validation = st.selectbox(
            "Full-BC deep validation", ["Kapalı", "İlk 3 aday", "İlk 5 aday"],
            help="BETA · Büyük ağlarda aday başına birkaç dakika sürebilir.",
        )
else:
    influence_mode, candidate_limit, deep_validation = "Ağ Aracılı Etki", 5, "Kapalı"
    st.info("BETA · Sonuçlar hesaplamalı telafi adaylarıdır; deneysel kompanzasyon kanıtı değildir.")

if st.button("Analizi çalıştır", type="primary", disabled=chosen is None, key="run_discovery"):
    target = chosen.rsplit(" · ", 1)[-1]
    try:
        with st.status("Aktif ağ ve baseline metrikleri hazırlanıyor…", expanded=True) as status:
            expected_context = {"species": species, "tissue": tissue}
            if st.session_state.get("workspace_graph") is not None and st.session_state.get("workspace_context") == expected_context:
                graph, scores = st.session_state["workspace_graph"], st.session_state["workspace_scores"]
                status.write("Mevcut doku grafı ve baseline skorları yeniden kullanıldı.")
            else:
                prepared = prepare_active_network(runtime, forced_genes=(target,), bc_sample_sources=st.session_state.get("bc_sample_sources", 1200))
                graph, scores = prepared.graph, prepared.scores
                st.session_state.update(workspace_graph=graph, workspace_scores=scores, workspace_context=expected_context)
                status.write(f"Ağ hazır: {graph.vcount():,} düğüm / {graph.ecount():,} kenar")
            if is_target_search:
                status.write("Adaylar ön seçimden geçiriliyor ve ayrı perturbasyonlarla doğrulanıyor…")
                deep_n = {"Kapalı": 0, "İlk 3 aday": 3, "İlk 5 aday": 5}[deep_validation]
                result = run_target_stress_search(
                    context=runtime, graph=graph, scores=scores, target_gene=target,
                    regulators=regulators_for_selection(is_mouse), gene_to_symbol=mapping,
                    mode=influence_mode, limit=candidate_limit,
                    block_strength=float(st.session_state.get("block_strength", .001)),
                    damping=float(st.session_state.get("damping", .85)), ghost_filter=_ghost_filter,
                    deep_validation_top_n=deep_n,
                )
                st.session_state["target_stress_results"] = result
                st.session_state["target_stress_context"] = {"target": target, "label": chosen.split(" · ")[0], "species": species, "tissue": tissue}
            else:
                status.write("Perturbasyon sonrası rol kazanan adaylar hesaplanıyor…")
                result = run_compensation_analysis(
                    context=runtime, graph=graph, scores=scores, perturbed_gene=target,
                    block_strength=float(st.session_state.get("block_strength", .001)),
                    damping=float(st.session_state.get("damping", .85)), ghost_filter=_ghost_filter,
                )
                st.session_state["compensation_results"] = result
                st.session_state["compensation_context"] = {"target": target, "label": chosen.split(" · ")[0], "species": species, "tissue": tissue}
            status.update(label="Keşif analizi tamamlandı", state="complete", expanded=False)
    except Exception as exc:
        user_error("Keşif analizi bu koşulda tamamlanamadı. Hedefi ve aktif doku bağlamını kontrol edin.", exception=exc)

if is_target_search:
    result = st.session_state.get("target_stress_results")
    context = st.session_state.get("target_stress_context") or {}
    if isinstance(result, pd.DataFrame) and context.get("species") == species and context.get("tissue") == tissue:
        render_target_stress_results(result, target_label=context.get("label", context.get("target", "")), tissue=tissue, species=species)
    else:
        st.info("Henüz Target Stress Search çalıştırmadınız. TP53 veya CFTR ile başlayabilirsiniz.")
else:
    result = st.session_state.get("compensation_results")
    context = st.session_state.get("compensation_context") or {}
    if isinstance(result, pd.DataFrame) and context.get("species") == species and context.get("tissue") == tissue:
        render_compensation_results(
            result, target_label=context.get("label", context.get("target", "")),
            tissue=tissue, species=species, gene_to_symbol=mapping,
        )
    else:
        st.info("Henüz Compensation Analysis çalıştırmadınız. Önce bir hedef seçin.")
