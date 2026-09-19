import streamlit as st

from src.ui.tables import render_dataframe
from src.scientific.perturbation import NULL_MODEL_PRESETS, SEVERITY_FRACTIONS

st.title("Ayarlar")
st.caption("Bilimsel senaryo parametreleri, görselleştirme sınırları ve geliştirici kontrolleri.")
with st.expander("Gelişmiş simülasyon ayarları", expanded=False):
    st.session_state["damping"] = st.slider("Ağ önemi yayılımı", .5, .99, float(st.session_state.get("damping", .85)))
    severity = st.selectbox(
        "Perturbasyon şiddeti / edge attenuation senaryosu",
        list(SEVERITY_FRACTIONS),
        index=list(SEVERITY_FRACTIONS).index(st.session_state.get("perturbation_severity", "Near-complete")),
    )
    st.session_state["perturbation_severity"] = severity
    st.session_state["block_strength"] = SEVERITY_FRACTIONS[severity]
    st.caption(f"Kalan edge strength fraksiyonu: {SEVERITY_FRACTIONS[severity]:.3f}")
    st.session_state["bc_sample_sources"] = st.select_slider("Betweenness örnek kaynağı", [300, 600, 1200], value=st.session_state.get("bc_sample_sources", 1200))
    st.session_state["hill_n"] = st.slider("Doz-yanıt Hill katsayısı", 0.5, 5.0, float(st.session_state.get("hill_n", 2.0)), 0.1)
    st.session_state["paralog_boost"] = st.slider("Paralog desteği (%)", 0, 100, int(st.session_state.get("paralog_boost", 10)), 5)
    st.session_state["efficiency_validation_mode"] = st.selectbox(
        "Global/local efficiency modu", ["Fast Approximate", "Full Validation"],
        index=["Fast Approximate", "Full Validation"].index(st.session_state.get("efficiency_validation_mode", "Fast Approximate")),
        help="Full Validation tüm kaynak/düğümleri kullanır ve büyük ağda çok uzun sürebilir.",
    )
    st.session_state["efficiency_sample_sources"] = st.select_slider("Global efficiency örnek kaynağı", [64, 128, 256, 512], value=int(st.session_state.get("efficiency_sample_sources", 256)))
    st.session_state["local_efficiency_sample_nodes"] = st.select_slider("Local efficiency örnek düğümü", [32, 64, 128, 256], value=int(st.session_state.get("local_efficiency_sample_nodes", 128)))
    st.session_state["redistribution_mode"] = st.selectbox(
        "Redistribution seçim modu", ["null_fdr", "percentile", "top_n"],
        index=["null_fdr", "percentile", "top_n"].index(st.session_state.get("redistribution_mode", "null_fdr")),
        help="Varsayılan null_fdr degree-matched empirical null ve BH-FDR kullanır.",
    )
    preset = st.selectbox("Null model hassasiyeti", list(NULL_MODEL_PRESETS), index=1)
    st.session_state["null_iterations"] = NULL_MODEL_PRESETS[preset]
    st.session_state["redistribution_top_n"] = st.number_input("Keşifsel Top-N", 10, 1000, int(st.session_state.get("redistribution_top_n", 250)), 10)
    st.session_state["redistribution_percentile"] = st.slider("|Ağ önemi değişimi| persentili", 90.0, 100.0, float(st.session_state.get("redistribution_percentile", 99.0)), .1)
    st.caption("Negatif yeniden dağılım görünümü yalnız sonuç seçimini etkiler; ağ önemi hesabını veya pozitif tabloyu değiştirmez.")
    st.session_state["negative_min_abs_delta_pct"] = st.number_input("Minimum ağ önemi kaybı (%)", 0.0, 100.0, float(st.session_state.get("negative_min_abs_delta_pct", .05)), .01)
    st.session_state["negative_redistribution_top_n"] = st.number_input("Gösterilecek negatif aday sayısı", 10, 2000, int(st.session_state.get("negative_redistribution_top_n", 100)), 10)
    st.session_state["tissue_normalization_mode"] = st.selectbox(
        "Doku normalizasyonu", ["within_tissue", "cross_tissue_comparable"],
        index=["within_tissue", "cross_tissue_comparable"].index(st.session_state.get("tissue_normalization_mode", "within_tissue")),
        help="Dokular Arası Karşılaştırma her zaman ortak tür referansını kullanır.",
    )
    st.session_state["edge_evidence_weighting_mode"] = st.selectbox(
        "Evidence → edge weighting modu", ["structural_only", "legacy_edge_modifiers"],
        index=["structural_only", "legacy_edge_modifiers"].index(st.session_state.get("edge_evidence_weighting_mode", "structural_only")),
        help="structural_only: STRING omurgası değişmez; TRRUST/OmniPath ayrı evidence overlay olarak kalır.",
    )
    st.session_state["scientific_random_seed"] = st.number_input("Random seed", 0, 2_147_483_647, int(st.session_state.get("scientific_random_seed", 42)))
with st.expander("Ağ haritası ayarları", expanded=False):
    st.session_state["harita_max_node"] = st.slider("Maksimum görünür düğüm", 50, 800, int(st.session_state.get("harita_max_node", 350)), 50)
    st.session_state["propagation_max_nodes"] = st.slider("Propagation Trace düğümü", 20, 100, int(st.session_state.get("propagation_max_nodes", 40)), 5)
    st.session_state["propagation_max_routes"] = st.slider("Propagation Trace rota sayısı", 3, 30, int(st.session_state.get("propagation_max_routes", 12)), 1)
st.subheader("İleri bilimsel doğrulama")
st.caption("Bu kontroller analizi değiştirmez; hassasiyet ve yaklaşık hesap kararlılığını ölçer.")
graph = st.session_state.get("workspace_graph")
if graph is not None:
    if st.button("Ağ önemi hassasiyet panelini çalıştır", type="secondary"):
        from src.scientific.validation import pagerank_robustness
        st.session_state["pagerank_robustness"] = pagerank_robustness(graph)
    robustness = st.session_state.get("pagerank_robustness")
    if robustness:
        st.write(robustness[1])
        render_dataframe(robustness[0], hide_index=True, width="stretch")
    st.warning("Full BC doğrulaması büyük ağlarda uzun sürebilir; yalnızca ileri validasyonda çalıştırın.")
    if st.button("Full BC deep validation çalıştır (BETA)", type="secondary"):
        from src.scientific.validation import deep_bc_validation
        st.session_state["bc_deep_validation"] = deep_bc_validation(
            graph, sample_sources=int(st.session_state.get("bc_sample_sources", 1200))
        )
    if st.session_state.get("bc_deep_validation"):
        st.json(st.session_state["bc_deep_validation"])

developer = st.toggle("Geliştirici modu", value=bool(st.session_state.get("ui_developer_mode", False)), key="ui_developer_mode")
if developer:
    st.markdown("### Simülasyon geçmişi")
    if st.session_state.get("analysis_snapshots"):
        render_dataframe(st.session_state["analysis_snapshots"], hide_index=True, width="stretch")
    else:
        st.info("Bu oturumda kayıtlı analiz yok.")
    st.write("Aktif session anahtarları", sorted(st.session_state.keys()))
