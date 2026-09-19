"""Target Stress Search result view; it does not mutate simulation state."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from src.ui.guides import render_mode_guide
from src.ui.insight_narrator import en_onemli_bulguyu_ozetle
from src.ui.workspace import render_insight_box
from src.ui.components import render_scientific_note, render_section_header
from src.ui.tables import render_dataframe


def render_target_stress_results(results: pd.DataFrame, *, target_label: str, tissue: str, species: str) -> bool:
    st.markdown("## Target Stress Search")
    st.caption(
        f"Hedef: {target_label} · {species} · {tissue}. Her satır ayrı bir forward perturbation "
        "ile doğrulandı. Ağ-aracılı tahmin, doğrudan biyolojik düzenleme kanıtı değildir."
    )
    st.markdown("Perturbation sensitivity analysis")
    render_mode_guide("Target Stress Search")
    if results is None or results.empty:
        st.warning("Aktif doku ağında seçilen kanıt filtresini geçen aday bulunamadı.")
        return st.button("Analize dön", key="target_stress_close_empty")
    insight, notable = en_onemli_bulguyu_ozetle(results, mode="Target Stress Search")
    render_insight_box(insight, is_notable=notable)
    display_columns = [
        "Aday", "Öngörülen Hedef Stresi", "ΔPageRank (%)", "Δ Yerel BC (2-hop projeksiyonu)", "Ağ Mesafesi",
        "Full ΔBC", "BC Doğrulama Modu", "Yönlü Kanıt", "Doku Uygunluğu", "Güven / Kanıt Gücü",
    ]
    available = [column for column in display_columns if column in results.columns]
    render_section_header("Candidate ranking", "Her aday ayrı forward perturbation ile değerlendirilir.", label="PERTURBATION SENSITIVITY")
    render_dataframe(
        results[available], width="stretch", hide_index=True,
        column_config={
            "Öngörülen Hedef Stresi": st.column_config.NumberColumn(format="%.4f"),
            "ΔPageRank (%)": st.column_config.NumberColumn("Ağ önemi değişimi (%)", format="%+.4f"),
            "Δ Yerel BC (2-hop projeksiyonu)": st.column_config.NumberColumn(format="%+.4f"),
            "Full ΔBC": st.column_config.NumberColumn(format="%+.4f"),
            "Ağ Mesafesi": st.column_config.NumberColumn(format="%.3f"),
        },
    )
    render_scientific_note("Ağ önemi değişimi pertürbasyon sonucudur. Local ΔBC ve Full ΔBC ağ-topolojisi ölçümleridir; doğrudan düzenleyici aktivite kanıtı değildir.")
    st.markdown("### Aday kanıtı")
    for _, row in results.iterrows():
        with st.expander(f"{row['Aday']} · öngörülen hedef stresi {float(row['Öngörülen Hedef Stresi']):.4f}"):
            st.write(row["Mekanizma / Açıklama"])
            st.caption(
                f"Sistemik kayma: {float(row['Sistemik Kayma (%)']):.4f}% · "
                f"En güçlü pozitif ağ önemi yeniden-dağılımı: {float(row['En Güçlü Pozitif PageRank Yeniden-Dağılımı (%)']):.4f}% · "
                f"Pareto katmanı: {int(row['Pareto Katmanı'])}"
            )
            if row.get("BC Doğrulama Modu") == "Full Validation":
                st.info("Bu adayda 2-hop projeksiyonuna ek olarak full weighted BC doğrulaması çalıştırıldı.")
            else:
                st.info("BC değeri açıkça 2-hop ağ projeksiyonudur; tam-ağ congestion metriği değildir.")
    with st.expander("Export", expanded=False):
        render_dataframe(results, label="target_stress_candidates", key="target_stress_export", row_limit=500)
    return st.button("Analize dön", key="target_stress_close")
