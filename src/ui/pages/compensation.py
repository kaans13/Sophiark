"""Computational compensation result view."""

from __future__ import annotations

import pandas as pd
import streamlit as st
from src.ui.guides import render_mode_guide
from src.ui.insight_narrator import en_onemli_bulguyu_ozetle
from src.ui.workspace import render_insight_box
from src.ui.components import render_scientific_note, render_section_header
from src.ui.tables import render_dataframe


def render_compensation_results(
    results: pd.DataFrame, *, target_label: str, tissue: str, species: str,
    gene_to_symbol: dict[str, str] | None = None,
) -> bool:
    st.markdown("## Network Redistribution Candidates")
    st.caption("Predicted compensation · BETA")
    st.caption(
        f"Perturbasyon hedefi: {target_label} · {species} · {tissue}. "
        "Adaylar pertürbasyon sonrası pozitif ağ önemi artışına göre sıralanır."
    )
    render_mode_guide("Compensation Analysis")
    if results is None or results.empty:
        st.warning("Ölçülen ağ önemi artışına dayalı kompanzasyon adayı bulunamadı.")
        return st.button("Analize dön", key="compensation_close_empty")
    display = results.copy()
    protein_ids = display.get("Protein ID", display.get("gen", pd.Series("—", index=display.index))).astype(str)
    display["Protein ID"] = protein_ids
    fallback = display.get("Gene", display.get("Aday", protein_ids)).astype(str)
    display["Gene"] = [
        (gene_to_symbol or {}).get(protein_id) or candidate
        for protein_id, candidate in zip(protein_ids, fallback)
    ]
    display["Tissue"] = tissue
    display["Baseline BC"] = display.get("BC Skoru (baseline)", pd.Series(float("nan"), index=display.index))
    display["Compartment Bottleneck"] = display.get(
        "Compartment Bottleneck", pd.Series("—", index=display.index)
    )
    display = display.sort_values("ΔPageRank (%)", ascending=False, kind="mergesort")
    insight, notable = en_onemli_bulguyu_ozetle(display, mode="Compensation Analysis", symbol_col="Gene")
    render_insight_box(insight, is_notable=notable)
    render_section_header("Candidate ranking", "Gene perturbation sonrasında göreli ağ rolü artan adaylar.", label="PREDICTED COMPENSATION")
    primary_columns = ["Gene", "Protein ID", "ΔPageRank (%)", "Baseline BC", "Compartment Bottleneck", "Tissue"]
    render_dataframe(display[primary_columns], label="compensation_primary", key="compensation_primary_table", column_config={
        "Gene": st.column_config.TextColumn("Gen"),
        "Protein ID": st.column_config.TextColumn("ENSP"),
        "ΔPageRank (%)": st.column_config.NumberColumn("Ağ önemi değişimi (%)", format="+%.2f%%"),
        "Baseline BC": st.column_config.NumberColumn("Baseline BC", format="%.2f"),
    })
    render_scientific_note("Computational network prediction; validated biological compensation değildir. Ağ önemi değişimi pertürbasyon sonucu, Baseline BC ise pertürbasyon öncesi ağ bağlamıdır.")
    with st.expander("İkincil ağ bağlamı ve kanıt", expanded=False):
        detail_columns = [column for column in [
            "Gene", "Protein ID", "ΔPageRank (%)", "Baseline BC", "Compartment Bottleneck", "Tissue",
            "Açıklama", "Kanıt Türü", "Doku Uygunluğu",
        ] if column in display.columns]
        render_dataframe(display[detail_columns], label="compensation_detail", key="compensation_detail_table")
    with st.expander("Export", expanded=False):
        render_dataframe(display, label="network_redistribution_candidates", key="compensation_export", row_limit=500)
    return st.button("Analize dön", key="compensation_close")
