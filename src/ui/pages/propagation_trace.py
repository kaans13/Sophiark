"""Propagation Trace sunumu; tüm değerler servis sonucundan gelir."""

from __future__ import annotations

import pandas as pd
import streamlit as st

from src.models import PropagationTrace
from src.services.propagation_trace import filter_trace_nodes
from src.ui.design_tokens import ACCENT, BG_SURFACE_2, BORDER
from src.ui.workspace import render_insight_box
from src.ui.components import render_scientific_note, render_section_header
from src.ui.tables import render_dataframe


def _narrative(trace: PropagationTrace, target: str) -> str:
    affected = trace.nodes[trace.nodes["Propagation_Layer"] > 0]
    if affected.empty:
        return f"{target} perturbasyonu sonrası yayılım izi için anlamlı ikinci bir ağ düğümü seçilemedi."
    early = affected.sort_values(
        ["Propagation_Layer", "Significant_Redistribution", "Delta_PageRank"],
        ascending=[True, False, False],
    ).head(2)
    labels = ", ".join(early["Gene_Symbol"].astype(str))
    directed = int(trace.edges.get("Directed_Evidence", pd.Series(dtype=bool)).sum()) if not trace.edges.empty else 0
    return (
        f"{target} perturbasyonu sonrası ölçülen ağ yeniden-dağılımı erken katmanlarda {labels} üzerinde öne çıkıyor. "
        f"Seçili yapısal rotalarda {directed} yönlü kanıt destekli kenar bulundu. Bu, predicted network propagation'dır; biyokimyasal sinyal iletimi kanıtı değildir."
    )


def _render_graph(trace: PropagationTrace, visible: pd.DataFrame) -> None:
    if trace.edges.empty or visible.empty:
        st.info("Bu filtre için gösterilecek rota yok.")
        return
    allowed = set(visible["Gene_Symbol"].astype(str))
    edges = trace.edges[trace.edges["Source"].isin(allowed) & trace.edges["Target"].isin(allowed)]
    if edges.empty:
        st.info("Bu filtre için gösterilecek rota yok.")
        return
    lines = ["digraph trace {", "rankdir=LR;", f'node [shape=ellipse, style=filled, fillcolor="{BG_SURFACE_2}", fontcolor="#E6ECF7", color="{BORDER}"];']
    target = trace.nodes.loc[trace.nodes["Propagation_Layer"] == 0, "Gene_Symbol"].iloc[0]
    lines.append(f'"{target}" [shape=doublecircle, fillcolor="{BG_SURFACE_2}", color="{ACCENT}"];')
    for _, edge in edges.iterrows():
        conflict = "conflicting_evidence" in str(edge["Evidence_Status"])
        attrs = f'dir=forward, color="{ACCENT}"' if bool(edge["Directed_Arrow"]) else f'dir=none, color="{BORDER}"'
        if conflict:
            attrs += ', style=dashed, label="çelişkili kanıt"'
        lines.append(f'"{edge["Source"]}" -> "{edge["Target"]}" [{attrs}];')
    lines.append("}")
    st.graphviz_chart("\n".join(lines), width="stretch")
    st.caption("Düz kenar: yapısal STRING PPI. Ok: belirtilen yönde regulatory evidence. Kesikli kenar: çelişkili yönlü kanıt. Hiçbiri tek başına biyolojik nedensellik kanıtı değildir.")


def render_propagation_trace(trace: PropagationTrace, *, target_label: str, tissue: str, species: str) -> bool:
    st.markdown("## Propagation Trace · BETA")
    st.caption(f"Hedef: {target_label} · {species} · {tissue}. Trace mevcut dokuya özgü grafı ve forward simülasyon ölçümlerini kullanır.")
    if trace.nodes.empty:
        st.warning("Yayılım izi için yeterli ölçülmüş ağ düğümü yok.")
        return False
    render_insight_box(_narrative(trace, target_label))
    option = st.radio(
        "Filtre", ["Layer 1", "Layer 1–2", "Layer 1–3+", "Yalnızca anlamlı", "Yönlü kanıt", "Yalnızca kritik"],
        horizontal=True, key="propagation_filter",
    )
    kwargs = {"layer_limit": 3}
    if option == "Layer 1": kwargs["layer_limit"] = 1
    elif option == "Layer 1–2": kwargs["layer_limit"] = 2
    elif option == "Yalnızca anlamlı": kwargs["significant_only"] = True
    elif option == "Yönlü kanıt": kwargs["directed_only"] = True
    elif option == "Yalnızca kritik": kwargs["critical_only"] = True
    visible = filter_trace_nodes(trace.nodes, **kwargs)
    metrics = st.columns(4)
    metrics[0].metric("Etkilenen düğüm", len(visible))
    metrics[1].metric("Anlamlı yeniden-dağılım", int(visible["Significant_Redistribution"].sum()) if not visible.empty else 0)
    metrics[2].metric("Maksimum katman", int(visible["Propagation_Layer"].max()) if not visible.empty else 0)
    metrics[3].metric("Yönlü kanıtlı kenar", int(trace.edges.get("Directed_Evidence", pd.Series(dtype=bool)).sum()))
    _render_graph(trace, visible)
    render_section_header("Propagation graph", "Structural route ve directed evidence ayrı görsel semantiklerle gösterilir.", label="NETWORK TRACE")
    st.markdown("### Düğüm ayrıntıları")
    render_dataframe(visible, label="propagation_trace_visible", key="trace_visible_table")
    with st.expander("Rotalar ve kanıt", expanded=False):
        st.markdown("#### Predicted Network Propagation Routes")
        render_dataframe(trace.routes, label="propagation_trace_routes", key="trace_routes_table")
        st.markdown("#### Yönlü kanıt")
        render_dataframe(trace.evidence, label="propagation_trace_evidence", key="trace_evidence_table")
    render_scientific_note("Propagation trace predicted network propagation gösterir; biyokimyasal sinyal iletimi veya nedensellik kanıtı değildir.")
    with st.expander("Export", expanded=False):
        render_dataframe(trace.nodes, label="propagation_trace_nodes", key="trace_nodes_export", row_limit=500)
        render_dataframe(trace.routes, label="propagation_trace_routes", key="trace_routes_export", row_limit=500)
        render_dataframe(trace.evidence, label="propagation_trace_evidence", key="trace_evidence_export", row_limit=500)
    return st.button("Simülasyon sonucuna dön", key="propagation_trace_close")
