"""Bounded 2D/3D views over an already-computed graph; no scientific mutation."""

from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import streamlit as st

from src.ui.design_tokens import ACCENT, BORDER, BG_SURFACE, TEXT_MUTED, TIER_DUSUK, TIER_KRITIK


def _subgraph(graph: Any, report: pd.DataFrame, targets: list[str], max_nodes: int):
    names = set(map(str, graph.vs["name"]))
    target_set = set(map(str, targets)) & names
    ranked = report.copy()
    metric = "Abs_Delta_PageRank" if "Abs_Delta_PageRank" in ranked else "Hinterland_Skoru"
    if metric in ranked:
        ranked[metric] = pd.to_numeric(ranked[metric], errors="coerce").fillna(0)
        ranked = ranked.sort_values(metric, ascending=False)
    core = [str(value) for value in ranked.get("gene", []) if str(value) in names][:max_nodes]
    selected = set(core) | target_set
    frontier = set(selected)
    while frontier and len(selected) < max_nodes:
        neighbours: set[int] = set()
        for name in frontier:
            neighbours.update(graph.neighbors(graph.vs.find(name=name).index))
        candidates = [str(graph.vs[index]["name"]) for index in neighbours if str(graph.vs[index]["name"]) not in selected]
        add = candidates[: max_nodes - len(selected)]
        selected.update(add)
        frontier = set(add)
        if len(selected) >= min(max_nodes, len(core) + 40):
            break
    indices = [graph.vs.find(name=name).index for name in selected]
    return graph.induced_subgraph(indices)


def build_network_figure(
    graph: Any, report: pd.DataFrame, targets: list[str], gene_to_symbol: dict[str, str],
    *, dimensions: int, max_nodes: int = 350, show_labels: bool = False,
):
    import plotly.graph_objects as go

    sub = _subgraph(graph, report, targets, max_nodes)
    if sub.vcount() == 0:
        return None
    layout_name = "fr3d" if dimensions == 3 else "fr"
    try:
        coords = np.asarray(sub.layout(layout_name).coords, dtype=float)
    except Exception:
        coords = np.asarray(sub.layout("random_3d" if dimensions == 3 else "random").coords, dtype=float)
    if "gene" in report:
        report_index = report.copy()
        report_index["gene"] = report_index["gene"].astype(str)
        report_index = report_index.drop_duplicates("gene").set_index("gene")
    else:
        report_index = pd.DataFrame()
    target_set = set(map(str, targets))
    x_edges: list[float | None] = []; y_edges: list[float | None] = []; z_edges: list[float | None] = []
    for edge in sub.es:
        source, target = edge.tuple
        x_edges.extend([coords[source, 0], coords[target, 0], None])
        y_edges.extend([coords[source, 1], coords[target, 1], None])
        if dimensions == 3:
            z_edges.extend([coords[source, 2], coords[target, 2], None])
    edge_kw = dict(x=x_edges, y=y_edges, mode="lines", line=dict(color=BORDER, width=1), hoverinfo="skip", showlegend=False)
    edge_trace = go.Scatter3d(z=z_edges, **edge_kw) if dimensions == 3 else go.Scatter(**edge_kw)
    colors, sizes, hover, labels = [], [], [], []
    for vertex in sub.vs:
        gene = str(vertex["name"]); symbol = gene_to_symbol.get(gene) or "Alias bulunamadı"
        row = report_index.loc[gene] if not report_index.empty and gene in report_index.index else pd.Series(dtype=object)
        delta = float(row.get("Delta_PageRank", 0) or 0)
        significant = bool(row.get("significant_redistribution", False) or row.get("Strong_Redistribution", False))
        colors.append(ACCENT if gene in target_set else TIER_DUSUK if delta > 0 and significant else TIER_KRITIK if delta < 0 and significant else TEXT_MUTED)
        sizes.append(18 if gene in target_set else 13 if significant else 5)
        hover.append(
            f"<b>{symbol}</b><br>Protein kimliği: {gene}<br>"
            f"Ağdaki Önemi: {float(row.get('Hinterland_Skoru', 0) or 0):.3f}<br>"
            f"Network importance change: {delta:+.6f}<br>Geçiş Merkeziliği: {float(row.get('BC_Skoru', 0) or 0):.3f}"
        )
        labels.append(symbol if show_labels and symbol != "Alias bulunamadı" else "")
    node_kw = dict(
        x=coords[:, 0], y=coords[:, 1], mode="markers+text" if show_labels else "markers",
        marker=dict(size=sizes, color=colors, line=dict(color=BG_SURFACE, width=.8)), text=labels,
        textposition="top center", hovertext=hover, hoverinfo="text", showlegend=False,
    )
    node_trace = go.Scatter3d(z=coords[:, 2], **node_kw) if dimensions == 3 else go.Scatter(**node_kw)
    scene = dict(xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False), bgcolor="rgba(0,0,0,0)")
    layout = go.Layout(
        height=640, margin=dict(l=0, r=0, t=20, b=0), paper_bgcolor="rgba(0,0,0,0)", plot_bgcolor="rgba(0,0,0,0)",
        hovermode="closest", scene=scene if dimensions == 3 else None,
        xaxis=dict(visible=False) if dimensions == 2 else None, yaxis=dict(visible=False) if dimensions == 2 else None,
    )
    return go.Figure(data=[edge_trace, node_trace], layout=layout)


def render_network_views(graph: Any, report: pd.DataFrame, targets: list[str], mapping: dict[str, str]) -> None:
    st.markdown("### Network visualization")
    st.caption("Mavi: pertürbasyon hedefi · Yeşil: pozitif ağ önemi değişimi · Pembe: negatif ağ önemi değişimi · Gri-mavi: başlangıç ağ bağlamı.")
    controls = st.columns([2, 2, 3])
    max_nodes = controls[0].slider("Gösterilecek düğüm", 50, 800, int(st.session_state.get("harita_max_node", 350)), 50)
    show_labels = controls[1].toggle("Gen sembollerini göster", value=False)
    controls[2].caption(f"Aktif ağ: {graph.vcount():,} düğüm / {graph.ecount():,} kenar")
    tab2d, tab3d = st.tabs(["2B ağ", "3B ağ"])
    with tab2d:
        figure = build_network_figure(graph, report, targets, mapping, dimensions=2, max_nodes=max_nodes, show_labels=show_labels)
        st.plotly_chart(figure, width="stretch", key="forward_network_2d") if figure is not None else st.info("Gösterilecek eşleşen düğüm bulunamadı.")
    with tab3d:
        figure = build_network_figure(graph, report, targets, mapping, dimensions=3, max_nodes=max_nodes, show_labels=show_labels)
        st.plotly_chart(figure, width="stretch", key="forward_network_3d") if figure is not None else st.info("Gösterilecek eşleşen düğüm bulunamadı.")
