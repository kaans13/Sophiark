"""Read-only presentation helpers for the forward network response.

This module never recalculates, filters, or mutates scientific results.  It
only aligns copies of already-produced result frames for a compact UI view.
"""

from __future__ import annotations

from collections.abc import Iterable
from html import escape

import numpy as np
import pandas as pd
import streamlit as st

from src.ui.components import render_scientific_note, render_section_header
from src.ui.localization import t
from src.ui.metric_presentation import (
    format_metric,
    is_missing as _is_missing,
    metric_help,
    metric_label,
    metric_number_column,
    metric_printf,
    metric_spec,
    metric_text_column,
    truthy as _truthy,
)
from src.ui.tables import render_dataframe


ROLE_TARGET = "Pertürbasyon hedefi"
ROLE_LOSS = "Ağ önemi kaybı"
ROLE_REDISTRIBUTION = "Yeniden dağılım"
GATEWAY_COLUMN = "Gümrük_Kapisi"
COMPARTMENT_BOTTLENECK_LABEL = "Compartment Bottleneck"
BIOLOGICAL_CONTEXT_COLUMNS = (
    "Protein_Adi",
    "MyGene Adı",
    "GO Biyolojik Süreç",
    "GO Moleküler İşlev",
    "GO Hücresel Bileşen",
    "GO_CC_Terimleri",
    "GO_MF_Terimleri",
    "Lokalizasyon",
    "Düzenleyici_TFler",
    "Essentiality",
    "Hedef Gen Etkisi",
)


def _as_frame(value: object) -> pd.DataFrame:
    if isinstance(value, pd.DataFrame):
        return value.copy(deep=True)
    return pd.DataFrame()


def _gene_key(frame: pd.DataFrame) -> pd.Series:
    if "gene" in frame.columns:
        return frame["gene"].astype(str)
    return pd.Series(dtype=str, index=frame.index)


def _first_present(rows: Iterable[pd.Series], column: str) -> object:
    for row in rows:
        if column not in row.index:
            continue
        value = row.get(column)
        if isinstance(value, (list, tuple, set, dict)):
            if value:
                return value
            continue
        try:
            if pd.notna(value) and str(value).strip().lower() not in {"", "nan", "none"}:
                return value
        except (TypeError, ValueError):
            continue
    return np.nan


def _display_text(value: object) -> str:
    if isinstance(value, dict):
        return " | ".join(f"{key}: {item}" for key, item in value.items())
    if isinstance(value, (list, tuple, set)):
        return " | ".join(map(str, value))
    return str(value)


def _display_location(value: object) -> object:
    """Remove storage separators only in the UI; source values stay untouched."""
    try:
        if pd.isna(value):
            return value
    except (TypeError, ValueError):
        pass
    return " ".join(_display_text(value).replace("_", " ").split())


def _display_annotation(value: object) -> str:
    """Present a recorded annotation without turning missing data into a claim."""

    return "—" if _is_missing(value) else _display_text(value)


def _response_table_columns(frame: pd.DataFrame, requested: Iterable[str] | None = None) -> list[str]:
    """Keep recorded biological context and place Gateway last in every export."""

    base = list(requested or frame.columns)
    visible = [
        column for column in base
        if column in frame.columns and column not in {"Kategori", GATEWAY_COLUMN}
    ]
    for column in BIOLOGICAL_CONTEXT_COLUMNS:
        if column in frame.columns and column not in visible:
            visible.append(column)
    if GATEWAY_COLUMN in frame.columns:
        visible.append(GATEWAY_COLUMN)
    return visible


def _display_table(frame: pd.DataFrame, columns: list[str]) -> tuple[pd.DataFrame, pd.DataFrame]:
    source_frame = frame[_response_table_columns(frame, columns)].copy(deep=True)
    display_frame = source_frame.copy(deep=True)
    if "Lokalizasyon" in display_frame.columns:
        display_frame["Lokalizasyon"] = display_frame["Lokalizasyon"].map(_display_location)
    # The underlying scientific key stays unchanged in memory. Only the
    # user-facing download header uses the established product term.
    export_frame = source_frame.rename(columns={GATEWAY_COLUMN: COMPARTMENT_BOTTLENECK_LABEL})
    return display_frame, export_frame


def _nonzero_responses(losses: pd.DataFrame, redistribution: pd.DataFrame) -> pd.DataFrame:
    """Combine existing signed display sets and exclude exact zero without re-ranking."""
    available = [frame for frame in (losses, redistribution) if not frame.empty]
    if not available:
        return pd.DataFrame()
    combined = pd.concat(available, ignore_index=True, sort=False)
    if "Delta_PageRank_Pct" not in combined.columns:
        return combined.copy(deep=True)
    delta = pd.to_numeric(combined["Delta_PageRank_Pct"], errors="coerce")
    combined = combined.loc[delta.notna() & delta.ne(0)].copy(deep=True)
    if "gene" in combined.columns:
        combined = combined.loc[~combined["gene"].astype(str).duplicated(keep="first")].copy()
    return combined


def build_network_response_matrix(
    targets: pd.DataFrame | None,
    losses: pd.DataFrame | None,
    redistribution: pd.DataFrame | None,
) -> pd.DataFrame:
    """Align existing rows by gene ID without changing source order or values."""
    source_frames = [
        (ROLE_TARGET, _as_frame(targets)),
        (ROLE_LOSS, _as_frame(losses)),
        (ROLE_REDISTRIBUTION, _as_frame(redistribution)),
    ]
    order: list[str] = []
    rows_by_gene: dict[str, list[tuple[str, pd.Series]]] = {}
    for role, frame in source_frames:
        if frame.empty or "gene" not in frame.columns:
            continue
        for (_, row), gene in zip(frame.iterrows(), _gene_key(frame)):
            if not gene or gene.casefold() in {"nan", "none"}:
                continue
            if gene not in rows_by_gene:
                rows_by_gene[gene] = []
                order.append(gene)
            rows_by_gene[gene].append((role, row))

    display_rows: list[dict[str, object]] = []
    for gene in order:
        tagged_rows = rows_by_gene[gene]
        roles = [role for role, _ in tagged_rows]
        rows = [row for _, row in tagged_rows]
        symbol = _first_present(rows, "Symbol")
        if pd.isna(symbol):
            symbol = gene
        gateway = _first_present(rows, GATEWAY_COLUMN)
        significant = _first_present(rows, "significant_redistribution")
        delta = pd.to_numeric(
            pd.Series([_first_present(rows, "Delta_PageRank_Pct")]), errors="coerce"
        ).iloc[0]
        response_direction = (
            "↑ Ağ önemi arttı" if pd.notna(delta) and delta > 0
            else "↓ Ağ önemi azaldı" if pd.notna(delta) and delta < 0
            else "—"
        )
        regulatory_relation = _first_present(rows, "Hedef Gen Etkisi")
        if _is_missing(regulatory_relation):
            regulatory_relation = _first_present(rows, "Sign_Context")
        regulatory_relation = (
            "Yön kanıtı yok" if _is_missing(regulatory_relation)
            else _display_text(regulatory_relation)
        )
        response_roles = list(dict.fromkeys(roles))
        if _truthy(significant):
            response_roles.append("İstatistiksel destek")
        display_rows.append({
            "Gen / protein": str(symbol),
            "Protein kimliği": gene,
            "Ağ yanıtı": " · ".join(response_roles),
            "Yanıt yönü": response_direction,
            "Düzenleyici ilişki": regulatory_relation,
            "Ağ önemi değişimi (%)": delta,
            "Hinterland Skoru": pd.to_numeric(
                pd.Series([_first_present(rows, "Hinterland_Skoru")]), errors="coerce"
            ).iloc[0],
            "BC": pd.to_numeric(
                pd.Series([_first_present(rows, "BC_Skoru")]), errors="coerce"
            ).iloc[0],
            "İstatistiksel destek": "✓" if _truthy(significant) else "—",
            "Hücresel konum": _display_location(_first_present(rows, "Lokalizasyon")),
            "Protein adı": _display_annotation(_first_present(rows, "Protein_Adi")),
            "MyGene Adı": _display_annotation(_first_present(rows, "MyGene Adı")),
            "GO Biyolojik Süreç": _display_annotation(_first_present(rows, "GO Biyolojik Süreç")),
            "GO Moleküler İşlev": _display_annotation(_first_present(rows, "GO Moleküler İşlev")),
            "GO Hücresel Bileşen": _display_annotation(_first_present(rows, "GO Hücresel Bileşen")),
            "PubMed araması": _display_annotation(_first_present(rows, "PubMed_Makale")),
            COMPARTMENT_BOTTLENECK_LABEL: "✓" if _truthy(gateway) else "—",
        })
    return pd.DataFrame(display_rows)


def _progress_config(frame: pd.DataFrame, column: str, source_field: str) -> object:
    values = pd.to_numeric(frame.get(column, pd.Series(dtype=float)), errors="coerce")
    finite = values[np.isfinite(values)]
    upper = float(finite.max()) if not finite.empty else 1.0
    lower = float(finite.min()) if not finite.empty else 0.0
    if lower == upper:
        lower = min(0.0, lower)
        upper = max(1.0, upper)
    return st.column_config.ProgressColumn(
        metric_label(source_field),
        help=f"{metric_help(source_field)} Çubuk yalnızca mevcut ham değerin bu görünümdeki görsel temsilidir.",
        width="small",
        format=metric_printf(source_field),
        min_value=lower,
        max_value=upper,
    )


def _render_response_matrix(matrix: pd.DataFrame, *, key: str) -> None:
    matrix_config = {
        "Gen / protein": st.column_config.TextColumn("Gen", width="small", pinned=True),
        "Protein kimliği": st.column_config.TextColumn("ENSP", width="medium"),
        "Ağ yanıtı": st.column_config.TextColumn("Ağ yanıtı", width="medium"),
        "Yanıt yönü": st.column_config.TextColumn("Yanıt yönü", width="medium"),
        "Düzenleyici ilişki": st.column_config.TextColumn(
            "Hedefle ilişki", width="medium",
            help="Yalnız kayıtlı yönlü kanıt. Ağ önemi değişiminin işareti aktivasyon/baskılama olarak yorumlanmaz.",
        ),
        "Ağ önemi değişimi (%)": _progress_config(matrix, "Ağ önemi değişimi (%)", "Delta_PageRank_Pct"),
        "Hinterland Skoru": _progress_config(matrix, "Hinterland Skoru", "Hinterland_Skoru"),
        "BC": _progress_config(matrix, "BC", "BC_Skoru"),
        "İstatistiksel destek": st.column_config.TextColumn("İstatistiksel destek", width="small"),
        "Hücresel konum": st.column_config.TextColumn("Konum", width="medium"),
        "Protein adı": st.column_config.TextColumn("Protein adı", width="medium"),
        "MyGene Adı": st.column_config.TextColumn("MyGene adı", width="medium"),
        "GO Biyolojik Süreç": st.column_config.TextColumn("GO biyolojik süreç", width="large"),
        "GO Moleküler İşlev": st.column_config.TextColumn("GO moleküler işlev", width="large"),
        "GO Hücresel Bileşen": st.column_config.TextColumn("GO hücresel bileşen", width="large"),
        "PubMed araması": st.column_config.LinkColumn(
            "PubMed", help="Seçilen hedef gen ve bu stresli gen için PubMed araması.", display_text="🔗 Ara",
        ),
        COMPARTMENT_BOTTLENECK_LABEL: metric_text_column("Gümrük_Kapisi", width="small"),
    }
    matrix_columns = ["Gen / protein", "Ağ yanıtı", "Yanıt yönü", "Düzenleyici ilişki"]
    for numeric_column in ("Ağ önemi değişimi (%)", "Hinterland Skoru", "BC"):
        if pd.to_numeric(matrix[numeric_column], errors="coerce").notna().any():
            matrix_columns.append(numeric_column)
    for annotation_column in (
        "İstatistiksel destek",
        "Hücresel konum",
        "Protein adı",
        "MyGene Adı",
        "GO Biyolojik Süreç",
        "GO Moleküler İşlev",
        "GO Hücresel Bileşen",
        "PubMed araması",
    ):
        if annotation_column not in matrix.columns:
            continue
        has_recorded_value = matrix[annotation_column].map(
            lambda value: not _is_missing(value) and str(value).strip() != "—"
        ).any()
        if has_recorded_value:
            matrix_columns.append(annotation_column)
    # The requested reading order applies to both the visible matrix and both exports.
    if COMPARTMENT_BOTTLENECK_LABEL in matrix.columns:
        matrix_columns.append(COMPARTMENT_BOTTLENECK_LABEL)
    matrix_columns.append("Protein kimliği")
    export_matrix = matrix[matrix_columns].copy(deep=True)
    render_dataframe(
        export_matrix,
        label=f"{key}_tablosu",
        key=key,
        height=min(620, 80 + len(matrix) * 34),
        column_config=matrix_config,
        export_df=export_matrix,
        row_limit=max(500, len(export_matrix)),
    )


def _analytical_rows_html(matrix: pd.DataFrame, *, limit: int = 18) -> str:
    """Responsive object rows for reading; the exact grid/export remains separate."""
    if matrix.empty:
        return ""
    delta_values = pd.to_numeric(
        matrix.get("Ağ önemi değişimi (%)", pd.Series(index=matrix.index, dtype=float)),
        errors="coerce",
    )
    finite = delta_values[np.isfinite(delta_values)]
    scale = float(finite.abs().max()) if not finite.empty else 1.0
    scale = max(scale, 1e-12)
    visible = matrix.head(limit)
    balance_note = ""
    has_gains = delta_values.gt(0).any()
    has_losses = delta_values.lt(0).any()
    if len(matrix) > limit and has_gains and has_losses:
        neutral = matrix.loc[delta_values.isna() | delta_values.eq(0)].head(2)
        remaining = max(0, limit - len(neutral))
        loss_budget = remaining // 2
        gain_budget = remaining - loss_budget
        losses = matrix.loc[delta_values.lt(0)].head(loss_budget)
        gains = matrix.loc[delta_values.gt(0)].head(gain_budget)
        if len(losses) < loss_budget:
            gains = matrix.loc[delta_values.gt(0)].head(remaining - len(losses))
        elif len(gains) < gain_budget:
            losses = matrix.loc[delta_values.lt(0)].head(remaining - len(gains))
        visible = pd.concat([neutral, losses, gains], axis=0)
        balance_note = f"{len(losses)} azalan + {len(gains)} artan yanıt"
    rows: list[str] = []
    for index, (_, row) in enumerate(visible.iterrows(), start=1):
        delta = pd.to_numeric(pd.Series([row.get("Ağ önemi değişimi (%)")]), errors="coerce").iloc[0]
        tone = "gain" if pd.notna(delta) and delta > 0 else "loss" if pd.notna(delta) and delta < 0 else "neutral"
        delta_label = "—" if pd.isna(delta) else f"{float(delta):+.2f}%"
        bar_width = 0 if pd.isna(delta) else max(4.0, min(100.0, abs(float(delta)) / scale * 100.0))
        importance = pd.to_numeric(pd.Series([row.get("Hinterland Skoru")]), errors="coerce").iloc[0]
        importance_label = "—" if pd.isna(importance) else f"{float(importance):.1f}"
        relation = str(row.get("Düzenleyici ilişki") or "Yön kanıtı yok")
        rows.append(
            f'<div class="sk-object-row sk-object-{tone}">'
            f'<div class="sk-object-rank">{index:02d}</div>'
            '<div class="sk-object-id">'
            f'<b>{escape(str(row.get("Gen / protein", "—")))}</b>'
            f'<small>{escape(str(row.get("Protein kimliği", "—")))}</small></div>'
            '<div class="sk-object-role">'
            f'<span>{escape(str(row.get("Ağ yanıtı", "—")))}</span>'
            f'<small>{escape(str(row.get("Yanıt yönü", "—")))}</small></div>'
            f'<div class="sk-object-delta"><b>{escape(delta_label)}</b>'
            f'<i><em style="width:{bar_width:.1f}%"></em></i></div>'
            f'<div class="sk-object-relation"><span>{escape(relation)}</span></div>'
            f'<div class="sk-object-importance"><small>AĞDAKİ ÖNEMİ</small><b>{escape(importance_label)}</b></div>'
            f'<div class="sk-object-location">{escape(str(row.get("Hücresel konum", "—")))}</div>'
            '</div>'
        )
    note = ""
    if len(matrix) > len(visible):
        summary = f"Özet görünümde {len(visible)} / {len(matrix)} kayıt"
        if balance_note:
            summary += f" · {balance_note}"
        note = f'<div class="sk-object-limit">{summary}; tam veri gridde korunuyor.</div>'
    return (
        '<div class="sk-object-table">'
        '<div class="sk-object-head"><span>#</span><span>GEN / PROTEİN</span><span>ROL VE YANIT</span>'
        '<span>DEĞİŞİM</span><span>HEDEFLE İLİŞKİ</span><span>BAŞLANGIÇ</span><span>HÜCRESEL KONUM</span></div>'
        + "".join(rows) + note + '</div>'
    )


def _render_analytical_rows(matrix: pd.DataFrame) -> None:
    html = _analytical_rows_html(matrix)
    if html:
        st.markdown(html, unsafe_allow_html=True)


def _status(label: str, frame: pd.DataFrame | None, *, not_run: bool = False,
            empty_message: str | None = None) -> None:
    if not_run:
        st.info(f"{label}: analiz çalıştırılmadı veya bu oturumda sonuç objesi bulunmuyor.")
    elif frame is None or frame.empty:
        st.info(empty_message or f"{label}: sonuç bulunamadı.")
    else:
        st.caption(f"{label}: {len(frame):,} kayıt mevcut.")


@st.fragment
def _render_gene_detail(
    matrix: pd.DataFrame,
    source_frames: list[pd.DataFrame],
    *,
    key_prefix: str,
) -> None:
    if matrix.empty:
        return
    options = matrix["Protein kimliği"].astype(str).tolist()
    labels = dict(zip(options, matrix["Gen / protein"].astype(str)))
    selected = st.selectbox(
        "Gen ayrıntısı",
        options,
        format_func=lambda gene: f"{labels.get(gene, gene)} · {gene}",
        key=f"{key_prefix}_gene_detail",
    )
    selected_rows: list[pd.Series] = []
    for frame in source_frames:
        if frame.empty or "gene" not in frame.columns:
            continue
        matches = frame[frame["gene"].astype(str) == selected]
        if not matches.empty:
            selected_rows.append(matches.iloc[0])
    if not selected_rows:
        return

    groups = {
        "Ağ": ["Hinterland_Skoru", "BC_Skoru", "Gümrük_Kapisi", "Topluluk_ID"],
        "Pertürbasyon": ["Delta_PageRank_Pct", "PageRank_Baseline", "PageRank_Perturbed", "Response_Direction", "Selection_Reason"],
        "Doğrulama": ["empirical_p", "q_value", "significant_redistribution"],
        "Biyolojik bağlam": ["Lokalizasyon", "Düzenleyici_TFler", "Essentiality", "Protein_Adi", "MyGene Adı", "GO Biyolojik Süreç", "GO Moleküler İşlev", "GO Hücresel Bileşen"],
    }
    available_groups: list[tuple[str, list[dict[str, str]]]] = []
    for group_name, columns in groups.items():
        values: list[dict[str, str]] = []
        for column in columns:
            value = _first_present(selected_rows, column)
            if not _is_missing(value):
                if column == "Lokalizasyon":
                    display_value = _display_location(value)
                else:
                    display_value = format_metric(column, value) if metric_spec(column) else _display_text(value)
                values.append({"Alan": metric_label(column), "Değer": str(display_value)})
        if values:
            available_groups.append((group_name, values))

    if not available_groups:
        return
    tabs = st.tabs([group_name for group_name, _ in available_groups])
    for group_index, (tab, (group_name, values)) in enumerate(zip(tabs, available_groups)):
        with tab:
            render_dataframe(
                pd.DataFrame(values),
                label=f"{key_prefix}_gen_ayrintisi_{group_index}",
                key=f"{key_prefix}_detail_{group_index}",
                row_limit=50,
            )


def render_network_response(
    *,
    targets: pd.DataFrame | None,
    losses: pd.DataFrame | None,
    redistribution: pd.DataFrame | None,
    target_columns: list[str] | None = None,
    loss_columns: list[str] | None = None,
    redistribution_columns: list[str] | None = None,
    target_column_config: dict | None = None,
    loss_column_config: dict | None = None,
    redistribution_column_config: dict | None = None,
    loss_empty_message: str | None = None,
) -> pd.DataFrame:
    """Render one coherent view over three existing scientific result frames."""
    target_frame = _as_frame(targets)
    loss_frame = _as_frame(losses)
    redistribution_frame = _as_frame(redistribution)
    nonzero_response_frame = _nonzero_responses(loss_frame, redistribution_frame)
    matrix = build_network_response_matrix(target_frame, loss_frame, redistribution_frame)
    loss_matrix = build_network_response_matrix(None, loss_frame, None)
    redistribution_matrix = build_network_response_matrix(None, None, nonzero_response_frame)

    render_section_header(
        "Pertürbasyon ve ağ yanıtı",
        "Hedefler, ağ önemi kayıpları ve yeniden-dağılım aynı ağ cevabının bağlantılı görünümleridir.",
        label=t("network_response"),
    )
    overview_tab, target_tab, loss_tab, redistribution_tab = st.tabs(
        ["Genel görünüm", "Hedefler", "Ağ önemi kaybı", "Yeniden dağılım"]
    )
    with overview_tab:
        if matrix.empty:
            st.info("Bu simülasyonda gösterilecek ağ yanıtı kaydı bulunamadı.")
        else:
            _render_analytical_rows(matrix)
            st.caption("Odak görünümü yeni skor veya sıralama üretmez; artan ve azalan yanıtları birlikte gösterir, her grubun kaynak sırasını korur.")
            with st.expander("Tam matris, sütunlar ve dışa aktarma", expanded=False):
                _render_response_matrix(matrix, key="network_response_overview_matrix")
            with st.expander("Seçili genin ayrıntıları", expanded=False):
                _render_gene_detail(
                    matrix,
                    [target_frame, loss_frame, redistribution_frame],
                    key_prefix="network_response_overview",
                )
    with target_tab:
        _status("Pertürbasyon hedefleri", target_frame)
        if not target_frame.empty:
            visible = [
                column for column in (target_columns or list(target_frame.columns))
                if column in target_frame.columns and column != "Kategori"
            ]
            display, export = _display_table(target_frame, visible)
            render_dataframe(
                display,
                label="birincil_hedefler",
                key="network_response_targets",
                column_config=target_column_config,
                export_df=export,
            )
    with loss_tab:
        _status("Ağ önemi kaybı", loss_frame, empty_message=loss_empty_message)
        if not loss_frame.empty:
            _render_analytical_rows(loss_matrix)
            with st.expander("Seçili kaybın ayrıntıları", expanded=False):
                _render_gene_detail(
                    loss_matrix,
                    [loss_frame],
                    key_prefix="network_response_loss",
                )
            with st.expander("Kaynak ağ önemi kaybı tablosu ve dışa aktarma", expanded=False):
                _render_response_matrix(loss_matrix, key="network_response_loss_matrix")
                visible = [
                    column for column in (loss_columns or list(loss_frame.columns))
                    if column in loss_frame.columns and column != "Kategori"
                ]
                display, export = _display_table(loss_frame, visible)
                render_dataframe(
                    display,
                    label="redistribution_losses",
                    key="network_response_losses",
                    column_config=loss_column_config,
                    export_df=export,
                )
        render_scientific_note("Ağ önemi kaybı ekspresyon, protein miktarı, aktivite veya biyolojik işlev kaybı anlamına gelmez.")
    with redistribution_tab:
        _status("Sıfırdan farklı ağ değişimi", nonzero_response_frame)
        if not nonzero_response_frame.empty:
            _render_analytical_rows(redistribution_matrix)
            if "Delta_PageRank_Pct" in nonzero_response_frame.columns:
                st.caption("Pozitif ve negatif kayıtlar birlikte gösterilir; yalnızca Ağ önemi değişimi (%) değeri tam olarak 0 olan kayıtlar bu görünümün dışındadır.")
            else:
                st.caption("Bu eski raporda Ağ önemi değişimi (%) alanı bulunmadığından mevcut pozitif ve negatif sonuç kayıtları yeniden hesaplanmadan gösterilir.")
            with st.expander("Seçili yanıtın ayrıntıları", expanded=False):
                _render_gene_detail(
                    redistribution_matrix,
                    [nonzero_response_frame],
                    key_prefix="network_response_redistribution",
                )
            with st.expander("Kaynak yeniden dağılım tablosu ve dışa aktarma", expanded=False):
                _render_response_matrix(redistribution_matrix, key="network_response_redistribution_matrix")
                requested_columns = list(dict.fromkeys([
                    *(redistribution_columns or list(nonzero_response_frame.columns)),
                    *(loss_columns or []),
                ]))
                visible = [
                    column for column in requested_columns
                    if column in nonzero_response_frame.columns and column != "Kategori"
                ]
                display, export = _display_table(nonzero_response_frame, visible)
                combined_config = dict(loss_column_config or {})
                combined_config.update(redistribution_column_config or {})
                render_dataframe(
                    display,
                    label="signed_redistribution",
                    key="network_response_redistribution",
                    column_config=combined_config or None,
                    export_df=export,
                )

    return matrix


def render_statistical_support(frame: pd.DataFrame | None, *, columns: list[str]) -> None:
    render_section_header(
        "İstatistiksel destek",
        "Mevcut p-değeri ve FDR alanları ağ yanıtıyla bağlantılı olarak gösterilir; yeniden hesaplama yapılmaz.",
        label="İSTATİSTİKSEL DESTEK",
    )
    if frame is None:
        st.info("İstatistiksel doğrulama bu oturumda çalıştırılmadı.")
        return
    if frame.empty:
        st.info("Mevcut eşiklerde istatistiksel olarak desteklenen aday bulunamadı.")
        return
    visible = [column for column in columns if column in frame.columns]
    if not visible:
        st.info("İstatistiksel doğrulama objesinde gösterilebilir alan bulunamadı.")
        return
    evidence_config = {
        field: (metric_text_column(field) if metric_spec(field).kind != "number" else metric_number_column(field))
        for field in visible
        if metric_spec(field) is not None
    }
    render_dataframe(
        frame[visible],
        label="statistical_validation",
        key="network_response_statistics",
        height=300,
        column_config=evidence_config or None,
    )


def render_functional_context(
    enrichment: pd.DataFrame | None,
    *,
    notices: list[str] | None = None,
    analysis_ran: bool | None = None,
) -> None:
    render_section_header(
        "Fonksiyonel yorum",
        "GO ve KEGG mevcut zenginleştirme sonuçlarının iki görünümüdür; yolak aktivasyonu veya inhibisyonu çıkarımı yapılmaz.",
        label="Biyolojik yorum",
    )
    notice_list = list(notices or [])
    notice_text = " ".join(notice_list).casefold()
    has_error_notice = any(
        marker in notice_text
        for marker in ("tamamlanamadı", "alınamadı", "hazırlanamadı", "erişilemedi", "error", "hata")
    )
    if enrichment is None and analysis_ran is False:
        st.info("GO / KEGG analizi çalıştırılmadı.")
    elif enrichment is None and has_error_notice:
        st.warning("GO / KEGG analizi tamamlanamadı; ağ simülasyonu sonucu korunuyor.")
    elif enrichment is None:
        st.info("GO / KEGG analizinde gösterilecek sonuç bulunamadı.")
    elif enrichment.empty:
        st.info("GO / KEGG analizinde gösterilecek sonuç bulunamadı.")
    else:
        columns = [column for column in ["Term", "Overlap", "Genes"] if column in enrichment.columns]
        labels = {
            "Term": "Terim",
            "Overlap": "Örtüşen genler",
            "Genes": "Genler",
        }
        source = enrichment.get("Kaynak", pd.Series("", index=enrichment.index)).astype(str)
        go_frame = enrichment[source.str.contains("GO", case=False, na=False)].copy()
        kegg_frame = enrichment[source.str.contains("KEGG", case=False, na=False)].copy()
        go_tab, kegg_tab = st.tabs(["GO", "KEGG"])
        with go_tab:
            _status("GO", go_frame)
            if not go_frame.empty:
                render_dataframe(
                    go_frame[columns].rename(columns=labels), label="go_zenginlestirme",
                    key="functional_go", export_df=go_frame,
                )
        with kegg_tab:
            _status("KEGG", kegg_frame)
            if not kegg_frame.empty:
                render_dataframe(
                    kegg_frame[columns].rename(columns=labels), label="kegg_zenginlestirme",
                    key="functional_kegg", export_df=kegg_frame,
                )
    for notice in notice_list:
        st.caption(f"Not: {notice}")


@st.fragment
def render_raw_data(report: pd.DataFrame) -> None:
    render_section_header(
        "Detaylı veri ve dışa aktarma",
        "Ham bilimsel rapor tüm satır ve sütunlarıyla CSV veya Excel olarak indirilebilir.",
        label="DATA & AUDIT",
    )
    if report is None or report.empty:
        st.info("Ham rapor verisi bulunamadı.")
        return
    st.caption(f"{len(report):,} satır · {len(report.columns):,} sütun")
    from src.ui.research_workspace import search_results
    query = st.text_input("Tüm sonuçlarda ara", key="complete_results_search", placeholder="Gen, protein kimliği, anotasyon veya kaynak değer")
    filtered = search_results(report, query)
    st.caption("Arama yalnızca görünümü süzer; indirme her zaman tam raporu içerir. Sütun başlığından sıralayabilir, tablo menüsünden sütunları inceleyebilirsiniz.")
    with st.expander("Tam veri tablosunu aç", expanded=True):
        if filtered.empty:
            st.info("Aramanızla eşleşen kayıt yok. Tam rapora dönmek için arama metnini temizleyin.")
        else:
            render_dataframe(filtered, export_df=report, label="ham_rapor", key="raw_report_data", height=560, paginate=True)
