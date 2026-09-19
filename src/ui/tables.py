"""Streamlit tablo ve güvenli dışa aktarma bileşenleri."""

from __future__ import annotations

import inspect
import io
from dataclasses import dataclass
from math import ceil
import re
from pathlib import Path
from typing import Optional

import pandas as pd
import streamlit as st

from src.ui.metric_presentation import metric_column_config, prepare_display_frame
from src.ui.localization import localize_text


@dataclass(frozen=True)
class DataFramePage:
    """A canonical, order-preserving slice of a result table."""

    frame: pd.DataFrame
    page: int
    page_count: int
    start: int
    end: int
    total: int


def dataframe_page(df: pd.DataFrame, *, page: int, page_size: int) -> DataFramePage:
    """Return one positional page without sorting, filtering, or changing values.

    Page numbers are one-based for the UI.  ``iloc`` keeps the canonical result
    order, so ranks and percentiles remain defined against their original full
    universe rather than against a separately sorted subset.
    """
    if page_size < 1:
        raise ValueError("page_size must be positive")
    total = len(df)
    page_count = max(1, ceil(total / page_size))
    selected = min(max(int(page), 1), page_count)
    start = (selected - 1) * page_size
    end = min(start + page_size, total)
    return DataFramePage(df.iloc[start:end].copy(), selected, page_count, start, end, total)


def _caller_widget_key() -> str:
    """Build a stable fallback from the call site, never from a temporary object id."""
    frame = inspect.currentframe()
    caller = frame.f_back if frame else None
    # Fragment wrappers add Streamlit frames; use the actual application callsite.
    while caller is not None:
        filename = Path(caller.f_code.co_filename)
        if filename != Path(__file__) and "streamlit" not in filename.parts:
            break
        caller = caller.f_back
    if caller is None:
        return "unknown_table"
    stem = Path(caller.f_code.co_filename).stem
    raw = f"{stem}_{caller.f_code.co_name}_{caller.f_lineno}"
    return re.sub(r"[^0-9A-Za-z_]+", "_", raw)


def export_bytes(df: pd.DataFrame, fmt: str) -> tuple[bytes, str, str]:
    """Raporu, karmaşık hücre tiplerini Excel'e uygunlaştırarak dışa aktar."""
    if fmt == "Excel (.xlsx)":
        export_df = df.copy()
        for column in export_df.columns:
            export_df[column] = export_df[column].map(
                lambda value: " | ".join(map(str, value))
                if isinstance(value, (list, tuple, set))
                else str(value) if isinstance(value, dict) else value
            )
        for engine in ("openpyxl", "xlsxwriter"):
            try:
                buffer = io.BytesIO()
                with pd.ExcelWriter(buffer, engine=engine) as writer:
                    export_df.to_excel(writer, index=False, sheet_name="rapor")
                return buffer.getvalue(), (
                    "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                ), "xlsx"
            except Exception:
                continue
        st.warning("Excel oluşturulamadı; CSV indirilebilir durumda.")
    return df.to_csv(index=False).encode("utf-8-sig"), "text/csv", "csv"


@st.fragment
def render_dataframe(
    df: pd.DataFrame,
    *,
    label: str = "",
    height: Optional[int] = None,
    key: Optional[str] = None,
    column_config: Optional[dict] = None,
    export_options: bool = True,
    export_secenekleri: Optional[bool] = None,
    export_df: Optional[pd.DataFrame] = None,
    row_limit: int = 500,
    paginate: bool = False,
    **display_kwargs,
) -> None:
    """Render a table and, when requested, a complete order-preserving inspector."""
    if df is None:
        st.info("Tablo boş.")
        return
    if not isinstance(df, pd.DataFrame):
        df = pd.DataFrame(df)
    if df.empty:
        st.info("Tablo boş.")
        return
    uid = f"{label or 'table'}_{key or _caller_widget_key()}"
    total = len(df)
    if paginate:
        page_size = max(1, int(row_limit))
        page_count = max(1, ceil(total / page_size))
        controls = st.columns([2, 3])
        with controls[0]:
            selected_page = int(st.number_input(
                "Sayfa", min_value=1, max_value=page_count, value=1, step=1,
                key=f"page_{uid}", help="Sayfa konumu canonical sonuç sırasını değiştirmez.",
            ))
        with controls[1]:
            st.caption(f"Tam veri gezgini · {total:,} satır · {page_count:,} sayfa · {page_size:,} satır/sayfa")
        current = dataframe_page(df, page=selected_page, page_size=page_size)
        source_slice = current.frame
        caption = f"Satırlar {current.start + 1:,}–{current.end:,} / {current.total:,} · Sayfa {current.page:,}/{current.page_count:,}"
    else:
        source_slice = df.head(row_limit)
        caption = f"{total:,} satır" + (f"; ilk {row_limit:,} satır gösteriliyor." if total > row_limit else "")
    shown = prepare_display_frame(source_slice)
    if isinstance(shown, pd.DataFrame):
        shown = shown.rename(columns=lambda column: localize_text(str(column)))
    kwargs = {"width": "stretch", "hide_index": True}
    kwargs.update({name: value for name, value in display_kwargs.items() if name in {"width", "hide_index", "use_container_width"}})
    if height:
        kwargs["height"] = min(height, row_limit * 35 + 60)
    kwargs["key"] = key or f"dataframe_{uid}"
    automatic_config = {
        column: config for column in shown.columns
        if (config := metric_column_config(column)) is not None
    }
    if automatic_config or column_config:
        kwargs["column_config"] = {**automatic_config, **(column_config or {})}
    st.dataframe(shown, **kwargs)
    st.caption(caption)
    if export_secenekleri is not None:
        export_options = export_secenekleri
    if not export_options:
        return
    export_source = export_df if isinstance(export_df, pd.DataFrame) else df
    export_total = len(export_source)
    format_col, download_col = st.columns([1, 2])
    with format_col:
        fmt = st.radio(
            "Dışa aktarma formatı", ["CSV (.csv)", "Excel (.xlsx)"],
            key=f"export_format_{uid}", horizontal=True, label_visibility="collapsed",
        )
    payload, mime, extension = export_bytes(export_source, fmt)
    with download_col:
        st.download_button(
            f"{extension.upper()} indir ({export_total:,} satır)", payload,
            f"{label or 'rapor'}.{extension}", mime, key=f"dl_{uid}", width="stretch",
            on_click="ignore",
        )
