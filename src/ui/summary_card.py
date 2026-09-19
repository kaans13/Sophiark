"""Pure HTML builders for the analysis summary card."""

from __future__ import annotations

from collections.abc import Iterable


def summary_cell_html(label: str, value: str, value_class: str = "") -> str:
    return (
        f'<div class="sk-readout-cell"><div class="sk-readout-label">{label}</div>'
        f'<div class="sk-readout-value {value_class}">{value}</div></div>'
    )


def analysis_summary_html(
    *, title: str, target_html: str, context_label: str, cells: Iterable[str],
) -> str:
    """Return compact markup so Markdown never treats nested HTML as a code block."""
    cell_html = "".join(cells)
    target_markup = f'<div class="sk-readout-target">{target_html}</div>' if target_html else ""
    context_markup = f'<div class="sk-eyebrow"><span class="dot"></span>{context_label}</div>' if context_label else ""
    return (
        '<div class="sk-readout-wrap"><div class="sk-readout-head"><div>'
        f'<div class="sk-readout-title">{title}</div>'
        f'{target_markup}</div>{context_markup}'
        f'</div><div class="sk-readout-grid">{cell_html}</div></div>'
    )
