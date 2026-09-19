"""Shared presentation-only layout primitives for analysis pages."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Callable
import streamlit as st

from src.ui.design_tokens import INSIGHT_BORDER, INSIGHT_HIGHLIGHT, render_badge


@dataclass(frozen=True)
class AnalysisPageConfig:
    title: str
    description: str
    beta: bool = False


def render_insight_box(text: str, *, tier: str | None = None, is_notable: bool = False) -> None:
    badges = render_badge(tier, tier) if tier else ""
    if is_notable:
        badges += render_badge("Dikkat çekici", "dikkat_cekici")
    st.markdown(f'<div class="sk-insight"><div class="sk-insight-title">Öne Çıkan Bulgu {badges}</div>'
                f'<div class="sk-insight-text">{text}</div></div>', unsafe_allow_html=True)


def render_analysis_page(config: AnalysisPageConfig, *, setup: Callable[[], None], results: Callable[[], None]) -> None:
    """Shared shell only; each analysis retains its own scientific data flow."""
    badge = render_badge("BETA", "beta") if config.beta else ""
    st.markdown(f'<h1 class="sk-title">{config.title} {badge}</h1><p class="sk-subtitle">{config.description}</p>', unsafe_allow_html=True)
    setup()
    results()
