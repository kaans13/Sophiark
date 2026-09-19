"""Minimum-risk workspace navigasyonu; hesaplama state'ini değiştirmez."""

from __future__ import annotations

import streamlit as st


WORKSPACE_PAGES = ("Özet", "Ağ", "Kritik Değişimler", "Yolaklar", "Kanıt", "Rapor")


def selected_page() -> str:
    return st.radio("Çalışma Alanı", WORKSPACE_PAGES, key="workspace_page", horizontal=True, label_visibility="collapsed")
