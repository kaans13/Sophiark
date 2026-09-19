"""Single UI component for persistent species and tissue context."""

from __future__ import annotations

import streamlit as st

from src import state as app_state


SPECIES_OPTIONS = ("İnsan", "Fare")


def _apply_context_selection() -> None:
    """Widget callback runs before the next script pass, avoiding key mutation after render."""
    new_species = st.session_state["context_species"]
    new_tissue = st.session_state["context_tissue"]
    changed = (
        st.session_state.get("ui_species") != new_species
        or st.session_state.get("ui_tissue") != new_tissue
    )
    if changed:
        app_state.clear_active_workspace()
    st.session_state["ui_species"] = new_species
    st.session_state["ui_tissue"] = new_tissue
    st.session_state["tur_secim_radio"] = "🐭 Fare (Mus musculus)" if new_species == "Fare" else "🧑 İnsan (Homo sapiens)"
    st.session_state["hedef_doku"] = new_tissue


def render_context_bar(*, species: str, tissue: str, tissue_options: list[str]) -> tuple[str, str]:
    """Render and update context state without invoking a scientific engine."""
    active_species = st.session_state.setdefault("ui_species", species)
    active_tissue = st.session_state.setdefault("ui_tissue", tissue)
    left, right = st.columns([5, 1], vertical_alignment="center")
    with left:
        st.markdown(f"<div class='sk-context-title'>Sophiark <span>{active_species} · {active_tissue}</span></div>", unsafe_allow_html=True)
        st.caption("Aktif analiz bağlamı")
    with right:
        with st.popover("Değiştir", width="stretch"):
            new_species = st.radio("Tür", SPECIES_OPTIONS, index=SPECIES_OPTIONS.index(active_species), key="context_species")
            new_tissue = st.selectbox("Doku", tissue_options, index=max(0, tissue_options.index(active_tissue)) if active_tissue in tissue_options else 0, key="context_tissue")
            st.button("Bağlamı uygula", type="primary", width="stretch", on_click=_apply_context_selection)
    return active_species, active_tissue
