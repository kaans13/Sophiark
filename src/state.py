"""Streamlit state anahtarları ve tür ayrımı için tek sorumluluk noktası."""

from __future__ import annotations

from collections.abc import Callable
from typing import TypeVar

import streamlit as st


_SavedUnifiedResult = TypeVar("_SavedUnifiedResult")


RESEARCH_SESSION_DEFAULTS = {
    "research_current_simulation_id": None,
    "research_current_snapshot_id": None,
    "research_external_snapshot_id": None,
    "research_simulation_scope": None,
    "research_snapshot": None,
    "research_bundle_snapshot_id": None,
    "research_offline_bundle": None,
    "research_explorer_open": False,
    "disease_bank_search_result": None,
    "disease_bank_reference": None,
}
RESEARCH_SESSION_KEYS = tuple(RESEARCH_SESSION_DEFAULTS)


DEFAULTS = {
    "rapor_df": None, "fonksiyon_cache": {}, "fonksiyon_sonuc": None,
    "sim_tamam": False, "sim_hedef": "", "debug_log": [],
    "hayalet_genler": [], "sim_gecmis": [], "enrichment_results": None,
    "enrichment_engine": None, "enrichment_candidate_fingerprint": None,
    "directed_run_bundle": None, "directed_engine_mode": "Classic",
    "signed_redistribution_result": None,
    "simulation_result_context": None,
    "_mygene_query_hash": None, "_last_csv_mtime": {},
    "damping": 0.85, "block_strength": 0.001, "hill_n": 2.0,
    "paralog_boost": 10, "bc_sample_sources": 1200,
    "harita_max_node": 600, "harita_max_etiket": 25,
    "target_stress_results": None, "target_stress_context": None,
    "compensation_results": None, "compensation_context": None,
    "propagation_trace": None, "propagation_context": None,
    "tissue_differential_result": None, "tissue_differential_context": None,
    "enrichment_notices": [],
    "ui_species": "İnsan", "ui_tissue": "None",
    "workspace_report": None, "workspace_graph": None, "workspace_scores": None,
    "unified_research_report": None, "unified_research_context": None,
    "saved_unified_result": None, "unified_history_entry": None, "unified_history_error": None,
    "workspace_targets": [], "workspace_context": None,
    "reports_by_species": {}, "analysis_snapshots": [],
    "target_stress_results": None, "compensation_results": None,
    "target_stress_context": None, "compensation_context": None,
    "disease_results": [], "disease_targets": [],
    **RESEARCH_SESSION_DEFAULTS,
}


def _fresh_default(value):
    return value.copy() if isinstance(value, dict) else list(value) if isinstance(value, list) else value


def build_simulation_result_context(*, species: str, tissue: str, engine: str,
                                    targets=(), localization: str = "",
                                    block_strength: float, damping: float,
                                    candidate_limit: int,
                                    tissue_normalization_mode: str,
                                    bc_sample_sources: int | None = 1200) -> dict:
    """Freeze the UI inputs that identify one completed forward result.

    This is presentation provenance only.  It deliberately contains no
    scientific result values and does not participate in a calculation.
    """
    return {
        "species": str(species),
        "tissue": str(tissue),
        "engine": str(engine),
        "targets": tuple(str(target) for target in targets),
        "localization": str(localization),
        "block_strength": float(block_strength),
        "damping": float(damping),
        "candidate_limit": int(candidate_limit),
        "tissue_normalization_mode": str(tissue_normalization_mode),
        "bc_sample_sources": bc_sample_sources,
    }


def simulation_result_context_matches(context: object, **current) -> bool:
    """Return true only when a stored result belongs to the live UI inputs."""
    if not isinstance(context, dict):
        return False
    return context == build_simulation_result_context(**current)


def bind_simulation_result(*, report, context: dict, ghosts=(), target_label: str = "") -> None:
    """Store a completed result and its provenance as one state transition.

    Both the in-memory engine return and the same-run disk fallback must use this
    boundary.  A report without its context would be indistinguishable from a
    legacy/orphaned CSV and is intentionally rejected here.
    """
    if report is None:
        raise ValueError("completed simulation result requires a report")
    if not isinstance(context, dict) or not context:
        raise ValueError("completed simulation result requires a context")
    st.session_state["rapor_df"] = report
    st.session_state["simulation_result_context"] = dict(context)
    st.session_state["hayalet_genler"] = list(ghosts or ())
    st.session_state["sim_tamam"] = True
    st.session_state["sim_hedef"] = str(target_label)


def initialize() -> None:
    for key, value in DEFAULTS.items():
        if key not in st.session_state:
            st.session_state[key] = _fresh_default(value)


def clear_research_session() -> None:
    """Clear transient Research Context state without touching its SQLite cache."""

    for key, value in RESEARCH_SESSION_DEFAULTS.items():
        st.session_state[key] = _fresh_default(value)


def store_research_simulation_scope(
    *,
    species: str,
    taxon_id: int,
    tissue: str,
    targets,
    attenuation: float | None,
    selection_mode: str,
    test_limit: int | None,
    tested_count: int | None,
    returned_count: int,
    top_n: int | None,
    threshold_parameters: dict | None = None,
) -> None:
    """Freeze the exact transient scope used by a successful simulation run."""

    species_value = str(species).strip()
    selection_value = str(selection_mode).strip()
    if not species_value or int(taxon_id) <= 0 or not selection_value:
        raise ValueError("research simulation scope requires species, taxon and selection mode")
    if returned_count < 0 or (tested_count is not None and tested_count < 0):
        raise ValueError("research simulation scope counts cannot be negative")
    if test_limit is not None and test_limit < 0:
        raise ValueError("research simulation test_limit cannot be negative")
    if top_n is not None and top_n < 0:
        raise ValueError("research simulation top_n cannot be negative")
    st.session_state["research_simulation_scope"] = {
        "species": species_value,
        "taxon_id": int(taxon_id),
        "tissue": str(tissue).strip() or "None",
        "targets": tuple(str(target).strip() for target in targets if str(target).strip()),
        "attenuation": float(attenuation) if attenuation is not None else None,
        "selection_mode": selection_value,
        "test_limit": test_limit,
        "tested_count": tested_count,
        "returned_count": int(returned_count),
        "top_n": int(top_n) if top_n is not None else None,
        "threshold_parameters": dict(threshold_parameters or {}),
    }


def research_simulation_scope() -> dict | None:
    """Return a defensive copy of the simulation-time Research scope, if any."""

    value = st.session_state.get("research_simulation_scope")
    if not isinstance(value, dict):
        return None
    copied = dict(value)
    copied["targets"] = tuple(value.get("targets") or ())
    copied["threshold_parameters"] = dict(value.get("threshold_parameters") or {})
    return copied


def bind_research_snapshot(*, snapshot_id: str, simulation_id: str) -> bool:
    """Bind the current snapshot and evict only stale derived research state.

    Returns ``True`` when an earlier, different snapshot/simulation was
    present.  Persistent provider/evidence caches are deliberately outside
    Streamlit session state and are never changed here.
    """

    snapshot_key = str(snapshot_id).strip()
    simulation_key = str(simulation_id).strip()
    if not snapshot_key or not simulation_key:
        raise ValueError("research snapshot and simulation IDs are required")
    previous_snapshot = st.session_state.get("research_current_snapshot_id")
    previous_simulation = st.session_state.get("research_current_simulation_id")
    stale = (
        previous_snapshot is not None
        and (previous_snapshot != snapshot_key or previous_simulation != simulation_key)
    )
    if previous_snapshot != snapshot_key or previous_simulation != simulation_key:
        for key in (
            "research_external_snapshot_id",
            "research_snapshot",
            "research_bundle_snapshot_id",
            "research_offline_bundle",
            "research_explorer_open",
            "disease_bank_search_result",
            "disease_bank_reference",
        ):
            st.session_state[key] = _fresh_default(RESEARCH_SESSION_DEFAULTS[key])
    st.session_state["research_current_snapshot_id"] = snapshot_key
    st.session_state["research_current_simulation_id"] = simulation_key
    return stale


def research_bundle_for(snapshot_id: str):
    """Return only the offline bundle belonging to the requested snapshot."""

    snapshot_key = str(snapshot_id).strip()
    if (
        snapshot_key
        and st.session_state.get("research_current_snapshot_id") == snapshot_key
        and st.session_state.get("research_bundle_snapshot_id") == snapshot_key
    ):
        return st.session_state.get("research_offline_bundle")
    return None


def store_research_bundle(bundle) -> None:
    """Store an offline bundle only when it matches the active snapshot."""

    snapshot_id = str(getattr(bundle, "snapshot_id", "")).strip()
    simulation_id = str(getattr(bundle, "simulation_id", "")).strip()
    if (
        not snapshot_id
        or snapshot_id != st.session_state.get("research_current_snapshot_id")
        or simulation_id != st.session_state.get("research_current_simulation_id")
    ):
        raise ValueError("research bundle does not belong to the active snapshot")
    st.session_state["research_offline_bundle"] = bundle
    st.session_state["research_bundle_snapshot_id"] = snapshot_id
    st.session_state["research_explorer_open"] = True


def clear_active_workspace() -> None:
    """Bağlam değiştiğinde eski grafın yeni doku/türde kullanılmasını engeller."""
    for key in (
        "workspace_report", "workspace_graph", "workspace_scores", "workspace_targets",
        "unified_research_report", "unified_research_context", "saved_unified_result",
        "unified_history_entry", "unified_history_error",
        "workspace_context", "propagation_trace", "target_stress_results",
        "target_stress_context", "compensation_results", "compensation_context",
        "enrichment_results", "enrichment_notices", "enrichment_engine",
        "enrichment_candidate_fingerprint", "signed_redistribution_result",
        "directed_run_bundle", "directed_engine_mode", "simulation_result_context",
        "dose_df", "sweep_df",
    ):
        if key in DEFAULTS:
            value = DEFAULTS[key]
            st.session_state[key] = _fresh_default(value)
        else:
            st.session_state.pop(key, None)
    clear_research_session()


def invalidate_saved_unified_result() -> None:
    """Clear only the transient bundle-open display state.

    Durable ``.sophiark`` files and their history index are deliberately not
    touched.  Call this at the start of a user-requested bundle operation so a
    failed request cannot be represented by an earlier successful report.
    """
    st.session_state["saved_unified_result"] = None


def load_saved_unified_result(loader: Callable[[], _SavedUnifiedResult]) -> _SavedUnifiedResult:
    """Replace transient display state only after a bundle loader succeeds.

    Invalidating first is intentional: if ``loader`` raises (for example for a
    corrupt bundle), no older report can be rendered as the requested result.
    """
    invalidate_saved_unified_result()
    saved = loader()
    st.session_state["saved_unified_result"] = saved
    return saved


def invalidate_live_unified_result() -> None:
    """Clear the current live Unified presentation before a fresh run starts."""
    st.session_state["unified_research_report"] = None
    st.session_state["unified_research_context"] = None


def reset_simulation_result() -> None:
    """Yeni simülasyondan önce sonuç state'ini temizler; ayarları ve grafı korur."""
    for key in (
        "rapor_df", "fonksiyon_sonuc", "hayalet_genler", "_mygene_query_hash",
        "enrichment_results", "enrichment_notices", "enrichment_engine",
        "enrichment_candidate_fingerprint", "signed_redistribution_result",
        "directed_run_bundle", "directed_engine_mode", "simulation_result_context",
    ):
        if key == "enrichment_notices":
            st.session_state[key] = []
        elif key == "directed_engine_mode":
            st.session_state[key] = "Classic"
        else:
            st.session_state[key] = None
    for key in ("dose_df", "sweep_df"):
        st.session_state.pop(key, None)
    clear_research_session()


def species_history(species: str) -> list[dict]:
    return [item for item in st.session_state.get("sim_gecmis", []) if item.get("tur", species) == species]


def append_history(*, target: str, tissue: str, species: str, stressed_count: int, loss_pct: float, timestamp: str) -> None:
    history = st.session_state["sim_gecmis"]
    history.append({"hedef": target, "doku": tissue, "tur": species, "stresli_sayisi": stressed_count,
                    "kayip_pct": loss_pct, "tarih": timestamp, "stresli_liste": []})
    if len(history) > 10:
        history.pop(0)
