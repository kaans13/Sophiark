from pathlib import Path
import gc
import hashlib
import io
import logging
import os
import pickle
import re
import sys
import tempfile
import time
import traceback
from typing import Optional
from urllib.parse import quote_plus

from src.runtime_logging import ensure_utf8_standard_streams

# Matplotlib'in kullanıcı profil klasörüne yazması bazı kurulumlarda engellenir.
# Proje içinde yazılabilir bir cache kullanmak Streamlit başlangıcını kararlı tutar.
os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".cache" / "matplotlib"))
ensure_utf8_standard_streams()

import matplotlib.pyplot as plt
import altair as alt  # grafikler için (opsiyonel, ama daha modern görünür)
from src.bio_interpreter import interpret_pathways, quick_summary
import igraph as ig
import pandas as pd
import numpy as np
from src.data_loader import get_similar_protein_ligands, query_structural_db

import requests
import streamlit as st

from src import main as human_motor
from src.evidence import main as evidence_motor
from src.metrics import compute_druggable_gate_score

try:
    from src.mouse import main as mouse_main
    MOUSE_MOTOR_OK = True
    _mouse_motor_hata = ""
except Exception as _e:
    mouse_main = None
    MOUSE_MOTOR_OK = False
    _mouse_motor_hata = str(_e)

try:
    from src.mouse import metrics as mouse_metrics  # noqa: F401
except Exception:
    mouse_metrics = None

try:
    from src.mouse import config as mouse_config
except Exception:
    mouse_config = None

motor = human_motor  # varsayılan; TÜR SEÇİMİ bölümünde gerekirse mouse_main'e bağlanır

from src.core.cache_manager import (
    auto_invalidate_on_change,
    force_gc,
    gc_after_dataload,
    scorched_earth_reset,
)
from src import state as app_state
from src.services.analysis_service import execute_simulation, prepare_network
from src.services.directed_analysis import directed_export_tables, run_optional_directed_engine
from src.services.engine_presentation import (
    base_engine_for_mode, candidate_set_fingerprint, compare_ui_frame, directed_ui_frame,
    select_active_result_flow, UNRESOLVED_SYMBOL,
)
from src.services.target_stress import preselect_candidates, validate_candidates
from src.services.compensation import predicted_compensation_candidates
from src.services.enrichment import run_post_simulation_enrichment
from src.services.propagation_trace import build_propagation_trace
from src.services.tissue_differential import run_tissue_differential
from src.cross_species import (
    ComparisonConfig,
    OrthologyIndex,
    check_mouse_suitability,
    compare_ortholog_responses,
    discover_functional_proxies,
    resolve_mouse_entity,
)
from src.ui.pages.target_stress import render_target_stress_results
from src.ui.pages.compensation import render_compensation_results
from src.ui.pages.propagation_trace import render_propagation_trace
from src.ui.pages.tissue_differential import render_tissue_differential
from src.ui.ui_styles import apply_theme
from src.ui.components import (
    candidate_stack_html,
    highlight_card_html,
    render_analysis_hero,
    render_command_identity,
    render_empty_workspace,
    render_landing,
    render_section_header,
    render_scientific_note,
)
from src.ui.app_constants import ENRICHMENT_SOURCE_COUNT, SUPPORTED_SPECIES_COUNT, SessionKeys
from src.ui.gene_cards import render_gene_card
from src.ui.mygene_display import mygene_display_fields
from src.ui.feedback import user_error, user_warning
from src.ui.guides import render_mode_guide, render_results_guide
from src.ui.insight_narrator import en_onemli_bulguyu_ozetle
from src.ui.workspace import render_insight_box
from src.ui.research_workspace import result_destinations, render_candidate_lens, render_reading_guide, align_candidate_evidence
from src.ui.summary_card import analysis_summary_html, summary_cell_html
from src.interpretation import interpret_forward_results, select_redistribution_losses
from src.ui.biological_interpretation import render_biological_interpretation
from src.ui.localization import install_streamlit_localization, locale as ui_locale, t as ui_t
from src.ui.metric_presentation import metric_number_column, metric_text_column
from src.ui.network_response import (
    render_functional_context,
    render_network_response,
    render_raw_data,
    render_statistical_support,
)
from src.research.facade import (
    DiseaseBankService,
    ProvenanceKind,
    ProvenanceRecord,
    LocalResearchContext,
    ResearchConfig,
    ResearchContextService,
    build_research_snapshot,
    render_disease_context,
    render_research_explorer,
)

log = logging.getLogger(__name__)


def _pipeline_log(stage: str, message: str, *args) -> None:
    """Emit grep-friendly progress records to the launch terminal."""
    log.info("[SOPHIARK][%s] " + message, stage, *args)


_pipeline_log("STARTUP", "Uygulama modülleri yüklendi; oturum hazırlanıyor")


@st.cache_data(show_spinner=False)
def _cached_directed_evidence(is_mouse: bool):
    """Load existing directed evidence once; it remains separate from the PPI engine."""
    from src.core.analysis_runtime import directed_evidence_for_selection
    return directed_evidence_for_selection(is_mouse)


@st.cache_resource(show_spinner=False)
def _cached_human_mouse_orthology(project_root: str) -> OrthologyIndex:
    root = Path(project_root)
    return OrthologyIndex.from_snapshot(
        root / "data" / "orthology" / "ensembl_human_mouse_orthology.tsv",
        root / "data" / "orthology" / "ensembl_human_mouse_orthology.metadata.json",
    )


@st.cache_data(show_spinner=False)
def _mouse_suitability_for_human_symbol(project_root: str, human_symbol: str, tissue: str | None) -> dict:
    """Perform a DB-only availability check; never starts a Mouse simulation."""

    import sqlite3

    root = Path(project_root)
    index = _cached_human_mouse_orthology(project_root)
    from src.mouse.config import MOUSE_SYMBOL_MAP

    symbol_to_protein = {
        str(symbol).casefold(): str(protein)
        for protein, symbol in MOUSE_SYMBOL_MAP.items() if symbol
    }
    records = tuple(record for record in index.lookup(human_symbol) if record.mouse_gene)
    candidate_ids = {
        symbol_to_protein.get(str(record.mouse_gene).casefold()) for record in records
    } - {None}
    graph_nodes: set[str] = set()
    context_nodes: set[str] | None = None
    db_path = root / "data" / "raw" / "mouse_hinterland_core.db"
    if db_path.exists() and candidate_ids:
        con = sqlite3.connect(str(db_path))
        try:
            for protein in candidate_ids:
                exists = con.execute(
                    "SELECT 1 FROM mouse_interactions WHERE combined_score >= 500 "
                    "AND (protein1 = ? OR protein2 = ?) LIMIT 1", (protein, protein),
                ).fetchone()
                if exists:
                    graph_nodes.add(str(protein))
            if tissue and str(tissue).casefold() != "none":
                context_nodes = set()
                for protein in candidate_ids:
                    exists = con.execute(
                        "SELECT 1 FROM mouse_tissue_expression WHERE protein_id = ? "
                        "AND lower(tissue) = lower(?) AND expression_level > 0 LIMIT 1",
                        (protein, tissue),
                    ).fetchone()
                    if exists:
                        context_nodes.add(str(protein))
        finally:
            con.close()
    result = check_mouse_suitability(
        human_symbol, index, symbol_to_protein,
        mouse_graph_nodes=graph_nodes, context_nodes=context_nodes, context=tissue,
    )
    return {
        "status": result.status.value,
        "mouse_genes": [record.mouse_gene for record in result.orthologs if record.mouse_gene],
        "mouse_protein_ids": list(result.mouse_protein_ids),
        "reason": result.reason,
        "source": index.metadata.get("source"),
    }


_ORTHOLOGY_STATUS_KEYS = {
    "DIRECT_ORTHOLOG_AVAILABLE": "orthology_status_direct_ortholog_available",
    "ORTHOLOG_EXISTS_NOT_IN_GRAPH": "orthology_status_ortholog_exists_not_in_graph",
    "ORTHOLOG_EXISTS_CONTEXT_UNAVAILABLE": "orthology_status_ortholog_exists_context_unavailable",
    "MULTIPLE_ORTHOLOGS": "orthology_status_multiple_orthologs",
    "NO_DIRECT_ORTHOLOG": "orthology_status_no_direct_ortholog",
    "FUNCTIONAL_PROXY_AVAILABLE": "orthology_status_functional_proxy_available",
    "NO_DEFENSIBLE_PROXY": "orthology_status_no_defensible_proxy",
    "UNRESOLVED": "orthology_status_unresolved",
}
_DIRECTED_PRESENTATION_KEYS = {
    "Relation": {
        "bidirectional": "directed_relation_bidirectional",
        "target_to_entity": "directed_relation_target_to_entity",
        "entity_to_target": "directed_relation_entity_to_target",
        "no_direct_relation": "directed_relation_no_direct_relation",
    },
    "Support": {
        "OmniPath": "directed_support_omnipath",
        "undirected_fallback": "directed_support_undirected_fallback",
        "none": "directed_support_none",
    },
    "Direction": {
        "Influence Gain": "directed_response_influence_gain",
        "Influence Loss": "directed_response_influence_loss",
        "No Change": "directed_response_no_change",
        "Perturbation Target": "directed_response_perturbation_target",
    },
}
def _localized_orthology_status(status: object) -> str:
    """Present a known suitability enum without altering the cached raw result."""
    key = _ORTHOLOGY_STATUS_KEYS.get(str(status))
    return ui_t(key) if key else str(status)


def _directed_presentation_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Localize only known Directed UI values; exports retain raw enums."""
    display = frame.copy(deep=True)
    for column, values in _DIRECTED_PRESENTATION_KEYS.items():
        if column in display:
            display[column] = display[column].map(
                lambda value, mapping=values: ui_t(mapping[str(value)]) if str(value) in mapping else value
            )
    return display


@st.cache_data(show_spinner=False)
def _functional_proxy_candidates(project_root: str, human_protein_id: str, limit: int = 5) -> list[dict]:
    """Read local annotations only; candidates remain opt-in and non-equivalent."""

    root = Path(project_root)
    try:
        with (root / "data" / "processed" / "mygene_cache.pkl").open("rb") as handle:
            human_cache = pickle.load(handle)
        with (root / "data" / "processed" / "mouse_mygene_cache.pkl").open("rb") as handle:
            mouse_cache = pickle.load(handle)
    except Exception:
        return []
    human = human_cache.get(str(human_protein_id), {})
    if not isinstance(human, dict):
        return []
    human_annotations = {
        "go_mf": human.get("go_mf", ()), "go_bp": human.get("go_bp", ()),
        "compartment": human.get("go_cc", ()), "family": (), "pathway": (),
    }
    mouse_annotations = {}
    for protein_id, record in mouse_cache.items():
        if not isinstance(record, dict) or not record.get("symbol"):
            continue
        mouse_annotations[str(record["symbol"])] = {
            "go_mf": record.get("go_mf", ()), "go_bp": record.get("go_bp", ()),
            "compartment": record.get("go_cc", ()), "family": (), "pathway": (),
            "tissue_supported": False,
            "evidence_source": f"local MyGene GO snapshot ({protein_id})",
        }
    return [candidate.__dict__ for candidate in discover_functional_proxies(
        human_annotations, mouse_annotations, limit=limit,
    )]


@st.cache_resource(show_spinner=False)
def _cached_research_service(project_root: str) -> ResearchContextService:
    """Load immutable local research indexes only after explicit UI opt-in."""

    config = ResearchConfig(
        research_context_enabled=True,
        external_context_enabled=False,
        directed_signaling_enabled=True,
        directed_signaling_max_depth=4,
        provider_uniprot_enabled=False,
        provider_interpro_enabled=False,
        provider_quickgo_enabled=False,
        provider_reactome_enabled=False,
        provider_europepmc_enabled=False,
    )
    return ResearchContextService.from_project_root(project_root, config=config)


@st.cache_resource(show_spinner=False)
def _cached_disease_bank_service(project_root: str) -> DiseaseBankService:
    """Construct the optional external-reference service without calling it."""

    return DiseaseBankService.from_project_root(project_root)


_RESEARCH_TOP_N_REASON = re.compile(r"\bTop-(\d+)\b", re.IGNORECASE)


def _research_selection_mode_from_report(report: pd.DataFrame) -> str:
    """Read an exact report mode without inferring scientific behavior."""

    if not isinstance(report, pd.DataFrame):
        return "legacy_unknown"
    modes = tuple(dict.fromkeys(
        value.strip()
        for value in report.get(
            "Redistribution_Selection_Mode", pd.Series(dtype=str)
        ).dropna().astype(str)
        if value.strip()
    ))
    if len(modes) == 1:
        return modes[0]
    return "legacy_unknown" if not modes else "mixed_source"


def _research_top_n_from_report(report: pd.DataFrame) -> int | None:
    """Conservatively recover legacy Top-N only from consistent report reasons."""

    if not isinstance(report, pd.DataFrame) or "Selection_Reason" not in report.columns:
        return None
    reasons = tuple(dict.fromkeys(
        value.strip()
        for value in report["Selection_Reason"].dropna().astype(str)
        if value.strip()
    ))
    if not reasons:
        return None
    matches = tuple(_RESEARCH_TOP_N_REASON.search(reason) for reason in reasons)
    if any(match is None for match in matches):
        return None
    values = {int(match.group(1)) for match in matches if match is not None}
    if len(values) != 1:
        return None
    value = values.pop()
    return value if value > 0 else None


def _research_scientific_threshold_parameters(selection_mode: str) -> dict:
    """Describe the existing strict Top-N rule; never invent other thresholds."""

    if str(selection_mode).strip().casefold() != "top_n":
        return {}
    return {
        "threshold": 0.05,
        "scientific_top_n_metric": "Delta_PageRank_Pct",
        "scientific_top_n_operator": ">",
        "scientific_top_n_strict": True,
        "scientific_top_n_positive_response_only": True,
        "scientific_top_n_unit": "percent",
    }


def _research_unique_gene_count(frame: pd.DataFrame | None) -> int | None:
    """Count the recorded entity universe without treating duplicate rows as entities."""

    if not isinstance(frame, pd.DataFrame):
        return None
    if "gene" not in frame.columns:
        return len(frame)
    identifiers = frame["gene"].dropna().astype(str).str.strip()
    return int(identifiers[identifiers != ""].nunique())


@st.fragment
def _render_research_surface(snapshot) -> None:
    """Automatically render one cached, snapshot-bound offline research bundle."""

    bundle = app_state.research_bundle_for(snapshot.snapshot_id)
    _render_config = ResearchConfig(
        research_context_enabled=True,
        external_context_enabled=False,
        directed_signaling_enabled=True,
        directed_signaling_max_depth=4,
    )
    if bundle is None:
        try:
            with st.spinner("Araştırma bağlamı hazırlanıyor..."):
                service = _cached_research_service(str(BASE_DIR))
                # Streamlit keeps cache_resource values across source hot
                # reloads. An object created from the previous class identity
                # rejects a fresh SimulationResultSnapshot even though its
                # data contract is unchanged, so rebuild only that stale
                # presentation service and leave the simulation untouched.
                if not isinstance(service, ResearchContextService):
                    log.info("Cached Research Context service is stale; rebuilding it")
                    service = ResearchContextService.from_project_root(
                        str(BASE_DIR), config=_render_config,
                    )
                bundle = service.build_offline_bundle(snapshot)
                app_state.store_research_bundle(bundle)
        except Exception as error:
            log.exception("Offline Research Context could not be prepared: %s", error)
            st.warning("Araştırma bağlamı hazırlanamadı; simülasyon sonucu değişmedi.")
            bundle = None
    _bundle_context = getattr(bundle, "local_context", None) if bundle is not None else None
    if bundle is not None and _bundle_context is not None and not isinstance(_bundle_context, LocalResearchContext):
        try:
            log.info("Research display bundle is stale; rebuilding it before rendering")
            service = ResearchContextService.from_project_root(str(BASE_DIR), config=_render_config)
            bundle = service.build_offline_bundle(snapshot)
            app_state.store_research_bundle(bundle)
        except Exception as error:
            log.warning("Research display bundle could not be refreshed: %s", error)
            bundle = None
    if bundle is not None:
        try:
            _bundle_config = getattr(bundle, "config", None)
            if isinstance(_bundle_config, ResearchConfig):
                _render_config = _bundle_config
            render_research_explorer(bundle, ui=st, config=_render_config, include_download=False)
            from src.ui.research_workspace import render_research_download
            render_research_download(bundle, config=_render_config)
        except Exception as error:
            log.warning("Research Explorer render failed safely: %s", error)
            st.warning("Araştırma bağlamı görüntülenemedi; simülasyon sonucu değişmedi.")
    try:
        render_disease_context(
            snapshot,
            service_factory=lambda: _cached_disease_bank_service(str(BASE_DIR)),
            ui=st,
            research_bundle=bundle,
        )
    except Exception as error:
        log.warning("Disease Context render failed safely: %s", error)
        st.info("Hastalık bağlamı kullanılamıyor; simülasyon sonucu değişmedi.")

try:
    from pyvis.network import Network
    PYVIS_OK = True
except ImportError:
    PYVIS_OK = False

try:
    import plotly.graph_objects as go
    PLOTLY_OK = True
except ImportError:
    PLOTLY_OK = False

try:
    import psutil
    PSUTIL_OK = True
except ImportError:
    PSUTIL_OK = False

st.set_page_config(
    page_title="Sophiark · Network Biology",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="collapsed",
)

apply_theme()
install_streamlit_localization()
st.markdown(
    """<style>
    .sk-form-label { display:block; margin:0.55rem 0 0.32rem; min-height:1rem;
        color:#cbd5e1; font-size:0.74rem; font-weight:700; line-height:1.2; }
    .sk-config-title { display:block; margin:0 0 0.65rem; font-size:1rem;
        font-weight:700; line-height:1.25; }
    </style>""",
    unsafe_allow_html=True,
)

BASE_DIR  = Path(__file__).resolve().parent
# One canonical, unfiltered brand asset is used in the header, sidebar and landing page.
LOGO_PATH = BASE_DIR / "assets" / "logo.png"
if not LOGO_PATH.exists():
    LOGO_PATH = BASE_DIR / "assets" / "logo-transparan.png"
MYGENE_URL        = "https://mygene.info/v3/query"
MYGENE_CACHE_PATH = BASE_DIR / "data" / "processed" / "mygene_cache.pkl"

OT_API_URL = "https://api.platform.opentargets.org/api/v4/graphql"

DOKU_SECENEKLERI_INSAN = [
    "None",
    "Adipose Tissue", "Adrenal Gland", "Appendix", "Bone Marrow",
    "Breast", "Bronchus", "Cartilage", "Caudate", "Cerebellum",
    "Cerebral Cortex", "Cervix", "Choroid Plexus", "Colon",
    "Dorsal Raphe", "Duodenum", "Endometrium", "Endometrium 1",
    "Endometrium 2", "Epididymis", "Esophagus", "Eye",
    "Fallopian Tube", "Gallbladder", "Hair", "Heart Muscle",
    "Hippocampus", "Hypothalamus", "Kidney", "Lactating Breast",
    "Liver", "Lung", "Lymph Node", "Nasopharynx", "Oral Mucosa",
    "Ovary", "Pancreas", "Parathyroid Gland", "Pituitary Gland",
    "Placenta", "Prostate", "Rectum", "Retina", "Salivary Gland",
    "Seminal Vesicle", "Skeletal Muscle", "Skin", "Skin 1", "Skin 2",
    "Small Intestine", "Smooth Muscle", "Soft Tissue 1", "Soft Tissue 2",
    "Sole Of Foot", "Spleen", "Stomach 1", "Stomach 2", "Substantia Nigra",
    "Testis", "Thymus", "Thyroid Gland", "Tonsil", "Urinary Bladder", "Vagina",
]

DOM_ROW_LIMIT = 500

# One result-bound command bar owns application identity, context, status and
# the transient configuration action. It deliberately reserves no side rail.
_shell_result = st.session_state.get("simulation_result_context")
_shell_has_report = st.session_state.get("rapor_df") is not None
_shell_live_targets = tuple(
    st.session_state.get("_selected_target_ids")
    or ((_shell_result or {}).get("targets", ()) if isinstance(_shell_result, dict) else ())
)
_shell_live_context = {
    "species": "Mouse" if str(st.session_state.get("tur_secim_radio", "")).startswith("🐭") else "Human",
    "tissue": st.session_state.get("analysis_tissue", (_shell_result or {}).get("tissue", "None")),
    "engine": st.session_state.get("calculation_engine", (_shell_result or {}).get("engine", "Classic")),
    "targets": _shell_live_targets,
    "localization": (_shell_result or {}).get("localization", "SNIPER_MODU"),
    "block_strength": float(st.session_state.get("block_strength", (_shell_result or {}).get("block_strength", 0.001))),
    "damping": float(st.session_state.get("damping", (_shell_result or {}).get("damping", 0.85))),
    "candidate_limit": int(st.session_state.get("redistribution_candidate_limit", (_shell_result or {}).get("candidate_limit", 100))),
    "tissue_normalization_mode": st.session_state.get("tissue_normalization_mode", (_shell_result or {}).get("tissue_normalization_mode", "within_tissue")),
    "bc_sample_sources": st.session_state.get("bc_sample_sources", (_shell_result or {}).get("bc_sample_sources", 1200)),
}
_shell_ready = bool(
    _shell_has_report
    and isinstance(_shell_result, dict)
    and app_state.simulation_result_context_matches(_shell_result, **_shell_live_context)
)
_shell_context = ""
if _shell_ready:
    _shell_selected_targets = st.session_state.get("gen_sec") or ()
    _shell_target_label = " + ".join(
        str(item).rsplit(" · ", 1)[0].strip() for item in _shell_selected_targets
    ) or str(st.session_state.get("_selected_target_display") or st.session_state.get("sim_hedef") or "").strip()
    _shell_report = st.session_state.get("rapor_df")
    if isinstance(_shell_report, pd.DataFrame) and "gene" in _shell_report:
        _shell_symbols = []
        for _shell_target in _shell_result.get("targets", ()):
            _shell_rows = _shell_report.loc[_shell_report["gene"].astype(str).eq(str(_shell_target))]
            _shell_symbol_field = next((field for field in ("Symbol", "gene_symbol") if field in _shell_report), None)
            if not _shell_rows.empty and _shell_symbol_field:
                _shell_symbol = _shell_rows.iloc[0].get(_shell_symbol_field)
                if pd.notna(_shell_symbol) and str(_shell_symbol).strip():
                    _shell_symbols.append(str(_shell_symbol).strip())
        if _shell_symbols:
            _shell_target_label = " + ".join(_shell_symbols)
    _shell_context = " · ".join(filter(None, (
        _shell_target_label,
        str(_shell_result.get("species", "")),
        str(_shell_result.get("tissue", "")),
        str(_shell_result.get("engine", "")),
    )))
_shell_identity, _shell_action = st.columns([8.3, 1.7], vertical_alignment="center")
with _shell_identity:
    render_command_identity(context=_shell_context, ready=_shell_ready)
with _shell_action:
    _config_label = ui_t("change_analysis") if st.session_state.get("rapor_df") is not None else ui_t("new_analysis")
    if st.session_state.get("rapor_df") is None:
        _config_label = ui_t("start_analysis_short")
    _config_panel = st.popover(
        _config_label,
        width="stretch",
    )

_config_panel.markdown(f'<div class="sk-config-title">{ui_t("analysis_settings")}</div>', unsafe_allow_html=True)
_config_workspace_surface = _config_panel.container()
_config_target_surface = _config_panel.container()
_config_context_surface = _config_panel.container()
_config_organism_surface, _config_tissue_surface = _config_context_surface.columns(2, gap="small")
_config_engine_surface = _config_panel.container()
_config_run_surface = _config_panel.container()
_config_advanced_surface = _config_panel.container()
_config_tools_surface = _config_panel.container()
_config_system_surface = _config_panel.container()

workspace_route = _config_workspace_surface.segmented_control(
    ui_t("workspace"), ["Perturbation", "Unified & records"],
    default="Perturbation", key="desktop_workspace_route",
    format_func=lambda value: ui_t("perturbation") if value == "Perturbation" else ui_t("unified_and_records"),
)

_config_organism_surface.markdown(f'<div class="sk-form-label">{ui_t("organism")}</div>', unsafe_allow_html=True)
_tur_secim = _config_organism_surface.segmented_control(
    ui_t("organism"), ["🧑 Human (Homo sapiens)", "🐭 Mouse (Mus musculus)"],
    default="🧑 Human (Homo sapiens)", label_visibility="collapsed", key="tur_secim_radio",
    format_func=lambda value: f"🧑 {ui_t('human')}" if value.startswith("🧑") else f"🐭 {ui_t('mouse')}",
    help="Tür değiştirmek motoru, sembol haritasını, doku listesini ve çıktı klasörünü baştan yükler.",
)
IS_MOUSE = _tur_secim.startswith("🐭")

if IS_MOUSE and mouse_config is not None and hasattr(mouse_config, "TARGET_TO_REGULATORS"):
    TARGET_TO_REGULATORS = mouse_config.TARGET_TO_REGULATORS
    def format_regulators(symbol, max_show=None):
        return mouse_config.format_regulators(symbol, max_show)
    def format_essentiality(symbol, min_test=1):
        return "—"  # farede essentiality verisi yok
else:
    from src.config import TARGET_TO_REGULATORS, format_regulators, format_essentiality

if IS_MOUSE:
    if not MOUSE_MOTOR_OK:
        _config_panel.error(f"❌ mouse_main.py içe aktarılamadı: {_mouse_motor_hata}")
        st.error("Fare motoru (mouse_main.py) yüklenemediği için Fare modu kullanılamıyor.")
        st.stop()
    motor = mouse_main
    _config_organism_surface.caption(ui_t("mouse_protein_network"))
else:
    motor = human_motor
    _config_organism_surface.caption(ui_t("human_protein_network"))

OUTPUTS_DIR = BASE_DIR / ("outputs_mouse" if IS_MOUSE else "outputs")
def _rapor_yolu(rel_adi: str):
    """Fare modunda rapor dosyalarının başına 'mouse_' ekler."""
    return OUTPUTS_DIR / "reports" / (f"mouse_{rel_adi}" if IS_MOUSE else rel_adi)

SOK_RAPORU        = _rapor_yolu("enfeksiyon_sok_dalgasi_raporu.csv")
GENOM_RAPORU      = _rapor_yolu("genom_analiz_sonuc_v5.csv")
DOZ_RAPORU        = _rapor_yolu("farmakolojik_doz_yanit_raporu.csv")
NULL_MODEL_RAPORU = _rapor_yolu("null_model_z_raporu.csv")
ENRICHMENT_RAPORU = _rapor_yolu("validation_enrichment.csv")

TAXON_PREFIX        = "10090." if IS_MOUSE else "9606."
PROTEIN_ID_PREFIX   = "ENSMUSP" if IS_MOUSE else "ENSP"
PROTEIN_ID_ORNEK    = "ENSMUSP00000000001 / Cftr" if IS_MOUSE else "ENSP00000003084 / CFTR"


@st.cache_data(show_spinner=False)
def _fare_dokularini_yukle() -> list:
    import sqlite3
    db_path = BASE_DIR / "data" / "raw" / "mouse_hinterland_core.db"
    if not db_path.exists():
        return ["None"]
    for tablo_adi in ("mouse_tissue_expression", "tissue_expression"):
        try:
            con = sqlite3.connect(str(db_path))
            try:
                rows = pd.read_sql_query(f"SELECT DISTINCT tissue FROM {tablo_adi} ORDER BY tissue", con)
                dokular = ["None"] + [str(t) for t in rows["tissue"].dropna().tolist()]
                return dokular
            finally:
                con.close()
        except Exception:
            continue
    _config_panel.warning("⚠ Fare doku listesi mouse_hinterland_core.db'den okunamadı (tablo adı kontrol edilmeli).")
    return ["None"]


DOKU_SECENEKLERI = _fare_dokularini_yukle() if IS_MOUSE else DOKU_SECENEKLERI_INSAN


@st.cache_data(show_spinner=False)
def _symbol_ensp_tablosunu_yukle() -> pd.DataFrame:
    ENSP_SYMBOLS = BASE_DIR / "data" / "processed" / "ensp_with_symbols.csv"
    if not ENSP_SYMBOLS.exists():
        return pd.DataFrame(columns=["Symbol", "gene"])
    try:
        cols_needed = ["gene", "Symbol", "Hinterland_Skoru"]
        df = pd.read_csv(str(ENSP_SYMBOLS), encoding="utf-8-sig",
                         usecols=lambda c: c in cols_needed)
        if "Symbol" not in df.columns or "gene" not in df.columns:
            return pd.DataFrame(columns=["Symbol", "gene"])
        df = df.dropna(subset=["Symbol"])
        df = df[df["Symbol"].str.strip() != ""]
        df = df.sort_values("Symbol")
        return df.reset_index(drop=True)
    except Exception:
        return pd.DataFrame(columns=["Symbol", "gene"])


@st.cache_data(show_spinner=False)
def _fare_sembol_haritasi_yukle() -> pd.DataFrame:
    # --- fare MyGene kütüphanesini yükle (local) ---
    mouse_cache_path = BASE_DIR / "data" / "processed" / "mouse_mygene_cache.pkl"
    _mygene_kutuphane = {}
    if mouse_cache_path.exists():
        try:
            with open(mouse_cache_path, "rb") as f:
                _mygene_kutuphane = pickle.load(f)
        except Exception:
            _mygene_kutuphane = {}

    if mouse_config is None or not hasattr(mouse_config, "MOUSE_SYMBOL_MAP"):
        return pd.DataFrame(columns=["Symbol", "gene"])
    harita = dict(mouse_config.MOUSE_SYMBOL_MAP)

    # MyGene kütüphanesinden eksik sembolleri tamamla
    if IS_MOUSE:
        for gene_id, kayit in _mygene_kutuphane.items():
            if gene_id in harita and harita[gene_id]:
                continue
            symbol = kayit.get("symbol")
            if symbol and symbol != "—":
                harita[gene_id] = symbol

    df = pd.DataFrame({"gene": list(harita.keys()), "Symbol": list(harita.values())})
    df = df.dropna(subset=["Symbol"])
    df = df[df["Symbol"].astype(str).str.strip() != ""]
    return df.sort_values("Symbol").reset_index(drop=True)


def symbol_den_ensp_coz(symbol: str, tablo: pd.DataFrame) -> Optional[str]:
    if tablo.empty or not symbol:
        return None
    eslesme = tablo[tablo["Symbol"].str.upper() == symbol.upper()]
    if eslesme.empty:
        if IS_MOUSE:
            try:
                mouse_mapping = dict(zip(tablo["Symbol"].astype(str), tablo["gene"].astype(str)))
                resolution = resolve_mouse_entity(
                    symbol, mouse_mapping, _cached_human_mouse_orthology(str(BASE_DIR)),
                )
                return resolution.protein_id if resolution.status == "RESOLVED" else None
            except Exception as error:
                log.warning("Mouse entity resolution unavailable: %s", error)
        return None
    if "Hinterland_Skoru" in eslesme.columns:
        eslesme = eslesme.sort_values("Hinterland_Skoru", ascending=False)
    return str(eslesme.iloc[0]["gene"])


def query_gene_expression(gene_id: str) -> pd.DataFrame:
    import sqlite3
    if IS_MOUSE:
        db_path = str(BASE_DIR / "data" / "raw" / "mouse_hinterland_core.db")
        tablo_adaylari = ("mouse_tissue_expression", "tissue_expression")
    else:
        db_path = str(BASE_DIR / "data" / "raw" / "hinterland_core.db")
        tablo_adaylari = ("tissue_expression",)
    con = sqlite3.connect(db_path)
    try:
        clean_id = gene_id.replace(TAXON_PREFIX, "")
        for tablo_adi in tablo_adaylari:
            try:
                df = pd.read_sql_query(
                    f"SELECT tissue AS 'Doku', MAX(expression_level) AS 'İfade Seviyesi' "
                    f"FROM {tablo_adi} WHERE protein_id = ? AND expression_level > 0 "
                    f"GROUP BY tissue ORDER BY \"İfade Seviyesi\" DESC",
                    con, params=(clean_id,)
                )
                return df
            except Exception:
                continue
        return pd.DataFrame(columns=["Doku", "İfade Seviyesi"])
    finally:
        con.close()


@gc_after_dataload
def csv_yukle_debug(dosya_yolu: Path) -> Optional[pd.DataFrame]:
    yol_str = str(dosya_yolu)
    for enc in ("utf-8-sig", "utf-8", "latin-1"):
        try:
            df = pd.read_csv(yol_str, encoding=enc)
            if not df.empty:
                return df
        except FileNotFoundError:
            return None
        except UnicodeDecodeError:
            continue
        except Exception as e:
            st.error(f"CSV okuma hatası ({enc}): {type(e).__name__}: {e}")
            return None
    return None


def _export_bytes(df: pd.DataFrame, fmt: str) -> tuple:
    if fmt == "Excel (.xlsx)":
        # Excel hücreleri dict/list/set veya timezone'lu değerleri kabul etmez.
        # Raporlardaki tüm hücreleri güvenli, tekil Excel değerlerine dönüştür.
        excel_df = df.copy()
        for col in excel_df.columns:
            excel_df[col] = excel_df[col].map(
                lambda value: (
                    " | ".join(map(str, value)) if isinstance(value, (list, tuple, set))
                    else str(value) if isinstance(value, dict)
                    else value
                )
            )
        errors = []
        for engine in ("openpyxl", "xlsxwriter"):
            try:
                buf = io.BytesIO()
                with pd.ExcelWriter(buf, engine=engine) as writer:
                    excel_df.to_excel(writer, index=False, sheet_name="rapor")
                return buf.getvalue(), "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet", "xlsx"
            except Exception as exc:
                errors.append(f"{engine}: {exc}")
        st.warning("Excel dışa aktarma başarısız; CSV oluşturuldu. " + " | ".join(errors))
    return df.to_csv(index=False).encode("utf-8"), "text/csv", "csv"


def safe_dataframe(
    df: pd.DataFrame,
    label: str = "",
    height: Optional[int] = None,
    key: Optional[str] = None,
    column_config: Optional[dict] = None,
    export_secenekleri: bool = True,
) -> None:
    if df is None or df.empty:
        st.info("Tablo boş.")
        return
    total   = len(df)
    trimmed = df.head(DOM_ROW_LIMIT)
    clipped = total > DOM_ROW_LIMIT
    kwargs  = dict(width="stretch", hide_index=True)
    if height:
        kwargs["height"] = min(height, DOM_ROW_LIMIT * 35 + 60)
    if key:
        kwargs["key"] = key
    if column_config:
        kwargs["column_config"] = column_config
    st.dataframe(trimmed, **kwargs)
    note = (f"▸ {total:,} satırdan ilk {DOM_ROW_LIMIT:,} gösteriliyor."
            if clipped else f"▸ {total:,} satır")
    st.markdown(f'<div class="table-limit-note">{note}</div>', unsafe_allow_html=True)

    _uid = f"{label}_{key or id(df)}"
    if export_secenekleri:
        col_fmt, col_btn = st.columns([1, 2])
        with col_fmt:
            fmt = st.radio(
                "Format", ["CSV (.csv)", "Excel (.xlsx)"],
                key=f"fmt_{_uid}", horizontal=True, label_visibility="collapsed",
            )
        payload, mime, ext = _export_bytes(df, fmt)
        with col_btn:
            dl_label = (f"⬇ Tamamını İndir ({total:,} satır · .{ext})"
                        if clipped else f"⬇ İndir ({total:,} satır · .{ext})")
            st.download_button(dl_label, payload, f"{label or 'rapor'}.{ext}", mime, key=f"dl_{_uid}")
    else:
        dl_label = (f"⬇ Tamamını İndir ({total:,} satır)"
                    if clipped else f"⬇ CSV İndir ({total:,} satır)")
        st.download_button(
            dl_label,
            df.to_csv(index=False).encode("utf-8"),
            f"{label or 'rapor'}.csv",
            "text/csv",
            key=f"dl_{_uid}",
        )


# Tablo/export arayüzü artık uygulama kabuğundan ayrıdır. Eski fonksiyonlar
# geçiş dönemi için yukarıda tutulur; aşağıdaki bağlama tüm ekranlarda modülü kullanır.
from src.ui.tables import render_dataframe as safe_dataframe


def _liste_hash(genler: list) -> str:
    joined = "|".join(sorted(str(g) for g in genler))
    return hashlib.md5(joined.encode()).hexdigest()


def mygene_sorgula(gen_id: str) -> dict:
    try:
        clean_id = str(gen_id).replace("9606.", "")
        payload  = {"q": clean_id, "scopes": "ensembl.protein,ensembl.gene,symbol",
                    "species": "human", "fields": "symbol,name,go.BP"}
        r = requests.post(MYGENE_URL, data=payload,
                          headers={"Content-Type": "application/x-www-form-urlencoded"},
                          timeout=8)
        r.raise_for_status()
        data = r.json()
        if isinstance(data, list) and data and not data[0].get("notfound"):
            hit    = data[0]
            go_raw = hit.get("go", {}).get("BP", [])
            if isinstance(go_raw, dict):
                go_raw = [go_raw]
            go_bp = [e.get("term", "") for e in go_raw[:3] if isinstance(e, dict)]
            return {"symbol": hit.get("symbol", "—"), "name": hit.get("name", "—"), "go_bp": go_bp}
    except requests.exceptions.Timeout:
        return {"symbol": "TIMEOUT", "name": "API yanıt vermedi", "go_bp": []}
    except Exception as e:
        return {"symbol": "HATA", "name": str(e), "go_bp": []}
    return {"symbol": "—", "name": "Bulunamadı", "go_bp": []}


def mygene_guvenli_toplu_sorgu(sorgu_listesi, max_sorgu=50, progress_placeholder=None):
    liste_h = _liste_hash(sorgu_listesi[:max_sorgu])
    if (st.session_state.get("_mygene_query_hash") == liste_h
            and st.session_state.get("fonksiyon_sonuc")):
        return st.session_state["fonksiyon_sonuc"]
    sorgulanan    = sorgu_listesi[:max_sorgu]
    sonuc_listesi = []
    cache         = st.session_state.setdefault("fonksiyon_cache", {})
    for i, gen_id in enumerate(sorgulanan):
        gid  = str(gen_id)
        veri = cache.get(gid)
        if veri is None:
            veri       = mygene_sorgula(gid)
            cache[gid] = veri
            time.sleep(0.15)
        sonuc_listesi.append({
            "Gen ID": gen_id,
            "Symbol": veri["symbol"],
            "Gene Name": veri["name"],
            "GO Biological Process": " | ".join(veri["go_bp"]) or "—",
        })
        if progress_placeholder:
            progress_placeholder.progress((i + 1) / max_sorgu, text=f"{gen_id} → {veri['symbol']}")
    st.session_state["_mygene_query_hash"] = liste_h
    st.session_state["fonksiyon_sonuc"]    = sonuc_listesi
    st.session_state["fonksiyon_cache"]    = cache
    force_gc()
    return sonuc_listesi


@st.cache_data(show_spinner=False)
def _mygene_kutuphanesi_yukle() -> dict:
    if not MYGENE_CACHE_PATH.exists():
        return {}
    try:
        with open(str(MYGENE_CACHE_PATH), "rb") as f:
            raw = pickle.load(f)
    except Exception as e:
        _config_panel.warning(f"⚠ mygene_cache.pkl okunamadı: {e}")
        return {}

    kutuphane: dict = {}
    try:
        if isinstance(raw, dict):
            kutuphane = raw
        elif isinstance(raw, pd.DataFrame):
            idx_col = "gene" if "gene" in raw.columns else raw.columns[0]
            for _, row in raw.iterrows():
                kutuphane[str(row[idx_col])] = row.to_dict()
        elif isinstance(raw, list):
            for item in raw:
                if isinstance(item, dict):
                    anahtar = item.get("query") or item.get("gene") or item.get("_id") or item.get("Gen ID")
                    if anahtar:
                        kutuphane[str(anahtar)] = item
    except Exception as e:
        _config_panel.warning(f"⚠ mygene_cache.pkl formatı tanınamadı: {e}")
        return {}
    return kutuphane


def mygene_kutuphane_coz(gen_id: str, kutuphane: dict) -> dict:
    if not gen_id or not kutuphane:
        return {}
    if gen_id in kutuphane:
        return kutuphane[gen_id]
    clean = str(gen_id).replace("9606.", "")
    if clean in kutuphane:
        return kutuphane[clean]
    prefixed = f"9606.{gen_id}"
    if prefixed in kutuphane:
        return kutuphane[prefixed]
    upper = clean.upper().strip()
    if upper in kutuphane:
        return kutuphane[upper]
    return {}


def mygene_alanlarini_ekle(df: pd.DataFrame, kutuphane: dict) -> pd.DataFrame:
    """Kart ve export için gerçek MyGene annotation alanlarını ekler; yorum üretmez."""
    if df.empty or "gene" not in df.columns:
        return df
    alanlar = df["gene"].astype(str).map(
        lambda gene: mygene_display_fields(mygene_kutuphane_coz(gene, kutuphane))
    )
    return pd.concat([df, pd.DataFrame(alanlar.tolist(), index=df.index)], axis=1)


def _go_terimlerini_cikar(kayit: dict) -> list:
    if not kayit:
        return []
    for anahtar in ("go_bp", "GO Biological Process", "go.BP", "goBP", "go_terms"):
        val = kayit.get(anahtar)
        if not val:
            continue
        if isinstance(val, str):
            return [t.strip() for t in val.replace("|", ";").split(";") if t.strip() and t.strip() != "—"]
        if isinstance(val, (list, tuple)):
            out = []
            for v in val:
                if isinstance(v, dict):
                    out.append(v.get("term", ""))
                else:
                    out.append(str(v))
            return [o for o in out if o]
    return []


_GO_TEMA_SOZLUGU = {
    "İyon/Klor Taşınımı":        ["chloride", "ion transport", "anion transport", "cation transport"],
    "Transmembran Taşıma":       ["transmembrane transport", "membrane transport"],
    "Transkripsiyon Düzenleme":  ["transcription", "dna-binding", "rna polymerase"],
    "Kinaz/Fosforilasyon":       ["kinase", "phosphorylation"],
    "Apoptoz/Hücre Ölümü":       ["apoptotic", "apoptosis", "programmed cell death"],
    "İmmün Yanıt":               ["immune", "inflammatory", "cytokine"],
    "Mitokondriyal Metabolizma": ["mitochondri", "oxidative phosphorylation", "respiratory chain"],
    "Hücre İskeleti":            ["cytoskeleton", "actin", "microtubule"],
    "DNA Onarımı":               ["dna repair", "damage response"],
    "Hücre Döngüsü":             ["cell cycle", "mitotic"],
    "Sinyal İletimi":            ["signaling pathway", "signal transduction"],
    "Protein Katlanması/Degradasyon": ["proteasom", "ubiquitin", "protein folding"],
}


def _ortak_go_temalarini_bul(go_terim_gruplari: list, min_gen: int = 2) -> list:
    tema_sayaci = {tema: 0 for tema in _GO_TEMA_SOZLUGU}
    for terimler in go_terim_gruplari:
        blok = " | ".join(terimler).lower()
        if not blok:
            continue
        for tema, anahtar_kelimeler in _GO_TEMA_SOZLUGU.items():
            if any(kw in blok for kw in anahtar_kelimeler):
                tema_sayaci[tema] += 1
    sonuc = [(tema, n) for tema, n in tema_sayaci.items() if n >= min_gen]
    sonuc.sort(key=lambda x: x[1], reverse=True)
    return sonuc


_TRANSMEMBRAN_ANAHTARLARI = ("cell_membrane", "membrane", "transmembrane", "plasma membrane")


def gen_satiri_yorumla(row: pd.Series, kutuphane_kaydi: dict) -> str:
    parcalar = []
    lok = str(row.get("Lokalizasyon", "") or "").lower()
    tf = str(row.get("Düzenleyici_TFler", "—") or "—")
    ess = str(row.get("Essentiality", "—") or "—")
    gate = str(row.get("Gümrük_Kapisi", "—") or "—")

    is_transmembran = any(k in lok for k in _TRANSMEMBRAN_ANAHTARLARI)

    if is_transmembran and tf.strip() in ("—", "", "nan"):
        parcalar.append(
            "Transmembran/zar lokalizasyonlu bir protein olduğu için doğrudan "
            "transkripsiyon faktörü (TF) ilişkisi beklenmez — bu bir veri "
            "eksikliği değil, biyolojik olarak tutarlı bir sonuçtur."
        )
    elif tf.strip() not in ("—", "", "nan"):
        parcalar.append(f"{tf} tarafından düzenleniyor.")

    go_terimleri = _go_terimlerini_cikar(kutuphane_kaydi)
    if go_terimleri:
        blok = " | ".join(go_terimleri).lower()
        eslesen_temalar = [tema for tema, kws in _GO_TEMA_SOZLUGU.items() if any(kw in blok for kw in kws)]
        if eslesen_temalar:
            parcalar.append(f"GO fonksiyon teması: {', '.join(eslesen_temalar[:2])}.")

    if "yüksek" in ess.lower() or "kritik" in ess.lower():
        parcalar.append("Essentiality verisine göre hücre canlılığı için kritik.")

    if gate.strip() == "✓":
        parcalar.append("Compartment Bottleneck olarak işaretli — ilaçlanabilirlik açısından öncelikli.")

    if not parcalar:
        return "Ek bir MyGene/domain kaydı bulunamadı."
    return " ".join(parcalar)


def yolak_satirini_yorumla(term_row: pd.Series, stresli_semboller_seti: set,
                            symbol_to_gene: dict, kutuphane: dict) -> str:
    overlap_raw = str(term_row.get("Overlap", ""))
    genes_raw   = str(term_row.get("Genes", ""))
    if "/" not in overlap_raw or not genes_raw:
        return ""
    try:
        n_tikanan, n_toplam = overlap_raw.split("/")
        n_tikanan, n_toplam = int(n_tikanan), int(n_toplam)
    except ValueError:
        return ""
    yuzde = (n_tikanan / n_toplam * 100) if n_toplam else 0.0

    yolak_genleri = [g.strip() for g in re.split("[;,]", genes_raw) if g.strip()]
    tikanan_bu_yolakta = [g for g in yolak_genleri if g in stresli_semboller_seti]

    cumle = (
        f"Bu yolakta toplam **{n_toplam}** gen var, bunlardan **{n_tikanan}** tanesi "
        f"(%{yuzde:.0f}) stres simülasyonunda tıkandı."
    )

    if tikanan_bu_yolakta:
        go_gruplari = []
        for sym in tikanan_bu_yolakta:
            ensp = symbol_to_gene.get(sym)
            kayit = mygene_kutuphane_coz(ensp, kutuphane) if ensp else {}
            go_gruplari.append(_go_terimlerini_cikar(kayit))
        temalar = _ortak_go_temalarini_bul(go_gruplari, min_gen=max(2, len(tikanan_bu_yolakta) // 3 or 1))
        if temalar:
            tema_metni = ", ".join(f"{t} ({n} gen)" for t, n in temalar[:2])
            cumle += f" Tıkanan genler ağırlıklı olarak **{tema_metni}** fonksiyonuyla ilişkili."
        else:
            cumle += f" Tıkanan genler: {', '.join(tikanan_bu_yolakta[:8])}{'…' if len(tikanan_bu_yolakta) > 8 else ''}."

    return cumle


def _gen_karti_html(row: "pd.Series") -> str:
    symbol   = row.get("Symbol", row.get("gene", "—"))
    gene_id  = row.get("gene", "—")
    drug_score = row.get("Drug_Score", None)
    hinterland = row.get("Hinterland_Skoru", None)
    bc         = row.get("BC_Skoru", None)
    ess  = str(row.get("Essentiality", "—") or "—")
    gate = str(row.get("Gümrük_Kapisi", "—") or "—")
    mygene_adi = str(row.get("MyGene Adı", "") or "")
    etki  = str(row.get("Hedef Gen Etkisi", "—") or "—")
    pubmed = row.get("PubMed_Makale", None)

    # Pandas/NumPy NaN değerlerini kullanıcıya "nan" olarak göstermeyin.
    if pd.isna(symbol) or str(symbol).strip().lower() in ("", "nan", "none"):
        symbol = gene_id

    def _score(value, precision: int) -> str:
        try:
            value = float(value)
            return f"{value:.{precision}f}" if np.isfinite(value) else "—"
        except (TypeError, ValueError):
            return "—"

    is_gate   = gate.strip() == "✓"
    score_val = _score(drug_score, 1)
    hint_val  = _score(hinterland, 1)
    bc_val    = _score(bc, 3)

    try:
        score_number = float(drug_score)
    except (TypeError, ValueError):
        score_number = float("nan")
    tier = "tier-hi" if np.isfinite(score_number) and score_number >= 65 else \
           "tier-mid" if np.isfinite(score_number) and score_number >= 40 else "tier-lo"

    etiketler = []
    if is_gate:
        etiketler.append('<span class="sk-gene-card-tag gate">✓ Compartment Bottleneck</span>')
    if etki == "Aktivasyon":
        etiketler.append('<span class="sk-gene-card-tag act">↑ Aktivasyon</span>')
    elif etki == "Baskılama":
        etiketler.append('<span class="sk-gene-card-tag rep">↓ Baskılama</span>')
    quote_html = f'<div class="sk-gene-card-quote">MyGene adı: {mygene_adi}</div>' if mygene_adi and mygene_adi != "—" else ""
    pubmed_html = f'<a href="{pubmed}" target="_blank">🔗 PubMed</a>' if pubmed else ""

    parts = [
        f'<div class="sk-gene-card {tier}{" gate" if is_gate else ""}">',
        '<div class="sk-gene-card-head">',
        f'<div class="sk-gene-card-symbol">{symbol}</div>',
        '<div class="sk-gene-card-score">',
        f'<div class="sk-gene-card-score-val">{score_val}</div>',
        '<div class="sk-gene-card-score-label">İlaç Skoru</div>',
        '</div></div>',
        '<div class="sk-gene-card-stats">',
        f'<div class="sk-gene-card-stat"><div class="sk-gene-card-stat-label">Ağdaki Önemi</div><div class="sk-gene-card-stat-val">{hint_val}</div></div>',
        f'<div class="sk-gene-card-stat"><div class="sk-gene-card-stat-label">Betweenness</div><div class="sk-gene-card-stat-val">{bc_val}</div></div>',
        f'<div class="sk-gene-card-stat"><div class="sk-gene-card-stat-label">Essentiality</div><div class="sk-gene-card-stat-val">{ess}</div></div>',
        '</div>',
        f'<div class="sk-gene-card-tags">{"".join(etiketler)}</div>',
        quote_html,
        '<div class="sk-gene-card-footer">',
        f'<span class="sk-gene-card-id">{gene_id}</span>',
        pubmed_html,
        '</div></div>',
    ]
    return "".join(parts)


def uygula_symbol_map(df: pd.DataFrame, cache: dict, col: str = "gene") -> pd.DataFrame:
    if df is None or df.empty or col not in df.columns:
        return df
    df = df.copy()
    df["Symbol"] = df[col].apply(lambda g: cache.get(str(g), {}).get("symbol", "—"))
    return df


def uygula_symbol_map_lokal(df: pd.DataFrame, symbol_map: dict, col: str = "gene") -> pd.DataFrame:
    if df is None or df.empty or col not in df.columns:
        return df
    df = df.copy()
    df["Symbol"] = df[col].map(symbol_map)
    directed_rows = (
        df.get("Engine", pd.Series("classic", index=df.index)).astype(str).str.casefold() == "directed"
    )
    df.loc[directed_rows & df["Symbol"].isna(), "Symbol"] = UNRESOLVED_SYMBOL
    df.loc[~directed_rows & df["Symbol"].isna(), "Symbol"] = df.loc[~directed_rows & df["Symbol"].isna(), col]
    return df


def uygula_trrust_regulators(df: pd.DataFrame, symbol_col: str = "Symbol") -> pd.DataFrame:
    if df is None or df.empty or symbol_col not in df.columns:
        return df
    df = df.copy()
    df["Düzenleyici_TFler"] = df[symbol_col].apply(format_regulators)
    return df


def uygula_essentiality(df: pd.DataFrame, symbol_col: str = "Symbol") -> pd.DataFrame:
    if df is None or df.empty or symbol_col not in df.columns:
        return df
    df = df.copy()
    df["Essentiality"] = df[symbol_col].apply(format_essentiality)
    return df


def render_gate_column(df: pd.DataFrame) -> pd.DataFrame:
    if "Gümrük_Kapisi" not in df.columns:
        return df
    df = df.copy()
    df["Gümrük_Kapisi"] = df["Gümrük_Kapisi"].apply(
        lambda v: "✓" if v is True or str(v).strip().lower() in ("true", "1", "yes") else "—"
    )
    return df


def _sembol_coz(nid: str, smap: dict) -> str:
    if not nid:
        return nid
    if nid in smap:
        return smap[nid]
    clean = nid.replace(TAXON_PREFIX, "")
    if clean in smap:
        return smap[clean]
    prefixed = f"{TAXON_PREFIX}{nid}"
    if prefixed in smap:
        return smap[prefixed]
    nid_upper = clean.upper().strip()
    if nid_upper in smap:
        return smap[nid_upper]
    return nid


def _harita_aliasi(nid: str, smap: dict) -> str:
    """Harita etiketi için yalnızca doğrulanmış sembol döndür; protein ID'sini yazma."""
    sembol = _sembol_coz(nid, smap)
    return str(sembol) if sembol and str(sembol) != str(nid) else ""


def _gercek_ag_alt_grafi(
    g_canli: ig.Graph,
    vuran_genler: list,
    stresli_genler: list,
    hinterland_skorlari: Optional[dict] = None,
    max_nodes: int = 600,
    hop: int = 1,
) -> ig.Graph:
    isim_seti = set(g_canli.vs["name"])
    cekirdek = [n for n in set(vuran_genler + stresli_genler) if n in isim_seti]
    if not cekirdek:
        return g_canli.induced_subgraph([])

    cekirdek_idx = [g_canli.vs.find(name=n).index for n in cekirdek]
    dahil_idx = set(cekirdek_idx)
    skor_map = hinterland_skorlari or {}

    if hop > 0 and len(dahil_idx) < max_nodes:
        komsu_idx: set = set()
        sinir = dahil_idx
        for _ in range(hop):
            yeni = set()
            for i in sinir:
                yeni.update(g_canli.neighbors(i))
            komsu_idx.update(yeni)
            sinir = yeni
            if len(dahil_idx | komsu_idx) >= max_nodes:
                break
        komsu_sirali = sorted(
            komsu_idx - dahil_idx,
            key=lambda i: skor_map.get(g_canli.vs[i]["name"], 0.0),
            reverse=True,
        )
        kalan_yer = max(0, max_nodes - len(dahil_idx))
        dahil_idx.update(komsu_sirali[:kalan_yer])

    if len(dahil_idx) < max_nodes:
        tum_sirali = sorted(
            range(g_canli.vcount()),
            key=lambda i: skor_map.get(g_canli.vs[i]["name"], 0.0),
            reverse=True,
        )
        for i in tum_sirali:
            if len(dahil_idx) >= max_nodes:
                break
            dahil_idx.add(i)

    return g_canli.induced_subgraph(sorted(dahil_idx))


def _koprulu_ilgili_altgraf(g_canli: ig.Graph, ilgili_isimler: list) -> tuple:
    isim_seti = set(g_canli.vs["name"])
    ilgili_isimler = [n for n in ilgili_isimler if n in isim_seti]
    if not ilgili_isimler:
        return g_canli.induced_subgraph([]), set()

    idx_map = {v["name"]: v.index for v in g_canli.vs}
    ilgili_idx = [idx_map[n] for n in ilgili_isimler]
    core_set = set(ilgili_isimler)

    sub_core = g_canli.induced_subgraph(ilgili_idx)
    izole_local = [i for i in range(sub_core.vcount()) if sub_core.degree(i) == 0]
    izole_isimler = [sub_core.vs[i]["name"] for i in izole_local]

    kopru_idx: set = set()
    for isim in izole_isimler:
        gi = idx_map[isim]
        komsu_idx = g_canli.neighbors(gi)
        adaylar = []
        for kj in komsu_idx:
            kj_name = g_canli.vs[kj]["name"]
            if kj_name in core_set:
                continue
            try:
                w = g_canli.es[g_canli.get_eid(gi, kj)]["weight"]
            except Exception:
                w = 0.5
            kj_komsulari = set(g_canli.vs[n]["name"] for n in g_canli.neighbors(kj))
            koprulu_mu = bool(kj_komsulari & (core_set - {isim}))
            adaylar.append((kj, float(w), koprulu_mu))
        adaylar.sort(key=lambda t: (t[2], t[1]), reverse=True)
        for kj, w, koprulu_mu in adaylar[:3]:
            if koprulu_mu:
                kopru_idx.add(kj)

    tum_idx = sorted(set(ilgili_idx) | kopru_idx)
    kopru_isimleri = set(g_canli.vs[i]["name"] for i in kopru_idx)
    return g_canli.induced_subgraph(tum_idx), kopru_isimleri


def _cakismayi_azalt(coords: "np.ndarray", boyutlar: "np.ndarray",
                      iterations: int = 140, pad: float = 1.55) -> "np.ndarray":
    from scipy.spatial import cKDTree
    coords = coords.copy()
    r = (boyutlar / 2.0) * pad
    max_r = float(r.max()) if len(r) else 0.0
    if max_r <= 0:
        return coords
    for _ in range(iterations):
        tree = cKDTree(coords)
        pairs = tree.query_pairs(r=2 * max_r, output_type="ndarray")
        if len(pairs) == 0:
            break
        diff = coords[pairs[:, 0]] - coords[pairs[:, 1]]
        dist = np.linalg.norm(diff, axis=1)
        dist[dist == 0] = 1e-3
        min_dist = r[pairs[:, 0]] + r[pairs[:, 1]]
        overlap = min_dist - dist
        mask = overlap > 0
        if not mask.any():
            break
        push = (diff[mask] / dist[mask, None]) * (overlap[mask, None] * 0.5)
        np.add.at(coords, pairs[mask, 0], push)
        np.add.at(coords, pairs[mask, 1], -push)
    return coords


def draw_incident_scene(
    g_canli: ig.Graph,
    vuran_genler: list,
    stresli_genler: list,
    hinterland_skorlari: Optional[dict] = None,
    symbol_map: Optional[dict] = None,
    tam_ag: bool = False,
    max_nodes: int = 600,
    max_etiket: int = 40,
) -> Optional[str]:
    if not PYVIS_OK:
        return None

    if tam_ag:
        sub_g = _gercek_ag_alt_grafi(
            g_canli, vuran_genler, stresli_genler, hinterland_skorlari,
            max_nodes=max_nodes, hop=1,
        )
        if sub_g.vcount() == 0:
            return None
        try:
            ham_layout = sub_g.layout("fr") if sub_g.vcount() <= 1500 else sub_g.layout("drl")
        except Exception:
            ham_layout = sub_g.layout("random")
        coords = np.array(ham_layout.coords, dtype=float)
        kopru_isimleri = set()
    else:
        ilgili_isimler = list(set(vuran_genler + stresli_genler) & set(g_canli.vs["name"]))
        if not ilgili_isimler:
            return None
        sub_g, kopru_isimleri = _koprulu_ilgili_altgraf(g_canli, ilgili_isimler)
        if sub_g.vcount() == 0:
            return None
        coords = None

    net = Network(height="560px", width="100%", bgcolor="#101012",
                  font_color="#d0d0d5", directed=False)
    vuran_set   = set(vuran_genler)
    stresli_set = set(stresli_genler)
    skor_map    = hinterland_skorlari or {}
    skor_degerleri = [skor_map.get(v["name"], 1.0) for v in sub_g.vs]
    max_skor = max(skor_degerleri) if skor_degerleri and max(skor_degerleri) > 0 else 1.0
    smap = symbol_map or {}

    cekirdek_isimleri_tum = (vuran_set | stresli_set) & set(sub_g.vs["name"])
    cekirdek_sirali = sorted(
        cekirdek_isimleri_tum,
        key=lambda n: skor_map.get(n, 0.0),
        reverse=True,
    )
    etiketli_core = set(cekirdek_sirali[:max_etiket]) if tam_ag else cekirdek_isimleri_tum

    coords_final = None
    boyutlar = None
    if tam_ag:
        boyutlar = np.empty(sub_g.vcount(), dtype=float)
        for i, v in enumerate(sub_g.vs):
            oran = max(0.0, skor_map.get(v["name"], 1.0)) / max_skor
            if v["name"] in vuran_set or v["name"] in stresli_set:
                boyutlar[i] = 12 + (oran ** 0.5) * 22
            else:
                boyutlar[i] = 4 + (oran ** 0.5) * 7

        genislik = (coords[:, 0].max() - coords[:, 0].min()) or 1.0
        yukseklik = (coords[:, 1].max() - coords[:, 1].min()) or 1.0
        hedef_boyut = max(1400.0, 105.0 * (sub_g.vcount() ** 0.5))
        katsayi = hedef_boyut / max(genislik, yukseklik)
        coords_olcekli = coords * katsayi

        coords_final = _cakismayi_azalt(coords_olcekli, boyutlar)

    for idx, v in enumerate(sub_g.vs):
        nid    = v["name"]
        sembol = _harita_aliasi(nid, smap)
        raw    = skor_map.get(nid, 1.0)

        if tam_ag:
            boyut = float(boyutlar[idx])
            # Gerçek ağ görünümünde etiket çizilmez; kimlik hover kartında kalır.
            label  = ""
            px, py = coords_final[idx]
            konum_kw = {"x": float(px), "y": float(py), "physics": False, "fixed": True}
        else:
            boyut = 10 + (raw / max_skor) * 40
            label = str(sembol)[:15] if sembol else ""
            konum_kw = {}

        if nid in vuran_set:
            net.add_node(nid, label=label, shape="diamond", color={"background": "#589bff", "border": "#a6caff",
                                                    "highlight": {"background": "#589bff", "border": "#ffffff"}},
                         size=boyut,
                         title=f"{sembol or 'Alias bulunamadı'} | Kilitlenen Kavşak | Skor: {raw:.2f}",
                         font={"color": "#ffffff", "size": 14, "bold": True, "strokeWidth": 3, "strokeColor": "#000000"},
                         borderWidth=2, borderWidthSelected=3,
                         shadow={"enabled": True, "color": "rgba(255,255,255,0.35)", "size": 12},
                         **konum_kw)
        elif nid in stresli_set:
            net.add_node(nid, label=label, color={"background": "#34bc6e", "border": "#8ce3b0",
                                                    "highlight": {"background": "#34bc6e", "border": "#ffffff"}},
                         size=boyut,
                         title=f"{sembol or 'Alias bulunamadı'} | Stresli Gen | Skor: {raw:.2f}",
                         font={"color": "#f0f0f2", "size": 13, "strokeWidth": 2, "strokeColor": "#000000"},
                         borderWidth=1.5,
                         **konum_kw)
        elif not tam_ag and nid in kopru_isimleri:
            net.add_node(nid, label=str(sembol)[:15] if sembol else "",
                         color={"background": "#3a2f0d", "border": "#fbbf24",
                                "highlight": {"background": "#4a3a10", "border": "#fde68a"}},
                         size=max(9, boyut * 0.6),
                         title=f"{sembol or 'Alias bulunamadı'} | Köprü Düğüm (ara istasyon)",
                         font={"color": "#fbbf24", "size": 11, "strokeWidth": 2, "strokeColor": "#000000"},
                         borderWidth=1.5,
                         **konum_kw)
        else:
            context_bg = "#232326" if tam_ag else "#38383e"
            net.add_node(nid, label=label, color={"background": context_bg, "border": "#45454a",
                                                    "highlight": {"background": "#55555c", "border": "#888890"}},
                         size=boyut if tam_ag else max(9, boyut * 0.55),
                         title=sembol or "Alias bulunamadı",
                         font={"color": "#9a9aa0", "size": 11, "strokeWidth": 2, "strokeColor": "#000000"},
                         borderWidth=1, borderWidthSelected=2,
                         **konum_kw)

    cekirdek_isimleri = vuran_set | stresli_set
    for e in sub_g.es:
        try:
            w = e["weight"]
        except (KeyError, TypeError):
            w = 0.5
        kaynak = sub_g.vs[e.source]["name"]
        hedef_ad = sub_g.vs[e.target]["name"]
        if tam_ag:
            dokunuyor = (kaynak in cekirdek_isimleri) or (hedef_ad in cekirdek_isimleri)
            ikisi_de_cekirdek = (kaynak in cekirdek_isimleri) and (hedef_ad in cekirdek_isimleri)
            if ikisi_de_cekirdek and w < 0.08:
                continue
            opaklik = 0.55 if dokunuyor else 0.08
            genislik_kenar = (0.4 + w * 3.5) if dokunuyor else 0.25
        else:
            opaklik = 0.7
            genislik_kenar = 0.5 + w * 5
        net.add_edge(
            kaynak, hedef_ad,
            color={"color": "#3a3a42", "highlight": "#c0c0c8", "opacity": opaklik},
            width=genislik_kenar,
        )

    if tam_ag:
        net.toggle_physics(False)
        net.set_options("""
        {
          "physics": { "enabled": false },
          "edges": { "smooth": false },
          "interaction": {
            "hover": true,
            "tooltipDelay": 100,
            "dragNodes": true,
            "zoomView": true,
            "dragView": true
          }
        }
        """)
    else:
        net.toggle_physics(True)
        net.set_options("""
        {
          "layout": { "randomSeed": 42 },
          "physics": {
            "barnesHut": {
              "gravitationalConstant": -3000,
              "centralGravity": 0.3,
              "springLength": 120,
              "springConstant": 0.04,
              "damping": 0.09
            },
            "stabilization": { "iterations": 150 }
          },
          "edges": {
            "smooth": { "type": "continuous" }
          },
          "interaction": {
            "hover": true,
            "tooltipDelay": 100
          }
        }
        """)

    # mktemp dosya adını ayırmadığı için yarış koşuluna açıktır; PyVis'in
    # üzerine yazacağı dosya yolunu atomik olarak rezerve ediyoruz.
    with tempfile.NamedTemporaryFile(suffix=".html", delete=False) as temp_html:
        html_path = temp_html.name
    net.save_graph(html_path)
    return html_path


def draw_incident_scene_3d(
    g_canli: ig.Graph,
    vuran_genler: list,
    stresli_genler: list,
    hinterland_skorlari: Optional[dict] = None,
    symbol_map: Optional[dict] = None,
    max_nodes: int = 600,
    max_etiket: int = 40,
):
    if not PLOTLY_OK:
        return None

    sub_g = _gercek_ag_alt_grafi(
        g_canli, vuran_genler, stresli_genler, hinterland_skorlari,
        max_nodes=max_nodes, hop=1,
    )
    if sub_g.vcount() == 0:
        return None

    try:
        layout3d = sub_g.layout("fr3d") if sub_g.vcount() <= 1500 else sub_g.layout("random3d")
    except Exception:
        layout3d = sub_g.layout("random3d")
    coords = np.array(layout3d.coords, dtype=float)

    vuran_set   = set(vuran_genler)
    stresli_set = set(stresli_genler)
    skor_map    = hinterland_skorlari or {}
    smap        = symbol_map or {}
    skor_degerleri = [skor_map.get(v["name"], 1.0) for v in sub_g.vs]
    max_skor = max(skor_degerleri) if skor_degerleri and max(skor_degerleri) > 0 else 1.0

    cekirdek_isimleri_tum = (vuran_set | stresli_set) & set(sub_g.vs["name"])
    cekirdek_sirali = sorted(cekirdek_isimleri_tum, key=lambda n: skor_map.get(n, 0.0), reverse=True)
    etiketli_core = set(cekirdek_sirali[:max_etiket])

    edge_x, edge_y, edge_z = [], [], []
    cekirdek_isimleri = vuran_set | stresli_set
    for e in sub_g.es:
        try:
            w = e["weight"]
        except (KeyError, TypeError):
            w = 0.5
        kaynak = sub_g.vs[e.source]["name"]
        hedef_ad = sub_g.vs[e.target]["name"]
        ikisi_de_cekirdek = (kaynak in cekirdek_isimleri) and (hedef_ad in cekirdek_isimleri)
        if ikisi_de_cekirdek and w < 0.08:
            continue
        x0, y0, z0 = coords[e.source]
        x1, y1, z1 = coords[e.target]
        edge_x += [x0, x1, None]
        edge_y += [y0, y1, None]
        edge_z += [z0, z1, None]

    edge_trace = go.Scatter3d(
        x=edge_x, y=edge_y, z=edge_z, mode="lines",
        line=dict(color="rgba(150,150,160,0.18)", width=1.5),
        hoverinfo="none", showlegend=False,
    )

    def _node_group(pred, size_range, color, opacity, is_core):
        xs, ys, zs, sizes, labels, hovers = [], [], [], [], [], []
        for i, v in enumerate(sub_g.vs):
            nid = v["name"]
            if not pred(nid):
                continue
            oran = max(0.0, skor_map.get(nid, 1.0)) / max_skor
            size = size_range[0] + (oran ** 0.5) * (size_range[1] - size_range[0])
            xs.append(coords[i, 0]); ys.append(coords[i, 1]); zs.append(coords[i, 2])
            sizes.append(size)
            sembol = _harita_aliasi(nid, smap)
            labels.append(str(sembol)[:14] if (is_core and nid in etiketli_core and sembol) else "")
            hovers.append(
                f"{sembol or 'Alias bulunamadı'}<br>Skor: {skor_map.get(nid, 1.0):.2f}"
            )
        return go.Scatter3d(
            x=xs, y=ys, z=zs, mode="markers+text" if is_core else "markers",
            marker=dict(size=sizes, color=color, opacity=opacity, line=dict(width=0.5, color="#000")),
            text=labels, textposition="top center",
            textfont=dict(color="#f5f5f5", size=10),
            hovertext=hovers, hoverinfo="text", showlegend=False,
        )

    vuran_trace   = _node_group(lambda n: n in vuran_set, (10, 26), "#589bff", 1.0, True)
    vuran_trace.marker.symbol = "diamond"
    stresli_trace = _node_group(lambda n: n in stresli_set, (8, 20), "#34bc6e", 0.95, True)
    baglam_trace  = _node_group(lambda n: n not in vuran_set and n not in stresli_set, (2, 6), "#3a3a42", 0.35, False)

    fig = go.Figure(data=[edge_trace, baglam_trace, stresli_trace, vuran_trace])
    fig.update_layout(
        paper_bgcolor="#101012", plot_bgcolor="#101012",
        scene=dict(
            xaxis=dict(visible=False), yaxis=dict(visible=False), zaxis=dict(visible=False),
            bgcolor="#101012",
        ),
        margin=dict(l=0, r=0, t=10, b=0),
        height=620,
        showlegend=False,
    )
    return fig


def ghost_gen_filtrele(df: pd.DataFrame, g_canli: ig.Graph) -> tuple:
    if df is None or df.empty or "gene" not in df.columns:
        return df, []
    all_names = set(g_canli.vs["name"])
    maske     = df["gene"].isin(all_names)
    hayalet   = df.loc[~maske, "gene"].tolist()
    return df[maske].copy(), hayalet




def build_pubmed_link(hedef_semboller, stresli_sembol):
    """Return a safe PubMed co-mention query; this never asserts a relationship."""

    target_symbols = [str(symbol).strip() for symbol in hedef_semboller if str(symbol).strip()]
    stressed_symbol = str(stresli_sembol).strip()
    if not target_symbols or not stressed_symbol:
        return ""
    target_query = " AND ".join(target_symbols)
    query = f"({target_query}) AND {stressed_symbol}"
    return f"https://pubmed.ncbi.nlm.nih.gov/?term={quote_plus(query)}"


@st.cache_data(show_spinner=False)
def _turler_arasi_karsilastirma() -> pd.DataFrame:
    """Son insan ve fare şok raporlarını insan→fare ortologlarıyla eşleştir."""
    human_path = BASE_DIR / "outputs" / "reports" / "enfeksiyon_sok_dalgasi_raporu.csv"
    mouse_path = BASE_DIR / "outputs_mouse" / "reports" / "mouse_enfeksiyon_sok_dalgasi_raporu.csv"
    ortholog_path = BASE_DIR / "data" / "processed" / "mouse_ortholog_map.pkl"
    human_symbol_path = BASE_DIR / "data" / "processed" / "ensp_with_symbols.csv"
    if not all(path.exists() for path in (human_path, mouse_path, ortholog_path, human_symbol_path)):
        return pd.DataFrame()
    try:
        human = pd.read_csv(human_path, encoding="utf-8-sig")
        mouse = pd.read_csv(mouse_path, encoding="utf-8-sig")
        human_symbols = pd.read_csv(human_symbol_path, encoding="utf-8-sig")
        with open(ortholog_path, "rb") as handle:
            human_to_mouse = pickle.load(handle)
    except Exception:
        return pd.DataFrame()
    if "gene" not in human or "gene" not in mouse:
        return pd.DataFrame()

    human_symbol_map = dict(zip(human_symbols["gene"].astype(str), human_symbols["Symbol"].astype(str)))
    mouse_symbol_map = getattr(mouse_config, "MOUSE_SYMBOL_MAP", {}) if mouse_config else {}
    human = human.copy()
    mouse = mouse.copy()
    human["İnsan_Sembol"] = human["gene"].astype(str).map(human_symbol_map)
    human["İnsan_Sembol"] = human["İnsan_Sembol"].fillna(human.get("Symbol")).replace("nan", np.nan)
    mouse_symbol_col = mouse["Symbol"] if "Symbol" in mouse.columns else pd.Series(index=mouse.index, dtype=object)
    mouse["Fare_Sembol"] = mouse_symbol_col.fillna(mouse["gene"].astype(str).map(mouse_symbol_map))
    mouse["Fare_Sembol"] = mouse["Fare_Sembol"].replace(["", "nan", "None"], np.nan)
    human["Fare_Sembol"] = human["İnsan_Sembol"].map(human_to_mouse)
    human = human[human["Fare_Sembol"].notna()]

    keep_h = [column for column in ["İnsan_Sembol", "Fare_Sembol", "Hasar_Tipi", "Hinterland_Skoru", "BC_Skoru", "Efficiency_Kayip_Pct"] if column in human]
    keep_m = [column for column in ["Fare_Sembol", "Hasar_Tipi", "Hinterland_Skoru", "BC_Skoru", "Efficiency_Kayip_Pct"] if column in mouse]
    human = human[keep_h].rename(columns={"Hasar_Tipi": "İnsan_Hasar_Tipi", "Hinterland_Skoru": "İnsan_Hinterland", "BC_Skoru": "İnsan_BC", "Efficiency_Kayip_Pct": "İnsan_Kayıp_Pct"})
    mouse = mouse[keep_m].rename(columns={"Hasar_Tipi": "Fare_Hasar_Tipi", "Hinterland_Skoru": "Fare_Hinterland", "BC_Skoru": "Fare_BC", "Efficiency_Kayip_Pct": "Fare_Kayıp_Pct"})
    # Outer merge, raporların ikisi de varken ortak ortolog yoksa bunu görünür
    # kılar; önceki inner merge yanlışlıkla "rapor yok" mesajına düşüyordu.
    merged = human.merge(mouse, on="Fare_Sembol", how="outer", indicator=True)
    for column in ["İnsan_Hinterland", "Fare_Hinterland", "İnsan_BC", "Fare_BC"]:
        if column in merged:
            merged[column] = pd.to_numeric(merged[column], errors="coerce")
    merged["Hinterland_Farkı"] = merged["İnsan_Hinterland"] - merged["Fare_Hinterland"]
    merged["Korunmuş_Etki"] = (
        merged["İnsan_Hasar_Tipi"].notna()
        & merged["Fare_Hasar_Tipi"].notna()
        & (merged["İnsan_Hasar_Tipi"] == merged["Fare_Hasar_Tipi"])
    )
    merged["Karşılaştırma_Durumu"] = merged["_merge"].map({
        "both": "Her iki türde", "left_only": "Yalnızca insan", "right_only": "Yalnızca fare",
    })
    return merged.drop(columns="_merge").sort_values(
        "Hinterland_Farkı", key=lambda values: values.abs().fillna(-1), ascending=False
    )



def search_diseases(query: str) -> list:
    graphql_query = {
        "query": """
        query($queryString: String!) {
            search(queryString: $queryString, entityNames: ["disease"]) {
                hits {
                    id
                    name
                }
            }
        }
        """,
        "variables": {"queryString": query}
    }
    try:
        r = requests.post(OT_API_URL, json=graphql_query, timeout=10)
        r.raise_for_status()
        data = r.json()
        if "errors" in data:
            error_msg = data["errors"][0]["message"]
            st.error(f"📛 Open Targets API GraphQL hatası: {error_msg}")
            return []
        hits = data.get("data", {}).get("search", {}).get("hits", [])
        return [{"id": hit["id"], "name": hit["name"]} for hit in hits]
    except requests.exceptions.Timeout:
        st.error("⏱️ Open Targets API zaman aşımına uğradı. Lütfen tekrar deneyin.")
        return []
    except requests.exceptions.ConnectionError:
        st.error("🌐 Open Targets API'ye bağlanılamadı. İnternet bağlantınızı kontrol edin.")
        return []
    except requests.exceptions.HTTPError as e:
        status = e.response.status_code if e.response is not None else "?"
        body = e.response.text[:200] if e.response is not None else ""
        st.error(f"🌐 API HTTP hatası: {status} — {body}")
        return []
    except ValueError:
        st.error("❌ Open Targets API beklenmeyen (JSON olmayan) bir yanıt döndürdü.")
        return []
    except Exception as e:
        st.error(f"❌ Beklenmeyen hata (search_diseases): {type(e).__name__}: {e}")
        return []

def get_genes_for_disease(efo_id: str) -> list:
    graphql_query = {
        "query": """
        query($efoId: String!, $page: Pagination!) {
            disease(efoId: $efoId) {
                associatedTargets(page: $page) {
                    count
                    rows {
                        target {
                            approvedSymbol
                        }
                        score
                    }
                }
            }
        }
        """,
        "variables": {"efoId": efo_id, "page": {"index": 0, "size": 500}}
    }
    try:
        r = requests.post(OT_API_URL, json=graphql_query, timeout=15)
        r.raise_for_status()
        data = r.json()
        if "errors" in data:
            error_msg = data["errors"][0]["message"]
            st.error(f"📛 API GraphQL hatası: {error_msg}")
            return []
        disease_node = data.get("data", {}).get("disease")
        if disease_node is None:
            st.warning("⚠️ Bu EFO ID için hastalık verisi bulunamadı.")
            return []
        rows = disease_node.get("associatedTargets", {}).get("rows", [])
        genes = []
        for row in rows:
            symbol = row.get("target", {}).get("approvedSymbol")
            if symbol and symbol not in genes:
                genes.append(symbol)
        return genes
    except requests.exceptions.Timeout:
        st.error("⏱️ Open Targets API zaman aşımına uğradı.")
        return []
    except requests.exceptions.ConnectionError:
        st.error("🌐 Open Targets API'ye bağlanılamadı. İnternet bağlantınızı kontrol edin.")
        return []
    except requests.exceptions.HTTPError as e:
        status = e.response.status_code if e.response is not None else "?"
        st.error(f"🌐 API HTTP hatası: {status}")
        return []
    except ValueError:
        st.error("❌ Open Targets API beklenmeyen (JSON olmayan) bir yanıt döndürdü.")
        return []
    except Exception as e:
        st.error(f"❌ Beklenmeyen hata (get_genes_for_disease): {type(e).__name__}: {e}")
        return []




def get_target_effect_on_stressed(target_symbols: list, stressed_symbol: str) -> str:
    if not target_symbols or not stressed_symbol:
        return "—"

    regulators = TARGET_TO_REGULATORS.get(stressed_symbol, [])
    if not regulators:
        return "—"

    for tf, effect in regulators:
        if tf in target_symbols:
            return effect
    return "—"
app_state.initialize()
_pipeline_log("SESSION", "Oturum durumu hazırlandı")

_stale_detected = auto_invalidate_on_change(SOK_RAPORU)

symbol_tablo = _fare_sembol_haritasi_yukle() if IS_MOUSE else _symbol_ensp_tablosunu_yukle()
if not symbol_tablo.empty:
    _secenekler = [
        f"{row['Symbol']}  ·  {row['gene']}"
        for _, row in symbol_tablo.iterrows()
    ]
    _gene_to_symbol = dict(zip(symbol_tablo["gene"], symbol_tablo["Symbol"]))
else:
    _secenekler = []
    _gene_to_symbol = {}

_symbol_to_gene = {v: k for k, v in _gene_to_symbol.items()}
_pipeline_log("IDENTIFIERS", "Protein→gen sembol haritası hazır (n=%d)", len(_gene_to_symbol))

if IS_MOUSE:
    MOUSE_MYGENE_CACHE_PATH = BASE_DIR / "data" / "processed" / "mouse_mygene_cache.pkl"
    @st.cache_data(show_spinner=False)
    def _mygene_kutuphanesi_yukle_mouse() -> dict:
        if not MOUSE_MYGENE_CACHE_PATH.exists():
            return {}
        try:
            with open(str(MOUSE_MYGENE_CACHE_PATH), "rb") as f:
                return pickle.load(f)
        except Exception as e:
            return {}
    _mygene_kutuphane = _mygene_kutuphanesi_yukle_mouse()
else:
    _mygene_kutuphane = _mygene_kutuphanesi_yukle()


with _config_panel:
    _config_tissue_surface.markdown(f'<div class="sk-form-label">{ui_t("tissue")}</div>', unsafe_allow_html=True)
    hedef_doku = _config_tissue_surface.selectbox(
        ui_t("tissue"), DOKU_SECENEKLERI,
        index=None if "analysis_tissue" in st.session_state else 1,
        label_visibility="collapsed",
        key="analysis_tissue",
        help="Seçilen dokuda üretilmeyen genler analizden çıkarılır.",
    )
    doku_aktif = hedef_doku != "None"
    if doku_aktif:
        _config_tissue_surface.caption(f"{hedef_doku} · {ui_t('tissue_filter_active')}")
    else:
        _config_tissue_surface.caption(ui_t("tissue_filter_off"))

    _config_engine_surface.markdown(f'<div class="sk-form-label">{ui_t("calculation_engine")}</div>', unsafe_allow_html=True)
    _engine_options = ["Classic"] if IS_MOUSE else ["Classic", "Directed", "Evidence", "Compare"]
    calculation_engine = _config_engine_surface.segmented_control(
        ui_t("calculation_engine"), _engine_options, default="Classic",
        key="calculation_engine", label_visibility="collapsed",
        format_func=lambda mode: "Directed · BETA" if mode == "Directed" else ("Evidence · BETA" if mode == "Evidence" else mode),
        help=ui_t("calculation_engine_help"),
    )
    if calculation_engine == "Directed":
        _config_engine_surface.caption("BETA · Directional result layer; it does not reproduce Classic-only metrics." if ui_locale() == "en" else "BETA · Yönlü sonuç katmanı; Classic-only metrikleri yeniden üretmez.")
    elif calculation_engine == "Evidence":
        _config_engine_surface.caption("BETA · E1 non-text weighting mode; text mining is limited validation only." if ui_locale() == "en" else "BETA · E1 non-text ağırlık modu; metin madenciliği yalnız sınırlı doğrulamadır.")
    elif calculation_engine == "Compare":
        _config_engine_surface.caption("Compares Classic and Directed results for the same input." if ui_locale() == "en" else "Classic ve Directed sonuçlarını aynı girdide karşılaştırır.")
    if IS_MOUSE:
        _config_engine_surface.caption("Mouse Directed mapping gate: UNAVAILABLE (%66.30 atomic edge coverage; gerekli eşik %70).")
    elif calculation_engine == "Evidence":
        motor = evidence_motor
    else:
        motor = human_motor

    _config_target_surface.markdown(f'<div class="sk-form-label">{ui_t("target")}</div>', unsafe_allow_html=True)

    hedef_genler_sniper: Optional[list] = None
    lokalizasyon = "Mitochondrion"

    MOD_SECENEKLERI = ["gene", "disease", "localization"]
    mod = _config_target_surface.segmented_control(
        ui_t("target_type"), MOD_SECENEKLERI, default="gene",
        label_visibility="collapsed", key="mod",
        format_func=lambda value: {"gene": ui_t("search_gene"), "disease": ui_t("select_disease"), "localization": ui_t("manual_localization")}[value],
    )

    if mod == "gene":
        if _secenekler:
            secilen_genler = _config_target_surface.multiselect(
                "gen_sec",
                options=_secenekler,
                max_selections=10,
                key="gen_sec",
                label_visibility="collapsed",
                help="Bir veya birden fazla gen seçin (en çok 10). Kombine hasarda ortak kenarlar yalnızca bir kez susturulur.",
            )
            if secilen_genler:
                hedef_genler_sniper = []
                for secilen_str in secilen_genler:
                    ensp_parsed = secilen_str.rsplit(" · ", 1)[-1].strip()
                    if ensp_parsed not in hedef_genler_sniper:
                        hedef_genler_sniper.append(ensp_parsed)
                lokalizasyon = "SNIPER_MODU"
                selected_labels = ", ".join(
                    item.rsplit(" · ", 1)[0].strip() for item in secilen_genler
                )
                if (
                    st.session_state.get("_selected_target_display") != selected_labels
                    or tuple(st.session_state.get("_selected_target_ids", ())) != tuple(hedef_genler_sniper)
                ):
                    st.session_state["_selected_target_display"] = selected_labels
                    st.session_state["_selected_target_ids"] = tuple(hedef_genler_sniper)
                    st.rerun()
                _config_target_surface.markdown(
                    f'<div class="ensp-resolved">▸ {len(hedef_genler_sniper)} hedef: {selected_labels}'
                    f'<br>▸ <span style="color:#85858c">{" | ".join(hedef_genler_sniper)}</span></div>',
                    unsafe_allow_html=True,
                )
                if not IS_MOUSE:
                    for _human_target_id in hedef_genler_sniper:
                        _human_target_symbol = _gene_to_symbol.get(_human_target_id, _human_target_id)
                        try:
                            _mouse_check = _mouse_suitability_for_human_symbol(
                                str(BASE_DIR), str(_human_target_symbol),
                                hedef_doku if doku_aktif else None,
                            )
                            _mouse_genes = ", ".join(_mouse_check["mouse_genes"]) or "—"
                            _config_target_surface.caption(
                                f"{ui_t('mouse_test')} · {_human_target_symbol}: {_localized_orthology_status(_mouse_check['status'])} "
                                f"· {_mouse_genes}"
                            )
                            if _mouse_check["status"] in {"NO_DIRECT_ORTHOLOG", "UNRESOLVED"}:
                                if _config_target_surface.button(
                                    f"View Functional Proxy Candidates · {_human_target_symbol}",
                                    key=f"mouse_proxy_{_human_target_id}",
                                ):
                                    _proxy_rows = _functional_proxy_candidates(
                                        str(BASE_DIR), str(_human_target_id), limit=5,
                                    )
                                    if _proxy_rows:
                                        st.warning(
                                            "Functional proxy ≠ ortholog/equivalent gene. "
                                            "No candidate is automatically used as a perturbation target."
                                        )
                                        _proxy_frame = pd.DataFrame(_proxy_rows)
                                        st.dataframe(_proxy_frame, hide_index=True, width="stretch")
                                        st.download_button(
                                            "Proxy evidence CSV",
                                            _proxy_frame.to_csv(index=False).encode("utf-8-sig"),
                                            file_name=f"{_human_target_symbol}_mouse_functional_proxy_candidates.csv",
                                            mime="text/csv",
                                            key=f"mouse_proxy_export_{_human_target_id}",
                                        )
                                    else:
                                        _config_target_surface.caption("NO_DEFENSIBLE_PROXY · Yerel anotasyon kanıtı yetersiz.")
                        except Exception as _mouse_check_error:
                            log.warning("Mouse suitability check unavailable: %s", _mouse_check_error)
                            _config_target_surface.caption(
                                f"{ui_t('mouse_test')} · {_human_target_symbol}: "
                                f"{_localized_orthology_status('UNRESOLVED')}"
                            )
            else:
                st.session_state.pop("_selected_target_display", None)
                st.session_state.pop("_selected_target_ids", None)
                _config_target_surface.caption(ui_t("select_one_or_more_genes"))
        else:
            if IS_MOUSE:
                st.warning("⚠️ mouse_config.MOUSE_SYMBOL_MAP boş ya da bulunamadı.")
            else:
                st.warning("⚠️ ensp_with_symbols.csv bulunamadı.")

        _tissue_query_panel = _config_advanced_surface.expander("Doku ifadesi sorgusu", expanded=False)
        gene_query = _tissue_query_panel.text_input(f"{PROTEIN_ID_PREFIX} veya Sembol", key="gene_query",
                                   placeholder=PROTEIN_ID_ORNEK,
                                   label_visibility="collapsed")
        if _tissue_query_panel.button("Dokuları sorgula", key="query_btn", width="stretch"):
            if gene_query:
                query_id = gene_query.strip()
                if not query_id.startswith(PROTEIN_ID_PREFIX):
                    resolved = None
                    if not symbol_tablo.empty:
                        match = symbol_tablo[symbol_tablo["Symbol"].str.upper() == query_id.upper()]
                        if not match.empty:
                            resolved = match.iloc[0]["gene"]
                    if resolved:
                        query_id = resolved
                        st.markdown(f'<div class="status-ok">▸ {query_id}</div>', unsafe_allow_html=True)
                    else:
                        st.warning("Sembol çevrilemedi.")
                        query_id = None
                if query_id and query_id.startswith(PROTEIN_ID_PREFIX):
                    result = query_gene_expression(query_id)
                    if not result.empty:
                        st.markdown(f'<div class="status-ok">▸ {len(result)} dokuda ifade</div>', unsafe_allow_html=True)
                        safe_dataframe(result, label="doku_ifadesi", key="tissue_expression_table", height=220)
                    else:
                        st.info("İfade verisi bulunamadı.")


    elif mod == "disease":
        st.markdown('<div class="sidebar-label">🧬 Hastalık Seç</div>', unsafe_allow_html=True)
        disease_query = _config_target_surface.text_input("Hastalık ara", key="disease_query",
                                      placeholder="Örn: Alzheimer, Parkinson, Diabetes...",
                                      label_visibility="collapsed")
        if disease_query:
            with st.spinner("Hastalık aranıyor..."):
                hastalik_listesi = search_diseases(disease_query)
            if hastalik_listesi:
                disease_options = {f"{h['name']} ({h['id']})": h for h in hastalik_listesi}
                secilen = _config_target_surface.selectbox("Hastalık seçin", list(disease_options.keys()),
                                       key="disease_select")
                if secilen:
                    secilen_hastalik = disease_options[secilen]
                    st.markdown(f'<div class="status-ok">▸ {secilen_hastalik["name"]}</div>',
                                unsafe_allow_html=True)
                    with st.spinner("Genler getiriliyor..."):
                        gen_sembolleri = get_genes_for_disease(secilen_hastalik["id"])
                    if gen_sembolleri:
                        st.markdown(f'<div class="status-ok">▸ {len(gen_sembolleri)} ilişkili gen bulundu</div>',
                                    unsafe_allow_html=True)
                        hedef_genler_sniper = []
                        cevrilemeyen = []
                        for sym in gen_sembolleri:
                            ensp = symbol_den_ensp_coz(sym, symbol_tablo)
                            if ensp:
                                hedef_genler_sniper.append(ensp)
                            else:
                                cevrilemeyen.append(sym)
                        lokalizasyon = f"Hastalık: {secilen_hastalik['name']}"
                        if cevrilemeyen:
                            st.warning(f"⚠ {len(cevrilemeyen)} sembol ENSP'ye çevrilemedi: {', '.join(cevrilemeyen[:10])}")
                        if hedef_genler_sniper:
                            st.success(f"✅ {len(hedef_genler_sniper)} gen hedef olarak seçildi.")
                        else:
                            st.error("Hiçbir gen seçilemedi.")
                    else:
                        st.warning("Bu hastalık için hedef gen bulunamadı.")
            else:
                st.warning("Hastalık bulunamadı ya da API'den sonuç dönmedi (üstteki hata mesajına bakın).")

    elif mod == "localization":
        lokalizasyon = _config_target_surface.selectbox(
            "lok",
            ["Mitochondrion", "Nucleus", "Cell_Membrane", "Endoplasmic_Reticulum",
             "Golgi_Apparatus", "Lysosome", "Cytoplasm", "Peroxisome"],
            label_visibility="collapsed",
            format_func=lambda value: str(value).replace("_", " "),
        )
        st.markdown(f'<div class="status-ok">▸ {str(lokalizasyon).replace("_", " ")}</div>', unsafe_allow_html=True)

    st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 0.8rem 0;"></div>', unsafe_allow_html=True)

    with _config_advanced_surface.expander(ui_t("advanced_settings"), expanded=False):
        st.selectbox(
            ui_t("number_genes_to_test"), options=[20, 30, 50, 100, 200], index=[20, 30, 50, 100, 200].index(int(st.session_state.get("redistribution_candidate_limit", 100))),
            key="redistribution_candidate_limit", help=ui_t("gene_count_help"),
        )
        if "damping" not in st.session_state:
            st.session_state.damping = 0.85
        if "block_strength" not in st.session_state:
            st.session_state.block_strength = 0.001
        if "hill_n" not in st.session_state:
            st.session_state.hill_n = 2.0
        if "paralog_boost" not in st.session_state:
            st.session_state.paralog_boost = 10
        if "bc_sample_sources" not in st.session_state:
            st.session_state.bc_sample_sources = 1200
        if "harita_max_node" not in st.session_state:
            st.session_state.harita_max_node = 600
        if "harita_max_etiket" not in st.session_state:
            st.session_state.harita_max_etiket = 25

        st.session_state.damping = st.slider(
            ui_t("network_importance_spread"),
            min_value=0.70, max_value=0.99, value=st.session_state.damping, step=0.01,
            help="Sinyalin ağda ne kadar yayılacağını belirler. Yüksek değer = daha global etki."
        )
        st.session_state.block_strength = st.slider(
            ui_t("remaining_edge_weight"),
            min_value=0.0001, max_value=0.1, value=st.session_state.block_strength, step=0.0001,
            format="%.4f",
            help="Susturma sonrası kalan kenar ağırlığı. 0.0001 = çok güçlü susturma, 0.1 = daha hafif/kısmi inhibisyon."
        )
        st.session_state.hill_n = st.slider(
            ui_t("hill_coefficient"),
            min_value=0.5, max_value=4.0, value=st.session_state.hill_n, step=0.5,
            help="Doz-yanıt eğrisinin dikliğini kontrol eder."
        )
        st.session_state.paralog_boost = st.slider(
            ui_t("paralog_boost"),
            min_value=0, max_value=50, value=st.session_state.paralog_boost, step=5,
            help="Hedef inhibe oldukça paralog genlerin aktivitesini artırma yüzdesi."
        )

        bc_accuracy = st.selectbox(
            ui_t("betweenness_accuracy"),
            options=["fast", "balanced", "full"],
            index=(2 if st.session_state.bc_sample_sources is None
                   else 1 if st.session_state.bc_sample_sources == 3000 else 0),
            format_func=lambda value: {"fast": ui_t("fast"), "balanced": ui_t("balanced"), "full": ui_t("full")}[value],
            help="Hızlı/Dengeli modlar sabit tohumlu kaynak örneklemesi kullanır. Tam mod, nihai raporlar için tüm düğümlerle hesaplar.",
        )
        st.session_state.bc_sample_sources = {"fast": 1200, "balanced": 3000, "full": None}[bc_accuracy]

        st.session_state.harita_max_node = st.slider(
            f"🌐 {ui_t('max_network_nodes')}",
            min_value=150, max_value=1500, value=st.session_state.harita_max_node, step=50,
            help="CSI Ağ Haritası'nda 'Gerçek Ağ' modu seçildiğinde tarayıcıya gönderilecek "
                 "en fazla düğüm sayısı. Kilitlenen ve stresli genler bu tavandan bağımsız "
                 "olarak HER ZAMAN dahil edilir; tavan yalnızca eklenecek bağlam (komşu/hub) "
                 "düğümlerini sınırlar. Daha düşük değer daha okunaklı, daha yüksek değer "
                 "daha kapsamlı ama daha yoğun bir harita üretir."
        )
        st.session_state.harita_max_etiket = st.slider(
            f"🏷️ {ui_t('max_network_labels')}",
            min_value=10, max_value=150, value=st.session_state.harita_max_etiket, step=5,
            help="Yalnızca bu kadar en yüksek skorlu kilitlenen/stresli gene metin "
                 "etiketi basılır; geri kalanlar üzerine gelince (hover) görünür. "
                 "Bu, yüzlerce üst üste yazının okunamaz bir yığın oluşturmasını önler."
        )

        if st.button(ui_t("restore_defaults")):
            st.session_state.damping = 0.85
            st.session_state.block_strength = 0.001
            st.session_state.hill_n = 2.0
            st.session_state.paralog_boost = 10
            st.session_state.bc_sample_sources = 1200
            st.session_state.harita_max_node = 600
            st.session_state.harita_max_etiket = 40
            st.rerun()

            st.divider()
        with st.expander("🧭 Parametre Kılavuzu (Nasıl Kullanılır?)", expanded=False):
            st.markdown(
                '<div style="font-size:0.75rem; color:#b8b8bd; line-height:1.6;">'
                '<p><strong style="color:#ffffff;">Ağ önemi yayılımı (α)</strong><br>'
                'Sinyalin ağda ne kadar uzağa yayılacağını belirler. <br>'
                '• <strong style="color:#4ade80;">0.85 (Varsayılan):</strong> Standart davranış.<br>'
                '• <strong style="color:#60a5fa;">0.95 - 0.99:</strong> Sinyal tüm ağa yayılır, global etkileri görmek için idealdir.<br>'
                '• <strong style="color:#fbbf24;">0.70 - 0.80:</strong> Sinyal hedef gene yakın kalır, sadece yerel komşuların stresini test etmek için kullanılır.</p>'
                '<p><strong style="color:#ffffff;">Susturma Şiddeti</strong><br>'
                'Hedef genin bağlantılarının ne kadar zayıflatılacağını belirler.<br>'
                '• <strong style="color:#4ade80;">0.0010 (Varsayılan):</strong> Başlangıç kenar ağırlığının %0.1’i kalır; deneysel knockout değildir.<br>'
                '• <strong style="color:#60a5fa;">0.0100 - 0.0500:</strong> Kısmi inhibisyon (ilaç etkisi), doz-yanıt testleri için uygundur.<br>'
                '• <strong style="color:#fbbf24;">0.1000:</strong> Başlangıç kenar ağırlığının %10’u kalır; biyolojik aktivite ölçümü değildir.</p>'
                '<p><strong style="color:#ffffff;">Hill Katsayısı (n)</strong><br>'
                'Doz-yanıt eğrisinin dikliğini kontrol eder.<br>'
                '• <strong style="color:#4ade80;">2.0 (Varsayılan):</strong> Klasik sigmoid eğri, çoğu biyolojik sistem için uygundur.<br>'
                '• <strong style="color:#60a5fa;">0.5 - 1.5:</strong> Daha yumuşak, kademeli etki (ilaç etkisi yavaş başlar).<br>'
                '• <strong style="color:#fbbf24;">3.0 - 4.0:</strong> Çok dik eğri (aşırı hassas sistemler veya işbirlikçi bağlanma).</p>'
                '<p><strong style="color:#ffffff;">Paralog Boost (%)</strong><br>'
                'Hedef inhibe olduğunda, benzer işlevli yedek genlerin aktivitesini ne kadar artıracağını simüle eder.<br>'
                '• <strong style="color:#4ade80;">%10 (Varsayılan):</strong> Ilımlı bir telafi mekanizması.<br>'
                '• <strong style="color:#60a5fa;">%0:</strong> Modelde paralog katkısının artırılmadığı senaryo.<br>'
                '• <strong style="color:#fbbf24;">%25 - 50:</strong> Modelde daha yüksek paralog katkısı; ilaç direncini kanıtlamaz.</p>'
                '</div>',
                unsafe_allow_html=True,
            )

        st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 0.8rem 0 0.4rem 0;"></div>', unsafe_allow_html=True)

        _aktif_tur = "Fare" if IS_MOUSE else "İnsan"
        _tur_gecmisi = [k for k in st.session_state.get("sim_gecmis", []) if k.get("tur", _aktif_tur) == _aktif_tur]
        if _tur_gecmisi:
            st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 0.8rem 0 0.4rem 0;"></div>', unsafe_allow_html=True)
            with st.expander("📜 Simülasyon Geçmişi (son 10)", expanded=False):
                for kayit in reversed(_tur_gecmisi):
                    st.markdown(f"""
                    <div style="background: #151517; border: 1px solid #2a2a2e; border-radius: 6px; padding: 0.5rem; margin-bottom: 0.3rem;">
                        <span style="color: #f5f5f5; font-weight: 700;">{kayit['hedef']}</span><br>
                        <span style="color: #85858c; font-size: 0.7rem;">
                            {kayit['doku']} · {kayit['stresli_sayisi']} stresli gen · Sistemik ağ kayması: %{kayit['kayip_pct']:.4f} · {kayit['tarih']}
                        </span>
                    </div>
                    """, unsafe_allow_html=True)
        if st.session_state.get("ui_developer_mode", False):
            st.caption(str(BASE_DIR))

    ates_buton = _config_run_surface.button(ui_t("start_analysis"), type="primary", width="stretch", disabled=workspace_route != "Perturbation" or (not bool(hedef_genler_sniper) and mod != "localization"), help="Runs the analysis with the selected target, tissue, and parameters." if ui_locale() == "en" else "Seçili hedef/doku ve yukarıdaki parametrelerle hesaplar.")

    with _config_tools_surface.expander(ui_t("target_tools"), expanded=False):
        st.markdown('<div class="sidebar-label">Tools · Target Stress</div>', unsafe_allow_html=True)
        with st.expander("Target Stress Search nasıl okunur?", expanded=False):
            st.markdown("""
    **Bu araç neyi yanıtlar?** Seçtiğiniz tek ana genin ağ içindeki konumunu hangi aday müdahalelerin en fazla değiştirebileceğini araştırır. Araç, aday genleri tek tek perturb eder; ardından ana genin ve yakın ağ çevresinin nasıl yanıt verdiğini karşılaştırır. Bu nedenle sonuç, “hangi gen ana hedefle ilişkilidir?” listesinden daha fazlasıdır: model içindeki duyarlılık sırasıdır.

    **Etki modu seçimi:** **Ağ-Aracılı Etki**, protein ağının bağlantı yapısından gelen sinyali ve forward-simulation doğrulamasını kullanır. **Doğrudan Etki**, yüklü TRRUST yönlü düzenleyici kanıtı varsa buna odaklanır. **Birleşik**, iki kanıt tipini yan yana sunar. Doğrudan kanıt bulunmaması, biyolojik ilişkinin olmadığı anlamına gelmez; sadece bu veri kaynağında yönlü kayıt olmadığı anlamına gelir.

    - **Ağ-Aracılı Etki:** Aday ile ana hedef arasındaki ağ-topolojisi üzerinden ortaya çıkan hesaplamalı etkiyi inceler. Doğrudan düzenleyici ilişki gerektirmez.
    - **Doğrudan Etki:** Yalnızca yüklü yönlü TRRUST düzenleyici kanıtıyla desteklenen aday-hedef ilişkilerini öne çıkarır. Kanıtın olmaması ilişkinin biyolojik olarak bulunmadığını göstermez.
    - **Birleşik:** Ağ-aracılı sonuçları ve mevcut doğrudan yönlü kanıtı aynı değerlendirmede sunar; iki kanıt türünü birbirinin yerine geçecek biçimde yorumlamaz.

    **Aday limiti:** Her aday için gerçek bir forward simülasyon çalıştırıldığı için limit arttıkça süre uzar. Küçük limit, ilk inceleme için; daha büyük limit, daha geniş aday taraması içindir. Limit yalnızca kaç adayın test edildiğini değiştirir, hiçbir ağ formülünü değiştirmez.

    **Sonuç tablosunu okuma:** **Ağ önemi değişimi**, aday pertürbasyonu sonrasında ilgili genin ağdaki göreli yapısal önemindeki değişimdir. Pozitif değer göreli ağ payının arttığını, negatif değer azaldığını söyler. Bu değer ifade düzeyi, protein miktarı, ilaç etkisi veya biyolojik aktivasyon değildir. **Yerel ΔBC** ana genin iki-hop çevresindeki geçiş noktası rolünün değişimidir; **Tam ΔBC** daha geniş ağ doğrulamasıdır. **Mesafe**, adayın ana gene ağdaki en kısa yol uzaklığıdır. Bu alanların hiçbiri tek başına nedensel biyolojik kanıt veya klinik öneri oluşturmaz.

    **Ağdaki Önemi:** Başlangıç ağındaki yapısal önem bağlamını verir. Yüksek olması Target Stress yanıtının otomatik olarak yüksek olacağını göstermez.
            """)
        target_stress_mode = st.selectbox(
            ui_t("effect_mode"), ["Ağ-Aracılı Etki", "Doğrudan Etki", "Birleşik"],
            key="target_stress_mode", label_visibility="collapsed",
            format_func=lambda value: {"Ağ-Aracılı Etki": ui_t("network_mediated_effect"), "Doğrudan Etki": ui_t("direct_effect"), "Birleşik": ui_t("combined")}[value],
            help="Doğrudan Etki yalnızca yüklü TRRUST yönlü kanıtını kullanır. Ağ-Aracılı Etki, ağ topolojisi ve forward simülasyon doğrulamasına dayanır.",
        )
        target_stress_limit = st.select_slider(
            ui_t("candidate_limit"), options=[5, 10, 15], value=5, key="target_stress_limit",
            help="Her aday gerçek forward simülasyonundan geçirilir; düşük değerler daha hızlı tamamlanır.",
        )
        target_stress_buton = st.button(
            "Target Stress Search", width="stretch", key="target_stress_start",
            disabled=not bool(hedef_genler_sniper and len(hedef_genler_sniper) == 1),
            help="Tam olarak bir hedef gen seçin. Bu gen doğrudan perturb edilmez; aday müdahaleler sonrası gözlenir.",
        )

        st.markdown('<div class="sidebar-label">Tools · Predicted Compensation · BETA</div>', unsafe_allow_html=True)
        compensation_buton = st.button(
            "Redistribution Candidates", width="stretch", key="compensation_start",
            disabled=not bool(hedef_genler_sniper and len(hedef_genler_sniper) == 1),
            help="Seçilen gen pertürbe edilir; ağ önemi artan genler hesaplamalı telafi adayları olarak raporlanır.",
        )

    hard_reset_buton = False
    sweep_ates = False
    with _config_system_surface.expander(ui_t("system"), expanded=False):
        st.selectbox(
            ui_t("language"), options=["en", "tr"],
            format_func=lambda value: "Türkçe" if value == "tr" else "English",
            key="ui_locale",
        )
        with st.expander(ui_t("local_resource_usage"), expanded=False):
            if PSUTIL_OK:
                proc = psutil.Process(os.getpid())
                ram_mb = proc.memory_info().rss / 1_048_576
                sys_ram = psutil.virtual_memory()
                pct = ram_mb / (sys_ram.total / 1_048_576) * 100
                st.caption(f"RAM · {ram_mb:.0f} MB · sistemin %{pct:.1f}'i")
        with st.expander(ui_t("maintenance_tools"), expanded=False):
            sweep_ates = st.button("Eşik taramasını başlat", width="stretch")
            hard_reset_buton = st.button(
                "Hard Reset", width="stretch",
                help="Tüm cache, session state ve RAM'i temizler.", key="hard_reset_btn",
            )


if hard_reset_buton:
    purged, collected = scorched_earth_reset()
    st.toast(f"Hard Reset tamamlandı — {purged} state, {collected} GC nesnesi.", icon="☢")
    time.sleep(0.8)
    st.rerun()


doku_label = hedef_doku if doku_aktif else ui_t("all_proteome")
_current_result_context = {
    "species": "Fare" if IS_MOUSE else "İnsan",
    "tissue": hedef_doku,
    "engine": calculation_engine,
    "targets": hedef_genler_sniper,
    "localization": lokalizasyon,
    "block_strength": float(st.session_state.block_strength),
    "damping": float(st.session_state.damping),
    "candidate_limit": int(st.session_state.get("redistribution_candidate_limit", 100)),
    "bc_sample_sources": st.session_state.bc_sample_sources,
    "tissue_normalization_mode": st.session_state.get("tissue_normalization_mode", "within_tissue"),
}
_has_matching_completed_result = (
    st.session_state.get("rapor_df") is not None
    and app_state.simulation_result_context_matches(
        st.session_state.get("simulation_result_context"), **_current_result_context,
    )
)
if hedef_genler_sniper:
    _analysis_context_label = ui_t("target_gene_mode")
    _analysis_context_value = " + ".join(
        _gene_to_symbol.get(gene, gene) for gene in hedef_genler_sniper
    )
else:
    if mod == "localization":
        _analysis_context_label = ui_t("localization_mode")
        _analysis_context_value = str(lokalizasyon).replace("_", " ")
    else:
        _analysis_context_label = "Hedef seçimi"
        _analysis_context_value = "Henüz hedef seçilmedi"

if _stale_detected:
    st.markdown(
        '<div class="stale-banner">⚡ Yeni simülasyon verisi tespit edildi — cache yenilendi.</div>',
        unsafe_allow_html=True,
    )

cache_key = tuple(hedef_genler_sniper) if hedef_genler_sniper else None

if workspace_route == "Unified & records":
    from src.ui.unified_workspace import render_unified_workspace
    render_unified_workspace(targets=hedef_genler_sniper, tissue=hedef_doku, is_mouse=IS_MOUSE, output_dir=OUTPUTS_DIR)
    st.stop()

if not hedef_genler_sniper and mod != "localization":
    render_empty_workspace()
    _example_option = next((item for item in _secenekler if item.startswith("CFTR  ·")), None)
    if _example_option:
        def _prepare_cftr_example() -> None:
            st.session_state["gen_sec"] = [_example_option]
            st.session_state["analysis_tissue"] = "Lung"

        _example_col, _guide_col = st.columns([.35, .65], vertical_alignment="center")
        with _example_col:
            st.button(
                "Prepare CFTR · Lung example" if ui_locale() == "en" else "CFTR · Lung örneğini hazırla", type="primary", width="stretch",
                on_click=_prepare_cftr_example, key="prepare_cftr_example",
            )
        with _guide_col:
            st.caption("The example prepares selections only; you start the calculation. Engine and output rules do not change." if ui_locale() == "en" else "Örnek yalnız seçimleri hazırlar; hesaplamayı siz başlatırsınız. Motor ve çıktı kuralları değişmez.")
    with st.expander("Çalışma alanını nasıl kullanırım?", expanded=False):
        render_reading_guide()
    st.stop()

# Invalidate a completed result before any potentially long network rebuild.
# The run button remains usable: on its click rerun we continue below with the
# newly selected context, while every ordinary selection-change rerun stops
# here so the previous result can never remain visible under the new header.
_result_is_stale = (
    st.session_state.get("rapor_df") is not None
    and not app_state.simulation_result_context_matches(
        st.session_state.get("simulation_result_context"),
        **_current_result_context,
    )
)
if _result_is_stale:
    st.markdown("## Seçimler değişti")
    st.caption("Mevcut sonuç eski bağlama ait olduğu için gizlendi.")
    if st.button("Yeni analizi çalıştır", type="primary", key="run_stale_analysis"):
        ates_buton = True
    if not ates_buton:
        st.stop()

_network_status_slot = st.empty()
_network_loading_slot = st.empty()
_network_loading_slot.markdown(
    ('<section class="sk-loading-workspace"><div><span class="sk-loading-pulse"></span><b>Preparing tissue network</b>'
     '<p>STRING interactions, the tissue filter, and baseline topology are being combined in one analysis context.</p></div>'
     '<ol><li class="active">Network data</li><li>Tissue context</li><li>Centrality</li><li>Ready to analyse</li></ol></section>')
    if ui_locale() == "en" else
    ('<section class="sk-loading-workspace"><div><span class="sk-loading-pulse"></span><b>Doku ağı hazırlanıyor</b>'
     '<p>STRING etkileşimleri, doku filtresi ve başlangıç topolojisi aynı analiz bağlamında birleştiriliyor.</p></div>'
     '<ol><li class="active">Ağ verisi</li><li>Doku bağlamı</li><li>Merkezilik</li><li>Analize hazır</li></ol></section>'),
    unsafe_allow_html=True,
)
with _network_status_slot.status("Preparing graph..." if ui_locale() == "en" else "Graf hazırlanıyor...", expanded=False) as _motor_status:
    _pipeline_log("NETWORK", "Ağ hazırlığı başladı (tissue=%s)", hedef_doku if doku_aktif else "None")
    _motor_status.write("▸ Reading raw network data (STRING)..." if ui_locale() == "en" else "▸ Ham ağ verisi (StringDB) okunuyor...")
    time.sleep(0.05)
    _motor_status.write("▸ Applying the tissue-expression filter..." if doku_aktif and ui_locale() == "en" else ("▸ Doku ifade filtresi uygulanıyor..." if doku_aktif else ("▸ Tissue filter off — using the full proteome" if ui_locale() == "en" else "▸ Doku filtresi kapalı — tüm proteom kullanılıyor")))
    _motor_status.write("▸ Calculating Network Importance / betweenness centrality (skipped when cached)..." if ui_locale() == "en" else "▸ Ağdaki Önemi / Geçiş Merkeziliği hesaplanıyor (cache varsa atlanır)...")
    _prepared = prepare_network(
        motor, cache_key, hedef_doku, st.session_state.bc_sample_sources,
        cache_identity=(evidence_motor.cache_identity() if calculation_engine == "Evidence" else None),
    )
    g_canli, df_canli = _prepared.graph, _prepared.scores
    _pipeline_log(
        "NETWORK", "Ağ hazır (nodes=%d edges=%d)",
        g_canli.vcount() if g_canli is not None else 0,
        g_canli.ecount() if g_canli is not None else 0,
    )
    _motor_status.write((f"▸ Complete — {g_canli.vcount() if g_canli is not None else 0:,} nodes, {g_canli.ecount() if g_canli is not None else 0:,} edges ready.") if ui_locale() == "en" else (f"▸ Tamamlandı — {g_canli.vcount() if g_canli is not None else 0:,} düğüm, {g_canli.ecount() if g_canli is not None else 0:,} kenar hazır."))
    _motor_status.update(label="✓ Graph ready" if ui_locale() == "en" else "✓ Graf hazır", state="complete", expanded=False)
_network_loading_slot.empty()

if _has_matching_completed_result:
    _network_status_slot.empty()

st.session_state["G_canli"] = g_canli

if hedef_genler_sniper and g_canli is not None:
    eksik = [gn for gn in hedef_genler_sniper if gn not in set(g_canli.vs["name"])]
    if eksik:
        _config_panel.warning(
            f"⚠ {len(eksik)} hedef seçili doku/ağ içinde bulunamadı; "
            "simülasyona dahil edilmedi."
        )


if sweep_ates and not IS_MOUSE:
    with st.spinner("Threshold Sweep çalışıyor..."):
        try:
            motor.run_threshold_sweep(str(motor.DB_PATH))
            st.success("Threshold Sweep tamamlandı.")
        except Exception as e:
            st.error(f"Hata: {e}")
elif sweep_ates and IS_MOUSE:
    st.warning("⚠ Fare modunda Threshold Sweep henüz desteklenmiyor.")


_tool_inputs = {
    "targets": tuple(hedef_genler_sniper or ()), "tissue": hedef_doku,
    "species": "Fare" if IS_MOUSE else "İnsan", "engine": calculation_engine,
    "block_strength": st.session_state.block_strength, "damping": st.session_state.damping,
    "bc_sample_sources": st.session_state.bc_sample_sources,
}
# Tool results are never presented under a changed scientific configuration.
for _tool_name in ("compensation", "target_stress"):
    _saved_tool = st.session_state.get(f"{_tool_name}_context") or {}
    _expected_inputs = dict(_tool_inputs)
    if _tool_name == "target_stress":
        _expected_inputs.update(mode=target_stress_mode, limit=target_stress_limit)
    if st.session_state.get(f"{_tool_name}_results") is not None and _saved_tool.get("inputs") != _expected_inputs:
        st.session_state[f"{_tool_name}_results"] = None
        st.session_state[f"{_tool_name}_context"] = None
        st.info("Hedef veya analiz parametreleri değişti. Önceki keşif sonucu gizlendi; aracı yeniden çalıştırın.")

if compensation_buton:
    perturbed_gene = hedef_genler_sniper[0]
    perturbed_symbol = _gene_to_symbol.get(perturbed_gene, perturbed_gene)
    with st.status(f"Compensation Analysis çalışıyor: {perturbed_symbol}", expanded=True) as _comp_status:
        try:
            _compensation_base_engine = base_engine_for_mode(calculation_engine)
            _pipeline_log(
                _compensation_base_engine.upper(),
                "%s temel simülasyon başladı (target=%s)",
                _compensation_base_engine.title(),
                perturbed_symbol,
            )
            _comp_status.write("Forward pertürbasyon çalıştırılıyor; ağ önemi artışları ölçülüyor...")
            _comp_results = predicted_compensation_candidates(
                motor_module=motor, graph=g_canli, scores=df_canli, perturbed_gene=perturbed_gene,
                output_root=OUTPUTS_DIR, block_strength=st.session_state.block_strength,
                damping=st.session_state.damping, ghost_filter=ghost_gen_filtrele,
                gene_to_symbol=_gene_to_symbol,
            )
            st.session_state["compensation_results"] = _comp_results
            st.session_state["compensation_context"] = {
                "target": perturbed_symbol, "tissue": doku_label, "species": "Fare" if IS_MOUSE else "İnsan",
                "inputs": dict(_tool_inputs)
            }
            _comp_status.update(label=f"Compensation Analysis tamamlandı: {len(_comp_results)} aday", state="complete", expanded=False)
        except Exception as error:
            _comp_status.update(label="Compensation Analysis hatası", state="error", expanded=True)
            user_error("Compensation Analysis tamamlanamadı. Lütfen hedefi ve doku seçimini kontrol edip yeniden deneyin.", exception=error)

if target_stress_buton:
    st.session_state["compensation_results"] = None
    st.session_state["compensation_context"] = None
    target_gene = hedef_genler_sniper[0]
    target_symbol = _gene_to_symbol.get(target_gene, target_gene)
    with st.status(f"Target Stress Search çalışıyor: {target_symbol}", expanded=True) as _target_status:
        try:
            _target_status.write("Aday ön seçimi: topoloji, topluluk, merkezilik, WBI ve yüklü yönlü kanıtlar değerlendiriliyor...")
            candidates = preselect_candidates(
                graph=g_canli, scores=df_canli, target_gene=target_gene,
                regulators=TARGET_TO_REGULATORS, gene_to_symbol=_gene_to_symbol,
                mode=target_stress_mode, limit=int(target_stress_limit),
            )
            if candidates.empty:
                st.session_state["target_stress_results"] = pd.DataFrame()
                st.session_state["target_stress_context"] = {
                    "inputs": {**_tool_inputs, "mode": target_stress_mode, "limit": target_stress_limit},
                    "target": target_symbol, "tissue": doku_label, "species": "Fare" if IS_MOUSE else "İnsan"
                }
                _target_status.update(label="Filtreyi geçen aday bulunamadı", state="complete", expanded=False)
            else:
                _target_status.write(f"{len(candidates)} aday mevcut forward simülasyonu ile doğrulanıyor...")
                results = validate_candidates(
                    candidates=candidates, motor_module=motor, graph=g_canli, scores=df_canli,
                    target_gene=target_gene, output_root=OUTPUTS_DIR,
                    block_strength=st.session_state.block_strength, damping=st.session_state.damping,
                    ghost_filter=ghost_gen_filtrele,
                )
                st.session_state["target_stress_results"] = results
                st.session_state["target_stress_context"] = {
                    "inputs": {**_tool_inputs, "mode": target_stress_mode, "limit": target_stress_limit},
                    "target": target_symbol, "tissue": doku_label, "species": "Fare" if IS_MOUSE else "İnsan",
                    "mode": target_stress_mode,
                }
                _target_status.update(label=f"Target Stress Search tamamlandı: {len(results)} doğrulanmış aday", state="complete", expanded=False)
        except Exception as error:
            _target_status.update(label="Target Stress Search hatası", state="error", expanded=True)
            user_error("Target Stress Search tamamlanamadı. Lütfen hedefi ve doku seçimini kontrol edip yeniden deneyin.", exception=error)

ates_buton = bool(ates_buton or st.session_state.pop("_run_from_ready", False))
if ates_buton:
    st.session_state["target_stress_results"] = None
    st.session_state["target_stress_context"] = None
    st.session_state["compensation_results"] = None
    st.session_state["compensation_context"] = None
    app_state.reset_simulation_result()
    st.cache_data.clear()
    force_gc()
    # Freeze the exploratory N before the scientific call. Later sidebar
    # changes must not rewrite the scope of this completed simulation.
    _simulation_exploratory_top_n = int(
        st.session_state.get("redistribution_candidate_limit", 100)
    )
    _simulation_result_context = app_state.build_simulation_result_context(
        species="Fare" if IS_MOUSE else "İnsan",
        tissue=hedef_doku,
        engine=calculation_engine,
        targets=hedef_genler_sniper,
        localization=lokalizasyon,
        block_strength=float(st.session_state.block_strength),
        damping=float(st.session_state.damping),
        candidate_limit=_simulation_exploratory_top_n,
        bc_sample_sources=st.session_state.bc_sample_sources,
        tissue_normalization_mode=st.session_state.get(
            "tissue_normalization_mode", "within_tissue"
        ),
    )

    hedef_label = (
        lokalizasyon if not hedef_genler_sniper
        else " + ".join(_gene_to_symbol.get(gene, gene) for gene in hedef_genler_sniper)
    )
    log_satirlari = []


    if hedef_genler_sniper:
        hedef_ensp_list = hedef_genler_sniper
        hedef_symbols = [_gene_to_symbol.get(g, g) for g in hedef_ensp_list]

        hedef_info_df = df_canli[df_canli["gene"].isin(hedef_ensp_list)].copy()
        if not hedef_info_df.empty:
            with st.container():
                st.markdown("### 🎯 Simulated target" if ui_locale() == "en" else "### 🎯 Simüle Edilen Hedef")
                cols = st.columns(min(len(hedef_info_df), 3))
                for idx, (_, row) in enumerate(hedef_info_df.iterrows()):
                    col = cols[idx % len(cols)]
                    with col:
                        sym = _gene_to_symbol.get(str(row["gene"])) or _harita_aliasi(str(row["gene"]), _gene_to_symbol) or "Alias bulunamadı"
                        skor = row.get("Hinterland_Skoru", "—")
                        lok = str(row.get("Lokalizasyon", "—")).replace("_", " ")
                        st.markdown(f"""
                        <div class="sk-running-target">
                            <div class="sk-running-target-alias">{sym}</div>
                            <div class="sk-running-target-id">{row["gene"]}</div>
                            <div class="sk-running-target-meta">
                                <span>{'Network importance' if ui_locale() == 'en' else 'Ağdaki Önemi'} <b>{skor:.1f}</b></span>
                                <span>{'Cellular location' if ui_locale() == 'en' else 'Hücresel konum'} <b>{lok}</b></span>
                            </div>
                        </div>
                        """, unsafe_allow_html=True)
                st.markdown(f"**{'Tissue' if ui_locale() == 'en' else 'Doku'}:** {doku_label}")
                st.markdown(f"**{'Target' if ui_locale() == 'en' else 'Hedef'}:** {len(hedef_ensp_list)}")
        else:
            st.info(f"Hedef genler ({', '.join(hedef_symbols)}) veri tablosunda bulunamadı.")
    else:
        st.markdown("### 🎯 Simulation target" if ui_locale() == "en" else "### 🎯 Simülasyon Hedefi")
        st.markdown(f"**{'Localization' if ui_locale() == 'en' else 'Lokalizasyon'}:** {str(lokalizasyon).replace('_', ' ')}")
        st.markdown(f"**{'Tissue' if ui_locale() == 'en' else 'Doku'}:** {doku_label}")
        if "Lokalizasyon" in df_canli.columns:
            lok_count = df_canli["Lokalizasyon"].str.contains(lokalizasyon, na=False).sum()
            st.markdown(
                f"Approximately **{lok_count}** genes will be targeted in this region."
                if ui_locale() == "en" else f"Bu bölgede yaklaşık **{lok_count}** gen hedeflenecek."
            )

    with st.status(f"{'Simulation running' if ui_locale() == 'en' else 'Simülasyon çalışıyor'} — {hedef_label}", expanded=True) as _sim_status:
        try:
            _sim_status.write("▸ Applying the edge-attenuation scenario to the network..." if ui_locale() == "en" else "▸ Susturma (knockout) senaryosu ağa uygulanıyor...")
            kw = (dict(spesifik_hedefler=hedef_genler_sniper)
                  if hedef_genler_sniper else dict(hedef_lokalizasyon=lokalizasyon))
            log_satirlari.append(f"run_infection_simulation çağrılıyor, kw={kw}")
            _sim_status.write("▸ Propagating the modeled network response (Network Importance + betweenness calculation)..." if ui_locale() == "en" else "▸ Enfeksiyon / şok dalgası ağda yayılıyor (ağ önemi + Geçiş Merkeziliği hesabı)...")
            sonuc_df, hayalet = execute_simulation(
                motor_module=motor, graph=g_canli, scores=df_canli, output_dir=OUTPUTS_DIR,
                targets=hedef_genler_sniper, localization=lokalizasyon,
                block_strength=st.session_state.block_strength, damping=st.session_state.damping,
                scientific_options={
                    "exploratory_top_n": _simulation_exploratory_top_n,
                    # This is provenance supplied to the existing motor; it
                    # neither changes graph construction nor any score.
                    "tissue": hedef_doku if doku_aktif else "None",
                    "tissue_normalization_mode": st.session_state.get(
                        "tissue_normalization_mode", "within_tissue"
                    ),
                },
                apply_druggability=(lambda result: result if IS_MOUSE else compute_druggable_gate_score(result, target_type="stressed")),
                ghost_filter=ghost_gen_filtrele,
                # The signed frame is a separate, already-computed scientific
                # result. Keep it in the session alongside the presentation
                # report so the loss view never falls back to a missing file.
                signed_result_sink=lambda frame: st.session_state.__setitem__(
                    "signed_redistribution_result", frame.copy(deep=True)
                ),
            )
            _base_engine = base_engine_for_mode(calculation_engine)
            _pipeline_log(
                _base_engine.upper(),
                "%s temel simülasyon tamamlandı (rows=%s)",
                _base_engine.title(),
                getattr(sonuc_df, "shape", None),
            )
            log_satirlari.append(f"motor döndü → boyut={getattr(sonuc_df, 'shape', 'N/A')}")
            if IS_MOUSE:
                _sim_status.write("▸ Fare modunda ilaçlanabilirlik gate-skoru atlanıyor (yapısal veri yok).")
            else:
                _sim_status.write("▸ Stresli/kilitlenen genler sınıflandırıldı, ilaçlanabilirlik skoru hesaplandı.")
            if sonuc_df is not None and isinstance(sonuc_df, pd.DataFrame) and not sonuc_df.empty:
                sonuc_df_filtered = sonuc_df
                _simulation_selection_mode = _research_selection_mode_from_report(sonuc_df_filtered)
                _simulation_signed = st.session_state.get("signed_redistribution_result")
                if isinstance(_simulation_signed, pd.DataFrame) and not _simulation_signed.empty:
                    _cross_species_signed = _simulation_signed.copy(deep=True)
                    if "Symbol" not in _cross_species_signed.columns and "gene" in _cross_species_signed.columns:
                        _cross_species_signed["Symbol"] = (
                            _cross_species_signed["gene"].astype(str).map(_gene_to_symbol)
                        )
                    st.session_state[
                        "cross_species_mouse_snapshot" if IS_MOUSE else "cross_species_human_snapshot"
                    ] = {
                        "response": _cross_species_signed,
                        "targets": tuple(map(str, hedef_genler_sniper or ())),
                        "tissue": hedef_doku if doku_aktif else "None",
                        "config": ComparisonConfig(
                            string_threshold=int(motor.HIGH_CONF_THRESHOLD),
                            damping=float(st.session_state.damping),
                            attenuation_fraction=float(st.session_state.block_strength),
                            selection_mode=_simulation_selection_mode,
                            algorithm_version="sophiark-classic-v5.3",
                        ),
                    }
                _simulation_targets = tuple(
                    sonuc_df_filtered.loc[
                        sonuc_df_filtered["Hasar_Tipi"] == motor.HASAR_BIRINCIL, "gene"
                    ].dropna().astype(str)
                )
                if calculation_engine in {"Directed", "Compare"}:
                    _sim_status.write("▸ OmniPath yönleri STRING topolojisine uygulanıyor; yönlü motor ayrı çalışıyor...")
                    _directed_bundle = run_optional_directed_engine(
                        graph=g_canli,
                        targets=_simulation_targets,
                        tissue=hedef_doku if doku_aktif else "None",
                        classic_report=(
                            _simulation_signed if isinstance(_simulation_signed, pd.DataFrame)
                            else sonuc_df_filtered
                        ),
                        project_root=Path(__file__).resolve().parent,
                        mode=calculation_engine,
                        block_weight_fraction=float(st.session_state.block_strength),
                        damping=float(st.session_state.damping),
                        bc_sample_sources=min(int(st.session_state.get("bc_sample_sources", 1200) or 4), 4),
                        taxon_id=10090 if IS_MOUSE else 9606,
                        scores=df_canli,
                        gene_to_symbol=_gene_to_symbol,
                        candidate_limit=_simulation_exploratory_top_n,
                    )
                    st.session_state["directed_run_bundle"] = _directed_bundle
                    st.session_state["directed_engine_mode"] = calculation_engine
                else:
                    st.session_state.pop("directed_run_bundle", None)
                    st.session_state["directed_engine_mode"] = calculation_engine
                _simulation_returned = sonuc_df_filtered.loc[
                    sonuc_df_filtered["Hasar_Tipi"] == motor.HASAR_UCUNCUL
                ]
                _simulation_returned_count = _research_unique_gene_count(
                    _simulation_returned
                )
                app_state.store_research_simulation_scope(
                    species="Mus musculus" if IS_MOUSE else "Homo sapiens",
                    taxon_id=10090 if IS_MOUSE else 9606,
                    tissue=doku_label,
                    targets=_simulation_targets,
                    attenuation=float(st.session_state.block_strength),
                    selection_mode=_simulation_selection_mode,
                    test_limit=None,
                    tested_count=(
                        _research_unique_gene_count(_simulation_signed)
                        if _simulation_selection_mode.casefold() == "top_n"
                        else None
                    ),
                    returned_count=int(_simulation_returned_count or 0),
                    top_n=(
                        _simulation_exploratory_top_n
                        if _simulation_selection_mode.casefold() == "top_n"
                        else None
                    ),
                    threshold_parameters=_research_scientific_threshold_parameters(
                        _simulation_selection_mode
                    ),
                )
                app_state.bind_simulation_result(
                    report=sonuc_df_filtered,
                    context=_simulation_result_context,
                    ghosts=hayalet,
                    target_label=hedef_label,
                )
                _sim_status.write(f"▸ {len(sonuc_df_filtered)} satırlık rapor derlendi.")
                _sim_status.update(label=f"✓ Simülasyon tamamlandı — {len(sonuc_df_filtered)} satır", state="complete", expanded=False)
                st.success(f"✓ {len(sonuc_df_filtered)} satır rapor hazır.")

                n_stresli = len(sonuc_df_filtered[sonuc_df_filtered["Hasar_Tipi"] == motor.HASAR_UCUNCUL])
                kayip_pct = float(sonuc_df_filtered["Efficiency_Kayip_Pct"].iloc[0])
                app_state.append_history(
                    target=hedef_label, tissue=doku_label, species="Fare" if IS_MOUSE else "İnsan",
                    stressed_count=n_stresli, loss_pct=kayip_pct, timestamp=time.strftime("%H:%M:%S"),
                )
            else:
                _sim_status.write("▸ Bellekten sonuç alınamadı, diskteki rapor deneniyor...")
                time.sleep(5.5)
                disk_df = csv_yukle_debug(SOK_RAPORU)
                if disk_df is not None and not disk_df.empty:
                    disk_df_filtered, hayalet = ghost_gen_filtrele(disk_df, g_canli)
                    app_state.bind_simulation_result(
                        report=disk_df_filtered,
                        context=_simulation_result_context,
                        ghosts=hayalet,
                        target_label=hedef_label,
                    )
                    _sim_status.update(label=f"✓ Diskten {len(disk_df_filtered)} satır yüklendi", state="complete", expanded=False)
                    st.success(f"✓ {len(disk_df_filtered)} satır diskten yüklendi.")
                else:
                    _sim_status.update(label="✗ Veri alınamadı", state="error", expanded=True)
                    st.error("Veri alınamadı.")
        except Exception as e:
            log_satirlari.append(f"❌ Exception: {e}")
            _sim_status.update(label="✗ Simülasyon hatası", state="error", expanded=True)
            user_error("Simülasyon tamamlanamadı. Hedefin seçili ağda bulunduğunu kontrol edip yeniden deneyin.", exception=e)

# Forward simülasyon tamamlandığında KEGG + GO-BP otomatik üretilir. Bu adım
# sonuç ekranını bloklayan bir kullanıcı seçimi değildir; uzak kaynak sorunları
# yalnızca yolak bölümünde açıklanır.
if ates_buton and st.session_state.get("sim_tamam"):
    _auto_report = st.session_state.get("rapor_df")
    if isinstance(_auto_report, pd.DataFrame) and not _auto_report.empty:
        _auto_engine = base_engine_for_mode(calculation_engine)
        _auto_fingerprint = candidate_set_fingerprint(
            _auto_engine, _auto_report[_auto_report["Hasar_Tipi"] == motor.HASAR_UCUNCUL]["gene"].astype(str).tolist(),
        )
        _auto_symbols = [
            _gene_to_symbol.get(gene, "")
            for gene in _auto_report[_auto_report["Hasar_Tipi"] == motor.HASAR_UCUNCUL]["gene"].astype(str)
        ]
        _auto_bundle = st.session_state.get("directed_run_bundle")
        if calculation_engine == "Directed" and getattr(_auto_bundle, "presentation", None) is not None:
            _auto_engine = "directed"
            _auto_fingerprint = _auto_bundle.presentation.candidate_fingerprint
            _auto_symbols = _auto_bundle.presentation.candidates["gene_symbol"].astype(str).tolist()
        _auto_symbols = [symbol for symbol in _auto_symbols if symbol and symbol != UNRESOLVED_SYMBOL]
        # ORA arka planı, varsayılan genom değil, o anda analiz edilen aktif
        # (doku filtreli ise filtrelenmiş) ağın sembol evrenidir.
        _auto_network_ids = list(g_canli.vs["name"]) if g_canli is not None and "name" in g_canli.vs.attributes() else []
        _auto_background = [_gene_to_symbol.get(gene, "") for gene in _auto_network_ids]
        st.session_state["enrichment_engine"] = _auto_engine
        st.session_state["enrichment_candidate_fingerprint"] = _auto_fingerprint
        with st.spinner("KEGG ve GO Biological Process yolak analizi otomatik hazırlanıyor…"):
            try:
                _pipeline_log("ENRICHMENT", "%s aday kümesi analiz ediliyor (n=%d)", _auto_engine, len(_auto_symbols))
                _auto_enrichment, _auto_notices = run_post_simulation_enrichment(
                    _auto_symbols, is_mouse=IS_MOUSE, background_symbols=_auto_background,
                )
                st.session_state["enrichment_results"] = _auto_enrichment if not _auto_enrichment.empty else None
                st.session_state["enrichment_notices"] = _auto_notices
            except Exception as _auto_error:
                # Simülasyon sonucu geçerlidir; enrichment hatası onu geçersiz kılmaz.
                st.session_state["enrichment_results"] = None
                st.session_state["enrichment_notices"] = ["Yolak analizi şu anda tamamlanamadı; simülasyon sonucu korunuyor."]
                log.warning("Otomatik enrichment hatası: %s", _auto_error)
            # Teknik traceback kullanıcı arayüzünde gösterilmez; log tarafında kalır.
            force_gc()

    st.session_state["debug_log"] = log_satirlari

    st.rerun()


# Never auto-promote a legacy CSV to the active result. It has no trustworthy
# target/tissue/engine fingerprint and would immediately (and correctly) fail
# the stale-result guard. Same-run disk fallback is handled above, where the
# exact live context is still available and is bound atomically with the report.


if st.session_state.get("ui_developer_mode", False):
    with st.expander("🛠 Debug", expanded=False):
        st.markdown("**Dosya**")
        st.write(f"SOK_RAPORU: `{SOK_RAPORU.exists()}` | `{os.path.getsize(str(SOK_RAPORU)) if SOK_RAPORU.exists() else 'N/A'}` byte")
        st.markdown("**Symbol Tablosu**")
        st.write(f"Yüklenen sembol: `{len(symbol_tablo)}`")
        st.markdown("**Graf**")
        _g = st.session_state.get("G_canli")
        if _g is not None and isinstance(_g, ig.Graph):
            st.write(f"Düğüm: `{_g.vcount()}` | Kenar: `{_g.ecount()}` | Doku: `{hedef_doku}`")
        else:
            st.write("Graf yüklenmedi.")
        st.markdown("**Session**")
        st.write(f"rapor_df: `{getattr(st.session_state.get('rapor_df'), 'shape', 'None')}` | sim_tamam: `{st.session_state.get('sim_tamam', False)}`")
        hayalet_d = st.session_state.get("hayalet_genler", [])
        if hayalet_d:
            st.write(f"Hayalet: `{hayalet_d[:10]}{'…' if len(hayalet_d)>10 else ''}`")
        if st.session_state["debug_log"]:
            for satir in st.session_state["debug_log"]:
                st.markdown(f"- `{satir}`")
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            if st.button("CSV Zorla Yükle", key="force_reload"):
                zorlanan = csv_yukle_debug(SOK_RAPORU)
                if zorlanan is not None:
                    if g_canli is not None:
                        zorlanan, _ = ghost_gen_filtrele(zorlanan, g_canli)
                    st.session_state["rapor_df"] = zorlanan
                    st.session_state["simulation_result_context"] = None
                    st.session_state["research_simulation_scope"] = None
                    st.success(f"✓ {len(zorlanan)} satır")
                    st.rerun()
                else:
                    st.error("Okunamadı.")
        with col_d2:
            if st.button("GC Çalıştır", key="force_gc_btn"):
                n = force_gc()
                st.success(f"✓ {n} nesne")


df = st.session_state.get("rapor_df")
# 23:37 öncesi sürüm signed sonucu DataFrame.attrs içine koyuyordu. Açık
# oturumlarda kalan bu eski metadata Pandas'ın nlargest/concat attrs
# karşılaştırmasını bozabilir; sonuç değerlerine dokunmadan bir kez temizle.
if isinstance(df, pd.DataFrame) and "sophiark_signed_redistribution" in df.attrs:
    df = df.copy(deep=True)
    df.attrs.pop("sophiark_signed_redistribution", None)
    st.session_state["rapor_df"] = df
_target_stress_results = st.session_state.get("target_stress_results")
_compensation_results = st.session_state.get("compensation_results")
_propagation_trace = st.session_state.get("propagation_trace")
_tissue_differential = st.session_state.get("tissue_differential_result")
if _compensation_results is not None:
    _comp_context = st.session_state.get("compensation_context") or {}
    _close_compensation = render_compensation_results(
        _compensation_results, target_label=_comp_context.get("target", "Seçilen hedef"),
        tissue=_comp_context.get("tissue", doku_label),
        species=_comp_context.get("species", "Fare" if IS_MOUSE else "İnsan"),
        gene_to_symbol=_gene_to_symbol,
    )
    if _close_compensation:
        st.session_state["compensation_results"] = None
        st.session_state["compensation_context"] = None
        st.rerun()
    st.stop()
if _target_stress_results is not None:
    _target_context = st.session_state.get("target_stress_context") or {}
    _close_target_stress = render_target_stress_results(
        _target_stress_results,
        target_label=_target_context.get("target", "Selected target"),
        tissue=_target_context.get("tissue", doku_label),
        species=_target_context.get("species", "Fare" if IS_MOUSE else "İnsan"),
    )
    if _close_target_stress:
        st.session_state["target_stress_results"] = None
        st.session_state["target_stress_context"] = None
        st.rerun()
    st.stop()
if _propagation_trace is not None:
    _propagation_context = st.session_state.get("propagation_context") or {}
    _close_propagation = render_propagation_trace(
        _propagation_trace,
        target_label=_propagation_context.get("target", "Seçilen hedef"),
        tissue=_propagation_context.get("tissue", doku_label),
        species=_propagation_context.get("species", "Fare" if IS_MOUSE else "İnsan"),
    )
    if _close_propagation:
        st.session_state["propagation_trace"] = None
        st.session_state["propagation_context"] = None
        st.rerun()
    st.stop()
if _tissue_differential is not None:
    _tissue_context = st.session_state.get("tissue_differential_context") or {}
    _close_tissue_differential = render_tissue_differential(
        _tissue_differential,
        target_label=_tissue_context.get("target", "Seçilen hedef"),
        species=_tissue_context.get("species", "Fare" if IS_MOUSE else "İnsan"),
    )
    if _close_tissue_differential:
        st.session_state["tissue_differential_result"] = None
        st.session_state["tissue_differential_context"] = None
        st.rerun()
    st.stop()

hayalet_genler = st.session_state.get("hayalet_genler", [])
if hayalet_genler:
    st.markdown(
        f'<div class="ghost-banner">👻 {len(hayalet_genler)} gen "{hedef_doku}" '
        f'dokusunda üretilmediği için çıkarıldı.</div>',
        unsafe_allow_html=True,
    )

if df is None or (isinstance(df, pd.DataFrame) and df.empty):
    st.markdown("""
    <div class="sk-stepper">
        <div class="sk-step done"><span class="sk-step-num">✓</span><span class="sk-step-label">Yapılandır</span></div>
        <div class="sk-step-connector"></div>
        <div class="sk-step active"><span class="sk-step-num">2</span><span class="sk-step-label">Çalıştır</span></div>
        <div class="sk-step-connector"></div>
        <div class="sk-step todo"><span class="sk-step-num">3</span><span class="sk-step-label">İncele</span></div>
    </div>
    """, unsafe_allow_html=True)
    _ready_target = " + ".join(_gene_to_symbol.get(gene, gene) for gene in (hedef_genler_sniper or []))
    _ready_species = "Mus musculus" if IS_MOUSE else "Homo sapiens"
    st.markdown(
        '<section class="sk-ready-workspace">'
        '<div class="sk-ready-main"><span>ÇALIŞMAYA HAZIR</span>'
        f'<h2>{_ready_target or str(lokalizasyon).replace("_", " ")}</h2>'
        '<p>Hedef ve ağ bağlamı doğrulandı. Simülasyon henüz çalıştırılmadı.</p></div>'
        '<div class="sk-ready-context">'
        f'<div><small>ORGANİZMA</small><b>{_ready_species}</b></div>'
        f'<div><small>DOKU</small><b>{doku_label}</b></div>'
        f'<div><small>MOTOR</small><b>{calculation_engine}</b></div>'
        f'<div><small>HEDEF</small><b>{len(hedef_genler_sniper or []) or "Lokalizasyon"}</b></div>'
        '</div></section>',
        unsafe_allow_html=True,
    )

    def _queue_ready_analysis() -> None:
        st.session_state["_run_from_ready"] = True

    _ready_action, _ready_note = st.columns([.34, .66], vertical_alignment="center")
    with _ready_action:
        st.button(
            "Analizi şimdi çalıştır  →", type="primary", width="stretch",
            on_click=_queue_ready_analysis, key="run_ready_analysis",
        )
    with _ready_note:
        st.caption("Aynı Classic / Directed / Evidence çağrısı kullanılır; yalnızca eylem sayfa içine taşındı.")
    with st.expander("Nasıl çalışır?", expanded=False):
        render_reading_guide()
        render_mode_guide("Perturbasyon Analizi")
    st.stop()

    _selected_mode = st.session_state.get(SessionKeys.ACTIVE_MODE)
    if _selected_mode:
        st.markdown(f"## {_selected_mode}")
        render_mode_guide(_selected_mode)
        if _selected_mode == "Target Stress Search":
            st.info("Bir hedef gen seçin; ardından sol panelden **Target Stress Search** başlatın. CFTR, GLP1R veya PSEN1 ile deneyebilirsiniz.")
        elif _selected_mode == "Compensation Analysis":
            st.info("Bir hedef gen seçin; ardından sol panelden **Kompanzasyon Adaylarını Bul** seçeneğini çalıştırın.")
            st.caption("BETA: Sonuçlar deneysel kanıt değil, hesaplamalı ağ tahminidir.")
        else:
            st.info("Bir hedef gen veya hedef seti seçin, ardından simülasyonu başlatın. CFTR, GLP1R veya PSEN1 ile deneyebilirsiniz.")
        if st.button("Mod seçimine dön", key="mode_landing_return"):
            st.session_state[SessionKeys.ACTIVE_MODE] = None
            st.rerun()
        st.stop()

    _landing_mode = render_landing(
        proteins=g_canli.vcount() if g_canli is not None else 0,
        edges=g_canli.ecount() if g_canli is not None else 0,
        enrichment_sources=ENRICHMENT_SOURCE_COUNT,
        species_count=SUPPORTED_SPECIES_COUNT,
    )
    if _landing_mode:
        st.session_state[SessionKeys.ACTIVE_MODE] = _landing_mode
        st.rerun()
    st.stop()

    st.markdown(f"""
    <div class="sk-hero">
        <span class="sk-eyebrow"><span class="dot"></span>SİSTEM HAZIR · {g_canli.vcount() if g_canli is not None else 0:,} PROTEİN DÜĞÜMÜ YÜKLÜ</span>
        <div class="sk-hero-title">Bir geni sustur, ağın nasıl<br>tepki verdiğini <span class="accent">gör.</span></div>
        <div class="sk-hero-sub">
            Sophiark, seçtiğin hedefi StringDB tabanlı gerçek bir insan protein-protein
            etkileşim ağında "susturur" ve şokun ağ boyunca nasıl yayıldığını simüle eder.
            Sonuçta hangi genlerin stres altına girdiğini, hangi bölgelerin en çok
            etkilendiğini ve hangi hedefin ilaç geliştirme açısından en uygun
            olduğunu gösteren bir hedef-önceliklendirme raporu üretir.
        </div>
        <div class="sk-capability-grid">
            <div class="sk-capability">
                <div class="sk-capability-num">01 · YAPILANDIR</div>
                <div class="sk-capability-title">Hedefini seç</div>
                <div class="sk-capability-desc">Gen adıyla ara, bir hastalığa bağlı gen setini içe aktar ya da hücresel bir bölgeyi manuel işaretle.</div>
            </div>
            <div class="sk-capability">
                <div class="sk-capability-num">02 · ÇALIŞTIR</div>
                <div class="sk-capability-title">Şoku yay</div>
                <div class="sk-capability-desc">Ağ önemi + Betweenness tabanlı bir difüzyon modeli, susturmanın ağdaki gerçek yayılımını hesaplar.</div>
            </div>
            <div class="sk-capability">
                <div class="sk-capability-num">03 · İNCELE</div>
                <div class="sk-capability-title">Hedefi önceliklendir</div>
                <div class="sk-capability-desc">Ağ haritası, doz-yanıt eğrisi ve yolak zenginleştirmesiyle desteklenen bir öncelik raporu al.</div>
            </div>
        </div>
        <div class="sk-hero-cta"><span class="arrow">←</span> Sol kenar çubuğundan bir hedef seç, sonra <strong style="color:#f5f5f5">Simülasyonu Başlat</strong>'a bas.</div>
    </div>
    """, unsafe_allow_html=True)

    st.markdown("""
    <div class="sk-stepper">
        <div class="sk-step active"><span class="sk-step-num">1</span><span class="sk-step-label">Yapılandır</span></div>
        <div class="sk-step-connector"></div>
        <div class="sk-step todo"><span class="sk-step-num">2</span><span class="sk-step-label">Çalıştır</span></div>
        <div class="sk-step-connector"></div>
        <div class="sk-step todo"><span class="sk-step-num">3</span><span class="sk-step-label">İncele</span></div>
    </div>
    """, unsafe_allow_html=True)

    _aktif_tur = "Fare" if IS_MOUSE else "İnsan"
    _tur_gecmisi = [k for k in st.session_state.get("sim_gecmis", []) if k.get("tur", _aktif_tur) == _aktif_tur]
    if _tur_gecmisi:
        st.markdown('<div class="section-header"><div class="section-header-icon">📜</div><div><div class="section-header-text">Son Simülasyonlar</div><div class="section-header-sub">Kaldığın yerden devam etmek için sol menüden aynı hedefi tekrar seç</div></div></div>', unsafe_allow_html=True)
        hist_cols = st.columns(min(len(_tur_gecmisi[-3:]), 3) or 1)
        for i, kayit in enumerate(reversed(_tur_gecmisi[-3:])):
            with hist_cols[i % len(hist_cols)]:
                st.markdown(f"""
                <div class="sk-capability">
                    <div class="sk-capability-num">{kayit['tarih']}</div>
                    <div class="sk-capability-title">{kayit['hedef']}</div>
                    <div class="sk-capability-desc">{kayit['doku']} · {kayit['stresli_sayisi']} stresli gen · Sistemik ağ kayması %{kayit['kayip_pct']:.2f}</div>
                </div>
                """, unsafe_allow_html=True)
    st.stop()

if not app_state.simulation_result_context_matches(
    st.session_state.get("simulation_result_context"),
    species="Fare" if IS_MOUSE else "İnsan",
    tissue=hedef_doku,
    engine=calculation_engine,
    targets=hedef_genler_sniper,
    localization=lokalizasyon,
    block_strength=float(st.session_state.block_strength),
    damping=float(st.session_state.damping),
    candidate_limit=int(st.session_state.get("redistribution_candidate_limit", 100)),
    bc_sample_sources=st.session_state.bc_sample_sources,
    tissue_normalization_mode=st.session_state.get("tissue_normalization_mode", "within_tissue"),
):
    st.info(
        "Ekrandaki seçimler tamamlanmış sonuçla eşleşmiyor. Eski target, tissue, "
        "engine veya parametre sonucunu göstermemek için yeni bir analiz çalıştırın."
    )
    st.stop()

if not isinstance(df, pd.DataFrame):
    st.error("Veri formatı beklenmeyen tipte.")
    st.stop()

# One canonical source feeds every downstream surface. Compare intentionally
# keeps Classic as its presentation result and renders the engine comparison
# separately; it never blends candidate sets.
_directed_bundle = st.session_state.get("directed_run_bundle")
_active_flow = select_active_result_flow(
    classic_report=df,
    classic_signed_response=st.session_state.get("signed_redistribution_result"),
    requested_mode=st.session_state.get("directed_engine_mode", "Classic"),
    directed_bundle=_directed_bundle,
)
df = _active_flow.report
_active_engine = _active_flow.engine
_active_presentation = _active_flow.presentation
_active_candidate_fingerprint = (
    _active_presentation.candidate_fingerprint if _active_presentation is not None
    else candidate_set_fingerprint(
        _active_engine, df.loc[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL, "gene"].astype(str).tolist(),
    )
)
_active_enrichment_matches = (
    st.session_state.get("enrichment_engine") == _active_engine
    and st.session_state.get("enrichment_candidate_fingerprint") == _active_candidate_fingerprint
)
_active_enrichment = st.session_state.get("enrichment_results") if _active_enrichment_matches else None
_active_enrichment_notices = st.session_state.get("enrichment_notices", []) if _active_enrichment_matches else []
_pipeline_log("PRESENTATION", "%s sonucu downstream yüzeylere bağlandı (rows=%d)", _active_engine, len(df))
_completed_engine_mode = st.session_state.get("directed_engine_mode", "Classic")
if calculation_engine != _completed_engine_mode:
    st.warning(
        f"Ekrandaki sonuç {_completed_engine_mode} motoruyla üretildi. "
        f"{calculation_engine} seçimini uygulamak için Analizi Başlat düğmesine yeniden basın."
    )

_context_target = " + ".join(_gene_to_symbol.get(gene, gene) for gene in (hedef_genler_sniper or [])) or lokalizasyon
render_analysis_hero(
    target=_context_target,
    species="Mus musculus" if IS_MOUSE else "Homo sapiens",
    tissue=doku_label,
    engine=_completed_engine_mode,
)
_findings_tab, _network_tab, _biology_tab, _research_tab, _data_tab, _tools_tab = result_destinations()
with _data_tab:
    _index_candidates = int(df["Hasar_Tipi"].astype(str).eq(motor.HASAR_UCUNCUL).sum()) if "Hasar_Tipi" in df else 0
    _index_enrichment = len(_active_enrichment) if isinstance(_active_enrichment, pd.DataFrame) else 0
    _index_signed = len(_active_flow.signed_response) if isinstance(_active_flow.signed_response, pd.DataFrame) else 0
    st.markdown("### Sonuç dizini")
    _index_cols = st.columns(4)
    _index_cols[0].metric("Ağ sonuçları", f"{len(df):,}")
    _index_cols[1].metric("Aday sonuçları", f"{_index_candidates:,}")
    _index_cols[2].metric("Biyolojik sonuçlar", f"{_index_enrichment:,}")
    _index_cols[3].metric("İşaretli ağ yanıtı", f"{_index_signed:,}")
    _complete_summary_surface = st.container()
    with st.expander("Analiz bağlamı ve kullanılan parametreler", expanded=False):
        _parameter_rows = [
            ("Hedef", _context_target),
            ("Organizma", "Mus musculus" if IS_MOUSE else "Homo sapiens"),
            ("Doku", doku_label),
            ("Hesaplama motoru", _completed_engine_mode),
            ("Analiz modu", mod),
            ("Kalan kenar ağırlığı", st.session_state.get("block_strength")),
            ("PageRank damping", st.session_state.get("damping")),
            ("BC örnek kaynak sayısı", st.session_state.get("bc_sample_sources")),
            ("Yeniden-dağılım aday sınırı", st.session_state.get("redistribution_candidate_limit")),
            ("Doku normalizasyonu", st.session_state.get("tissue_normalization_mode")),
            ("Ağ görünümü düğüm tavanı", st.session_state.get("harita_max_node")),
            ("Ağ görünümü etiket tavanı", st.session_state.get("harita_max_etiket")),
        ]
        safe_dataframe(
            pd.DataFrame(
                [(name, "—" if value is None else str(value)) for name, value in _parameter_rows],
                columns=["Parametre", "Kaydedilen değer"],
            ),
            label="analysis_context_and_parameters", key="analysis_context_and_parameters",
            row_limit=100,
        )
        st.caption("Bu tablo sunum amaçlıdır; hesaplama parametrelerini değiştirmez.")
    if _directed_bundle is not None and st.session_state.get("directed_engine_mode") != "Classic":
        _directed_result = _directed_bundle.calculation
        st.markdown(
            f"### {ui_t('classic_vs_directed_comparison')}"
            if _completed_engine_mode == "Compare" else
            f"### {ui_t('directed_calculation_engine')} · BETA"
        )
        if _directed_result.status.value == "AVAILABLE":
            _coverage = float(_directed_result.metadata.get("direction_coverage", 0.0)) * 100.0
            _directed_tables = directed_export_tables(_directed_bundle, gene_to_symbol=_gene_to_symbol)
            _directed_column_config = {
                "Directed % değişim": st.column_config.NumberColumn(ui_t("directed_percent_change"), format="%+.4f%%"),
                "Direction": st.column_config.TextColumn(ui_t("directed_direction")),
                "Rank": st.column_config.NumberColumn("Rank", format="%d"),
                "Directed BC": st.column_config.NumberColumn("Directed BC", format="%.6f"),
                "Relation": st.column_config.TextColumn(ui_t("directed_relation")),
                "Support": st.column_config.TextColumn(ui_t("directed_support")),
                "Coverage": st.column_config.NumberColumn(ui_t("directed_coverage"), format="%.2f"),
                "Localization": st.column_config.TextColumn(ui_t("directed_localization")),
            }
            _comparison_column_config = {
                "Classic % değişim": st.column_config.NumberColumn(ui_t("classic_percent_change"), format="%+.4f%%"),
                "Directed % değişim": st.column_config.NumberColumn(ui_t("directed_percent_change"), format="%+.4f%%"),
                "Fark": st.column_config.NumberColumn(ui_t("difference"), format="%+.4f%%"),
                "Classic Rank": st.column_config.NumberColumn(ui_t("classic_rank"), format="%d"),
                "Directed Rank": st.column_config.NumberColumn(ui_t("directed_rank"), format="%d"),
                "Rank Shift": st.column_config.NumberColumn(ui_t("rank_shift"), format="%+d"),
                "Directed BC": st.column_config.NumberColumn("Directed BC", format="%.6f"),
                "Relation": st.column_config.TextColumn(ui_t("directed_relation")),
                "Support": st.column_config.TextColumn(ui_t("directed_support")),
                "Coverage": st.column_config.NumberColumn(ui_t("directed_coverage"), format="%.2f"),
                "Localization": st.column_config.TextColumn(ui_t("directed_localization")),
            }
            st.caption(
                f"{st.session_state.get('directed_engine_mode')} · yön çözümleme kapsamı %{_coverage:.2f} · "
                f"pertürbasyon: {_directed_result.metadata.get('perturbation_strategy', '—')} · "
                "işaret bilgisi hesap ağırlığına katılmamıştır."
            )
            if _directed_bundle.comparison is not None:
                st.info(
                    "Compare, bu analiz çalıştırmasındaki aynı hedef/doku/parametrelerle üretilen "
                    "Classic ve Directed yeniden-dağılımlarını karşılaştırır."
                )
                _comparison_display = compare_ui_frame(
                    _directed_bundle.comparison.table, _directed_bundle.presentation,
                ) if _directed_bundle.presentation is not None else directed_export_tables(
                    _directed_bundle, gene_to_symbol=_gene_to_symbol,
                )["Classic vs Directed"]
                safe_dataframe(
                    _directed_presentation_frame(_comparison_display).head(300),
                    label="Classic vs Directed Comparison", key="classic_directed_comparison", height=520,
                    column_config=_comparison_column_config,
                    export_df=_directed_tables["Classic vs Directed"],
                )
                _validation = dict(_directed_bundle.comparison.validation)
                st.caption(
                    f"Spearman ρ={_validation.get('spearman_rank_correlation', 0):.4f} · "
                    f"overlap@50={_validation.get('overlap_at_50', 0)} · "
                    f"overlap@100={_validation.get('overlap_at_100', 0)} · "
                    f"medyan |rank shift|={_validation.get('median_absolute_rank_shift', 0):.1f}"
                )
            elif _completed_engine_mode == "Compare" and _directed_bundle.comparison_error:
                st.warning(
                    "Classic vs Directed comparison tamamlanamadı: bilimsel karşılaştırma için "
                    "full signed Classic response gereklidir. Eksik kayıtlar sıfır kabul edilmedi. "
                    + _directed_bundle.comparison_error
                )
            else:
                _directed_display = (
                    directed_ui_frame(_directed_bundle.presentation)
                    if _directed_bundle.presentation is not None else _directed_result.report
                )
                safe_dataframe(
                    _directed_presentation_frame(_directed_display).head(300),
                    label="Directed Redistribution · BETA", key="directed_redistribution", height=520,
                    column_config=_directed_column_config,
                    export_df=_directed_tables["Directed Redistribution"],
                )
        else:
            st.warning(
                "Directed engine kullanılamadı; classic analiz değişmeden korunuyor. "
                + str(_directed_result.error or "Yönlü veri mevcut değil.")
            )
    if _completed_engine_mode == "Evidence":
        from src.evidence.engine import latest_result as _latest_evidence_result
        _evidence_result = _latest_evidence_result()
        if _evidence_result is not None:
            st.markdown("### Evidence Calculation Engine · BETA")
            st.caption(
                "STRING v12 non-text confidence · bounded text corroboration · physical context only · "
                "CORUM complex response is post-calculation and cannot change ranking."
            )
            _evidence_columns = [c for c in (
                "gene", "Delta_PageRank_Pct", "BC_Skoru", "Non_Text_Evidence",
                "Text_Dependency", "Applied_Text_Modifier", "Physical_Supported_Edges",
            ) if c in _evidence_result.signed_response]
            _evidence_column_config = {
                "Non_Text_Evidence": metric_number_column("Evidence_s_nontext"),
                "Text_Dependency": metric_number_column("Evidence_text_dependency"),
                "Applied_Text_Modifier": st.column_config.NumberColumn("Uygulanan metin değiştiricisi", format="%.3f"),
                "Physical_Supported_Edges": st.column_config.NumberColumn("Fiziksel destekli kenarlar", format="%d"),
            }
            safe_dataframe(
                _evidence_result.signed_response[_evidence_columns].head(300),
                label="Evidence Redistribution · BETA", key="evidence_redistribution", height=520,
                column_config=_evidence_column_config,
            )
            if not _evidence_result.complexes.empty:
                safe_dataframe(
                    _evidence_result.complexes.drop(columns=["Evidence_Run_Provenance_JSON"], errors="ignore"),
                    label="Target-associated Complex Response", key="evidence_complex_response", height=360,
                )
            if not _evidence_result.comparison.empty:
                safe_dataframe(
                    _evidence_result.comparison.drop(columns=["Evidence_Run_Provenance_JSON"], errors="ignore").head(300),
                    label="Classic vs Evidence", key="classic_evidence_comparison", height=520,
                )
            st.download_button(
                "Evidence sonucu CSV indir",
                data=_evidence_result.signed_response.to_csv(index=False).encode("utf-8-sig"),
                file_name="sophiark_evidence_calculation.csv", mime="text/csv",
                key="download_evidence_calculation",
            )
with _findings_tab:
    _insight_frame = (
        _active_presentation.candidates.copy()
        if _active_engine == "directed" and _active_presentation is not None
        else df.copy()
    )
    if "gene" in _insight_frame.columns:
        _insight_frame = uygula_symbol_map_lokal(_insight_frame, _gene_to_symbol)
        _insight_frame["MyGene Adı"] = _insight_frame["gene"].astype(str).map(
            lambda gene: mygene_display_fields(mygene_kutuphane_coz(gene, _mygene_kutuphane)).get("MyGene Adı", "—")
        )
    _forward_insight, _forward_notable = en_onemli_bulguyu_ozetle(_insight_frame, mode="Perturbasyon Analizi")

    if "Hasar_Tipi" not in df.columns or "gene" not in df.columns:
        st.error("CSV formatı hatalı. Sütunlar: " + str(list(df.columns)))
        with st.expander("Ham CSV"):
            safe_dataframe(df.head(20), label="hata_debug")
        st.stop()


    hedef_genler     = df[df["Hasar_Tipi"] == motor.HASAR_BIRINCIL]["gene"].tolist()
    stresli_genler   = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL]["gene"].tolist()
    dokulen_genler   = df[df["Hasar_Tipi"] == motor.HASAR_IKINCIL]["gene"].tolist()
    toplam_etkilenen = len(hedef_genler) + len(stresli_genler)


    def _rapor_sembolu(row: pd.Series) -> str:
        """Rapor metninde hiçbir zaman NaN ya da ham protein kimliği döndürme."""
        value = row.get("Symbol")
        if value is not None and not pd.isna(value) and str(value).strip().lower() not in ("", "nan", "none"):
            return str(value).strip()
        gene_id = str(row.get("gene", ""))
        return _harita_aliasi(gene_id, _gene_to_symbol) or "Alias bulunamadı"


    def _hedef_fonksiyonu(gene_id: str) -> str:
        record = mygene_kutuphane_coz(gene_id, _mygene_kutuphane)
        name = str(record.get("name", "") or "").strip()
        if name and name.lower() not in ("nan", "none"):
            return name
        go_terms = _go_terimlerini_cikar(record)
        return go_terms[0] if go_terms else "MyGene işlev kaydı bulunamadı"


    _hedef_ozetleri = []
    for _gene_id in hedef_genler[:6]:
        _target_rows = df[df["gene"].astype(str) == str(_gene_id)]
        _target_symbol = _rapor_sembolu(_target_rows.iloc[0]) if not _target_rows.empty else _harita_aliasi(str(_gene_id), _gene_to_symbol)
        _hedef_ozetleri.append(
            f'<div style="margin-top:0.22rem"><b>{_target_symbol or "Alias bulunamadı"}</b> '
            f'<span style="font-size:0.68rem;color:#85858c">{_gene_id}</span><br>'
            f'<span style="font-size:0.72rem;color:#b8b8bd">{_hedef_fonksiyonu(str(_gene_id))}</span></div>'
        )
    _hedef_ozet_html = "".join(_hedef_ozetleri) or '<span style="color:#85858c">Hedef bilgisi yok</span>'

    hinterland_skor_map: dict = {}
    if "Hinterland_Skoru" in df.columns:
        hinterland_skor_map = dict(zip(df["gene"], df["Hinterland_Skoru"].fillna(1.0)))


    def _safe_float(dataframe, col):
        try:
            value = float(dataframe[col].iloc[0]) if col in dataframe.columns and not dataframe.empty else float("nan")
            return value if np.isfinite(value) else None
        except (ValueError, TypeError):
            return None

    global_efficiency_change = _safe_float(df, "Global_Efficiency_Change_Pct")
    mean_local_efficiency_change = _safe_float(df, "Mean_Local_Efficiency_Change_Pct")
    target_local_efficiency_change = _safe_float(df, "Target_Local_Efficiency_Change_Pct")

    def _efficiency_severity(value):
        magnitude = abs(value) if value is not None else 0.0
        return "crit" if magnitude > 15 else "warn" if magnitude > 5 else "ok"


    def _efficiency_value(value):
        return "—" if value is None else f"{value:.2f}<span class=\"sk-readout-unit\">%</span>"


    _global_cls = _efficiency_severity(global_efficiency_change)
    _mean_local_cls = _efficiency_severity(mean_local_efficiency_change)
    _target_local_cls = _efficiency_severity(target_local_efficiency_change)


    if _active_engine == "directed" and _active_presentation is not None:
        _directed_meta = _active_presentation.metadata
        _directed_delta = pd.to_numeric(_active_presentation.full_response["Delta_PageRank_Pct"], errors="coerce")
        _strongest_gain = _directed_delta.max() if not _directed_delta.empty else np.nan
        _strongest_loss = _directed_delta.min() if not _directed_delta.empty else np.nan
        _complete_summary_cells = "".join((
            summary_cell_html("Pertürbasyon hedefi", str(len(hedef_genler))),
            summary_cell_html("Directed aday", str(len(stresli_genler))),
            summary_cell_html("Etkilenen kayıt", str(_directed_meta["affected_record_count"])),
            summary_cell_html("Pozitif yanıt", str(_directed_meta["positive_response_count"])),
            summary_cell_html("Negatif yanıt", str(_directed_meta["negative_response_count"])),
            summary_cell_html("En güçlü artış", f'{_strongest_gain:.2f}<span class="sk-readout-unit">%</span>'),
            summary_cell_html("En güçlü kayıp", f'{_strongest_loss:.2f}<span class="sk-readout-unit">%</span>'),
            summary_cell_html("Motor / anlamlılık", "DIRECTED BETA · N/A"),
        ))
        _summary_cells = "".join((
            summary_cell_html("Pertürbasyon hedefi", str(len(hedef_genler))),
            summary_cell_html("Directed aday", str(len(stresli_genler))),
            summary_cell_html("Etkilenen kayıt", str(_directed_meta["affected_record_count"])),
            summary_cell_html("En güçlü artış", f'{_strongest_gain:.2f}<span class="sk-readout-unit">%</span>'),
        ))
    else:
        _complete_summary_cells = "".join((
            summary_cell_html(ui_t("perturbation_target"), str(len(hedef_genler))),
            summary_cell_html(ui_t("redistribution_candidate"), str(len(stresli_genler))),
            summary_cell_html(ui_t("affected_records"), str(toplam_etkilenen)),
            summary_cell_html(ui_t("mean_local_efficiency_change"), _efficiency_value(mean_local_efficiency_change), _mean_local_cls),
            summary_cell_html(ui_t("global_efficiency_change"), _efficiency_value(global_efficiency_change), _global_cls),
            summary_cell_html(ui_t("target_local_efficiency_change"), _efficiency_value(target_local_efficiency_change), _target_local_cls),
        ))
        _summary_cells = "".join((
            summary_cell_html(ui_t("perturbation_target"), str(len(hedef_genler))),
            summary_cell_html(ui_t("redistribution_candidate"), str(len(stresli_genler))),
            summary_cell_html(ui_t("affected_records"), str(toplam_etkilenen)),
            summary_cell_html(ui_t("global_efficiency_change"), _efficiency_value(global_efficiency_change), _global_cls),
        ))

    _summary_mode_label = (
        "COMPARE · CLASSIC BASE"
        if st.session_state.get("directed_engine_mode") == "Compare"
        else _active_engine.upper()
    )
    _summary_html = analysis_summary_html(
        title=ui_t("analysis_summary"), target_html="",
        context_label="", cells=(_summary_cells,),
    )
    st.markdown(_summary_html, unsafe_allow_html=True)
    with _complete_summary_surface:
        st.markdown("### Ayrıntılı çıktı özeti")
        st.markdown(
            analysis_summary_html(
                title="Kararlı sonuç ölçümleri",
                target_html=_hedef_ozet_html,
                context_label=f"{doku_label} · {_summary_mode_label}",
                cells=(_complete_summary_cells,),
            ),
            unsafe_allow_html=True,
        )

    _overview_candidates = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL].copy()
    _overview_candidates = uygula_symbol_map_lokal(_overview_candidates, _gene_to_symbol)
    _overview_signed = _active_flow.signed_response.copy(deep=True)
    _overview_signed = uygula_symbol_map_lokal(_overview_signed, _gene_to_symbol)
    _overview_delta = pd.to_numeric(
        _overview_signed.get("Delta_PageRank_Pct", pd.Series(index=_overview_signed.index, dtype=float)),
        errors="coerce",
    )
    _overview_losses = _overview_signed.loc[_overview_delta.lt(0)].copy()
    if not _overview_losses.empty:
        _overview_losses = _overview_losses.assign(_display_delta=_overview_delta.loc[_overview_losses.index])
        _overview_losses = _overview_losses.sort_values("_display_delta", ascending=True, kind="mergesort")

    _target_symbols_for_relation = [_gene_to_symbol.get(gene, gene) for gene in (hedef_genler or [])]

    def _overview_rows(source: pd.DataFrame, *, tone: str, limit: int = 4) -> list[dict[str, str]]:
        rows: list[dict[str, str]] = []
        for _, _candidate in source.head(limit).iterrows():
            _candidate_delta = pd.to_numeric(pd.Series([_candidate.get("Delta_PageRank_Pct")]), errors="coerce").iloc[0]
            _candidate_importance = pd.to_numeric(pd.Series([_candidate.get("Hinterland_Skoru")]), errors="coerce").iloc[0]
            _candidate_location = str(_candidate.get("Lokalizasyon", "—")).replace("_", " ").replace("|", " · ")
            _candidate_symbol = str(_candidate.get("Symbol") or _candidate.get("gene") or "—")
            _recorded_relation = _candidate.get("Hedef Gen Etkisi")
            if pd.isna(_recorded_relation) if not isinstance(_recorded_relation, (list, dict)) else False:
                _recorded_relation = None
            if not _recorded_relation or str(_recorded_relation).strip() in {"—", "nan", "None"}:
                _recorded_relation = get_target_effect_on_stressed(_target_symbols_for_relation, _candidate_symbol)
            _relation = str(_recorded_relation or "—")
            if _relation == "—":
                _relation = "Yön kanıtı yok"
            rows.append({
                "symbol": str(_candidate.get("Symbol") or _candidate.get("gene") or "—"),
                "delta": "—" if pd.isna(_candidate_delta) else f"{_candidate_delta:+.2f}%",
                "importance": "—" if pd.isna(_candidate_importance) else f"{_candidate_importance:.1f}",
                "localization": _candidate_location,
                "tone": tone,
                "response": "↑ Ağ önemi arttı" if tone == "gain" else "↓ Ağ önemi azaldı",
                "relation": _relation,
            })
        return rows

    _overview_gain_rows = _overview_rows(_overview_candidates, tone="gain")
    _overview_loss_rows = _overview_rows(_overview_losses, tone="loss")

    st.markdown("### Yanıt dengesi")
    st.caption("Pozitif ve negatif ağ yanıtları birlikte gösterilir. Ağ yanıtının işareti, biyolojik aktivasyon/baskılama anlamına gelmez.")
    _overview_left, _overview_right = st.columns(2, gap="large")
    with _overview_left:
        st.markdown("#### Ağ önemi artanlar")
        if _overview_gain_rows:
            st.markdown(candidate_stack_html(_overview_gain_rows), unsafe_allow_html=True)
        else:
            st.caption("Bu sonuçta pozitif yeniden dağılım adayı yok.")
    with _overview_right:
        st.markdown("#### Ağ önemi azalanlar")
        if _overview_loss_rows:
            st.markdown(candidate_stack_html(_overview_loss_rows), unsafe_allow_html=True)
        else:
            st.caption("Bu sonuçta negatif yeniden dağılım adayı yok.")
    st.markdown('<div class="sk-inline-action">Tüm pozitif ve negatif kayıtlar <span>→</span> <b>2 · Ağı Keşfet › Ağ Yanıtı</b></div>', unsafe_allow_html=True)

    _overview_network, _overview_finding = st.columns([.85, 1.15], gap="large")
    with _overview_network:
        st.markdown("### Ağ yanıtı")
        _network_value = "—" if global_efficiency_change is None else f"{global_efficiency_change:.2f}%"
        st.markdown(highlight_card_html(
            title="Global verimlilik değişimi",
            body=f'<strong>{_network_value}</strong><p>{toplam_etkilenen} etkilenmiş kayıt</p>',
            meta="Hedef kenar zayıflatması sonrası hesaplanan ağ sonucu",
            tone="blue",
        ), unsafe_allow_html=True)
    with _overview_finding:
        st.markdown("### Öne çıkan bulgu")
        render_insight_box(_forward_insight, is_notable=_forward_notable)

    _overview_biology, _overview_context = st.columns(2, gap="large")
    with _overview_biology:
        st.markdown("### Biyolojik vurgu")
        if isinstance(_active_enrichment, pd.DataFrame) and not _active_enrichment.empty:
            _bio_row = _active_enrichment.iloc[0]
            _bio_term_field = next((field for field in ("Term", "Pathway") if field in _active_enrichment), None)
            _bio_term = str(_bio_row.get(_bio_term_field, "Kayıtlı zenginleştirme sonucu")) if _bio_term_field else "Kayıtlı zenginleştirme sonucu"
            st.markdown(highlight_card_html(title=_bio_term, body="Mevcut zenginleştirme sonuçlarında ilk kayıt.", meta="Mevcut biyolojik bağlam", tone="green"), unsafe_allow_html=True)
        else:
            st.markdown(highlight_card_html(title="Biyolojik bağlam", body="Bu sonuç için kayıtlı zenginleştirme sonucu bulunmuyor.", tone="neutral"), unsafe_allow_html=True)
    with _overview_context:
        st.markdown("### Önemli bağlam")
        st.markdown(highlight_card_html(
            title="Sonuç sınırı",
            body="Ağ yeniden dağılımı biyolojik nedensellik veya tedavi önerisi değildir.",
            meta=f"{doku_label} · {_completed_engine_mode}", tone="amber",
        ), unsafe_allow_html=True)

    if stresli_genler:
        _verdict_rank = "Delta_PageRank_Pct" if _active_engine == "directed" else "Hinterland_Skoru"
        top3_df = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL].nlargest(3, _verdict_rank)
        top3_list = []
        for _, row in top3_df.iterrows():
            sym = _rapor_sembolu(row)
            top3_list.append(sym)

        lok_df = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL].copy()
        lok_df = render_gate_column(lok_df)
        lok_list = []
        for loc in lok_df["Lokalizasyon"]:
            if isinstance(loc, str):
                lok_list.extend([l.strip().replace("_", " ") for l in loc.split("|")])
        en_cok_lok = max(set(lok_list), key=lok_list.count) if lok_list else "Bilinmiyor"

        gate_df = lok_df[lok_df["Gümrük_Kapisi"] == "✓"] if "Gümrük_Kapisi" in lok_df else pd.DataFrame()
        if not gate_df.empty:
            best_gate = gate_df.nlargest(1, _verdict_rank).iloc[0]
            best_gate_sym = _rapor_sembolu(best_gate)
        else:
            best_gate_sym = "Yok"
        if _active_engine == "directed":
            verdict = (
                f"Directed hesapta en yüksek pozitif yeniden-dağılım {' · '.join(f'<b>{item}</b>' for item in top3_list)} "
                f"kayıtlarında gözlendi. Kayıtlar en sık <b>{en_cok_lok}</b> hücresel bağlamında görülüyor. "
                "Classic verimlilik ve Compartment Bottleneck ölçümleri Directed sonuç gibi sunulmamıştır."
            )
        elif ui_locale() == "en":
            verdict = (
                f"The highest relative network redistribution was observed for {' · '.join(f'<b>{item}</b>' for item in top3_list)}. "
                f"These records most often occur in the <b>{en_cok_lok}</b> cellular context. "
                f"<b>{best_gate_sym}</b> was identified as a Compartment Bottleneck candidate. This is a network-model result, not a treatment or biological-causality recommendation."
            )
        else:
            verdict = (
                f"En yüksek relatif ağ yeniden-dağılımı {' · '.join(f'<b>{item}</b>' for item in top3_list)} kayıtlarında gözlendi. "
                f"Kayıtlar en sık <b>{en_cok_lok}</b> hücresel bağlamında görülüyor. "
                f"<b>{best_gate_sym}</b> Compartment Bottleneck adayı olarak işaretlendi. Bu özet ağ-modeli sonucudur; tedavi veya biyolojik nedensellik önerisi değildir."
            )
        st.markdown(f"""
        <div class="sk-verdict">
            <div class="sk-verdict-tag">{ui_t("network_redistribution")}</div>
            <div class="sk-verdict-text">
                {verdict}
            </div>
        </div>
        """, unsafe_allow_html=True)

    render_scientific_note(
        "Directed mode: values are directed PageRank redistribution after incident-edge attenuation; significance is unavailable and Classic efficiency metrics are not reused."
        if _active_engine == "directed" else
        "Ağ önemi ve verimlilik değişimleri hedef kenar zayıflatması sonrasında hesaplanan pertürbasyon sonuçlarıdır. Geçiş Merkeziliği, Ağdaki Önemi ve topluluk alanları başlangıç ağ-topolojisi bağlamını temsil eder."
    )
    with st.expander("Nasıl yorumlanır?", expanded=False):
        render_reading_guide()
        render_mode_guide("Perturbasyon Analizi")
        render_results_guide()

# Result containers establish the desktop reading order while the existing
# render blocks keep their scientific inputs and component implementations.
with _network_tab:
    _candidate_surface = st.container()
    _cards_surface = st.container()
    _system_response_surface = st.container()
with _biology_tab:
    _biological_context_surface = st.container()
    _evidence_surface = st.container()
    _functional_surface = st.container()
    _deep_interpretation_surface = st.container()
with _research_tab:
    _research_surface = st.container()
with _data_tab:
    _statistical_audit_surface = st.container()
    _data_audit_surface = st.container()
with _tools_tab:
    _advanced_surface = st.container()

with _system_response_surface:
    render_section_header(
        "Sistem yanıtı",
        "Gerçek ağ görünümü, ortak Network Response matrisi ve kaynak sonuç tabloları aynı pertürbasyonun katmanlarıdır.",
        label="Ağ yanıtı",
    )

with _biological_context_surface:
    render_section_header(
        "Biyolojik öne çıkanlar",
        "",
        label="Biyolojik bağlam",
    )
    if isinstance(_active_enrichment, pd.DataFrame) and not _active_enrichment.empty:
        _bio_highlight_fields = [field for field in (
            "Term", "Pathway", "Overlap", "Genes",
        ) if field in _active_enrichment]
        if _bio_highlight_fields:
            _bio_highlight_labels = {
                "Term": "Terim",
                "Pathway": "Yolak",
                "Overlap": "Örtüşme",
                "Genes": "Genler",
            }
            st.dataframe(
                _active_enrichment.loc[:, _bio_highlight_fields].head(5).rename(columns=_bio_highlight_labels),
                width="stretch", hide_index=True, height=220,
            )
    _bio_expression_panel, _bio_localization_panel = st.columns(2, gap="large")

with _deep_interpretation_surface:
    render_section_header(
        "Derinlemesine yorum",
        "Mevcut biyolojik yorumlama ve anotasyon ayrıntıları; ağ hesaplaması veya yeni skor üretmez.",
        label="Yorum",
    )

with _advanced_surface:
    render_section_header(
        "İleri analizler",
        "Ana sonuç akışının dışındaki mevcut deneysel ve karşılaştırmalı araçlar.",
        label="İleri araçlar",
    )
    if _active_engine == "classic" and _completed_engine_mode != "Evidence":
        from src.ui.advanced_workflows import render_advanced_workflows
        render_advanced_workflows(motor=motor, graph=g_canli, scores=df_canli, report=df,
            targets=hedef_genler, tissue=hedef_doku, is_mouse=IS_MOUSE, output_dir=OUTPUTS_DIR,
            mapping=_gene_to_symbol, tissues=DOKU_SECENEKLERI, ghost_filter=ghost_gen_filtrele,
            context=st.session_state["simulation_result_context"])
    else:
        st.info("Doz–yanıt, şiddet, yayılım ve doku araçları Classic modeline aittir. Bu araçlar için Classic analizi seçin.")


with _system_response_surface:
    @st.fragment
    def _render_local_network_map():
        st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 1.5rem 0;"></div>', unsafe_allow_html=True)


        st.markdown(f"""
        <div class="section-header">
            <div class="section-header-icon">🔬</div>
            <div>
                    <div class="section-header-text">{ui_t("network_visualization")}</div>
                <div class="section-header-sub">
                    Mavi elmas = pertürbasyon hedefi · Yeşil = yeniden-dağılım adayı · Gri = ağ bağlamı · Sarı = 2D köprü düğümü
                </div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        _g = st.session_state.get("G_canli")
        _harita_max_node = st.session_state.get("harita_max_node", 600)
        _harita_max_etiket = st.session_state.get("harita_max_etiket", 40)

        if _g is not None and isinstance(_g, ig.Graph) and _g.vcount() > 0 and (hedef_genler or stresli_genler):
            col_mod, col_bos = st.columns([2, 3])
            with col_mod:
                harita_modu = st.radio(
                    "Harita Modu",
                    ["🎯 Sadece İlgili Genler (Hızlı)", "🌐 Gerçek Ağ — 2B", "🧊 Gerçek Ağ — 3B"],
                    horizontal=True,
                    key="harita_modu_secim",
                    label_visibility="collapsed",
                    help=(
                        "Sadece İlgili: yalnızca kilitlenen+stresli genler ve aralarındaki "
                        "bağlantılar (küçük, akıcı, fizik simülasyonlu).\n\n"
                        "Gerçek Ağ — 2B: bu genlerin ağdaki gerçek komşuluğunu (hub'lar dahil) "
                        "gösterir; yalnızca en yüksek skorlu çekirdek genlere etiket basılır.\n\n"
                        "Gerçek Ağ — 3B: aynı alt-graf, döndürülebilir/yakınlaştırılabilir 3 "
                        "boyutlu görünümde. Çok sayıda çekirdek düğüm varsa (500+) 2B görünüm "
                        "kaçınılmaz olarak sıkışık kalabilir — 3B, düğümleri Z ekseninde de "
                        "ayırarak okunabilirliği artırır."
                    ),
                )
            tam_ag_modu = harita_modu.startswith("🌐")
            tam_ag_3d   = harita_modu.startswith("🧊")

            st.markdown(
                f'<div class="table-limit-note">✓ Bu, StringDB tabanlı GERÇEK ağdan çizilmiştir: '
                f'toplam ağ {_g.vcount():,} düğüm / {_g.ecount():,} kenar içeriyor. '
                f'"Sadece İlgili Genler" modunda, çekirdek genler arasında doğrudan bağlantı yoksa '
                f'ama ara bir gen (Y) üzerinden bağlıysa (X—Y—Z), o ara gen sarı "Köprü Düğüm" '
                f'olarak gösterilir — X ile Z arasında asla sahte/doğrudan bir çizgi çizilmez.</div>',
                unsafe_allow_html=True,
            )
            _tum_etiket_goster = st.checkbox(
                "🏷️ Tüm çekirdek genlere etiket göster (yavaşlatabilir)",
                value=False, key="tum_etiket_goster_chk",
                help="Varsayılan olarak yalnızca en yüksek skorlu genlere metin etiketi basılır, "
                     "geri kalanı üzerine gelince (hover) görünür. Bunu işaretlerseniz TÜM "
                     "kilitlenen/stresli genlere etiket basılır — çok sayıda gen varsa harita "
                     "kalabalıklaşabilir.",
            )

            core_sayisi = len(set(hedef_genler) | set(stresli_genler))
            if (tam_ag_modu or tam_ag_3d) and core_sayisi > 600:
                st.markdown(
                    f'<div class="stale-banner">⚠ Çekirdek (kilitlenen+stresli) düğüm sayısı '
                    f'{core_sayisi:,} — bu tavandan bağımsız olarak HER ZAMAN dahil edilir, '
                    f'bu yüzden çok yoğun olabilir. "3B" görünümü veya "⚙️ Gelişmiş Simülasyon '
                    f'Ayarları" içindeki etiket/düğüm tavanlarını düşürmek okunabilirliği artırır.</div>',
                    unsafe_allow_html=True,
                )
            elif (tam_ag_modu or tam_ag_3d) and _g.vcount() > _harita_max_node:
                st.markdown(
                    f'<div class="stale-banner">🌐 Tüm ağ {_g.vcount():,} düğüm içeriyor — '
                    f'tarayıcı performansı için hedef+stresli genler ve en yüksek skorlu '
                    f'{_harita_max_node:,} düğüme kadar bağlam gösteriliyor. Tavanı '
                    f'"⚙️ Gelişmiş Simülasyon Ayarları" içinden değiştirebilirsin.</div>',
                    unsafe_allow_html=True,
                )

            _etiket_limiti = 100000 if _tum_etiket_goster else _harita_max_etiket

            with st.spinner("Harita çiziliyor..."):
                try:
                    if tam_ag_3d:
                        fig3d = draw_incident_scene_3d(
                            _g, hedef_genler, stresli_genler, hinterland_skor_map, _gene_to_symbol,
                            max_nodes=_harita_max_node, max_etiket=_etiket_limiti,
                        )
                        if fig3d is not None:
                            st.plotly_chart(fig3d, width="stretch")
                        elif not PLOTLY_OK:
                            st.warning("3B görünüm için `plotly` kurulu değil. `pip install plotly` ile ekleyebilirsiniz.")
                        else:
                            st.info("Grafta eşleşen düğüm yok.")
                    else:
                        html_path = draw_incident_scene(
                            _g, hedef_genler, stresli_genler, hinterland_skor_map, _gene_to_symbol,
                            tam_ag=tam_ag_modu, max_nodes=_harita_max_node, max_etiket=_etiket_limiti,
                        )
                        if html_path:
                            with open(html_path, "r", encoding="utf-8") as f:
                                st.iframe(f.read(), height=560)
                            os.unlink(html_path)
                        elif not PYVIS_OK:
                            st.warning("2B ağ görünümü için `pyvis` kurulu değil. `requirements.txt` bağımlılıklarını kurun.")
                        else:
                            st.info("Grafta eşleşen düğüm yok.")
                except Exception as e:
                    st.warning(f"Harita hatası: {e}")
        elif _g is not None and _g.vcount() == 0:
            st.warning("Graf boş — seçilen doku için düğüm bulunamadı.")
        else:
            st.info("Harita için önce simülasyonu çalıştırın.")

        st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 1.5rem 0;"></div>', unsafe_allow_html=True)

    _render_local_network_map()

with _bio_expression_panel:
    if hedef_genler:
        st.markdown("""
        <div class="section-header">
            <div class="section-header-icon">🫁</div>
            <div>
                <div class="section-header-text">Başlangıç ifade bağlamı</div>
                <div class="section-header-sub">Seçili hedefin kayıtlı başlangıç ifadesi; pertürbasyon etkisini veya yan etkiyi ölçmez</div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        try:
            ilk_hedef = hedef_genler[0]
            doku_df = query_gene_expression(ilk_hedef)
            if not doku_df.empty:
                top10 = doku_df.head(10)
                _expression_series = top10.set_index("Doku")["İfade Seviyesi"].copy()
                _expression_series.name = "Expression level" if ui_locale() == "en" else "İfade Seviyesi"
                st.bar_chart(_expression_series, width="stretch")
            else:
                st.info("Bu gen için doku ifade verisi bulunamadı.")
        except Exception as e:
            st.warning(f"Doku ifade grafiği çizilemedi: {e}")
    if not hedef_genler:
        st.info("Başlangıç ifade bağlamı için pertürbasyon hedefi bulunamadı.")
    render_scientific_note("Başlangıç doku ifadesi bir anotasyon bağlamıdır; pertürbasyon sonucu veya hedef zayıflatması sonrası ifade değişimi değildir.")

with _bio_localization_panel:
    if stresli_genler:
        st.markdown("""
        <div class="section-header">
            <div class="section-header-icon">📍</div>
            <div>
                <div class="section-header-text">Yeniden dağılım adayları · hücresel bağlam</div>
                <div class="section-header-sub">Hangi bölgede ne kadar stres birikmiş?</div>
            </div>
        </div>
        """, unsafe_allow_html=True)

        lok_df = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL].copy()
        lok_df = render_gate_column(lok_df)
        lok_list = []
        for loc in lok_df["Lokalizasyon"]:
            if isinstance(loc, str):
                lok_list.extend([l.strip().replace("_", " ") for l in loc.split("|")])

        if lok_list:
            from collections import Counter
            lok_say = Counter(lok_list)

            sirali_lok = lok_say.most_common(8)
            kalan = sum(lok_say.values()) - sum(n for _, n in sirali_lok)
            if kalan:
                sirali_lok.append(("Diğer", kalan))
            etiketler, degerler = zip(*sirali_lok)
            renkler = ['#4ade80', '#60a5fa', '#fbbf24', '#f87171', '#a78bfa',
                       '#34d399', '#f472b6', '#818cf8', '#64748b']

            fig, ax = plt.subplots(figsize=(7.2, 4.8), facecolor='#0c0c0d')
            wedges, texts, autotexts = ax.pie(
                degerler,
                labels=None,
                autopct=lambda pct: f"%{pct:.0f}" if pct >= 5 else "",
                colors=renkler[:len(degerler)],
                startangle=90,
                pctdistance=0.74,
                wedgeprops={'width': 0.52, 'edgecolor': '#0c0c0d', 'linewidth': 2},
                textprops={'color': '#f5f5f5', 'fontsize': 9}
            )
            for t in autotexts:
                t.set_color('#0c0c0d')
                t.set_fontweight('bold')
                t.set_fontsize(8)

            ax.text(0, 0.08, "Stresli", ha="center", va="center", color="#b8b8bd", fontsize=10)
            ax.text(0, -0.10, f"{sum(degerler)}", ha="center", va="center", color="#f5f5f5", fontsize=17, fontweight="bold")
            ax.legend(wedges, [f"{name} ({count})" for name, count in zip(etiketler, degerler)],
                      loc="center left", bbox_to_anchor=(1.0, 0.5), frameon=False,
                      labelcolor="#c8c8ce", fontsize=8)
            ax.set_facecolor('#0c0c0d')
            fig.tight_layout()
            st.pyplot(fig)
        else:
            st.info("Lokalizasyon verisi yok.")
    if not stresli_genler:
        st.info("Hücresel dağılım bağlamı için yeniden-dağılım adayı bulunamadı.")
    render_scientific_note("Lokalizasyon kayıtlı hücresel konum anotasyonudur; pertürbasyonun o organelde gerçekleştiğinin kanıtı değildir.")

with _cards_surface:
    st.markdown(f"""
    <div class="section-header">
        <div class="section-header-icon">🎯</div>
        <div>
            <div class="section-header-text">{ui_t("perturbation_targets")}</div>
            <div class="section-header-sub">{("Targets with applied attenuation" if ui_locale() == "en" else "Zayıflatma uygulanmış hedefler")}</div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.markdown(f'<span class="sk-eyebrow" style="margin-bottom:0.8rem;"><span class="dot"></span>{ui_t("network_redistribution")} · {ui_t("perturbation_targets")}</span>', unsafe_allow_html=True)

hedef_df = df[df["Hasar_Tipi"] == motor.HASAR_BIRINCIL].copy()
if not hedef_df.empty:
    hedef_df = uygula_symbol_map_lokal(hedef_df, _gene_to_symbol)
    hedef_df = uygula_trrust_regulators(hedef_df)
    hedef_df = uygula_essentiality(hedef_df)
    hedef_df = render_gate_column(hedef_df)
    hedef_df = mygene_alanlarini_ekle(hedef_df, _mygene_kutuphane)

    with _cards_surface:
        for _idx, _row in hedef_df.head(3).iterrows():
            with st.container(border=False):
                render_gene_card(_row, key=f"primary_{_idx}")

    goster = [c for c in ["Symbol", "gene", "PageRank_Baseline", "PageRank_Perturbed", "Delta_PageRank_Pct",
                           "Hinterland_Skoru", "Lokalizasyon",
                           "Gümrük_Kapisi", "BC_Skoru", "Efficiency_Kayip_Pct",
                           "Local_Efficiency_Kayip_Pct", "Drug_Score",
                           "Düzenleyici_TFler", "Essentiality", "Protein_Adi",
                           "GO_CC_Terimleri", "GO_MF_Terimleri", "Yapisal_Kanit",
                           "MyGene Adı", "GO Biyolojik Süreç", "GO Moleküler İşlev",
                           "GO Hücresel Bileşen"] if c in hedef_df.columns]
    col_cfg = {
        field: metric_number_column(field)
        for field in (
            "PageRank_Baseline", "PageRank_Perturbed", "Delta_PageRank_Pct",
            "Hinterland_Skoru", "BC_Skoru", "Efficiency_Kayip_Pct",
            "Local_Efficiency_Kayip_Pct", "Drug_Score",
        )
        if field in hedef_df.columns
    }
    if "Gümrük_Kapisi" in hedef_df.columns:
        col_cfg["Gümrük_Kapisi"] = metric_text_column("Gümrük_Kapisi")
    if "Essentiality" in hedef_df.columns:
        col_cfg["Essentiality"] = metric_text_column("Essentiality")
    col_cfg = col_cfg or None
else:
    goster = []
    col_cfg = None
    with _cards_surface:
        st.info("Birincil hedef yok.")

with _cards_surface:
    st.markdown(f"""
    <div class="section-header">
        <div class="section-header-icon">⚡</div>
        <div>
            <div class="section-header-text">{ui_t("network_redistribution")}</div>
            <div class="section-header-sub">{("Relative network role after target attenuation" if ui_locale() == "en" else "Hedef zayıflatması sonrasında değişen göreli ağ rolü")}</div>
        </div>
    </div>
    """, unsafe_allow_html=True)
    st.markdown(f'<span class="sk-eyebrow" style="margin-bottom:0.8rem;"><span class="dot"></span>{"PERTURBATION RESULT" if ui_locale() == "en" else "PERTÜRBASYON SONUCU"} · {ui_t("network_importance_change").upper()}</span>', unsafe_allow_html=True)

stresli_df = df[df["Hasar_Tipi"] == motor.HASAR_UCUNCUL].copy()
_negative_loss_df = pd.DataFrame()
_negative_columns: list[str] = []
_negative_config: dict = {}
_loss_empty_message: str | None = None
_interpretation_result = None
_signed_redistribution = _active_flow.signed_response
# Eski/diskten yüklenmiş raporlarda iki-yönlü sonuç ayrıca session state'e
# yazılmamış olabilir. Signed sütunlar raporun kendisindeyse yalnızca sunum
# için bu mevcut veriyi kullan; yeni bir ağ hesabı veya sıralama yapılmaz.
if not isinstance(_signed_redistribution, pd.DataFrame) and "Delta_PageRank_Pct" in df.columns:
    _signed_redistribution = df
if not isinstance(_signed_redistribution, pd.DataFrame):
    _loss_empty_message = (
        "Bu oturumda iki-yönlü signed sonuç bellekte yok. Eski/diskten yüklenmiş "
        "raporlar bu ayrı sonuç objesini içermez; yeni bir analiz çalıştırın."
    )
elif _signed_redistribution.empty:
    _loss_empty_message = "Motor iki-yönlü signed sonuç objesini üretti, ancak sonuç boş."
elif "Delta_PageRank_Pct" not in _signed_redistribution.columns:
    _loss_empty_message = (
        "Bu analiz sonucu Ağ önemi değişimi (%) alanını içermiyor; kayıp görünümü "
        "yeniden hesaplanmadan gösterilemez."
    )
with _system_response_surface:
    with st.expander("Ağ önemi kaybı · görünüm ayarları", expanded=False):
        st.caption("This list shows genes whose relative structural importance decreases after perturbation. Values describe only the relative change in the network model; they do not indicate loss of expression, protein amount, activity, or biological function." if ui_locale() == "en" else "Bu liste, pertürbasyon sonrasında ağ içindeki göreli yapısal önemi azalan genleri gösterir. Değerler yalnızca ağ modelindeki göreli değişimi ifade eder; ekspresyon, protein miktarı, aktivite veya biyolojik işlev kaybı anlamına gelmez.")
        _negative_threshold = st.number_input(
            ui_t("minimum_network_importance_loss"), 0.0, 100.0,
            float(st.session_state.get("negative_min_abs_delta_pct", .05)), .01,
            key="negative_min_abs_delta_pct",
        )
        _negative_limit = st.number_input(
            "Gösterilecek negatif aday sayısı", 10, 2000,
            int(st.session_state.get("negative_redistribution_top_n", 100)), 10,
            key="negative_redistribution_top_n",
        )
        _negative_loss_df = select_redistribution_losses(
            _signed_redistribution, min_abs_delta_pct=float(_negative_threshold), limit=int(_negative_limit),
        )
        if _negative_loss_df.empty:
            st.info(
                _loss_empty_message
                or "Signed iki-yönlü sonuç mevcut, ancak seçili eşikte negatif aday yok."
            )
        else:
            _negative_loss_df = uygula_symbol_map_lokal(_negative_loss_df, _gene_to_symbol)
            _negative_loss_df = uygula_trrust_regulators(_negative_loss_df)
            _negative_loss_df = uygula_essentiality(_negative_loss_df)
            _negative_loss_df = render_gate_column(_negative_loss_df)
            _negative_loss_df = mygene_alanlarini_ekle(_negative_loss_df, _mygene_kutuphane)
            _negative_columns = [column for column in [
                "Symbol", "gene", "Delta_PageRank_Pct", "PageRank_Baseline", "PageRank_Perturbed", "Directed_Rank",
                "Direction_Relation", "Direction_Support_Type", "Sign_Context", "Direction_Source_Coverage",
                "Hinterland_Skoru", "Gümrük_Kapisi", "BC_Skoru",
                "Derece", "Topluluk_ID", "Lokalizasyon", "Düzenleyici_TFler", "Essentiality",
            ] if column in _negative_loss_df.columns]
            _negative_config = {
                "Symbol": st.column_config.TextColumn("Gen"),
                "gene": st.column_config.TextColumn("ENSP"),
            "Delta_PageRank_Pct": metric_number_column("Delta_PageRank_Pct"),
            "PageRank_Baseline": metric_number_column("PageRank_Baseline"),
            "Hinterland_Skoru": metric_number_column("Hinterland_Skoru"),
            "Gümrük_Kapisi": metric_text_column("Gümrük_Kapisi"),
            "BC_Skoru": metric_number_column("BC_Skoru"),
            "Lokalizasyon": st.column_config.TextColumn("Hücresel konum"),
            "Essentiality": metric_text_column("Essentiality"),
            }
            _negative_config = {key: value for key, value in _negative_config.items() if key in _negative_columns}
if not stresli_df.empty:
    stresli_df = uygula_symbol_map_lokal(stresli_df, _gene_to_symbol)
    stresli_df = uygula_trrust_regulators(stresli_df)
    stresli_df = uygula_essentiality(stresli_df)
    stresli_df = render_gate_column(stresli_df)
    stresli_df = mygene_alanlarini_ekle(stresli_df, _mygene_kutuphane)

    if hedef_genler:
        hedef_semboller = [_gene_to_symbol.get(g, g) for g in hedef_genler]
        stresli_df["PubMed_Makale"] = [
            build_pubmed_link(hedef_semboller, s) for s in stresli_df["Symbol"]
        ]
        stresli_df["Hedef Gen Etkisi"] = stresli_df["Symbol"].apply(
            lambda s: get_target_effect_on_stressed(hedef_semboller, s)
        )

    # Network redistribution is ranked by the measured perturbation response;
    # Hinterland remains visible as baseline topology context.
    _siralama_kolonu = (
        "Delta_PageRank_Pct" if "Delta_PageRank_Pct" in stresli_df.columns
        else "Hinterland_skoru" if "Hinterland_skoru" in stresli_df.columns
        else "Hinterland_Skoru"
    )
    stresli_df = stresli_df.sort_values(_siralama_kolonu, ascending=False, kind="mergesort")
    _kart_sayisi = min(3, len(stresli_df))
    _kart_df = stresli_df.sort_values(_siralama_kolonu, ascending=False).head(_kart_sayisi)

    with _cards_surface:
        with st.expander("What does this table show?" if ui_locale() == "en" else "Bu tablo ne anlatıyor? · Basit okuma rehberi", expanded=False):
            st.markdown("""
    Bu tablo, seçtiğiniz hedefin bağlantıları zayıflatıldıktan sonra ağdaki hangi genlerin **göreli ağ öneminin değiştiğini** gösterir. En üstte görünen kayıt, “en önemli gen” olmak zorunda değildir; bu deneyde **en yüksek göreli ağ yeniden-dağılımını** gösteren kayıttır.

    **Önce şuna bakın: Ağ önemi değişimi (%).** Bu ana sonuç sütunudur. Pozitif değer, hedef pertürbasyonundan sonra o genin ağ içindeki göreli payının arttığını; negatif değer ise azaldığını söyler. Bu bir ifade artışı, protein miktarı artışı veya biyolojik aktivasyon ölçümü değildir. Yalnızca bu ağ modelindeki göreli ağ öneminin değişimidir.

    **Ağdaki Önemi** ve **Geçiş Merkeziliği** farklı şeylerdir. İkisi de simülasyon başlamadan önceki ağın fotoğrafıdır:

    - **Ağdaki Önemi:** Genin başlangıç ağındaki yapısal önemini özetler; pertürbasyona en güçlü yanıtı vereceğini tek başına söylemez.
    - **Geçiş Merkeziliği:** Genin farklı ağ bölgeleri arasında ne kadar sık geçiş noktası olduğunu gösterir; ağ önemi değişimi sonucu değildir.
    - **Lokalizasyon:** İlgili proteinin kayıtlı hücresel konum bilgisidir. Sonucun o bölgede fiziksel olarak gerçekleştiğini kanıtlamaz.

    Tablodaki diğer alanlar ek bağlam sağlar:

    - **Compartment Bottleneck:** Ağ bölgeleri veya hücresel bağlamlar arasındaki geçişte köprü rolü gösterebilen bir adaydır; ilaç hedefi veya biyolojik zorunluluk kanıtı değildir.
    - **Düzenleyici TFler / Hedef Gen Etkisi:** Yüklü TRRUST yönlü kanıt eşleşmeleridir. Kanıt yoksa bu biyolojik ilişkinin olmadığı anlamına gelmez.
    - **Essentiality:** Mevcut external essentiality annotation’ıdır; bu forward simülasyonun ürettiği bir sonuç değildir.
    - **Protein adı, MyGene ve GO alanları:** Kimlik ve biyolojik açıklama içindir. Ağda görülen değişimin doğrudan moleküler mekanizmasını ispatlamaz.
    - **PubMed:** Seçilen hedef ile ilgili aday için arama bağlantısıdır; otomatik literatür doğrulaması değildir.

    **q-value neden ana tabloda yok?** q-value yalnızca null-model/FDR analizi gerçekten çalıştırıldıysa “İstatistiksel Doğrulama” bölümünde gösterilir. Bu koşul yoksa q-value görmemek normaldir.

    Pratik okuma sırası: önce ağ önemi değişimi (%), sonra Ağdaki Önemi ve Geçiş Merkeziliği ile bunun zaten yüksek başlangıç merkeziliğine sahip bir gen mi yoksa düşük başlangıç öneminden belirgin yanıt veren bir gen mi olduğuna bakın. Son olarak kanıt ve annotation alanlarını bağımsız doğrulama için kullanın.
            """)

        st.markdown(
            f'<span class="sk-eyebrow" style="margin-bottom:0.8rem;"><span class="dot"></span>'
            f'TOP {_kart_sayisi} {("REDISTRIBUTION CANDIDATES" if ui_locale() == "en" else "YENİDEN DAĞILIM ADAYI")} · {ui_t("network_importance_change").upper()} {("ORDER" if ui_locale() == "en" else "SIRASI")}</span>',
            unsafe_allow_html=True,
        )
        for _idx, _row in _kart_df.iterrows():
            with st.container(border=True):
                render_gene_card(_row, key=f"stressed_{_idx}")

    goster_s = [c for c in ["Symbol", "gene", "Delta_PageRank_Pct", "PageRank_Baseline", "PageRank_Perturbed",
                              "Directed_Rank", "Direction_Relation", "Direction_Support_Type", "Sign_Context", "Direction_Source_Coverage",
                              "Hinterland_Skoru", "BC_Skoru", "Lokalizasyon",
                              "Gümrük_Kapisi", "Drug_Score",
                              "Düzenleyici_TFler", "Essentiality", "Protein_Adi",
                              "GO_CC_Terimleri", "GO_MF_Terimleri", "Yapisal_Kanit",
                              "MyGene Adı", "GO Biyolojik Süreç", "GO Moleküler İşlev",
                              "GO Hücresel Bileşen",
                              "PubMed_Makale", "Hedef Gen Etkisi"] if c in stresli_df.columns]
    col_cfg_s = {}
    _turkce_basliklar = {
        "Symbol": "Gen", "gene": "ENSP",
        "Lokalizasyon": "Hücresel konum", "Drug_Score": "İlaçlanabilirlik skoru",
        "Düzenleyici_TFler": "Düzenleyici TF'ler", "Essentiality": "Yaşamsallık",
        "Protein_Adi": "Protein adı", "GO_CC_Terimleri": "GO hücresel bileşen terimleri",
        "GO_MF_Terimleri": "GO moleküler işlev terimleri", "Yapisal_Kanit": "Yapısal kanıt",
        "MyGene Adı": "Gen açıklaması", "GO Biyolojik Süreç": "GO biyolojik süreç",
        "GO Moleküler İşlev": "GO moleküler işlev", "GO Hücresel Bileşen": "GO hücresel bileşen",
    }
    for _key, _label in _turkce_basliklar.items():
        if _key in stresli_df.columns:
            col_cfg_s[_key] = st.column_config.TextColumn(_label)
    if "Gümrük_Kapisi" in stresli_df.columns:
        col_cfg_s["Gümrük_Kapisi"] = metric_text_column("Gümrük_Kapisi")
    if "Essentiality" in stresli_df.columns:
        col_cfg_s["Essentiality"] = metric_text_column("Essentiality")
    if "PubMed_Makale" in stresli_df.columns:
        col_cfg_s["PubMed_Makale"] = st.column_config.LinkColumn(
            "PubMed", help="Hedef gen ve bu stresli geni birlikte PubMed'de ara", display_text="🔗 Ara",
        )
    if "Hedef Gen Etkisi" in stresli_df.columns:
        col_cfg_s["Hedef Gen Etkisi"] = st.column_config.TextColumn(
            "Hedef Gen Etkisi",
            help="Ana hedef genin bu stresli gen üzerindeki düzenleyici etkisi (TRRUST verisine göre)",
        )
    if "Delta_PageRank_Pct" in stresli_df.columns:
        col_cfg_s["Delta_PageRank_Pct"] = metric_number_column("Delta_PageRank_Pct")
    if "Hinterland_Skoru" in stresli_df.columns:
        col_cfg_s["Hinterland_Skoru"] = metric_number_column("Hinterland_Skoru")
    if "BC_Skoru" in stresli_df.columns:
        col_cfg_s["BC_Skoru"] = metric_number_column("BC_Skoru")
    if "Drug_Score" in stresli_df.columns:
        col_cfg_s["Drug_Score"] = metric_number_column("Drug_Score")
    col_cfg_s = col_cfg_s or None

    with _system_response_surface:
        render_network_response(
            targets=hedef_df,
            losses=_negative_loss_df,
            redistribution=stresli_df,
            target_columns=goster,
            loss_columns=_negative_columns,
            redistribution_columns=goster_s,
            target_column_config=col_cfg,
            loss_column_config=_negative_config,
            redistribution_column_config=col_cfg_s,
            loss_empty_message=_loss_empty_message,
        )

    _selection_modes = set(df.get("Redistribution_Selection_Mode", pd.Series(dtype=str)).dropna().astype(str))
    _has_null_fdr = _selection_modes == {"null_fdr"}
    _validation_columns = [column for column in ["Symbol", "gene", "q_value", "empirical_p", "significant_redistribution", "Selection_Reason"] if column in stresli_df.columns]
    _validation_frame = None
    if _has_null_fdr:
        _validation_frame = stresli_df[_validation_columns].copy()
        if "significant_redistribution" in _validation_frame.columns:
            _validation_frame = _validation_frame[
                _validation_frame["significant_redistribution"].fillna(False).astype(bool)
            ].copy()
    with _statistical_audit_surface:
        render_statistical_support(_validation_frame, columns=_validation_columns)

    try:
        if IS_MOUSE:
            raise RuntimeError("MOUSE_BIOLOGICAL_INTERPRETATION_DISABLED")
        _pipeline_log("INTERPRETATION", "%s adaylarından biyolojik yorum hazırlanıyor (n=%d)", _active_engine, len(stresli_df))
        _directed_evidence, _directed_metadata = _cached_directed_evidence(IS_MOUSE)
        _interpretation = interpret_forward_results(
            report=df,
            candidates=stresli_df,
            targets=list(hedef_genler or []),
            gene_to_symbol=_gene_to_symbol,
            tissue=doku_label,
            organism="Mus musculus" if IS_MOUSE else "Homo sapiens",
            edge_attenuation=float(st.session_state.get("block_strength", 0.001)),
            network_nodes=g_canli.vcount() if g_canli is not None else None,
            network_edges=g_canli.ecount() if g_canli is not None else None,
            directed_evidence=_directed_evidence,
            enrichment=_active_enrichment,
            loss_candidates=_negative_loss_df,
        )
        _interpretation_result = _interpretation
        _pipeline_log("INTERPRETATION", "Biyolojik yorum hazır (engine=%s)", _active_engine)
    except Exception as _interpretation_error:
        # Interpretation is intentionally non-blocking: no report or scientific output is invalidated.
        if not IS_MOUSE:
            log.exception("Biyolojik yorumlama katmanı render edilemedi: %s", _interpretation_error)
            with _deep_interpretation_surface:
                st.warning("Biyolojik yorumlama bu rapor için görüntülenemedi; hesaplanmış ağ sonuçları korunuyor.")
else:
    with _system_response_surface:
        render_network_response(
            targets=hedef_df,
            losses=_negative_loss_df,
            redistribution=stresli_df,
            target_columns=goster,
            loss_columns=_negative_columns,
            target_column_config=col_cfg,
            loss_column_config=_negative_config,
            loss_empty_message=_loss_empty_message,
        )
    with _statistical_audit_surface:
        render_statistical_support(None, columns=[])
    with _cards_surface:
        st.info("Stresli gen yok.")

with _advanced_surface:
    with st.expander("İleri teknik analizler · Doz-yanıt ve eşik taraması", expanded=False):
        st.caption("Doz–yanıt araçları yukarıdaki sekmelerde yer alır. Aşağıdaki eşik grafiği ayrı STRING eşik taramasıdır; hedef pertürbasyonuyla aynı deney değildir.")
        sweep_csv = OUTPUTS_DIR / "reports" / "threshold_sweep_raporu.csv"
        sweep_df = st.session_state.get("sweep_df")
        if sweep_df is None and sweep_csv.exists():
            sweep_df = pd.read_csv(sweep_csv)

        if isinstance(sweep_df, pd.DataFrame) and not sweep_df.empty:
            st.markdown("### 📈 Ağ Direnci ve LCC Oranı")
            fig, ax1 = plt.subplots(figsize=(8, 4), facecolor="#141929")
            ax1.set_facecolor("#111522")
            resilience = pd.to_numeric(sweep_df["Resilience"], errors="coerce")
            lcc_ratio = pd.to_numeric(sweep_df["LCC_Orani"], errors="coerce")
            thresholds = sweep_df["Esik"]
            ax1.plot(thresholds, resilience, color="#34BC6E", marker="o", linewidth=2, label="Ağ direnci")
            ax1.set_xlabel("Threshold value" if ui_locale() == "en" else "Eşik değeri", color="#BCC6D9")
            ax1.set_ylabel("Network resilience" if ui_locale() == "en" else "Ağ direnci", color="#34BC6E")
            ax1.tick_params(axis="both", colors="#BCC6D9")
            ax2 = ax1.twinx()
            ax2.plot(thresholds, lcc_ratio, color="#589BFF", marker="s", linewidth=2, label="LCC oranı")
            ax2.set_ylabel("LCC ratio" if ui_locale() == "en" else "LCC oranı", color="#589BFF")
            ax2.tick_params(axis="y", colors="#BCC6D9")
            for axis in (ax1, ax2):
                for spine in axis.spines.values():
                    spine.set_color("#2B3655")
            ax1.grid(axis="y", color="#2B3655", alpha=0.65, linewidth=0.8)
            for threshold, value in zip(thresholds, resilience):
                if pd.notna(value):
                    ax1.annotate(f"{value:.3f}", (threshold, value), xytext=(0, 8), textcoords="offset points", ha="center", color="#62D990", fontsize=8)
            for threshold, value in zip(thresholds, lcc_ratio):
                if pd.notna(value):
                    ax2.annotate(f"{value:.3f}", (threshold, value), xytext=(0, -14), textcoords="offset points", ha="center", color="#78B1FF", fontsize=8)
            lines = [ax1.get_lines()[0], ax2.get_lines()[0]]
            ax1.legend(lines, [line.get_label() for line in lines], loc="upper right", facecolor="#141929", edgecolor="#2B3655", labelcolor="#E6ECF7")
            fig.tight_layout()
            st.pyplot(fig, clear_figure=True)
        else:
            st.caption("Görüntülenecek eşik taraması verisi yok.")


with _functional_surface:
    if not stresli_genler:
        st.info("Önce bir simülasyon çalıştırın; KEGG ve GO sonuçları simülasyonla birlikte otomatik hazırlanır.")
    else:
        st.caption("KEGG ve GO Biological Process sonuçları, seçili doku ağındaki gen evreni kullanılarak otomatik hesaplanır.")
        if st.button("Yolak analizini yeniden dene", key="enrich_retry_btn", help="Yalnızca uzak kaynak erişimi geçici olarak kesildiyse kullanın."):
            _retry_symbols = [
                symbol for symbol in stresli_df.get("gene_symbol", stresli_df.get("Symbol", pd.Series(dtype=str))).astype(str)
                if symbol and symbol != UNRESOLVED_SYMBOL
            ]
            _retry_ids = list(g_canli.vs["name"]) if g_canli is not None and "name" in g_canli.vs.attributes() else []
            _retry_background = [_gene_to_symbol.get(gene, "") for gene in _retry_ids]
            with st.spinner("KEGG ve GO sonuçları yeniden hazırlanıyor…"):
                try:
                    _retry_result, _retry_notices = run_post_simulation_enrichment(
                        _retry_symbols, is_mouse=IS_MOUSE, background_symbols=_retry_background,
                    )
                    st.session_state["enrichment_results"] = _retry_result if not _retry_result.empty else None
                    st.session_state["enrichment_notices"] = _retry_notices
                    st.session_state["enrichment_engine"] = _active_engine
                    st.session_state["enrichment_candidate_fingerprint"] = _active_candidate_fingerprint
                    _active_enrichment = st.session_state["enrichment_results"]
                    _active_enrichment_notices = _retry_notices
                except Exception as error:
                    log.warning("Yolak analizi yeniden denemesi başarısız: %s", error)
                    st.session_state["enrichment_notices"] = ["Yolak sonuçları şu anda hazırlanamadı; simülasyon sonucu korunuyor."]

    _enrichment_notices = _active_enrichment_notices
    _auto_enrichment_display = _active_enrichment
    render_functional_context(
        _auto_enrichment_display if isinstance(_auto_enrichment_display, pd.DataFrame) else None,
        notices=list(_enrichment_notices or []),
        analysis_ran=bool(stresli_genler),
    )

if _interpretation_result is not None:
    with _deep_interpretation_surface:
        with st.expander("Ek biyolojik yorum ve anotasyon ayrıntıları", expanded=False):
            render_biological_interpretation(_interpretation_result)

with _research_surface:
    if not IS_MOUSE:
        render_section_header(
            "Araştırma bağlamı",
            "Mevcut analiz sonucuna bağlı yerel, deterministik ve salt okunur kanıt.",
        )
    _research_snapshot = None
    try:
        if IS_MOUSE:
            raise RuntimeError("MOUSE_RESEARCH_EXPLORER_DISABLED")
        _pipeline_log("RESEARCH", "%s sonucu için snapshot hazırlanıyor", _active_engine)
        _research_selection_mode = _research_selection_mode_from_report(df)
        _research_run_scope = app_state.research_simulation_scope()
        _research_scope_matches_report = bool(
            _research_run_scope
            and _research_run_scope.get("selection_mode") == _research_selection_mode
            and _active_engine == "classic"
        )
        _research_top_n = (
            _research_run_scope.get("top_n")
            if _research_scope_matches_report
            else _research_top_n_from_report(df)
        )
        if _research_selection_mode.casefold() != "top_n":
            _research_top_n = None
        _research_signed_source = _active_flow.signed_response
        _research_tested_count = (
            _research_run_scope.get("tested_count")
            if _research_scope_matches_report
            else _research_unique_gene_count(_research_signed_source)
            if _research_selection_mode.casefold() == "top_n"
            else None
        )
        _research_threshold_parameters = (
            dict(_research_run_scope.get("threshold_parameters") or {})
            if _research_scope_matches_report
            else _research_scientific_threshold_parameters(_research_selection_mode)
        )
        _research_threshold_parameters.update({
            "engine": _active_engine,
            "loss_presentation_metric": "Delta_PageRank_Pct",
            "loss_presentation_operator": "< -abs(threshold)",
            "loss_presentation_min_abs_delta_pct": float(_negative_threshold),
            "loss_presentation_limit": int(_negative_limit),
            "loss_presentation_only": True,
        })
        _research_enrichment = _active_enrichment
        if not isinstance(_research_enrichment, pd.DataFrame):
            _research_enrichment = None
        _research_provenance = []
        _research_directed_bundle = st.session_state.get("directed_run_bundle")
        if _active_engine == "directed" and _research_directed_bundle is not None:
            _research_directed = _research_directed_bundle.calculation
            _research_provenance.append(ProvenanceRecord(
                source="Sophiark Directed Calculation Engine",
                kind=ProvenanceKind.SOPHIARK_COMPUTED,
                version="v1",
                snapshot_id=str(_research_directed.metadata.get("direction_dataset_version", "unknown")),
                details={
                    "status": _research_directed.status.value,
                    "engine": _research_directed.metadata.get("engine", "directed"),
                    "graph_mode": _research_directed.metadata.get("graph_mode"),
                    "directed_perturbation_strategy": _research_directed.metadata.get("directed_perturbation_strategy"),
                    "suppression_factor": _research_directed.metadata.get("suppression_factor"),
                    "direction_coverage": _research_directed.metadata.get("direction_coverage"),
                    "fallback_edge_count": _research_directed.metadata.get("fallback_edge_count"),
                    "comparison_validation": (
                        dict(_research_directed_bundle.comparison.validation)
                        if _research_directed_bundle.comparison is not None else {}
                    ),
                    "scientific_boundary": "direction affects topology; activation/inhibition is context only",
                },
            ))
        _research_signed_input = (
            _active_flow.signed_response
            if _active_engine == "directed"
            else _negative_loss_df
        )
        _research_snapshot = build_research_snapshot(
            df,
            _research_signed_input,
            _research_enrichment,
            species=(
                _research_run_scope["species"]
                if _research_scope_matches_report
                else "Mus musculus" if IS_MOUSE else "Homo sapiens"
            ),
            taxon_id=(
                _research_run_scope["taxon_id"]
                if _research_scope_matches_report
                else 10090 if IS_MOUSE else 9606
            ),
            tissue=(
                _research_run_scope["tissue"]
                if _research_scope_matches_report
                else doku_label
            ),
            targets=(
                _research_run_scope["targets"]
                if _research_scope_matches_report
                else tuple(map(str, hedef_genler or ()))
            ),
            attenuation=(
                _research_run_scope["attenuation"]
                if _research_scope_matches_report
                else float(st.session_state.get("block_strength", 0.001))
            ),
            selection_mode=_research_selection_mode,
            test_limit=None,
            tested_count=_research_tested_count,
            returned_count=int(_research_unique_gene_count(stresli_df) or 0),
            top_n=_research_top_n,
            threshold_parameters=_research_threshold_parameters,
            provenance=_research_provenance,
            simulation_id=(
                st.session_state.get("research_current_simulation_id")
                if _research_scope_matches_report
                else None
            ),
        )
        _pipeline_log("RESEARCH", "Research snapshot hazır (engine=%s)", _active_engine)
    except Exception as error:
        if not IS_MOUSE:
            log.warning("Research snapshot could not be prepared: %s", error)
            st.caption("Bu eski sonuç biçimi için araştırma bağlamı kullanılamıyor; ham sonuçlar erişilebilir durumda.")

    if _research_snapshot is not None:
        _research_stale = app_state.bind_research_snapshot(
            snapshot_id=_research_snapshot.snapshot_id,
            simulation_id=_research_snapshot.simulation_id,
        )
        st.session_state["research_snapshot"] = _research_snapshot
        if _research_stale:
            st.caption("Aktif analiz sonucu değiştiği için araştırma bağlamı yenilendi.")
        _render_research_surface(_research_snapshot)

with _candidate_surface:
    _candidate_frames = [part for part in (hedef_df, stresli_df, _negative_loss_df) if isinstance(part, pd.DataFrame) and not part.empty]
    if _candidate_frames:
        _candidate_evidence = pd.concat(_candidate_frames, ignore_index=True, sort=False).drop_duplicates("gene", keep="first")
        _candidate_evidence = uygula_symbol_map_lokal(_candidate_evidence, _gene_to_symbol)
        _candidate_evidence = align_candidate_evidence(_candidate_evidence, _active_flow.signed_response)
        render_candidate_lens(_candidate_evidence, _active_enrichment)

df_goster = uygula_symbol_map_lokal(df, _gene_to_symbol)
with _data_audit_surface:
    render_raw_data(df_goster)

with _advanced_surface:
    with st.expander("Human–Mouse Comparison", expanded=False):
        _human_cross = st.session_state.get("cross_species_human_snapshot")
        _mouse_cross = st.session_state.get("cross_species_mouse_snapshot")
        if not _human_cross or not _mouse_cross:
            st.caption(
                "Karşılaştırma için Human ve Mouse Classic analizlerini ayrı ayrı çalıştırın. "
                "Tek tür sonucu diğer tür başarısız olduğunda korunur."
            )
        else:
            try:
                _orthology_index = _cached_human_mouse_orthology(str(BASE_DIR))
                _cross_result = compare_ortholog_responses(
                    _human_cross["response"], _mouse_cross["response"], _orthology_index,
                    human_config=_human_cross["config"], mouse_config=_mouse_cross["config"],
                )
                _cross_metrics = _cross_result["metrics"]
                _cross_cols = st.columns(4)
                _cross_cols[0].metric("Comparable 1:1", _cross_metrics["comparable_one_to_one_count"])
                _cross_cols[1].metric("Shared positive", _cross_metrics["concordant_positive"])
                _cross_cols[2].metric("Shared negative", _cross_metrics["concordant_negative"])
                _cross_cols[3].metric("Discordant", _cross_metrics["discordant"])
                st.caption(
                    "Primary comparison uses Ensembl Compara ONE_TO_ONE orthologs and normalized "
                    "rank percentiles. Raw ΔPageRank values are descriptive and are not interpreted "
                    "as directly comparable biological effect sizes."
                )
                st.dataframe(_cross_result["matched_responses"], width="stretch", hide_index=True)
                st.download_button(
                    "Ortholog-aware comparison CSV",
                    _cross_result["matched_responses"].to_csv(index=False).encode("utf-8-sig"),
                    file_name="sophiark_human_mouse_ortholog_comparison.csv",
                    mime="text/csv",
                )
            except ValueError as _cross_error:
                st.warning(str(_cross_error))
            except Exception as _cross_error:
                log.warning("Human–Mouse comparison unavailable: %s", _cross_error)
                st.warning("Human–Mouse comparison is safely unavailable; both source results are preserved.")

st.markdown('<div style="border-top:1px solid #2a2a2e; margin: 1.5rem 0 0.5rem 0;"></div>', unsafe_allow_html=True)
st.markdown("""
<div style="text-align:center; color:#5c5c62; font-family:'JetBrains Mono',monospace;
            font-size:0.62rem; letter-spacing:0.15em; padding-bottom:1.5rem;">
    SOPHIARK · NETWORK BIOLOGY WORKSPACE · v5.3
</div>
""", unsafe_allow_html=True)
