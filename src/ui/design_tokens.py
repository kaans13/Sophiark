"""Sophiark dark scientific workstation design system."""

from __future__ import annotations

import streamlit as st
from functools import lru_cache
from pathlib import Path

BG_PAGE = "#0B0D11"
BG_SIDEBAR = "#0E1116"
BG_SURFACE = "#12161C"
BG_SURFACE_2 = "#0F1318"
BG_SURFACE_HOVER = "#181D25"
BG_INPUT = "#151A21"
BG_ELEVATED = "#1A2029"
BORDER = "#242B35"
BORDER_STRONG = "#343E4B"
BORDER_FOCUS = "#607FA8"
TEXT_PRIMARY = "#F0F2F5"
TEXT_SECONDARY = "#C1C7D0"
TEXT_MUTED = "#89929F"
ACCENT = "#78A9E6"
ACCENT_BRIGHT = "#9BC1F0"
ACCENT_DARK = "#527DB5"
ACCENT_HOVER = "#8AB5EA"
ACCENT_SOFT = "rgba(120,169,230,.10)"
TIER_KRITIK = "#D97888"
TIER_ORTA = "#D5A55E"
TIER_DUSUK = "#72B99A"
BETA_BADGE = "#8996AD"
KANIT_GUCLU = ACCENT
KANIT_ZAYIF = TEXT_MUTED
INSIGHT_HIGHLIGHT = "rgba(120,169,230,.055)"
INSIGHT_BORDER = ACCENT
FONT_XS, FONT_SM, FONT_BASE, FONT_MD, FONT_LG, FONT_XL = "0.75rem", "0.8125rem", "0.875rem", "1rem", "1.25rem", "1.875rem"


@lru_cache(maxsize=2)
def _stylesheet(name: str) -> str:
    return Path(__file__).with_name(name).read_text(encoding="utf-8")


def apply_theme() -> None:
    """Apply one restrained workstation language to all Streamlit surfaces."""
    st.markdown("<style>" + _stylesheet("legacy_components.css") + "</style>", unsafe_allow_html=True)
    st.markdown(f"""
    <style>
    :root {{
      --bg-main:{BG_PAGE}; --bg-sidebar:{BG_SIDEBAR}; --bg-surface:{BG_SURFACE}; --bg-surface-2:{BG_SURFACE_2}; --bg-surface-hover:{BG_SURFACE_HOVER}; --bg-input:{BG_INPUT}; --bg-elevated:{BG_ELEVATED};
      --border-subtle:{BORDER}; --border-default:{BORDER_STRONG}; --border-focus:{BORDER_FOCUS}; --text-title:{TEXT_PRIMARY}; --text-main:{TEXT_SECONDARY}; --text-secondary:{TEXT_MUTED};
      --accent-primary:{ACCENT}; --accent-primary-hover:{ACCENT_HOVER}; --accent-primary-soft:{ACCENT_SOFT}; --accent-positive:{TIER_DUSUK}; --accent-negative:{TIER_KRITIK}; --accent-warning:{TIER_ORTA};
      --graph-target:{ACCENT}; --graph-positive:{TIER_DUSUK}; --graph-negative:{TIER_KRITIK}; --graph-neutral:#7d8794; --graph-edge:#3a4654; --accent-violet:#9a91c7;
      --bg-base:{BG_PAGE}; --bg-panel:{BG_SURFACE}; --bg-panel-alt:{BG_SURFACE_2}; --border-soft:{BORDER}; --border-strong:{BORDER_STRONG}; --text-primary:{TEXT_PRIMARY}; --text-tertiary:{TEXT_MUTED}; --accent:{ACCENT}; --accent-white:{ACCENT}; --accent-green:{TIER_DUSUK}; --accent-red:{TIER_KRITIK}; --accent-amber:{TIER_ORTA};
    }}
    /* Typography is inherited by UI text. Do not target every descendant: Streamlit
       renders Material Symbols as ligature text inside span nodes. */
    .stApp {{ font-family:Inter,"Segoe UI",system-ui,sans-serif!important; }}
    [data-testid="stIconMaterial"],.material-symbols-rounded,.material-symbols-outlined {{
      font-family:"Material Symbols Rounded","Material Symbols Outlined"!important;
      font-weight:normal!important;font-style:normal!important;letter-spacing:normal!important;
      text-transform:none!important;white-space:nowrap!important;word-wrap:normal!important;
      direction:ltr!important;-webkit-font-feature-settings:"liga"!important;font-feature-settings:"liga"!important;
    }}
    .stApp {{ background:{BG_PAGE}!important;color:{TEXT_SECONDARY}!important; }}
    .stApp h1,.stApp h2,.stApp h3,.stApp h4 {{ color:{TEXT_PRIMARY}!important;font-weight:600!important;letter-spacing:-.015em; }}
    .stApp h1 {{ font-size:30px!important; }} .stApp h2 {{ font-size:20px!important; }} .stApp p,.stApp label,.stApp span {{ color:{TEXT_SECONDARY}; }}
    [data-testid="stHeader"] {{ background:rgba(9,11,17,.94)!important; }} [data-testid="stDecoration"] {{ background:{ACCENT}!important;height:2px!important; }}
    [data-testid="stSidebar"],section[data-testid="stSidebar"],[data-testid="stSidebarCollapsedControl"] {{ display:none!important; }}
    .main .block-container {{ max-width:1680px!important;padding-left:clamp(1.25rem,3vw,3.5rem)!important;padding-right:clamp(1.25rem,3vw,3.5rem)!important; }}
    [data-testid="stSidebar"] p,[data-testid="stSidebar"] label,[data-testid="stSidebar"] h1,[data-testid="stSidebar"] h2,[data-testid="stSidebar"] h3,[data-testid="stSidebar"] .sidebar-label {{ color:{TEXT_PRIMARY}!important; }}
    .sidebar-label,.sk-section-label {{ color:{TEXT_MUTED}!important;font-size:12px!important;font-weight:600!important;letter-spacing:.06em!important;text-transform:uppercase; }}
    .sk-shell-brand {{ color:{TEXT_PRIMARY};font-size:22px;font-weight:650;letter-spacing:-.02em;line-height:1.15; }}
    .sk-active-context {{ display:flex;flex-wrap:wrap;gap:8px 20px;margin:.7rem 0 1rem;color:{TEXT_MUTED};font-size:13px; }}
    .sk-title {{ color:{TEXT_PRIMARY};font-size:{FONT_XL};font-weight:650;letter-spacing:-.025em;margin:0; }} .sk-subtitle {{ color:{TEXT_MUTED};font-size:{FONT_BASE};margin:6px 0 18px;max-width:760px; }}
    .sk-context {{ display:flex;align-items:center;gap:6px;margin:0 0 14px;color:{TEXT_MUTED};font-size:12px; }} .sk-context b {{ color:{TEXT_PRIMARY};font-weight:600;letter-spacing:.04em; }}
    .sk-context-header {{ background:linear-gradient(180deg,rgba(88,155,255,.07),rgba(88,155,255,0));border:1px solid {BORDER};border-radius:8px;padding:16px 18px;margin:10px 0 14px; }} .sk-context-header h1 {{ margin:0!important;font-size:24px!important; }} .sk-context-meta {{ color:{TEXT_MUTED}!important;font-size:13px!important;margin-top:5px!important; }}
    .sk-section {{ border-top:1px solid {BORDER};margin:24px 0 12px;padding-top:18px; }} .sk-section h2 {{ margin:3px 0!important; }} .sk-section p {{ color:{TEXT_MUTED}!important;margin:2px 0 0!important;font-size:13px!important; }}
    .sk-insight {{ background:{INSIGHT_HIGHLIGHT};border:1px solid {BORDER};border-left:2px solid {INSIGHT_BORDER};border-radius:6px;padding:12px 14px;margin:12px 0 16px; }} .sk-insight-title {{ color:{TEXT_PRIMARY}!important;font-size:12px!important;font-weight:600!important;letter-spacing:.04em;text-transform:uppercase;margin-bottom:4px; }} .sk-insight-text {{ color:{TEXT_SECONDARY}!important;font-size:13px!important;line-height:1.55; }}
    .sk-scientific-note {{ background:rgba(88,155,255,.06);border-left:2px solid {ACCENT};color:{TEXT_SECONDARY};padding:9px 11px;font-size:13px;margin:10px 0; }}
    .sk-badge {{ display:inline-block;border-radius:4px;padding:2px 6px;font-size:11px!important;font-weight:600!important;margin-right:5px;white-space:nowrap; }}
    [data-testid="stMetric"] {{ background:transparent!important;border:0!important;border-right:1px solid {BORDER}!important;border-radius:0!important;padding:4px 14px!important;box-shadow:none!important; }} [data-testid="stMetric"]:last-child {{ border-right:0!important; }} [data-testid="stMetricLabel"] {{ color:{TEXT_MUTED}!important;font-size:12px!important; }}
    [data-testid="stMetricValue"],[data-testid="stMetricDelta"],[data-testid="stDataFrame"] {{ font-family:"JetBrains Mono","Cascadia Code",monospace!important;font-variant-numeric:tabular-nums; }} [data-testid="stMetricValue"] {{ color:{TEXT_PRIMARY}!important;font-size:20px!important; }}
    .stButton>button,button[data-baseweb="tab"], [data-testid="stSidebar"] label {{ font-family:Inter,"Segoe UI",system-ui,sans-serif!important; }}
    .stButton>button {{ background:{BG_SURFACE}!important;color:{TEXT_SECONDARY}!important;border:1px solid {BORDER_STRONG}!important;border-radius:6px!important;box-shadow:none!important;font-size:13px!important;font-weight:500!important;min-height:2.25rem!important; }} .stButton>button:hover {{ background:{BG_SURFACE_HOVER}!important;color:{TEXT_PRIMARY}!important;border-color:{BORDER_FOCUS}!important; }}
    .stButton>button[kind="primary"] {{ background:linear-gradient(135deg,{ACCENT} 0%,{ACCENT_DARK} 100%)!important;color:#fff!important;border-color:{ACCENT}!important;font-weight:600!important; }} .stButton>button[kind="primary"]:hover {{ background:linear-gradient(135deg,{ACCENT_HOVER} 0%,#518BE1 100%)!important;border-color:{ACCENT_HOVER}!important; }}
    [data-testid="stSidebar"] .stButton>button {{ background:{BG_SURFACE_2}!important;border-color:{BORDER_STRONG}!important; }} [data-testid="stSidebar"] .stButton>button[kind="primary"] {{ min-height:2.65rem!important;border-radius:6px!important;background:linear-gradient(135deg,{ACCENT} 0%,{ACCENT_DARK} 100%)!important;color:white!important;box-shadow:0 1px 2px rgba(0,0,0,.30),0 8px 24px rgba(0,0,0,.12)!important; }}
    [data-baseweb="select"] > div,[data-baseweb="input"] > div,[data-baseweb="base-input"],[data-testid="stTextInput"] input,[data-testid="stTextArea"] textarea {{ background:{BG_INPUT}!important;color:{TEXT_PRIMARY}!important;border-color:{BORDER_STRONG}!important;border-radius:6px!important; }} [data-baseweb="radio"] div[role="radio"][aria-checked="true"],[data-testid="stSlider"] [role="slider"] {{ background:{ACCENT}!important;border-color:{ACCENT}!important; }}
    button[data-baseweb="tab"] {{ color:{TEXT_MUTED}!important;font-size:13px!important;font-weight:500!important; }} button[data-baseweb="tab"][aria-selected="true"] {{ color:{ACCENT_BRIGHT}!important;border-bottom-color:{ACCENT}!important; }}
    div[data-testid="stExpander"] {{ background:{BG_SURFACE_2}!important;border:1px solid {BORDER}!important;border-radius:6px!important; }} [data-testid="stDataFrame"],[data-testid="stDataEditor"] {{ border:1px solid {BORDER}!important;border-radius:6px!important; }} [data-testid="stDataFrame"] [role="columnheader"] {{ background:#111624!important;color:{TEXT_SECONDARY}!important;font-size:12px!important; }} [data-testid="stDataFrame"] [role="gridcell"] {{ color:{TEXT_SECONDARY}!important;font-size:12px!important;border-color:#202A40!important; }}
    [data-testid="stAlert"] {{ background:{BG_SURFACE_2}!important;border:1px solid {BORDER}!important;border-left:2px solid {ACCENT}!important;border-radius:6px!important; }} .section-header,.validation-card,.sk-verdict,.analysis-card,.metric-card {{ background:{BG_SURFACE}!important;border-color:{BORDER}!important;box-shadow:none!important; }} .section-header-text {{ color:{TEXT_PRIMARY}!important; }} .section-header-sub,.table-limit-note,.stale-banner {{ color:{TEXT_MUTED}!important;background:{BG_SURFACE_2}!important;border-color:{BORDER}!important; }} .sk-readout-wrap,.sk-pipeline {{ background:linear-gradient(145deg,{BG_SURFACE} 0%,{BG_SURFACE_2} 100%)!important;border-color:{BORDER}!important;box-shadow:none!important; }} .sk-readout-title,.sk-readout-value,.sk-verdict-text {{ color:{TEXT_PRIMARY}!important; }} .sk-readout-label,.sk-eyebrow,.sk-pipeline-stage-label {{ color:{TEXT_MUTED}!important; }}
    /* Brand assets stay source-faithful and independent of the application theme. */
    img[alt="Sophiark logo"],.top-header-logo-wrap img {{ filter:none!important;opacity:1!important;mix-blend-mode:normal!important;object-fit:contain!important;box-shadow:none!important; }}
    </style>""", unsafe_allow_html=True)
    st.markdown("<style>" + _stylesheet("research_workspace.css") + "</style>", unsafe_allow_html=True)


def render_badge(text: str, kind: str = "kanit_zayif") -> str:
    colors = {"tier_kritik": TIER_KRITIK, "tier_orta": TIER_ORTA, "tier_dusuk": TIER_DUSUK, "beta": BETA_BADGE, "kanit_guclu": KANIT_GUCLU, "kanit_zayif": KANIT_ZAYIF, "dikkat_cekici": TIER_ORTA}
    color = colors.get(kind, KANIT_ZAYIF)
    return f'<span class="sk-badge" style="color:{color};border:1px solid {color};">{text}</span>'
