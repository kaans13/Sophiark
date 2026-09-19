"""Presentation-layer constants and centralized Streamlit state keys."""

from __future__ import annotations


class SessionKeys:
    ACTIVE_MODE = "ui_active_mode"
    LANDING_MODE = "ui_landing_mode"
    NAVIGATION = "ui_navigation"
    DEVELOPER_MODE = "ui_developer_mode"


ANALYSIS_MODES = ("Perturbasyon Analizi", "Target Stress Search", "Compensation Analysis · BETA")
SUPPORTED_SPECIES_COUNT = 2
ENRICHMENT_SOURCE_COUNT = 2  # KEGG and GO are the currently rendered source families.
# Yerel STRING snapshot'ında combined_score >= 500 için doğrulanmış sayılar.
# Landing ekranında 6.8M satırlı veritabanını her açılışta yeniden taramamak için
# veri snapshot metadata'sı olarak merkezi tutulur.
HUMAN_ACTIVE_PROTEIN_COUNT = 19_220
HUMAN_ACTIVE_EDGE_COUNT = 569_217
