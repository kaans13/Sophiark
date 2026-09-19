"""Data-driven, conservative biological-theme vocabulary.

Every match represents annotation overlap only.  No theme is a claim of a
pathway's activation, inhibition, causal mechanism, or clinical relevance.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class ThemeDefinition:
    identifier: str
    label: str
    mechanism_class: str
    keywords: tuple[str, ...]


THEME_ONTOLOGY: tuple[ThemeDefinition, ...] = (
    ThemeDefinition("chloride_anion_transport", "Klorür / anyon taşınması", "iyon taşınması", ("chloride", "anion transport", "anion channel")),
    ThemeDefinition("bicarbonate_transport", "Bikarbonat taşınması", "iyon taşınması", ("bicarbonate",)),
    ThemeDefinition("water_transport", "Su taşınması", "membran taşınması", ("water transport", "aquaporin")),
    ThemeDefinition("cell_volume_osmotic", "Hücre hacmi / ozmotik düzenleme", "homeostaz", ("osmotic", "cell volume", "cellular volume", "volume regulation")),
    ThemeDefinition("ph_ion_homeostasis", "pH / iyon homeostazı", "homeostaz", ("ph homeostasis", "regulation of ph", "ion homeostasis", "proton transport")),
    ThemeDefinition("membrane_transport", "Membran taşınması", "membran taşınması", ("membrane transport", "transmembrane transport", "ion transport", "transporter activity", "channel activity")),
    ThemeDefinition("calcium_handling", "Kalsiyum işlenmesi", "iyon sinyalleşmesi", ("calcium ion", "calcium channel", "calcium release", "calcium homeostasis")),
    ThemeDefinition("excitation_contraction", "Uyarılma-kasılma eşleşmesi", "kas fizyolojisi", ("excitation-contraction", "excitation contraction", "cardiac muscle contraction")),
    ThemeDefinition("sarcoplasmic_reticulum_calcium", "Sarkoplazmik retikulum kalsiyum bağlamı", "hücresel bölme", ("sarcoplasmic reticulum", "junctional sarcoplasmic reticulum")),
    ThemeDefinition("cholesterol_handling", "Kolesterol işlenmesi", "lipit homeostazı", ("cholesterol", "sterol")),
    ThemeDefinition("lipoprotein_handling", "Lipoprotein işlenmesi", "lipit taşınması", ("lipoprotein", "low-density lipoprotein", "ldl receptor")),
)
