"""Annotation- and symbol-based pattern detection for read-only interpretation."""

from __future__ import annotations

from collections import defaultdict
import re

import pandas as pd


THEME_RULES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Klorür / anyon taşınması", ("chloride", "anion transport", "anion channel")),
    ("Bikarbonat taşınması", ("bicarbonate",)),
    ("Su taşınması", ("water transport", "aquaporin")),
    ("Hücre hacmi / ozmotik düzenleme", ("osmotic", "cell volume", "cellular volume", "volume regulation")),
    ("pH / iyon homeostazı", ("ph homeostasis", "regulation of ph", "ion homeostasis", "proton transport")),
    ("Membran taşınması", ("membrane transport", "transmembrane transport", "ion transport", "transporter activity", "channel activity")),
)

FAMILY_RULES: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("ANO", re.compile(r"^ANO\d", re.I)),
    ("AQP", re.compile(r"^AQP\d", re.I)),
    ("SLC4", re.compile(r"^SLC4", re.I)),
    ("SLC9", re.compile(r"^SLC9", re.I)),
    ("SLC13", re.compile(r"^SLC13", re.I)),
    ("SLC26", re.compile(r"^SLC26", re.I)),
    ("LRRC8", re.compile(r"^LRRC8", re.I)),
    ("CLCN", re.compile(r"^CLCN", re.I)),
    ("BEST", re.compile(r"^BEST\d", re.I)),
    ("CA", re.compile(r"^CA\d", re.I)),
)

ANNOTATION_COLUMNS = (
    "GO Biyolojik Süreç", "GO Moleküler İşlev", "GO Hücresel Bileşen",
    "GO_CC_Terimleri", "GO_MF_Terimleri", "Protein_Adi", "MyGene Adı",
)


def _symbols(frame: pd.DataFrame) -> pd.Series:
    column = "Symbol" if "Symbol" in frame.columns else "gene"
    return frame.get(column, pd.Series(dtype=str)).fillna("").astype(str).str.strip()


def annotation_text(row: pd.Series) -> str:
    return " | ".join(str(row.get(column, "")) for column in ANNOTATION_COLUMNS if pd.notna(row.get(column))) .casefold()


def detect_functional_themes(candidates: pd.DataFrame, *, max_genes: int = 10) -> list[dict[str, object]]:
    """Detect repeated terms only where an existing candidate annotation supports it."""
    found: dict[str, list[str]] = defaultdict(list)
    for _, row in candidates.iterrows():
        symbol = str(row.get("Symbol") or row.get("gene") or "").strip()
        if not symbol:
            continue
        text = annotation_text(row)
        for name, keywords in THEME_RULES:
            if any(keyword in text for keyword in keywords) and symbol not in found[name]:
                found[name].append(symbol)
    return [
        {"theme": name, "genes": genes[:max_genes], "count": len(genes)}
        for name, _ in THEME_RULES
        if (genes := found.get(name, []))
    ]


def detect_family_patterns(candidates: pd.DataFrame, *, min_members: int = 2) -> list[dict[str, object]]:
    """Report repeated symbol families; this does not infer a pathway state."""
    symbols = list(dict.fromkeys(symbol for symbol in _symbols(candidates) if symbol))
    patterns: list[dict[str, object]] = []
    for family, matcher in FAMILY_RULES:
        members = [symbol for symbol in symbols if matcher.match(symbol)]
        if len(members) >= min_members:
            patterns.append({"family": family, "genes": members, "count": len(members)})
    return patterns
