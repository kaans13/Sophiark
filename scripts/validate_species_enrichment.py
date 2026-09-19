"""Network-backed human/mouse enrichment smoke test; writes only validation JSON."""

from __future__ import annotations

import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.services.enrichment import run_post_simulation_enrichment


human_mapping = pd.read_csv("data/processed/ensp_with_symbols.csv")
human_universe = human_mapping["Symbol"].dropna().astype(str).tolist()
from src.mouse.config import MOUSE_SYMBOL_MAP
mouse_universe = [str(value) for value in MOUSE_SYMBOL_MAP.values() if value]

result = {}
for key, selected, universe, is_mouse in (
    ("human", ["TP53", "MDM2", "ATM", "CHEK2"], human_universe, False),
    ("mouse", ["Trp53", "Mdm2", "Atm", "Chek2"], mouse_universe, True),
):
    frame, notices = run_post_simulation_enrichment(
        selected, is_mouse=is_mouse, background_symbols=universe,
    )
    result[key] = {
        "rows": len(frame),
        "notices": notices,
        "columns": list(frame.columns),
        "sources": sorted(frame["Kaynak"].unique().tolist()) if not frame.empty else [],
        "species": sorted(frame["Species"].unique().tolist()) if not frame.empty else [],
        "background_sizes": sorted(frame["Background_Universe_Size"].unique().tolist()) if not frame.empty else [],
    }

output = Path("outputs/scientific_validation/species_enrichment.json")
output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
