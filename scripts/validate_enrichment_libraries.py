"""Manual smoke check for currently advertised Enrichr libraries."""

from __future__ import annotations

import json
from pathlib import Path

import gseapy as gp


result = {}
for organism in ("Human", "Mouse"):
    try:
        libraries = gp.get_library_name(organism=organism)
        result[organism] = {
            "count": len(libraries),
            "relevant": [name for name in libraries if "KEGG" in name or "GO_Biological" in name],
        }
    except Exception as error:
        result[organism] = {"error": type(error).__name__, "message": str(error)}

output = Path("outputs") / "scientific_validation" / "enrichment_libraries.json"
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
