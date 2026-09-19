"""Build one real cross-tissue graph and persist normalization metadata."""

from __future__ import annotations

import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.config import DB_PATH, HIGH_CONF_THRESHOLD
from src.graph_engine import build_high_conf_subgraph


graph = build_high_conf_subgraph(
    DB_PATH, HIGH_CONF_THRESHOLD, hedef_doku="Liver",
    tissue_normalization_mode="cross_tissue_comparable",
)
result = {
    "tissue": "Liver", "nodes": graph.vcount(), "edges": graph.ecount(),
    "normalization": graph["tissue_normalization"],
}
path = Path("outputs/scientific_validation/cross_tissue_normalization.json")
path.write_text(json.dumps(result, ensure_ascii=False, indent=2), encoding="utf-8")
