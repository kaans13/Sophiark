"""One-off, file-backed UI refactor validation for the human Liver network."""
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.core.analysis_runtime import prepare_active_network, runtime_for_selection

out = Path("outputs/ui_refactor_validation/liver_network.json")
out.parent.mkdir(parents=True, exist_ok=True)
context = runtime_for_selection(is_mouse=False, tissue="Liver")
prepared = prepare_active_network(context, forced_genes=None, bc_sample_sources=300)
out.write_text(json.dumps({"tissue": "Liver", "nodes": prepared.graph.vcount(), "edges": prepared.graph.ecount()}), encoding="utf-8")
