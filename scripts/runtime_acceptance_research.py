"""Build the active offline Research Explorer from an already executed result."""
from __future__ import annotations

import json
from pathlib import Path
import sys
import time

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.research.facade import ResearchConfig, ResearchContextService, build_research_snapshot
from src.research.ui import build_research_export_tables


def _count(value) -> int:
    try:
        return len(value)
    except TypeError:
        return 0


def main() -> None:
    report = pd.read_csv(ROOT / "outputs" / "runtime_acceptance" / "secondary_smoke" / "multi_target" / "enfeksiyon_sok_dalgasi_raporu.csv")
    snapshot = build_research_snapshot(
        report, None, None, species="Homo sapiens", taxon_id=9606, tissue="Lung",
        targets=("ENSP00000003084", "ENSP00000447297"), attenuation=0.001,
        selection_mode="top_n", test_limit=None, tested_count=None, returned_count=len(report), top_n=20,
        threshold_parameters={"engine": "classic", "scientific_top_n_metric": "Delta_PageRank_Pct"},
    )
    config = ResearchConfig(
        research_context_enabled=True, external_context_enabled=False, directed_signaling_enabled=True,
        directed_signaling_max_depth=4, provider_uniprot_enabled=False, provider_interpro_enabled=False,
        provider_quickgo_enabled=False, provider_reactome_enabled=False, provider_europepmc_enabled=False,
    )
    started = time.perf_counter()
    service = ResearchContextService.from_project_root(ROOT, config=config)
    bundle = service.build_offline_bundle(snapshot)
    tables = build_research_export_tables(bundle, config=config)
    payload = {
        "seconds": time.perf_counter() - started,
        "snapshot_id": snapshot.snapshot_id,
        "errors": list(bundle.errors), "notices": list(bundle.notices),
        "observations": _count(bundle.observation_result.observations),
        "relationships": _count(bundle.relationships.assertions),
        "semantic_records": sum(_count(getattr(bundle.semantic_result, field, ())) for field in ("target_records", "positive_redistribution", "network_losses", "fdr_supported")),
        "local_entities": _count(bundle.local_context.entities),
        "table_counts": {name: len(rows) for name, rows in tables.items()},
        "identity_fields_in_general_view": sorted({key for row in tables.get("General View", ()) for key in row if key in {"Gen", "ENSP", "Protein ID"}}),
    }
    target = ROOT / "outputs" / "runtime_acceptance" / "research_validation.json"
    target.write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    print(json.dumps(payload, indent=2, default=str))


if __name__ == "__main__":
    main()
