from __future__ import annotations

import json

import pandas as pd

from src.unified.service import UnifiedResearchReport, UnifiedStatus, export_unified_report


def test_unified_json_export_sanitizes_machine_paths_without_changing_machine_schema(tmp_path):
    report = UnifiedResearchReport(
        target="ENSP00000000001", tissue="Lung", status=UnifiedStatus.COMPLETE,
        classic_status="COMPLETE", directed_status="AVAILABLE",
        candidates=pd.DataFrame({"symbol": ["TP53"], "entity_id": ["ENSP00000000001"]}),
        classic_raw=pd.DataFrame({"gene": ["ENSP00000000001"]}),
        directed_raw=pd.DataFrame({"entity": ["ENSP00000000001"]}),
        complex_context=pd.DataFrame(), agreement={},
        provenance={"source": "C:/Users/example/private/input.tsv"},
        context_availability={"evidence_provenance": {"index_path": "C:/Users/example/private/index.sqlite"}},
        errors={},
    )

    paths = export_unified_report(report, tmp_path)
    payload = json.loads(paths["report"].read_text(encoding="utf-8"))

    assert payload["target"] == "ENSP00000000001"
    assert payload["provenance"]["source"] == "input.tsv"
    assert payload["context_availability"]["evidence_provenance"]["index_path"] == "index.sqlite"
