"""Validate export-to-bundle-to-history integrity from a completed Unified export."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from src.product.history import HistoryService
from src.unified.service import UnifiedResearchReport, UnifiedStatus


def main() -> None:
    source = ROOT / "outputs" / "runtime_acceptance" / "unified_panel" / "CFTR"
    payload = json.loads((source / "unified_report.json").read_text(encoding="utf-8"))
    report = UnifiedResearchReport(
        target=payload["target"], tissue=payload["tissue"], status=UnifiedStatus(payload["status"]),
        classic_status=payload["classic_status"], directed_status=payload["directed_status"],
        candidates=pd.read_csv(source / "unified_candidates.csv"),
        classic_raw=pd.read_csv(source / "unified_classic_raw.csv"),
        directed_raw=pd.read_csv(source / "unified_directed_raw.csv"),
        complex_context=pd.read_csv(source / "unified_complex_context.csv"),
        agreement=payload["agreement"], provenance=payload["provenance"],
        context_availability=payload["context_availability"], errors=payload["errors"],
    )
    history_dir = ROOT / "outputs" / "runtime_acceptance" / "bundle_history"
    service = HistoryService(history_dir=history_dir)
    entry = service.save(report)
    saved = service.open_entry(entry.analysis_id)
    rebuild = service.rebuild_history_index()
    checks = {
        "candidate_identity_equal": report.candidates["entity_id"].astype(str).equals(saved.report.candidates["entity_id"].astype(str)),
        "candidate_rows_equal": len(report.candidates) == len(saved.report.candidates),
        "classic_rows_equal": len(report.classic_raw) == len(saved.report.classic_raw),
        "directed_rows_equal": len(report.directed_raw) == len(saved.report.directed_raw),
        "complex_rows_equal": len(report.complex_context) == len(saved.report.complex_context),
        "statuses_equal": (report.status == saved.report.status and report.classic_status == saved.report.classic_status and report.directed_status == saved.report.directed_status),
        "provenance_scientific_fields_equal": all(
            report.provenance.get(key) == saved.report.provenance.get(key)
            for key in ("target", "tissue", "classic_config_fingerprint", "directed_config_fingerprint", "run_manifest")
        ),
        "evidence_context_status_equal": all(
            report.context_availability.get(key, {}).get("status") == saved.report.context_availability.get(key, {}).get("status")
            for key in ("evidence_provenance", "corum_complex_context")
        ),
        "portable_paths_sanitized": (
            saved.report.context_availability.get("evidence_provenance", {}).get("index_path") == "string_evidence_v12_9606.sqlite"
        ),
    }
    result = {
        "entry": asdict(entry), "bundle_path": str(history_dir / entry.bundle_filename),
        "checks": checks, "rebuild": asdict(rebuild),
    }
    destination = ROOT / "outputs" / "runtime_acceptance" / "bundle_history_check.json"
    destination.write_text(json.dumps(result, indent=2, default=str), encoding="utf-8")
    print(json.dumps(result, indent=2, default=str))


if __name__ == "__main__":
    main()
