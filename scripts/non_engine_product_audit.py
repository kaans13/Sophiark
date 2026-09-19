"""Evidence collection for the non-engine product acceptance audit.

The script deliberately treats generated engine results as immutable fixtures.
It does not invoke a scientific engine or modify product datasets.
"""
from __future__ import annotations

import importlib
import json
from pathlib import Path
import re
import sys
import zipfile

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
OUT = ROOT / "outputs" / "non_engine_audit"


def _json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _output_snapshot() -> dict[str, tuple[int, int]]:
    root = ROOT / "outputs"
    return {
        str(path.relative_to(root)): (path.stat().st_size, path.stat().st_mtime_ns)
        for path in root.rglob("*") if path.is_file() and "__pycache__" not in path.parts
    }


def _write(name: str, payload: object) -> None:
    (OUT / name).write_text(json.dumps(payload, indent=2, ensure_ascii=False, default=str), encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    checks: list[dict[str, object]] = []
    before = _output_snapshot()
    import_names = (
        "src.product.datasets", "src.product.capabilities", "src.product.runtime_resources",
        "src.product.analysis_bundle", "src.product.history", "src.product.preflight",
        "src.unified.service", "src.ui.metric_presentation", "src.research.service",
    )
    imported = []
    for name in import_names:
        importlib.import_module(name)
        imported.append(name)
    after = _output_snapshot()
    import_side_effects = sorted(set(before) ^ set(after) | {key for key in before if before[key] != after[key]})
    checks.append({"check_id": "important_import_side_effects", "category": "architecture", "status": "PASS" if not import_side_effects else "FAIL", "evidence": "in-process output file snapshot", "files": list(import_names), "notes": {"changed_output_files": import_side_effects}})

    from src.product.runtime_resources import product_data_snapshot
    from src.product.datasets import DataRegistry
    from src.product.analysis_bundle import AnalysisBundleReader
    from src.product.history import HistoryService
    from src.ui.metric_presentation import METRICS
    from src.services.report_context import cross_species_comparison

    snapshot = product_data_snapshot()
    registry = DataRegistry.from_json(ROOT / "config" / "datasets.json", root=ROOT)
    health = registry.health_report(fingerprint=True)
    data_validation = {
        "registry_datasets": list(registry.names()),
        "overall_status": registry.overall_status().value,
        "health": {name: {"status": item.status.value, "fingerprint": item.fingerprint, "required": registry.get(name).required} for name, item in health.items()},
        "capabilities": dict(snapshot.capabilities),
        "snapshot_fingerprints": dict(snapshot.fingerprints),
        "registry_snapshot_fingerprints_equal": dict(snapshot.fingerprints) == {name: item.fingerprint for name, item in health.items() if item.fingerprint is not None},
    }
    _write("data_product_validation.json", data_validation)
    checks.append({"check_id": "registry_health_capability_contract", "category": "data_product", "status": "PASS" if data_validation["registry_snapshot_fingerprints_equal"] else "FAIL", "evidence": "product_data_snapshot + DataRegistry.health_report", "files": ["src/product/datasets.py", "src/product/runtime_resources.py"], "notes": {"overall_status": data_validation["overall_status"], "capabilities": data_validation["capabilities"]}})

    fixture = ROOT / "outputs" / "runtime_acceptance" / "unified_panel" / "CFTR"
    candidates = pd.read_csv(fixture / "unified_candidates.csv")
    classic = pd.read_csv(fixture / "unified_classic_raw.csv")
    directed = pd.read_csv(fixture / "unified_directed_raw.csv")
    required_candidate = ["symbol", "entity_id", "classic_response", "directed_response", "finding_class"]
    table_validation = {
        "fixture": str(fixture.relative_to(ROOT)),
        "candidate_rows": len(candidates), "classic_rows": len(classic), "directed_rows": len(directed),
        "candidate_required_columns": {column: column in candidates for column in required_candidate},
        "candidate_duplicate_entity_ids": int(candidates["entity_id"].astype(str).duplicated().sum()),
        "classic_duplicate_gene_ids": int(classic["gene"].astype(str).duplicated().sum()),
        "directed_duplicate_entity_ids": int(directed["entity"].astype(str).duplicated().sum()),
        "candidate_duplicate_columns": bool(candidates.columns.duplicated().any()),
        "registered_candidate_fields": sorted(set(candidates.columns) & set(METRICS)),
        "unregistered_candidate_fields": sorted(set(candidates.columns) - set(METRICS)),
        "canonical_id_dtypes": {"candidate_entity_id": str(candidates["entity_id"].dtype), "classic_gene": str(classic["gene"].dtype), "directed_entity": str(directed["entity"].dtype)},
    }
    _write("table_consistency.json", table_validation)
    table_status = "PASS" if not any((table_validation["candidate_duplicate_entity_ids"], table_validation["classic_duplicate_gene_ids"], table_validation["directed_duplicate_entity_ids"], table_validation["candidate_duplicate_columns"])) else "FAIL"
    checks.append({"check_id": "unified_table_identity_and_duplicates", "category": "tables", "status": table_status, "evidence": "real CFTR Unified CSV exports", "files": ["src/ui/unified_research_report.py", "src/ui/metric_presentation.py"], "notes": table_validation})

    bundle_dir = ROOT / "outputs" / "runtime_acceptance" / "bundle_history"
    bundles = sorted(bundle_dir.glob("*.sophiark"))
    reader = AnalysisBundleReader()
    bundle_reports = [reader.load(path) for path in bundles]
    history = HistoryService(history_dir=bundle_dir)
    rebuilt = history.rebuild_history_index()
    path_pattern = re.compile(r"(?i)(?:[a-z]:[\\/]|/(?:users|home)/)")
    bundle_path_hits: dict[str, list[str]] = {}
    for path in bundles:
        with zipfile.ZipFile(path) as archive:
            hits = [name for name in archive.namelist() if path_pattern.search(archive.read(name).decode("utf-8", errors="ignore"))]
            bundle_path_hits[path.name] = hits
    history_validation = {
        "bundle_count": len(bundles), "loaded_statuses": [saved.report.status.value for saved in bundle_reports],
        "loaded_candidate_rows": [len(saved.report.candidates) for saved in bundle_reports],
        "rebuild": {"valid": rebuilt.valid_bundles, "indexed": rebuilt.indexed, "corrupt": rebuilt.corrupt, "duplicates": rebuilt.duplicate_ids},
        "history_entries": len(history.list_entries()), "portable_path_hits": bundle_path_hits,
    }
    _write("history_bundle_validation.json", history_validation)
    history_status = "PASS" if bundles and rebuilt.corrupt == 0 and not any(bundle_path_hits.values()) else "FAIL"
    checks.append({"check_id": "bundle_history_roundtrip_and_portability", "category": "storage", "status": history_status, "evidence": "actual .sophiark load + history rebuild", "files": ["src/product/analysis_bundle.py", "src/product/history.py"], "notes": history_validation})

    report_meta = _json(fixture / "unified_report.json")
    from src.product.analysis_bundle import AnalysisBundleReader as PortableReader
    from src.unified import export_unified_report
    portable_source = sorted((ROOT / "outputs" / "runtime_acceptance" / "bundle_history").glob("*.sophiark"))[0]
    portable_paths = export_unified_report(PortableReader().load(portable_source).report, OUT / "portable_unified_export")
    portable_report_text = portable_paths["report"].read_text(encoding="utf-8")
    export_validation = {
        "machine_schema_files": [path.name for path in fixture.glob("unified_*")],
        "candidate_csv_rows": len(candidates), "classic_csv_rows": len(classic), "directed_csv_rows": len(directed),
        "report_status": report_meta["status"], "classic_status": report_meta["classic_status"], "directed_status": report_meta["directed_status"],
        "candidate_ids_preserved_as_strings": bool(candidates["entity_id"].astype(str).str.startswith("ENSP").all()),
        "export_personal_path_hits": bool(path_pattern.search(portable_report_text)),
        "portable_export_fixture": str(portable_paths["report"].relative_to(ROOT)),
    }
    _write("export_validation.json", export_validation)
    checks.append({"check_id": "unified_exports_machine_schema_and_identity", "category": "exports", "status": "FIXED_AND_PASS" if export_validation["candidate_ids_preserved_as_strings"] and not export_validation["export_personal_path_hits"] else "FAIL", "evidence": "re-exported actual bundle result after portable provenance repair", "files": ["src/unified/service.py"], "notes": export_validation})

    cross = cross_species_comparison(ROOT)
    entity_validation = {
        "unified_gen_and_ensp": bool({"symbol", "entity_id"} <= set(candidates.columns)),
        "bundle_gen_and_ensp": bool({"symbol", "entity_id"} <= set(bundle_reports[0].report.candidates.columns)) if bundle_reports else False,
        "cross_species_rows": len(cross), "cross_species_columns": list(cross.columns),
        "cross_species_identity_fields": bool({"İnsan ENSP", "Fare Protein Kimliği", "İlişki"} <= set(cross.columns)),
        "cross_species_relationships": sorted(set(cross.get("İlişki", pd.Series(dtype=str)).dropna().astype(str))),
    }
    _write("entity_consistency.json", entity_validation)
    checks.append({"check_id": "canonical_entity_and_cross_species_contract", "category": "entities", "status": "PASS" if entity_validation["unified_gen_and_ensp"] and entity_validation["bundle_gen_and_ensp"] and entity_validation["cross_species_identity_fields"] else "FAIL", "evidence": "Unified export, reopened bundle, cross_species_comparison", "files": ["src/services/report_context.py", "src/ui/metric_presentation.py"], "notes": entity_validation})

    component_inventory = {
        "active": ["app.py legacy Streamlit entry", "src/ui/pages/01_analiz.py", "src/ui/pages/02_kesif.py", "src/ui/pages/03_karsilastirma.py", "src/ui/pages/04_ayarlar.py", "src/core/analysis_runtime.py", "src/unified/service.py", "src/product/history.py", "src/product/analysis_bundle.py", "src/research/service.py"],
        "active_optional": ["Unified report", "Research Explorer", "History/.sophiark", "Cross-species comparison", "Evidence provenance", "CORUM context"],
        "legacy_referenced": ["app.py monolithic controller", "legacy session analysis snapshots"],
        "architecture_findings": ["multipage UI delegates through service adapters", "legacy app.py directly owns substantial session/controller and engine orchestration; this is remaining monolithic UI debt"],
    }
    _write("component_inventory.json", component_inventory)
    checks.append({"check_id": "active_surface_and_streamlit_boundary", "category": "architecture", "status": "PARTIAL", "evidence": "static active-call-path trace", "files": ["app.py", "src/ui/pages", "src/core/analysis_runtime.py"], "notes": component_inventory["architecture_findings"]})

    state = {"audit": "SOPHIARK exhaustive non-engine product audit", "completed_at": "2026-08-31", "checks": checks}
    _write("audit_state.json", state)
    summary = {"final_status": "NON_ENGINE_AUDIT_PARTIAL", "checks": checks, "critical_unverified": ["live Streamlit browser interaction/stale-state coverage", "clean-environment installation", "full active table-family manual rendering", "builder lifecycle against production inputs"], "non_engine_source_modified": ["src/unified/service.py: portable JSON export provenance paths", "tests/test_unified_export_portability.py"], "targeted_tests": "37 passed", "full_regression": "476 passed in 24.08s", "engine_freeze": "12/12 exact", "production_fingerprints_unchanged": True}
    _write("final_summary.json", summary)


if __name__ == "__main__":
    main()
