"""Sparse, UI-framework-neutral rendering for the optional Disease Context."""

from __future__ import annotations

import io
import zipfile
from typing import Any, Callable

import pandas as pd

from .benchmark import benchmark_disease_reference
from .service import DiseaseBankService, DiseaseLoad


def _table_downloads(ui: Any, tables: dict[str, pd.DataFrame], *, stem: str) -> None:
    """Expose every Disease Context table in XLSX and a UTF-8 CSV ZIP."""

    buffer = io.BytesIO()
    with pd.ExcelWriter(buffer, engine="openpyxl") as writer:
        for name, frame in tables.items():
            frame.to_excel(writer, index=False, sheet_name=name[:31])
    ui.download_button("Tüm hastalık bağlamı tablolarını indir (Excel)", data=buffer.getvalue(), file_name=f"{stem}.xlsx", mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
    archive = io.BytesIO()
    with zipfile.ZipFile(archive, "w", zipfile.ZIP_DEFLATED) as zipped:
        for name, frame in tables.items():
            slug = "_".join(name.casefold().split())
            zipped.writestr(f"{slug}.csv", frame.to_csv(index=False).encode("utf-8-sig"))
    ui.download_button("Tüm hastalık bağlamı tablolarını indir (CSV ZIP)", data=archive.getvalue(), file_name=f"{stem}.zip", mime="application/zip")


def _relationship_rows(research_bundle: Any, reference) -> pd.DataFrame:
    """Read existing Explorer assertions; never create relationship evidence."""

    reference_ids = reference.canonical_proteins
    rows: list[dict[str, object]] = []
    for assertion in tuple(getattr(getattr(research_bundle, "relationships", None), "assertions", ()) or ()):
        subject = getattr(assertion, "subject", None)
        object_entity = getattr(assertion, "object", None)
        subject_id = str(getattr(subject, "canonical_id", "")).strip()
        object_id = str(getattr(object_entity, "canonical_id", "")).strip()
        if subject_id not in reference_ids and object_id not in reference_ids:
            continue
        rows.append({
            "Relationship ID": getattr(assertion, "assertion_id", ""),
            "Relationship": str(getattr(getattr(assertion, "relationship_type", None), "value", "")),
            "Subject": subject_id, "Object": object_id,
        })
    return pd.DataFrame(rows, columns=["Relationship ID", "Relationship", "Subject", "Object"])


def render_disease_context(snapshot, *, service_factory: Callable[[], DiseaseBankService], ui: Any, research_bundle: Any = None) -> None:
    """Render only after analysis; the external source is explicit and optional."""

    key = f"disease_bank_{snapshot.snapshot_id}"
    with ui.expander("Hastalık bağlamı / Benchmark", expanded=False):
        ui.caption("Yalnızca dış referans karşılaştırmasıdır; hastalık verisi ağ hesabını değiştirmez.")
        if not ui.session_state.get(f"{key}_enabled", False):
            if ui.button("Hastalık bağlamını aç", key=f"{key}_open"):
                ui.session_state[f"{key}_enabled"] = True
            return
        # Resolver and cache are constructed only after the user explicitly
        # opens this optional context; no external request occurs here.
        service = service_factory()
        term = ui.text_input("Hastalık ara", key=f"{key}_term", placeholder="Parkinson hastalığı")
        if ui.button("Hastalık ara", key=f"{key}_search"):
            result = service.search(term)
            ui.session_state[f"{key}_search_result"] = result
            ui.session_state.pop(f"{key}_reference", None)
        search_result = ui.session_state.get(f"{key}_search_result")
        if isinstance(search_result, DiseaseLoad):
            if search_result.message:
                ui.info(f"Disease Bank: {search_result.status} — {search_result.message}")
            return
        if not search_result:
            return
        options = {f"{item.disease_name} ({item.disease_id})": item.disease_id for item in search_result}
        selected_label = ui.selectbox("Hastalık", tuple(options), key=f"{key}_selection")
        if ui.button("Hastalık referansını yükle", key=f"{key}_load"):
            ui.session_state[f"{key}_reference"] = service.load_reference(options[selected_label])
        load = ui.session_state.get(f"{key}_reference")
        if not isinstance(load, DiseaseLoad):
            return
        if load.reference is None:
            ui.info(f"Disease Bank: {load.status}" + (f" — {load.message}" if load.message else ""))
            return
        reference = load.reference
        mapping = reference.mapping_summary
        ui.caption(f"{reference.disease.disease_name} · {reference.disease.disease_id} · source {reference.source_status}")
        ui.dataframe(pd.DataFrame([{
            "Known associated genes": mapping.total, "Successfully mapped": mapping.mapped,
            "Unmapped": mapping.unmapped, "Ambiguous": mapping.ambiguous,
        }]), hide_index=True, width="stretch")
        report = snapshot.report.to_frame()
        redistribution = snapshot.signed_redistribution.to_frame() if snapshot.signed_redistribution is not None else None
        benchmark = benchmark_disease_reference(
            reference, report=report, redistribution=redistribution,
            perturbed_proteins=service.resolve_perturbed_proteins(snapshot.targets), population_size=snapshot.tested_count,
        )
        ui.dataframe(pd.DataFrame(benchmark.summary_rows()), hide_index=True, width="stretch")
        recalls = pd.DataFrame([{
            "Metric": f"Recall@{metric.k}", "Returned": metric.returned_count, "Overlap": metric.overlap_count,
            "Evaluation reference": metric.evaluation_reference_count, "Recall": metric.recall,
            "Expected random overlap": metric.expected_random_overlap, "Enrichment ratio": metric.enrichment_ratio,
        } for metric in benchmark.recalls])
        ui.dataframe(recalls, hide_index=True, width="stretch")
        comparison = pd.DataFrame({
            "Known Disease-Associated Response": pd.Series(benchmark.known_response),
            "Disease-Associated Network Loss": pd.Series(benchmark.known_network_loss),
            "Disease-Associated Redistribution": pd.Series(benchmark.known_redistribution),
            "Network candidates outside the reference set": pd.Series(benchmark.network_candidates_outside_reference),
        })
        ui.dataframe(comparison, hide_index=True, width="stretch")
        mappings = pd.DataFrame([{
            "Gene symbol": item.association.gene_symbol, "Ensembl gene ID": item.association.ensembl_gene_id,
            "Resolved protein": item.resolution.entity.canonical_id if item.resolution.entity is not None else "",
            "Resolution status": item.resolution.status.value, "Query used": item.query_used,
            "Association score": item.association.association_score,
        } for item in reference.resolved])
        relationships = _relationship_rows(research_bundle, reference)
        if not relationships.empty:
            ui.dataframe(relationships, hide_index=True, width="stretch")
        _table_downloads(ui, {
            "Summary": pd.DataFrame(benchmark.summary_rows()), "Recall": recalls,
            "Comparison": comparison, "Mapping": mappings, "Research Relationships": relationships,
        }, stem=f"disease_benchmark_{reference.disease.disease_id}")
        ui.caption(benchmark.limitation)


__all__ = ["render_disease_context"]
