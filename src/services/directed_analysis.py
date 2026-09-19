"""Failure-isolated application adapter for the optional directed engine."""

from __future__ import annotations

from dataclasses import dataclass
from collections import OrderedDict
from functools import lru_cache
import logging
from pathlib import Path
from threading import Lock

import pandas as pd

from src.directed import (
    DirectedCalculationResult, DirectedEngineStatus, EngineComparisonResult,
    build_hybrid_directed_graph, compare_classic_and_directed,
    prepare_effective_directed_source_graph, safe_run_directed_calculation,
)
from src.signaling.storage import load_signaling_cache
from src.services.engine_presentation import (
    EnginePresentationResult, build_directed_presentation, directed_export_frame,
)


log = logging.getLogger(__name__)


@dataclass(slots=True)
class DirectedRunBundle:
    calculation: DirectedCalculationResult
    comparison: EngineComparisonResult | None = None
    presentation: EnginePresentationResult | None = None
    comparison_error: str | None = None


_HYBRID_GRAPH_CACHE: OrderedDict[tuple, object] = OrderedDict()
_HYBRID_GRAPH_CACHE_LOCK = Lock()
_HYBRID_GRAPH_CACHE_MAXSIZE = 1


@lru_cache(maxsize=2)
def _cached_signaling_dataset(path: str, size: int, mtime_ns: int):
    # SignalingDataset is a frozen tuple-based contract and is safe to share.
    return load_signaling_cache(path)


def _load_immutable_signaling_dataset(path: Path):
    stat = path.stat()
    return _cached_signaling_dataset(str(path.resolve()), stat.st_size, stat.st_mtime_ns)


def _cached_hybrid_graph(
    *, graph, dataset, tissue: str, string_dataset_version: str,
    effective_graph_fingerprint: str, taxon_id: int,
    edge_evidence_weighting_mode: str,
):
    key = (
        int(taxon_id), tissue, string_dataset_version, effective_graph_fingerprint,
        dataset.fingerprint, "omnipath_direction_when_available", "bidirectional_equal_weight",
        edge_evidence_weighting_mode,
    )
    with _HYBRID_GRAPH_CACHE_LOCK:
        cached = _HYBRID_GRAPH_CACHE.get(key)
        if cached is not None:
            _HYBRID_GRAPH_CACHE.move_to_end(key)
            return cached
    built = build_hybrid_directed_graph(
        graph, dataset, tissue_context=tissue,
        string_dataset_version=string_dataset_version,
        effective_graph_fingerprint=effective_graph_fingerprint,
        taxon_id=taxon_id,
        edge_evidence_weighting_mode=edge_evidence_weighting_mode,
    )
    with _HYBRID_GRAPH_CACHE_LOCK:
        _HYBRID_GRAPH_CACHE[key] = built
        _HYBRID_GRAPH_CACHE.move_to_end(key)
        while len(_HYBRID_GRAPH_CACHE) > _HYBRID_GRAPH_CACHE_MAXSIZE:
            _HYBRID_GRAPH_CACHE.popitem(last=False)
    return built


def run_optional_directed_engine(
    *, graph, targets: list[str] | tuple[str, ...], tissue: str,
    classic_report: pd.DataFrame | None, project_root: str | Path,
    mode: str, block_weight_fraction: float, damping: float,
    bc_sample_sources: int | None = 4, taxon_id: int = 9606,
    string_dataset_version: str = "local_string",
    scores: pd.DataFrame | None = None,
    gene_to_symbol: dict[str, str] | None = None,
    candidate_limit: int = 20,
    edge_evidence_weighting_mode: str = "legacy_edge_modifiers",
) -> DirectedRunBundle:
    """Run directed calculation without allowing failure to affect classic output."""
    if mode == "Classic":
        return DirectedRunBundle(DirectedCalculationResult(
            status=DirectedEngineStatus.UNAVAILABLE,
            metadata={"engine": "classic", "reason": "directed engine not selected"},
        ))
    if taxon_id != 9606:
        return DirectedRunBundle(DirectedCalculationResult(
            status=DirectedEngineStatus.UNAVAILABLE,
            metadata={"engine": "directed", "species_status": "experimental"},
            error="Mouse directed calculation is unavailable until exact mapping reaches production threshold",
        ))
    try:
        log.info("[DIRECTED][1/7] OmniPath signaling cache yükleniyor")
        root = Path(project_root)
        dataset = _load_immutable_signaling_dataset(
            root / "data" / "processed" / "omnipath_signaling_9606.sqlite"
        )
        log.info("[DIRECTED][2/7] Classic effective biological graph copy hazırlanıyor")
        source = prepare_effective_directed_source_graph(
            graph, targets, edge_evidence_weighting_mode=edge_evidence_weighting_mode,
        )
        log.info("[DIRECTED][3/7] Hybrid directed STRING grafı hazırlanıyor (tissue=%s)", tissue)
        hybrid = _cached_hybrid_graph(
            graph=source.graph, dataset=dataset, tissue=tissue,
            string_dataset_version=string_dataset_version,
            effective_graph_fingerprint=source.effective_fingerprint,
            taxon_id=taxon_id,
            edge_evidence_weighting_mode=edge_evidence_weighting_mode,
        )
        log.info("[DIRECTED][4/7] Directed PageRank ve incident-edge perturbation çalışıyor")
        calculation = safe_run_directed_calculation(
            hybrid, targets=targets, block_weight_fraction=block_weight_fraction,
            damping=damping, bc_sample_sources=bc_sample_sources,
            source_graph_metadata={
                "classic_effective_graph_fingerprint": source.effective_fingerprint,
                "classic_source_fingerprint_before": source.original_fingerprint_before,
                "classic_source_fingerprint_after": source.original_fingerprint_after,
                "classic_source_graph_unchanged": (
                    source.original_fingerprint_before == source.original_fingerprint_after
                ),
                "edge_evidence_weighting_mode": source.modifier_mode,
            },
        )
        presentation = None
        if (
            calculation.status is DirectedEngineStatus.AVAILABLE
            and isinstance(scores, pd.DataFrame) and gene_to_symbol is not None
        ):
            log.info("[DIRECTED][5/7] Canonical Directed presentation sonucu oluşturuluyor")
            presentation = build_directed_presentation(
                calculation, targets=targets, scores=scores,
                gene_to_symbol=gene_to_symbol, candidate_limit=candidate_limit,
            )
        comparison = None
        comparison_error = None
        if (
            mode == "Compare" and calculation.status is DirectedEngineStatus.AVAILABLE
            and isinstance(classic_report, pd.DataFrame) and not classic_report.empty
        ):
            log.info("[DIRECTED][6/7] Classic vs Directed karşılaştırması hazırlanıyor")
            try:
                comparison = compare_classic_and_directed(classic_report, calculation)
            except ValueError as error:
                comparison_error = str(error)
                log.warning("Directed comparison safely unavailable: %s", comparison_error)
        log.info(
            "[DIRECTED][7/7] Directed pipeline hazır: status=%s candidates=%s",
            calculation.status.value,
            len(presentation.candidates) if presentation is not None else "not-built",
        )
        return DirectedRunBundle(calculation, comparison, presentation, comparison_error)
    except Exception as error:
        return DirectedRunBundle(DirectedCalculationResult(
            status=DirectedEngineStatus.UNAVAILABLE,
            metadata={"engine": "directed", "classic_engine_affected": False},
            error=f"{type(error).__name__}: {error}",
        ))


def directed_export_tables(
    bundle: DirectedRunBundle, *, gene_to_symbol: dict[str, str] | None = None,
) -> dict[str, pd.DataFrame]:
    if bundle.presentation is not None:
        directed_table = bundle.presentation.export.copy(deep=True)
    elif gene_to_symbol is not None:
        directed_table = directed_export_frame(bundle.calculation.report, gene_to_symbol=gene_to_symbol)
    else:
        directed_table = bundle.calculation.report.copy()
    tables = {"Directed Redistribution": directed_table}
    if bundle.comparison is not None:
        comparison = bundle.comparison.table.copy(deep=True)
        if gene_to_symbol is not None and not comparison.empty:
            comparison["protein_id"] = comparison["entity"].astype(str)
            comparison["gene"] = comparison["protein_id"].map(gene_to_symbol).fillna("Unresolved")
            comparison["Symbol_Mapping_Status"] = comparison["protein_id"].map(gene_to_symbol).notna().map(
                {True: "EXACT_LOCAL", False: "UNRESOLVED"}
            )
            first = ["gene", "protein_id", "Symbol_Mapping_Status"]
            comparison = comparison[[*first, *[column for column in comparison if column not in {*first, "entity"}]]]
        tables["Classic vs Directed"] = comparison
        tables["Directed Validation"] = pd.DataFrame([
            {key: value for key, value in bundle.comparison.validation.items() if not isinstance(value, tuple)}
        ])
    provenance = {
        **dict(bundle.calculation.metadata), **dict(bundle.calculation.timings),
        "status": bundle.calculation.status.value, "error": bundle.calculation.error,
        "comparison_status": "AVAILABLE" if bundle.comparison is not None else (
            "COMPARISON_INCOMPLETE" if bundle.comparison_error else "NOT_REQUESTED"
        ),
        "comparison_error": bundle.comparison_error,
    }
    tables["Directed Engine Provenance"] = pd.DataFrame([
        {key: str(value) if isinstance(value, (dict, tuple, list)) else value for key, value in provenance.items()}
    ])
    return tables
