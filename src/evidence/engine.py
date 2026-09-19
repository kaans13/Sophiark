from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from types import MappingProxyType
from typing import Any, Mapping
import tempfile

import numpy as np
import pandas as pd

from src import biology_logic as classic_biology
from src.config import TF_TARGETS, _ensp_to_symbol
from src.graph_engine import compute_pagerank, run_community_detection, tag_graph_nodes
from src.metrics import build_hinterland_score, classify_genes, compute_druggable_gate_score, detect_customs_gates
from src.data_loader import fetch_degree_summary
from src.config import OUT_CSV, PAGERANK_DAMPING, PAGERANK_MAX_ITER

from .config import CLASSIC_DB_PATH, EvidenceConfig, EvidenceMode
from .complex_response import target_complex_response
from .graph import build_evidence_graph
from .provenance import attach_provenance, build_run_provenance
from .provenance import provenance_frame
from .audit import audit_modifier_bias
from .comparison import compare_classic_evidence


@dataclass(slots=True)
class EvidenceRunResult:
    report: pd.DataFrame
    signed_response: pd.DataFrame
    complexes: pd.DataFrame
    edge_evidence: pd.DataFrame
    provenance: Mapping[str, Any]
    comparison: pd.DataFrame = field(default_factory=pd.DataFrame)


_LATEST_RESULT: EvidenceRunResult | None = None


def latest_result() -> EvidenceRunResult | None:
    return _LATEST_RESULT


def _score_graph(graph, bc_sample_sources: int | None) -> pd.DataFrame:
    annotations = pd.read_csv(OUT_CSV)
    amap = annotations[[c for c in ("gene", "Lokalizasyon") if c in annotations]].drop_duplicates("gene").set_index("gene")
    classic_degree = fetch_degree_summary(str(CLASSIC_DB_PATH), 500).set_index("gene")
    degree = pd.DataFrame({"gene": graph.vs["name"], "k_i": graph.degree()})
    degree["max_cs"] = degree["gene"].map(classic_degree.get("max_cs", pd.Series(dtype=float))).fillna(0).astype(int)
    degree["Classic_Degree"] = degree["gene"].map(classic_degree.get("k_i", pd.Series(dtype=float))).fillna(0).astype(int)
    pr = compute_pagerank(graph, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
    if "classic_weight" in graph.es.attribute_names():
        classic_graph = graph.copy()
        classic_graph.es["weight"] = classic_graph.es["classic_weight"]
        classic_pr = compute_pagerank(classic_graph, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
    else:
        classic_pr = pr
    degree["Classic_PageRank"] = degree["gene"].map(classic_pr).fillna(0.0)
    degree = build_hinterland_score(degree, pr)
    degree = classify_genes(degree)
    degree, _ = run_community_detection(graph, degree)
    degree["Lokalizasyon"] = degree["gene"].map(amap.get("Lokalizasyon", pd.Series(dtype=str))).fillna("Unknown")
    tag_graph_nodes(graph, dict(zip(degree["gene"], degree["Lokalizasyon"])))
    degree = detect_customs_gates(graph, degree, bc_sample_sources=bc_sample_sources)
    return compute_druggable_gate_score(degree)


def prepare_evidence_network(*, config: EvidenceConfig, forced_genes: list[str] | None,
                             tissue: str | None, bc_sample_sources: int | None = None):
    graph, profiles = build_evidence_graph(config=config, forced_genes=forced_genes, tissue=tissue)
    modifier_audit: dict[str, dict[str, object]] = {}
    if config.mode is not EvidenceMode.PARITY:
        edge_context = _edge_frame(graph)
        node_context = _node_evidence(edge_context)
        classic_graph = graph.copy()
        classic_graph.es["weight"] = classic_graph.es["classic_weight"]
        classic_pr = compute_pagerank(classic_graph, PAGERANK_DAMPING, PAGERANK_MAX_ITER)
        controls = pd.DataFrame({
            "gene": graph.vs["name"], "Classic_Degree": graph.degree(),
            "Classic_PageRank": [classic_pr[name] for name in graph.vs["name"]],
        }).merge(node_context, on="gene", how="left")
        for column, edge_attribute in (
            ("Applied_Text_Modifier", "text_bonus"),
            ("Physical_Coverage", "physical_bonus"),
        ):
            decision = audit_modifier_bias(controls, modifier_column=column)
            modifier_audit[column] = {
                "usable_as_multiplier": decision.usable_as_multiplier,
                "reason": decision.reason,
                "spearman_classic_degree": decision.spearman_classic_degree,
                "spearman_classic_pagerank": decision.spearman_classic_pagerank,
                "spearman_annotation_availability": decision.spearman_annotation_availability,
            }
            configured = (column == "Applied_Text_Modifier" and config.lambda_text > 0) or (
                column == "Physical_Coverage" and config.lambda_physical > 0
            )
            if configured and not decision.usable_as_multiplier:
                values = np.asarray(graph.es[edge_attribute], dtype=float)
                graph.es["weight"] = (np.asarray(graph.es["weight"], dtype=float) / values).tolist()
                graph.es["distance"] = (1.0 / np.asarray(graph.es["weight"], dtype=float)).tolist()
                graph.es[edge_attribute] = [1.0] * graph.ecount()
        graph["evidence_modifier_audit"] = modifier_audit
    scores = _score_graph(graph, bc_sample_sources)
    graph["evidence_config"] = config.to_dict()
    graph["evidence_tissue"] = tissue or "ALL"
    return graph, scores, profiles


def _apply_tf_only(targets: list[str], graph, factor: float = .20) -> dict[int, tuple[float, float]]:
    snapshot: dict[int, tuple[float, float]] = {}
    names = set(graph.vs["name"])
    for target in targets:
        if target not in names:
            continue
        symbol = _ensp_to_symbol(target)
        known = set(TF_TARGETS.get(symbol, ())) if symbol else set()
        if not known:
            continue
        tidx = graph.vs.find(name=target).index
        for nidx in graph.neighbors(tidx):
            if _ensp_to_symbol(graph.vs[nidx]["name"]) not in known:
                continue
            eid = graph.get_eid(tidx, nidx)
            snapshot[eid] = (graph.es[eid]["weight"], graph.es[eid]["distance"])
            graph.es[eid]["weight"] *= 1.0 + factor
            graph.es[eid]["distance"] = 1.0 / graph.es[eid]["weight"]
    return snapshot


def _restore(graph, snapshot: dict[int, tuple[float, float]]) -> None:
    for eid, (weight, distance) in snapshot.items():
        graph.es[eid]["weight"] = weight
        graph.es[eid]["distance"] = distance


def _edge_frame(graph) -> pd.DataFrame:
    rows = []
    attrs = set(graph.es.attribute_names())
    wanted = [c for c in (
        "combined_score", "reconstructed_full", "s_nontext", "textmining",
        "textmining_transferred", "text_dependency", "proposed_text_bonus", "text_bonus",
        "physical_score", "physical_status", "proposed_physical_bonus", "physical_bonus",
        "experiments", "experiments_transferred", "database", "database_transferred",
        "coexpression", "coexpression_transferred", "neighborhood",
        "neighborhood_transferred", "fusion", "cooccurence", "homology",
    ) if c in attrs]
    for edge in graph.es:
        row = {"protein1": graph.vs[edge.source]["name"], "protein2": graph.vs[edge.target]["name"]}
        row.update({column: edge[column] for column in wanted})
        rows.append(row)
    result = pd.DataFrame(rows)
    channel_columns = [c for c in (
        "neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
        "coexpression", "coexpression_transferred", "experiments",
        "experiments_transferred", "database", "database_transferred",
        "textmining", "textmining_transferred",
    ) if c in result]
    if channel_columns:
        result["Evidence_Channel_Count"] = result[channel_columns].apply(
            pd.to_numeric, errors="coerce"
        ).gt(0).sum(axis=1)
    return result


def _node_evidence(edge_frame: pd.DataFrame) -> pd.DataFrame:
    if edge_frame.empty:
        return pd.DataFrame(columns=["gene"])
    left = edge_frame.rename(columns={"protein1": "gene"})
    right = edge_frame.rename(columns={"protein2": "gene"})
    both = pd.concat([left, right], ignore_index=True)
    result = both.groupby("gene", as_index=False).agg(
        Non_Text_Evidence=("s_nontext", "median"),
        Text_Dependency=("text_dependency", "median"),
        Proposed_Text_Modifier=("proposed_text_bonus", "mean"),
        Applied_Text_Modifier=("text_bonus", "mean"),
        Applied_Physical_Modifier=("physical_bonus", "mean"),
        Physical_Supported_Edges=("physical_status", lambda s: int((s == "SUPPORTED").sum())),
        Physical_Available_Edges=("physical_status", lambda s: int(s.isin(["SUPPORTED", "NOT_SUPPORTED"]).sum())),
        Annotation_Availability=("Evidence_Channel_Count", "median"),
    )
    result["Physical_Coverage"] = (
        result["Physical_Supported_Edges"]
        / result["Physical_Available_Edges"].replace(0, np.nan)
    )
    return result


def run_evidence_simulation(*, graph, scores: pd.DataFrame, out_dir: Path,
                            targets: list[str], config: EvidenceConfig,
                            tissue: str | None, block_weight_fraction: float,
                            damping: float, scientific_options: dict | None = None) -> EvidenceRunResult:
    global _LATEST_RESULT
    options = dict(scientific_options or {})
    classic_signed = pd.DataFrame()
    if config.mode is not EvidenceMode.E2_TOPOLOGY_WEIGHT and "classic_weight" in graph.es.attribute_names():
        classic_graph = graph.copy()
        classic_graph.es["weight"] = classic_graph.es["classic_weight"]
        classic_graph.es["distance"] = (1.0 / np.asarray(classic_graph.es["weight"], dtype=float)).tolist()
        with tempfile.TemporaryDirectory(prefix="sophiark_evidence_classic_control_") as control_dir:
            classic_biology.run_infection_simulation(
                classic_graph, scores, Path(control_dir), spesifik_hedefler=targets,
                block_weight_fraction=block_weight_fraction, damping=damping,
                edge_evidence_weighting_mode="legacy_edge_modifiers", tissue=tissue,
                **{k: v for k, v in options.items() if k not in {"edge_evidence_weighting_mode", "tissue"}},
            )
            latest_classic = classic_biology.latest_signed_redistribution()
            classic_signed = latest_classic.copy(deep=True) if isinstance(latest_classic, pd.DataFrame) else pd.DataFrame()
    snapshot: dict[int, tuple[float, float]] = {}
    # Parity deliberately invokes the complete unchanged Classic modifier path.
    modifier_mode = "legacy_edge_modifiers" if config.mode is EvidenceMode.PARITY else "structural_only"
    if config.mode is not EvidenceMode.PARITY:
        snapshot = _apply_tf_only(targets, graph)
        classic_biology.apply_functional_penalty(targets, graph, snapshot=snapshot)
    try:
        report = classic_biology.run_infection_simulation(
            graph, scores, out_dir, spesifik_hedefler=targets,
            block_weight_fraction=block_weight_fraction, damping=damping,
            edge_evidence_weighting_mode=modifier_mode, tissue=tissue,
            **{k: v for k, v in options.items() if k not in {"edge_evidence_weighting_mode", "tissue"}},
        )
        signed = classic_biology.latest_signed_redistribution()
        signed = signed.copy(deep=True) if isinstance(signed, pd.DataFrame) else pd.DataFrame()
    finally:
        _restore(graph, snapshot)
    edges = _edge_frame(graph)
    node_context = _node_evidence(edges)
    provenance = build_run_provenance(
        config=config, tissue=tissue, targets=targets,
        perturbation={
            "edge_attenuation_fraction": block_weight_fraction, "pagerank_damping": damping,
            "modifier_bias_audit": graph["evidence_modifier_audit"] if "evidence_modifier_audit" in graph.attributes() else {},
            "actual_text_multiplier_enabled": bool(any(abs(float(v)-1.0) > 1e-12 for v in graph.es["text_bonus"])),
            "actual_physical_multiplier_enabled": bool(any(abs(float(v)-1.0) > 1e-12 for v in graph.es["physical_bonus"])),
            **options,
        },
    )
    report = attach_provenance(report.merge(node_context, on="gene", how="left"), provenance)
    signed = attach_provenance(signed.merge(node_context, on="gene", how="left"), provenance)
    symbols = pd.read_csv(Path(__file__).resolve().parents[2] / "data" / "processed" / "ensp_with_symbols.csv")
    ensp_to_symbol = dict(zip(symbols["gene"].astype(str), symbols["Symbol"].astype(str)))
    complexes = target_complex_response([ensp_to_symbol.get(t, t) for t in targets], signed, ensp_to_symbol)
    complexes = attach_provenance(complexes, provenance)
    edges = attach_provenance(edges, provenance)
    comparison = (
        attach_provenance(compare_classic_evidence(classic_signed, signed), provenance)
        if not classic_signed.empty else pd.DataFrame()
    )
    reports_dir = Path(out_dir) / "reports"
    reports_dir.mkdir(parents=True, exist_ok=True)
    # Replace the compatibility CSV written by the reused perturbation routine:
    # an Evidence artifact may never exist without its actual run provenance.
    report.to_csv(Path(out_dir) / "enfeksiyon_sok_dalgasi_raporu.csv", index=False, encoding="utf-8-sig")
    report.to_csv(reports_dir / "evidence_redistribution.csv", index=False, encoding="utf-8-sig")
    complexes.to_csv(reports_dir / "evidence_complex_response.csv", index=False, encoding="utf-8-sig")
    comparison.to_csv(reports_dir / "classic_evidence_comparison.csv", index=False, encoding="utf-8-sig")
    provenance_frame(provenance).to_csv(
        reports_dir / "evidence_run_provenance.csv", index=False, encoding="utf-8-sig"
    )
    _LATEST_RESULT = EvidenceRunResult(report, signed, complexes, edges, MappingProxyType(provenance), comparison)
    return _LATEST_RESULT
