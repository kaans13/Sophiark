"""Canonical, engine-specific result flow for UI, interpretation and research."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
from types import MappingProxyType
from typing import Mapping, Sequence

import numpy as np
import pandas as pd

from src.directed.models import DirectedCalculationResult, DirectedEngineStatus


PRIMARY_DAMAGE = "Birincil_Hasar_Hedef"
REDISTRIBUTION_DAMAGE = "Üçüncül_Hasar_Stresli"
UNRESOLVED_SYMBOL = "Unresolved"

_DIRECTED_RENAME = {
    "Directed_PageRank_Baseline": "PageRank_Baseline",
    "Directed_PageRank_Perturbed": "PageRank_Perturbed",
    "Directed_Delta_PageRank": "Delta_PageRank",
    "Directed_Redistribution_Pct": "Delta_PageRank_Pct",
    "Directed_Abs_Delta_PageRank": "Abs_Delta_PageRank",
    "Directed_Response_Direction": "Response_Direction",
    "Directed_BC": "BC_Skoru",
}
_ANNOTATION_COLUMNS = (
    "Lokalizasyon", "Protein_Adi", "GO Biyolojik Süreç", "GO Moleküler İşlev",
    "GO Hücresel Bileşen", "GO_CC_Terimleri", "GO_MF_Terimleri",
)


@dataclass(frozen=True, slots=True)
class EnginePresentationResult:
    engine: str
    report: pd.DataFrame
    candidates: pd.DataFrame
    losses: pd.DataFrame
    full_response: pd.DataFrame
    export: pd.DataFrame
    candidate_fingerprint: str
    metadata: Mapping[str, object]


@dataclass(frozen=True, slots=True)
class ActiveResultFlow:
    engine: str
    report: pd.DataFrame
    signed_response: pd.DataFrame | None
    presentation: EnginePresentationResult | None = None


DIRECTED_UI_COLUMNS = (
    "Gen", "ENSP", "Directed % değişim", "Direction", "Rank", "Directed BC",
    "Relation", "Support", "Coverage", "Localization",
)
COMPARE_UI_COLUMNS = (
    "Gen", "ENSP", "Classic % değişim", "Directed % değişim", "Fark",
    "Classic Rank", "Directed Rank", "Rank Shift", "Directed BC",
    "Relation", "Support", "Coverage", "Localization",
)


def base_engine_for_mode(mode: str) -> str:
    """Return the scientific motor that produces the base report for a UI mode.

    Directed and Compare deliberately start from the Classic report. Evidence
    has its own drop-in motor and must retain a distinct provenance identity.
    """
    return "evidence" if str(mode).strip().casefold() == "evidence" else "classic"


def select_active_result_flow(
    *, classic_report: pd.DataFrame, classic_signed_response: pd.DataFrame | None,
    requested_mode: str, directed_bundle: object | None,
) -> ActiveResultFlow:
    """Select one authoritative downstream result without merging engines."""
    presentation = getattr(directed_bundle, "presentation", None)
    if requested_mode == "Directed" and isinstance(presentation, EnginePresentationResult):
        return ActiveResultFlow(
            engine="directed", report=presentation.report.copy(deep=True),
            signed_response=presentation.full_response.copy(deep=True),
            presentation=presentation,
        )
    if requested_mode == "Evidence":
        return ActiveResultFlow(
            engine="evidence", report=classic_report.copy(deep=True),
            signed_response=(classic_signed_response.copy(deep=True) if isinstance(classic_signed_response, pd.DataFrame) else None),
        )
    return ActiveResultFlow(
        engine="classic", report=classic_report.copy(deep=True),
        signed_response=(
            classic_signed_response.copy(deep=True)
            if isinstance(classic_signed_response, pd.DataFrame) else None
        ),
    )


def directed_ui_frame(presentation: EnginePresentationResult) -> pd.DataFrame:
    """Stable, user-facing Directed result projection in the requested order."""
    frame = presentation.full_response.copy(deep=True)
    projected = pd.DataFrame({
        "Gen": frame.get("Symbol", pd.Series("—", index=frame.index)).fillna("—"),
        "ENSP": frame["protein_id"],
        "Directed % değişim": pd.to_numeric(frame["Delta_PageRank_Pct"], errors="coerce"),
        "Direction": frame.get("Response_Direction", "—"),
        "Rank": frame.get("Directed_Rank", pd.Series(pd.NA, index=frame.index, dtype="Int64")),
        "Directed BC": pd.to_numeric(frame.get("BC_Skoru"), errors="coerce"),
        "Relation": frame.get("Direction_Relation", "—"),
        "Support": frame.get("Direction_Support_Type", "—"),
        "Coverage": frame.get("Direction_Source_Coverage", "—"),
        "Localization": frame.get("Lokalizasyon", "Unknown").fillna("Unknown")
        if isinstance(frame.get("Lokalizasyon"), pd.Series) else "Unknown",
    })
    return projected.loc[:, DIRECTED_UI_COLUMNS].sort_values("Rank", kind="mergesort").reset_index(drop=True)


def compare_ui_frame(comparison_table: pd.DataFrame, presentation: EnginePresentationResult) -> pd.DataFrame:
    """Join comparison metrics to Directed aliases/localization without blending results."""
    comparison = comparison_table.copy(deep=True)
    source = presentation.full_response.copy(deep=True)
    if "Lokalizasyon" not in source:
        source["Lokalizasyon"] = "Unknown"
    context = source[["protein_id", "Symbol", "Lokalizasyon"]].copy()
    context = context.drop_duplicates("protein_id", keep="first")
    comparison = comparison.merge(
        context, left_on="entity", right_on="protein_id", how="left", validate="one_to_one",
    )
    projected = pd.DataFrame({
        "Gen": comparison["Symbol"].fillna("—"),
        "ENSP": comparison["protein_id"],
        "Classic % değişim": comparison["Classic_Redistribution_Pct"],
        "Directed % değişim": comparison["Directed_Redistribution_Pct"],
        "Fark": comparison["Response_Delta"],
        "Classic Rank": comparison["Classic_Rank"],
        "Directed Rank": comparison["Directed_Rank"],
        "Rank Shift": comparison["Rank_Shift"],
        "Directed BC": comparison["Directed_BC"],
        "Relation": comparison["Direction_Relation"],
        "Support": comparison["Direction_Support_Type"],
        "Coverage": comparison["Direction_Source_Coverage"],
        "Localization": comparison["Lokalizasyon"].fillna("Unknown"),
    })
    return projected.loc[:, COMPARE_UI_COLUMNS].sort_values("Directed Rank", kind="mergesort").reset_index(drop=True)


def _symbol_for(protein_id: str, mapping: Mapping[str, str]) -> tuple[str | None, str]:
    value = mapping.get(str(protein_id))
    if value is None or pd.isna(value) or not str(value).strip():
        return None, "UNRESOLVED"
    return str(value).strip(), "EXACT_LOCAL"


def candidate_set_fingerprint(engine: str, protein_ids: Sequence[str]) -> str:
    """Engine-sensitive identity for enrichment/interpretation candidate inputs."""
    payload = {
        "engine": str(engine).strip().casefold(),
        "protein_ids": sorted({str(value) for value in protein_ids if str(value).strip()}),
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _fingerprint(engine: str, candidates: pd.DataFrame) -> str:
    return candidate_set_fingerprint(
        engine, candidates.get("protein_id", pd.Series(dtype=str)).dropna().astype(str).tolist(),
    )


def _attach_symbols(frame: pd.DataFrame, mapping: Mapping[str, str]) -> pd.DataFrame:
    result = frame.copy(deep=True)
    resolved = [
        _symbol_for(protein_id, mapping)
        for protein_id in result["protein_id"].astype(str)
    ]
    result["Symbol"] = [item[0] for item in resolved]
    result["gene_symbol"] = [item[0] or UNRESOLVED_SYMBOL for item in resolved]
    result["Symbol_Mapping_Status"] = [item[1] for item in resolved]
    return result


def _annotation_frame(scores: pd.DataFrame) -> pd.DataFrame:
    if not isinstance(scores, pd.DataFrame) or "gene" not in scores:
        return pd.DataFrame(columns=["protein_id"])
    columns = ["gene", *[column for column in _ANNOTATION_COLUMNS if column in scores]]
    result = scores[columns].drop_duplicates("gene", keep="first").copy(deep=True)
    return result.rename(columns={"gene": "protein_id"})


def _canonical_responses(
    calculation: DirectedCalculationResult,
    *, scores: pd.DataFrame, gene_to_symbol: Mapping[str, str],
) -> pd.DataFrame:
    frame = calculation.report.copy(deep=True)
    if frame.empty:
        return frame
    frame["protein_id"] = frame["entity"].astype(str)
    frame["gene"] = frame["protein_id"]  # authoritative internal network identifier
    for source, target in _DIRECTED_RENAME.items():
        frame[target] = frame[source]
    frame["Directed_Rank"] = pd.to_numeric(
        frame["Delta_PageRank_Pct"], errors="coerce",
    ).abs().rank(method="min", ascending=False).astype("Int64")
    frame["Hasar_Tipi"] = REDISTRIBUTION_DAMAGE
    frame["Engine"] = "directed"
    frame["Graph_Mode"] = calculation.metadata.get("graph_mode", "hybrid_string_omnipath")
    frame["Perturbation_Strategy"] = calculation.metadata.get("directed_perturbation_strategy")
    frame["Suppression_Factor"] = calculation.metadata.get("suppression_factor")
    frame["Direction_Dataset_Version"] = calculation.metadata.get("direction_dataset_version")
    frame["Direction_Coverage_Global"] = calculation.metadata.get("direction_coverage")
    frame["Directed_Significance_Status"] = "unavailable"
    frame["Redistribution_Selection_Mode"] = "top_n"
    annotations = _annotation_frame(scores)
    if not annotations.empty:
        frame = frame.merge(annotations, on="protein_id", how="left", validate="one_to_one")
    return _attach_symbols(frame, gene_to_symbol)


def _target_rows(
    calculation: DirectedCalculationResult,
    *, targets: Sequence[str], scores: pd.DataFrame, gene_to_symbol: Mapping[str, str],
) -> pd.DataFrame:
    annotations = _annotation_frame(scores).set_index("protein_id", drop=False)
    target_pr = dict(calculation.metadata.get("target_pagerank", {}) or {})
    rows: list[dict[str, object]] = []
    for target in dict.fromkeys(map(str, targets)):
        values = dict(target_pr.get(target, {}) or {})
        before = float(values.get("baseline", np.nan))
        after = float(values.get("perturbed", np.nan))
        delta = after - before
        row: dict[str, object] = {
            "entity": target, "protein_id": target, "gene": target,
            "PageRank_Baseline": before, "PageRank_Perturbed": after,
            "Delta_PageRank": delta,
            "Delta_PageRank_Pct": (delta / before * 100.0 if np.isfinite(before) and before > 0 else np.nan),
            "Abs_Delta_PageRank": abs(delta), "Response_Direction": "Perturbation Target",
            "Hasar_Tipi": PRIMARY_DAMAGE, "Engine": "directed",
            "Graph_Mode": calculation.metadata.get("graph_mode", "hybrid_string_omnipath"),
            "Perturbation_Strategy": calculation.metadata.get("directed_perturbation_strategy"),
            "Suppression_Factor": calculation.metadata.get("suppression_factor"),
            "Direction_Dataset_Version": calculation.metadata.get("direction_dataset_version"),
            "Direction_Coverage_Global": calculation.metadata.get("direction_coverage"),
            "Directed_Significance_Status": "unavailable",
            "Redistribution_Selection_Mode": "top_n",
        }
        if target in annotations.index:
            for column in annotations.columns:
                if column != "protein_id":
                    row[column] = annotations.at[target, column]
        rows.append(row)
    return _attach_symbols(pd.DataFrame(rows), gene_to_symbol) if rows else pd.DataFrame()


def directed_export_frame(
    frame: pd.DataFrame, *, gene_to_symbol: Mapping[str, str],
) -> pd.DataFrame:
    """Return a user-facing export where ``gene`` is never a disguised ENSP ID."""
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return pd.DataFrame()
    result = frame.copy(deep=True)
    if "protein_id" not in result:
        source = "entity" if "entity" in result else "gene"
        result["protein_id"] = result[source].astype(str)
    result = _attach_symbols(result.drop(columns=[c for c in ("Symbol", "gene_symbol", "Symbol_Mapping_Status") if c in result]), gene_to_symbol)
    result["gene"] = result["gene_symbol"]
    preferred = [
        "gene", "protein_id", "Symbol_Mapping_Status",
        "Directed_PageRank_Baseline", "Directed_PageRank_Perturbed",
        "Directed_Delta_PageRank", "Directed_Redistribution_Pct", "Directed_BC",
        "Direction_Relation", "Direction_Support_Type", "Directed_Hop_Distance",
        "Directed_Distance",
        "Sign_Context", "Direction_Source_Coverage", "Alternative_Route_Status",
        "Engine", "Graph_Mode", "Perturbation_Strategy", "Suppression_Factor",
        "Direction_Dataset_Version", "Direction_Coverage_Global",
        "Directed_Significance_Status",
    ]
    return result[[*[column for column in preferred if column in result], *[column for column in result if column not in preferred and column not in {"entity", "gene_symbol", "Symbol"}]]]


def build_directed_presentation(
    calculation: DirectedCalculationResult,
    *, targets: Sequence[str], scores: pd.DataFrame,
    gene_to_symbol: Mapping[str, str], candidate_limit: int,
    minimum_abs_response_pct: float = 0.05,
) -> EnginePresentationResult:
    if calculation.status is not DirectedEngineStatus.AVAILABLE:
        raise ValueError("Directed presentation requires an AVAILABLE calculation")
    full = _canonical_responses(calculation, scores=scores, gene_to_symbol=gene_to_symbol)
    delta = pd.to_numeric(full.get("Delta_PageRank_Pct"), errors="coerce")
    candidates = full.loc[delta > abs(float(minimum_abs_response_pct))].sort_values(
        "Delta_PageRank_Pct", ascending=False, kind="mergesort",
    ).head(max(0, int(candidate_limit))).copy()
    losses = full.loc[delta < -abs(float(minimum_abs_response_pct))].sort_values(
        "Delta_PageRank_Pct", ascending=True, kind="mergesort",
    ).copy()
    targets_frame = _target_rows(
        calculation, targets=targets, scores=scores, gene_to_symbol=gene_to_symbol,
    )
    report = pd.concat([targets_frame, candidates], ignore_index=True, sort=False)
    export = directed_export_frame(full, gene_to_symbol=gene_to_symbol)
    metadata = MappingProxyType({
        **dict(calculation.metadata),
        "engine": "directed",
        "candidate_count": len(candidates),
        "positive_response_count": int((delta > abs(float(minimum_abs_response_pct))).sum()),
        "negative_response_count": int((delta < -abs(float(minimum_abs_response_pct))).sum()),
        "affected_record_count": len(targets_frame) + len(candidates),
        "strongest_gain_pct": float(delta.max()) if delta.notna().any() else None,
        "strongest_loss_pct": float(delta.min()) if delta.notna().any() else None,
        "candidate_limit": int(candidate_limit),
        "minimum_abs_response_pct": float(minimum_abs_response_pct),
    })
    return EnginePresentationResult(
        engine="directed", report=report, candidates=candidates, losses=losses,
        full_response=full, export=export,
        candidate_fingerprint=_fingerprint("directed", candidates), metadata=metadata,
    )
