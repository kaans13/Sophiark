"""Classic + Directed research-report orchestration.

This module deliberately owns no graph construction, edge weighting or
PageRank mathematics.  It calls the existing engines on isolated graph copies,
then adds descriptive synthesis, STRING provenance and CORUM context.
"""
from __future__ import annotations

from collections import OrderedDict
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from enum import Enum
import hashlib
import json
from pathlib import Path
import sqlite3
import tempfile
from threading import Lock
from typing import Any, Callable, Mapping

import numpy as np
import pandas as pd

from src.product.runtime_resources import (
    DetachedPreparedContextCache, file_identity, product_data_snapshot,
    symbol_mapping,
)
from src.product.fingerprints import stable_fingerprint


UNIFIED_ENGINE_VERSION = "unified-research-report-1.0"
EVIDENCE_WEIGHTING_STATUS = "BIAS_WARNING"
_CACHE_MAXSIZE = 8
_CACHE: OrderedDict[str, "UnifiedResearchReport"] = OrderedDict()
_CACHE_LOCK = Lock()
_PREPARED_CONTEXT_CACHE = DetachedPreparedContextCache(maxsize=2)


class UnifiedStatus(str, Enum):
    COMPLETE = "COMPLETE"
    ENGINE_PARTIAL = "ENGINE_PARTIAL"
    FAILED = "FAILED"


@dataclass(frozen=True, slots=True)
class PresentationConvention:
    """Explicit presentation policy, separate from both calculation engines."""

    top_n: int = 50
    meaningful_response_pct: float = 0.05
    direction_sensitivity_percentile: float = 0.25


@dataclass(slots=True)
class UnifiedResearchReport:
    target: str
    tissue: str
    status: UnifiedStatus
    classic_status: str
    directed_status: str
    candidates: pd.DataFrame = field(default_factory=pd.DataFrame)
    classic_raw: pd.DataFrame = field(default_factory=pd.DataFrame)
    directed_raw: pd.DataFrame = field(default_factory=pd.DataFrame)
    complex_context: pd.DataFrame = field(default_factory=pd.DataFrame)
    agreement: Mapping[str, object] = field(default_factory=dict)
    provenance: Mapping[str, object] = field(default_factory=dict)
    context_availability: Mapping[str, object] = field(default_factory=dict)
    errors: Mapping[str, str] = field(default_factory=dict)

    def copy(self) -> "UnifiedResearchReport":
        return UnifiedResearchReport(
            target=self.target, tissue=self.tissue, status=self.status,
            classic_status=self.classic_status, directed_status=self.directed_status,
            candidates=self.candidates.copy(deep=True), classic_raw=self.classic_raw.copy(deep=True),
            directed_raw=self.directed_raw.copy(deep=True), complex_context=self.complex_context.copy(deep=True),
            agreement=dict(self.agreement), provenance=dict(self.provenance),
            context_availability=dict(self.context_availability), errors=dict(self.errors),
        )


def _json_hash(value: Mapping[str, object]) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, default=str, separators=(",", ":")).encode()).hexdigest()


def _safe_identity(path: Path) -> dict[str, object] | None:
    if not path.exists():
        return None
    stat = path.stat()
    return {"path": str(path), "bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}


def _symbol_mapping() -> dict[str, str]:
    path = Path(__file__).resolve().parents[2] / "data" / "processed" / "ensp_with_symbols.csv"
    return symbol_mapping(path)


def _rank_percentile(rank: pd.Series) -> pd.Series:
    result = pd.Series(pd.NA, index=rank.index, dtype="Float64")
    usable = pd.to_numeric(rank, errors="coerce").dropna()
    count = len(usable)
    if count == 1:
        result.loc[usable.index] = 1.0
    elif count > 1:
        result.loc[usable.index] = 1.0 - (usable - 1.0) / float(count - 1)
    return result


def _engine_response_frame(raw: pd.DataFrame, *, engine: str, symbols: Mapping[str, str]) -> pd.DataFrame:
    if not isinstance(raw, pd.DataFrame) or raw.empty:
        return pd.DataFrame({
            "entity_id": pd.Series(dtype="string"), "symbol": pd.Series(dtype="string"),
            f"{engine}_response": pd.Series(dtype=float), f"{engine}_baseline": pd.Series(dtype=float),
            f"{engine}_perturbed": pd.Series(dtype=float),
        })
    entity_col = "gene" if engine == "classic" else "entity"
    response_col = "Delta_PageRank_Pct" if engine == "classic" else "Directed_Redistribution_Pct"
    baseline_col = "PageRank_Baseline" if engine == "classic" else "Directed_PageRank_Baseline"
    perturbed_col = "PageRank_Perturbed" if engine == "classic" else "Directed_PageRank_Perturbed"
    required = {entity_col, response_col}
    missing = required - set(raw.columns)
    if missing:
        raise ValueError(f"{engine} response universe missing required columns: {sorted(missing)}")
    source = raw.copy(deep=True)
    source["entity_id"] = source[entity_col].astype(str)
    if source["entity_id"].duplicated().any():
        raise ValueError(f"{engine} response universe contains duplicate canonical entities")
    response = pd.to_numeric(source[response_col], errors="coerce")
    if response.isna().any():
        raise ValueError(f"{engine} response universe has missing responses; missing is not zero")
    prefix = engine.capitalize()
    result = pd.DataFrame({
        "entity_id": source["entity_id"], "symbol": source["entity_id"].map(symbols),
        f"{engine}_response": response,
        f"{engine}_baseline": pd.to_numeric(source.get(baseline_col), errors="coerce"),
        f"{engine}_perturbed": pd.to_numeric(source.get(perturbed_col), errors="coerce"),
    })
    absolute_rank = response.abs().rank(method="min", ascending=False)
    positive_rank = response.where(response > 0).rank(method="min", ascending=False)
    negative_rank = response.abs().where(response < 0).rank(method="min", ascending=False)
    result[f"{engine}_absolute_response_rank"] = absolute_rank.astype("Int64")
    result[f"{engine}_positive_rank"] = positive_rank.astype("Int64")
    result[f"{engine}_negative_rank"] = negative_rank.astype("Int64")
    result[f"{engine}_absolute_rank_percentile"] = _rank_percentile(absolute_rank)
    result[f"{engine}_positive_rank_percentile"] = _rank_percentile(positive_rank)
    result[f"{engine}_negative_rank_percentile"] = _rank_percentile(negative_rank)
    result[f"{engine}_raw_columns_available"] = ";".join(source.columns)
    return result


def _sign(values: pd.Series, threshold: float) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    return pd.Series(np.select([numeric > threshold, numeric < -threshold], ["POSITIVE", "NEGATIVE"], default="NEUTRAL"), index=values.index)


def _top_set(values: pd.Series, *, positive: bool | None, count: int) -> set[str]:
    numeric = pd.to_numeric(values, errors="coerce")
    if positive is True:
        numeric = numeric.loc[numeric > 0]
        order = numeric.sort_values(ascending=False, kind="mergesort")
    elif positive is False:
        numeric = numeric.loc[numeric < 0].abs()
        order = numeric.sort_values(ascending=False, kind="mergesort")
    else:
        order = numeric.abs().sort_values(ascending=False, kind="mergesort")
    return set(order.head(min(int(count), len(order))).index.astype(str))


def _agreement_metrics(classic: pd.DataFrame, directed: pd.DataFrame, aligned: pd.DataFrame) -> dict[str, object]:
    if classic.empty or directed.empty:
        return {"status": "UNAVAILABLE", "reason": "one calculation engine has no full response universe"}
    joined = aligned.loc[aligned["classic_available"] & aligned["directed_available"]].copy()
    if joined.empty:
        return {"status": "UNAVAILABLE", "reason": "no canonical entity overlap"}
    left = joined["classic_response"].astype(float)
    right = joined["directed_response"].astype(float)
    top = {}
    for label, positive, sizes in (("positive", True, (20, 50)), ("negative", False, (20, 50)), ("absolute", None, (50, 100))):
        for size in sizes:
            a = _top_set(left.set_axis(joined["entity_id"]), positive=positive, count=size)
            b = _top_set(right.set_axis(joined["entity_id"]), positive=positive, count=size)
            top[f"{label}_top{size}_overlap"] = len(a & b)
    same = joined["same_sign"].fillna(False).astype(bool)
    return {
        "status": "AVAILABLE", "aligned_entities": len(joined),
        "signed_spearman": float(left.corr(right, method="spearman")) if left.nunique() > 1 and right.nunique() > 1 else None,
        "signed_pearson": float(left.corr(right, method="pearson")) if left.nunique() > 1 and right.nunique() > 1 else None,
        "sign_agreement_count": int(same.sum()), "sign_agreement_rate": float(same.mean()),
        **top,
    }


def _synthesize(*, target: str, tissue: str, classic_raw: pd.DataFrame, directed_raw: pd.DataFrame,
                symbols: Mapping[str, str], convention: PresentationConvention) -> tuple[pd.DataFrame, dict[str, object]]:
    classic = _engine_response_frame(classic_raw, engine="classic", symbols=symbols)
    directed = _engine_response_frame(directed_raw, engine="directed", symbols=symbols)
    classic = classic.rename(columns={"symbol": "classic_symbol"})
    directed = directed.rename(columns={"symbol": "directed_symbol"})
    aligned = classic.merge(directed, on="entity_id", how="outer", validate="one_to_one")
    aligned["symbol"] = aligned["classic_symbol"].combine_first(aligned["directed_symbol"])
    aligned["classic_available"] = aligned["classic_response"].notna()
    aligned["directed_available"] = aligned["directed_response"].notna()
    both = aligned["classic_available"] & aligned["directed_available"]
    aligned["classic_sign"] = pd.Series(pd.NA, index=aligned.index, dtype="string")
    aligned["directed_sign"] = pd.Series(pd.NA, index=aligned.index, dtype="string")
    aligned.loc[aligned["classic_available"], "classic_sign"] = _sign(aligned.loc[aligned["classic_available"], "classic_response"], convention.meaningful_response_pct)
    aligned.loc[aligned["directed_available"], "directed_sign"] = _sign(aligned.loc[aligned["directed_available"], "directed_response"], convention.meaningful_response_pct)
    aligned["same_sign"] = both & aligned["classic_sign"].eq(aligned["directed_sign"])
    aligned["sign_flip"] = both & aligned["classic_sign"].isin(["POSITIVE", "NEGATIVE"]) & aligned["directed_sign"].isin(["POSITIVE", "NEGATIVE"]) & ~aligned["same_sign"]
    aligned["percentile_difference"] = aligned["directed_absolute_rank_percentile"] - aligned["classic_absolute_rank_percentile"]
    aligned["direction_sensitivity"] = aligned["percentile_difference"].abs()
    classic_prominent = aligned["classic_absolute_response_rank"].le(convention.top_n)
    directed_prominent = aligned["directed_absolute_response_rank"].le(convention.top_n)
    aligned["robust_cross_model"] = both & aligned["same_sign"] & classic_prominent & directed_prominent
    aligned["direction_sensitive"] = both & (
        aligned["sign_flip"] | aligned["direction_sensitivity"].ge(convention.direction_sensitivity_percentile)
    )
    aligned["classic_dominant"] = both & classic_prominent & ~directed_prominent
    aligned["directed_emergent"] = both & directed_prominent & ~classic_prominent
    aligned["finding_class"] = np.select(
        [aligned["robust_cross_model"], aligned["direction_sensitive"], aligned["classic_dominant"], aligned["directed_emergent"], aligned["classic_available"] & ~aligned["directed_available"], aligned["directed_available"] & ~aligned["classic_available"]],
        ["ROBUST_CROSS_MODEL", "DIRECTION_SENSITIVE", "CLASSIC_DOMINANT", "DIRECTED_EMERGENT", "CLASSIC_ONLY_AVAILABLE", "DIRECTED_ONLY_AVAILABLE"],
        default="LOW_AGREEMENT",
    )
    aligned["target"] = target
    aligned["tissue"] = tissue
    agreement = _agreement_metrics(classic, directed, aligned)
    return aligned.sort_values(["robust_cross_model", "direction_sensitivity", "classic_absolute_rank_percentile", "directed_absolute_rank_percentile"], ascending=[False, False, False, False], kind="mergesort").reset_index(drop=True), agreement


def _target_edge_provenance(candidates: pd.DataFrame, *, target: str) -> tuple[pd.DataFrame, dict[str, object]]:
    """Attach target-edge profiles from the existing local STRING index only."""
    from src.evidence.config import EVIDENCE_DB_PATH
    result = candidates.copy(deep=True)
    if result.empty:
        return result, {"status": "AVAILABLE", "indexed_candidates": 0}
    status = pd.Series("NO_DIRECT_INDEXED_TARGET_EDGE", index=result.index, dtype="string")
    result["Evidence_Provenance_Available"] = False
    result["Evidence_Provenance_Status"] = status
    if not EVIDENCE_DB_PATH.exists():
        result["Evidence_Provenance_Status"] = "INDEX_UNAVAILABLE"
        return result, {"status": "UNAVAILABLE", "reason": "local STRING evidence index missing"}
    columns = ("combined_score", "s_nontext", "text_dependency", "textmining", "textmining_transferred", "experiments", "experiments_transferred", "database", "database_transferred", "coexpression", "coexpression_transferred", "neighborhood", "neighborhood_transferred", "fusion", "cooccurence", "physical_score", "physical_status")
    lookup: dict[str, dict[str, object]] = {}
    clean_target = str(target).replace("9606.", "")
    con = sqlite3.connect(EVIDENCE_DB_PATH)
    try:
        for entity in result["entity_id"].astype(str):
            clean_entity = entity.replace("9606.", "")
            if clean_entity == clean_target:
                continue
            first, second = sorted((clean_target, clean_entity))
            row = con.execute(f"SELECT {','.join(columns)} FROM evidence WHERE protein1=? AND protein2=?", (first, second)).fetchone()
            if row is not None:
                lookup[entity] = dict(zip(columns, row))
    finally:
        con.close()
    for source in columns:
        destination = "cooccurrence" if source == "cooccurence" else source
        result[f"Evidence_{destination}"] = pd.NA
    for index, entity in result["entity_id"].astype(str).items():
        values = lookup.get(entity)
        if values is None:
            continue
        result.at[index, "Evidence_Provenance_Available"] = True
        result.at[index, "Evidence_Provenance_Status"] = "INDEXED_TARGET_EDGE"
        for source, value in values.items():
            destination = "cooccurrence" if source == "cooccurence" else source
            if source not in {"physical_status"} and value is not None:
                value = float(value) / 1000.0
            result.at[index, f"Evidence_{destination}"] = value
    experimental = result.get("Evidence_experiments", pd.Series(pd.NA, index=result.index)).fillna(0).gt(0) | result.get("Evidence_experiments_transferred", pd.Series(pd.NA, index=result.index)).fillna(0).gt(0)
    database = result.get("Evidence_database", pd.Series(pd.NA, index=result.index)).fillna(0).gt(0) | result.get("Evidence_database_transferred", pd.Series(pd.NA, index=result.index)).fillna(0).gt(0)
    transferred = sum((pd.to_numeric(result.get(f"Evidence_{column}", pd.Series(pd.NA, index=result.index)), errors="coerce").fillna(0) for column in ("textmining_transferred", "experiments_transferred", "database_transferred", "coexpression_transferred", "neighborhood_transferred")), start=pd.Series(0.0, index=result.index))
    result["Evidence_Experimental_Support_Present"] = experimental.where(result["Evidence_Provenance_Available"], pd.NA)
    result["Evidence_Database_Support_Present"] = database.where(result["Evidence_Provenance_Available"], pd.NA)
    result["Evidence_Transferred_Channel_Sum"] = transferred.where(result["Evidence_Provenance_Available"], pd.NA)
    return result, {"status": "AVAILABLE", "indexed_candidates": len(lookup), "index_path": str(EVIDENCE_DB_PATH)}


def _complex_label(row: pd.Series) -> str:
    if int(row.get("Response_Members", 0) or 0) == 0:
        return "LOW RESPONSE / STABLE"
    coherence = float(row.get("Coherence", 0.0) or 0.0)
    balance = float(row.get("Signed_Balance", 0.0) or 0.0)
    if coherence >= 0.5 and balance > 0:
        return "COORDINATED POSITIVE REDISTRIBUTION"
    if coherence >= 0.5 and balance < 0:
        return "COORDINATED NEGATIVE REDISTRIBUTION"
    return "MIXED REDISTRIBUTION"


def _complex_context(*, target: str, symbols: Mapping[str, str], classic_raw: pd.DataFrame, directed_raw: pd.DataFrame,
                     candidates: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, object]]:
    try:
        from src.evidence.complex_response import target_complex_response
        target_symbol = symbols.get(str(target), str(target))
        classic_input = classic_raw[["gene", "Delta_PageRank_Pct"]].copy(deep=True) if not classic_raw.empty else pd.DataFrame(columns=["gene", "Delta_PageRank_Pct"])
        directed_input = (directed_raw[["entity", "Directed_Redistribution_Pct"]].rename(columns={"entity": "gene", "Directed_Redistribution_Pct": "Delta_PageRank_Pct"}).copy(deep=True) if not directed_raw.empty else pd.DataFrame(columns=["gene", "Delta_PageRank_Pct"]))
        classic = target_complex_response([target_symbol], classic_input, dict(symbols))
        directed = target_complex_response([target_symbol], directed_input, dict(symbols))
        for frame in (classic, directed):
            if not frame.empty:
                frame["Response_Label"] = frame.apply(_complex_label, axis=1)
        if classic.empty:
            classic = pd.DataFrame(columns=["Complex_ID", "Complex", "Members", "Response_Label"])
        if directed.empty:
            directed = pd.DataFrame(columns=["Complex_ID", "Complex", "Members", "Response_Label"])
        left = classic.add_prefix("Classic_").rename(columns={"Classic_Complex_ID": "Complex_ID", "Classic_Complex": "Complex"})
        right = directed.add_prefix("Directed_").rename(columns={"Directed_Complex_ID": "Complex_ID", "Directed_Complex": "Complex"})
        context = left.merge(right, on=["Complex_ID", "Complex"], how="outer")
        if not context.empty:
            def cross(row: pd.Series) -> str:
                a, b = row.get("Classic_Response_Label"), row.get("Directed_Response_Label")
                if pd.isna(a) or pd.isna(b):
                    return "MIXED / UNCERTAIN"
                if a == b and a != "MIXED REDISTRIBUTION":
                    return "CROSS-MODEL CONSISTENT COMPLEX RESPONSE"
                if a != b:
                    return "DIRECTION-SENSITIVE COMPLEX RESPONSE"
                return "MIXED / UNCERTAIN"
            context["Cross_Engine_Complex_Class"] = context.apply(cross, axis=1)
        result = candidates.copy(deep=True)
        memberships: dict[str, list[str]] = {}
        for _, row in context.iterrows():
            members = str(row.get("Classic_Members", row.get("Directed_Members", ""))).split(";")
            for member in (value.upper() for value in members if value):
                memberships.setdefault(member, []).append(str(row.get("Complex") or row.get("Complex_ID")))
        if "symbol" in result:
            result["CORUM_Target_Complex_Context"] = result["symbol"].astype(str).str.upper().map(
                lambda value: ";".join(sorted(set(memberships.get(value, []))) )
            )
        else:
            result["CORUM_Target_Complex_Context"] = pd.NA
        result["CORUM_Target_Complex_Member"] = result["CORUM_Target_Complex_Context"].astype(str).ne("") if "CORUM_Target_Complex_Context" in result else pd.NA
        return result, context, {"status": "AVAILABLE", "target_complexes": len(context)}
    except Exception as error:
        return candidates.copy(deep=True), pd.DataFrame(), {"status": "UNAVAILABLE", "reason": f"{type(error).__name__}: {error}"}


def _provenance(*, target: str, tissue: str, block_strength: float, damping: float, directed_metadata: Mapping[str, object], cache_key: str) -> dict[str, object]:
    from src import config as classic_config
    from src.evidence.config import CORUM_PATH, EVIDENCE_DB_PATH, STRING_VERSION
    try:
        from src.evidence.scores import index_metadata
        evidence_index = index_metadata() if EVIDENCE_DB_PATH.exists() else None
    except Exception:
        evidence_index = None
    classic_configuration = {"threshold": classic_config.HIGH_CONF_THRESHOLD, "pagerank_damping": damping, "block_weight_fraction": block_strength, "edge_modifier_mode": "legacy_edge_modifiers"}
    directed_configuration = {
        "direction_policy": directed_metadata.get("direction_policy", "omnipath_direction_when_available"),
        "fallback_policy": directed_metadata.get("fallback_policy", "bidirectional_equal_weight"),
        "direction_dataset_version": directed_metadata.get("direction_dataset_version"),
        "effective_graph_fingerprint": directed_metadata.get("classic_effective_graph_fingerprint"),
    }
    product_manifest: dict[str, object]
    try:
        from src.product.provenance import build_run_manifest

        snapshot = product_data_snapshot()
        product_manifest = build_run_manifest(
            target=target,
            tissue=tissue,
            dataset_versions=snapshot.versions,
            dataset_fingerprints=snapshot.fingerprints,
            capabilities=snapshot.capability_report,
            engine_config_fingerprints={
                "classic": _json_hash(classic_configuration),
                "directed": _json_hash(directed_configuration),
            },
        )
        product_manifest["status"] = "READY"
    except Exception as error:
        product_manifest = {
            "status": "UNAVAILABLE",
            "reason": f"{type(error).__name__}: {error}",
        }
    return {
        "unified_engine_version": UNIFIED_ENGINE_VERSION, "timestamp": datetime.now(timezone.utc).isoformat(), "target": target, "tissue": tissue,
        "classic_engine_version": "Classic existing API", "directed_engine_version": "Directed BETA existing API",
        "classic_config_fingerprint": _json_hash(classic_configuration), "directed_config_fingerprint": _json_hash(directed_configuration),
        "classic_configuration": classic_configuration, "directed_configuration": directed_configuration,
        "string_version": STRING_VERSION, "evidence_channel_index": evidence_index,
        "omnipath": {key: value for key, value in directed_metadata.items() if "direction" in str(key).lower() or "omnipath" in str(key).lower()},
        "evidence_weighting_status": EVIDENCE_WEIGHTING_STATUS, "text_multiplier": "OFF", "physical_multiplier": "OFF",
        "corum_ranking_modifier": False, "corum_context_source": _safe_identity(CORUM_PATH), "cache_key": cache_key,
        "run_manifest": product_manifest,
    }


def _default_engine_runs(*, target: str, tissue: str, block_strength: float, damping: float, bc_sample_sources: int | None) -> tuple[pd.DataFrame, pd.DataFrame, Mapping[str, object], dict[str, str], dict[str, str]]:
    """Call existing engines on independent graph copies; no formula duplication."""
    from src import biology_logic, main as classic_main
    from src.services.directed_analysis import run_optional_directed_engine
    from src import config as classic_config

    preparation_configuration = {
        name: value for name, value in vars(classic_config).items()
        if name.isupper() and isinstance(value, (bool, int, float, str, tuple, list))
    }

    prepared_key = (
        "classic-prepared-context-v1", tissue, bc_sample_sources, "within_tissue",
        file_identity(classic_config.DB_PATH), file_identity(classic_config.OUT_CSV),
        stable_fingerprint(preparation_configuration),
    )
    graph, scores = _PREPARED_CONTEXT_CACHE.get_or_build(
        prepared_key,
        lambda: classic_main.arayuz_icin_motoru_hazirla(
            forced_genes=None,
            hedef_doku=None if tissue == "None" else tissue,
            bc_sample_sources=bc_sample_sources,
        ),
    )
    if target not in set(graph.vs["name"]):
        # Preserve the frozen engine's low-confidence target rescue semantics.
        graph, scores = classic_main.arayuz_icin_motoru_hazirla(
            forced_genes=[target],
            hedef_doku=None if tissue == "None" else tissue,
            bc_sample_sources=bc_sample_sources,
        )
    symbols = _symbol_mapping()
    errors: dict[str, str] = {}
    classic = pd.DataFrame()
    try:
        with tempfile.TemporaryDirectory(prefix="sophiark-unified-") as directory:
            biology_logic.run_infection_simulation(
                graph.copy(), scores.copy(deep=True), Path(directory), spesifik_hedefler=[target],
                block_weight_fraction=block_strength, damping=damping, tissue=None if tissue == "None" else tissue,
                edge_evidence_weighting_mode="legacy_edge_modifiers",
            )
            classic = biology_logic.latest_signed_redistribution()
        if not isinstance(classic, pd.DataFrame) or classic.empty:
            errors["classic"] = "Classic engine returned no signed response universe"
            classic = pd.DataFrame()
    except Exception as error:
        errors["classic"] = f"{type(error).__name__}: {error}"
        classic = pd.DataFrame()
    bundle = run_optional_directed_engine(
        graph=graph.copy(), targets=[target], tissue=tissue, classic_report=classic if not classic.empty else None,
        project_root=Path(__file__).resolve().parents[2], mode="Compare" if not classic.empty else "Directed", block_weight_fraction=block_strength,
        damping=damping, bc_sample_sources=bc_sample_sources, scores=scores.copy(deep=True), gene_to_symbol=symbols,
        string_dataset_version="local_string", edge_evidence_weighting_mode="legacy_edge_modifiers",
    )
    directed = bundle.calculation.report.copy(deep=True) if not bundle.calculation.report.empty else pd.DataFrame()
    metadata = {**dict(bundle.calculation.metadata), "status": bundle.calculation.status.value, "error": bundle.calculation.error}
    if bundle.calculation.error:
        errors["directed"] = bundle.calculation.error
    return classic.copy(deep=True), directed, metadata, symbols, errors


def _cache_key(*, target: str, tissue: str, block_strength: float, damping: float, convention: PresentationConvention) -> str:
    return _json_hash({"namespace": UNIFIED_ENGINE_VERSION, "target": target, "tissue": tissue, "block": block_strength, "damping": damping, "convention": asdict(convention)})


def clear_unified_cache() -> None:
    with _CACHE_LOCK:
        _CACHE.clear()
    _PREPARED_CONTEXT_CACHE.clear()


def run_unified_analysis(*, target: str, tissue: str, block_strength: float = 0.001, damping: float = 0.85,
                         bc_sample_sources: int | None = 4, convention: PresentationConvention = PresentationConvention(),
                         use_cache: bool = True,
                         engine_runner: Callable[..., tuple[pd.DataFrame, pd.DataFrame, Mapping[str, object], dict[str, str], dict[str, str]]] | None = None) -> UnifiedResearchReport:
    """Run Classic and Directed separately, then add non-calculating context layers."""
    target = str(target).replace("9606.", "")
    tissue = str(tissue or "None")
    key = _cache_key(target=target, tissue=tissue, block_strength=float(block_strength), damping=float(damping), convention=convention)
    if use_cache and engine_runner is None:
        with _CACHE_LOCK:
            cached = _CACHE.get(key)
            if cached is not None:
                _CACHE.move_to_end(key)
                return cached.copy()
    runner = engine_runner or _default_engine_runs
    classic_raw = pd.DataFrame(); directed_raw = pd.DataFrame(); directed_metadata: Mapping[str, object] = {}; errors: dict[str, str] = {}
    symbols = _symbol_mapping()
    try:
        classic_raw, directed_raw, directed_metadata, runner_symbols, runner_errors = runner(
            target=target, tissue=tissue, block_strength=float(block_strength), damping=float(damping), bc_sample_sources=bc_sample_sources,
        )
        symbols = runner_symbols or symbols
        errors.update({key: value for key, value in dict(runner_errors).items() if value})
    except Exception as error:
        errors["engine_orchestration"] = f"{type(error).__name__}: {error}"
    classic_ok = isinstance(classic_raw, pd.DataFrame) and not classic_raw.empty
    directed_ok = isinstance(directed_raw, pd.DataFrame) and not directed_raw.empty
    if not classic_ok and not directed_ok:
        status = UnifiedStatus.FAILED
    elif classic_ok and directed_ok:
        status = UnifiedStatus.COMPLETE
    else:
        status = UnifiedStatus.ENGINE_PARTIAL
    try:
        aligned, agreement = _synthesize(target=target, tissue=tissue, classic_raw=classic_raw, directed_raw=directed_raw, symbols=symbols, convention=convention)
    except Exception as error:
        errors["synthesis"] = f"{type(error).__name__}: {error}"
        aligned = pd.DataFrame(); agreement = {"status": "UNAVAILABLE", "reason": errors["synthesis"]}
        if classic_ok or directed_ok:
            status = UnifiedStatus.ENGINE_PARTIAL
        else:
            status = UnifiedStatus.FAILED
    important = aligned.loc[aligned.get("finding_class", pd.Series(dtype=str)).ne("LOW_AGREEMENT")].copy() if not aligned.empty else pd.DataFrame()
    if len(important) > 100:
        important = important.head(100).copy()
    provenance_candidates, evidence_status = _target_edge_provenance(important, target=target)
    candidates, complex_context, corum_status = _complex_context(target=target, symbols=symbols, classic_raw=classic_raw, directed_raw=directed_raw, candidates=provenance_candidates)
    provenance = _provenance(target=target, tissue=tissue, block_strength=float(block_strength), damping=float(damping), directed_metadata=directed_metadata, cache_key=key)
    result = UnifiedResearchReport(
        target=target, tissue=tissue, status=status,
        classic_status="COMPLETE" if classic_ok else "FAILED",
        directed_status=str(directed_metadata.get("status", "COMPLETE" if directed_ok else "FAILED")),
        candidates=candidates, classic_raw=classic_raw.copy(deep=True), directed_raw=directed_raw.copy(deep=True), complex_context=complex_context,
        agreement=agreement, provenance=provenance,
        context_availability={"evidence_provenance": evidence_status, "corum_complex_context": corum_status}, errors=errors,
    )
    if use_cache and engine_runner is None and result.status is not UnifiedStatus.FAILED:
        with _CACHE_LOCK:
            _CACHE[key] = result.copy(); _CACHE.move_to_end(key)
            while len(_CACHE) > _CACHE_MAXSIZE:
                _CACHE.popitem(last=False)
    return result


def export_unified_report(report: UnifiedResearchReport, destination: str | Path) -> dict[str, Path]:
    """Write machine-readable report/context exports without changing engine outputs."""
    # Match `.sophiark` portability: machine-specific source paths are
    # provenance metadata, not part of the portable export contract.
    from src.product.analysis_bundle import _portable

    folder = Path(destination)
    folder.mkdir(parents=True, exist_ok=True)
    paths = {
        "report": folder / "unified_report.json", "candidates": folder / "unified_candidates.csv",
        "complex_context": folder / "unified_complex_context.csv", "classic_raw": folder / "unified_classic_raw.csv",
        "directed_raw": folder / "unified_directed_raw.csv",
    }
    payload = _portable({
        "target": report.target, "tissue": report.tissue, "status": report.status.value,
        "classic_status": report.classic_status, "directed_status": report.directed_status,
        "agreement": dict(report.agreement), "provenance": dict(report.provenance),
        "context_availability": dict(report.context_availability), "errors": dict(report.errors),
    })
    paths["report"].write_text(json.dumps(payload, indent=2, default=str), encoding="utf-8")
    report.candidates.to_csv(paths["candidates"], index=False, encoding="utf-8-sig")
    report.complex_context.to_csv(paths["complex_context"], index=False, encoding="utf-8-sig")
    report.classic_raw.to_csv(paths["classic_raw"], index=False, encoding="utf-8-sig")
    report.directed_raw.to_csv(paths["directed_raw"], index=False, encoding="utf-8-sig")
    return paths
