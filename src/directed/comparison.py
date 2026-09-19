"""Classic versus directed redistribution comparison."""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .models import DirectedCalculationResult, EngineComparisonResult


def _classic_columns(frame: pd.DataFrame) -> tuple[str, str]:
    entity = "gene" if "gene" in frame else ("entity" if "entity" in frame else "")
    response = "Delta_PageRank_Pct" if "Delta_PageRank_Pct" in frame else ""
    if not entity or not response:
        raise ValueError("Classic report requires gene/entity and Delta_PageRank_Pct")
    return entity, response


def compare_classic_and_directed(
    classic_report: pd.DataFrame, directed: DirectedCalculationResult,
    *, top_ks: tuple[int, ...] = (50, 100, 300),
) -> EngineComparisonResult:
    entity_col, response_col = _classic_columns(classic_report)
    classic = classic_report[[entity_col, response_col]].copy()
    classic.columns = ["entity", "Classic_Redistribution_Pct"]
    classic = classic.drop_duplicates("entity", keep="first")
    selected = directed.report[[
        "entity", "Directed_Redistribution_Pct", "Directed_BC",
        "Direction_Relation", "Direction_Support_Type", "Directed_Hop_Distance",
        "Directed_Distance",
        "Sign_Context", "Direction_Source_Coverage", "Alternative_Route_Status",
    ]].copy()
    classic_entities = set(classic["entity"].astype(str))
    directed_entities = set(selected["entity"].astype(str))
    missing_classic = directed_entities - classic_entities
    missing_directed = classic_entities - directed_entities
    if missing_classic or missing_directed:
        raise ValueError(
            "COMPARISON_INCOMPLETE: full Classic and Directed response universes are required; "
            f"missing_classic={len(missing_classic)}, missing_directed={len(missing_directed)}"
        )
    table = classic.merge(selected, on="entity", how="inner", validate="one_to_one")
    table["Classic_Redistribution_Pct"] = pd.to_numeric(
        table["Classic_Redistribution_Pct"], errors="raise"
    )
    table["Directed_Redistribution_Pct"] = pd.to_numeric(
        table["Directed_Redistribution_Pct"], errors="raise"
    )
    if table[["Classic_Redistribution_Pct", "Directed_Redistribution_Pct"]].isna().any().any():
        raise ValueError("COMPARISON_INCOMPLETE: missing response values cannot be interpreted as zero")
    table["Response_Delta"] = table["Directed_Redistribution_Pct"] - table["Classic_Redistribution_Pct"]
    table["Classic_Rank"] = table["Classic_Redistribution_Pct"].abs().rank(method="min", ascending=False).astype(int)
    table["Directed_Rank"] = table["Directed_Redistribution_Pct"].abs().rank(method="min", ascending=False).astype(int)
    table["Rank_Shift"] = table["Classic_Rank"] - table["Directed_Rank"]
    table = table.sort_values("Directed_Rank", kind="mergesort").reset_index(drop=True)
    if table["Classic_Rank"].nunique() <= 1 or table["Directed_Rank"].nunique() <= 1:
        rho = 1.0 if table["Classic_Rank"].equals(table["Directed_Rank"]) else 0.0
    else:
        rho = spearmanr(table["Classic_Rank"], table["Directed_Rank"]).statistic
    overlaps = {}
    for k in top_ks:
        classic_top = set(table.nsmallest(min(k, len(table)), "Classic_Rank")["entity"])
        directed_top = set(table.nsmallest(min(k, len(table)), "Directed_Rank")["entity"])
        overlaps[f"overlap_at_{k}"] = len(classic_top & directed_top)
        overlaps[f"overlap_fraction_at_{k}"] = len(classic_top & directed_top) / min(k, len(table)) if len(table) else 1.0
        overlaps[f"jaccard_at_{k}"] = (
            len(classic_top & directed_top) / len(classic_top | directed_top)
            if classic_top | directed_top else 1.0
        )
    if (
        table["Classic_Redistribution_Pct"].nunique() <= 1
        or table["Directed_Redistribution_Pct"].nunique() <= 1
    ):
        pearson = (
            1.0 if table["Classic_Redistribution_Pct"].equals(table["Directed_Redistribution_Pct"])
            else 0.0
        )
    else:
        pearson = table["Classic_Redistribution_Pct"].corr(
            table["Directed_Redistribution_Pct"], method="pearson",
        )
    validation = {
        "spearman_rank_correlation": float(rho) if np.isfinite(rho) else 1.0,
        "pearson_response_correlation": float(pearson) if np.isfinite(pearson) else 1.0,
        **overlaps,
        "median_absolute_rank_shift": float(table["Rank_Shift"].abs().median()) if len(table) else 0.0,
        "mean_absolute_rank_shift": float(table["Rank_Shift"].abs().mean()) if len(table) else 0.0,
        "largest_rank_increases": tuple(table.nlargest(min(10, len(table)), "Rank_Shift")["entity"]),
        "largest_rank_decreases": tuple(table.nsmallest(min(10, len(table)), "Rank_Shift")["entity"]),
        "ranking_metric": "absolute redistribution percent",
    }
    return EngineComparisonResult(table=table, validation=validation)
