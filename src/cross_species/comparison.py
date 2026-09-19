from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from .orthology import OrthologyIndex


@dataclass(frozen=True)
class ComparisonConfig:
    string_threshold: int
    damping: float
    attenuation_fraction: float
    selection_mode: str
    algorithm_version: str


def _ranked(frame: pd.DataFrame, symbol_column: str, response_column: str) -> pd.DataFrame:
    result = frame[[symbol_column, response_column]].dropna().copy()
    result[symbol_column] = result[symbol_column].astype(str)
    result[response_column] = pd.to_numeric(result[response_column], errors="coerce")
    result = result.dropna(subset=[response_column]).sort_values(
        response_column, key=lambda values: values.abs(), ascending=False, kind="mergesort"
    )
    result = result.drop_duplicates(symbol_column, keep="first").reset_index(drop=True)
    result["raw_rank"] = np.arange(1, len(result) + 1)
    result["rank_percentile"] = result["raw_rank"] / max(1, len(result))
    return result


def compare_ortholog_responses(
    human: pd.DataFrame,
    mouse: pd.DataFrame,
    orthology: OrthologyIndex,
    *,
    human_config: ComparisonConfig,
    mouse_config: ComparisonConfig,
    human_symbol_column: str = "Symbol",
    mouse_symbol_column: str = "Symbol",
    response_column: str = "Delta_PageRank_Pct",
    top_ks: tuple[int, ...] = (20, 50, 100),
) -> dict:
    if human_config != mouse_config:
        raise ValueError("CONFIG_MISMATCH: Human and Mouse calculation configs must match")
    human_ranked = _ranked(human, human_symbol_column, response_column)
    mouse_ranked = _ranked(mouse, mouse_symbol_column, response_column)
    human_ranked["Mouse_Symbol"] = human_ranked[human_symbol_column].map(
        lambda symbol: (orthology.one_to_one(symbol).mouse_gene if orthology.one_to_one(symbol) else None)
    )
    projected_human = human_ranked.dropna(subset=["Mouse_Symbol"]).copy()
    non_bijective = projected_human["Mouse_Symbol"].duplicated(keep=False)
    excluded_non_bijective = int(non_bijective.sum())
    projected_human = projected_human.loc[~non_bijective].copy()
    projected = projected_human.merge(
        mouse_ranked, left_on="Mouse_Symbol", right_on=mouse_symbol_column,
        how="inner", suffixes=("_Human", "_Mouse"), validate="one_to_one",
    )
    h_response = f"{response_column}_Human"
    m_response = f"{response_column}_Mouse"
    projected["Direction_Consistency"] = np.select(
        [
            (projected[h_response] > 0) & (projected[m_response] > 0),
            (projected[h_response] < 0) & (projected[m_response] < 0),
            np.sign(projected[h_response]) != np.sign(projected[m_response]),
        ],
        ["CONCORDANT_POSITIVE", "CONCORDANT_NEGATIVE", "DISCORDANT"],
        default="SPECIES_SPECIFIC_OR_ZERO",
    )
    metrics: dict[str, float | int] = {
        "comparable_one_to_one_count": len(projected),
        "excluded_non_bijective_symbol_mappings": excluded_non_bijective,
        "concordant_positive": int((projected["Direction_Consistency"] == "CONCORDANT_POSITIVE").sum()),
        "concordant_negative": int((projected["Direction_Consistency"] == "CONCORDANT_NEGATIVE").sum()),
        "discordant": int((projected["Direction_Consistency"] == "DISCORDANT").sum()),
    }
    metrics["response_direction_concordance"] = (
        float((projected["Direction_Consistency"].str.startswith("CONCORDANT")).mean())
        if len(projected) else float("nan")
    )
    metrics["rank_percentile_spearman"] = (
        float(spearmanr(projected["rank_percentile_Human"], projected["rank_percentile_Mouse"]).statistic)
        if len(projected) >= 2 else float("nan")
    )
    mouse_symbols = set(mouse_ranked[mouse_symbol_column])
    for k in top_ks:
        human_top = set(human_ranked.head(k)["Mouse_Symbol"].dropna())
        mouse_top = set(mouse_ranked.head(k)[mouse_symbol_column])
        intersection = human_top & mouse_top
        union = human_top | mouse_top
        metrics[f"top_{k}_overlap"] = len(intersection)
        metrics[f"top_{k}_jaccard"] = len(intersection) / len(union) if union else float("nan")
        metrics[f"human_top_{k}_mapped_into_mouse"] = len(human_top & mouse_symbols)
    return {"status": "AVAILABLE", "config_match": True, "metrics": metrics, "matched_responses": projected}
