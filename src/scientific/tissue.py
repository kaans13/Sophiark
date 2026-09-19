"""Tür çapında ortak referanslı doku ifade normalizasyonu."""

from __future__ import annotations

from functools import lru_cache
from contextlib import closing
import sqlite3

import numpy as np
import pandas as pd


def _weighted_quantile(values: np.ndarray, counts: np.ndarray, quantile: float) -> float:
    order = np.argsort(values, kind="mergesort")
    values, counts = values[order], counts[order]
    cumulative = np.cumsum(counts)
    position = quantile * cumulative[-1]
    return float(values[np.searchsorted(cumulative, position, side="left")])


@lru_cache(maxsize=8)
def pooled_expression_reference(db_path: str, table: str) -> dict[str, float | str]:
    if table not in {"tissue_expression", "mouse_tissue_expression"}:
        raise ValueError("İzin verilmeyen expression tablosu")
    with closing(sqlite3.connect(db_path)) as connection:
        rows = connection.execute(
            f"SELECT expression_level, COUNT(*) FROM {table} "
            "WHERE expression_level > 0 GROUP BY expression_level ORDER BY expression_level"
        ).fetchall()
    if not rows:
        return {"q01": 0.0, "q99": 1.0, "method": "log1p_pooled_q01_q99"}
    values = np.asarray([float(row[0]) for row in rows], dtype=float)
    counts = np.asarray([int(row[1]) for row in rows], dtype=np.int64)
    logged = np.log1p(values)
    q01 = _weighted_quantile(logged, counts, 0.01)
    q99 = _weighted_quantile(logged, counts, 0.99)
    if q99 <= q01:
        q99 = q01 + 1.0
    return {"q01": q01, "q99": q99, "method": "log1p_pooled_q01_q99"}


def tissue_multiplier(
    expression_1: pd.Series,
    expression_2: pd.Series | None = None,
    *,
    mode: str,
    db_path: str,
    expression_table: str,
) -> tuple[pd.Series, dict[str, object]]:
    first = expression_1.astype(float)
    second = expression_2.astype(float) if expression_2 is not None else first
    values = np.sqrt(first.clip(lower=0.0) * second.clip(lower=0.0))
    if mode == "within_tissue":
        minimum, maximum = float(values.min()), float(values.max())
        if maximum > minimum:
            scaled = (values - minimum) / (maximum - minimum)
        else:
            scaled = pd.Series(1.0, index=values.index)
        return (0.1 + 0.9 * scaled), {
            "mode": mode, "reference_min": minimum, "reference_max": maximum,
            "cross_tissue_comparable": False,
        }
    if mode != "cross_tissue_comparable":
        raise ValueError("tissue normalization mode within_tissue veya cross_tissue_comparable olmalıdır")
    reference = pooled_expression_reference(str(db_path), expression_table)
    denominator = float(reference["q99"]) - float(reference["q01"])
    scaled_first = ((np.log1p(first.clip(lower=0.0)) - float(reference["q01"])) / denominator).clip(0.0, 1.0)
    scaled_second = ((np.log1p(second.clip(lower=0.0)) - float(reference["q01"])) / denominator).clip(0.0, 1.0)
    edge_scaled = np.sqrt(scaled_first * scaled_second)
    return (0.1 + 0.9 * edge_scaled), {
        "mode": mode, **reference, "cross_tissue_comparable": True,
        "reference_scope": f"species-wide pooled {expression_table}",
        "edge_aggregation": "geometric_mean_of_endpoint_shared_scaled_expression",
    }
