"""İleri doğrulama panelleri: BC ve PageRank hassasiyeti."""

from __future__ import annotations

import igraph as ig
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def sampled_betweenness(
    graph: ig.Graph, *, sample_sources: int | None, random_seed: int = 42
) -> tuple[np.ndarray, dict[str, object]]:
    n = graph.vcount()
    if sample_sources is None or sample_sources >= n:
        values = np.asarray(graph.betweenness(weights="distance", directed=False), dtype=float)
        return values, {"mode": "Full Validation", "sample_size": n, "random_seed": random_seed}
    sources = sorted(np.random.default_rng(random_seed).choice(n, size=sample_sources, replace=False).tolist())
    values = np.asarray(graph.betweenness(weights="distance", directed=False, sources=sources), dtype=float)
    values *= n / sample_sources
    return values, {"mode": "Fast Approximate", "sample_size": sample_sources, "random_seed": random_seed}


def deep_bc_validation(
    graph: ig.Graph, *, sample_sources: int = 1200, top_k: int = 100, random_seed: int = 42
) -> dict[str, object]:
    approximate, approx_meta = sampled_betweenness(graph, sample_sources=sample_sources, random_seed=random_seed)
    exact, exact_meta = sampled_betweenness(graph, sample_sources=None, random_seed=random_seed)
    rho = float(spearmanr(approximate, exact).statistic)
    approx_top = set(np.argsort(approximate)[-top_k:])
    exact_top = set(np.argsort(exact)[-top_k:])
    jaccard = len(approx_top & exact_top) / len(approx_top | exact_top) if approx_top | exact_top else 1.0
    return {
        "spearman_rank_correlation": rho,
        "top_k": top_k,
        "top_k_jaccard": jaccard,
        "approximate": approx_meta,
        "full": exact_meta,
        "graph_size": {"nodes": graph.vcount(), "edges": graph.ecount()},
        "weighted": True,
    }


def pagerank_robustness(
    graph: ig.Graph,
    *,
    damping_grid: tuple[float, ...] = (0.75, 0.85, 0.95),
    reference_damping: float = 0.85,
    top_k: int = 100,
) -> tuple[pd.DataFrame, dict[str, object]]:
    ranks = {
        damping: np.asarray(graph.pagerank(damping=damping, weights="weight", directed=False), dtype=float)
        for damping in damping_grid
    }
    reference = ranks[min(damping_grid, key=lambda value: abs(value - reference_damping))]
    reference_top = set(np.argsort(reference)[-top_k:])
    rows = []
    reference_ranks = np.argsort(np.argsort(reference))
    for damping, values in ranks.items():
        top = set(np.argsort(values)[-top_k:])
        rows.append({
            "Damping": damping,
            "Spearman": float(spearmanr(reference, values).statistic),
            "Top_K_Jaccard": len(reference_top & top) / len(reference_top | top) if reference_top | top else 1.0,
            "Critical_Gene_Stability": len(reference_top & top) / len(reference_top | top) if reference_top | top else 1.0,
            "Median_Absolute_Rank_Shift": float(np.median(np.abs(reference_ranks - np.argsort(np.argsort(values))))),
            "P95_Absolute_Rank_Shift": float(np.percentile(np.abs(reference_ranks - np.argsort(np.argsort(values))), 95)),
            "Max_Absolute_Rank_Shift": int(np.max(np.abs(reference_ranks - np.argsort(np.argsort(values))))),
        })
    frame = pd.DataFrame(rows)
    min_rho, min_jaccard = float(frame["Spearman"].min()), float(frame["Top_K_Jaccard"].min())
    if min_rho >= 0.95 and min_jaccard >= 0.75:
        classification = "Stable"
    elif min_rho >= 0.85 and min_jaccard >= 0.50:
        classification = "Moderately Sensitive"
    else:
        classification = "Sensitive"
    return frame, {
        "classification": classification,
        "thresholds": {"stable": {"spearman": 0.95, "top_k_jaccard": 0.75}, "moderate": {"spearman": 0.85, "top_k_jaccard": 0.50}},
        "reference_damping": reference_damping,
        "top_k": top_k,
    }
