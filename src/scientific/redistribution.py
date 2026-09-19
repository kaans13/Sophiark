"""İki yönlü PageRank redistribution sınıflandırması."""

from __future__ import annotations

from collections.abc import Mapping, Sequence

import numpy as np
import pandas as pd

from .statistics import benjamini_hochberg


def classify_redistribution(
    baseline: Mapping[str, float],
    perturbed: Mapping[str, float],
    *,
    empirical_p: Mapping[str, float] | None,
    target_names: Sequence[str],
    selection_mode: str,
    exploratory_top_n: int,
    percentile: float,
    fdr_alpha: float,
) -> pd.DataFrame:
    targets = set(map(str, target_names))
    rows: list[dict[str, object]] = []
    for gene in baseline:
        if gene in targets:
            continue
        before = float(baseline.get(gene, 0.0))
        after = float(perturbed.get(gene, 0.0))
        delta = after - before
        pct = delta / before * 100.0 if before > 0 else np.nan
        rows.append({
            "gene": gene,
            "PageRank_Baseline": before,
            "PageRank_Perturbed": after,
            "Delta_PageRank": delta,
            "Delta_PageRank_Pct": pct,
            "Abs_Delta_PageRank": abs(delta),
            "Response_Direction": "Influence Gain" if delta > 0 else ("Influence Loss" if delta < 0 else "No Change"),
            "empirical_p": float((empirical_p or {}).get(gene, 1.0)),
        })
    frame = pd.DataFrame(rows)
    if frame.empty:
        return frame
    frame["q_value"] = benjamini_hochberg(frame["empirical_p"].to_numpy())
    frame["significant_redistribution"] = frame["q_value"] <= float(fdr_alpha)
    cutoff = float(np.percentile(frame["Abs_Delta_PageRank"], percentile))
    frame["Strong_Redistribution"] = (
        (frame["Abs_Delta_PageRank"] >= cutoff) & (frame["Abs_Delta_PageRank"] > 0)
    )
    ordered = frame.sort_values(
        ["significant_redistribution", "Abs_Delta_PageRank"], ascending=[False, False], kind="mergesort"
    )
    if selection_mode == "percentile":
        selected = ordered[ordered["Strong_Redistribution"]].copy()
        selected["Selection_Reason"] = f"≥ {percentile:g}. |ΔPageRank| persentili"
    elif selection_mode == "top_n":
        selected = ordered.head(exploratory_top_n).copy()
        selected["Selection_Reason"] = f"Keşifsel Top-{exploratory_top_n} |ΔPageRank|"
    else:
        significant = ordered[ordered["significant_redistribution"]].copy()
        significant["Selection_Reason"] = f"Empirical null + BH-FDR ≤ {fdr_alpha:g}"
        exploratory = ordered[~ordered.index.isin(significant.index)].head(exploratory_top_n).copy()
        exploratory["Selection_Reason"] = f"FDR dışı keşifsel Top-{exploratory_top_n}"
        selected = pd.concat([significant, exploratory]).drop_duplicates("gene")
    return selected.reset_index(drop=True)


def top_positive_mean_pct(frame: pd.DataFrame, top_n: int = 100) -> float:
    if frame.empty or "Delta_PageRank_Pct" not in frame:
        return 0.0
    values = frame.loc[
        frame["Response_Direction"] == "Influence Gain", "Delta_PageRank_Pct"
    ].replace([np.inf, -np.inf], np.nan).dropna().nlargest(top_n)
    return float(values.mean()) if not values.empty else 0.0
