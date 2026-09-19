from __future__ import annotations

from dataclasses import dataclass
import numpy as np
import pandas as pd

from .config import EvidenceConfig
from .graph import text_bonus


@dataclass(frozen=True, slots=True)
class ModifierBiasDecision:
    modifier: str
    usable_as_multiplier: bool
    reason: str
    spearman_classic_degree: float | None
    spearman_classic_pagerank: float | None
    spearman_annotation_availability: float | None


def _rho(left: pd.Series, right: pd.Series) -> float | None:
    pair = pd.concat([pd.to_numeric(left, errors="coerce"), pd.to_numeric(right, errors="coerce")], axis=1).dropna()
    if len(pair) < 20 or pair.iloc[:, 0].nunique() < 2 or pair.iloc[:, 1].nunique() < 2:
        return None
    return float(pair.iloc[:, 0].corr(pair.iloc[:, 1], method="spearman"))


def audit_modifier_bias(nodes: pd.DataFrame, *, modifier_column: str,
                        degree_column: str = "Classic_Degree",
                        pagerank_column: str = "Classic_PageRank",
                        annotation_column: str = "Annotation_Availability",
                        systematic_threshold: float = .50) -> ModifierBiasDecision:
    """Gate a modifier: centrality/annotation association keeps it context-only."""
    required = {modifier_column, degree_column, pagerank_column, annotation_column}
    missing = required.difference(nodes.columns)
    if missing:
        return ModifierBiasDecision(modifier_column, False, f"audit unavailable; missing {sorted(missing)}", None, None, None)
    correlations = (
        _rho(nodes[modifier_column], nodes[degree_column]),
        _rho(nodes[modifier_column], nodes[pagerank_column]),
        _rho(nodes[modifier_column], nodes[annotation_column]),
    )
    strong = [value for value in correlations if value is not None and value >= systematic_threshold]
    usable = not strong
    reason = (
        "no strong systematic centrality/annotation association detected"
        if usable else "modifier mainly rewards central or well-annotated genes; retain as context only"
    )
    return ModifierBiasDecision(modifier_column, usable, reason, *correlations)


def promotion_audit(classic: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Classic-controlled degree/PageRank decile promotion table; positive means promoted."""
    base = classic[["gene", "Classic_Degree", "Classic_PageRank", "Classic_Rank"]].copy()
    merged = base.merge(evidence[["gene", "Evidence_Rank"]], on="gene", how="inner", validate="one_to_one")
    merged["Promotion"] = merged["Classic_Rank"] - merged["Evidence_Rank"]
    for source, target in (("Classic_Degree", "Degree_Decile"), ("Classic_PageRank", "PageRank_Decile")):
        merged[target] = pd.qcut(merged[source].rank(method="first"), 10, labels=[f"D{i}" for i in range(1, 11)])
    return merged


def mid_evidence_audit(edges: pd.DataFrame, meaningful_bonus: float = 1.001) -> pd.DataFrame:
    bonus_column = "proposed_text_bonus" if "proposed_text_bonus" in edges else "text_bonus"
    values = edges[["s_nontext", bonus_column]].dropna().rename(columns={bonus_column: "text_bonus"})
    values["Evidence_Band"] = pd.qcut(values["s_nontext"].rank(method="first"), 3, labels=["LOW", "MID", "HIGH"])
    return values.groupby("Evidence_Band", observed=True).agg(
        edges=("text_bonus", "size"),
        meaningful_bonus_fraction=("text_bonus", lambda s: float((s >= meaningful_bonus).mean())),
        median_bonus=("text_bonus", "median"),
        p95_bonus=("text_bonus", lambda s: float(np.quantile(s, .95))),
    ).reset_index()


def hill_sensitivity(s_nontext: pd.Series, text: pd.Series) -> pd.DataFrame:
    """Small predefined diagnostic sweep; it does not optimize against known biology."""
    e = pd.to_numeric(s_nontext, errors="coerce").fillna(0).to_numpy()
    t = pd.to_numeric(text, errors="coerce").fillna(0).to_numpy()
    rows = []
    bands = pd.qcut(pd.Series(e).rank(method="first"), 3, labels=["LOW", "MID", "HIGH"])
    for k_text, k_evidence in ((.4, .4), (.6, .6), (.8, .8)):
        for lam in (0.0, .025, .05, .075, .10):
            config = EvidenceConfig(hill_k_text=k_text, hill_k_evidence=k_evidence, lambda_text=lam)
            bonus = text_bonus(e, t, config)
            for band in ("LOW", "MID", "HIGH"):
                selected = bonus[np.asarray(bands == band)]
                rows.append({
                    "k_text": k_text, "k_evidence": k_evidence, "lambda_text": lam,
                    "band": band, "median_bonus": float(np.median(selected)),
                    "p95_bonus": float(np.quantile(selected, .95)), "max_bonus": float(selected.max()),
                    "meaningful_fraction": float((selected >= 1.001).mean()),
                })
    return pd.DataFrame(rows)
