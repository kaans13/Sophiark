from __future__ import annotations

import pandas as pd


def compare_classic_evidence(classic: pd.DataFrame, evidence: pd.DataFrame) -> pd.DataFrame:
    """Outer comparison; missing responses remain missing and are never coerced to zero."""
    c = classic[["gene", "Delta_PageRank_Pct"]].drop_duplicates("gene").rename(
        columns={"Delta_PageRank_Pct": "Classic_Response"}
    )
    e_columns = [c for c in (
        "gene", "Delta_PageRank_Pct", "Non_Text_Evidence", "Text_Dependency",
        "Physical_Supported_Edges", "Physical_Coverage",
    ) if c in evidence]
    e = evidence[e_columns].drop_duplicates("gene").rename(
        columns={"Delta_PageRank_Pct": "Evidence_Response"}
    )
    c["Classic_Rank"] = pd.to_numeric(c["Classic_Response"], errors="coerce").abs().rank(method="min", ascending=False).astype("Int64")
    e["Evidence_Rank"] = pd.to_numeric(e["Evidence_Response"], errors="coerce").abs().rank(method="min", ascending=False).astype("Int64")
    result = c.merge(e, on="gene", how="outer", validate="one_to_one")
    result["Response_Delta"] = result["Evidence_Response"] - result["Classic_Response"]
    result["Rank_Shift"] = result["Classic_Rank"] - result["Evidence_Rank"]
    return result
