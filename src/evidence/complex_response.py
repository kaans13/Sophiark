from __future__ import annotations

import pandas as pd
import numpy as np

from .config import CORUM_PATH


def target_complex_response(target_symbols: list[str], redistribution: pd.DataFrame,
                            ensp_to_symbol: dict[str, str], zero_tolerance: float = 0.05) -> pd.DataFrame:
    corum = pd.read_csv(CORUM_PATH, sep="\t", dtype=str)
    target_set = {s.upper() for s in target_symbols if s}
    response = redistribution.copy()
    response["symbol"] = response["gene"].astype(str).map(ensp_to_symbol)
    response_map = dict(zip(response["symbol"].str.upper(), pd.to_numeric(response["Delta_PageRank_Pct"], errors="coerce")))
    mapped_symbols = {str(value).upper() for value in ensp_to_symbol.values()}
    rows = []
    for record in corum.to_dict("records"):
        members = [x.strip().upper() for x in str(record.get("subunits(Gene name)", "")).split(";") if x.strip() and x != "None"]
        if not target_set.intersection(members):
            continue
        values = np.asarray([response_map[m] for m in members if m in response_map and pd.notna(response_map[m])], dtype=float)
        denom = float(np.abs(values).sum()) if len(values) else 0.0
        rows.append({
            "Complex_ID": record.get("ComplexID"), "Complex": record.get("ComplexName"),
            "Total_Members": len(members), "Mapped_Members": sum(m in mapped_symbols for m in members),
            "Response_Members": len(values), "Perturbed_Members": len(target_set.intersection(members)),
            "Positive_Response_Members": int((values > zero_tolerance).sum()),
            "Negative_Response_Members": int((values < -zero_tolerance).sum()),
            "Low_Response_Members": int((np.abs(values) <= zero_tolerance).sum()),
            "Unresolved_Members": len(members) - len(values),
            "Median_Abs_Response": float(np.median(np.abs(values))) if len(values) else 0.0,
            "Mean_Abs_Response": float(np.mean(np.abs(values))) if len(values) else 0.0,
            "Max_Abs_Response": float(np.max(np.abs(values))) if len(values) else 0.0,
            "Coherence": abs(float(values.sum())) / denom if denom else 0.0,
            "Signed_Balance": float(values.sum()) / denom if denom else 0.0,
            "Members": ";".join(members),
        })
    return pd.DataFrame(rows)
