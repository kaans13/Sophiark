from __future__ import annotations

import sqlite3
from typing import Iterable

import numpy as np
import pandas as pd

from src.graph_engine import build_high_conf_subgraph, _mask_inf_distances
from .config import CLASSIC_DB_PATH, EVIDENCE_DB_PATH, EvidenceConfig, EvidenceMode


PROFILE_COLUMNS = (
    "neighborhood", "neighborhood_transferred", "fusion", "cooccurence", "homology",
    "coexpression", "coexpression_transferred", "experiments", "experiments_transferred",
    "database", "database_transferred", "textmining", "textmining_transferred",
    "combined_score", "reconstructed_full", "s_nontext", "text_dependency",
    "physical_score", "physical_status",
)


def hill(value: np.ndarray, k: float, exponent: float) -> np.ndarray:
    powered = np.power(np.clip(value, 0.0, 1.0), exponent)
    return powered / (float(k) ** exponent + powered)


def text_bonus(s_nontext: np.ndarray, text: np.ndarray, config: EvidenceConfig) -> np.ndarray:
    ht = hill(text, config.hill_k_text, config.hill_n)
    he = hill(s_nontext, config.hill_k_evidence, config.hill_m)
    if config.text_gate == "G2":
        gate = ht * np.sqrt(he)
    elif config.text_gate == "G3":
        gate = ht * (2.0 * he / (1.0 + he))
    else:
        gate = ht * he
    return 1.0 + config.lambda_text * gate


def _profiles_for_graph(graph) -> pd.DataFrame:
    pairs = []
    for edge in graph.es:
        a, b = sorted((graph.vs[edge.source]["name"], graph.vs[edge.target]["name"]))
        pairs.append((edge.index, a, b))
    con = sqlite3.connect(EVIDENCE_DB_PATH)
    try:
        con.execute("CREATE TEMP TABLE wanted (eid INTEGER, protein1 TEXT, protein2 TEXT, PRIMARY KEY(protein1,protein2))")
        con.executemany("INSERT OR IGNORE INTO wanted VALUES (?,?,?)", pairs)
        select = ",".join(f"e.{c}" for c in PROFILE_COLUMNS)
        return pd.read_sql_query(
            f"SELECT w.eid,w.protein1,w.protein2,{select} FROM wanted w LEFT JOIN evidence e USING(protein1,protein2) ORDER BY w.eid",
            con,
        )
    finally:
        con.close()


def build_evidence_graph(*, config: EvidenceConfig, forced_genes: list[str] | None = None,
                         tissue: str | None = None):
    """Build on a fresh Classic graph, then alter only Evidence-owned copies."""
    if not EVIDENCE_DB_PATH.exists():
        raise FileNotFoundError(f"Evidence index missing: {EVIDENCE_DB_PATH}")
    graph = build_high_conf_subgraph(
        str(CLASSIC_DB_PATH), config.string_threshold, forced_genes=forced_genes, hedef_doku=tissue,
    )
    profiles = _profiles_for_graph(graph)
    missing = profiles["combined_score"].isna()
    if missing.any() and config.mode is not EvidenceMode.PARITY:
        raise RuntimeError(
            f"{int(missing.sum())} rescued edge(s) fall below the indexed production threshold and lack "
            "validated channel profiles; Evidence refuses to infer S_nontext or borrow Classic weights"
        )
    for column in PROFILE_COLUMNS:
        values = profiles[column].tolist()
        graph.es[column] = values
    if config.mode is EvidenceMode.PARITY:
        graph.es["proposed_text_bonus"] = [1.0] * graph.ecount()
        graph.es["text_bonus"] = [1.0] * graph.ecount()
        graph.es["proposed_physical_bonus"] = [1.0] * graph.ecount()
        graph.es["physical_bonus"] = [1.0] * graph.ecount()
        graph["evidence_mode"] = config.mode.value
        return graph, profiles

    official = pd.to_numeric(profiles["combined_score"], errors="coerce") / 1000.0
    snt = pd.to_numeric(profiles["s_nontext"], errors="coerce")
    text_direct = pd.to_numeric(profiles["textmining"], errors="coerce").fillna(0) / 1000.0
    text_transfer = pd.to_numeric(profiles["textmining_transferred"], errors="coerce").fillna(0) / 1000.0
    # Direct and transferred text are corroboration inputs, combined with the same validated combiner.
    prior = config.string_prior
    direct_np, transfer_np = text_direct.to_numpy(), text_transfer.to_numpy()
    d = np.where(direct_np > 0, np.maximum(0, (direct_np-prior)/(1-prior)), 0)
    t = np.where(transfer_np > 0, np.maximum(0, (transfer_np-prior)/(1-prior)), 0)
    text = prior + (1-prior) * (1-(1-d)*(1-t))
    bonuses = text_bonus(snt.fillna(official).to_numpy(), text, config)
    physical_supported = profiles["physical_status"].eq("SUPPORTED").to_numpy(dtype=float)
    physical_bonus = 1.0 + config.lambda_physical * physical_supported
    old_weight = np.asarray(graph.es["weight"], dtype=float)
    official_values = official.to_numpy(dtype=float)
    official_values = np.where(np.isnan(official_values), np.sqrt(old_weight), official_values)
    snt_values = snt.to_numpy(dtype=float)
    snt_values = np.where(np.isnan(snt_values), official_values, snt_values)
    tissue_multiplier = old_weight / np.maximum(np.square(official_values), 1e-12)
    new_weight = np.square(snt_values) * tissue_multiplier * bonuses * physical_bonus
    new_weight = np.maximum(new_weight, 1e-9)
    graph.es["classic_weight"] = old_weight.tolist()
    graph.es["proposed_text_bonus"] = bonuses.tolist()
    graph.es["text_bonus"] = bonuses.tolist()
    graph.es["proposed_physical_bonus"] = physical_bonus.tolist()
    graph.es["physical_bonus"] = physical_bonus.tolist()
    graph.es["weight"] = new_weight.tolist()
    graph.es["distance"] = (1.0 / new_weight).tolist()
    if config.mode is EvidenceMode.E2_TOPOLOGY_WEIGHT:
        remove = [i for i, value in enumerate(snt.fillna(-1).to_numpy()) if value * 1000 < config.string_threshold]
        graph.delete_edges(remove)
        isolates = [v.index for v in graph.vs if v.degree() == 0]
        if isolates:
            graph.delete_vertices(isolates)
    _mask_inf_distances(graph)
    graph["evidence_mode"] = config.mode.value
    return graph, profiles
