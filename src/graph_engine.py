# ---------------------------------------------------------------------------
# graph_engine.py
# igraph tabanlı ağ inşası, PageRank, topluluk tespiti, congestion ve
# global efficiency hesapları.
# ---------------------------------------------------------------------------
import gc
import math
import sqlite3
from pathlib import Path
from typing import Optional

import igraph as ig
import numpy as np
import pandas as pd

from .config import (
    DB_FIX_MAP,
    LOUVAIN_SEED,
    RESCUE_MIN_PARTNERS,
    RESCUE_THRESHOLDS,
    log,
)
from .data_loader import _attach_tissue_weights, _load_tissue_pivot, filter_edges_by_tissue


# ===========================================================================
# YARDIMCI FONKSİYONLAR
# ===========================================================================

def _edge_attrs(score: int) -> dict:
    """Basit ağırlık formülü — rescue kenarlarda kullanılır."""
    w    = float(score / 1000.0) ** 2
    dist = 1000.0 / float(score) if score > 0 else float("inf")
    return {"weight": w, "distance": dist}


def _mask_inf_distances(g: ig.Graph) -> None:
    dists = np.array(g.es["distance"], dtype=float)
    finite_mask = np.isfinite(dists)
    if np.any(finite_mask):
        max_finite = float(np.max(dists[finite_mask]))
        replacement = max_finite * 2.0
    else:
        replacement = 1e9  # Sabit büyük değer
    g.es["distance"] = [
        float(d) if np.isfinite(d) else replacement
        for d in dists
    ]


def _lcc_node_set(g: ig.Graph) -> set:
    if g.vcount() == 0:
        return set()
    clusters = g.clusters()
    giant    = clusters.giant()
    return set(giant.vs["name"])


def _name_to_idx(g: ig.Graph, name: str) -> int:
    try:
        return g.vs.find(name=name).index
    except ValueError:
        return -1


def _restore_edges(g: ig.Graph, snapshot: dict) -> None:
    for eid, (ow, od) in snapshot.items():
        g.es[eid]["weight"]   = ow
        g.es[eid]["distance"] = od
    gc.collect()


def _get_paralogs(g: ig.Graph, target_idx: int, blast_min: int) -> list:
    paralogs = []
    for nb_idx in g.neighbors(target_idx):
        eid = g.get_eid(target_idx, nb_idx, error=False)
        if eid == -1:
            continue
        w        = g.es[eid]["weight"]
        proxy_sc = 1000.0 * math.sqrt(max(w, 0))
        if proxy_sc >= blast_min:
            paralogs.append(nb_idx)
    return paralogs


# ===========================================================================
# ADIM 2a: Graf İnşası + Rescue
# ===========================================================================

def build_high_conf_subgraph(
    db_path: str,
    threshold: int,
    forced_genes: list = None,
    hedef_doku: str = "Cerebral Cortex",
) -> ig.Graph:
    """
    igraph Graf inşası — v5.0

    K2 düzeltmesi: rescue kenarları simplify() SONRASI ekleniyor;
    bu sayede combine_edges ile rescue_bridge özniteliği çakışmıyor.
    """
    log.info(
        "Adım 2a: igraph alt-graf v5.0 (eşik=%d, doku=%s)…",
        threshold, hedef_doku or "TÜM",
    )

    con = sqlite3.connect(db_path)
    try:
        log.info("    [1/3] Etkileşimler çekiliyor (SQL)...")
        edges_df = pd.read_sql_query(
            "SELECT protein1, protein2, combined_score FROM interactions WHERE combined_score >= ?",
            con, params=(threshold,)
        )

        if hedef_doku:
            log.info("    [2/3] Doku pivot SQL'de hesaplanıyor (MAX agregasyonu)...")
            tissue_pivot = _load_tissue_pivot(con, hedef_doku)
        else:
            tissue_pivot = {}
    finally:
        con.close()

    if edges_df.empty:
        log.warning("⚠ Sorgu sonuç döndürmedi. Eşik çok kısıtlayıcı?")
        g = ig.Graph(directed=False)
        g.vs["name"] = []
        return g

    edges_df["protein1"] = edges_df["protein1"].astype(str).str.replace("9606.", "", regex=False)
    edges_df["protein2"] = edges_df["protein2"].astype(str).str.replace("9606.", "", regex=False)
    edges_df["combined_score"] = edges_df["combined_score"].astype(np.int32)
    log.info("  Toplam kenar (ham): %d", len(edges_df))

    if hedef_doku:
        edges_df = filter_edges_by_tissue(edges_df, tissue_pivot, hedef_doku)

    if edges_df.empty:
        log.warning("Seçilen doku için pozitif ifadeli protein çifti bulunamadı.")
        g = ig.Graph(directed=False)
        g.vs["name"] = []
        return g

    log.info("    [3/3] Doku ağırlıkları LEFT JOIN ile ekleniyor...")
    edges_df = _attach_tissue_weights(edges_df, tissue_pivot)
    del tissue_pivot
    gc.collect()

    deg_series   = pd.concat(
        [edges_df["protein1"], edges_df["protein2"]]
    ).value_counts()
    edges_df["ka"] = edges_df["protein1"].map(deg_series).fillna(1).astype(np.int32)
    edges_df["kb"] = edges_df["protein2"].map(deg_series).fillna(1).astype(np.int32)

    edges_df["w_raw"] = (edges_df["combined_score"] / 1000.0) ** 2
    # Hub penalty devre dışı (w_hub = w_raw); aktif etmek için
    # _edge_attrs_hub_penalty mantığını burada uygulayın
    edges_df["w_hub"] = edges_df["w_raw"]

    # fillna(1.0) sayesinde tissue_raw hiçbir zaman sıfır olmaz
    edges_df["tissue_raw"] = np.sqrt(
        edges_df["exp1"].clip(lower=0.01) * edges_df["exp2"].clip(lower=0.01)
    )

    t_min = edges_df["tissue_raw"].min()
    t_max = edges_df["tissue_raw"].max()
    if t_max > t_min:
        edges_df["tissue_multiplier"] = (
            0.1 + 0.9 * (edges_df["tissue_raw"] - t_min) / (t_max - t_min)
        )
    else:
        edges_df["tissue_multiplier"] = 1.0

    edges_df["final_weight"]   = (edges_df["w_hub"] * edges_df["tissue_multiplier"]).clip(lower=1e-9)
    edges_df["final_distance"] = (1.0 / edges_df["final_weight"]).replace(
        [np.inf, -np.inf], 1e9
    )

    all_proteins = pd.unique(
        pd.concat([edges_df["protein1"], edges_df["protein2"]])
    )
    name_to_idx  = {name: idx for idx, name in enumerate(all_proteins)}
    vertex_names = list(all_proteins)

    edges_df["u_idx"] = edges_df["protein1"].map(name_to_idx)
    edges_df["v_idx"] = edges_df["protein2"].map(name_to_idx)

    edge_list = list(zip(edges_df["u_idx"], edges_df["v_idx"]))
    weights   = edges_df["final_weight"].tolist()
    distances = edges_df["final_distance"].tolist()

    del edges_df
    gc.collect()

    g = ig.Graph(n=len(vertex_names), edges=edge_list, directed=False)
    g.vs["name"]         = vertex_names
    g.vs["lokalizasyon"] = ["Unknown"] * len(vertex_names)
    g.es["weight"]       = weights
    g.es["distance"]     = distances
    # K2: rescue_bridge başlangıçta False — simplify ÖNCE çalışıyor,
    # rescue kenarları SONRA ekleniyor; bu sayede combine_edges çakışması yok.
    g.es["rescue_bridge"] = [False] * len(weights)

    g = g.simplify(combine_edges={"weight": "max", "distance": "min", "rescue_bridge": "first"})
    name_to_idx = {v["name"]: v.index for v in g.vs}

    log.info("    [Inf Mask] Ana ağda sonsuz mesafeler maskeleniyor...")
    _mask_inf_distances(g)

    log.info(
        "  Ağ inşası tamamlandı: %d düğüm, %d kenar (simplify + inf-mask sonrası).",
        g.vcount(), g.ecount(),
    )

    # Rescue — simplify() SONRASI; K2 düzeltmesi
    if forced_genes:
        log.info(
            "  [forced_genes] %d gen için rescue (eşik listesi: %s)…",
            len(forced_genes), [threshold] + RESCUE_THRESHOLDS,
        )
        rescue_query = """
            SELECT protein1, protein2, combined_score
            FROM interactions
            WHERE (protein1 IN (?, ?) OR protein2 IN (?, ?))
              AND combined_score >= ?
            ORDER BY combined_score DESC LIMIT 50
        """
        con = sqlite3.connect(db_path)
        try:
            for gene in forced_genes:
                raw_id    = str(gene).replace("9606.", "")
                prefix_id = "9606." + raw_id
                alias_id  = DB_FIX_MAP.get(raw_id, raw_id)

                if raw_id in name_to_idx or prefix_id in name_to_idx or alias_id in name_to_idx:
                    log.info("    ✓ Zaten mevcut: %s", raw_id)
                    continue

                partners_added  = 0
                rescue_thr_used = None

                for thr_try in [threshold] + RESCUE_THRESHOLDS:
                    rows = con.execute(
                        rescue_query,
                        (raw_id, alias_id, raw_id, alias_id, thr_try),
                    ).fetchall()
                    if not rows:
                        continue
                    for r in rows[: RESCUE_MIN_PARTNERS * 3]:
                        src, dst, sc = r[0], r[1], r[2]
                        partner = dst if src in (raw_id, alias_id) else src
                        if raw_id not in name_to_idx:
                            g.add_vertex(name=raw_id, lokalizasyon="Unknown")
                            name_to_idx[raw_id] = g.vs[g.vcount() - 1].index
                        if partner not in name_to_idx:
                            g.add_vertex(name=partner, lokalizasyon="Unknown")
                            name_to_idx[partner] = g.vs[g.vcount() - 1].index
                        ui    = name_to_idx[raw_id]
                        vi    = name_to_idx[partner]
                        attrs = _edge_attrs(sc)
                        attrs["rescue_bridge"] = False
                        eid   = g.get_eid(ui, vi, error=False)
                        if eid == -1:
                            g.add_edge(ui, vi, **attrs)
                        partners_added += 1
                        if partners_added >= RESCUE_MIN_PARTNERS:
                            break
                    if partners_added >= RESCUE_MIN_PARTNERS:
                        rescue_thr_used = thr_try
                        break

                if partners_added == 0:
                    log.warning(
                        "    ⚠ %s izole kaldı. LCC yapay köprüsü atılıyor "
                        "(rescue_bridge=True — BC analizinde bu etiket izlenebilir)…",
                        raw_id,
                    )
                    lcc_names = _lcc_node_set(g)
                    if lcc_names:
                        bridge_target = max(lcc_names, key=lambda n: g.vs.find(name=n).degree())
                        if raw_id not in name_to_idx:
                            g.add_vertex(name=raw_id, lokalizasyon="Unknown")
                            name_to_idx[raw_id] = g.vs[g.vcount() - 1].index
                        ui = name_to_idx[raw_id]
                        vi = name_to_idx[bridge_target]
                        bridge_attrs = _edge_attrs(threshold)
                        bridge_attrs["rescue_bridge"] = True
                        g.add_edge(ui, vi, **bridge_attrs)
                    else:
                        if raw_id not in name_to_idx:
                            g.add_vertex(name=raw_id, lokalizasyon="Unknown")
                            name_to_idx[raw_id] = g.vs[g.vcount() - 1].index
                else:
                    log.info(
                        "    ✓ %s kurtarıldı: %d partner, eşik=%d",
                        raw_id, partners_added, rescue_thr_used,
                    )

                lcc_names = _lcc_node_set(g)
                gene_idx  = name_to_idx.get(raw_id, -1)
                if gene_idx != -1 and raw_id not in lcc_names and lcc_names:
                    lcc_neighbors = [
                        nb for nb in g.neighbors(gene_idx)
                        if g.vs[nb]["name"] in lcc_names
                    ]
                    if not lcc_neighbors:
                        bridge_target = max(lcc_names, key=lambda n: g.vs.find(name=n).degree())
                        vi = name_to_idx[bridge_target]
                        bridge_attrs = _edge_attrs(threshold)
                        bridge_attrs["rescue_bridge"] = True
                        g.add_edge(gene_idx, vi, **bridge_attrs)
                        log.info(
                            "    🔗 LCC yapay köprüsü (rescue_bridge=True): %s → %s",
                            raw_id, bridge_target,
                        )
        finally:
            con.close()

    _mask_inf_distances(g)

    rescue_bridge_count = sum(1 for e in g.es if e["rescue_bridge"])
    if rescue_bridge_count > 0:
        log.warning(
            "  ⚠ %d yapay köprü kenarı mevcut (rescue_bridge=True). "
            "BC sıralamasını hafifçe bozabilir — raporda etiketlendi.",
            rescue_bridge_count,
        )

    log.info("  → igraph Graf (final): %d düğüm, %d kenar.", g.vcount(), g.ecount())
    return g


# ===========================================================================
# ADIM 2b: PageRank
# ===========================================================================

def compute_pagerank(g: ig.Graph, damping: float, max_iter: int) -> dict:
    log.info("Adım 2b: igraph PageRank (d=%.2f, deterministic)…", damping)

    # Her hesaplama öncesi rastgeleliği sıfırla
    np.random.seed(42)

    try:
        # PRPACK deterministiktir
        pr_values = g.pagerank(
            damping=damping, 
            weights="weight", 
            directed=False,
            implementation="prpack"
        )
    except Exception:
        # PRPACK yoksa ARPACK kullan ama seed'i sabitlediğimiz için
        # daha kararlı olacaktır
        pr_values = g.pagerank(
            damping=damping, 
            weights="weight", 
            directed=False,
            implementation="arpack"
        )

    pr = {g.vs[i]["name"]: pr_values[i] for i in range(g.vcount())}
    log.info("  → %d düğüm.", len(pr))
    return pr



# ===========================================================================
# ADIM 3b: Topluluk Tespiti
# ===========================================================================

def run_community_detection(g: ig.Graph, df: pd.DataFrame) -> tuple:
    log.info("Adım 3b: igraph Louvain (seed=%d)…", LOUVAIN_SEED)

    # SEED UYUMLULUK DÜZELTMESİ:
    # igraph < 0.10: random.seed() ile
    # igraph >= 0.10: seed= parametresi ile
    try:
        communities = g.community_multilevel(weights="weight", seed=LOUVAIN_SEED)
    except TypeError:
        import random as rnd
        rnd.seed(LOUVAIN_SEED)
        communities = g.community_multilevel(weights="weight")

    node_comm: dict[str, int] = {}
    for cid, members in enumerate(communities):
        for vidx in members:
            node_comm[g.vs[vidx]["name"]] = cid

    df = df.copy()
    df["Topluluk_ID"] = df["gene"].map(node_comm).fillna(-1).astype(int)
    log.info("  → %d topluluk.", len(communities))
    return df, communities


# ===========================================================================
# GÖREV 2: Graf Düğüm Etiketleme
# ===========================================================================

def tag_graph_nodes(g: ig.Graph, loc_map: dict) -> ig.Graph:
    log.info("Görev 2: Mekansal etiketleme…")
    tagged = 0
    for v in g.vs:
        loc = loc_map.get(v["name"], "Unknown")
        v["lokalizasyon"] = loc
        if loc != "Unknown":
            tagged += 1
    log.info("  → %d / %d etiketlendi.", tagged, g.vcount())
    return g


# ===========================================================================
# Congestion Hesabı — Lokal Betweenness
# ===========================================================================

def _calculate_congestion(
    g_pre: ig.Graph,
    g_post: ig.Graph,
    gate_genes: Optional[list] = None,
    max_vertices: int = 3000,  # Artırıldı: 1000 → 3000
    hop_distance: int = 3,     # YENİ: 2 → 3 hop
) -> np.ndarray:
    """
    Lokal Betweenness Congestion — 1., 2. ve 3. derece komşular.

    v5.1: hop_distance=3 ile daha geniş sinyal yayılımı yakalanır.
    """
    pre_copy  = g_pre.copy()
    post_copy = g_post.copy()
    _mask_inf_distances(pre_copy)
    _mask_inf_distances(post_copy)

    n = pre_copy.vcount()
    delta_bc = np.zeros(n, dtype=float)

    name_to_idx_pre = {pre_copy.vs[i]["name"]: i for i in range(n)}

    if gate_genes:
        seed_indices = [
            name_to_idx_pre[gn] for gn in gate_genes if gn in name_to_idx_pre
        ]
    else:
        seed_indices = []

    if not seed_indices:
        degrees = pre_copy.degree()
        seed_indices = sorted(range(n), key=lambda i: degrees[i], reverse=True)[
            : max_vertices // 2
        ]

    # Genişletilmiş komşuluk (hop_distance kadar)
    target_set: set = set(seed_indices)
    current_frontier = set(seed_indices)

    hop = 0
    for hop in range(hop_distance):
        next_frontier = set()
        for node in current_frontier:
            neighbors = set(pre_copy.neighbors(node))
            new_nodes = neighbors - target_set
            next_frontier.update(new_nodes)
            target_set.update(new_nodes)

            if len(target_set) >= max_vertices:
                break

        current_frontier = next_frontier

        if len(target_set) >= max_vertices:
            log.info(
                "  Congestion: %d düğüm toplandı (%d hop'ta max_vertices aşıldı)",
                len(target_set), hop + 1
            )
            break

    target_list = sorted(target_set)[:max_vertices]
    if not target_list:
        return delta_bc

    log.info(
        "  Congestion BC: %d hedef düğüm üzerinde hesaplanıyor (hop=%d)…",
        len(target_list), hop_distance
    )

    bc_pre_local  = pre_copy.betweenness(
        vertices=target_list, weights="distance", directed=False
    )
    bc_post_local = post_copy.betweenness(
        vertices=target_list, weights="distance", directed=False
    )

    for rank, vidx in enumerate(target_list):
        delta_bc[vidx] = bc_post_local[rank] - bc_pre_local[rank]

    return delta_bc


# ===========================================================================
# Global Efficiency
# ===========================================================================

def compute_global_efficiency(g: ig.Graph, sample_nodes: list):
    n = len(sample_nodes)
    if n < 2:
        return 0.0, {}

    name_to_idx = {v["name"]: v.index for v in g.vs}
    src_indices = [name_to_idx[s] for s in sample_nodes if s in name_to_idx]
    if len(src_indices) < 2:
        return 0.0, {}

    dist_matrix = g.distances(source=src_indices, target=src_indices, weights="distance")

    total_inv  = 0.0
    pair_count = 0
    node_effs: dict[str, float] = {}

    for i, src_name in enumerate(sample_nodes):
        if src_name not in name_to_idx:
            node_effs[src_name] = 0.0
            continue
        row_sum = 0.0
        row_cnt = 0
        for j, tgt_name in enumerate(sample_nodes):
            if i == j:
                continue
            d = dist_matrix[i][j]
            if d != float("inf") and d > 0:
                inv        = 1.0 / d
                total_inv += inv
                row_sum   += inv
                row_cnt   += 1
                pair_count += 1
        node_effs[src_name] = row_sum / row_cnt if row_cnt > 0 else 0.0

    possible = n * (n - 1)
    return float(total_inv / possible) if possible > 0 else 0.0, node_effs
