# ---------------------------------------------------------------------------
# mouse_graph_engine.py
# Fare (Mus musculus) için ayrı graf motoru.
# İnsan graph_engine.py ile aynı matematiksel prensipleri kullanır,
# ancak VERİTABANI ve PREFIX olarak fare (10090.) kullanır.
# İnsan motoruna hiçbir bağımlılığı yoktur.
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
    LOUVAIN_SEED,
    RESCUE_MIN_PARTNERS,
    RESCUE_THRESHOLDS,
    DB_FIX_MAP,
    log,
)
from .data_loader import (
    _attach_tissue_weights,
    _load_tissue_pivot,
    filter_edges_by_tissue,
)
from ..scientific.efficiency import strengths_to_distances


# ===========================================================================
# YARDIMCI FONKSİYONLAR
# ===========================================================================

def _edge_attrs(score: int) -> dict:
    """Basit ağırlık formülü — rescue kenarlarda kullanılır."""
    w    = float(score / 1000.0) ** 2
    dist = 1000.0 / float(score) if score > 0 else float("inf")
    return {"weight": w, "distance": dist}


def _mask_inf_distances(g: ig.Graph) -> None:
    """Sonsuz mesafeleri, ağdaki en büyük sonlu mesafenin 2 katıyla değiştir."""
    dists = np.array(g.es["distance"], dtype=float)
    finite_mask = np.isfinite(dists)
    if np.any(finite_mask):
        max_finite = float(np.max(dists[finite_mask]))
        replacement = max_finite * 2.0
    else:
        replacement = 1e9
    g.es["distance"] = [
        float(d) if np.isfinite(d) else replacement
        for d in dists
    ]


def _lcc_node_set(g: ig.Graph) -> set:
    """En büyük bağlı bileşendeki (LCC) düğüm isimlerini döndür."""
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
    """Kenar ağırlıklarını ve mesafeleri snapshot'tan geri yükle."""
    for eid, (ow, od) in snapshot.items():
        g.es[eid]["weight"]   = ow
        g.es[eid]["distance"] = od
    gc.collect()


def _get_paralogs(g: ig.Graph, target_idx: int, blast_min: int) -> list:
    """Hedef düğümün yüksek ağırlıklı komşularını paralog olarak döndür."""
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
# ADIM 2a: Graf İnşası + Rescue (Fare Versiyonu)
# ===========================================================================

def build_high_conf_subgraph(
    db_path: str,
    threshold: int,
    forced_genes: list = None,
    hedef_doku: str = "Cerebral Cortex",
    tissue_normalization_mode: str = "within_tissue",
) -> ig.Graph:
    """
    igraph Graf inşası — Fare v1.0

    İnsan v5.0 ile aynı akış:
    - mouse_interactions tablosundan okur
    - 10090. prefix'ini temizler
    - Doku filtresi uygular
    - Doku ağırlıklandırması yapar
    - simplify + inf-maskeleme + rescue
    """
    log.info(
        "Mouse Adım 2a: igraph alt-graf (eşik=%d, doku=%s)…",
        threshold, hedef_doku or "TÜM",
    )

    con = sqlite3.connect(db_path)
    try:
        log.info("    [1/3] Fare etkileşimleri çekiliyor (SQL)...")
        edges_df = pd.read_sql_query(
            "SELECT protein1, protein2, combined_score FROM mouse_interactions WHERE combined_score >= ?",
            con, params=(threshold,)
        )

        if hedef_doku:
            log.info("    [2/3] Fare doku pivotu SQL'de hesaplanıyor (MAX agregasyonu)...")
            tissue_pivot = _load_tissue_pivot(con, hedef_doku)
        else:
            tissue_pivot = {}
    finally:
        con.close()

    if edges_df.empty:
        log.warning("⚠ Fare sorgu sonuç döndürmedi. Eşik çok kısıtlayıcı?")
        g = ig.Graph(directed=False)
        g.vs["name"] = []
        return g

    # Fare STRING verisindeki 10090. prefix'ini temizle
    edges_df["protein1"] = (
        edges_df["protein1"].astype(str).str.replace("10090.", "", regex=False)
    )
    edges_df["protein2"] = (
        edges_df["protein2"].astype(str).str.replace("10090.", "", regex=False)
    )
    edges_df["combined_score"] = edges_df["combined_score"].astype(np.int32)
    log.info("  Toplam fare kenarı (ham): %d", len(edges_df))

    if hedef_doku:
        edges_df = filter_edges_by_tissue(edges_df, tissue_pivot, hedef_doku)

    if edges_df.empty:
        log.warning(
            "⚠ Seçilen doku için pozitif ifadeli fare protein çifti bulunamadı."
        )
        g = ig.Graph(directed=False)
        g.vs["name"] = []
        return g

    log.info("    [3/3] Fare doku ağırlıkları LEFT JOIN ile ekleniyor...")
    edges_df = _attach_tissue_weights(edges_df, tissue_pivot)
    gc.collect()

    # Derece hesapla
    deg_series   = pd.concat(
        [edges_df["protein1"], edges_df["protein2"]]
    ).value_counts()
    edges_df["ka"] = edges_df["protein1"].map(deg_series).fillna(1).astype(np.int32)
    edges_df["kb"] = edges_df["protein2"].map(deg_series).fillna(1).astype(np.int32)

    # Temel ağırlık
    edges_df["w_raw"] = (edges_df["combined_score"] / 1000.0) ** 2
    edges_df["w_hub"] = edges_df["w_raw"]

    # Doku çarpanı
    edges_df["tissue_raw"] = np.sqrt(
        edges_df["exp1"].clip(lower=0.01) * edges_df["exp2"].clip(lower=0.01)
    )

    from ..scientific.tissue import tissue_multiplier
    edges_df["tissue_multiplier"], normalization_metadata = tissue_multiplier(
        edges_df["exp1"], edges_df["exp2"], mode=tissue_normalization_mode,
        db_path=db_path, expression_table="mouse_tissue_expression",
    )

    edges_df["final_weight"]   = (edges_df["w_hub"] * edges_df["tissue_multiplier"]).clip(lower=1e-9)
    edges_df["final_distance"] = strengths_to_distances(edges_df["final_weight"])

    # Düğüm isimleri
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
    g.es["rescue_bridge"] = [False] * len(weights)
    g["tissue_normalization"] = normalization_metadata

    g = g.simplify(combine_edges={"weight": "max", "distance": "min", "rescue_bridge": "first"})
    name_to_idx = {v["name"]: v.index for v in g.vs}

    log.info("    [Mouse Inf Mask] Sonsuz mesafeler maskeleniyor...")
    _mask_inf_distances(g)

    log.info(
        "  Fare ağı tamamlandı: %d düğüm, %d kenar (simplify + inf-mask sonrası).",
        g.vcount(), g.ecount(),
    )

    # Rescue — simplify() SONRASI
    if forced_genes:
        log.info(
            "  [Mouse forced_genes] %d gen için rescue (eşik listesi: %s)…",
            len(forced_genes), [threshold] + RESCUE_THRESHOLDS,
        )
        rescue_query = """
            SELECT protein1, protein2, combined_score
            FROM mouse_interactions
            WHERE (protein1 IN (?, ?) OR protein2 IN (?, ?))
              AND combined_score >= ?
            ORDER BY combined_score DESC LIMIT 50
        """
        con = sqlite3.connect(db_path)
        try:
            for gene in forced_genes:
                raw_id    = str(gene).replace("10090.", "")
                prefix_id = "10090." + raw_id
                alias_id  = DB_FIX_MAP.get(raw_id, raw_id)

                if hedef_doku and raw_id not in tissue_pivot:
                    log.info(
                        "    [Mouse Doku filtresi] Rescue atlandı: %s, %s için "
                        "pozitif ifade kaydı yok.", raw_id, hedef_doku,
                    )
                    continue

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
                        "    ⚠ %s izole kaldı. LCC yapay köprüsü atılıyor...",
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
            "  ⚠ %d yapay köprü kenarı mevcut (rescue_bridge=True).",
            rescue_bridge_count,
        )

    log.info("  → Mouse igraph Graf (final): %d düğüm, %d kenar.", g.vcount(), g.ecount())
    return g


# ===========================================================================
# ADIM 2b: PageRank
# ===========================================================================

def compute_pagerank(g: ig.Graph, damping: float, max_iter: int) -> dict:
    log.info("Mouse Adım 2b: PageRank (d=%.2f, deterministic)…", damping)
    np.random.seed(42)

    try:
        pr_values = g.pagerank(
            damping=damping,
            weights="weight",
            directed=False,
            implementation="prpack"
        )
    except Exception:
        pr_values = g.pagerank(
            damping=damping,
            weights="weight",
            directed=False,
            implementation="arpack"
        )

    pr = {g.vs[i]["name"]: pr_values[i] for i in range(g.vcount())}
    log.info("  → Mouse %d düğüm.", len(pr))
    return pr


# ===========================================================================
# ADIM 3b: Topluluk Tespiti
# ===========================================================================

def run_community_detection(g: ig.Graph, df: pd.DataFrame) -> tuple:
    log.info("Mouse Adım 3b: Louvain (seed=%d)…", LOUVAIN_SEED)

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
    log.info("  → Mouse %d topluluk.", len(communities))
    return df, communities


# ===========================================================================
# GÖREV 2: Graf Düğüm Etiketleme
# ===========================================================================

def tag_graph_nodes(g: ig.Graph, loc_map: dict) -> ig.Graph:
    log.info("Mouse Görev 2: Mekansal etiketleme…")
    tagged = 0
    for v in g.vs:
        loc = loc_map.get(v["name"], "Unknown")
        v["lokalizasyon"] = loc
        if loc != "Unknown":
            tagged += 1
    log.info("  → Mouse %d / %d etiketlendi.", tagged, g.vcount())
    return g


# ===========================================================================
# Congestion Hesabı — Lokal Betweenness
# ===========================================================================

def _calculate_congestion(
    g_pre: ig.Graph,
    g_post: ig.Graph,
    gate_genes: Optional[list] = None,
    max_vertices: int = 3000,
    hop_distance: int = 3,
) -> np.ndarray:
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

    target_set: set = set(seed_indices)
    current_frontier = set(seed_indices)

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
                "  Mouse Congestion: %d düğüm toplandı (%d hop'ta max_vertices aşıldı)",
                len(target_set), hop + 1
            )
            break

    target_list = sorted(target_set)[:max_vertices]
    if not target_list:
        return delta_bc

    log.info(
        "  Mouse Congestion BC: %d hedef düğüm üzerinde hesaplanıyor (hop=%d)…",
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
