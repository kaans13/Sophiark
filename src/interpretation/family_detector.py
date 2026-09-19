"""Deterministic family, paralog and complex-context detection.

This module is a read-only consumer of result symbols.  It produces explicit,
auditable clusters and does not feed any score back into Sophiark's engine.
"""

from __future__ import annotations

from collections import Counter
import re
from typing import Iterable, Mapping, Sequence

import pandas as pd

from .family_reference import complex_portal_available, hgnc_groups, hgnc_groups_by_symbol
from .models import InterpretationCluster


_SYMBOL_COLUMN_CANDIDATES = ("Symbol", "gene", "Gene")
_DELTA_COLUMNS = ("Delta_PageRank_Pct", "delta_pagerank_pct")
_HINTERLAND_COLUMNS = ("Hinterland_Skoru", "hinterland_baseline")
_COMMUNITY_COLUMNS = ("Community_Context", "Community", "community")

# These patterns cover HGNC omissions caused by unstable duplicated-locus
# nomenclature. They are deliberately narrow; a shared prefix is never a
# general family assertion.  METTL is intentionally absent.
_NARROW_HEURISTICS: tuple[tuple[str, re.Pattern[str]], ...] = (
    ("NOTCH2NL paralogları", re.compile(r"^NOTCH2NL[ABC]$", re.I)),
)


def _symbols(values: Iterable[str]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(
        value.strip().upper() for value in values if isinstance(value, str) and value.strip()
    ))


def _cluster_type(group_name: str) -> str:
    lowered = group_name.casefold()
    if "complex" in lowered or "complexes" in lowered:
        return "protein_complex"
    if "subfamily" in lowered or "subunits" in lowered:
        return "subfamily"
    if "members" in lowered or "paralog" in lowered:
        return "paralog"
    return "family"


def _number(row: pd.Series, names: Sequence[str]) -> float | None:
    for name in names:
        if name in row.index:
            value = pd.to_numeric(pd.Series([row[name]]), errors="coerce").iloc[0]
            if pd.notna(value):
                return float(value)
    return None


def _frame_symbol_column(frame: pd.DataFrame) -> str | None:
    return next((column for column in _SYMBOL_COLUMN_CANDIDATES if column in frame.columns), None)


def _metrics(
    members: tuple[str, ...],
    candidates: pd.DataFrame | None,
    direct_neighbors: Iterable[str] | None,
) -> dict[str, float | int | None]:
    """Summarize a cluster without touching its candidate values."""
    base: dict[str, float | int | None] = {
        "member_count": len(members), "top100_fraction": None,
        "mean_delta_pagerank_pct": None, "median_delta_pagerank_pct": None,
        "mean_hinterland": None, "median_hinterland": None,
        "direct_target_neighbor_fraction": None, "same_community_fraction": None,
    }
    if candidates is None or candidates.empty:
        return base
    symbol_column = _frame_symbol_column(candidates)
    if symbol_column is None:
        return base
    frame = candidates.copy(deep=True)
    frame["_symbol"] = frame[symbol_column].fillna("").astype(str).str.strip().str.upper()
    for delta_column in _DELTA_COLUMNS:
        if delta_column in frame.columns:
            frame["_delta"] = pd.to_numeric(frame[delta_column], errors="coerce")
            frame = frame.sort_values("_delta", ascending=False, kind="mergesort")
            break
    member_frame = frame[frame["_symbol"].isin(members)].copy()
    top_symbols = set(frame.head(100)["_symbol"])
    base["top100_fraction"] = len(set(members) & top_symbols) / 100.0
    if "_delta" in member_frame:
        base["mean_delta_pagerank_pct"] = float(member_frame["_delta"].mean()) if member_frame["_delta"].notna().any() else None
        base["median_delta_pagerank_pct"] = float(member_frame["_delta"].median()) if member_frame["_delta"].notna().any() else None
    hinterland_column = next((column for column in _HINTERLAND_COLUMNS if column in member_frame.columns), None)
    if hinterland_column:
        values = pd.to_numeric(member_frame[hinterland_column], errors="coerce")
        base["mean_hinterland"] = float(values.mean()) if values.notna().any() else None
        base["median_hinterland"] = float(values.median()) if values.notna().any() else None
    if direct_neighbors is not None:
        neighbors = set(_symbols(direct_neighbors))
        base["direct_target_neighbor_fraction"] = len(set(members) & neighbors) / len(members)
    community_column = next((column for column in _COMMUNITY_COLUMNS if column in member_frame.columns), None)
    if community_column and len(member_frame) > 0:
        communities = [str(value) for value in member_frame[community_column] if pd.notna(value) and str(value).strip()]
        if communities:
            base["same_community_fraction"] = Counter(communities).most_common(1)[0][1] / len(members)
    return base


def _risk(metrics: Mapping[str, float | int | None]) -> str:
    count = int(metrics["member_count"] or 0)
    top_fraction = float(metrics["top100_fraction"] or 0.0)
    same_community = float(metrics["same_community_fraction"] or 0.0)
    if count >= 5 or top_fraction >= 0.15 or (count >= 3 and same_community >= 0.8):
        return "yüksek"
    if count >= 3 or top_fraction >= 0.05:
        return "orta"
    return "düşük"


def _curated_cluster(group_id: str, group: Mapping[str, object], members: tuple[str, ...],
                     candidates: pd.DataFrame | None, direct_neighbors: Iterable[str] | None) -> InterpretationCluster:
    metrics = _metrics(members, candidates, direct_neighbors)
    risk = _risk(metrics)
    group_name = str(group["name"])
    source = "HGNC Gene Group"
    if _cluster_type(group_name) == "protein_complex" and complex_portal_available():
        source += "; Complex Portal referansı yerelde mevcut"
    return InterpretationCluster(
        cluster_id=f"hgnc:{group_id}", cluster_type=_cluster_type(group_name), name=group_name,
        members=members, member_count=len(members), confidence_score=0.97,
        confidence_label="yüksek", evidence=f"HGNC Group ID {group_id}", source=source,
        heuristic_used=False, possible_topology_amplification=risk != "düşük",
        warning=("Aynı aile/kompleks üyeleri bağımsız biyolojik kanıtlar gibi toplanmamalıdır."
                 if risk != "düşük" else None), metrics=metrics,
        topology_amplification_risk=risk,
    )


def _heuristic_clusters(symbols: tuple[str, ...], candidates: pd.DataFrame | None,
                        direct_neighbors: Iterable[str] | None) -> list[InterpretationCluster]:
    clusters: list[InterpretationCluster] = []
    for name, pattern in _NARROW_HEURISTICS:
        members = tuple(symbol for symbol in symbols if pattern.fullmatch(symbol))
        if len(members) < 2:
            continue
        metrics = _metrics(members, candidates, direct_neighbors)
        risk = _risk(metrics)
        clusters.append(InterpretationCluster(
            cluster_id=f"heuristic:{name.casefold().replace(' ', '-')}", cluster_type="heuristic", name=name,
            members=members, member_count=len(members), confidence_score=0.35,
            confidence_label="düşük", evidence="Daraltılmış sembol deseni; küratörlü grup kaydı bulunamadı.",
            source="Sembol sezgiseli", heuristic_used=True,
            possible_topology_amplification=risk != "düşük",
            warning="Bu grup sembol deseniyle bulundu; aile/kompleks kanıtı olarak bağımsız doğrulama gerektirir.",
            metrics=metrics, topology_amplification_risk=risk,
        ))
    return clusters


def _calcium_functional_system(symbols: tuple[str, ...], candidates: pd.DataFrame | None,
                               direct_neighbors: Iterable[str] | None) -> list[InterpretationCluster]:
    """Keep CACNA/CACNB/CACNG separate while optionally reporting their system."""
    classes = {"CACNA": [], "CACNB": [], "CACNG": [], "CACNA2D": []}
    for symbol in symbols:
        for prefix in classes:
            if symbol.startswith(prefix):
                classes[prefix].append(symbol)
                break
    present_classes = sum(bool(members) for members in classes.values())
    members = tuple(member for group in classes.values() for member in group)
    if present_classes < 2 or len(members) < 2:
        return []
    metrics = _metrics(members, candidates, direct_neighbors)
    risk = _risk(metrics)
    return [InterpretationCluster(
        cluster_id="functional-system:voltage-gated-calcium-channel-subunits",
        cluster_type="functional_system", name="Voltaj kapılı kalsiyum kanalı alt birimleri",
        members=members, member_count=len(members), confidence_score=0.88, confidence_label="yüksek",
        evidence="CACNA, CACNB, CACNG ve CACNA2D sınıfları aynı işlevsel sistem altında; aynı aile değildir.",
        source="Küratörlü sınıf ayrımı", heuristic_used=False,
        possible_topology_amplification=risk != "düşük",
        warning="Bu kayıt üst işlevsel sistemdir; alt birim sınıflarını tek paralog ailesi gibi birleştirmez.",
        metrics=metrics, topology_amplification_risk=risk,
    )]


def detect_family_complex_clusters(
    symbols: Iterable[str], *, candidates: pd.DataFrame | None = None,
    direct_neighbors: Iterable[str] | None = None,
) -> list[dict[str, object]]:
    """Return transparent, deterministic clusters for the supplied symbols.

    `symbols` may be all network candidates or only a selected subset. The
    function is intentionally side-effect-free and has no scientific-engine
    imports.
    """
    normalized = _symbols(symbols)
    groups = hgnc_groups()
    by_symbol = hgnc_groups_by_symbol()
    group_ids = sorted({group_id for symbol in normalized for group_id in by_symbol.get(symbol, ())})
    clusters: list[InterpretationCluster] = []
    for group_id in group_ids:
        group = groups[group_id]
        members = tuple(symbol for symbol in normalized if symbol in group["members"])
        if len(members) >= 2:
            clusters.append(_curated_cluster(group_id, group, members, candidates, direct_neighbors))
    clusters.extend(_calcium_functional_system(normalized, candidates, direct_neighbors))
    clusters.extend(_heuristic_clusters(normalized, candidates, direct_neighbors))
    # A stable order keeps exported and rendered output deterministic.
    clusters.sort(key=lambda cluster: (-cluster.confidence_score, -cluster.member_count, cluster.name, cluster.cluster_id))
    return [cluster.to_dict() for cluster in clusters]
