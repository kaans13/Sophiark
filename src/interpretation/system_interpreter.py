"""Deterministic system-level interpretation from existing efficiency fields."""

from __future__ import annotations

from typing import Iterable, Mapping
import pandas as pd

from .biological_rules import first_number, fmt_percent
from .sentence_bank import SENTENCES


def compose_system_response(report: pd.DataFrame) -> str:
    """Describe graph metrics without treating them as cellular efficiency."""
    systemic = first_number(report, ("Systemic_Network_Shift_Pct",))
    global_change = first_number(report, ("Global_Efficiency_Change_Pct",))
    mean_local = first_number(report, ("Mean_Local_Efficiency_Change_Pct",))
    target_local = first_number(report, ("Target_Local_Efficiency_Change_Pct", "Local_Efficiency_Kayip_Pct"))
    available = [value for value in (systemic, global_change, mean_local, target_local) if value is not None]
    if not available:
        return "Sistem düzeyi ağ metrikleri bu raporda yok; yorum aday bazlı ΔPageRank verisiyle sınırlıdır."
    text = f"Sistemik ağ kayması {fmt_percent(systemic)}; Grafik verimliliği (global) değişimi {fmt_percent(global_change)}."
    if mean_local is not None:
        text += f" Ortalama yerel grafik verimliliği değişimi {fmt_percent(mean_local)}."
    if target_local is not None:
        text += f" Hedef-yerel grafik verimliliği değişimi {fmt_percent(target_local)}."
    if global_change is not None and target_local is not None and abs(target_local) > abs(global_change) * 1.5:
        text += " " + SENTENCES["system_localized"]
    else:
        text += " " + SENTENCES["system_neutral"]
    return text


def _independent_support(
    candidates: Iterable[Mapping[str, object]], theme_genes: set[str],
    clusters: Iterable[Mapping[str, object]], limit: int = 5,
) -> list[str]:
    """Choose representatives without counting one repeated family repeatedly."""
    membership: dict[str, list[str]] = {}
    for cluster in clusters:
        cluster_id = str(cluster.get("cluster_id") or "")
        for member in cluster.get("members", []):
            membership.setdefault(str(member).upper(), []).append(cluster_id)
    selected: list[str] = []
    represented_clusters: set[str] = set()
    for candidate in candidates:
        gene = str(candidate.get("gene") or "").upper()
        if not gene or gene not in theme_genes:
            continue
        groups = membership.get(gene, [])
        if groups and any(group in represented_clusters for group in groups):
            continue
        selected.append(gene)
        represented_clusters.update(groups)
        if len(selected) >= limit:
            break
    return selected


def compose_analysis_summary(
    *, target: str, tissue: str, themes: list[Mapping[str, object]],
    clusters: list[Mapping[str, object]], candidates: list[Mapping[str, object]],
    loss_themes: list[Mapping[str, object]] | None = None,
    loss_candidates: list[Mapping[str, object]] | None = None,
) -> dict[str, object]:
    """Compose a cohesive, source-bounded primary/secondary interpretation."""
    axes = [theme for theme in themes if int(theme.get("gene_count") or 0) > 0]
    if not axes:
        return {
            "primary": "Mevcut aday anotasyonları baskın bir işlevsel eksen oluşturmak için yetersizdir; bu nedenle yorum perturbasyon sonrası topolojik değişimle sınırlıdır.",
            "secondary": "", "supporting_candidates": [], "family_warning": "",
            "validation": [SENTENCES["validation_transcript_protein"], SENTENCES["validation_perturbation"]],
            "primary_theme": None, "secondary_theme": None,
        }
    primary_theme = axes[0]
    primary_genes = {str(gene).upper() for gene in primary_theme.get("genes", [])}
    support = _independent_support(candidates, primary_genes, clusters)
    primary = SENTENCES["summary_primary"].format(target=target or "Seçilen hedef", theme=primary_theme["theme"])
    if support:
        primary += " " + SENTENCES["summary_support"].format(genes=", ".join(support), theme=primary_theme["theme"])
    secondary_theme = next((theme for theme in axes[1:] if theme.get("theme_id") != "membrane_transport"), None)
    secondary = ""
    if secondary_theme:
        secondary_support = _independent_support(
            candidates, {str(gene).upper() for gene in secondary_theme.get("genes", [])}, clusters, limit=4,
        )
        if secondary_support:
            secondary = SENTENCES["summary_secondary"].format(
                theme=secondary_theme["theme"], genes=", ".join(secondary_support),
            )
    family_candidates = [
        cluster for cluster in clusters
        if cluster.get("possible_topology_amplification") and set(str(item).upper() for item in cluster.get("members", [])) & primary_genes
    ]
    family_candidates.sort(
        key=lambda cluster: (
            -len(set(str(item).upper() for item in cluster.get("members", [])) & primary_genes),
            -float(cluster.get("metrics", {}).get("mean_delta_pagerank_pct") or 0.0),
            str(cluster.get("name") or ""),
        )
    )
    family_names = [str(cluster["name"]) for cluster in family_candidates[:2]]
    family_warning = SENTENCES["summary_family_caveat"].format(families=", ".join(family_names)) if family_names else ""
    tissue_note = SENTENCES["summary_tissue"].format(tissue=tissue) if tissue and tissue.casefold() not in {"none", "nan"} else ""
    validation = [
        SENTENCES["validation_theme"].format(theme=primary_theme["theme"]),
        SENTENCES["validation_transcript_protein"],
        SENTENCES["validation_perturbation"],
    ]
    loss_themes = list(loss_themes or [])
    loss_candidates = list(loss_candidates or [])
    loss_summary = ""
    if loss_candidates:
        strongest = loss_candidates[0]
        loss_axes = [item for item in loss_themes if int(item.get("gene_count") or 0) > 0]
        if loss_axes:
            loss_summary = (
                f"Negatif yeniden dağılım tarafında en güçlü göreli topolojik pay kaybı "
                f"{strongest.get('gene')} üzerinde {fmt_percent(float(strongest.get('delta_pagerank_pct')), signed=True)} ile görüldü; "
                f"baskın anotasyon ekseni {loss_axes[0].get('theme')} olarak belirdi. "
                "Bu bulgu ekspresyon, aktivite veya biyolojik işlev kaybı anlamına gelmez."
            )
        else:
            loss_summary = (
                f"Negatif tarafta {strongest.get('gene')} için {fmt_percent(float(strongest.get('delta_pagerank_pct')), signed=True)} "
                "göreli topolojik pay kaybı görüldü; mevcut anotasyonlar baskın bir işlevsel eksen oluşturmak için yeterli değildir."
            )
    return {
        "primary": primary, "secondary": secondary, "supporting_candidates": support,
        "family_warning": family_warning, "tissue_note": tissue_note, "validation": validation,
        "primary_theme": dict(primary_theme), "secondary_theme": dict(secondary_theme) if secondary_theme else None,
        "loss_summary": loss_summary,
    }
