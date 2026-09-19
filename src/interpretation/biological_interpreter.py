"""Structured, read-only interpretation of an existing forward-analysis report."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from typing import Mapping, Sequence

import pandas as pd

from .biological_rules import LIMITATIONS, directional_effect, first_number, first_text, fmt_percent
from .biological_facts import extract_biological_facts
from .candidate_interpreter import compose_candidate_interpretation
from .system_interpreter import compose_analysis_summary, compose_system_response
from .theme_detection import annotation_text, detect_family_patterns, detect_functional_themes
from .theme_detector_v2 import detect_themes_from_facts


@dataclass(frozen=True)
class BiologicalInterpretation:
    overview: str
    system_response: str
    analysis_summary: dict[str, object]
    top_candidates: list[dict[str, object]]
    loss_candidates: list[dict[str, object]]
    baseline_vs_perturbation: dict[str, list[dict[str, object]]]
    functional_themes: list[dict[str, object]]
    loss_functional_themes: list[dict[str, object]]
    family_patterns: list[dict[str, object]]
    family_clusters: list[dict[str, object]]
    loss_family_clusters: list[dict[str, object]]
    redistribution_balance: dict[str, object]
    canonical_component_check: list[dict[str, object]]
    directed_evidence: list[dict[str, object]]
    enrichment: dict[str, object]
    statistical_validation: dict[str, object] | None
    debug: dict[str, object]
    limitations: str = LIMITATIONS

    def to_dict(self) -> dict[str, object]:
        return asdict(self)


def _symbol(row: pd.Series) -> str:
    value = row.get("Symbol")
    return str(value).strip() if pd.notna(value) and str(value).strip() else "—"


def _numeric(row: pd.Series, *names: str) -> float | None:
    for name in names:
        value = pd.to_numeric(pd.Series([row.get(name)]), errors="coerce").iloc[0]
        if pd.notna(value):
            return float(value)
    return None


def _candidate_sentence(row: pd.Series, hinterland_median: float | None) -> str:
    symbol = _symbol(row)
    delta = _numeric(row, "Delta_PageRank_Pct")
    hinterland = _numeric(row, "Hinterland_Skoru")
    direction = "pozitif göreli PageRank yeniden dağılımı" if (delta or 0) > 0 else "göreli PageRank kaybı"
    sentence = f"{symbol}, {fmt_percent(delta, signed=True)} ile {direction} gösterdi."
    if hinterland is not None and hinterland_median is not None:
        if hinterland < hinterland_median:
            sentence += " Başlangıç Ağdaki Önemi değeri aday kümesinin medyanının altında olduğundan, yanıt başlangıçta baskın olmayan bir ağ konumunda görülüyor."
        else:
            sentence += " Ağdaki Önemi başlangıç topolojik bağlamıdır; bu değer perturbasyon sonucunun kendisi değildir."
    annotation = annotation_text(row)
    if annotation:
        sentence += " Mevcut anotasyonlar bu aday için ek işlevsel bağlam sağlar; mekanizma veya aktivasyon kanıtı oluşturmaz."
    return sentence


def _build_overview(report: pd.DataFrame, *, targets: Sequence[str], gene_to_symbol: Mapping[str, str], tissue: str,
                    organism: str, edge_attenuation: float | None, network_nodes: int | None,
                    network_edges: int | None) -> str:
    labels = ", ".join(gene_to_symbol.get(target, target) for target in targets) or "seçilen lokalizasyon"
    parts = [f"{labels}, {tissue}-özgül ağırlıklandırılmış {organism} PPI ağında analiz edildi."]
    if edge_attenuation is not None:
        parts.append(f"Hedefe komşu kenarlar başlangıç ağırlıklarının yaklaşık %{edge_attenuation * 100:.3g}'üne attenüe edildi ve göreli PageRank yeniden dağılımı yeniden hesaplandı.")
    else:
        parts.append("Hedefe komşu kenarlar attenüe edildi ve göreli PageRank yeniden dağılımı yeniden hesaplandı.")
    if network_nodes is not None and network_edges is not None:
        parts.append(f"Ağ bağlamı: {network_nodes:,} düğüm ve {network_edges:,} kenar.")
    severity = first_text(report, ("Perturbation_Severity",))
    if severity:
        parts.append(f"Perturbasyon şiddeti kaydı: {severity}.")
    mode = first_text(report, ("Redistribution_Selection_Mode",))
    if mode:
        parts.append(f"Aday seçimi modu: {mode}.")
    return " ".join(parts)


def _build_system_response(report: pd.DataFrame) -> str:
    systemic = first_number(report, ("Systemic_Network_Shift_Pct",))
    global_change = first_number(report, ("Global_Efficiency_Change_Pct",))
    mean_local = first_number(report, ("Mean_Local_Efficiency_Change_Pct",))
    target_local = first_number(report, ("Target_Local_Efficiency_Change_Pct", "Local_Efficiency_Kayip_Pct"))
    available = [value for value in (systemic, global_change, mean_local, target_local) if value is not None]
    if not available:
        return "Sistem düzeyi verimlilik metrikleri bu raporda yok; yorum yalnızca aday bazlı PageRank yeniden dağılımıyla sınırlıdır."
    sentence = f"Sistemik ağ kayması {fmt_percent(systemic)}; global grafik verimliliği değişimi {fmt_percent(global_change)}."
    if mean_local is not None:
        sentence += f" Ortalama yerel grafik verimliliği değişimi {fmt_percent(mean_local)}."
    if target_local is not None:
        sentence += f" Hedef-yerel grafik verimliliği değişimi {fmt_percent(target_local)}."
    if global_change is not None and target_local is not None and abs(target_local) > abs(global_change) * 1.5:
        sentence += " Hedef-yerel değişimin global değişimden daha büyük olması, yapısal etkinin ağın tamamına eşdeğer bir çöküşten çok hedef çevresinde yoğunlaştığıyla uyumludur."
    sentence += " Grafik verimliliği, hücresel veya metabolik verimlilik ölçümü değildir."
    return sentence


def select_redistribution_losses(signed_result: pd.DataFrame | None, *, min_abs_delta_pct: float = .05,
                                 limit: int | None = 100) -> pd.DataFrame:
    """Select already-computed negative redistribution rows; never recompute PageRank."""
    if not isinstance(signed_result, pd.DataFrame) or "Delta_PageRank_Pct" not in signed_result:
        return pd.DataFrame()
    threshold = abs(float(min_abs_delta_pct))
    losses = signed_result.copy(deep=True)
    losses["_delta"] = pd.to_numeric(losses["Delta_PageRank_Pct"], errors="coerce")
    losses = losses.loc[losses["_delta"] < -threshold].sort_values("_delta", ascending=True, kind="mergesort")
    if limit is not None:
        losses = losses.head(max(0, int(limit)))
    return losses.drop(columns=["_delta"])


def _candidate_records(candidates: pd.DataFrame, *, limit: int = 15, direction: str = "gain") -> tuple[list[dict[str, object]], pd.DataFrame]:
    usable = candidates.copy(deep=True)
    if "Delta_PageRank_Pct" not in usable.columns:
        return [], usable.iloc[0:0]
    usable["_delta"] = pd.to_numeric(usable["Delta_PageRank_Pct"], errors="coerce")
    usable = usable.dropna(subset=["_delta"])
    usable = usable[usable["_delta"] < 0] if direction == "loss" else usable[usable["_delta"] > 0]
    usable = usable.sort_values("_delta", ascending=(direction == "loss"), kind="mergesort")
    top = usable.head(limit).copy()
    median = pd.to_numeric(usable.get("Hinterland_Skoru"), errors="coerce").median() if "Hinterland_Skoru" in usable else None
    records: list[dict[str, object]] = []
    for _, row in top.iterrows():
        raw_id = row.get("gene")
        entity_id = str(raw_id).strip() if pd.notna(raw_id) and str(raw_id).strip() else "—"
        records.append({
            "gene": _symbol(row), "entity_id": entity_id,
            "delta_pagerank_pct": _numeric(row, "Delta_PageRank_Pct"),
            "pagerank_baseline": _numeric(row, "PageRank_Baseline"),
            "pagerank_perturbed": _numeric(row, "PageRank_Perturbed"),
            "hinterland_baseline": _numeric(row, "Hinterland_Skoru"), "bc_baseline": _numeric(row, "BC_Skoru"),
            "localization": str(row.get("Lokalizasyon") or "veri yok"),
            "bridge_label": str(row.get("Gümrük_Kapisi") or "veri yok"),
            "drug_score": _numeric(row, "Drug_Score"),
            "annotation": str(row.get("GO Biyolojik Süreç") or row.get("GO Moleküler İşlev") or row.get("Protein_Adi") or "veri yok"),
            "external_evidence": str(row.get("Düzenleyici_TFler") or row.get("Essentiality") or "veri yok"),
            "interpretation": _candidate_sentence(row, median),
        })
    return records, usable


def _theme_maps(themes: list[dict[str, object]]) -> dict[str, list[str]]:
    mapped: dict[str, list[str]] = {}
    for theme in themes:
        for symbol in theme["genes"]:
            mapped.setdefault(str(symbol), []).append(str(theme["theme"]))
    return mapped


def _annotate_candidates(*, records: list[dict[str, object]], facts: list[dict[str, object]], themes: list[dict[str, object]],
                         clusters: list[dict[str, object]]) -> None:
    facts_by_symbol = {str(item.get("symbol") or ""): item for item in facts if isinstance(item, dict)}
    theme_names = _theme_maps(themes)
    cluster_names: dict[str, list[str]] = {}
    for cluster in clusters:
        for symbol in cluster.get("members", []):
            cluster_names.setdefault(str(symbol), []).append(str(cluster["name"]))
    primary_genes = tuple(themes[0]["genes"]) if themes else ()
    primary_set = set(map(str.upper, primary_genes))
    for candidate in records:
        symbol = str(candidate["gene"]).upper()
        fact = facts_by_symbol.get(symbol, {
            "symbol": candidate["gene"], "delta_pagerank_pct": candidate["delta_pagerank_pct"],
            "hinterland_baseline": candidate["hinterland_baseline"], "bc_baseline": candidate["bc_baseline"], "annotations": {},
        })
        composed = compose_candidate_interpretation(
            fact, theme_labels=theme_names.get(symbol, ()), cluster_names=cluster_names.get(symbol, ()),
            theme_peer_symbols=primary_genes if symbol in primary_set else (),
        )
        candidate.update({
            "interpretation": composed["text"], "interpretation_confidence_score": composed["confidence_score"],
            "interpretation_confidence_label": composed["confidence_label"], "validation_suggestions": composed["validation_suggestions"],
        })


def _redistribution_balance(gains: list[dict[str, object]], losses: list[dict[str, object]],
                            gain_themes: list[dict[str, object]], loss_themes: list[dict[str, object]]) -> dict[str, object]:
    gain_names = [str(item["theme"]) for item in gain_themes]
    loss_names = [str(item["theme"]) for item in loss_themes]
    shared = sorted(set(gain_names) & set(loss_names))
    return {
        "positive_candidates": len(gains), "negative_candidates": len(losses),
        "strongest_gain": gains[0] if gains else None, "strongest_loss": losses[0] if losses else None,
        "dominant_gain_themes": gain_names[:3], "dominant_loss_themes": loss_names[:3], "shared_themes": shared,
        "consistency_warning": "Aynı aday hem kazanç hem kayıp listesinde yer alıyor; işaret tutarlılığı incelenmeli."
        if set(str(item["gene"]).upper() for item in gains) & set(str(item["gene"]).upper() for item in losses) else "",
    }


def _canonical_component_check(components: Sequence[str] | None, gains: list[dict[str, object]],
                               losses: list[dict[str, object]]) -> list[dict[str, object]]:
    """Check supplied canonical components against both signed sides before any absence note."""
    if not components:
        return []
    gain_by_gene = {str(item["gene"]).upper(): item for item in gains}
    loss_by_gene = {str(item["gene"]).upper(): item for item in losses}
    records: list[dict[str, object]] = []
    for component in components:
        key = str(component).upper()
        if key in gain_by_gene:
            records.append({"gene": component, "side": "pozitif", "delta_pagerank_pct": gain_by_gene[key].get("delta_pagerank_pct"), "message": "Canonical bileşen göreli topolojik pay artışı tarafında yer alıyor."})
        elif key in loss_by_gene:
            records.append({"gene": component, "side": "negatif", "delta_pagerank_pct": loss_by_gene[key].get("delta_pagerank_pct"), "message": "Canonical bileşen göreli topolojik pay kaybı tarafında yer alıyor."})
        else:
            records.append({"gene": component, "side": "eşik dışı", "delta_pagerank_pct": None, "message": "Canonical bileşen mevcut signed sonuçta seçili eşik aralığında değil."})
    return records


def _baseline_groups(usable: pd.DataFrame) -> dict[str, list[dict[str, object]]]:
    if usable.empty or "Hinterland_Skoru" not in usable:
        return {"perturbation_sensitive_baseline_nondominant": [], "baseline_central_weakly_responsive": []}
    work = usable.copy(deep=True)
    work["_hinterland"] = pd.to_numeric(work["Hinterland_Skoru"], errors="coerce")
    work["_abs_delta"] = work["_delta"].abs()
    h_low, h_high = work["_hinterland"].quantile([.50, .75]).tolist()
    delta_high = work["_delta"].quantile(.75)
    abs_delta_low = work["_abs_delta"].quantile(.50)
    sensitive = work[(work["_delta"] >= delta_high) & (work["_hinterland"] <= h_low)].head(10)
    central = work[(work["_hinterland"] >= h_high) & (work["_abs_delta"] <= abs_delta_low)].sort_values("_hinterland", ascending=False).head(10)

    def records(frame: pd.DataFrame) -> list[dict[str, object]]:
        return [{
            "gene": _symbol(row),
            "entity_id": str(row.get("gene")).strip() if pd.notna(row.get("gene")) and str(row.get("gene")).strip() else "—",
            "delta_pagerank_pct": _numeric(row, "Delta_PageRank_Pct"),
            "hinterland_baseline": _numeric(row, "Hinterland_Skoru"),
        } for _, row in frame.iterrows()]
    return {
        "perturbation_sensitive_baseline_nondominant": records(sensitive),
        "baseline_central_weakly_responsive": records(central),
    }


def _directed_records(evidence: pd.DataFrame | None, relevant_symbols: set[str]) -> list[dict[str, object]]:
    if evidence is None or evidence.empty or not {"source", "target"}.issubset(evidence.columns):
        return []
    frame = evidence.copy(deep=True)
    source = frame["source"].fillna("").astype(str)
    target = frame["target"].fillna("").astype(str)
    directed = frame.get("consensus_direction", pd.Series(False, index=frame.index)).fillna(False).astype(bool)
    frame = frame[directed & (source.isin(relevant_symbols) | target.isin(relevant_symbols))]
    records: list[dict[str, object]] = []
    dedupe_columns = [column for column in ("source", "target", "evidence_status", "provenance") if column in frame.columns]
    for _, row in frame.drop_duplicates(subset=dedupe_columns).head(50).iterrows():
        records.append({
            "source": str(row["source"]), "target": str(row["target"]), "effect": directional_effect(row),
            "provenance": str(row.get("provenance") or "veri kaynağı belirtilmemiş"),
            "references": str(row.get("references") or "referans alanı yok"),
        })
    return records


def _enrichment_context(enrichment: pd.DataFrame | None) -> dict[str, object]:
    if enrichment is None or enrichment.empty:
        return {"available": False, "message": "Bu analiz için gösterilecek GO/KEGG zenginleştirme kaydı yok."}
    frame = enrichment.copy(deep=True)
    adjusted_column = next((column for column in ("Adjusted P-value", "adjusted_p_value", "FDR", "q_value") if column in frame), None)
    if adjusted_column:
        frame["_adjusted"] = pd.to_numeric(frame[adjusted_column], errors="coerce")
        frame = frame.sort_values("_adjusted", ascending=True, na_position="last")
        significant = frame[frame["_adjusted"] <= .05]
    else:
        significant = frame.iloc[0:0]
    term_column = next((column for column in ("Term", "term", "Yolak") if column in frame), None)
    genes_column = next((column for column in ("Genes", "genes", "Genler") if column in frame), None)
    rows = []
    for _, row in frame.head(20).iterrows():
        rows.append({"term": str(row.get(term_column) or "terim yok"), "adjusted_p": _numeric(row, adjusted_column) if adjusted_column else None, "genes": str(row.get(genes_column) or "gen alanı yok")})
    message = "Terimler, yeniden dağılım adayları arasında aşırı temsil edilen anotasyon bağlamını gösterir; yolak aktivasyonu veya inhibisyonu anlamına gelmez."
    if adjusted_column and significant.empty:
        message = "Çoklu test düzeltmesinden sonra hiçbir terim anlamlı kalmadı."
    return {"available": True, "message": message, "rows": rows, "adjusted_column": adjusted_column}


def _statistical_validation(report: pd.DataFrame, candidates: pd.DataFrame) -> dict[str, object] | None:
    modes = set(report.get("Redistribution_Selection_Mode", pd.Series(dtype=str)).dropna().astype(str))
    if modes != {"null_fdr"}:
        return None
    columns = [column for column in ("Symbol", "gene", "empirical_p", "q_value", "significant_redistribution") if column in candidates]
    if not columns:
        return {"message": "null_fdr seçilmiş ancak gösterilecek istatistik alanı yok.", "rows": []}
    rows = candidates[columns].head(25).copy(deep=True).to_dict("records")
    return {"message": "Bu bölüm yalnızca null-model/FDR seçimi gerçekten kullanıldığında gösterilir.", "rows": rows}


def interpret_forward_results(*, report: pd.DataFrame, candidates: pd.DataFrame, targets: Sequence[str],
                              gene_to_symbol: Mapping[str, str], tissue: str, organism: str,
                              edge_attenuation: float | None = None, network_nodes: int | None = None,
                              network_edges: int | None = None, directed_evidence: pd.DataFrame | None = None,
                              enrichment: pd.DataFrame | None = None,
                              loss_candidates: pd.DataFrame | None = None,
                              canonical_components: Sequence[str] | None = None) -> BiologicalInterpretation:
    """Create a deterministic interpretation without mutating the supplied outputs."""
    top_candidates, usable = _candidate_records(candidates)
    facts_bundle = extract_biological_facts(usable.head(100), organism=organism)
    themes_v2 = detect_themes_from_facts(facts_bundle["candidates"])
    clusters = list(facts_bundle["clusters"])
    _annotate_candidates(records=top_candidates, facts=facts_bundle["candidates"], themes=themes_v2, clusters=clusters)
    top_losses, usable_losses = _candidate_records(loss_candidates if isinstance(loss_candidates, pd.DataFrame) else pd.DataFrame(), direction="loss")
    loss_facts_bundle = extract_biological_facts(usable_losses.head(100), organism=organism)
    loss_themes = detect_themes_from_facts(loss_facts_bundle["candidates"], direction="loss")
    loss_clusters = list(loss_facts_bundle["clusters"])
    _annotate_candidates(records=top_losses, facts=loss_facts_bundle["candidates"], themes=loss_themes, clusters=loss_clusters)
    relevant = set(gene_to_symbol.get(target, target) for target in targets)
    relevant.update(item["gene"] for item in top_candidates)
    target_labels = [gene_to_symbol.get(target, target) for target in targets]
    analysis_summary = compose_analysis_summary(
        target=", ".join(target_labels), tissue=tissue, themes=themes_v2,
        clusters=clusters, candidates=top_candidates, loss_themes=loss_themes, loss_candidates=top_losses,
    )
    return BiologicalInterpretation(
        overview=_build_overview(report, targets=targets, gene_to_symbol=gene_to_symbol, tissue=tissue, organism=organism, edge_attenuation=edge_attenuation, network_nodes=network_nodes, network_edges=network_edges),
        system_response=compose_system_response(report),
        analysis_summary=analysis_summary,
        top_candidates=top_candidates,
        loss_candidates=top_losses,
        baseline_vs_perturbation=_baseline_groups(usable),
        functional_themes=[{
            **theme, "count": theme["gene_count"],
        } for theme in themes_v2],
        loss_functional_themes=[{**theme, "count": theme["gene_count"]} for theme in loss_themes],
        family_patterns=detect_family_patterns(usable.head(50)),
        family_clusters=clusters,
        loss_family_clusters=loss_clusters,
        redistribution_balance=_redistribution_balance(top_candidates, top_losses, themes_v2, loss_themes),
        canonical_component_check=_canonical_component_check(canonical_components, top_candidates, top_losses),
        directed_evidence=_directed_records(directed_evidence, relevant),
        enrichment=_enrichment_context(enrichment),
        statistical_validation=_statistical_validation(report, candidates),
        debug={
            **facts_bundle.get("debug", {}),
            "theme_scores": [{"theme": item["theme"], "priority": item["interpretation_priority"], "genes": item["genes"]} for item in themes_v2],
            "loss_theme_scores": [{"theme": item["theme"], "priority": item["interpretation_priority"], "genes": item["genes"]} for item in loss_themes],
            "family_clusters": clusters,
            "selected_primary_theme": analysis_summary.get("primary_theme", {}).get("theme") if analysis_summary.get("primary_theme") else None,
            "selected_secondary_theme": analysis_summary.get("secondary_theme", {}).get("theme") if analysis_summary.get("secondary_theme") else None,
            "sentence_rules": ["summary_primary", "summary_support", "summary_secondary", "summary_family_caveat", "candidate_interpretation"],
        },
    )
