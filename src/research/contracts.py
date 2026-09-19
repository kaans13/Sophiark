"""Machine-readable interpretation contracts for existing Sophiark fields.

This registry centralizes presentation semantics without importing Streamlit or
performing any scientific calculation.  Labels and numeric formatting for
existing UI fields mirror ``src.ui.metric_presentation``; limitations add the
scientific boundaries required by the Research Context specification.
"""

from __future__ import annotations

import math
from types import MappingProxyType
from typing import Any

import numpy as np
import pandas as pd

from .models import InterpretationContract, MetricKind, MetricOrigin


_NO_CAUSALITY = (
    "gene expression changed",
    "protein abundance changed",
    "activation or inhibition",
    "causality",
    "experimental validation",
)
_TOPOLOGY_ONLY = (
    "biological necessity",
    "therapeutic suitability",
    "causality",
    "experimental validation",
)


def _contract(
    metric_id: str,
    source_fields: tuple[str, ...],
    display_name: str,
    definition: str,
    *,
    origin: MetricOrigin = MetricOrigin.SOPHIARK_COMPUTED,
    unit: str | None = None,
    precision: int | None = None,
    signed: bool = False,
    scope: str = "network_context",
    supports: tuple[str, ...] = (),
    does_not_support: tuple[str, ...] = _TOPOLOGY_ONLY,
    kind: MetricKind = MetricKind.NUMBER,
    suffix: str = "",
) -> InterpretationContract:
    return InterpretationContract(
        metric_id=metric_id,
        source_fields=source_fields,
        display_name=display_name,
        origin=origin,
        unit=unit,
        precision=precision,
        signed=signed,
        scope=scope,
        definition=definition,
        supports=supports,
        does_not_support=does_not_support,
        kind=kind,
        suffix=suffix,
    )


_CONTRACT_LIST = (
    _contract(
        "delta_pagerank_pct", ("Delta_PageRank_Pct",), "Ağ önemi değişimi (%)",
        "Relative PageRank change after target attenuation in the current network model.",
        unit="percent", precision=4, signed=True, suffix="%", scope="perturbation_response",
        supports=("relative network importance changed",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "delta_pagerank", ("Delta_PageRank",), "Ağ önemi değişimi",
        "Absolute PageRank difference between perturbed and baseline model states.",
        precision=8, signed=True, scope="perturbation_response",
        supports=("absolute PageRank share changed",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "pagerank_baseline", ("PageRank_Baseline", "pagerank_raw"), "Başlangıç ağ önemi",
        "PageRank value in the tissue-specific graph before target attenuation.",
        precision=6, scope="baseline_topology",
        supports=("relative baseline network importance",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "pagerank_perturbed", ("PageRank_Perturbed",), "Pertürbasyon sonrası ağ önemi",
        "PageRank value calculated after target attenuation.",
        precision=6, scope="perturbation_response",
        supports=("relative network importance after attenuation",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "hinterland", ("Hinterland_Skoru",), "Ağdaki Önemi",
        "Baseline rank-based composite of PageRank and degree in the constructed graph.",
        unit="score_0_100", precision=2, scope="baseline_topology",
        supports=("baseline structural prominence in this graph",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "betweenness_centrality", ("BC_Skoru",), "Geçiş Merkeziliği",
        "Weighted betweenness centrality in the baseline graph, with the reported calculation mode.",
        precision=3, scope="baseline_topology",
        supports=("baseline shortest-path mediation in this graph",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "gateway", ("Gümrük_Kapisi",), "Compartment Bottleneck",
        "Existing hybrid gateway marker based on cross-compartment neighborhood ratio and BC threshold.",
        unit="boolean", precision=None, scope="baseline_topology", kind=MetricKind.MARKER,
        supports=("the existing Compartment Bottleneck rule matched",), does_not_support=(
            "therapeutic target status", "biological necessity", "causality", "experimental validation",
        ),
    ),
    _contract(
        "degree", ("k_i",), "Derece",
        "Number of incident high-confidence STRING edges represented in the constructed graph summary.",
        unit="edge_count", precision=0, scope="baseline_topology",
        supports=("observed graph connectivity",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "community", ("Topluluk_ID",), "Topluluk",
        "Community identifier assigned by the existing graph community detection step.",
        unit="identifier", precision=None, scope="baseline_topology", kind=MetricKind.TEXT,
        supports=("shared algorithmic graph partition",), does_not_support=(
            "shared protein complex", "shared pathway", "physical colocalization", "causality",
        ),
    ),
    _contract(
        "drug_score", ("Drug_Score",), "İlaçlanabilirlik skoru",
        "Existing Sophiark score derived from available local structural and ligand fields.",
        unit="score_0_100", precision=2, scope="local_structural_context",
        supports=("the existing Drug Score value",), does_not_support=(
            "clinical efficacy", "clinical safety", "approved drug status", "experimental validation",
        ),
    ),
    _contract(
        "essentiality", ("Essentiality",), "Yaşamsallık (Essentiality)",
        "Existing external essentiality annotation displayed alongside the simulation result.",
        origin=MetricOrigin.LOCAL_ANNOTATION, precision=None, scope="annotation_context",
        kind=MetricKind.TEXT, supports=("the stored essentiality annotation",),
        does_not_support=("essentiality in the selected tissue", "causal simulation evidence", "therapeutic suitability"),
    ),
    _contract(
        "empirical_p", ("empirical_p",), "Ampirik p-değeri",
        "Existing empirical p-value from the degree-matched perturbation null model.",
        unit="probability", precision=6, scope="statistical_support",
        supports=("extremeness under the reported null model",), does_not_support=(
            "effect size", "biological importance", "experimental validation", "causality",
        ),
    ),
    _contract(
        "q_value", ("q_value",), "q-değeri (FDR)",
        "Existing Benjamini-Hochberg adjusted empirical p-value.",
        unit="probability", precision=6, scope="statistical_support",
        supports=("multiple-testing-adjusted support under the reported null model",),
        does_not_support=("effect size", "biological importance", "experimental validation", "causality"),
    ),
    _contract(
        "significant_redistribution", ("significant_redistribution",), "İstatistiksel destek",
        "Existing boolean flag indicating that q_value met the configured FDR threshold.",
        unit="boolean", precision=None, scope="statistical_support", kind=MetricKind.MARKER,
        supports=("the configured FDR rule matched",), does_not_support=(
            "biological significance", "experimental validation", "causality",
        ),
    ),
    _contract(
        "global_efficiency_baseline", ("Global_Efficiency_Baseline",), "Başlangıç global ağ verimliliği",
        "Weighted global efficiency before attenuation.", precision=8, scope="structural_efficiency",
        supports=("baseline weighted graph efficiency",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "global_efficiency_perturbed", ("Global_Efficiency_Perturbed",), "Pertürbasyon sonrası global ağ verimliliği",
        "Weighted global efficiency after attenuation.", precision=8, scope="structural_efficiency",
        supports=("post-attenuation weighted graph efficiency",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "global_efficiency_change", ("Global_Efficiency_Change",), "Global ağ verimliliği değişimi",
        "Perturbed minus baseline weighted global efficiency.", precision=8, signed=True,
        scope="structural_efficiency", supports=("absolute graph-efficiency change",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "global_efficiency_change_pct", ("Global_Efficiency_Change_Pct",), "Global ağ verimliliği değişimi (%)",
        "Percent change in weighted global efficiency after attenuation.", unit="percent", precision=2,
        signed=True, suffix="%", scope="structural_efficiency",
        supports=("relative graph-efficiency change",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "mean_local_efficiency_baseline", ("Mean_Local_Efficiency_Baseline",), "Başlangıç ortalama yerel verimlilik",
        "Mean weighted local efficiency before attenuation, using the reported sampling mode.",
        precision=8, scope="structural_efficiency",
        supports=("baseline mean local graph efficiency",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "mean_local_efficiency_perturbed", ("Mean_Local_Efficiency_Perturbed",), "Pertürbasyon sonrası ortalama yerel verimlilik",
        "Mean weighted local efficiency after attenuation, using the reported sampling mode.",
        precision=8, scope="structural_efficiency",
        supports=("post-attenuation mean local graph efficiency",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "mean_local_efficiency_change_pct", ("Mean_Local_Efficiency_Change_Pct",), "Ortalama yerel verimlilik değişimi (%)",
        "Percent change in mean weighted local efficiency after attenuation.", unit="percent",
        precision=2, signed=True, suffix="%", scope="structural_efficiency",
        supports=("relative mean local graph-efficiency change",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "target_local_efficiency_baseline", ("Target_Local_Efficiency_Baseline",), "Başlangıç hedef-yerel verimlilik",
        "Weighted efficiency of the reported closed one-hop target neighborhood before attenuation.",
        precision=8, scope="structural_efficiency",
        supports=("baseline target-neighborhood graph efficiency",), does_not_support=_TOPOLOGY_ONLY,
    ),
    _contract(
        "target_local_efficiency_perturbed", ("Target_Local_Efficiency_Perturbed",), "Pertürbasyon sonrası hedef-yerel verimlilik",
        "Weighted efficiency of the reported closed one-hop target neighborhood after attenuation.",
        precision=8, scope="structural_efficiency",
        supports=("post-attenuation target-neighborhood graph efficiency",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "target_local_efficiency_change_pct", ("Target_Local_Efficiency_Change_Pct",), "Hedef-yerel verimlilik değişimi (%)",
        "Percent change in the reported target-neighborhood weighted efficiency.", unit="percent",
        precision=2, signed=True, suffix="%", scope="structural_efficiency",
        supports=("relative target-neighborhood graph-efficiency change",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "systemic_network_shift_pct", ("Systemic_Network_Shift_Pct",), "Sistemik ağ kayması (%)",
        "Existing hybrid of PageRank total-variation distance and affected-node fraction.",
        unit="percent", precision=2, signed=False, suffix="%", scope="perturbation_response",
        supports=("system-wide redistribution magnitude in this model",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "legacy_efficiency_kayip_pct", ("Efficiency_Kayip_Pct",), "Sistemik ağ kayması (%) · legacy",
        "Legacy display field that stores Systemic_Network_Shift_Pct; despite its historical name, it is not global efficiency.",
        unit="percent", precision=2, signed=True, suffix="%", scope="perturbation_response",
        supports=("the legacy systemic-network-shift value",), does_not_support=(
            "global graph efficiency", *_NO_CAUSALITY,
        ),
    ),
    _contract(
        "top_positive_pagerank_mean_pct", ("Top_Positive_PageRank_Mean_Pct",), "Pozitif ağ önemi yanıtı ortalaması (%)",
        "Mean positive Delta_PageRank_Pct among the selected top responses.",
        unit="percent", precision=2, signed=False, suffix="%", scope="perturbation_response",
        supports=("mean positive relative PageRank response",), does_not_support=_NO_CAUSALITY,
    ),
    _contract(
        "legacy_local_efficiency_kayip_pct", ("Local_Efficiency_Kayip_Pct",), "Pozitif ağ önemi yanıtı ortalaması (%) · legacy",
        "Legacy display field that stores Top_Positive_PageRank_Mean_Pct; despite its historical name, it is not local efficiency.",
        unit="percent", precision=2, signed=True, suffix="%", scope="perturbation_response",
        supports=("the legacy mean-positive-PageRank-response value",), does_not_support=(
            "local graph efficiency", *_NO_CAUSALITY,
        ),
    ),
    _contract(
        "go_enrichment", ("Term", "Adjusted P-value", "Genes", "Pathway_Family", "Kaynak"), "GO zenginleştirmesi",
        "Existing GO over-representation result for the returned candidate set and active-network background.",
        origin=MetricOrigin.ONTOLOGY, precision=None, scope="returned_set_enrichment", kind=MetricKind.TABLE,
        supports=("GO terms over-represented in the tested returned set",), does_not_support=(
            "pathway activation", "pathway inhibition", "causality", "experimental validation",
        ),
    ),
    _contract(
        "kegg_enrichment", ("Term", "Adjusted P-value", "Genes", "Pathway_Family", "Kaynak"), "KEGG zenginleştirmesi",
        "Existing KEGG over-representation result for the returned candidate set and active-network background.",
        origin=MetricOrigin.PATHWAY, precision=None, scope="returned_set_enrichment", kind=MetricKind.TABLE,
        supports=("KEGG terms over-represented in the tested returned set",), does_not_support=(
            "pathway activation", "pathway inhibition", "causality", "experimental validation",
        ),
    ),
    _contract(
        "localization", ("Lokalizasyon", "GO_CC_Terimleri", "GO Hücresel Bileşen"), "Hücresel konum",
        "Stored cellular-component/localization annotation shown as biological context.",
        origin=MetricOrigin.LOCAL_ANNOTATION, precision=None, scope="annotation_context", kind=MetricKind.TEXT,
        supports=("the recorded localization annotation",), does_not_support=(
            "where the perturbation physically occurred", "protein abundance", "causality",
        ),
    ),
    _contract(
        "tf_annotations", ("Düzenleyici_TFler", "Hedef Gen Etkisi"), "Düzenleyici TF kanıtı",
        "Stored directed regulatory annotation, kept separate from the undirected STRING PPI simulation.",
        origin=MetricOrigin.CURATED_DATABASE, precision=None, scope="directed_evidence_overlay", kind=MetricKind.TEXT,
        supports=("a stored directed-regulation record exists",), does_not_support=(
            "causal effect in the simulated tissue", "absence of regulation when no record exists", "experimental validation",
        ),
    ),
)


CONTRACTS = MappingProxyType({contract.metric_id: contract for contract in _CONTRACT_LIST})
_SOURCE_INDEX: dict[str, tuple[str, ...]] = {}
for _contract_item in _CONTRACT_LIST:
    for _source_field in _contract_item.source_fields:
        _SOURCE_INDEX[_source_field] = (*_SOURCE_INDEX.get(_source_field, ()), _contract_item.metric_id)
SOURCE_FIELD_INDEX = MappingProxyType(_SOURCE_INDEX)


def contracts_for_source_field(field: str) -> tuple[InterpretationContract, ...]:
    return tuple(CONTRACTS[metric_id] for metric_id in SOURCE_FIELD_INDEX.get(field, ()))


def get_contract(metric_or_field: str) -> InterpretationContract | None:
    """Resolve a metric ID or an unambiguous source field."""

    direct = CONTRACTS.get(metric_or_field)
    if direct is not None:
        return direct
    matches = contracts_for_source_field(metric_or_field)
    return matches[0] if len(matches) == 1 else None


def get_metric_label(metric_or_field: str) -> str:
    contract = get_contract(metric_or_field)
    return contract.display_name if contract else metric_or_field


def get_metric_definition(metric_or_field: str) -> str | None:
    contract = get_contract(metric_or_field)
    return contract.definition if contract else None


def get_metric_limitations(metric_or_field: str) -> tuple[str, ...]:
    contract = get_contract(metric_or_field)
    return contract.does_not_support if contract else ()


def _is_missing(value: object) -> bool:
    if value is None:
        return True
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(missing) if isinstance(missing, (bool, np.bool_)) else False


def _truthy(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().casefold() in {"true", "1", "yes", "evet", "✓"}


def format_metric(metric_or_field: str, value: object) -> str:
    """Format an existing value without calculating or transforming it."""

    if _is_missing(value):
        return "—"
    contract = get_contract(metric_or_field)
    if contract is None or contract.kind in {MetricKind.TEXT, MetricKind.TABLE}:
        return str(value)
    if contract.kind is MetricKind.MARKER:
        return "✓" if _truthy(value) else "—"
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(numeric):
        return "—"
    precision = contract.precision if contract.precision is not None else 6
    sign = "+" if contract.signed and numeric >= 0 else ""
    return f"{sign}{numeric:.{precision}f}{contract.suffix}"


# Familiar aliases for callers migrating from the existing UI presentation API.
metric_label = get_metric_label
metric_definition = get_metric_definition
metric_limitations = get_metric_limitations
