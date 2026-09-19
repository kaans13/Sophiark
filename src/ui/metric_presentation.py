"""Shared display conventions for existing scientific metrics.

The helpers in this module only define labels, tooltips, missing-value text,
and number formatting. They never transform or write scientific source data.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Literal

import numpy as np
import pandas as pd
from src.ui.localization import localize_text
@dataclass(frozen=True)
class MetricPresentation:
    label: str
    help: str
    digits: int | None = None
    signed: bool = False
    suffix: str = ""
    kind: Literal["number", "text", "boolean", "status"] = "number"
    short_label: str | None = None
    unit: str | None = None
    group: Literal["identity", "response", "rank", "comparison", "evidence", "context", "engine"] = "engine"
    order: int = 100
    missing: str = "—"
    fraction_as_percent: bool = False
    true_text: str = "Evet"
    false_text: str = "Hayır"


METRICS: dict[str, MetricPresentation] = {
    "Delta_PageRank_Pct": MetricPresentation(
        "Ağ önemi değişimi (%)",
        "Hedef zayıflatması sonrası göreli yapısal ağ önemi değişimi; ekspresyon veya aktivite değişimi değildir.",
        digits=4,
        signed=True,
        suffix="%",
        unit="percent", group="response", order=10,
    ),
    "Hinterland_Skoru": MetricPresentation(
        "Ağdaki Önemi",
        "Pertürbasyon öncesi başlangıç ağ-topolojisi bağlamı.",
        digits=2,
        group="rank", order=30,
    ),
    "BC_Skoru": MetricPresentation(
        "Geçiş Merkeziliği",
        "Pertürbasyon öncesi betweenness centrality bağlamı.",
        digits=3,
        group="engine", order=20,
    ),
    "PageRank_Baseline": MetricPresentation(
        "Başlangıç ağ önemi",
        "Pertürbasyon öncesi PageRank değeri.",
        digits=6,
        group="response", order=1,
    ),
    "PageRank_Perturbed": MetricPresentation(
        "Pertürbasyon sonrası ağ önemi",
        "Hedef zayıflatması sonrasında hesaplanmış PageRank değeri.",
        digits=6,
        group="response", order=2,
    ),
    "Gümrük_Kapisi": MetricPresentation(
        "Compartment Bottleneck",
        "Ağ bölgeleri veya hücresel bağlamlar arasında köprü rolü gösterebilen düğüm işareti.",
        kind="boolean", group="context", order=10, missing="Kullanılamıyor", true_text="✓",
    ),
    "Topluluk_ID": MetricPresentation(
        "Topluluk kimliği",
        "Başlangıç ağındaki kaydın mevcut topluluk kimliği.",
        kind="text", group="context", order=11,
    ),
    "Efficiency_Kayip_Pct": MetricPresentation(
        "Sistemik ağ kayması (%) · legacy",
        "Systemic_Network_Shift_Pct için geriye dönük alan adı; global ağ verimliliği değildir.",
        digits=2,
        signed=True,
        suffix="%",
    ),
    "Local_Efficiency_Kayip_Pct": MetricPresentation(
        "Pozitif ağ önemi yanıtı ortalaması (%) · legacy",
        "Top_Positive_PageRank_Mean_Pct için geriye dönük alan adı; yerel ağ verimliliği değildir.",
        digits=2,
        signed=True,
        suffix="%",
    ),
    "empirical_p": MetricPresentation(
        "Ampirik p-değeri",
        "Mevcut null-model analizinin ürettiği ampirik p-değeri; burada yeniden hesaplanmaz.",
        digits=6,
    ),
    "q_value": MetricPresentation(
        "q-değeri (FDR)",
        "Mevcut çoklu test düzeltmesinin ürettiği q-değeri; burada yeniden hesaplanmaz.",
        digits=6,
    ),
    "significant_redistribution": MetricPresentation(
        "İstatistiksel destek",
        "Mevcut FDR sonucundaki significant_redistribution işareti.",
        kind="boolean", group="comparison", order=30, missing="Kullanılamıyor", true_text="✓",
    ),
    "Essentiality": MetricPresentation(
        "Yaşamsallık (Essentiality)",
        "Mevcut dış essentiality anotasyonu; simülasyon tarafından üretilmez.",
        kind="text",
    ),
    "Drug_Score": MetricPresentation(
        "İlaçlanabilirlik skoru",
        "Mevcut anotasyon kaynağındaki Drug Score değeri; yeni bir birleşik skor değildir.",
        digits=2,
    ),
    "entity_id": MetricPresentation("ENSP", "Canonical protein/entity kimliği.", kind="text", group="identity", order=2),
    "symbol": MetricPresentation("Gen", "Bulunduğunda canonical gen sembolü.", kind="text", group="identity", order=1),
    "Symbol": MetricPresentation("Gen", "Bulunduğunda canonical gen sembolü.", kind="text", group="identity", order=1),
    "gene": MetricPresentation("ENSP", "Canonical protein kimliği.", kind="text", group="identity", order=2),
    "target": MetricPresentation("Hedef", "Analiz hedefi.", kind="text", group="identity", order=3),
    "tissue": MetricPresentation("Doku", "Analizin doku bağlamı.", kind="text", group="identity", order=4),
    "finding_class": MetricPresentation("Bulgu sınıfı", "Unified sınıflandırma kodu.", kind="status", group="comparison", order=1),
    "classic_response": MetricPresentation("Classic yanıtı (%)", "Classic işaretli yeniden dağılım yanıtı.", digits=4, signed=True, suffix="%", unit="percent", group="response", order=10),
    "directed_response": MetricPresentation("Directed yanıtı (%)", "Directed işaretli yeniden dağılım yanıtı.", digits=4, signed=True, suffix="%", unit="percent", group="response", order=11),
    "classic_baseline": MetricPresentation("Classic başlangıç PageRank", "Classic başlangıç PageRank değeri.", digits=8, group="response", order=1),
    "classic_perturbed": MetricPresentation("Classic pertürbe PageRank", "Classic pertürbasyon sonrası PageRank değeri.", digits=8, group="response", order=2),
    "directed_baseline": MetricPresentation("Directed başlangıç PageRank", "Directed başlangıç PageRank değeri.", digits=8, group="response", order=3),
    "directed_perturbed": MetricPresentation("Directed pertürbe PageRank", "Directed pertürbasyon sonrası PageRank değeri.", digits=8, group="response", order=4),
    "Directed_Redistribution_Pct": MetricPresentation("Directed yanıtı (%)", "Directed işaretli yeniden dağılım yanıtı.", digits=4, signed=True, suffix="%", unit="percent", group="response", order=11),
    "Directed_PageRank_Baseline": MetricPresentation("Directed başlangıç PageRank", "Directed başlangıç PageRank değeri.", digits=8, group="response", order=3),
    "Directed_PageRank_Perturbed": MetricPresentation("Directed pertürbe PageRank", "Directed pertürbasyon sonrası PageRank değeri.", digits=8, group="response", order=4),
    "Directed_BC": MetricPresentation("Directed BC", "Directed betweenness centrality değeri.", digits=6, group="engine", order=20),
    "classic_absolute_response_rank": MetricPresentation("Classic mutlak yanıt sırası", "Classic mutlak yanıt sırası.", digits=0, group="rank", order=1),
    "directed_absolute_response_rank": MetricPresentation("Directed mutlak yanıt sırası", "Directed mutlak yanıt sırası.", digits=0, group="rank", order=2),
    "classic_absolute_rank_percentile": MetricPresentation("Classic sıra yüzdeliği", "0–1 aralığında sıra fraksiyonu; 1 en yüksektir.", digits=4, unit="fraction", group="rank", order=10),
    "directed_absolute_rank_percentile": MetricPresentation("Directed sıra yüzdeliği", "0–1 aralığında sıra fraksiyonu; 1 en yüksektir.", digits=4, unit="fraction", group="rank", order=11),
    "direction_sensitivity": MetricPresentation("Yön duyarlılığı", "0–1 ölçeğinde mutlak yüzdelik farkı.", digits=4, unit="fraction", group="comparison", order=10),
    "percentile_difference": MetricPresentation("Yüzdelik farkı", "0–1 ölçeğinde Directed eksi Classic sıra yüzdeliği.", digits=4, signed=True, unit="fraction", group="comparison", order=11),
    "classic_available": MetricPresentation("Classic kullanılabilir", "Classic sonucunun kullanılabilirliği.", kind="boolean", group="comparison", order=20, missing="Kullanılamıyor"),
    "directed_available": MetricPresentation("Directed kullanılabilir", "Directed sonucunun kullanılabilirliği.", kind="boolean", group="comparison", order=21, missing="Kullanılamıyor"),
    "same_sign": MetricPresentation("Aynı işaret", "Classic ve Directed yanıt işareti uyumu.", kind="boolean", group="comparison", order=22, missing="Kullanılamıyor"),
    "sign_flip": MetricPresentation("İşaret değişimi", "Nötr olmayan karşıt yanıt işaretleri.", kind="boolean", group="comparison", order=23, missing="Kullanılamıyor"),
    "sign_agreement_rate": MetricPresentation("İşaret uyumu", "Uyum fraksiyonu yüzde olarak gösterilir.", digits=1, suffix="%", unit="fraction", group="comparison", order=24, fraction_as_percent=True),
    "Evidence_Provenance_Available": MetricPresentation("Evidence provenance kullanılabilir", "STRING kenar provenance sorgusunun kullanılabilirliği.", kind="boolean", group="evidence", order=1, missing="Kullanılamıyor"),
    "Evidence_Provenance_Status": MetricPresentation("Evidence provenance durumu", "Makine tarafından okunabilir provenance durumu.", kind="status", group="evidence", order=2),
    "Evidence_combined_score": MetricPresentation("STRING combined score", "Saklanan STRING skor fraksiyonu.", digits=3, unit="fraction", group="evidence", order=10),
    "Evidence_s_nontext": MetricPresentation("STRING non-text score", "Text-mining kanalları hariç STRING skoru.", digits=3, unit="fraction", group="evidence", order=11),
    "Evidence_text_dependency": MetricPresentation("Metin bağımlılığı", "Combined score içinde metin kanallarına atfedilen pay.", digits=3, unit="fraction", group="evidence", order=12),
    "Evidence_Experimental_Support_Present": MetricPresentation("Deneysel destek", "Evidence kaynağı kontrol edildi: destek var veya yok.", kind="boolean", group="evidence", order=20, missing="Kullanılamıyor"),
    "Evidence_Database_Support_Present": MetricPresentation("Veritabanı desteği", "Evidence kaynağı kontrol edildi: destek var veya yok.", kind="boolean", group="evidence", order=21, missing="Kullanılamıyor"),
    "CORUM_Target_Complex_Context": MetricPresentation("CORUM kompleks bağlamı", "Hedef ilişkili CORUM kompleks adları.", kind="text", group="context", order=1),
    "CORUM_Target_Complex_Member": MetricPresentation("CORUM kompleks üyesi", "Adayın hedef kompleks bağlamında olup olmadığı.", kind="boolean", group="context", order=2, missing="Kullanılamıyor"),
}


STATUS_LABELS = {
    "READY": "Hazır", "COMPLETE": "Tamamlandı", "PARTIAL": "Kısmi",
    "ENGINE_PARTIAL": "Motor kısmi", "FAILED": "Başarısız",
    "DATA_UNAVAILABLE": "Veri kullanılamıyor", "UNAVAILABLE": "Kullanılamıyor",
    "AVAILABLE": "Kullanılabilir", "BUNDLE_CORRUPT": "Bundle bozuk",
    "UNSUPPORTED_BUNDLE_VERSION": "Desteklenmeyen bundle sürümü",
    "ID_CONFLICT": "Kimlik çakışması",
    "INDEXED_TARGET_EDGE": "Hedef kenarı indekslendi",
    "NO_DIRECT_INDEXED_TARGET_EDGE": "Doğrudan indeksli hedef kenarı yok",
    "INDEX_UNAVAILABLE": "İndeks kullanılamıyor",
    "ROBUST_CROSS_MODEL": "Modeller arası sağlam",
    "DIRECTION_SENSITIVE": "Yön duyarlı", "CLASSIC_DOMINANT": "Classic baskın",
    "DIRECTED_EMERGENT": "Directed belirgin", "LOW_AGREEMENT": "Düşük uyum",
}

_GROUP_ORDER = {name: index for index, name in enumerate(
    ("identity", "response", "rank", "comparison", "evidence", "context", "engine")
)}


def metric_spec(field: str) -> MetricPresentation | None:
    """Return the display specification for a source field, if registered."""
    return METRICS.get(field)


def metric_label(field: str) -> str:
    """Return one stable UI label without renaming the source column."""
    spec = metric_spec(field)
    return localize_text(spec.label) if spec else localize_text(field)


def metric_help(field: str) -> str | None:
    spec = metric_spec(field)
    return localize_text(spec.help) if spec else None


def metric_printf(field: str) -> str | None:
    """Return a Streamlit-compatible printf display format."""
    spec = metric_spec(field)
    if spec is None or spec.digits is None:
        return None
    sign = "+" if spec.signed else ""
    suffix = spec.suffix.replace("%", "%%")
    return f"%{sign}.{spec.digits}f{suffix}"


def display_status(value: object) -> str:
    return localize_text(STATUS_LABELS.get(str(value), str(value)))


def ordered_columns(columns: Iterable[str]) -> list[str]:
    """Deterministic identity→response→rank→comparison→evidence ordering."""
    indexed = list(enumerate(columns))
    return [column for _, column in sorted(indexed, key=lambda item: (
        _GROUP_ORDER.get(METRICS[item[1]].group, 99) if item[1] in METRICS else 99,
        METRICS[item[1]].order if item[1] in METRICS else item[0],
        item[0],
    ))]


def prepare_display_frame(frame: pd.DataFrame) -> pd.DataFrame:
    """Build a presentation-only view; source/export values stay untouched."""
    display = frame.loc[:, ordered_columns(frame.columns)].copy()
    for field in display.columns:
        spec = metric_spec(field)
        if spec and (spec.kind in {"boolean", "status"} or spec.fraction_as_percent):
            display[field] = display[field].map(lambda value, name=field: format_metric(name, value))
    return display


def metric_column_config(field: str) -> object | None:
    spec = metric_spec(field)
    if spec is None:
        return None
    if spec.kind in {"text", "boolean", "status"} or spec.fraction_as_percent:
        return metric_text_column(field)
    return metric_number_column(field)


def is_missing(value: object) -> bool:
    if value is None:
        return True
    if isinstance(value, (list, tuple, set, dict)):
        return not bool(value)
    try:
        missing = pd.isna(value)
    except (TypeError, ValueError):
        return False
    return bool(missing) if isinstance(missing, (bool, np.bool_)) else False


def truthy(value: object) -> bool:
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    return str(value).strip().casefold() in {"true", "1", "yes", "evet", "✓"}


def format_metric(field: str, value: object) -> str:
    """Format a value for detail views without changing the value itself."""
    if is_missing(value):
        spec = metric_spec(field)
        return localize_text(spec.missing) if spec else "—"
    spec = metric_spec(field)
    if spec is None or spec.kind == "text":
        return localize_text(str(value))
    if spec.kind == "status":
        return display_status(value)
    if spec.kind == "boolean":
        return localize_text(spec.true_text if truthy(value) else spec.false_text)
    try:
        numeric = float(value)
    except (TypeError, ValueError):
        return localize_text(str(value))
    if spec.fraction_as_percent:
        numeric *= 100.0
    sign = "+" if spec.signed and numeric >= 0 else ""
    suffix = spec.suffix
    return f"{sign}{numeric:.{spec.digits}f}{suffix}"


def metric_number_column(field: str, *, width: str | int | None = None) -> object:
    """Build a numeric column config from the shared presentation convention."""
    import streamlit as st
    return st.column_config.NumberColumn(
        metric_label(field),
        help=metric_help(field),
        format=metric_printf(field),
        width=width,
    )


def metric_text_column(field: str, *, width: str | int | None = None) -> object:
    """Build a text/marker column config from the shared convention."""
    import streamlit as st
    return st.column_config.TextColumn(
        metric_label(field),
        help=metric_help(field),
        width=width,
    )
