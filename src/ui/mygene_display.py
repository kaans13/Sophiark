"""MyGene kayıtlarını yorum katmadan kullanıcıya gösterme yardımcıları."""

from __future__ import annotations

from collections.abc import Iterable, Mapping


def terms_from_record(record: Mapping[str, object] | None, field: str, limit: int = 4) -> str:
    """MyGene alanındaki gerçek terimleri, kaynak sırasını koruyarak döndürür."""
    if not record:
        return "—"
    raw = record.get(field)
    if not raw:
        return "—"
    values: Iterable[object] = raw if isinstance(raw, (list, tuple)) else (raw,)
    terms: list[str] = []
    seen: set[str] = set()
    for item in values:
        term = item.get("term") if isinstance(item, Mapping) else str(item)
        text = str(term or "").strip()
        if text and text.lower() not in {"nan", "none", "—"} and text not in seen:
            seen.add(text)
            terms.append(text)
        if len(terms) >= limit:
            break
    return "; ".join(terms) if terms else "—"


def mygene_display_fields(record: Mapping[str, object] | None) -> dict[str, str]:
    """Sentetik biyolojik yorum üretmeden MyGene annotation alanlarını hazırlar."""
    record = record or {}
    name = str(record.get("name") or record.get("fullname") or "—").strip() or "—"
    return {
        "MyGene Adı": name,
        "GO Biyolojik Süreç": terms_from_record(record, "go_bp"),
        "GO Moleküler İşlev": terms_from_record(record, "go_mf"),
        "GO Hücresel Bileşen": terms_from_record(record, "go_cc"),
    }
