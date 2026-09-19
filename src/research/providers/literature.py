"""On-demand Europe PMC adapter for PubMed-indexed literature context.

Only publication metadata and abstracts returned by the official search API
are used.  The adapter never downloads PDFs or full text and never turns a
query match into a validation, support, contradiction, or causal claim.
"""

from __future__ import annotations

from datetime import datetime, timezone
import html
import re
from types import MappingProxyType
from typing import Any, Mapping

from src.research.config import ProviderTimeout
from src.research.http_transport import JsonTransport
from src.research.provider_runtime import ProviderCapability, ProviderRequest, ProviderResult, ProviderStatus
from src.research.providers.external_biology import JsonEvidenceProvider


EUROPEPMC_SCHEMA_VERSION = "europepmc-core-json-v1"
_TAG = re.compile(r"<[^>]+>")


def _now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str | None:
    if value is None:
        return None
    value = " ".join(html.unescape(_TAG.sub(" ", str(value))).split())
    return value or None


def _year(value: Any) -> int | None:
    text = _text(value)
    if text is None:
        return None
    match = re.search(r"(?:19|20)[0-9]{2}", text)
    return int(match.group(0)) if match else None


class EuropePMCProvider(JsonEvidenceProvider):
    """Bounded `resultType=core` search over Europe PMC/PubMed metadata."""

    provider_name = "europepmc"
    operation = "literature_search"
    legal_metadata = MappingProxyType({
        "provider_id": "europepmc",
        "attribution": "Europe PMC",
        "license_note": "Publication content remains subject to its source copyright and license terms.",
        "cache_policy": "Metadata and abstracts only; seven-day default research cache.",
        "commercial_use_note": "Review Europe PMC and underlying publication terms before commercial reuse.",
    })

    def __init__(self, *, transport: JsonTransport | None = None, enabled: bool = False) -> None:
        super().__init__(
            self.provider_name,
            (ProviderCapability(self.operation, (9606, 10090), "Europe PMC publication metadata and abstracts"),),
            transport=transport,
            enabled=enabled,
        )

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        query = " ".join(request.canonical_query.split())
        if len(query) < 3:
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Literature query is empty or too broad")
        maximum = request.params.get("max_publications", 10)
        if isinstance(maximum, bool) or not isinstance(maximum, int):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="max_publications must be an integer")
        maximum = min(max(maximum, 1), 100)
        outcome = self._get(
            request,
            timeout,
            url="https://www.ebi.ac.uk/europepmc/webservices/rest/search",
            params={
                "query": query,
                "resultType": "core",
                "format": "json",
                "pageSize": maximum,
                "cursorMark": "*",
            },
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        payload = outcome.data
        result_list = payload.get("resultList") if isinstance(payload, Mapping) else None
        raw_records = result_list.get("result") if isinstance(result_list, Mapping) else None
        if not isinstance(raw_records, list):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="Europe PMC response has an invalid schema")
        retrieved_at = _now()
        source_version = _text(payload.get("version")) if isinstance(payload, Mapping) else None
        matched_entities = request.params.get("matched_entities", ())
        if not isinstance(matched_entities, (tuple, list)):
            matched_entities = ()
        reason_shown = _text(request.params.get("reason_shown")) or (
            "Shown for the current explicit literature query."
        )
        query_record = {
            "query_text": query,
            "query_version": request.query_version,
            "matched_entities": list(matched_entities),
            "reason_shown": reason_shown,
        }
        records = []
        seen: set[str] = set()
        for raw in raw_records[:maximum]:
            if not isinstance(raw, Mapping):
                continue
            source = _text(raw.get("source")) or "EPMC"
            external_id = _text(raw.get("id"))
            title = _text(raw.get("title"))
            if external_id is None or title is None:
                continue
            publication_id = f"{source}:{external_id}"
            if publication_id in seen:
                continue
            seen.add(publication_id)
            pmid = _text(raw.get("pmid"))
            if pmid is None and source == "MED" and external_id.isdigit():
                pmid = external_id
            records.append({
                "publication_id": publication_id,
                "pmid": pmid,
                "pmcid": _text(raw.get("pmcid")),
                "doi": _text(raw.get("doi")),
                "title": title,
                "year": _year(raw.get("pubYear") or raw.get("firstPublicationDate")),
                "journal": _text(raw.get("journalTitle")),
                "authors_summary": _text(raw.get("authorString")),
                "abstract": _text(raw.get("abstractText")),
                "source": source,
                "retrieved_at": retrieved_at,
                "is_open_access": str(raw.get("isOpenAccess", "")).upper() == "Y",
                "provenance": {
                    "source": "Europe PMC",
                    "source_identifier": publication_id,
                    "retrieved_at": retrieved_at,
                    "provider": self.provider_name,
                    "source_version": source_version,
                    "locator": f"https://europepmc.org/article/{source}/{external_id}",
                },
            })
        status = ProviderStatus.AVAILABLE if records else ProviderStatus.NO_MATCH
        return ProviderResult(
            request,
            status,
            data={
                "query": query_record,
                "records": records,
                "matched_publication_count": len(records),
                "provider_hit_count": payload.get("hitCount"),
            } if records else None,
            provider_version=source_version,
            metadata={
                "retrieved_at": retrieved_at,
                "page_followed": False,
                "full_text_downloaded": False,
                "legal": dict(self.legal_metadata),
            },
        )


__all__ = ["EUROPEPMC_SCHEMA_VERSION", "EuropePMCProvider"]
