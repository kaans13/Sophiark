"""Phase 8 acceptance tests for Europe PMC metadata/abstract retrieval."""

from __future__ import annotations

from typing import Any

from src.research.config import ProviderTimeout
from src.research.http_transport import JsonHttpResponse
from src.research.provider_runtime import ProviderRequest, ProviderStatus
from src.research.providers.literature import EUROPEPMC_SCHEMA_VERSION, EuropePMCProvider


TIMEOUT = ProviderTimeout(0.1, 0.2, 0.3)


class FakeTransport:
    def __init__(self, data: Any) -> None:
        self.data = data
        self.calls = []

    def get_json(self, url, *, params, timeout, headers=None):
        self.calls.append((url, params, timeout))
        return JsonHttpResponse(200, self.data, {}, "https://www.ebi.ac.uk/europepmc/webservices/rest/search")


def _request(**params) -> ProviderRequest:
    return ProviderRequest(
        provider="europepmc",
        operation="literature_search",
        taxon_id=9606,
        canonical_query='"GLP1R" AND "RAMP1" AND "Pancreas"',
        provider_schema_version=EUROPEPMC_SCHEMA_VERSION,
        query_version="literature-query-v1",
        params={
            "max_publications": 2,
            "matched_entities": ["9606:protein:ENSP1", "9606:protein:ENSP2"],
            "reason_shown": "Shown because of a current-result pattern.",
            **params,
        },
    )


def _payload() -> dict[str, Any]:
    return {
        "version": "6.9",
        "hitCount": 12,
        "nextCursorMark": "not-followed",
        "resultList": {
            "result": [
                {
                    "id": "12345678",
                    "source": "MED",
                    "pmid": "12345678",
                    "pmcid": "PMC123",
                    "doi": "10.1000/fixture",
                    "title": "GLP1R and RAMP1 in pancreatic cells",
                    "pubYear": "2024",
                    "journalTitle": "Fixture Journal",
                    "authorString": "A. Author et al.",
                    "abstractText": "<p>GLP1R and RAMP1 were measured in pancreatic cells.</p>",
                    "isOpenAccess": "Y",
                    "fullTextUrlList": {"fullTextUrl": [{"url": "https://pdf.example"}]},
                },
                {
                    "id": "AGR-2",
                    "source": "AGR",
                    "title": "A second metadata record",
                    "firstPublicationDate": "2023-01-02",
                },
                {"id": "ignored-without-title", "source": "MED"},
            ]
        },
    }


def test_europepmc_returns_bounded_normalized_publications_and_query_record() -> None:
    transport = FakeTransport(_payload())
    provider = EuropePMCProvider(transport=transport, enabled=True)

    result = provider.execute(_request(), TIMEOUT)

    assert result.status is ProviderStatus.AVAILABLE
    assert len(result.data["records"]) == 2
    assert result.data["matched_publication_count"] == 2
    assert result.data["provider_hit_count"] == 12
    assert result.data["query"] == {
        "query_text": '"GLP1R" AND "RAMP1" AND "Pancreas"',
        "query_version": "literature-query-v1",
        "matched_entities": ("9606:protein:ENSP1", "9606:protein:ENSP2"),
        "reason_shown": "Shown because of a current-result pattern.",
    }
    first = result.data["records"][0]
    assert first["publication_id"] == "MED:12345678"
    assert first["abstract"] == "GLP1R and RAMP1 were measured in pancreatic cells."
    assert first["year"] == 2024
    assert "full_text" not in first
    assert "pdf" not in first
    assert "score" not in first
    assert first["provenance"]["source_identifier"] == "MED:12345678"
    assert first["provenance"]["source_version"] == "6.9"
    assert result.metadata["full_text_downloaded"] is False
    assert result.metadata["page_followed"] is False
    assert result.metadata["legal"]["attribution"] == "Europe PMC"
    assert transport.calls[0][1] == {
        "query": '"GLP1R" AND "RAMP1" AND "Pancreas"',
        "resultType": "core",
        "format": "json",
        "pageSize": 2,
        "cursorMark": "*",
    }


def test_provider_no_results_uses_no_match_without_false_novelty() -> None:
    provider = EuropePMCProvider(
        transport=FakeTransport({"version": "6.9", "hitCount": 0, "resultList": {"result": []}}),
        enabled=True,
    )
    result = provider.execute(_request(), TIMEOUT)
    assert result.status is ProviderStatus.NO_MATCH
    assert result.data is None
    assert "novel" not in (result.message or "").casefold()


def test_provider_rejects_invalid_publication_limit_without_transport() -> None:
    transport = FakeTransport(_payload())
    provider = EuropePMCProvider(transport=transport, enabled=True)
    result = provider.execute(_request(max_publications="many"), TIMEOUT)
    assert result.status is ProviderStatus.PROVIDER_ERROR
    assert transport.calls == []

