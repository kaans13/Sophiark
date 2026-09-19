"""Phase 7 contract tests for bounded external biology adapters."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from typing import Any, Mapping

import pytest

from src.research.config import ProviderTimeout
from src.research.http_transport import JsonHttpResponse
from src.research.provider_runtime import ProviderRequest, ProviderStatus
from src.research.providers.external_biology import (
    INTERPRO_SCHEMA_VERSION,
    QUICKGO_SCHEMA_VERSION,
    REACTOME_SCHEMA_VERSION,
    UNIPROT_SCHEMA_VERSION,
    InterProProvider,
    QuickGOProvider,
    ReactomeProvider,
    UniProtProvider,
)


TIMEOUT = ProviderTimeout(0.1, 0.2, 0.3)


class FakeTransport:
    def __init__(self, response: JsonHttpResponse | BaseException) -> None:
        self.response = response
        self.calls: list[dict[str, Any]] = []

    def get_json(self, url, *, params, timeout, headers=None):
        self.calls.append({"url": url, "params": params, "timeout": timeout, "headers": headers})
        if isinstance(self.response, BaseException):
            raise self.response
        return self.response


def _response(data: Any, status: int = 200, headers: Mapping[str, str] | None = None) -> JsonHttpResponse:
    return JsonHttpResponse(status, data, headers or {}, "https://official.example/record")


def _request(provider: str, operation: str, schema: str, query: str = "P43220", **changes) -> ProviderRequest:
    values = dict(
        provider=provider,
        operation=operation,
        taxon_id=9606,
        canonical_query=query,
        provider_schema_version=schema,
        query_version="research-context-query-v1",
    )
    values.update(changes)
    return ProviderRequest(**values)


def _assert_provenance(record: Mapping[str, Any], source_identifier: str) -> None:
    provenance = record["provenance"]
    assert provenance["source"]
    assert provenance["source_identifier"] == source_identifier
    assert provenance["retrieved_at"].endswith("+00:00")
    assert provenance["provider"]
    assert "source_version" in provenance
    assert provenance["locator"] == "https://official.example/record"


def test_uniprot_parses_only_source_fields_and_rejects_cross_taxon_record() -> None:
    payload = {
        "primaryAccession": "P43220",
        "entryType": "UniProtKB reviewed (Swiss-Prot)",
        "organism": {"taxonId": 9606, "scientificName": "Homo sapiens"},
        "proteinDescription": {"recommendedName": {"fullName": {"value": "Glucagon-like peptide 1 receptor"}}},
        "genes": [{"geneName": {"value": "GLP1R"}}],
        "comments": [
            {"commentType": "FUNCTION", "texts": [{"value": "Source-provided function statement."}]},
            {"commentType": "SUBCELLULAR LOCATION", "subcellularLocations": [{"location": {"value": "Cell membrane"}}]},
        ],
        "uniProtKBCrossReferences": [{"database": "Ensembl", "id": "ENSP00000362353"}],
    }
    transport = FakeTransport(_response(payload, headers={"X-UniProt-Release": "2026_03"}))
    provider = UniProtProvider(transport=transport, enabled=True)
    request = _request("uniprot", provider.operation, UNIPROT_SCHEMA_VERSION)

    result = provider.execute(request, TIMEOUT)

    assert result.status is ProviderStatus.AVAILABLE
    record = result.data["records"][0]
    assert record["protein_name"] == "Glucagon-like peptide 1 receptor"
    assert record["function_statements"] == ("Source-provided function statement.",)
    assert record["subcellular_locations"] == ("Cell membrane",)
    assert "family" not in record
    assert "score" not in record
    _assert_provenance(record, "P43220")
    assert result.provider_version == "2026_03"
    with pytest.raises(TypeError):
        record["protein_name"] = "changed"

    cross_taxon = UniProtProvider(
        transport=FakeTransport(_response({**payload, "organism": {"taxonId": 10090}})),
        enabled=True,
    ).execute(request, TIMEOUT)
    assert cross_taxon.status is ProviderStatus.NO_MATCH
    assert cross_taxon.data is None


def test_interpro_is_bounded_and_describes_entries_without_family_inference() -> None:
    payload = {
        "count": 200,
        "next": "https://next-page.example",
        "results": [
            {"metadata": {"accession": f"IPR{i:06d}", "name": f"Entry {i}", "type": "domain", "source_database": "interpro"}}
            for i in range(25)
        ],
    }
    transport = FakeTransport(_response(payload))
    provider = InterProProvider(transport=transport, enabled=True)
    request = _request(
        "interpro",
        provider.operation,
        INTERPRO_SCHEMA_VERSION,
        params={"page_size": 3},
    )

    result = provider.execute(request, TIMEOUT)

    assert result.status is ProviderStatus.AVAILABLE
    assert len(result.data["records"]) == 3
    assert result.data["total_available"] == 200
    assert result.metadata["page_followed"] is False
    assert transport.calls[0]["params"] == {"page_size": 3}
    _assert_provenance(result.data["records"][0], "IPR000000")
    assert "family_membership" not in result.data["records"][0]


def test_quickgo_annotation_filters_cross_species_and_preserves_evidence() -> None:
    payload = {
        "numberOfHits": 2,
        "results": [
            {"geneProductId": "UniProtKB:P43220", "goId": "GO:0005886", "goName": "plasma membrane", "goAspect": "cellular_component", "evidenceCode": "ECO:0000269", "reference": "PMID:1", "taxonId": 9606},
            {"geneProductId": "UniProtKB:P43220", "goId": "GO:1234567", "taxonId": 10090},
        ],
    }
    provider = QuickGOProvider(transport=FakeTransport(_response(payload)), enabled=True)
    request = _request(
        "quickgo",
        provider.annotation_operation,
        QUICKGO_SCHEMA_VERSION,
        params={"limit": 25},
    )

    result = provider.execute(request, TIMEOUT)

    assert result.status is ProviderStatus.AVAILABLE
    assert len(result.data["records"]) == 1
    record = result.data["records"][0]
    assert record["taxon_id"] == 9606
    assert record["evidence_code"] == "ECO:0000269"
    assert result.metadata["cross_taxon_records_rejected"] == 1
    _assert_provenance(record, "P43220:GO:0005886")


def test_quickgo_exact_term_and_reactome_mapping_have_provenance() -> None:
    term_provider = QuickGOProvider(
        transport=FakeTransport(_response({"results": [{"id": "GO:0005886", "name": "plasma membrane", "aspect": "cellular_component", "definition": {"text": "A membrane."}}]})),
        enabled=True,
    )
    term_request = _request("quickgo", term_provider.term_operation, QUICKGO_SCHEMA_VERSION, query="GO:0005886")
    term = term_provider.execute(term_request, TIMEOUT)
    assert term.status is ProviderStatus.AVAILABLE
    _assert_provenance(term.data["records"][0], "GO:0005886")

    reactome_provider = ReactomeProvider(
        transport=FakeTransport(_response([{"stId": "R-HSA-123", "dbId": 123, "displayName": "Fixture pathway", "speciesName": "Homo sapiens", "isInferred": False, "isInDisease": False}])) ,
        enabled=True,
    )
    reactome_request = _request("reactome", reactome_provider.operation, REACTOME_SCHEMA_VERSION, params={"limit": 5})
    pathway = reactome_provider.execute(reactome_request, TIMEOUT)
    assert pathway.status is ProviderStatus.AVAILABLE
    assert pathway.data["records"][0]["requested_taxon_id"] == 9606
    assert "causal" not in pathway.data["records"][0]
    _assert_provenance(pathway.data["records"][0], "R-HSA-123")


def test_reactome_rejects_explicit_cross_species_records() -> None:
    provider = ReactomeProvider(
        transport=FakeTransport(_response([
            {"stId": "R-MMU-1", "displayName": "Mouse pathway", "speciesName": "Mus musculus"},
            {"stId": "R-HSA-1", "displayName": "Human pathway", "speciesName": "Homo sapiens"},
        ])),
        enabled=True,
    )
    request = _request("reactome", provider.operation, REACTOME_SCHEMA_VERSION)
    result = provider.execute(request, TIMEOUT)
    assert result.status is ProviderStatus.AVAILABLE
    assert [item["pathway_id"] for item in result.data["records"]] == ["R-HSA-1"]
    assert result.metadata["cross_taxon_records_rejected"] == 1


@pytest.mark.parametrize(
    ("status_code", "expected"),
    ((404, ProviderStatus.NO_MATCH), (408, ProviderStatus.TIMEOUT), (429, ProviderStatus.RATE_LIMIT), (500, ProviderStatus.PROVIDER_ERROR)),
)
def test_http_statuses_remain_distinct(status_code: int, expected: ProviderStatus) -> None:
    provider = UniProtProvider(transport=FakeTransport(_response({}, status=status_code)), enabled=True)
    result = provider.execute(_request("uniprot", provider.operation, UNIPROT_SCHEMA_VERSION), TIMEOUT)
    assert result.status is expected


@pytest.mark.parametrize(
    ("error", "expected"),
    ((TimeoutError("slow"), ProviderStatus.TIMEOUT), (ConnectionError("offline"), ProviderStatus.OFFLINE), (ValueError("bad json"), ProviderStatus.PROVIDER_ERROR)),
)
def test_transport_failures_are_fail_soft(error: BaseException, expected: ProviderStatus) -> None:
    provider = UniProtProvider(transport=FakeTransport(error), enabled=True)
    result = provider.execute(_request("uniprot", provider.operation, UNIPROT_SCHEMA_VERSION), TIMEOUT)
    assert result.status is expected


def test_non_exact_identifiers_never_reach_transport() -> None:
    transport = FakeTransport(_response({}))
    provider = UniProtProvider(transport=transport, enabled=True)
    request = _request("uniprot", provider.operation, UNIPROT_SCHEMA_VERSION, query="GLP1R")
    assert provider.execute(request, TIMEOUT).status is ProviderStatus.NO_MATCH
    assert transport.calls == []
