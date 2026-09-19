"""Verified, on-demand adapters for external biological context providers.

These adapters are read-only, disabled unless explicitly configured, bounded
to small entity-level responses, and never participate in the scientific
simulation path.  They preserve source statements as context; they do not
compute scores, significance, causality, or family membership heuristics.
"""

from __future__ import annotations

from datetime import datetime, timezone
import re
from typing import Any, Iterable, Mapping
from urllib.parse import quote

from src.research.config import ProviderTimeout
from src.research.http_transport import JsonHttpResponse, JsonTransport, RequestsJsonTransport
from src.research.provider_runtime import (
    EvidenceProvider,
    ProviderCapability,
    ProviderRequest,
    ProviderResult,
    ProviderStatus,
)


SUPPORTED_TAXA = (9606, 10090)
UNIPROT_SCHEMA_VERSION = "uniprot-json-v1"
INTERPRO_SCHEMA_VERSION = "interpro-json-v1"
QUICKGO_SCHEMA_VERSION = "quickgo-json-v1"
REACTOME_SCHEMA_VERSION = "reactome-json-v1"

_UNIPROT_ACCESSION = re.compile(r"^[A-Z0-9][A-Z0-9]{5,9}(?:-[1-9][0-9]*)?$")
_GO_ID = re.compile(r"^GO:[0-9]{7}$")


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _text(value: Any) -> str | None:
    if value is None:
        return None
    normalized = str(value).strip()
    return normalized or None


def _nested_text(value: Any, *path: str) -> str | None:
    current = value
    for key in path:
        if not isinstance(current, Mapping):
            return None
        current = current.get(key)
    return _text(current)


def _items(value: Any) -> tuple[Mapping[str, Any], ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(item for item in value if isinstance(item, Mapping))


def _bounded_unique(values: Iterable[str | None], limit: int) -> list[str]:
    result: list[str] = []
    seen: set[str] = set()
    for value in values:
        text = _text(value)
        if text is None or text in seen:
            continue
        seen.add(text)
        result.append(text)
        if len(result) >= limit:
            break
    return result


def _source_version(response: JsonHttpResponse) -> str | None:
    return (
        response.header("x-uniprot-release")
        or response.header("x-api-release")
        or response.header("x-reactome-release")
    )


def _provenance(
    *,
    provider: str,
    source: str,
    source_identifier: str,
    response: JsonHttpResponse,
    retrieved_at: str,
) -> dict[str, Any]:
    return {
        "source": source,
        "source_identifier": source_identifier,
        "retrieved_at": retrieved_at,
        "provider": provider,
        "source_version": _source_version(response),
        "locator": response.url,
    }


class JsonEvidenceProvider(EvidenceProvider):
    def __init__(
        self,
        provider_name: str,
        capabilities: Iterable[ProviderCapability],
        *,
        transport: JsonTransport | None,
        enabled: bool,
    ) -> None:
        super().__init__(provider_name, capabilities, enabled=enabled)
        self._transport = transport or RequestsJsonTransport()

    def _get(
        self,
        request: ProviderRequest,
        timeout: ProviderTimeout,
        *,
        url: str,
        params: Mapping[str, Any] | None = None,
    ) -> JsonHttpResponse | ProviderResult:
        try:
            response = self._transport.get_json(url, params=params, timeout=timeout)
        except TimeoutError as exc:
            return ProviderResult(request, ProviderStatus.TIMEOUT, message=str(exc) or "Provider timed out")
        except ConnectionError as exc:
            return ProviderResult(request, ProviderStatus.OFFLINE, message=str(exc) or "Provider is offline")
        except Exception as exc:
            return ProviderResult(
                request,
                ProviderStatus.PROVIDER_ERROR,
                message=str(exc) or "Provider transport failed",
                metadata={"exception_type": type(exc).__name__},
            )
        status = response.status_code
        if status == 404:
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Provider returned no exact match")
        if status in (408, 504):
            return ProviderResult(request, ProviderStatus.TIMEOUT, message=f"Provider returned HTTP {status}")
        if status == 429:
            return ProviderResult(
                request,
                ProviderStatus.RATE_LIMIT,
                message="Provider rate limit reached",
                metadata={"retry_after": response.header("retry-after")},
            )
        if status < 200 or status >= 300:
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message=f"Provider returned HTTP {status}")
        return response


class UniProtProvider(JsonEvidenceProvider):
    """Exact accession lookup from UniProtKB; no search or bulk endpoint."""

    provider_name = "uniprot"
    operation = "protein_record"

    def __init__(self, *, transport: JsonTransport | None = None, enabled: bool = False) -> None:
        super().__init__(
            self.provider_name,
            (ProviderCapability(self.operation, SUPPORTED_TAXA, "Exact UniProtKB protein record"),),
            transport=transport,
            enabled=enabled,
        )

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        accession = request.canonical_query.upper()
        if not _UNIPROT_ACCESSION.fullmatch(accession):
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Canonical query is not an exact UniProt accession")
        outcome = self._get(
            request,
            timeout,
            url=f"https://rest.uniprot.org/uniprotkb/{quote(accession, safe='-')}.json",
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        payload = outcome.data
        if not isinstance(payload, Mapping):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="UniProt response root is not an object")
        organism_taxon = payload.get("organism", {}).get("taxonId") if isinstance(payload.get("organism"), Mapping) else None
        returned_accession = _text(payload.get("primaryAccession"))
        requested_base = accession.split("-", 1)[0]
        if returned_accession != requested_base:
            return ProviderResult(
                request,
                ProviderStatus.NO_MATCH,
                message="UniProt response accession does not match the exact request",
                metadata={"returned_accession": returned_accession},
            )
        if organism_taxon != request.taxon_id:
            return ProviderResult(
                request,
                ProviderStatus.NO_MATCH,
                message="UniProt record taxon does not match the requested taxon",
                metadata={"returned_taxon_id": organism_taxon},
            )

        description = payload.get("proteinDescription")
        protein_name = _nested_text(description, "recommendedName", "fullName", "value")
        if protein_name is None and isinstance(description, Mapping):
            alternatives = _items(description.get("alternativeNames"))
            protein_name = next(
                (_nested_text(item, "fullName", "value") for item in alternatives if _nested_text(item, "fullName", "value")),
                None,
            )
        genes = _bounded_unique(
            (
                _nested_text(gene, "geneName", "value")
                for gene in _items(payload.get("genes"))
            ),
            20,
        )
        functions: list[str] = []
        locations: list[str] = []
        for comment in _items(payload.get("comments")):
            comment_type = _text(comment.get("commentType"))
            if comment_type == "FUNCTION":
                functions.extend(
                    _bounded_unique((_text(item.get("value")) for item in _items(comment.get("texts"))), 20)
                )
            elif comment_type == "SUBCELLULAR LOCATION":
                for location in _items(comment.get("subcellularLocations")):
                    label = _nested_text(location, "location", "value")
                    if label:
                        locations.append(label)
        cross_references = []
        for item in _items(payload.get("uniProtKBCrossReferences"))[:100]:
            database = _text(item.get("database"))
            identifier = _text(item.get("id"))
            if database and identifier:
                cross_references.append({"database": database, "identifier": identifier})
        retrieved_at = _utc_now()
        record = {
            "accession": _text(payload.get("primaryAccession")) or accession,
            "entry_type": _text(payload.get("entryType")),
            "protein_name": protein_name,
            "genes": genes,
            "organism": {
                "taxon_id": organism_taxon,
                "scientific_name": _nested_text(payload, "organism", "scientificName"),
            },
            "function_statements": _bounded_unique(functions, 20),
            "subcellular_locations": _bounded_unique(locations, 30),
            "cross_references": cross_references,
            "provenance": _provenance(
                provider=self.provider_name,
                source="UniProtKB",
                source_identifier=accession,
                response=outcome,
                retrieved_at=retrieved_at,
            ),
        }
        return ProviderResult(
            request,
            ProviderStatus.AVAILABLE,
            data={"records": [record]},
            provider_version=_source_version(outcome),
            metadata={"retrieved_at": retrieved_at, "record_count": 1},
        )


class InterProProvider(JsonEvidenceProvider):
    """Bounded InterPro entries matching one exact UniProt accession."""

    provider_name = "interpro"
    operation = "protein_entries"

    def __init__(self, *, transport: JsonTransport | None = None, enabled: bool = False) -> None:
        super().__init__(
            self.provider_name,
            (ProviderCapability(self.operation, SUPPORTED_TAXA, "InterPro entries for an exact UniProt protein"),),
            transport=transport,
            enabled=enabled,
        )

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        accession = request.canonical_query.upper()
        if not _UNIPROT_ACCESSION.fullmatch(accession):
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Canonical query is not an exact UniProt accession")
        limit = min(max(int(request.params.get("page_size", 20)), 1), 100)
        outcome = self._get(
            request,
            timeout,
            url=f"https://www.ebi.ac.uk/interpro/api/entry/interpro/protein/uniprot/{quote(accession, safe='-')}/",
            params={"page_size": limit},
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        payload = outcome.data
        if not isinstance(payload, Mapping) or not isinstance(payload.get("results"), list):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="InterPro response has an invalid schema")
        retrieved_at = _utc_now()
        records = []
        for item in _items(payload.get("results"))[:limit]:
            metadata = item.get("metadata") if isinstance(item.get("metadata"), Mapping) else item
            identifier = _text(metadata.get("accession"))
            if not identifier:
                continue
            source_database = metadata.get("source_database")
            if isinstance(source_database, Mapping):
                source_database = source_database.get("name")
            records.append({
                "entry_id": identifier,
                "name": _text(metadata.get("name")),
                "entry_type": _text(metadata.get("type")),
                "source_database": _text(source_database) or "InterPro",
                "protein_accession": accession,
                "provenance": _provenance(
                    provider=self.provider_name,
                    source="InterPro",
                    source_identifier=identifier,
                    response=outcome,
                    retrieved_at=retrieved_at,
                ),
            })
        status = ProviderStatus.AVAILABLE if records else ProviderStatus.NO_MATCH
        return ProviderResult(
            request,
            status,
            data={"records": records, "total_available": payload.get("count")} if records else None,
            provider_version=_source_version(outcome),
            metadata={"retrieved_at": retrieved_at, "record_count": len(records), "page_followed": False},
        )


class QuickGOProvider(JsonEvidenceProvider):
    """Exact GO term or bounded UniProt GO annotation lookup."""

    provider_name = "quickgo"
    annotation_operation = "protein_go_annotations"
    term_operation = "go_term"

    def __init__(self, *, transport: JsonTransport | None = None, enabled: bool = False) -> None:
        super().__init__(
            self.provider_name,
            (
                ProviderCapability(self.annotation_operation, SUPPORTED_TAXA, "GO annotations for an exact UniProt protein"),
                ProviderCapability(self.term_operation, SUPPORTED_TAXA, "Exact GO ontology term"),
            ),
            transport=transport,
            enabled=enabled,
        )

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        if request.operation == self.term_operation:
            return self._term(request, timeout)
        return self._annotations(request, timeout)

    def _term(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        go_id = request.canonical_query.upper()
        if not _GO_ID.fullmatch(go_id):
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Canonical query is not an exact GO identifier")
        outcome = self._get(
            request,
            timeout,
            url=f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{quote(go_id, safe=':')}",
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        payload = outcome.data
        results = _items(payload.get("results")) if isinstance(payload, Mapping) else ()
        exact = next((item for item in results if _text(item.get("id")) == go_id), None)
        if exact is None:
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="QuickGO returned no exact GO term")
        retrieved_at = _utc_now()
        record = {
            "go_id": go_id,
            "name": _text(exact.get("name")),
            "aspect": _text(exact.get("aspect")),
            "definition": _nested_text(exact, "definition", "text") or _text(exact.get("definition")),
            "obsolete": bool(exact.get("isObsolete", False)),
            "provenance": _provenance(
                provider=self.provider_name,
                source="Gene Ontology via QuickGO",
                source_identifier=go_id,
                response=outcome,
                retrieved_at=retrieved_at,
            ),
        }
        return ProviderResult(request, ProviderStatus.AVAILABLE, data={"records": [record]}, metadata={"retrieved_at": retrieved_at, "record_count": 1})

    def _annotations(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        accession = request.canonical_query.upper()
        if not _UNIPROT_ACCESSION.fullmatch(accession):
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Canonical query is not an exact UniProt accession")
        limit = min(max(int(request.params.get("limit", 25)), 1), 100)
        outcome = self._get(
            request,
            timeout,
            url="https://www.ebi.ac.uk/QuickGO/services/annotation/search",
            params={"geneProductId": f"UniProtKB:{accession}", "taxonId": request.taxon_id, "limit": limit},
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        payload = outcome.data
        if not isinstance(payload, Mapping) or not isinstance(payload.get("results"), list):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="QuickGO response has an invalid schema")
        retrieved_at = _utc_now()
        records = []
        rejected_taxa = 0
        for item in _items(payload.get("results"))[:limit]:
            gene_product_id = _text(item.get("geneProductId"))
            if gene_product_id and gene_product_id.casefold() != f"UniProtKB:{accession}".casefold():
                continue
            returned_taxon = item.get("taxonId")
            try:
                returned_taxon = int(returned_taxon)
            except (TypeError, ValueError):
                returned_taxon = None
            if returned_taxon != request.taxon_id:
                rejected_taxa += 1
                continue
            go_id = _text(item.get("goId"))
            if not go_id:
                continue
            records.append({
                "protein_accession": accession,
                "go_id": go_id,
                "go_name": _text(item.get("goName")),
                "go_aspect": _text(item.get("goAspect")),
                "evidence_code": _text(item.get("evidenceCode")),
                "qualifier": _text(item.get("qualifier")),
                "reference": _text(item.get("reference")),
                "assigned_by": _text(item.get("assignedBy")),
                "taxon_id": returned_taxon,
                "provenance": _provenance(
                    provider=self.provider_name,
                    source="Gene Ontology Annotation via QuickGO",
                    source_identifier=f"{accession}:{go_id}",
                    response=outcome,
                    retrieved_at=retrieved_at,
                ),
            })
        status = ProviderStatus.AVAILABLE if records else ProviderStatus.NO_MATCH
        return ProviderResult(
            request,
            status,
            data={"records": records, "total_available": payload.get("numberOfHits")} if records else None,
            metadata={"retrieved_at": retrieved_at, "record_count": len(records), "cross_taxon_records_rejected": rejected_taxa},
        )


class ReactomeProvider(JsonEvidenceProvider):
    """Bounded pathway mappings for one exact, species-specific accession."""

    provider_name = "reactome"
    operation = "protein_pathways"

    def __init__(self, *, transport: JsonTransport | None = None, enabled: bool = False) -> None:
        super().__init__(
            self.provider_name,
            (ProviderCapability(self.operation, SUPPORTED_TAXA, "Reactome pathways mapped to an exact UniProt protein"),),
            transport=transport,
            enabled=enabled,
        )

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        accession = request.canonical_query.upper()
        if not _UNIPROT_ACCESSION.fullmatch(accession):
            return ProviderResult(request, ProviderStatus.NO_MATCH, message="Canonical query is not an exact UniProt accession")
        limit = min(max(int(request.params.get("limit", 50)), 1), 100)
        outcome = self._get(
            request,
            timeout,
            url=f"https://reactome.org/ContentService/data/mapping/UniProt/{quote(accession, safe='-')}/pathways",
        )
        if isinstance(outcome, ProviderResult):
            return outcome
        if not isinstance(outcome.data, list):
            return ProviderResult(request, ProviderStatus.PROVIDER_ERROR, message="Reactome response root is not a list")
        retrieved_at = _utc_now()
        records = []
        rejected_taxa = 0
        expected_species = {9606: "Homo sapiens", 10090: "Mus musculus"}[request.taxon_id]
        for item in _items(outcome.data)[:limit]:
            stable_id = _text(item.get("stId")) or _text(item.get("stIdVersion"))
            if not stable_id:
                continue
            species_name = _text(item.get("speciesName"))
            if species_name is not None and species_name.casefold() != expected_species.casefold():
                rejected_taxa += 1
                continue
            records.append({
                "pathway_id": stable_id,
                "database_id": item.get("dbId"),
                "display_name": _text(item.get("displayName")),
                "species_name": species_name,
                "is_inferred": item.get("isInferred") if isinstance(item.get("isInferred"), bool) else None,
                "is_disease": item.get("isInDisease") if isinstance(item.get("isInDisease"), bool) else None,
                "protein_accession": accession,
                "requested_taxon_id": request.taxon_id,
                "provenance": _provenance(
                    provider=self.provider_name,
                    source="Reactome",
                    source_identifier=stable_id,
                    response=outcome,
                    retrieved_at=retrieved_at,
                ),
            })
        status = ProviderStatus.AVAILABLE if records else ProviderStatus.NO_MATCH
        return ProviderResult(
            request,
            status,
            data={"records": records} if records else None,
            provider_version=_source_version(outcome),
            metadata={
                "retrieved_at": retrieved_at,
                "record_count": len(records),
                "cross_taxon_records_rejected": rejected_taxa,
                "taxon_binding": "exact species-scoped UniProt accession supplied by the external biology router",
            },
        )


__all__ = [
    "INTERPRO_SCHEMA_VERSION",
    "JsonEvidenceProvider",
    "QUICKGO_SCHEMA_VERSION",
    "REACTOME_SCHEMA_VERSION",
    "UNIPROT_SCHEMA_VERSION",
    "InterProProvider",
    "QuickGOProvider",
    "ReactomeProvider",
    "UniProtProvider",
]
