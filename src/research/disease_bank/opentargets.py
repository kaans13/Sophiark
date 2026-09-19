"""Open Targets GraphQL adapter used only after an explicit Disease Bank request."""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Callable, Mapping

import requests

from ..config import ProviderTimeout
from ..provider_runtime import EvidenceProvider, ProviderCapability, ProviderRequest, ProviderResult, ProviderStatus
from .models import Disease, DiseaseGeneAssociation
from .source_registry import OPEN_TARGETS_SOURCE


OPEN_TARGETS_GRAPHQL_URL = "https://api.platform.opentargets.org/api/v4/graphql"
_SEARCH_QUERY = """query($term: String!) { search(queryString: $term, entityNames: [\"disease\"]) { hits { id name } } }"""
_ASSOCIATIONS_QUERY = """query($id: String!, $page: Pagination!) {
 disease(efoId: $id) {
   id name synonyms { relation terms } dbXRefs parents { id name }
   associatedTargets(page: $page) {
     count
     rows { score datatypeScores { id score } target { id approvedSymbol approvedName } }
   }
 }
}"""


def _utc_now() -> str:
    return datetime.now(timezone.utc).isoformat()


def _string_items(value: object) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    return tuple(dict.fromkeys(str(item).strip() for item in value if str(item).strip()))


def _mondo_id(disease_id: object, xrefs: object) -> str | None:
    for item in (str(disease_id).strip(), *_string_items(xrefs)):
        normalized = item.strip()
        if normalized.upper().startswith("MONDO_"):
            return f"MONDO:{normalized.split('_', 1)[1]}"
        if normalized.upper().startswith("MONDO:"):
            return normalized.upper()
    return None


def _synonyms(value: object) -> tuple[str, ...]:
    if not isinstance(value, (list, tuple)):
        return ()
    values: list[str] = []
    for item in value:
        if isinstance(item, Mapping):
            values.extend(_string_items(item.get("terms")))
        else:
            values.extend(_string_items((item,)))
    return tuple(dict.fromkeys(values))


class OpenTargetsDiseaseProvider(EvidenceProvider):
    """Small v4 adapter with no import-time network activity."""

    def __init__(
        self,
        *,
        enabled: bool = True,
        post: Callable[..., Any] | None = None,
        endpoint: str = OPEN_TARGETS_GRAPHQL_URL,
    ) -> None:
        super().__init__(
            "open_targets_disease",
            (
                ProviderCapability("disease_search", (9606,), "Search disease labels."),
                ProviderCapability("disease_associations", (9606,), "Retrieve disease-target associations."),
            ),
            provider_version=OPEN_TARGETS_SOURCE.source_version,
            enabled=enabled,
        )
        self._post = post or requests.post
        self._endpoint = str(endpoint)

    def _graphql(self, query: str, variables: Mapping[str, object], timeout: ProviderTimeout) -> Mapping[str, object]:
        response = self._post(
            self._endpoint,
            json={"query": query, "variables": dict(variables)},
            timeout=(timeout.connect_seconds, timeout.read_seconds),
        )
        response.raise_for_status()
        payload = response.json()
        if not isinstance(payload, Mapping):
            raise ValueError("Open Targets returned a non-object response")
        if payload.get("errors"):
            raise ValueError("Open Targets GraphQL returned errors")
        data = payload.get("data")
        if not isinstance(data, Mapping):
            raise ValueError("Open Targets response has no data object")
        return data

    def execute(self, request: ProviderRequest, timeout: ProviderTimeout) -> ProviderResult:
        try:
            if request.operation == "disease_search":
                data = self._graphql(_SEARCH_QUERY, {"term": request.canonical_query}, timeout)
                hits = data.get("search", {}).get("hits", []) if isinstance(data.get("search"), Mapping) else []
                records = tuple(
                    {"disease_id": str(item.get("id", "")).strip(), "disease_name": str(item.get("name", "")).strip()}
                    for item in hits if isinstance(item, Mapping) and item.get("id") and item.get("name")
                )
                return ProviderResult(
                    request=request,
                    status=ProviderStatus.AVAILABLE if records else ProviderStatus.NO_MATCH,
                    data={"diseases": records},
                    provider_version=self.provider_version,
                    metadata={"source_key": OPEN_TARGETS_SOURCE.source_key, "retrieved_at": _utc_now()},
                )

            if request.operation == "disease_associations":
                size = min(max(int(request.params.get("limit", 200)), 1), 500)
                data = self._graphql(_ASSOCIATIONS_QUERY, {"id": request.canonical_query, "page": {"index": 0, "size": size}}, timeout)
                disease = data.get("disease")
                if not isinstance(disease, Mapping) or not disease.get("id"):
                    return ProviderResult(request=request, status=ProviderStatus.NO_MATCH, data=None, provider_version=self.provider_version)
                rows = disease.get("associatedTargets", {}).get("rows", []) if isinstance(disease.get("associatedTargets"), Mapping) else []
                associations: list[dict[str, object]] = []
                for row in rows:
                    if not isinstance(row, Mapping) or not isinstance(row.get("target"), Mapping):
                        continue
                    target = row["target"]
                    ensembl = str(target.get("id", "")).strip()
                    symbol = str(target.get("approvedSymbol", "")).strip()
                    if not (ensembl or symbol):
                        continue
                    associations.append({
                        "disease_id": str(disease["id"]),
                        "gene_symbol": symbol,
                        "ensembl_gene_id": ensembl,
                        "target_name": str(target.get("approvedName", "")).strip(),
                        "association_score": row.get("score"),
                        "evidence": {"datatype_scores": row.get("datatypeScores", [])},
                    })
                parents = tuple(
                    {"id": str(parent.get("id", "")).strip(), "name": str(parent.get("name", "")).strip()}
                    for parent in disease.get("parents", []) if isinstance(parent, Mapping)
                )
                return ProviderResult(
                    request=request,
                    status=ProviderStatus.AVAILABLE,
                    data={
                        "disease": {
                            "disease_id": str(disease["id"]), "disease_name": str(disease.get("name", disease["id"])),
                            "synonyms": _synonyms(disease.get("synonyms")), "mondo_id": _mondo_id(disease["id"], disease.get("dbXRefs")),
                            "parents": parents,
                        },
                        "associations": tuple(associations),
                    },
                    provider_version=self.provider_version,
                    metadata={"source_key": OPEN_TARGETS_SOURCE.source_key, "retrieved_at": _utc_now()},
                )
            return ProviderResult(request=request, status=ProviderStatus.DISABLED, data=None, message="Unsupported disease operation.")
        except requests.Timeout as error:
            return ProviderResult(request=request, status=ProviderStatus.TIMEOUT, data=None, message=str(error), provider_version=self.provider_version)
        except requests.RequestException as error:
            return ProviderResult(request=request, status=ProviderStatus.OFFLINE, data=None, message=str(error), provider_version=self.provider_version)
        except (TypeError, ValueError, KeyError) as error:
            return ProviderResult(request=request, status=ProviderStatus.PROVIDER_ERROR, data=None, message=str(error), provider_version=self.provider_version)


def disease_from_payload(value: Mapping[str, object]) -> Disease:
    return Disease(
        disease_id=value["disease_id"], disease_name=value["disease_name"], synonyms=_string_items(value.get("synonyms")),
        parents=tuple(str(parent.get("name") or parent.get("id")) for parent in value.get("parents", ()) if isinstance(parent, Mapping)),
        mondo_id=value.get("mondo_id"),
    )


def associations_from_payload(values: object, *, retrieved_at: str) -> tuple[DiseaseGeneAssociation, ...]:
    if not isinstance(values, (list, tuple)):
        return ()
    return tuple(
        DiseaseGeneAssociation(
            disease_id=item["disease_id"], gene_symbol=item.get("gene_symbol", ""), ensembl_gene_id=item.get("ensembl_gene_id", ""),
            target_name=item.get("target_name", ""), association_score=item.get("association_score"),
            evidence=item.get("evidence", {}), retrieved_at=retrieved_at,
        ) for item in values if isinstance(item, Mapping)
    )


__all__ = ["OPEN_TARGETS_GRAPHQL_URL", "OpenTargetsDiseaseProvider", "associations_from_payload", "disease_from_payload"]
