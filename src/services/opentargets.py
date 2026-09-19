"""Small, UI-neutral Open Targets adapter used by the analysis setup page."""

from __future__ import annotations

import requests

API_URL = "https://api.platform.opentargets.org/api/v4/graphql"


class OpenTargetsUnavailable(RuntimeError):
    pass


def _post(query: str, variables: dict, timeout: int = 15) -> dict:
    try:
        response = requests.post(API_URL, json={"query": query, "variables": variables}, timeout=timeout)
        response.raise_for_status()
        payload = response.json()
    except (requests.RequestException, ValueError) as exc:
        raise OpenTargetsUnavailable("Open Targets hizmetine şu anda erişilemiyor.") from exc
    if payload.get("errors"):
        raise OpenTargetsUnavailable("Open Targets sorgusu tamamlanamadı.")
    return payload.get("data", {})


def search_diseases(term: str) -> list[dict[str, str]]:
    query = """query($term: String!) { search(queryString:$term, entityNames:[\"disease\"]) { hits { id name } } }"""
    hits = _post(query, {"term": term}, timeout=10).get("search", {}).get("hits", [])
    return [{"id": str(hit["id"]), "name": str(hit["name"])} for hit in hits if hit.get("id") and hit.get("name")]


def disease_genes(efo_id: str, limit: int = 200) -> list[str]:
    query = """query($id:String!,$page:Pagination!){ disease(efoId:$id){ associatedTargets(page:$page){ rows{ target{ approvedSymbol } score } } } }"""
    disease = _post(query, {"id": efo_id, "page": {"index": 0, "size": limit}}).get("disease") or {}
    rows = disease.get("associatedTargets", {}).get("rows", [])
    return list(dict.fromkeys(str(row.get("target", {}).get("approvedSymbol", "")).strip() for row in rows if row.get("target", {}).get("approvedSymbol")))
