from __future__ import annotations

from dataclasses import dataclass
from enum import Enum

from .orthology import OrthologRecord, OrthologyClass, OrthologyIndex


class SuitabilityStatus(str, Enum):
    DIRECT_ORTHOLOG_AVAILABLE = "DIRECT_ORTHOLOG_AVAILABLE"
    ORTHOLOG_EXISTS_NOT_IN_GRAPH = "ORTHOLOG_EXISTS_NOT_IN_GRAPH"
    ORTHOLOG_EXISTS_CONTEXT_UNAVAILABLE = "ORTHOLOG_EXISTS_CONTEXT_UNAVAILABLE"
    MULTIPLE_ORTHOLOGS = "MULTIPLE_ORTHOLOGS"
    NO_DIRECT_ORTHOLOG = "NO_DIRECT_ORTHOLOG"
    FUNCTIONAL_PROXY_AVAILABLE = "FUNCTIONAL_PROXY_AVAILABLE"
    NO_DEFENSIBLE_PROXY = "NO_DEFENSIBLE_PROXY"
    UNRESOLVED = "UNRESOLVED"


@dataclass(frozen=True)
class SuitabilityResult:
    human_gene: str
    status: SuitabilityStatus
    orthologs: tuple[OrthologRecord, ...]
    mouse_protein_ids: tuple[str, ...]
    context: str | None
    reason: str


def check_mouse_suitability(
    human_gene: str,
    orthology: OrthologyIndex,
    mouse_symbol_to_protein: dict[str, str],
    *,
    mouse_graph_nodes: set[str] | None = None,
    context_nodes: set[str] | None = None,
    context: str | None = None,
) -> SuitabilityResult:
    records = tuple(record for record in orthology.lookup(human_gene) if record.mouse_gene)
    if not records:
        known = orthology.lookup(human_gene)
        status = SuitabilityStatus.NO_DIRECT_ORTHOLOG if known else SuitabilityStatus.UNRESOLVED
        return SuitabilityResult(human_gene, status, tuple(known), (), context, "No authoritative direct Mouse ortholog")
    if len({record.mouse_gene for record in records}) > 1 or any(
        record.orthology_class in {OrthologyClass.ONE_TO_MANY, OrthologyClass.MANY_TO_MANY}
        for record in records
    ):
        return SuitabilityResult(human_gene, SuitabilityStatus.MULTIPLE_ORTHOLOGS, records, (), context, "Multiple orthologs require an explicit user choice")

    mouse_gene = records[0].mouse_gene or ""
    lookup = {str(key).casefold(): str(value) for key, value in mouse_symbol_to_protein.items()}
    protein = lookup.get(mouse_gene.casefold())
    if not protein or (mouse_graph_nodes is not None and protein not in mouse_graph_nodes):
        return SuitabilityResult(human_gene, SuitabilityStatus.ORTHOLOG_EXISTS_NOT_IN_GRAPH, records, tuple(filter(None, [protein])), context, "Ortholog exists but its protein is absent from the Mouse graph")
    if context_nodes is not None and protein not in context_nodes:
        return SuitabilityResult(human_gene, SuitabilityStatus.ORTHOLOG_EXISTS_CONTEXT_UNAVAILABLE, records, (protein,), context, "Ortholog exists but is unavailable in the selected Mouse context")
    return SuitabilityResult(human_gene, SuitabilityStatus.DIRECT_ORTHOLOG_AVAILABLE, records, (protein,), context, "Direct one-to-one ortholog and required Mouse calculation inputs are available")
