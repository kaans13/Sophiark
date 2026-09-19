from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class FunctionalProxyCandidate:
    mouse_gene: str
    same_family: bool
    pathway_overlap: tuple[str, ...]
    go_mf_overlap: tuple[str, ...]
    go_bp_overlap: tuple[str, ...]
    compartment_overlap: tuple[str, ...]
    tissue_supported: bool
    evidence_source: str
    automatically_selected: bool = False


def discover_functional_proxies(
    human_annotations: dict,
    mouse_annotations: dict[str, dict],
    *,
    limit: int = 10,
) -> tuple[FunctionalProxyCandidate, ...]:
    """Rank explanatory annotation overlaps; never substitute a target."""

    def values(record: dict, key: str) -> set[str]:
        return {str(item).strip() for item in record.get(key, ()) if str(item).strip()}

    rows: list[tuple[tuple[int, ...], FunctionalProxyCandidate]] = []
    human_family = values(human_annotations, "family")
    for mouse_gene, annotation in mouse_annotations.items():
        same_family = bool(human_family & values(annotation, "family"))
        pathway = tuple(sorted(values(human_annotations, "pathway") & values(annotation, "pathway")))
        mf = tuple(sorted(values(human_annotations, "go_mf") & values(annotation, "go_mf")))
        bp = tuple(sorted(values(human_annotations, "go_bp") & values(annotation, "go_bp")))
        compartment = tuple(sorted(values(human_annotations, "compartment") & values(annotation, "compartment")))
        tissue = bool(annotation.get("tissue_supported", False))
        evidence_count = int(same_family) + len(pathway) + len(mf) + len(bp) + len(compartment) + int(tissue)
        if evidence_count == 0:
            continue
        candidate = FunctionalProxyCandidate(
            mouse_gene=str(mouse_gene), same_family=same_family,
            pathway_overlap=pathway, go_mf_overlap=mf, go_bp_overlap=bp,
            compartment_overlap=compartment, tissue_supported=tissue,
            evidence_source=str(annotation.get("evidence_source", "unspecified")),
        )
        rank = (int(same_family), len(pathway), len(mf), len(bp), len(compartment), int(tissue))
        rows.append((rank, candidate))
    rows.sort(key=lambda item: (item[0], item[1].mouse_gene), reverse=True)
    return tuple(candidate for _, candidate in rows[: max(0, int(limit))])
