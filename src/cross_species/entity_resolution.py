from __future__ import annotations

from dataclasses import dataclass

from .orthology import OrthologyIndex


@dataclass(frozen=True)
class MouseEntityResolution:
    query: str
    official_symbol: str | None
    protein_id: str | None
    status: str
    mapping_type: str
    source: str
    ambiguity: tuple[str, ...] = ()


def resolve_mouse_entity(
    query: str,
    mouse_symbol_to_protein: dict[str, str],
    orthology: OrthologyIndex | None = None,
) -> MouseEntityResolution:
    raw = str(query or "").strip()
    if not raw:
        return MouseEntityResolution(raw, None, None, "UNRESOLVED", "EMPTY", "none")
    if raw.startswith("ENSP") and not raw.startswith("ENSMUSP"):
        return MouseEntityResolution(raw, None, None, "SPECIES_MISMATCH", "REJECTED_HUMAN_PROTEIN_ID", "species boundary")
    normalized = {str(symbol).casefold(): (str(symbol), str(protein)) for symbol, protein in mouse_symbol_to_protein.items()}
    if raw.startswith("ENSMUSP"):
        for symbol, protein in mouse_symbol_to_protein.items():
            if str(protein) == raw:
                return MouseEntityResolution(raw, str(symbol), raw, "RESOLVED", "EXACT_MOUSE_PROTEIN_ID", "STRING protein.info")
        return MouseEntityResolution(raw, None, None, "UNRESOLVED", "UNKNOWN_MOUSE_PROTEIN_ID", "STRING protein.info")
    direct = normalized.get(raw.casefold())
    if direct:
        return MouseEntityResolution(raw, direct[0], direct[1], "RESOLVED", "EXACT_MOUSE_SYMBOL", "STRING protein.info")
    if orthology is not None:
        records = tuple(record for record in orthology.lookup(raw) if record.mouse_gene)
        candidates = tuple(sorted({record.mouse_gene for record in records if record.mouse_gene}))
        if len(candidates) > 1:
            return MouseEntityResolution(raw, None, None, "AMBIGUOUS", "MULTIPLE_ORTHOLOGS", str(orthology.metadata.get("source", "Ensembl Compara")), candidates)
        if len(candidates) == 1:
            mapped = normalized.get(candidates[0].casefold())
            if mapped and orthology.one_to_one(raw) is not None:
                return MouseEntityResolution(raw, mapped[0], mapped[1], "RESOLVED", "AUTHORITATIVE_ONE_TO_ONE_ORTHOLOG_ALIAS", str(orthology.metadata.get("source", "Ensembl Compara")))
    return MouseEntityResolution(raw, None, None, "UNRESOLVED", "NO_EXACT_MAPPING", "STRING protein.info / Ensembl Compara")
