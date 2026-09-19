"""Species-aware, ambiguity-safe entity resolution for Research Context.

The resolver is a read-only consumer of project mappings.  Loaders are
explicit: importing this module never reads a project dataset and never makes
an external request.
"""

from __future__ import annotations

import csv
import gzip
import re
import unicodedata
from dataclasses import dataclass
from pathlib import Path
from types import MappingProxyType
from typing import Callable, Iterable, Mapping

from .models import (
    AliasType,
    Entity,
    EntityAlias,
    EntityResolution,
    EntityType,
    ResolutionStatus,
)


HUMAN_TAXON_ID = 9606
MOUSE_TAXON_ID = 10090
_PROJECT_ROOT = Path(__file__).resolve().parents[2]
_MULTIVALUE_SEPARATOR = re.compile(r"[|;,]")


def _lookup_key(value: object) -> str:
    return unicodedata.normalize("NFKC", str(value or "")).strip().casefold()


def _entity_type(value: EntityType | str) -> EntityType:
    if isinstance(value, EntityType):
        return value
    return EntityType(str(value).strip().casefold())


def _alias_type(value: AliasType | str) -> AliasType:
    if isinstance(value, AliasType):
        return value
    normalized = str(value).strip().casefold()
    try:
        return AliasType(normalized)
    except ValueError:
        return AliasType[normalized.upper()]


def _split_values(value: object) -> tuple[str, ...]:
    raw = str(value or "").strip()
    if not raw:
        return ()
    return tuple(item.strip() for item in _MULTIVALUE_SEPARATOR.split(raw) if item.strip())


@dataclass(frozen=True, slots=True)
class _IndexedMatch:
    entity: Entity
    alias: EntityAlias | None = None


@dataclass(frozen=True, slots=True)
class ProjectEntityIndex:
    """Immutable lookup indexes for one taxon and one entity type."""

    taxon_id: int
    entity_type: EntityType
    canonical: Mapping[str, tuple[_IndexedMatch, ...]]
    project_mappings: Mapping[str, tuple[_IndexedMatch, ...]]
    symbols: Mapping[str, tuple[_IndexedMatch, ...]]
    aliases: Mapping[str, tuple[_IndexedMatch, ...]]

    @property
    def entity_count(self) -> int:
        return len({match.entity.logical_key for matches in self.canonical.values() for match in matches})


@dataclass(frozen=True, slots=True)
class ExternalCandidateSet:
    """Explicit, provenance-bearing contract for cache/provider resolution.

    Zero candidates is an explicit negative result.  Multiple candidates stay
    ambiguous; the resolver never chooses one silently.
    """

    candidates: tuple[Entity, ...]
    source: str
    matched_alias: EntityAlias | None = None
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "candidates", tuple(self.candidates))
        if not str(self.source).strip():
            raise ValueError("ExternalCandidateSet.source is required for provenance.")


ExternalCacheKey = tuple[int, EntityType, str]
ExternalResolver = Callable[..., ExternalCandidateSet]


def external_cache_key(
    query: object, taxon_id: int, entity_type: EntityType | str = EntityType.PROTEIN
) -> ExternalCacheKey:
    """Build a species- and type-scoped cache key."""

    return int(taxon_id), _entity_type(entity_type), _lookup_key(query)


class _IndexBuilder:
    def __init__(self, taxon_id: int, entity_type: EntityType) -> None:
        self.taxon_id = int(taxon_id)
        self.entity_type = entity_type
        self._symbols: dict[str, set[str]] = {}
        self._display_names: dict[str, set[str]] = {}
        self._project: list[tuple[str, str, AliasType, str]] = []
        self._aliases: list[tuple[str, str, AliasType, str]] = []
        self._extra_symbols: list[tuple[str, str, str]] = []

    def add_entity(self, canonical_id: object, *, symbol: object = None, display_name: object = None) -> None:
        canonical = str(canonical_id or "").strip()
        if not canonical:
            return
        self._symbols.setdefault(canonical, set())
        self._display_names.setdefault(canonical, set())
        if str(symbol or "").strip():
            self._symbols[canonical].add(str(symbol).strip())
        if str(display_name or "").strip():
            self._display_names[canonical].add(str(display_name).strip())

    def add_project_mapping(
        self,
        lookup_value: object,
        canonical_id: object,
        *,
        alias_type: AliasType | str,
        source: str,
    ) -> None:
        value = str(lookup_value or "").strip()
        canonical = str(canonical_id or "").strip()
        if value and canonical:
            self._project.append((value, canonical, _alias_type(alias_type), source))

    def add_symbol(self, symbol: object, canonical_id: object, *, source: str) -> None:
        value = str(symbol or "").strip()
        canonical = str(canonical_id or "").strip()
        if value and canonical:
            self._extra_symbols.append((value, canonical, source))

    def add_alias(
        self,
        alias_value: object,
        canonical_id: object,
        *,
        alias_type: AliasType | str,
        source: str,
    ) -> None:
        value = str(alias_value or "").strip()
        canonical = str(canonical_id or "").strip()
        if value and canonical:
            self._aliases.append((value, canonical, _alias_type(alias_type), source))

    def build(self) -> ProjectEntityIndex:
        entities: dict[str, Entity] = {}
        for canonical, symbols in self._symbols.items():
            names = self._display_names.get(canonical, set())
            entities[canonical] = Entity(
                taxon_id=self.taxon_id,
                entity_type=self.entity_type,
                canonical_id=canonical,
                symbol=next(iter(symbols)) if len(symbols) == 1 else None,
                display_name=next(iter(names)) if len(names) == 1 else None,
            )

        canonical_map: dict[str, list[_IndexedMatch]] = {}
        project_map: dict[str, list[_IndexedMatch]] = {}
        symbol_map: dict[str, list[_IndexedMatch]] = {}
        alias_map: dict[str, list[_IndexedMatch]] = {}

        for canonical, entity in entities.items():
            canonical_map.setdefault(_lookup_key(canonical), []).append(_IndexedMatch(entity))
            for symbol in self._symbols.get(canonical, set()):
                alias = EntityAlias(_alias_type("gene_symbol"), symbol, "project_symbol_mapping")
                symbol_map.setdefault(_lookup_key(symbol), []).append(_IndexedMatch(entity, alias))

        for symbol, canonical, source in self._extra_symbols:
            entity = entities.get(canonical)
            if entity is not None:
                alias = EntityAlias(_alias_type("gene_symbol"), symbol, source)
                symbol_map.setdefault(_lookup_key(symbol), []).append(_IndexedMatch(entity, alias))

        for value, canonical, alias_type, source in self._project:
            entity = entities.get(canonical)
            if entity is not None:
                alias = EntityAlias(alias_type, value, source)
                project_map.setdefault(_lookup_key(value), []).append(_IndexedMatch(entity, alias))

        for value, canonical, alias_type, source in self._aliases:
            entity = entities.get(canonical)
            if entity is not None:
                alias = EntityAlias(alias_type, value, source)
                alias_map.setdefault(_lookup_key(value), []).append(_IndexedMatch(entity, alias))

        return ProjectEntityIndex(
            taxon_id=self.taxon_id,
            entity_type=self.entity_type,
            canonical=_freeze_match_map(canonical_map),
            project_mappings=_freeze_match_map(project_map),
            symbols=_freeze_match_map(symbol_map),
            aliases=_freeze_match_map(alias_map),
        )


def _freeze_match_map(
    source: Mapping[str, Iterable[_IndexedMatch]],
) -> Mapping[str, tuple[_IndexedMatch, ...]]:
    frozen: dict[str, tuple[_IndexedMatch, ...]] = {}
    for key, matches in source.items():
        unique: dict[tuple[str, str, str, str], _IndexedMatch] = {}
        for match in matches:
            alias = match.alias
            identity = (
                match.entity.logical_key,
                str(alias.alias_type.value) if alias is not None else "",
                alias.alias_value if alias is not None else "",
                alias.source if alias is not None else "",
            )
            unique[identity] = match
        frozen[key] = tuple(sorted(unique.values(), key=lambda item: item.entity.logical_key))
    return MappingProxyType(frozen)


def empty_project_index(
    taxon_id: int, entity_type: EntityType | str = EntityType.PROTEIN
) -> ProjectEntityIndex:
    """Return an immutable empty index without touching the filesystem."""

    return _IndexBuilder(int(taxon_id), _entity_type(entity_type)).build()


def _read_csv_rows(path: Path, *, delimiter: str = ",") -> Iterable[dict[str, str]]:
    if not path.is_file():
        return ()

    def rows() -> Iterable[dict[str, str]]:
        with path.open("r", encoding="utf-8-sig", newline="") as handle:
            yield from csv.DictReader(handle, delimiter=delimiter)

    return rows()


def load_human_project_index(
    project_root: str | Path = _PROJECT_ROOT,
    *,
    symbol_path: str | Path | None = None,
    ensp_to_ensg_path: str | Path | None = None,
    hgnc_path: str | Path | None = None,
) -> ProjectEntityIndex:
    """Load existing human ENSP/ENSG/symbol/HGNC mappings read-only."""

    root = Path(project_root)
    symbols_file = Path(symbol_path) if symbol_path is not None else root / "data" / "processed" / "ensp_with_symbols.csv"
    ensg_file = Path(ensp_to_ensg_path) if ensp_to_ensg_path is not None else root / "data" / "processed" / "ensp_to_ensg_map.csv"
    hgnc_file = Path(hgnc_path) if hgnc_path is not None else root / "data" / "interpretation_reference" / "hgnc" / "hgnc_gene_groups.tsv"

    builder = _IndexBuilder(HUMAN_TAXON_ID, EntityType.PROTEIN)
    for row in _read_csv_rows(symbols_file):
        canonical = str(row.get("gene", "")).strip()
        symbol = str(row.get("Symbol", "")).strip()
        if not canonical:
            continue
        builder.add_entity(canonical, symbol=symbol)
        builder.add_alias(
            f"{HUMAN_TAXON_ID}.{canonical}", canonical,
            alias_type="string_id", source="STRING protein identifier mapping",
        )

    ensg_to_canonical: dict[str, set[str]] = {}
    for row in _read_csv_rows(ensg_file):
        canonical = str(row.get("ENSP", "")).strip()
        ensg = str(row.get("ENSG", "")).strip()
        if not canonical or not ensg:
            continue
        ensg_to_canonical.setdefault(_lookup_key(ensg), set()).add(canonical)
        builder.add_project_mapping(
            ensg, canonical, alias_type="ensembl_gene", source="ensp_to_ensg_map.csv",
        )

    for row in _read_csv_rows(hgnc_file, delimiter="\t"):
        ensembl_ids = _split_values(row.get("Ensembl gene ID", ""))
        canonical_ids = {
            canonical
            for ensembl_id in ensembl_ids
            for canonical in ensg_to_canonical.get(_lookup_key(ensembl_id), set())
        }
        if not canonical_ids:
            continue
        approved_symbol = str(row.get("Approved symbol", "")).strip()
        approved_name = str(row.get("Approved name", "")).strip()
        for canonical in canonical_ids:
            builder.add_entity(canonical, symbol=approved_symbol, display_name=approved_name)
            builder.add_symbol(approved_symbol, canonical, source="HGNC approved symbol")
            for hgnc_id in _split_values(row.get("HGNC ID", "")):
                builder.add_project_mapping(
                    hgnc_id, canonical, alias_type="synonym", source="HGNC identifier",
                )
            for ncbi_id in _split_values(row.get("NCBI Gene ID", "")):
                builder.add_project_mapping(
                    ncbi_id, canonical, alias_type="ncbi_gene", source="HGNC NCBI Gene ID",
                )
            for field_name in ("Previous symbols", "Alias symbols"):
                for alias_value in _split_values(row.get(field_name, "")):
                    builder.add_alias(
                        alias_value, canonical, alias_type="synonym", source=f"HGNC {field_name}",
                    )

    return builder.build()


def _read_string_info(path: Path) -> Iterable[dict[str, str]]:
    if not path.is_file():
        return ()

    def rows() -> Iterable[dict[str, str]]:
        with path.open("rb") as binary:
            is_gzip = binary.read(2) == b"\x1f\x8b"
        opener = gzip.open if is_gzip else open
        with opener(path, "rt", encoding="utf-8-sig", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                yield {str(key).lstrip("#"): value for key, value in row.items()}

    return rows()


def load_mouse_project_index(
    project_root: str | Path = _PROJECT_ROOT,
    *,
    string_info_path: str | Path | None = None,
) -> ProjectEntityIndex:
    """Load the local mouse STRING protein-to-symbol mapping read-only."""

    root = Path(project_root)
    info_file = Path(string_info_path) if string_info_path is not None else root / "data" / "raw" / "10090.protein.info.v12.0.txt.gz"
    builder = _IndexBuilder(MOUSE_TAXON_ID, EntityType.PROTEIN)
    prefix = f"{MOUSE_TAXON_ID}."
    for row in _read_string_info(info_file):
        string_id = str(row.get("string_protein_id", "")).strip()
        canonical = string_id[len(prefix):] if string_id.startswith(prefix) else string_id
        symbol = str(row.get("preferred_name", "")).strip()
        if not canonical:
            continue
        builder.add_entity(canonical, symbol=symbol)
        if string_id and string_id != canonical:
            builder.add_alias(
                string_id, canonical, alias_type="string_id", source="STRING mouse protein.info",
            )
    return builder.build()


class EntityResolver:
    """Resolve identifiers in the mandated order without guessing."""

    def __init__(
        self,
        indexes: Iterable[ProjectEntityIndex] = (),
        *,
        local_external_cache: Mapping[ExternalCacheKey, ExternalCandidateSet] | None = None,
        external_fallback: ExternalResolver | None = None,
    ) -> None:
        self._indexes = {
            (int(index.taxon_id), index.entity_type): index
            for index in indexes
        }
        self._local_external_cache = MappingProxyType(dict(local_external_cache or {}))
        self._external_fallback = external_fallback

    @classmethod
    def from_project_data(
        cls,
        project_root: str | Path = _PROJECT_ROOT,
        *,
        include_human: bool = True,
        include_mouse: bool = True,
        local_external_cache: Mapping[ExternalCacheKey, ExternalCandidateSet] | None = None,
        external_fallback: ExternalResolver | None = None,
    ) -> "EntityResolver":
        """Explicitly load selected local indexes; never called at import time."""

        indexes: list[ProjectEntityIndex] = []
        if include_human:
            indexes.append(load_human_project_index(project_root))
        if include_mouse:
            indexes.append(load_mouse_project_index(project_root))
        return cls(
            indexes,
            local_external_cache=local_external_cache,
            external_fallback=external_fallback,
        )

    def resolve(
        self,
        query: object,
        *,
        taxon_id: int,
        entity_type: EntityType | str = EntityType.PROTEIN,
        allow_external: bool = False,
    ) -> EntityResolution:
        """Resolve canonical → project ID → symbol → alias → cache → fallback."""

        query_text = str(query or "").strip()
        requested_type = _entity_type(entity_type)
        requested_taxon = int(taxon_id)
        if not query_text:
            return self._unresolved("<empty>", requested_taxon, requested_type, "Empty entity query.")

        key = _lookup_key(query_text)
        index = self._indexes.get((requested_taxon, requested_type))
        if index is not None:
            ordered_stages = (
                (index.canonical, ResolutionStatus.EXACT, "project_canonical"),
                (index.project_mappings, ResolutionStatus.EXACT, "project_mapping"),
                (index.symbols, ResolutionStatus.EXACT, "project_symbol"),
                (index.aliases, ResolutionStatus.ALIAS, "project_alias"),
            )
            for lookup, status, source in ordered_stages:
                matches = lookup.get(key, ())
                if matches:
                    return self._from_index_matches(
                        query_text, requested_taxon, requested_type, matches, status, source,
                    )

        cache_key = external_cache_key(query_text, requested_taxon, requested_type)
        if cache_key in self._local_external_cache:
            return self._from_external_candidates(
                query_text,
                requested_taxon,
                requested_type,
                self._local_external_cache[cache_key],
                fallback_source="local_external_cache",
            )

        if allow_external and self._external_fallback is not None:
            try:
                result = self._external_fallback(
                    query=query_text,
                    taxon_id=requested_taxon,
                    entity_type=requested_type,
                )
            except Exception as exc:  # provider failures remain isolated from core
                return self._unresolved(
                    query_text,
                    requested_taxon,
                    requested_type,
                    f"External resolution failed: {type(exc).__name__}.",
                    source="external_fallback",
                )
            if not isinstance(result, ExternalCandidateSet):
                return self._unresolved(
                    query_text,
                    requested_taxon,
                    requested_type,
                    "External resolver returned an invalid candidate contract.",
                    source="external_fallback",
                )
            return self._from_external_candidates(
                query_text,
                requested_taxon,
                requested_type,
                result,
                fallback_source="external_fallback",
            )

        return self._unresolved(
            query_text,
            requested_taxon,
            requested_type,
            "No exact species-aware entity mapping was found.",
        )

    @staticmethod
    def _from_index_matches(
        query: str,
        taxon_id: int,
        entity_type: EntityType,
        matches: Iterable[_IndexedMatch],
        status: ResolutionStatus,
        source: str,
    ) -> EntityResolution:
        ordered = tuple(matches)
        candidates = _unique_entities(match.entity for match in ordered)
        if len(candidates) != 1:
            return EntityResolution(
                query=query,
                taxon_id=taxon_id,
                entity_type=entity_type,
                status=ResolutionStatus.AMBIGUOUS,
                candidates=candidates,
                source=source,
                message=f"{len(candidates)} exact species-aware candidates matched; no candidate was selected.",
            )
        aliases = {
            (match.alias.alias_type, match.alias.alias_value, match.alias.source): match.alias
            for match in ordered
            if match.alias is not None and match.entity.logical_key == candidates[0].logical_key
        }
        matched_alias = next(iter(aliases.values())) if aliases else None
        return EntityResolution(
            query=query,
            taxon_id=taxon_id,
            entity_type=entity_type,
            status=status,
            entity=candidates[0],
            matched_alias=matched_alias,
            source=source,
        )

    @staticmethod
    def _from_external_candidates(
        query: str,
        taxon_id: int,
        entity_type: EntityType,
        result: ExternalCandidateSet,
        *,
        fallback_source: str,
    ) -> EntityResolution:
        in_scope = _unique_entities(
            candidate
            for candidate in result.candidates
            if candidate.taxon_id == taxon_id and candidate.entity_type == entity_type
        )
        rejected = len(result.candidates) - len(in_scope)
        source = result.source or fallback_source
        if rejected:
            return EntityResolution(
                query=query,
                taxon_id=taxon_id,
                entity_type=entity_type,
                status=ResolutionStatus.UNRESOLVED,
                candidates=in_scope,
                source=source,
                message="External candidates crossed the requested taxon/entity boundary and were rejected.",
            )
        if len(in_scope) > 1:
            return EntityResolution(
                query=query,
                taxon_id=taxon_id,
                entity_type=entity_type,
                status=ResolutionStatus.AMBIGUOUS,
                candidates=in_scope,
                matched_alias=result.matched_alias,
                source=source,
                message=result.message or "External source returned multiple candidates; none was selected.",
            )
        if len(in_scope) == 1:
            return EntityResolution(
                query=query,
                taxon_id=taxon_id,
                entity_type=entity_type,
                status=ResolutionStatus.EXTERNAL_RESOLVED,
                entity=in_scope[0],
                matched_alias=result.matched_alias,
                source=source,
                message=result.message,
            )
        return EntityResolver._unresolved(
            query,
            taxon_id,
            entity_type,
            result.message or "External source has an explicit negative result for this entity.",
            source=source,
        )

    @staticmethod
    def _unresolved(
        query: str,
        taxon_id: int,
        entity_type: EntityType,
        message: str,
        *,
        source: str | None = None,
    ) -> EntityResolution:
        return EntityResolution(
            query=query,
            taxon_id=taxon_id,
            entity_type=entity_type,
            status=ResolutionStatus.UNRESOLVED,
            source=source,
            message=message,
        )


def _unique_entities(entities: Iterable[Entity]) -> tuple[Entity, ...]:
    unique = {entity.logical_key: entity for entity in entities}
    return tuple(unique[key] for key in sorted(unique))


__all__ = [
    "EntityResolver",
    "ExternalCacheKey",
    "ExternalCandidateSet",
    "ExternalResolver",
    "HUMAN_TAXON_ID",
    "MOUSE_TAXON_ID",
    "ProjectEntityIndex",
    "empty_project_index",
    "external_cache_key",
    "load_human_project_index",
    "load_mouse_project_index",
]
