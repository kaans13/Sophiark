"""Explicit, read-only adapters for biological data already shipped with Sophiark.

No provider performs I/O during import or construction.  Every lookup is local,
species-scoped, provenance-bearing, and failure-isolated; there is deliberately
no HTTP client or external fallback in this module.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from enum import Enum
import json
from pathlib import Path
import pickle
import re
import sqlite3
from typing import Any, Generic, Iterable, Mapping, TypeVar

from src.research.models import (
    ProvenanceKind,
    ProvenanceRecord,
    immutable_mapping,
)


HUMAN_TAXON_ID = 9606
MOUSE_TAXON_ID = 10090
SUPPORTED_TAXA = frozenset({HUMAN_TAXON_ID, MOUSE_TAXON_ID})


class LocalProviderStatus(str, Enum):
    AVAILABLE = "AVAILABLE"
    PARTIAL = "PARTIAL"
    NO_MATCH = "NO_MATCH"
    UNAVAILABLE = "UNAVAILABLE"
    UNSUPPORTED_TAXON = "UNSUPPORTED_TAXON"
    PROVIDER_ERROR = "PROVIDER_ERROR"


T = TypeVar("T")


@dataclass(frozen=True, slots=True)
class LocalProviderResult(Generic[T]):
    provider_id: str
    taxon_id: int
    status: LocalProviderStatus
    records: tuple[T, ...] = ()
    provenance: tuple[ProvenanceRecord, ...] = ()
    requested_count: int = 0
    matched_count: int = 0
    missing_sources: tuple[str, ...] = ()
    errors: tuple[str, ...] = ()
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "status", LocalProviderStatus(self.status))
        object.__setattr__(self, "records", tuple(self.records))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        object.__setattr__(self, "missing_sources", tuple(self.missing_sources))
        object.__setattr__(self, "errors", tuple(self.errors))

    @property
    def source_available(self) -> bool:
        return self.status in {
            LocalProviderStatus.AVAILABLE,
            LocalProviderStatus.PARTIAL,
            LocalProviderStatus.NO_MATCH,
        }


@dataclass(frozen=True, slots=True)
class GeneAnnotation:
    taxon_id: int
    canonical_id: str
    symbol: str | None
    name: str | None
    go_bp: tuple[Any, ...]
    go_cc: tuple[Any, ...]
    go_mf: tuple[Any, ...]
    provenance: ProvenanceRecord


@dataclass(frozen=True, slots=True)
class FamilyMembership:
    taxon_id: int
    query_identifier: str
    approved_symbol: str
    group_id: str
    group_name: str
    matched_via: str
    aliases: tuple[str, ...]
    ensembl_gene_id: str | None
    ncbi_gene_id: str | None
    provenance: ProvenanceRecord
    evidence_type: str = "AUTHORITATIVE_FAMILY_MEMBERSHIP"


@dataclass(frozen=True, slots=True)
class ComplexMembership:
    taxon_id: int
    member_identifier: str
    identifier_type: str
    complex_id: str
    complex_name: str
    provenance: ProvenanceRecord
    evidence_type: str = "COMPLEX_MEMBERSHIP"


@dataclass(frozen=True, slots=True)
class PartialUniProtContext:
    taxon_id: int
    ensembl_protein_id: str
    uniprot_accession: str
    structure: Mapping[str, Any]
    ligand_summary: Mapping[str, Any]
    available_fields: tuple[str, ...]
    unavailable_fields: tuple[str, ...]
    provenance: ProvenanceRecord

    def __post_init__(self) -> None:
        object.__setattr__(self, "structure", immutable_mapping(self.structure))
        object.__setattr__(self, "ligand_summary", immutable_mapping(self.ligand_summary))
        object.__setattr__(self, "available_fields", tuple(self.available_fields))
        object.__setattr__(self, "unavailable_fields", tuple(self.unavailable_fields))


@dataclass(frozen=True, slots=True)
class RegulatoryRelation:
    taxon_id: int
    source_symbol: str
    target_symbol: str
    database: str
    effect: str | None
    mechanism: str | None
    directed: bool | None
    stimulation: bool | None
    inhibition: bool | None
    reference_ids: tuple[str, ...]
    record_identifier: str | None
    qualifiers: Mapping[str, Any]
    provenance: ProvenanceRecord
    evidence_type: str = "DIRECTED_REGULATION"

    def __post_init__(self) -> None:
        object.__setattr__(self, "reference_ids", tuple(self.reference_ids))
        object.__setattr__(self, "qualifiers", immutable_mapping(self.qualifiers))


def _queries(values: Iterable[object]) -> tuple[str, ...]:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        text = str(value or "").strip()
        if text and text not in seen:
            seen.add(text)
            ordered.append(text)
    return tuple(ordered)


def _freeze_annotation(value: Any) -> Any:
    if isinstance(value, Mapping):
        return immutable_mapping(value)
    if isinstance(value, list):
        return tuple(_freeze_annotation(item) for item in value)
    return value


def _annotation_tuple(value: Any) -> tuple[Any, ...]:
    if value is None:
        return ()
    if not isinstance(value, (list, tuple)):
        value = (value,)
    return tuple(_freeze_annotation(item) for item in value)


def _split_aliases(value: object) -> tuple[str, ...]:
    text = str(value or "").strip()
    if not text or text.casefold() == "nan":
        return ()
    return tuple(
        part.strip()
        for part in re.split(r"[;,|]", text)
        if part.strip() and part.strip() != "-"
    )


def _bool(value: object) -> bool | None:
    if value is None or str(value).strip() == "":
        return None
    text = str(value).strip().casefold()
    if text in {"1", "true", "yes", "y"}:
        return True
    if text in {"0", "false", "no", "n"}:
        return False
    return None


def _source_provenance(
    *,
    source: str,
    path: Path,
    kind: ProvenanceKind,
    version: str | None = None,
    retrieved_at: str | None = None,
    details: Mapping[str, Any] | None = None,
) -> ProvenanceRecord:
    return ProvenanceRecord(
        source=source,
        kind=kind,
        version=version,
        retrieved_at=retrieved_at,
        locator=str(path),
        details=details or {},
    )


def _unsupported(provider_id: str, taxon_id: int, requested_count: int, message: str) -> LocalProviderResult[Any]:
    return LocalProviderResult(
        provider_id=provider_id,
        taxon_id=taxon_id,
        status=LocalProviderStatus.UNSUPPORTED_TAXON,
        requested_count=requested_count,
        message=message,
    )


def _missing(provider_id: str, taxon_id: int, requested_count: int, path: Path) -> LocalProviderResult[Any]:
    return LocalProviderResult(
        provider_id=provider_id,
        taxon_id=taxon_id,
        status=LocalProviderStatus.UNAVAILABLE,
        requested_count=requested_count,
        missing_sources=(str(path),),
        message="Required local source is unavailable; no external fallback was attempted.",
    )


class LocalMyGeneProvider:
    """Exact ENSP/ENSMUSP annotation lookup from existing pickle snapshots."""

    provider_id = "local_mygene"

    def __init__(
        self,
        *,
        human_cache_path: str | Path | None = None,
        mouse_cache_path: str | Path | None = None,
    ) -> None:
        self._paths = {
            HUMAN_TAXON_ID: Path(human_cache_path) if human_cache_path is not None else None,
            MOUSE_TAXON_ID: Path(mouse_cache_path) if mouse_cache_path is not None else None,
        }
        self._loaded: dict[int, Mapping[str, Any]] = {}

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalMyGeneProvider":
        root = Path(project_root)
        return cls(
            human_cache_path=root / "data" / "processed" / "mygene_cache.pkl",
            mouse_cache_path=root / "data" / "processed" / "mouse_mygene_cache.pkl",
        )

    def get_annotations(
        self, canonical_ids: Iterable[object], *, taxon_id: int
    ) -> LocalProviderResult[GeneAnnotation]:
        requested = _queries(canonical_ids)
        if taxon_id not in SUPPORTED_TAXA:
            return _unsupported(
                self.provider_id, taxon_id, len(requested), "MyGene local caches support only taxon 9606 and 10090."
            )
        path = self._paths[taxon_id]
        if path is None or not path.is_file():
            return _missing(self.provider_id, taxon_id, len(requested), path or Path("<not-configured>"))
        provenance = _source_provenance(
            source="MyGene local annotation cache",
            path=path,
            kind=ProvenanceKind.LOCAL_ANNOTATION,
            details={"taxon_id": taxon_id, "format": "pickle", "lookup": "exact canonical ID"},
        )
        try:
            cache = self._loaded.get(taxon_id)
            if cache is None:
                with path.open("rb") as handle:
                    loaded = pickle.load(handle)
                if not isinstance(loaded, Mapping):
                    raise ValueError("MyGene cache root must be a mapping")
                cache = loaded
                self._loaded[taxon_id] = cache
            records: list[GeneAnnotation] = []
            for canonical_id in requested:
                raw = cache.get(canonical_id)
                if not isinstance(raw, Mapping):
                    continue
                records.append(
                    GeneAnnotation(
                        taxon_id=taxon_id,
                        canonical_id=canonical_id,
                        symbol=str(raw.get("symbol") or "").strip() or None,
                        name=str(raw.get("name") or "").strip() or None,
                        go_bp=_annotation_tuple(raw.get("go_bp")),
                        go_cc=_annotation_tuple(raw.get("go_cc")),
                        go_mf=_annotation_tuple(raw.get("go_mf")),
                        provenance=provenance,
                    )
                )
        except Exception as exc:
            return LocalProviderResult(
                self.provider_id,
                taxon_id,
                LocalProviderStatus.PROVIDER_ERROR,
                provenance=(provenance,),
                requested_count=len(requested),
                errors=(f"{type(exc).__name__}: {exc}",),
                message="Local MyGene cache could not be read; no external fallback was attempted.",
            )
        if not records:
            status = LocalProviderStatus.NO_MATCH
        elif len(records) < len(requested):
            status = LocalProviderStatus.PARTIAL
        else:
            status = LocalProviderStatus.AVAILABLE
        return LocalProviderResult(
            self.provider_id,
            taxon_id,
            status,
            records=tuple(records),
            provenance=(provenance,),
            requested_count=len(requested),
            matched_count=len(records),
            message=(
                "Some exact canonical IDs were not present in the local cache."
                if status is LocalProviderStatus.PARTIAL
                else None
            ),
        )


class LocalHGNCFamilyProvider:
    """Authoritative human family/group lookup; never infers family from symbols."""

    provider_id = "local_hgnc_family"

    def __init__(self, path: str | Path, *, metadata_path: str | Path | None = None) -> None:
        self.path = Path(path)
        self.metadata_path = Path(metadata_path) if metadata_path is not None else None
        self._index: Mapping[str, tuple[Mapping[str, str], ...]] | None = None
        self._provenance: ProvenanceRecord | None = None

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalHGNCFamilyProvider":
        root = Path(project_root)
        return cls(
            root / "data" / "interpretation_reference" / "hgnc" / "hgnc_gene_groups.tsv",
            metadata_path=root / "data" / "interpretation_reference" / "metadata" / "hgnc_gene_groups.json",
        )

    def _load(self) -> tuple[Mapping[str, tuple[Mapping[str, str], ...]], ProvenanceRecord]:
        if self._index is not None and self._provenance is not None:
            return self._index, self._provenance
        metadata: dict[str, Any] = {}
        if self.metadata_path is not None and self.metadata_path.is_file():
            try:
                metadata = json.loads(self.metadata_path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                metadata = {}
        provenance = _source_provenance(
            source="HGNC Gene Groups",
            path=self.path,
            kind=ProvenanceKind.CURATED_DATABASE,
            version=metadata.get("version_or_release"),
            retrieved_at=metadata.get("retrieved_at"),
            details={"taxon_id": HUMAN_TAXON_ID, "authoritative_family_source": True},
        )
        mutable: dict[str, list[Mapping[str, str]]] = {}
        with self.path.open("r", encoding="utf-8-sig", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                approved = str(row.get("Approved symbol") or "").strip()
                group_id = str(row.get("Group ID") or "").strip()
                group_name = str(row.get("Group name") or "").strip()
                if not approved or not group_id or not group_name:
                    continue
                previous = _split_aliases(row.get("Previous symbols"))
                aliases = _split_aliases(row.get("Alias symbols"))
                normalized = immutable_mapping(
                    {
                        "approved": approved,
                        "group_id": group_id,
                        "group_name": group_name,
                        "previous": previous,
                        "aliases": aliases,
                        "ensembl_gene_id": str(row.get("Ensembl gene ID") or "").strip(),
                        "ncbi_gene_id": str(row.get("NCBI Gene ID") or "").strip(),
                    }
                )
                for identifier in (approved, *previous, *aliases):
                    mutable.setdefault(identifier.casefold(), []).append(normalized)
        self._index = {key: tuple(rows) for key, rows in mutable.items()}
        self._provenance = provenance
        return self._index, provenance

    def get_families(
        self, identifiers: Iterable[object], *, taxon_id: int
    ) -> LocalProviderResult[FamilyMembership]:
        requested = _queries(identifiers)
        if taxon_id != HUMAN_TAXON_ID:
            return _unsupported(
                self.provider_id,
                taxon_id,
                len(requested),
                "HGNC Gene Groups is human-only; no mouse family mapping was inferred.",
            )
        if not self.path.is_file():
            return _missing(self.provider_id, taxon_id, len(requested), self.path)
        try:
            index, provenance = self._load()
            records: list[FamilyMembership] = []
            seen: set[tuple[str, str, str]] = set()
            for query in requested:
                for row in index.get(query.casefold(), ()):
                    approved = str(row["approved"])
                    previous = tuple(row["previous"])
                    aliases = tuple(row["aliases"])
                    if query.casefold() == approved.casefold():
                        matched_via = "approved_symbol"
                    elif any(query.casefold() == value.casefold() for value in previous):
                        matched_via = "previous_symbol"
                    else:
                        matched_via = "alias_symbol"
                    identity = (query.casefold(), str(row["group_id"]), approved)
                    if identity in seen:
                        continue
                    seen.add(identity)
                    records.append(
                        FamilyMembership(
                            taxon_id=HUMAN_TAXON_ID,
                            query_identifier=query,
                            approved_symbol=approved,
                            group_id=str(row["group_id"]),
                            group_name=str(row["group_name"]),
                            matched_via=matched_via,
                            aliases=tuple(dict.fromkeys((*previous, *aliases))),
                            ensembl_gene_id=str(row["ensembl_gene_id"]) or None,
                            ncbi_gene_id=str(row["ncbi_gene_id"]) or None,
                            provenance=provenance,
                        )
                    )
        except Exception as exc:
            return LocalProviderResult(
                self.provider_id,
                taxon_id,
                LocalProviderStatus.PROVIDER_ERROR,
                requested_count=len(requested),
                errors=(f"{type(exc).__name__}: {exc}",),
                message="HGNC family source could not be read.",
            )
        status = LocalProviderStatus.AVAILABLE if records else LocalProviderStatus.NO_MATCH
        return LocalProviderResult(
            self.provider_id,
            taxon_id,
            status,
            records=tuple(records),
            provenance=(provenance,),
            requested_count=len(requested),
            matched_count=len(records),
            message=(
                None
                if records
                else "No exact authoritative HGNC family annotation matched; no symbol-prefix heuristic was used."
            ),
        )


class LocalComplexProvider:
    """Complex memberships from Complex Portal/CORUM-derived maps, never families."""

    provider_id = "local_complex_membership"

    def __init__(
        self,
        *,
        human_processed_path: str | Path | None = None,
        mouse_processed_path: str | Path | None = None,
        complex_portal_path: str | Path | None = None,
    ) -> None:
        self._processed_paths = {
            HUMAN_TAXON_ID: Path(human_processed_path) if human_processed_path is not None else None,
            MOUSE_TAXON_ID: Path(mouse_processed_path) if mouse_processed_path is not None else None,
        }
        self.complex_portal_path = Path(complex_portal_path) if complex_portal_path is not None else None

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalComplexProvider":
        root = Path(project_root)
        return cls(
            human_processed_path=root / "data" / "processed" / "complex_members.pkl",
            mouse_processed_path=root / "data" / "processed" / "mouse_complex_members.pkl",
            complex_portal_path=root / "data" / "interpretation_reference" / "complex_portal" / "complex_portal_9606.tsv",
        )

    def get_memberships(
        self, identifiers: Iterable[object], *, taxon_id: int
    ) -> LocalProviderResult[ComplexMembership]:
        requested = _queries(identifiers)
        if taxon_id not in SUPPORTED_TAXA:
            return _unsupported(
                self.provider_id, taxon_id, len(requested), "Local complex mappings support only taxon 9606 and 10090."
            )
        keys = {query.casefold(): query for query in requested}
        records: list[ComplexMembership] = []
        provenances: list[ProvenanceRecord] = []
        missing: list[str] = []
        errors: list[str] = []

        processed_path = self._processed_paths[taxon_id]
        if processed_path is None or not processed_path.is_file():
            missing.append(str(processed_path or Path("<processed-complex-map-not-configured>")))
        else:
            provenance = _source_provenance(
                source="Processed CORUM complex membership",
                path=processed_path,
                kind=ProvenanceKind.CURATED_DATABASE,
                details={
                    "taxon_id": taxon_id,
                    "evidence_type": "COMPLEX_MEMBERSHIP",
                    "not_family_annotation": True,
                    "derivation": "mouse-projected" if taxon_id == MOUSE_TAXON_ID else "human-local",
                },
            )
            provenances.append(provenance)
            try:
                with processed_path.open("rb") as handle:
                    loaded = pickle.load(handle)
                if not isinstance(loaded, Mapping):
                    raise ValueError("processed complex mapping must be a mapping")
                index = {str(key).casefold(): (str(key), value) for key, value in loaded.items()}
                for folded, original_query in keys.items():
                    hit = index.get(folded)
                    if hit is None:
                        continue
                    member, complexes = hit
                    if not isinstance(complexes, (list, tuple, set)):
                        complexes = (complexes,)
                    for name in complexes:
                        complex_name = str(name or "").strip()
                        if complex_name:
                            records.append(
                                ComplexMembership(
                                    taxon_id=taxon_id,
                                    member_identifier=member or original_query,
                                    identifier_type="gene_symbol",
                                    complex_id=f"processed:{complex_name}",
                                    complex_name=complex_name,
                                    provenance=provenance,
                                )
                            )
            except Exception as exc:
                errors.append(f"processed complex map: {type(exc).__name__}: {exc}")

        if taxon_id == HUMAN_TAXON_ID:
            portal_path = self.complex_portal_path
            if portal_path is None or not portal_path.is_file():
                missing.append(str(portal_path or Path("<complex-portal-not-configured>")))
            else:
                provenance = _source_provenance(
                    source="Complex Portal human ComplexTab",
                    path=portal_path,
                    kind=ProvenanceKind.CURATED_DATABASE,
                    details={
                        "taxon_id": HUMAN_TAXON_ID,
                        "evidence_type": "COMPLEX_MEMBERSHIP",
                        "not_family_annotation": True,
                    },
                )
                provenances.append(provenance)
                try:
                    with portal_path.open("r", encoding="utf-8-sig", newline="") as handle:
                        for row in csv.DictReader(handle, delimiter="\t"):
                            if str(row.get("Taxonomy identifier") or "").strip() != str(HUMAN_TAXON_ID):
                                continue
                            participants = str(
                                row.get("Expanded participant list")
                                or row.get("Identifiers (and stoichiometry) of molecules in complex")
                                or ""
                            )
                            for token in participants.split("|"):
                                accession = re.sub(r"\([^)]*\)$", "", token.strip())
                                if accession.casefold() not in keys:
                                    continue
                                complex_id = str(row.get("#Complex ac") or "").strip()
                                complex_name = str(row.get("Recommended name") or complex_id).strip()
                                if complex_id and complex_name:
                                    records.append(
                                        ComplexMembership(
                                            taxon_id=HUMAN_TAXON_ID,
                                            member_identifier=accession,
                                            identifier_type="uniprot",
                                            complex_id=complex_id,
                                            complex_name=complex_name,
                                            provenance=provenance,
                                        )
                                    )
                except Exception as exc:
                    errors.append(f"Complex Portal: {type(exc).__name__}: {exc}")

        unique = {
            (
                record.taxon_id,
                record.member_identifier.casefold(),
                record.complex_id,
                record.provenance.source,
            ): record
            for record in records
        }
        records = list(unique.values())
        available_sources = max(0, len(provenances) - len(errors))
        if records and (missing or errors):
            status = LocalProviderStatus.PARTIAL
        elif records:
            status = LocalProviderStatus.AVAILABLE
        elif available_sources:
            status = LocalProviderStatus.NO_MATCH
        elif errors:
            status = LocalProviderStatus.PROVIDER_ERROR
        else:
            status = LocalProviderStatus.UNAVAILABLE
        return LocalProviderResult(
            self.provider_id,
            taxon_id,
            status,
            records=tuple(records),
            provenance=tuple(provenances),
            requested_count=len(requested),
            matched_count=len(records),
            missing_sources=tuple(missing),
            errors=tuple(errors),
            message="Complex membership is not protein-family evidence.",
        )


class LocalUniProtProvider:
    """Human-only partial UniProt context from structural_binding.db, opened RO."""

    provider_id = "local_uniprot_partial"
    unavailable_fields = (
        "protein_name",
        "protein_family",
        "function",
        "subcellular_location",
        "cross_references",
    )

    def __init__(self, database_path: str | Path) -> None:
        self.database_path = Path(database_path)

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalUniProtProvider":
        return cls(Path(project_root) / "data" / "processed" / "structural_binding.db")

    def get_context(
        self, identifiers: Iterable[object], *, taxon_id: int
    ) -> LocalProviderResult[PartialUniProtContext]:
        requested = _queries(identifiers)
        if taxon_id != HUMAN_TAXON_ID:
            return _unsupported(
                self.provider_id,
                taxon_id,
                len(requested),
                "The local structural/BindingDB UniProt mapping is human-only.",
            )
        if not self.database_path.is_file():
            return _missing(self.provider_id, taxon_id, len(requested), self.database_path)
        provenance = _source_provenance(
            source="Sophiark structural/BindingDB UniProt mapping",
            path=self.database_path,
            kind=ProvenanceKind.LOCAL_ANNOTATION,
            details={
                "taxon_id": HUMAN_TAXON_ID,
                "scope": "accession, AlphaFold structural descriptors, BindingDB ligand summary only",
                "partial_uniprot": True,
                "sqlite_mode": "ro",
            },
        )
        records: list[PartialUniProtContext] = []
        connection: sqlite3.Connection | None = None
        try:
            connection = sqlite3.connect(
                self.database_path.resolve().as_uri() + "?mode=ro",
                uri=True,
                isolation_level=None,
            )
            connection.row_factory = sqlite3.Row
            connection.execute("PRAGMA query_only = ON")
            for offset in range(0, len(requested), 400):
                chunk = requested[offset : offset + 400]
                if not chunk:
                    continue
                placeholders = ",".join("?" for _ in chunk)
                rows = connection.execute(
                    f"""
                    SELECT
                        em.ensp, em.uniprot,
                        sd.length, sd.mean_plddt, sd.sasa,
                        sd.largest_pocket_volume, sd.hydrophobicity,
                        sd.radius_gyration,
                        bs.num_ligands, bs.avg_ki, bs.min_ki,
                        bs.avg_ic50, bs.min_ic50
                    FROM ensp_uniprot_map AS em
                    LEFT JOIN structural_descriptors AS sd ON sd.ensp = em.ensp
                    LEFT JOIN bindingdb_summary AS bs ON bs.uniprot = em.uniprot
                    WHERE em.ensp IN ({placeholders}) OR em.uniprot IN ({placeholders})
                    """,
                    (*chunk, *chunk),
                ).fetchall()
                for row in rows:
                    structure = {
                        key: row[key]
                        for key in (
                            "length",
                            "mean_plddt",
                            "sasa",
                            "largest_pocket_volume",
                            "hydrophobicity",
                            "radius_gyration",
                        )
                        if row[key] is not None
                    }
                    ligand = {
                        key: row[key]
                        for key in ("num_ligands", "avg_ki", "min_ki", "avg_ic50", "min_ic50")
                        if row[key] is not None
                    }
                    records.append(
                        PartialUniProtContext(
                            taxon_id=HUMAN_TAXON_ID,
                            ensembl_protein_id=str(row["ensp"]),
                            uniprot_accession=str(row["uniprot"]),
                            structure=structure,
                            ligand_summary=ligand,
                            available_fields=(
                                "uniprot_accession",
                                *(f"structure.{key}" for key in structure),
                                *(f"ligand_summary.{key}" for key in ligand),
                            ),
                            unavailable_fields=self.unavailable_fields,
                            provenance=provenance,
                        )
                    )
        except (OSError, sqlite3.Error) as exc:
            return LocalProviderResult(
                self.provider_id,
                taxon_id,
                LocalProviderStatus.PROVIDER_ERROR,
                provenance=(provenance,),
                requested_count=len(requested),
                errors=(f"{type(exc).__name__}: {exc}",),
                message="Partial local UniProt source could not be queried read-only.",
            )
        finally:
            if connection is not None:
                connection.close()
        unique = {(record.ensembl_protein_id, record.uniprot_accession): record for record in records}
        records = list(unique.values())
        status = LocalProviderStatus.PARTIAL if records else LocalProviderStatus.NO_MATCH
        return LocalProviderResult(
            self.provider_id,
            taxon_id,
            status,
            records=tuple(records),
            provenance=(provenance,),
            requested_count=len(requested),
            matched_count=len(records),
            message=(
                "Local UniProt is partial: family, function, location and cross-reference fields are unavailable."
                if records
                else "No exact ENSP or UniProt accession matched the partial local source."
            ),
        )


class LocalRegulatoryProvider:
    """Species-scoped TRRUST, raw SIGNOR and OmniPath regulatory evidence."""

    provider_id = "local_regulatory"

    def __init__(
        self,
        *,
        human_trrust_path: str | Path | None = None,
        mouse_trrust_path: str | Path | None = None,
        signor_raw_path: str | Path | None = None,
        human_omnipath_path: str | Path | None = None,
        mouse_omnipath_path: str | Path | None = None,
        human_omnipath_metadata_path: str | Path | None = None,
        mouse_omnipath_metadata_path: str | Path | None = None,
    ) -> None:
        self._trrust = {
            HUMAN_TAXON_ID: Path(human_trrust_path) if human_trrust_path is not None else None,
            MOUSE_TAXON_ID: Path(mouse_trrust_path) if mouse_trrust_path is not None else None,
        }
        self.signor_raw_path = Path(signor_raw_path) if signor_raw_path is not None else None
        self._omnipath = {
            HUMAN_TAXON_ID: Path(human_omnipath_path) if human_omnipath_path is not None else None,
            MOUSE_TAXON_ID: Path(mouse_omnipath_path) if mouse_omnipath_path is not None else None,
        }
        self._omnipath_metadata = {
            HUMAN_TAXON_ID: Path(human_omnipath_metadata_path) if human_omnipath_metadata_path is not None else None,
            MOUSE_TAXON_ID: Path(mouse_omnipath_metadata_path) if mouse_omnipath_metadata_path is not None else None,
        }

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalRegulatoryProvider":
        root = Path(project_root)
        omni = root / "data" / "regulatory" / "omnipath"
        return cls(
            human_trrust_path=root / "data" / "processed" / "target_to_regulators.pkl",
            mouse_trrust_path=root / "data" / "processed" / "mouse_target_to_regulators.pkl",
            signor_raw_path=root / "data" / "raw" / "signor.tsv",
            human_omnipath_path=omni / "omnipath_9606.tsv",
            mouse_omnipath_path=omni / "omnipath_10090.tsv",
            human_omnipath_metadata_path=omni / "omnipath_9606_metadata.json",
            mouse_omnipath_metadata_path=omni / "omnipath_10090_metadata.json",
        )

    def get_relations(
        self, symbols: Iterable[object], *, taxon_id: int
    ) -> LocalProviderResult[RegulatoryRelation]:
        requested = _queries(symbols)
        if taxon_id not in SUPPORTED_TAXA:
            return _unsupported(
                self.provider_id, taxon_id, len(requested), "Local regulatory sources support only taxon 9606 and 10090."
            )
        lookup = {symbol.casefold() for symbol in requested}
        records: list[RegulatoryRelation] = []
        provenances: list[ProvenanceRecord] = []
        missing: list[str] = []
        errors: list[str] = []

        trrust_path = self._trrust[taxon_id]
        if trrust_path is None or not trrust_path.is_file():
            missing.append(str(trrust_path or Path("<trrust-not-configured>")))
        else:
            provenance = _source_provenance(
                source="TRRUST local processed lookup",
                path=trrust_path,
                kind=ProvenanceKind.CURATED_DATABASE,
                details={
                    "taxon_id": taxon_id,
                    "direction": "source_to_target",
                    "derivation": "mouse-ortholog-projected" if taxon_id == MOUSE_TAXON_ID else "human-curated",
                },
            )
            provenances.append(provenance)
            try:
                with trrust_path.open("rb") as handle:
                    mapping = pickle.load(handle)
                if not isinstance(mapping, Mapping):
                    raise ValueError("TRRUST lookup must be a mapping")
                for target, regulators in mapping.items():
                    if not isinstance(regulators, (list, tuple, set)):
                        continue
                    for item in regulators:
                        if isinstance(item, (list, tuple)):
                            source = str(item[0] if item else "").strip()
                            effect = str(item[1] if len(item) > 1 else "").strip() or None
                        else:
                            source, effect = str(item).strip(), None
                        target_text = str(target).strip()
                        if not source or not target_text:
                            continue
                        if source.casefold() not in lookup and target_text.casefold() not in lookup:
                            continue
                        records.append(
                            RegulatoryRelation(
                                taxon_id=taxon_id,
                                source_symbol=source,
                                target_symbol=target_text,
                                database="TRRUST",
                                effect=effect,
                                mechanism=None,
                                directed=True,
                                stimulation=None,
                                inhibition=None,
                                reference_ids=(),
                                record_identifier=None,
                                qualifiers={},
                                provenance=provenance,
                            )
                        )
            except Exception as exc:
                errors.append(f"TRRUST: {type(exc).__name__}: {exc}")

        signor_path = self.signor_raw_path
        if signor_path is None or not signor_path.is_file():
            missing.append(str(signor_path or Path("<signor-not-configured>")))
        else:
            provenance = _source_provenance(
                source="SIGNOR raw local snapshot",
                path=signor_path,
                kind=ProvenanceKind.CURATED_DATABASE,
                details={"taxon_filter_required": True, "taxon_id": taxon_id},
            )
            provenances.append(provenance)
            try:
                with signor_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        row_taxa = {
                            int(value)
                            for value in re.findall(r"(?<!\d)\d+(?!\d)", str(row.get("TAX_ID") or ""))
                        }
                        if taxon_id not in row_taxa:
                            continue
                        source = str(row.get("ENTITYA") or "").strip()
                        target = str(row.get("ENTITYB") or "").strip()
                        if not source or not target:
                            continue
                        if source.casefold() not in lookup and target.casefold() not in lookup:
                            continue
                        pmids = _split_aliases(row.get("PMID"))
                        records.append(
                            RegulatoryRelation(
                                taxon_id=taxon_id,
                                source_symbol=source,
                                target_symbol=target,
                                database="SIGNOR",
                                effect=str(row.get("EFFECT") or "").strip() or None,
                                mechanism=str(row.get("MECHANISM") or "").strip() or None,
                                directed=_bool(row.get("DIRECT")),
                                stimulation=None,
                                inhibition=None,
                                reference_ids=pmids,
                                record_identifier=str(row.get("SIGNOR_ID") or "").strip() or None,
                                qualifiers={
                                    "source_uniprot": str(row.get("IDA") or "").strip() or None,
                                    "target_uniprot": str(row.get("IDB") or "").strip() or None,
                                    "score": str(row.get("SCORE") or "").strip() or None,
                                },
                                provenance=provenance,
                            )
                        )
            except Exception as exc:
                errors.append(f"SIGNOR: {type(exc).__name__}: {exc}")

        omni_path = self._omnipath[taxon_id]
        if omni_path is None or not omni_path.is_file():
            missing.append(str(omni_path or Path("<omnipath-not-configured>")))
        else:
            metadata_path = self._omnipath_metadata[taxon_id]
            metadata: dict[str, Any] = {}
            try:
                if metadata_path is not None and metadata_path.is_file():
                    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
                    declared_taxon = int(metadata.get("organism_taxid"))
                    if declared_taxon != taxon_id:
                        raise ValueError(
                            f"OmniPath metadata taxon {declared_taxon} does not match requested {taxon_id}"
                        )
                provenance = _source_provenance(
                    source="OmniPath local snapshot",
                    path=omni_path,
                    kind=ProvenanceKind.CURATED_DATABASE,
                    version=metadata.get("release"),
                    retrieved_at=metadata.get("downloaded_at"),
                    details={
                        "taxon_id": taxon_id,
                        "checksum_sha256": metadata.get("checksum_sha256"),
                        "metadata_available": bool(metadata),
                    },
                )
                provenances.append(provenance)
                with omni_path.open("r", encoding="utf-8-sig", newline="") as handle:
                    for row in csv.DictReader(handle, delimiter="\t"):
                        source = str(row.get("source_genesymbol") or row.get("source") or "").strip()
                        target = str(row.get("target_genesymbol") or row.get("target") or "").strip()
                        if not source or not target:
                            continue
                        if source.casefold() not in lookup and target.casefold() not in lookup:
                            continue
                        records.append(
                            RegulatoryRelation(
                                taxon_id=taxon_id,
                                source_symbol=source,
                                target_symbol=target,
                                database="OmniPath",
                                effect=None,
                                mechanism=None,
                                directed=_bool(row.get("consensus_direction") or row.get("is_directed")),
                                stimulation=_bool(row.get("consensus_stimulation") or row.get("is_stimulation")),
                                inhibition=_bool(row.get("consensus_inhibition") or row.get("is_inhibition")),
                                reference_ids=_split_aliases(row.get("references")),
                                record_identifier=None,
                                qualifiers={
                                    "sources": _split_aliases(row.get("sources")),
                                    "source_uniprot": str(row.get("source") or "").strip() or None,
                                    "target_uniprot": str(row.get("target") or "").strip() or None,
                                },
                                provenance=provenance,
                            )
                        )
            except Exception as exc:
                errors.append(f"OmniPath: {type(exc).__name__}: {exc}")

        unique = {
            (
                record.taxon_id,
                record.database,
                record.source_symbol.casefold(),
                record.target_symbol.casefold(),
                record.effect,
                record.mechanism,
                record.reference_ids,
                record.record_identifier,
            ): record
            for record in records
        }
        records = list(unique.values())
        # Metadata validation can fail before a provenance record is accepted.
        # Clamp at zero so an error-only lookup cannot be mistaken for an
        # available source merely because negative integers are truthy.
        successful_sources = max(0, len(provenances) - len(errors))
        if records and (missing or errors):
            status = LocalProviderStatus.PARTIAL
        elif records:
            status = LocalProviderStatus.AVAILABLE
        elif successful_sources:
            status = LocalProviderStatus.NO_MATCH
        elif errors:
            status = LocalProviderStatus.PROVIDER_ERROR
        else:
            status = LocalProviderStatus.UNAVAILABLE
        return LocalProviderResult(
            self.provider_id,
            taxon_id,
            status,
            records=tuple(records),
            provenance=tuple(provenances),
            requested_count=len(requested),
            matched_count=len(records),
            missing_sources=tuple(missing),
            errors=tuple(errors),
            message=(
                "Regulatory records are source-specific directed evidence; absence is not negative biological evidence."
            ),
        )


@dataclass(slots=True)
class LocalEvidenceAdapters:
    """Convenience facade; construction remains I/O-free."""

    mygene: LocalMyGeneProvider
    families: LocalHGNCFamilyProvider
    complexes: LocalComplexProvider
    uniprot: LocalUniProtProvider
    regulatory: LocalRegulatoryProvider

    @classmethod
    def from_project_root(cls, project_root: str | Path) -> "LocalEvidenceAdapters":
        return cls(
            mygene=LocalMyGeneProvider.from_project_root(project_root),
            families=LocalHGNCFamilyProvider.from_project_root(project_root),
            complexes=LocalComplexProvider.from_project_root(project_root),
            uniprot=LocalUniProtProvider.from_project_root(project_root),
            regulatory=LocalRegulatoryProvider.from_project_root(project_root),
        )


__all__ = [
    "ComplexMembership",
    "FamilyMembership",
    "GeneAnnotation",
    "HUMAN_TAXON_ID",
    "LocalComplexProvider",
    "LocalEvidenceAdapters",
    "LocalHGNCFamilyProvider",
    "LocalMyGeneProvider",
    "LocalProviderResult",
    "LocalProviderStatus",
    "LocalRegulatoryProvider",
    "LocalUniProtProvider",
    "MOUSE_TAXON_ID",
    "PartialUniProtContext",
    "RegulatoryRelation",
]
