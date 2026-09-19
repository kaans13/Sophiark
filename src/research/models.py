"""Immutable domain models for the read-only Research Context layer.

The models in this module deliberately contain no scientific algorithms and no
Streamlit state.  They describe already-computed Sophiark results, entity
resolution outcomes, provenance, evidence relationships, and deterministic
observations without allowing those objects to feed values back into the core.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import json
from types import MappingProxyType
from typing import Any, Mapping


class EntityType(str, Enum):
    GENE = "gene"
    PROTEIN = "protein"
    PROTEIN_FAMILY = "protein_family"
    GO_TERM = "go_term"
    PATHWAY = "pathway"
    COMPARTMENT = "compartment"
    COMMUNITY = "community"
    PUBLICATION = "publication"
    SIMULATION = "simulation"


class AliasType(str, Enum):
    GENE_SYMBOL = "gene_symbol"
    ENSEMBL_GENE = "ensembl_gene"
    ENSEMBL_PROTEIN = "ensembl_protein"
    UNIPROT = "uniprot"
    STRING_ID = "string_id"
    NCBI_GENE = "ncbi_gene"
    SYNONYM = "synonym"


class ResolutionStatus(str, Enum):
    EXACT = "EXACT"
    ALIAS = "ALIAS"
    EXTERNAL_RESOLVED = "EXTERNAL_RESOLVED"
    AMBIGUOUS = "AMBIGUOUS"
    UNRESOLVED = "UNRESOLVED"


class ProvenanceKind(str, Enum):
    SOPHIARK_COMPUTED = "SOPHIARK_COMPUTED"
    LOCAL_ANNOTATION = "LOCAL_ANNOTATION"
    CURATED_DATABASE = "CURATED_DATABASE"
    ONTOLOGY = "ONTOLOGY"
    PATHWAY = "PATHWAY"
    LITERATURE = "LITERATURE"


class RelationshipType(str, Enum):
    PPI_EDGE = "PPI_EDGE"
    DIRECTED_REGULATION = "DIRECTED_REGULATION"
    SHARED_PATHWAY = "SHARED_PATHWAY"
    SHARED_GO = "SHARED_GO"
    SHARED_FAMILY = "SHARED_FAMILY"
    SHARED_COMPARTMENT = "SHARED_COMPARTMENT"
    SHARED_COMMUNITY = "SHARED_COMMUNITY"
    NETWORK_PATH = "NETWORK_PATH"
    CO_MENTIONED_IN_PUBLICATION = "CO_MENTIONED_IN_PUBLICATION"


class ObservationType(str, Enum):
    FAMILY_COOCCURRENCE = "FAMILY_COOCCURRENCE"
    FUNCTIONAL_COOCCURRENCE = "FUNCTIONAL_COOCCURRENCE"
    LOCALIZATION_PATTERN = "LOCALIZATION_PATTERN"
    POSITIVE_NEGATIVE_CONTRAST = "POSITIVE_NEGATIVE_CONTRAST"
    GATEWAY_RESPONSE_OVERLAP = "GATEWAY_RESPONSE_OVERLAP"
    FDR_SUPPORTED_SUBSET = "FDR_SUPPORTED_SUBSET"
    TARGET_RESPONSE_RELATIONSHIP = "TARGET_RESPONSE_RELATIONSHIP"


class MetricOrigin(str, Enum):
    SOPHIARK_COMPUTED = "sophiark_computed"
    LOCAL_ANNOTATION = "local_annotation"
    CURATED_DATABASE = "curated_database"
    ONTOLOGY = "ontology"
    PATHWAY = "pathway"


class MetricKind(str, Enum):
    NUMBER = "number"
    TEXT = "text"
    MARKER = "marker"
    TABLE = "table"


def _freeze_value(value: Any) -> Any:
    """Recursively make caller-owned containers immutable."""

    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze_value(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze_value(item) for item in value)
    return value


def immutable_mapping(value: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    """Return a recursively immutable, detached mapping."""

    return MappingProxyType({
        str(key): _freeze_value(item)
        for key, item in (value or {}).items()
    })


@dataclass(frozen=True, slots=True)
class Entity:
    taxon_id: int
    entity_type: EntityType
    canonical_id: str
    symbol: str | None = None
    display_name: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "entity_type", EntityType(self.entity_type))
        object.__setattr__(self, "canonical_id", str(self.canonical_id).strip())
        if self.taxon_id <= 0:
            raise ValueError("taxon_id must be a positive NCBI taxonomy identifier")
        if not self.canonical_id:
            raise ValueError("canonical_id cannot be empty")
        if self.symbol is not None:
            object.__setattr__(self, "symbol", str(self.symbol).strip() or None)
        if self.display_name is not None:
            object.__setattr__(self, "display_name", str(self.display_name).strip() or None)

    @property
    def logical_key(self) -> str:
        return f"{self.taxon_id}:{self.entity_type.value}:{self.canonical_id}"


@dataclass(frozen=True, slots=True)
class EntityAlias:
    alias_type: AliasType
    alias_value: str
    source: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "alias_type", AliasType(self.alias_type))
        object.__setattr__(self, "alias_value", str(self.alias_value).strip())
        object.__setattr__(self, "source", str(self.source).strip())
        if not self.alias_value or not self.source:
            raise ValueError("alias_value and source cannot be empty")


@dataclass(frozen=True, slots=True)
class EntityResolution:
    query: str
    taxon_id: int
    entity_type: EntityType
    status: ResolutionStatus
    entity: Entity | None = None
    matched_alias: EntityAlias | None = None
    candidates: tuple[Entity, ...] = ()
    source: str | None = None
    message: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "query", str(self.query).strip())
        object.__setattr__(self, "entity_type", EntityType(self.entity_type))
        object.__setattr__(self, "status", ResolutionStatus(self.status))
        object.__setattr__(self, "candidates", tuple(self.candidates))
        if not self.query:
            raise ValueError("resolution query cannot be empty")
        if self.taxon_id <= 0:
            raise ValueError("taxon_id must be positive")
        if self.entity is not None and self.entity.taxon_id != self.taxon_id:
            raise ValueError("resolved entity taxon does not match the requested taxon")
        if any(candidate.taxon_id != self.taxon_id for candidate in self.candidates):
            raise ValueError("resolution candidates must share the requested taxon")
        resolved_statuses = {
            ResolutionStatus.EXACT,
            ResolutionStatus.ALIAS,
            ResolutionStatus.EXTERNAL_RESOLVED,
        }
        if self.status in resolved_statuses and self.entity is None:
            raise ValueError(f"{self.status.value} resolution requires an entity")
        if self.status in {ResolutionStatus.AMBIGUOUS, ResolutionStatus.UNRESOLVED} and self.entity is not None:
            raise ValueError(f"{self.status.value} resolution cannot auto-select an entity")
        if self.status is ResolutionStatus.ALIAS and self.matched_alias is None:
            raise ValueError("ALIAS resolution requires matched_alias provenance")


@dataclass(frozen=True, slots=True)
class ProvenanceRecord:
    source: str
    kind: ProvenanceKind
    version: str | None = None
    retrieved_at: str | None = None
    locator: str | None = None
    snapshot_id: str | None = None
    attribution: str | None = None
    details: Mapping[str, Any] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "source", str(self.source).strip())
        object.__setattr__(self, "kind", ProvenanceKind(self.kind))
        object.__setattr__(self, "details", immutable_mapping(self.details))
        if not self.source:
            raise ValueError("provenance source cannot be empty")


@dataclass(frozen=True, slots=True)
class EvidenceAssertion:
    assertion_id: str
    subject: Entity
    relationship_type: RelationshipType
    object: Entity | None
    statement: str
    provenance: tuple[ProvenanceRecord, ...]
    value: Any = None
    qualifiers: Mapping[str, Any] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "relationship_type", RelationshipType(self.relationship_type))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        object.__setattr__(self, "value", _freeze_value(self.value))
        object.__setattr__(self, "qualifiers", immutable_mapping(self.qualifiers))
        if not str(self.assertion_id).strip() or not str(self.statement).strip():
            raise ValueError("assertion_id and statement cannot be empty")
        if not self.provenance:
            raise ValueError("every evidence assertion requires provenance")


@dataclass(frozen=True, slots=True)
class ObservationScope:
    selection_mode: str
    tested_count: int | None
    returned_count: int
    threshold: float | None
    top_n: int | None
    tissue: str
    species: str
    taxon_id: int
    threshold_parameters: Mapping[str, Any] = field(default_factory=immutable_mapping)

    def __post_init__(self) -> None:
        object.__setattr__(self, "threshold_parameters", immutable_mapping(self.threshold_parameters))
        if self.returned_count < 0 or (self.tested_count is not None and self.tested_count < 0):
            raise ValueError("scope counts cannot be negative")
        if self.top_n is not None and self.top_n < 0:
            raise ValueError("top_n cannot be negative")
        if self.taxon_id <= 0 or not str(self.species).strip():
            raise ValueError("scope requires explicit species and taxon_id")


@dataclass(frozen=True, slots=True)
class Observation:
    observation_id: str
    simulation_id: str
    type: ObservationType
    members: tuple[Entity, ...]
    basis: str
    scope: ObservationScope
    source_fields: tuple[str, ...]
    created_at: str
    limitations: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "type", ObservationType(self.type))
        object.__setattr__(self, "members", tuple(self.members))
        object.__setattr__(self, "source_fields", tuple(map(str, self.source_fields)))
        if not self.observation_id or not self.simulation_id or not self.basis:
            raise ValueError("observation identity and basis cannot be empty")
        if any(member.taxon_id != self.scope.taxon_id for member in self.members):
            raise ValueError("observation members must match the scope taxon")
        if not self.limitations:
            raise ValueError("observation limitations must be explicit")


@dataclass(frozen=True, slots=True)
class InterpretationContract:
    metric_id: str
    source_fields: tuple[str, ...]
    display_name: str
    origin: MetricOrigin
    unit: str | None
    precision: int | None
    signed: bool
    scope: str
    definition: str
    supports: tuple[str, ...]
    does_not_support: tuple[str, ...]
    missing_policy: str = "unavailable"
    kind: MetricKind = MetricKind.NUMBER
    suffix: str = ""

    def __post_init__(self) -> None:
        object.__setattr__(self, "origin", MetricOrigin(self.origin))
        object.__setattr__(self, "kind", MetricKind(self.kind))
        object.__setattr__(self, "source_fields", tuple(map(str, self.source_fields)))
        object.__setattr__(self, "supports", tuple(map(str, self.supports)))
        object.__setattr__(self, "does_not_support", tuple(map(str, self.does_not_support)))
        if not self.metric_id or not self.source_fields or not self.display_name:
            raise ValueError("metric_id, source_fields and display_name are required")
        if self.precision is not None and self.precision < 0:
            raise ValueError("precision cannot be negative")


@dataclass(frozen=True, slots=True)
class FrozenResultTable:
    """Immutable serialized representation of one source DataFrame."""

    name: str
    columns: tuple[str, ...]
    dtypes: tuple[str, ...]
    rows_json: str
    row_count: int
    fingerprint: str

    def __post_init__(self) -> None:
        if len(self.columns) != len(self.dtypes):
            raise ValueError("columns and dtypes must have the same length")
        if self.row_count < 0:
            raise ValueError("row_count cannot be negative")
        if not self.name or not self.fingerprint:
            raise ValueError("table name and fingerprint are required")

    @property
    def empty(self) -> bool:
        return self.row_count == 0

    def records(self) -> tuple[dict[str, Any], ...]:
        """Return newly decoded records; callers can never mutate the snapshot."""

        decoded = json.loads(self.rows_json)
        return tuple(dict(record) for record in decoded)

    def to_frame(self):
        """Materialize a detached DataFrame on demand for read-only adapters."""

        import pandas as pd

        return pd.DataFrame(self.records(), columns=list(self.columns))


@dataclass(frozen=True, slots=True)
class SimulationResultSnapshot:
    simulation_id: str
    snapshot_id: str
    created_at: str
    species: str
    taxon_id: int
    tissue: str
    targets: tuple[str, ...]
    attenuation: float | None
    selection_mode: str
    test_limit: int | None
    threshold_parameters: Mapping[str, Any]
    tested_count: int | None
    returned_count: int
    top_n: int | None
    result_schema_version: int
    report: FrozenResultTable
    signed_redistribution: FrozenResultTable | None
    enrichment: FrozenResultTable | None
    provenance: tuple[ProvenanceRecord, ...]
    scientific_fingerprint: str
    snapshot_fingerprint: str

    def __post_init__(self) -> None:
        object.__setattr__(self, "targets", tuple(map(str, self.targets)))
        object.__setattr__(self, "threshold_parameters", immutable_mapping(self.threshold_parameters))
        object.__setattr__(self, "provenance", tuple(self.provenance))
        if not self.simulation_id or not self.snapshot_id:
            raise ValueError("snapshot and simulation IDs are required")
        if self.taxon_id <= 0 or not self.species.strip():
            raise ValueError("snapshot requires explicit species and taxon_id")
        if self.returned_count < 0 or (self.tested_count is not None and self.tested_count < 0):
            raise ValueError("snapshot counts cannot be negative")
        if self.test_limit is not None and self.test_limit < 0:
            raise ValueError("test_limit cannot be negative")
        if self.top_n is not None and self.top_n < 0:
            raise ValueError("top_n cannot be negative")
        if self.result_schema_version < 1:
            raise ValueError("result_schema_version must be >= 1")

    @property
    def scope(self) -> ObservationScope:
        threshold = self.threshold_parameters.get("threshold")
        try:
            numeric_threshold = float(threshold) if threshold is not None else None
        except (TypeError, ValueError):
            numeric_threshold = None
        return ObservationScope(
            selection_mode=self.selection_mode,
            tested_count=self.tested_count,
            returned_count=self.returned_count,
            threshold=numeric_threshold,
            top_n=self.top_n,
            tissue=self.tissue,
            species=self.species,
            taxon_id=self.taxon_id,
            threshold_parameters=self.threshold_parameters,
        )

    @property
    def fingerprint(self) -> str:
        return self.snapshot_fingerprint
