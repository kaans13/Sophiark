"""Pure, read-only semantic views over a simulation result snapshot.

This adapter only classifies and labels values that Sophiark has already
computed.  It deliberately has no access to the graph, scientific engines,
session state, files, caches, or providers.  Source order is presentation
order: no row is ranked, limited, thresholded, or statistically recalculated.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
import math
from types import MappingProxyType
from typing import Any, Mapping

from .contracts import CONTRACTS, contracts_for_source_field
from .models import InterpretationContract, ObservationScope, SimulationResultSnapshot


class ResultSource(str, Enum):
    REPORT = "report"
    SIGNED_REDISTRIBUTION = "signed_redistribution"
    ENRICHMENT = "enrichment"
    SNAPSHOT_TARGETS = "snapshot_targets"


class ResultRole(str, Enum):
    PERTURBATION_TARGET = "PERTURBATION_TARGET"
    POSITIVE_REDISTRIBUTION = "POSITIVE_REDISTRIBUTION"
    NETWORK_LOSS = "NETWORK_LOSS"
    FDR_SUPPORTED = "FDR_SUPPORTED"


class ResponseDirection(str, Enum):
    GAIN = "GAIN"
    LOSS = "LOSS"
    NO_CHANGE = "NO_CHANGE"
    UNKNOWN = "UNKNOWN"
    UNAVAILABLE = "UNAVAILABLE"


class FunctionalContextKind(str, Enum):
    GO = "GO"
    KEGG = "KEGG"
    UNCLASSIFIED = "UNCLASSIFIED"


def _freeze_value(value: Any) -> Any:
    if isinstance(value, Mapping):
        return MappingProxyType({str(key): _freeze_value(item) for key, item in value.items()})
    if isinstance(value, (list, tuple)):
        return tuple(_freeze_value(item) for item in value)
    if isinstance(value, (set, frozenset)):
        return frozenset(_freeze_value(item) for item in value)
    return value


def _immutable_values(value: Mapping[str, Any] | None = None) -> Mapping[str, Any]:
    return MappingProxyType({
        str(key): _freeze_value(item)
        for key, item in (value or {}).items()
    })


@dataclass(frozen=True, slots=True)
class MetricValue:
    """One unchanged source value bound to its interpretation contract."""

    contract: InterpretationContract
    source_field: str
    raw_value: Any

    def __post_init__(self) -> None:
        object.__setattr__(self, "source_field", str(self.source_field))
        object.__setattr__(self, "raw_value", _freeze_value(self.raw_value))
        if self.source_field not in self.contract.source_fields:
            raise ValueError("metric source_field is not declared by its contract")

    @property
    def metric_id(self) -> str:
        return self.contract.metric_id


@dataclass(frozen=True, slots=True)
class SemanticRecord:
    """Immutable presentation record with explicit source and result roles."""

    source: ResultSource
    source_position: int
    canonical_id: str | None
    symbol: str | None
    roles: tuple[ResultRole, ...]
    direction: ResponseDirection
    values: Mapping[str, Any] = field(default_factory=_immutable_values)
    metrics: tuple[MetricValue, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "source", ResultSource(self.source))
        object.__setattr__(self, "roles", tuple(ResultRole(role) for role in self.roles))
        object.__setattr__(self, "direction", ResponseDirection(self.direction))
        object.__setattr__(self, "values", _immutable_values(self.values))
        object.__setattr__(self, "metrics", tuple(self.metrics))
        if self.source_position < 0:
            raise ValueError("source_position cannot be negative")
        if self.canonical_id is not None:
            object.__setattr__(self, "canonical_id", str(self.canonical_id).strip() or None)
        if self.symbol is not None:
            object.__setattr__(self, "symbol", str(self.symbol).strip() or None)

    def metric(self, metric_id: str) -> MetricValue | None:
        return next((item for item in self.metrics if item.metric_id == metric_id), None)

    @property
    def contracts(self) -> tuple[InterpretationContract, ...]:
        return tuple(item.contract for item in self.metrics)


@dataclass(frozen=True, slots=True)
class FunctionalContextRow:
    """An existing enrichment row, partitioned without new enrichment logic."""

    kind: FunctionalContextKind
    source_position: int
    contract: InterpretationContract | None
    term: str | None
    source_library: str | None
    pathway_family: str | None
    members: tuple[str, ...]
    values: Mapping[str, Any] = field(default_factory=_immutable_values)

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", FunctionalContextKind(self.kind))
        object.__setattr__(self, "members", tuple(map(str, self.members)))
        object.__setattr__(self, "values", _immutable_values(self.values))
        if self.source_position < 0:
            raise ValueError("source_position cannot be negative")
        expected = {
            FunctionalContextKind.GO: "go_enrichment",
            FunctionalContextKind.KEGG: "kegg_enrichment",
        }.get(self.kind)
        if expected is not None and (self.contract is None or self.contract.metric_id != expected):
            raise ValueError(f"{self.kind.value} rows require the {expected} contract")
        if self.kind is FunctionalContextKind.UNCLASSIFIED and self.contract is not None:
            raise ValueError("unclassified enrichment cannot be assigned a GO/KEGG contract")


@dataclass(frozen=True, slots=True)
class SemanticResultView:
    snapshot_id: str
    simulation_id: str
    scope: ObservationScope
    target_ids: tuple[str, ...]
    target_records: tuple[SemanticRecord, ...]
    positive_redistribution: tuple[SemanticRecord, ...]
    network_losses: tuple[SemanticRecord, ...]
    fdr_supported: tuple[SemanticRecord, ...]
    go_terms: tuple[FunctionalContextRow, ...]
    kegg_pathways: tuple[FunctionalContextRow, ...]
    unclassified_enrichment: tuple[FunctionalContextRow, ...]
    notices: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        for field_name in (
            "target_ids", "target_records", "positive_redistribution",
            "network_losses", "fdr_supported", "go_terms", "kegg_pathways",
            "unclassified_enrichment", "notices",
        ):
            object.__setattr__(self, field_name, tuple(getattr(self, field_name)))

    @property
    def targets(self) -> tuple[SemanticRecord, ...]:
        return self.target_records

    @property
    def functional_rows(self) -> tuple[FunctionalContextRow, ...]:
        # Reconstitute the original enrichment order across the three typed
        # views.  This is source-position restoration, not a scientific rank.
        return tuple(sorted(
            (*self.go_terms, *self.kegg_pathways, *self.unclassified_enrichment),
            key=lambda item: item.source_position,
        ))


def _is_missing(value: Any) -> bool:
    if value is None:
        return True
    if isinstance(value, float):
        return math.isnan(value)
    if isinstance(value, str):
        return value.strip().casefold() in {"", "nan", "none", "<na>"}
    return False


def _text(value: Any) -> str | None:
    return None if _is_missing(value) else str(value).strip()


def _finite_number(value: Any) -> float | None:
    if isinstance(value, bool) or _is_missing(value):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number if math.isfinite(number) else None


def _strict_true(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)) and not isinstance(value, bool):
        return value == 1
    return isinstance(value, str) and value.strip().casefold() in {"true", "1", "yes", "evet", "✓"}


def _direction(row: Mapping[str, Any]) -> ResponseDirection:
    raw = _text(row.get("Response_Direction"))
    if raw is not None:
        normalized = raw.casefold()
        if normalized == "influence gain":
            return ResponseDirection.GAIN
        if normalized == "influence loss":
            return ResponseDirection.LOSS
        if normalized == "no change":
            return ResponseDirection.NO_CHANGE
        return ResponseDirection.UNKNOWN
    delta = _finite_number(row.get("Delta_PageRank_Pct"))
    if delta is None:
        return ResponseDirection.UNAVAILABLE
    if delta > 0:
        return ResponseDirection.GAIN
    if delta < 0:
        return ResponseDirection.LOSS
    return ResponseDirection.NO_CHANGE


def _ordered_values(table, record: Mapping[str, Any]) -> Mapping[str, Any]:
    # JSON object keys are canonicalized for fingerprinting; restore the original
    # DataFrame column order for deterministic presentation/contract traversal.
    return {column: record.get(column) for column in table.columns}


def _metric_values(values: Mapping[str, Any], columns: tuple[str, ...]) -> tuple[MetricValue, ...]:
    found: list[MetricValue] = []
    seen: set[str] = set()
    for source_field in columns:
        matches = contracts_for_source_field(source_field)
        if len(matches) != 1:
            continue
        contract = matches[0]
        if contract.metric_id in seen:
            continue
        seen.add(contract.metric_id)
        found.append(MetricValue(contract, source_field, values.get(source_field)))
    return tuple(found)


def _semantic_record(
    *, source: ResultSource, position: int, table, row: Mapping[str, Any],
    roles: tuple[ResultRole, ...], direction: ResponseDirection,
) -> SemanticRecord:
    values = _ordered_values(table, row)
    return SemanticRecord(
        source=source,
        source_position=position,
        canonical_id=_text(values.get("gene")),
        symbol=_text(values.get("Symbol")),
        roles=roles,
        direction=direction,
        values=values,
        metrics=_metric_values(values, table.columns),
    )


def _is_positive_report_row(row: Mapping[str, Any], direction: ResponseDirection) -> bool:
    damage_type = _text(row.get("Hasar_Tipi"))
    if damage_type is not None and damage_type != "Üçüncül_Hasar_Stresli":
        return False
    return direction is ResponseDirection.GAIN


def _functional_kind(row: Mapping[str, Any]) -> FunctionalContextKind:
    family = _text(row.get("Pathway_Family"))
    if family is not None:
        normalized = family.casefold()
        if normalized == "kegg":
            return FunctionalContextKind.KEGG
        if normalized == "go" or normalized.startswith("go:"):
            return FunctionalContextKind.GO
        return FunctionalContextKind.UNCLASSIFIED
    # Legacy schemas used library names instead of Pathway_Family.  This
    # fallback is intentionally prefix-only and never infers from Term text.
    legacy = _text(row.get("Gene_set")) or _text(row.get("Kaynak"))
    normalized = legacy.casefold() if legacy is not None else ""
    if normalized.startswith("kegg_") or normalized == "kegg":
        return FunctionalContextKind.KEGG
    if normalized.startswith("go_") or normalized.startswith("go:") or normalized == "go":
        return FunctionalContextKind.GO
    return FunctionalContextKind.UNCLASSIFIED


def _members(value: Any) -> tuple[str, ...]:
    if _is_missing(value):
        return ()
    if isinstance(value, str):
        return tuple(item.strip() for item in value.split(";") if item.strip())
    if isinstance(value, (list, tuple, set, frozenset)):
        return tuple(str(item).strip() for item in value if str(item).strip())
    return (str(value).strip(),)


def adapt_snapshot(snapshot: SimulationResultSnapshot) -> SemanticResultView:
    """Adapt an immutable snapshot without scientific work or any I/O."""

    if not isinstance(snapshot, SimulationResultSnapshot):
        raise TypeError("snapshot must be a SimulationResultSnapshot")

    target_ids = tuple(snapshot.targets)
    target_set = set(target_ids)
    target_records: list[SemanticRecord] = []
    positives: list[SemanticRecord] = []
    fdr_records: list[SemanticRecord] = []
    notices: list[str] = []
    fdr_mode = snapshot.selection_mode.strip().casefold() == "null_fdr"

    report_rows = snapshot.report.records()
    has_fdr_field = "significant_redistribution" in snapshot.report.columns
    for position, raw_row in enumerate(report_rows):
        row = _ordered_values(snapshot.report, raw_row)
        canonical_id = _text(row.get("gene"))
        direction = _direction(row)
        is_target = canonical_id is not None and canonical_id in target_set
        is_positive = not is_target and _is_positive_report_row(row, direction)
        is_fdr = (
            fdr_mode and not is_target and has_fdr_field
            and _strict_true(row.get("significant_redistribution"))
        )
        roles: list[ResultRole] = []
        if is_target:
            roles.append(ResultRole.PERTURBATION_TARGET)
        if is_positive:
            roles.append(ResultRole.POSITIVE_REDISTRIBUTION)
        if is_fdr:
            if direction is ResponseDirection.LOSS:
                roles.append(ResultRole.NETWORK_LOSS)
            elif direction is ResponseDirection.GAIN and ResultRole.POSITIVE_REDISTRIBUTION not in roles:
                roles.append(ResultRole.POSITIVE_REDISTRIBUTION)
            roles.append(ResultRole.FDR_SUPPORTED)
        if not roles:
            continue
        record = _semantic_record(
            source=ResultSource.REPORT, position=position, table=snapshot.report,
            row=row, roles=tuple(roles), direction=direction,
        )
        if is_target:
            target_records.append(record)
        if is_positive:
            positives.append(record)
        if is_fdr:
            fdr_records.append(record)

    missing_target_ids = target_set.difference(
        record.canonical_id for record in target_records if record.canonical_id is not None
    )
    if missing_target_ids:
        notices.append("Some snapshot target IDs have no exact report row; target IDs were retained without symbol inference.")

    network_losses: list[SemanticRecord] = []
    if snapshot.signed_redistribution is None:
        notices.append("Signed redistribution source is unavailable; network-loss rows were not inferred from the report.")
    else:
        signed_rows = snapshot.signed_redistribution.records()
        if not signed_rows:
            notices.append("Signed redistribution source is present but empty.")
        for position, raw_row in enumerate(signed_rows):
            row = _ordered_values(snapshot.signed_redistribution, raw_row)
            canonical_id = _text(row.get("gene"))
            direction = _direction(row)
            if canonical_id in target_set or direction is not ResponseDirection.LOSS:
                continue
            roles = [ResultRole.NETWORK_LOSS]
            if fdr_mode and _strict_true(row.get("significant_redistribution")):
                roles.append(ResultRole.FDR_SUPPORTED)
            network_losses.append(_semantic_record(
                source=ResultSource.SIGNED_REDISTRIBUTION,
                position=position,
                table=snapshot.signed_redistribution,
                row=row,
                roles=tuple(roles),
                direction=direction,
            ))

    if not fdr_mode:
        notices.append("FDR-supported subset is not exposed outside null_fdr selection mode.")
    elif not has_fdr_field:
        notices.append("The report has no stored significant_redistribution field; q-values were not reinterpreted.")

    go_rows: list[FunctionalContextRow] = []
    kegg_rows: list[FunctionalContextRow] = []
    unclassified_rows: list[FunctionalContextRow] = []
    if snapshot.enrichment is None:
        notices.append("Functional enrichment source is unavailable.")
    else:
        enrichment_rows = snapshot.enrichment.records()
        if not enrichment_rows:
            notices.append("Functional enrichment source is present but empty.")
        for position, raw_row in enumerate(enrichment_rows):
            row = _ordered_values(snapshot.enrichment, raw_row)
            kind = _functional_kind(row)
            contract = {
                FunctionalContextKind.GO: CONTRACTS["go_enrichment"],
                FunctionalContextKind.KEGG: CONTRACTS["kegg_enrichment"],
            }.get(kind)
            functional = FunctionalContextRow(
                kind=kind,
                source_position=position,
                contract=contract,
                term=_text(row.get("Term")),
                source_library=_text(row.get("Kaynak")) or _text(row.get("Gene_set")),
                pathway_family=_text(row.get("Pathway_Family")),
                members=_members(row.get("Genes")),
                values=row,
            )
            if kind is FunctionalContextKind.GO:
                go_rows.append(functional)
            elif kind is FunctionalContextKind.KEGG:
                kegg_rows.append(functional)
            else:
                unclassified_rows.append(functional)
        if unclassified_rows:
            notices.append("Some enrichment rows lack an authoritative GO/KEGG family and remain unclassified.")

    return SemanticResultView(
        snapshot_id=snapshot.snapshot_id,
        simulation_id=snapshot.simulation_id,
        scope=snapshot.scope,
        target_ids=target_ids,
        target_records=tuple(target_records),
        positive_redistribution=tuple(positives),
        network_losses=tuple(network_losses),
        fdr_supported=tuple(fdr_records),
        go_terms=tuple(go_rows),
        kegg_pathways=tuple(kegg_rows),
        unclassified_enrichment=tuple(unclassified_rows),
        notices=tuple(notices),
    )


# Explicit aliases for callers that prefer architecture terminology.
SemanticSimulationResult = SemanticResultView
build_semantic_result = adapt_snapshot
