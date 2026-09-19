"""Future interpreter boundary backed by a deterministic offline default.

Only bounded, structured evidence may cross this boundary. Raw scientific
tables, model clients, credentials, provider execution and network I/O are
deliberately absent.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass
from enum import Enum
import hashlib
import json
from typing import Any, Iterable

from .contracts import CONTRACTS
from .external_context import ExternalContextSnapshot
from .presenter import PresenterDeck, build_presenter_deck


class InterpreterSectionKind(str, Enum):
    SIMULATION_CONTEXT = "Simulation Context"
    SELECTED_OBSERVATIONS = "Selected Observations"
    INTERPRETATION_CONTRACTS = "Interpretation Contracts"
    EVIDENCE_GRAPH_SUBSET = "Evidence Graph subset"
    EXTERNAL_ASSERTIONS = "External Assertions"
    RELEVANT_PUBLICATIONS = "Relevant Publications"
    PROVENANCE = "Provenance"
    LIMITATIONS = "Limitations"


INTERPRETER_SECTION_ORDER = tuple(InterpreterSectionKind)


@dataclass(frozen=True, slots=True)
class InterpreterContextItem:
    identifier: str
    text: str
    source_badge: str
    source_identifier: str | None = None

    def __post_init__(self) -> None:
        for name in ("identifier", "text", "source_badge"):
            value = str(getattr(self, name)).strip()
            if not value:
                raise ValueError(f"{name} cannot be empty")
            object.__setattr__(self, name, value)
        if self.source_identifier is not None:
            object.__setattr__(
                self,
                "source_identifier",
                str(self.source_identifier).strip() or None,
            )
        forbidden = ("rows_json", "raw_report", "raw report")
        if any(token in self.identifier.casefold() for token in forbidden):
            raise ValueError("raw scientific tables are not valid interpreter context")


@dataclass(frozen=True, slots=True)
class InterpreterContext:
    context_id: str
    simulation_id: str
    simulation_snapshot_id: str
    external_context_snapshot_id: str | None
    simulation_context: tuple[InterpreterContextItem, ...]
    selected_observations: tuple[InterpreterContextItem, ...]
    interpretation_contracts: tuple[InterpreterContextItem, ...]
    evidence_graph_subset: tuple[InterpreterContextItem, ...]
    external_assertions: tuple[InterpreterContextItem, ...]
    relevant_publications: tuple[InterpreterContextItem, ...]
    provenance: tuple[InterpreterContextItem, ...]
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        item_fields = (
            "simulation_context",
            "selected_observations",
            "interpretation_contracts",
            "evidence_graph_subset",
            "external_assertions",
            "relevant_publications",
            "provenance",
        )
        for name in item_fields:
            values = tuple(getattr(self, name))
            if not all(isinstance(item, InterpreterContextItem) for item in values):
                raise TypeError(f"{name} must contain InterpreterContextItem values")
            object.__setattr__(self, name, values)
        object.__setattr__(self, "limitations", _unique(self.limitations))
        if not self.context_id.startswith("interpreter-context-"):
            raise ValueError("context_id must be a deterministic interpreter context ID")
        if not self.simulation_id or not self.simulation_snapshot_id:
            raise ValueError("interpreter context must be snapshot-bound")
        expected = _context_id(
            self.simulation_id,
            self.simulation_snapshot_id,
            self.external_context_snapshot_id,
            self.section_items(),
        )
        if self.context_id != expected:
            raise ValueError("interpreter context integrity check failed")

    def section_items(
        self,
    ) -> tuple[tuple[InterpreterSectionKind, tuple[InterpreterContextItem, ...]], ...]:
        return _sections(
            self.simulation_context,
            self.selected_observations,
            self.interpretation_contracts,
            self.evidence_graph_subset,
            self.external_assertions,
            self.relevant_publications,
            self.provenance,
            self.limitations,
        )


@dataclass(frozen=True, slots=True)
class InterpreterSection:
    kind: InterpreterSectionKind
    items: tuple[InterpreterContextItem, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", InterpreterSectionKind(self.kind))
        values = tuple(self.items)
        if not all(isinstance(item, InterpreterContextItem) for item in values):
            raise TypeError("items must contain InterpreterContextItem values")
        object.__setattr__(self, "items", values)


@dataclass(frozen=True, slots=True)
class InterpreterRender:
    interpreter: str
    context_id: str
    sections: tuple[InterpreterSection, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "sections", tuple(self.sections))
        if tuple(section.kind for section in self.sections) != INTERPRETER_SECTION_ORDER:
            raise ValueError("interpreter render must retain the exact context section order")
        if not self.interpreter or not self.context_id:
            raise ValueError("interpreter render identity cannot be empty")


class ResearchInterpreter(ABC):
    """Stable extension point; concrete model adapters are intentionally absent."""

    @abstractmethod
    def render(self, context: InterpreterContext) -> InterpreterRender:
        """Render only the bounded, provenance-preserving context."""


class DeterministicInterpreter(ResearchInterpreter):
    """Offline default that preserves context without generating new claims."""

    def render(self, context: InterpreterContext) -> InterpreterRender:
        if not isinstance(context, InterpreterContext):
            raise TypeError("context must be an InterpreterContext")
        return InterpreterRender(
            interpreter="deterministic",
            context_id=context.context_id,
            sections=tuple(
                InterpreterSection(kind, items)
                for kind, items in context.section_items()
            ),
        )


def _unique(values: Iterable[Any]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(
        str(value).strip() for value in values
        if value is not None and str(value).strip()
    ))


def _item(
    identifier: str,
    text: Any,
    source_badge: str,
    source_identifier: str | None = None,
) -> InterpreterContextItem:
    return InterpreterContextItem(
        identifier=identifier,
        text=str(text),
        source_badge=source_badge,
        source_identifier=source_identifier,
    )


def _sections(
    simulation_context: tuple[InterpreterContextItem, ...],
    selected_observations: tuple[InterpreterContextItem, ...],
    interpretation_contracts: tuple[InterpreterContextItem, ...],
    evidence_graph_subset: tuple[InterpreterContextItem, ...],
    external_assertions: tuple[InterpreterContextItem, ...],
    relevant_publications: tuple[InterpreterContextItem, ...],
    provenance: tuple[InterpreterContextItem, ...],
    limitations: tuple[str, ...],
) -> tuple[tuple[InterpreterSectionKind, tuple[InterpreterContextItem, ...]], ...]:
    limitation_items = tuple(
        _item(f"limitation:{index}", text, "INTERPRETATION BOUNDARY")
        for index, text in enumerate(limitations, start=1)
    )
    return (
        (InterpreterSectionKind.SIMULATION_CONTEXT, simulation_context),
        (InterpreterSectionKind.SELECTED_OBSERVATIONS, selected_observations),
        (InterpreterSectionKind.INTERPRETATION_CONTRACTS, interpretation_contracts),
        (InterpreterSectionKind.EVIDENCE_GRAPH_SUBSET, evidence_graph_subset),
        (InterpreterSectionKind.EXTERNAL_ASSERTIONS, external_assertions),
        (InterpreterSectionKind.RELEVANT_PUBLICATIONS, relevant_publications),
        (InterpreterSectionKind.PROVENANCE, provenance),
        (InterpreterSectionKind.LIMITATIONS, limitation_items),
    )


def _stable_item(item: InterpreterContextItem) -> dict[str, Any]:
    return {
        "identifier": item.identifier,
        "text": item.text,
        "source_badge": item.source_badge,
        "source_identifier": item.source_identifier,
    }


def _context_id(
    simulation_id: str,
    simulation_snapshot_id: str,
    external_context_snapshot_id: str | None,
    sections: tuple[tuple[InterpreterSectionKind, tuple[InterpreterContextItem, ...]], ...],
) -> str:
    payload = {
        "simulation_id": simulation_id,
        "simulation_snapshot_id": simulation_snapshot_id,
        "external_context_snapshot_id": external_context_snapshot_id,
        "sections": [
            {
                "kind": kind.value,
                "items": [_stable_item(item) for item in items],
            }
            for kind, items in sections
        ],
    }
    encoded = json.dumps(
        payload,
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "interpreter-context-" + hashlib.sha256(encoded).hexdigest()[:32]


def _presenter_items(deck: PresenterDeck, attribute: str) -> tuple[InterpreterContextItem, ...]:
    values: list[InterpreterContextItem] = []
    for section in deck.sections:
        for index, value in enumerate(getattr(section, attribute), start=1):
            values.append(_item(
                f"{section.kind.name.casefold()}:{attribute}:{index}",
                value.text,
                value.source_badge,
                value.source_identifier,
            ))
    return tuple(values)


def build_interpreter_context(
    bundle: Any,
    *,
    presenter_deck: PresenterDeck | None = None,
    external_snapshot: ExternalContextSnapshot | None = None,
) -> InterpreterContext:
    """Project completed evidence into the bounded future-adapter contract."""

    deck = presenter_deck or build_presenter_deck(
        bundle,
        external_snapshot=external_snapshot,
    )
    if deck.simulation_id != bundle.snapshot.simulation_id:
        raise ValueError("presenter deck does not belong to the current simulation")
    if deck.simulation_snapshot_id != bundle.snapshot.snapshot_id:
        raise ValueError("presenter deck is stale")
    if external_snapshot is not None:
        if external_snapshot.simulation_id != bundle.snapshot.simulation_id:
            raise ValueError("external context does not belong to the current simulation")
        if external_snapshot.simulation_snapshot_id != bundle.snapshot.snapshot_id:
            raise ValueError("external context is stale")
        if deck.external_context_snapshot_id != external_snapshot.context_snapshot_id:
            raise ValueError("presenter deck and external context snapshot do not match")

    snapshot = bundle.snapshot
    simulation_context = (
        _item("simulation_id", snapshot.simulation_id, "SOPHIARK COMPUTED"),
        _item("simulation_snapshot_id", snapshot.snapshot_id, "SOPHIARK COMPUTED"),
        _item("species", snapshot.species, "SOPHIARK COMPUTED"),
        _item("taxon_id", snapshot.taxon_id, "SOPHIARK COMPUTED"),
        _item("tissue", snapshot.tissue or "Not recorded", "SOPHIARK COMPUTED"),
        _item("targets", ", ".join(snapshot.targets) or "None recorded", "SOPHIARK COMPUTED"),
        _item("selection_mode", snapshot.selection_mode, "SOPHIARK COMPUTED"),
        _item("returned_count", snapshot.returned_count, "SOPHIARK COMPUTED"),
    )
    selected_observations = tuple(
        _item(
            detected.observation_id,
            (
                f"{detected.type.value}: {detected.basis} "
                f"Members: {', '.join(member.logical_key for member in detected.members)}"
            ),
            "SOPHIARK COMPUTED",
            detected.observation_id,
        )
        for detected in bundle.observation_result.observations[:12]
    )
    contract_ids = tuple(dict.fromkeys(
        fact.contract_id
        for section in deck.sections
        for fact in section.sophiark_data
        if fact.contract_id in CONTRACTS
    ))
    interpretation_contracts = tuple(
        _item(
            contract_id,
            (
                f"Definition: {CONTRACTS[contract_id].definition} "
                f"Supports: {'; '.join(CONTRACTS[contract_id].supports)}. "
                f"Does not support: {'; '.join(CONTRACTS[contract_id].does_not_support)}."
            ),
            "INTERPRETATION CONTRACT",
            contract_id,
        )
        for contract_id in contract_ids[:16]
    )
    evidence_graph_subset = tuple(
        _item(
            assertion.assertion_id,
            assertion.statement,
            (
                assertion.provenance[0].kind.value.replace("_", " ")
                if assertion.provenance
                else "SOPHIARK COMPUTED"
            ),
            assertion.assertion_id,
        )
        for assertion in bundle.relationships.assertions[:12]
    )
    external_assertions = _presenter_items(deck, "external_context")[:12]
    relevant_publications = _presenter_items(deck, "literature_context")[:12]
    provenance_items = [
        _item(
            f"simulation-provenance:{index}",
            record.source,
            record.kind.value.replace("_", " "),
            record.locator or record.snapshot_id,
        )
        for index, record in enumerate(snapshot.provenance, start=1)
    ]
    provenance_items.extend(
        _item(f"source-badge:{index}", badge, badge)
        for index, badge in enumerate(deck.source_badges, start=1)
    )
    if external_snapshot is not None:
        provenance_items.append(_item(
            "external_context_snapshot_id",
            external_snapshot.context_snapshot_id,
            "CACHED EXTERNAL CONTEXT",
            external_snapshot.context_snapshot_id,
        ))
    limitations = _unique((
        *(detected.limitations for detected in bundle.observation_result.observations[:12]),
        *(section.interpretation_boundary for section in deck.sections),
        *deck.limitations,
        *bundle.notices,
        *bundle.errors,
        "This context excludes raw scientific result tables and cannot alter scientific results.",
        "External assertions and publication co-mentions do not establish causality or evidence strength.",
    ))
    sections = _sections(
        simulation_context,
        selected_observations,
        interpretation_contracts,
        evidence_graph_subset,
        external_assertions,
        relevant_publications,
        tuple(provenance_items),
        limitations,
    )
    return InterpreterContext(
        context_id=_context_id(
            snapshot.simulation_id,
            snapshot.snapshot_id,
            deck.external_context_snapshot_id,
            sections,
        ),
        simulation_id=snapshot.simulation_id,
        simulation_snapshot_id=snapshot.snapshot_id,
        external_context_snapshot_id=deck.external_context_snapshot_id,
        simulation_context=simulation_context,
        selected_observations=selected_observations,
        interpretation_contracts=interpretation_contracts,
        evidence_graph_subset=evidence_graph_subset,
        external_assertions=external_assertions,
        relevant_publications=relevant_publications,
        provenance=tuple(provenance_items),
        limitations=limitations,
    )


__all__ = [
    "DeterministicInterpreter",
    "INTERPRETER_SECTION_ORDER",
    "InterpreterContext",
    "InterpreterContextItem",
    "InterpreterRender",
    "InterpreterSection",
    "InterpreterSectionKind",
    "ResearchInterpreter",
    "build_interpreter_context",
]
