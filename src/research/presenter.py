"""Deterministic, offline Presenter Mode for Research Context.

The presenter consumes an already-built offline bundle and an optional cached
external snapshot.  It performs no provider, network, file, graph, or session
work and does not generate scientific claims beyond Interpretation Contracts.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from enum import Enum
from typing import Any, Iterable, Mapping

from .contracts import CONTRACTS, format_metric
from .external_context import ExternalContextSnapshot
from .models import MetricOrigin, ObservationType, ProvenanceKind, RelationshipType
from .provider_runtime import ProviderResult, ProviderStatus
from .semantic_adapter import ResultRole


class PresenterSectionKind(str, Enum):
    PERTURBATION = "Perturbation Context"
    NETWORK_RESPONSE = "System / Network Response"
    STRUCTURAL_PATTERNS = "Detected Structural Patterns"
    FUNCTIONAL_CONTEXT = "Functional Context"
    RELATIONSHIPS = "Relationships"
    STATISTICAL_SUPPORT = "Statistical Support"
    LITERATURE = "Literature Context"
    LIMITATIONS = "Limitations"


PRESENTER_ORDER = tuple(PresenterSectionKind)


@dataclass(frozen=True, slots=True)
class PresenterFact:
    label: str
    value: str
    source_badge: str
    contract_id: str | None = None
    source_identifier: str | None = None


@dataclass(frozen=True, slots=True)
class PresenterContextItem:
    text: str
    source_badge: str
    source_identifier: str | None = None


@dataclass(frozen=True, slots=True)
class PresenterQuestion:
    question: str
    answer: str
    contract_id: str


@dataclass(frozen=True, slots=True)
class SpeakerNotes:
    say: tuple[str, ...]
    do_not_say: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "say", tuple(self.say))
        object.__setattr__(self, "do_not_say", tuple(self.do_not_say))


@dataclass(frozen=True, slots=True)
class PresenterSection:
    order: int
    kind: PresenterSectionKind
    title: str
    observation: str
    why_shown: str
    sophiark_data: tuple[PresenterFact, ...]
    biological_context: tuple[PresenterContextItem, ...]
    external_context: tuple[PresenterContextItem, ...]
    literature_context: tuple[PresenterContextItem, ...]
    interpretation_boundary: str
    possible_questions: tuple[PresenterQuestion, ...]
    speaker_notes: SpeakerNotes

    def __post_init__(self) -> None:
        object.__setattr__(self, "kind", PresenterSectionKind(self.kind))
        for name in (
            "sophiark_data", "biological_context", "external_context",
            "literature_context", "possible_questions",
        ):
            object.__setattr__(self, name, tuple(getattr(self, name)))
        if self.order != PRESENTER_ORDER.index(self.kind) + 1:
            raise ValueError("presenter section order is fixed and non-causal")


@dataclass(frozen=True, slots=True)
class PresenterDeck:
    simulation_id: str
    simulation_snapshot_id: str
    external_context_snapshot_id: str | None
    offline_ready: bool
    sections: tuple[PresenterSection, ...]
    source_badges: tuple[str, ...]
    limitations: tuple[str, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "sections", tuple(self.sections))
        object.__setattr__(self, "source_badges", tuple(self.source_badges))
        object.__setattr__(self, "limitations", tuple(self.limitations))
        if tuple(section.kind for section in self.sections) != PRESENTER_ORDER:
            raise ValueError("PresenterDeck must contain the exact eight-section order")


class ExternalDisplayState(str, Enum):
    AVAILABLE = "Available"
    CACHED = "Cached"
    NO_MATCH = "No match"
    OFFLINE = "Offline"
    TIMEOUT = "Timeout"
    RATE_LIMIT = "Rate limited"
    PROVIDER_ERROR = "Provider error"
    DISABLED = "Disabled"
    NOT_REQUESTED = "Not requested"


@dataclass(frozen=True, slots=True)
class ProviderDisplayStatus:
    provider: str
    state: ExternalDisplayState
    captured_statuses: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "state", ExternalDisplayState(self.state))
        object.__setattr__(self, "captured_statuses", tuple(self.captured_statuses))


@dataclass(frozen=True, slots=True)
class ExternalContextUIState:
    sophiark_analysis: ExternalDisplayState
    local_research_context: ExternalDisplayState
    cached_external_context: ExternalDisplayState
    live_external_context: ExternalDisplayState
    context_snapshot_id: str | None
    providers: tuple[ProviderDisplayStatus, ...]
    notices: tuple[str, ...]

    def __post_init__(self) -> None:
        for name in (
            "sophiark_analysis", "local_research_context",
            "cached_external_context", "live_external_context",
        ):
            object.__setattr__(self, name, ExternalDisplayState(getattr(self, name)))
        object.__setattr__(self, "providers", tuple(self.providers))
        object.__setattr__(self, "notices", tuple(self.notices))


def _unique(values: Iterable[Any]) -> tuple[str, ...]:
    return tuple(dict.fromkeys(
        str(value).strip() for value in values
        if value is not None and str(value).strip()
    ))


def _contract_answer(contract_id: str) -> str:
    contract = CONTRACTS[contract_id]
    supports = "; ".join(contract.supports) or contract.definition
    boundary = "; ".join(contract.does_not_support)
    return f"This field can describe: {supports}. It does not establish: {boundary}."


def _question(contract_id: str, question: str) -> PresenterQuestion:
    return PresenterQuestion(question, _contract_answer(contract_id), contract_id)


def _badge(origin: MetricOrigin) -> str:
    return {
        MetricOrigin.SOPHIARK_COMPUTED: "SOPHIARK COMPUTED",
        MetricOrigin.LOCAL_ANNOTATION: "LOCAL ANNOTATION",
        MetricOrigin.CURATED_DATABASE: "CURATED DATABASE",
        MetricOrigin.ONTOLOGY: "ONTOLOGY",
        MetricOrigin.PATHWAY: "PATHWAY",
    }[origin]


def _provenance_badge(kind: ProvenanceKind) -> str:
    return {
        ProvenanceKind.SOPHIARK_COMPUTED: "SOPHIARK COMPUTED",
        ProvenanceKind.LOCAL_ANNOTATION: "LOCAL ANNOTATION",
        ProvenanceKind.CURATED_DATABASE: "CURATED DATABASE",
        ProvenanceKind.ONTOLOGY: "ONTOLOGY",
        ProvenanceKind.PATHWAY: "PATHWAY",
        ProvenanceKind.LITERATURE: "LITERATURE",
    }[ProvenanceKind(kind)]


def _record_label(record: Any) -> str:
    return str(getattr(record, "symbol", None) or getattr(record, "canonical_id", "Unknown"))


def _value(record: Any, field: str) -> Any:
    values = getattr(record, "values", {})
    return values.get(field) if isinstance(values, Mapping) else None


def _metric_facts(records: Iterable[Any], field: str, contract_id: str, limit: int = 5) -> tuple[PresenterFact, ...]:
    contract = CONTRACTS[contract_id]
    facts = []
    for record in records:
        value = _value(record, field)
        if value is None:
            continue
        facts.append(PresenterFact(
            f"{_record_label(record)} · {contract.display_name}",
            format_metric(field, value),
            _badge(contract.origin),
            contract_id,
            getattr(record, "canonical_id", None),
        ))
        if len(facts) >= limit:
            break
    return tuple(facts)


def _external_items(snapshot: ExternalContextSnapshot | None) -> tuple[PresenterContextItem, ...]:
    if snapshot is None:
        return ()
    items = []
    for result in snapshot.provider_results():
        if result.provider == "europepmc":
            continue
        count = 0
        if isinstance(result.data, Mapping):
            records = result.data.get("records", ())
            count = len(records) if isinstance(records, (tuple, list)) else 0
        items.append(PresenterContextItem(
            f"{result.provider}: {result.status.value}; {count} cached record(s).",
            result.provider.upper(),
            result.request.canonical_query,
        ))
    return tuple(items)


def _literature_items(snapshot: ExternalContextSnapshot | None) -> tuple[PresenterContextItem, ...]:
    if snapshot is None:
        return ()
    matches = snapshot.payload.get("literature_matches", ())
    unique_publications: dict[str, Mapping[str, Any]] = {}
    passages: list[Mapping[str, Any]] = []
    for match in matches:
        if not isinstance(match, Mapping):
            continue
        publication = match.get("publication")
        if isinstance(publication, Mapping) and publication.get("publication_id"):
            unique_publications.setdefault(str(publication["publication_id"]), publication)
        for passage in match.get("relevant_passages", ()):
            if isinstance(passage, Mapping):
                passages.append(passage)
    items = [PresenterContextItem(
        f"{len(unique_publications)} matched publication(s).",
        "LITERATURE · EUROPE PMC",
    )]
    for passage in passages[:3]:
        items.append(PresenterContextItem(
            f"Relevant passage: {passage.get('text', '')}",
            "LITERATURE · EUROPE PMC",
            str(passage.get("publication_id") or "") or None,
        ))
    return tuple(items)


def _local_context_items(bundle: Any, limit: int = 8) -> tuple[PresenterContextItem, ...]:
    context = getattr(bundle, "local_context", None)
    if context is None:
        return ()
    items = []
    for entity in getattr(context, "entities", ()):
        label = entity.entity.symbol or entity.entity.canonical_id
        for annotation in entity.annotations:
            if annotation.annotation_type in {"gateway", "significant_redistribution"}:
                continue
            items.append(PresenterContextItem(
                f"{label}: {annotation.annotation_type}",
                str(annotation.source).upper(),
                entity.entity.canonical_id,
            ))
            if len(items) >= limit:
                return tuple(items)
    return tuple(items)


def _section(
    kind: PresenterSectionKind,
    *,
    observation: str,
    why: str,
    boundary: str,
    data: Iterable[PresenterFact] = (),
    biology: Iterable[PresenterContextItem] = (),
    external: Iterable[PresenterContextItem] = (),
    literature: Iterable[PresenterContextItem] = (),
    questions: Iterable[PresenterQuestion] = (),
    say: Iterable[str] = (),
    do_not_say: Iterable[str] = (),
) -> PresenterSection:
    return PresenterSection(
        order=PRESENTER_ORDER.index(kind) + 1,
        kind=kind,
        title=kind.value,
        observation=observation,
        why_shown=why,
        sophiark_data=tuple(data),
        biological_context=tuple(biology),
        external_context=tuple(external),
        literature_context=tuple(literature),
        interpretation_boundary=boundary,
        possible_questions=tuple(questions),
        speaker_notes=SpeakerNotes(tuple(say), tuple(do_not_say)),
    )


def build_presenter_deck(
    bundle: Any,
    *,
    external_snapshot: ExternalContextSnapshot | None = None,
) -> PresenterDeck:
    """Build the fixed narrative from cached/offline inputs only."""

    snapshot = bundle.snapshot
    if external_snapshot is not None and (
        external_snapshot.simulation_id != snapshot.simulation_id
        or external_snapshot.simulation_snapshot_id != snapshot.snapshot_id
    ):
        raise ValueError("external snapshot does not belong to the offline bundle")
    semantic = bundle.semantic_result
    observations = tuple(getattr(bundle.observation_result, "observations", ()))
    relationships = tuple(getattr(bundle.relationships, "assertions", ()))
    local_items = _local_context_items(bundle)
    external_items = _external_items(external_snapshot)
    literature_items = _literature_items(external_snapshot)
    observation_counts = Counter(item.type for item in observations)

    target_labels = _unique(_record_label(item) for item in semantic.target_records) or tuple(snapshot.targets)
    perturbation_facts = (
        PresenterFact("Species", snapshot.species, "SOPHIARK COMPUTED"),
        PresenterFact("Tissue", snapshot.tissue, "SOPHIARK COMPUTED"),
        PresenterFact("Target(s)", " · ".join(target_labels) or "—", "SOPHIARK COMPUTED"),
        PresenterFact("Attenuation", str(snapshot.attenuation) if snapshot.attenuation is not None else "—", "SOPHIARK COMPUTED"),
    )
    network_facts = (
        PresenterFact("Returned positive redistribution", str(len(semantic.positive_redistribution)), "SOPHIARK COMPUTED"),
        PresenterFact("Returned network losses", str(len(semantic.network_losses)), "SOPHIARK COMPUTED"),
        *_metric_facts(semantic.positive_redistribution, "Delta_PageRank_Pct", "delta_pagerank_pct"),
        *_metric_facts(semantic.network_losses, "Delta_PageRank_Pct", "delta_pagerank_pct"),
    )
    structural_facts = tuple(
        PresenterFact(
            observation.type.value.replace("_", " ").title(),
            str(len(observation.members)),
            "SOPHIARK COMPUTED + LOCAL ANNOTATION",
        )
        for observation in observations
        if observation.type in {
            ObservationType.FAMILY_COOCCURRENCE,
            ObservationType.GATEWAY_RESPONSE_OVERLAP,
            ObservationType.LOCALIZATION_PATTERN,
            ObservationType.POSITIVE_NEGATIVE_CONTRAST,
        }
    )[:8]
    functional_rows = tuple(getattr(bundle.local_context, "functional_context", ()))
    functional_facts = tuple(
        PresenterFact(
            str(getattr(row, "term", None) or getattr(row, "source_library", "Functional context")),
            "Recorded among the current returned set",
            "ONTOLOGY" if "GO" in str(getattr(row, "pathway_family", "")).upper() else "PATHWAY",
        )
        for row in functional_rows[:8]
    )
    relationship_counts = Counter(assertion.relationship_type for assertion in relationships)
    relationship_facts = tuple(
        PresenterFact(
            kind.value,
            str(count),
            (
                "LITERATURE"
                if kind is RelationshipType.CO_MENTIONED_IN_PUBLICATION
                else "SOPHIARK COMPUTED"
                if kind in {
                    RelationshipType.PPI_EDGE,
                    RelationshipType.NETWORK_PATH,
                    RelationshipType.SHARED_COMMUNITY,
                }
                else "CURATED DATABASE"
            ),
        )
        for kind, count in sorted(relationship_counts.items(), key=lambda item: item[0].value)
    )
    fdr_facts = _metric_facts(semantic.fdr_supported, "q_value", "q_value")
    if not fdr_facts:
        fdr_facts = (
            PresenterFact("FDR-supported returned records", str(len(semantic.fdr_supported)), "SOPHIARK COMPUTED", "q_value"),
        )

    sections = (
        _section(
            PresenterSectionKind.PERTURBATION,
            observation="The presentation is bound to one immutable Sophiark simulation snapshot.",
            why="This establishes the species, tissue, target and attenuation scope before any context is shown.",
            boundary="Attenuation is a model operation; it is not a measured expression, abundance or clinical intervention.",
            data=perturbation_facts,
            say=("State the exact simulation scope and target.",),
            do_not_say=("The target was experimentally knocked down.",),
        ),
        _section(
            PresenterSectionKind.NETWORK_RESPONSE,
            observation="The current result contains returned network-importance gains and losses.",
            why="These records summarize PageRank redistribution already computed by Sophiark.",
            boundary="PageRank redistribution does not establish expression change, activation, inhibition or causality.",
            data=network_facts,
            questions=(_question("delta_pagerank_pct", "Does this mean gene expression changed?"),),
            say=("Describe relative network-importance changes in this model.",),
            do_not_say=("These genes were activated or inhibited.",),
        ),
        _section(
            PresenterSectionKind.STRUCTURAL_PATTERNS,
            observation=f"{len(structural_facts)} structural/local annotation pattern(s) are available for presentation.",
            why="Only deterministic observations detected among the current returned records are shown.",
            boundary="Co-occurrence does not establish enrichment, coordinated activation or a shared causal mechanism.",
            data=structural_facts,
            biology=local_items,
            say=("Describe which annotated members co-occur in the returned set.",),
            do_not_say=("The family or compartment is activated.", "The pattern is experimentally validated."),
        ),
        _section(
            PresenterSectionKind.FUNCTIONAL_CONTEXT,
            observation=f"{len(functional_rows)} existing functional-context row(s) are attached to the current result.",
            why="Recorded GO/KEGG/Reactome context can help organize the returned entities.",
            boundary="Functional annotation or enrichment does not mean pathway activation, inhibition or causality.",
            data=functional_facts,
            biology=local_items,
            external=external_items,
            questions=(_question("go_enrichment", "Does enrichment mean pathway activation?"),) if functional_rows else (),
            say=("Name the recorded functional term and its source.",),
            do_not_say=("The pathway is activated.",),
        ),
        _section(
            PresenterSectionKind.RELATIONSHIPS,
            observation=f"{len(relationships)} typed relationship record(s) cover selected current-result pairs.",
            why="Typed records keep Sophiark network, curated annotation and literature co-mention separate.",
            boundary="A shared annotation or publication co-mention is not a physical interaction or causal edge.",
            data=relationship_facts,
            external=external_items,
            literature=literature_items,
            say=("Always name the visible relationship type.",),
            do_not_say=("Co-mentioned proteins interact.", "External context adds an edge to the Sophiark graph."),
        ),
        _section(
            PresenterSectionKind.STATISTICAL_SUPPORT,
            observation=f"{len(semantic.fdr_supported)} returned record(s) carry the existing FDR-supported role.",
            why="Stored q-values and flags disclose support under the configured null model.",
            boundary="FDR support is not effect size, biological importance, causality or experimental validation.",
            data=fdr_facts,
            questions=(_question("q_value", "Is this experimentally validated?"),),
            say=("State the null-model and multiple-testing boundary.",),
            do_not_say=("The result is biologically proven.",),
        ),
        _section(
            PresenterSectionKind.LITERATURE,
            observation=(literature_items[0].text if literature_items else "No cached literature context is attached."),
            why="Matched publications and source passages are shown only for the stored bounded queries.",
            boundary="A query match is a co-mention, not validation, proof, contradiction or causal evidence.",
            literature=literature_items,
            say=("Use the terms matched publications, literature context and relevant passages.",),
            do_not_say=("The literature validates the mechanism.", "No result means a novel discovery."),
        ),
        _section(
            PresenterSectionKind.LIMITATIONS,
            observation="Interpretation boundaries remain part of the presentation.",
            why="The research layer must remain visibly separate from Sophiark-computed results.",
            boundary="No presenter section changes the scientific result, adds causal claims or substitutes for experimental review.",
            external=external_items,
            literature=literature_items,
            say=("Close with provenance, unavailable context and explicit limitations.",),
            do_not_say=("External context improved the simulation score.",),
        ),
    )
    badges = _unique((
        "SOPHIARK COMPUTED",
        *(
            _provenance_badge(record.kind)
            for record in getattr(bundle.local_context, "provenance", ())
        ),
        *(item.source_badge for item in local_items),
        *(item.source_badge for item in external_items),
        *(item.source_badge for item in literature_items),
    ))
    return PresenterDeck(
        simulation_id=snapshot.simulation_id,
        simulation_snapshot_id=snapshot.snapshot_id,
        external_context_snapshot_id=(
            external_snapshot.context_snapshot_id if external_snapshot is not None else None
        ),
        offline_ready=True,
        sections=sections,
        source_badges=badges,
        limitations=tuple(section.interpretation_boundary for section in sections),
    )


_LIVE_STATE = {
    ProviderStatus.AVAILABLE: ExternalDisplayState.AVAILABLE,
    ProviderStatus.NO_MATCH: ExternalDisplayState.NO_MATCH,
    ProviderStatus.OFFLINE: ExternalDisplayState.OFFLINE,
    ProviderStatus.TIMEOUT: ExternalDisplayState.TIMEOUT,
    ProviderStatus.RATE_LIMIT: ExternalDisplayState.RATE_LIMIT,
    ProviderStatus.PROVIDER_ERROR: ExternalDisplayState.PROVIDER_ERROR,
    ProviderStatus.DISABLED: ExternalDisplayState.DISABLED,
}


def build_external_context_ui_state(
    bundle: Any,
    *,
    cached_snapshot: ExternalContextSnapshot | None = None,
    live_results: Iterable[ProviderResult] = (),
    live_requested: bool = False,
) -> ExternalContextUIState:
    if cached_snapshot is not None and (
        cached_snapshot.simulation_id != bundle.snapshot.simulation_id
        or cached_snapshot.simulation_snapshot_id != bundle.snapshot.snapshot_id
    ):
        raise ValueError("cached snapshot does not match the current result")
    live = tuple(live_results)
    provider_states: dict[str, ProviderDisplayStatus] = {}
    if cached_snapshot is not None:
        for provider in cached_snapshot.providers:
            provider_states[provider] = ProviderDisplayStatus(
                provider,
                ExternalDisplayState.CACHED,
                cached_snapshot.provider_statuses.get(provider, ()),
            )
    for result in live:
        previous = provider_states.get(result.provider)
        provider_states[result.provider] = ProviderDisplayStatus(
            result.provider,
            _LIVE_STATE[result.status],
            tuple((*previous.captured_statuses, result.status.value)) if previous else (result.status.value,),
        )
    if live:
        live_state = (
            ExternalDisplayState.AVAILABLE
            if any(result.status in {ProviderStatus.AVAILABLE, ProviderStatus.NO_MATCH} for result in live)
            else _LIVE_STATE[live[0].status]
        )
    else:
        live_state = ExternalDisplayState.OFFLINE if live_requested else ExternalDisplayState.NOT_REQUESTED
    local_errors = tuple(getattr(bundle.local_context, "errors", ()))
    notices = []
    if cached_snapshot is None:
        notices.append("No cached external context is attached to the current result.")
    if live_requested and not live:
        notices.append("Live external context is offline; Sophiark and local context remain available.")
    return ExternalContextUIState(
        sophiark_analysis=ExternalDisplayState.AVAILABLE,
        local_research_context=(
            ExternalDisplayState.AVAILABLE if not local_errors else ExternalDisplayState.PROVIDER_ERROR
        ),
        cached_external_context=(
            ExternalDisplayState.CACHED if cached_snapshot is not None else ExternalDisplayState.NOT_REQUESTED
        ),
        live_external_context=live_state,
        context_snapshot_id=(cached_snapshot.context_snapshot_id if cached_snapshot else None),
        providers=tuple(provider_states.values()),
        notices=tuple(notices),
    )


def render_presenter_mode(deck: PresenterDeck, *, ui: Any) -> None:
    """Render the prebuilt deck through an injected Streamlit-like object."""

    ui.subheader("Presenter Mode")
    ui.caption(f"Simulation: {deck.simulation_id}")
    for section in deck.sections:
        with ui.expander(f"{section.order}. {section.title}", expanded=section.order == 1):
            ui.markdown(f"**Observation**  \n{section.observation}")
            ui.markdown(f"**Why it is shown**  \n{section.why_shown}")
            for fact in section.sophiark_data:
                ui.write(f"{fact.label}: {fact.value} [{fact.source_badge}]")
            for item in section.biological_context:
                ui.write(f"{item.text} [{item.source_badge}]")
            for item in section.external_context:
                ui.write(f"{item.text} [{item.source_badge}]")
            for item in section.literature_context:
                ui.write(f"{item.text} [{item.source_badge}]")
            ui.warning(section.interpretation_boundary)
            for question in section.possible_questions:
                ui.markdown(f"**Possible question:** {question.question}  \n{question.answer}")
            if section.speaker_notes.say:
                ui.caption("SAY: " + " ".join(section.speaker_notes.say))
            if section.speaker_notes.do_not_say:
                ui.caption("DO NOT SAY: " + " ".join(section.speaker_notes.do_not_say))


def render_external_context_status(state: ExternalContextUIState, *, ui: Any) -> None:
    """Render safe high-level status without provider exception dumps."""

    ui.write(f"Sophiark analysis: {state.sophiark_analysis.value}")
    ui.write(f"Local Research Context: {state.local_research_context.value}")
    ui.write(f"Cached External Context: {state.cached_external_context.value}")
    ui.write(f"Live External Context: {state.live_external_context.value}")
    for provider in state.providers:
        ui.caption(f"{provider.provider}: {provider.state.value}")
    for notice in state.notices:
        ui.info(notice)


__all__ = [
    "ExternalContextUIState",
    "ExternalDisplayState",
    "PRESENTER_ORDER",
    "PresenterContextItem",
    "PresenterDeck",
    "PresenterFact",
    "PresenterQuestion",
    "PresenterSection",
    "PresenterSectionKind",
    "ProviderDisplayStatus",
    "SpeakerNotes",
    "build_external_context_ui_state",
    "build_presenter_deck",
    "render_external_context_status",
    "render_presenter_mode",
]
