"""Phase 10 acceptance tests for deterministic offline Presenter Mode."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from types import SimpleNamespace
from unittest.mock import patch

import pytest

from src.research.presenter import (
    PRESENTER_ORDER,
    ExternalDisplayState,
    PresenterDeck,
    PresenterSectionKind,
    build_external_context_ui_state,
    build_presenter_deck,
    render_external_context_status,
    render_presenter_mode,
)
from src.research.provider_runtime import ProviderRequest, ProviderResult, ProviderStatus
from tests.test_research_external_context_snapshot import _snapshot


class _Expander:
    def __init__(self, ui, label):
        self.ui = ui
        self.label = label

    def __enter__(self):
        self.ui.calls.append(("enter", (self.label,), {}))
        return self

    def __exit__(self, exc_type, exc, tb):
        self.ui.calls.append(("exit", (self.label,), {}))


class FakeUI:
    def __init__(self):
        self.calls = []

    def expander(self, label, **kwargs):
        self.calls.append(("expander", (label,), kwargs))
        return _Expander(self, label)

    def __getattr__(self, name):
        def method(*args, **kwargs):
            self.calls.append((name, args, kwargs))
        return method

    @property
    def text(self):
        return " ".join(
            str(value)
            for _, args, kwargs in self.calls
            for value in (*args, *kwargs.values())
        )


def test_presenter_has_exact_non_causal_order_and_complete_section_contract() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle, external_snapshot=external)
    assert isinstance(deck, PresenterDeck)
    assert deck.offline_ready is True
    assert tuple(section.kind for section in deck.sections) == PRESENTER_ORDER
    assert len(deck.sections) == 8
    for index, section in enumerate(deck.sections, 1):
        assert section.order == index
        assert section.title
        assert section.observation
        assert section.why_shown
        assert section.interpretation_boundary
        assert section.speaker_notes.say
        assert section.speaker_notes.do_not_say
    assert deck.external_context_snapshot_id == external.context_snapshot_id
    assert "SOPHIARK COMPUTED" in deck.source_badges
    assert "UNIPROT" in deck.source_badges
    assert "LITERATURE · EUROPE PMC" in deck.source_badges
    with pytest.raises(FrozenInstanceError):
        deck.offline_ready = False


def test_presenter_uses_contract_answers_and_separates_computed_external_literature() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle, external_snapshot=external)
    network = deck.sections[1]
    functional = deck.sections[3]
    statistical = deck.sections[5]
    literature = deck.sections[6]
    assert network.possible_questions[0].contract_id == "delta_pagerank_pct"
    assert "It does not establish:" in network.possible_questions[0].answer
    assert statistical.possible_questions[0].contract_id == "q_value"
    assert all(fact.source_badge == "SOPHIARK COMPUTED" for fact in network.sophiark_data)
    assert functional.external_context
    assert all(item.source_badge in {"UNIPROT"} for item in functional.external_context)
    assert literature.literature_context[0].text == "1 matched publication(s)."
    assert any("Relevant passage:" in item.text for item in literature.literature_context)
    assert "co-mention" in literature.interpretation_boundary.casefold()
    assert not hasattr(deck, "score")


def test_presenter_build_is_provider_free_and_works_without_external_snapshot() -> None:
    bundle, _, external = _snapshot()
    with patch("requests.Session.get", side_effect=AssertionError("provider access forbidden")):
        cached = build_presenter_deck(bundle, external_snapshot=external)
        local_only = build_presenter_deck(bundle)
    assert cached.offline_ready and local_only.offline_ready
    assert local_only.external_context_snapshot_id is None
    literature = local_only.sections[PRESENTER_ORDER.index(PresenterSectionKind.LITERATURE)]
    assert literature.observation == "No cached literature context is attached."
    assert literature.literature_context == ()


def test_presenter_rejects_external_snapshot_from_another_result() -> None:
    bundle, _, external = _snapshot()
    # Integrity validation rejects tampering even before the presenter can use it.
    with pytest.raises(ValueError, match="integrity"):
        replace(external, simulation_snapshot_id="other")
    other_bundle = SimpleNamespace(
        snapshot=SimpleNamespace(
            simulation_id=bundle.simulation_id,
            snapshot_id="other",
        )
    )
    with pytest.raises(ValueError, match="does not belong"):
        build_presenter_deck(other_bundle, external_snapshot=external)


def _result(provider: str, status: ProviderStatus) -> ProviderResult:
    request = ProviderRequest(provider, "fixture", 9606, "fixture", "1", "1")
    return ProviderResult(request, status, message="technical exception details must not render")


def test_external_ui_state_keeps_analysis_local_cached_and_live_states_separate() -> None:
    bundle, _, external = _snapshot()
    state = build_external_context_ui_state(
        bundle,
        cached_snapshot=external,
        live_results=(_result("uniprot", ProviderStatus.TIMEOUT), _result("reactome", ProviderStatus.OFFLINE)),
        live_requested=True,
    )
    assert state.sophiark_analysis is ExternalDisplayState.AVAILABLE
    assert state.local_research_context is ExternalDisplayState.AVAILABLE
    assert state.cached_external_context is ExternalDisplayState.CACHED
    assert state.live_external_context is ExternalDisplayState.TIMEOUT
    providers = {item.provider: item for item in state.providers}
    assert providers["uniprot"].state is ExternalDisplayState.TIMEOUT
    assert providers["uniprot"].captured_statuses == ("AVAILABLE", "TIMEOUT")
    assert providers["europepmc"].state is ExternalDisplayState.CACHED
    assert providers["reactome"].state is ExternalDisplayState.OFFLINE
    assert "exception" not in repr(state).casefold()


def test_offline_live_status_never_degrades_sophiark_or_local_context() -> None:
    bundle, _, _ = _snapshot()
    state = build_external_context_ui_state(bundle, live_requested=True)
    assert state.sophiark_analysis is ExternalDisplayState.AVAILABLE
    assert state.local_research_context is ExternalDisplayState.AVAILABLE
    assert state.cached_external_context is ExternalDisplayState.NOT_REQUESTED
    assert state.live_external_context is ExternalDisplayState.OFFLINE
    assert "remain available" in state.notices[-1]


def test_injected_renderers_render_only_prebuilt_models_without_tabs_or_new_kpis() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle, external_snapshot=external)
    state = build_external_context_ui_state(bundle, cached_snapshot=external)
    ui = FakeUI()
    render_presenter_mode(deck, ui=ui)
    render_external_context_status(state, ui=ui)
    methods = [name for name, _, __ in ui.calls]
    assert methods.count("expander") == 8
    assert "tabs" not in methods
    assert "metric" not in methods
    assert "SOPHIARK COMPUTED" in ui.text
    assert "Cached External Context: Cached" in ui.text
    assert "technical exception details" not in ui.text
