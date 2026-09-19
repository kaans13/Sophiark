"""Phase 12 acceptance tests for the future interpreter boundary only."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from pathlib import Path
from unittest.mock import patch

import pytest

from src.research.interpreter import (
    DeterministicInterpreter,
    INTERPRETER_SECTION_ORDER,
    InterpreterContext,
    InterpreterContextItem,
    ResearchInterpreter,
    build_interpreter_context,
)
from src.research.presenter import build_presenter_deck
from tests.test_research_external_context_snapshot import _snapshot


def test_interpreter_is_abstract_and_default_is_deterministic() -> None:
    with pytest.raises(TypeError):
        ResearchInterpreter()
    bundle, _, external = _snapshot()
    context = build_interpreter_context(bundle, external_snapshot=external)
    first = DeterministicInterpreter().render(context)
    second = DeterministicInterpreter().render(context)
    assert first == second
    assert first.context_id == context.context_id
    assert tuple(section.kind for section in first.sections) == INTERPRETER_SECTION_ORDER
    assert first.interpreter == "deterministic"


def test_context_has_exact_eight_bounded_sections_and_provenance() -> None:
    bundle, _, external = _snapshot()
    context = build_interpreter_context(bundle, external_snapshot=external)
    assert isinstance(context, InterpreterContext)
    assert context.context_id.startswith("interpreter-context-")
    assert context.simulation_id == bundle.simulation_id
    assert context.simulation_snapshot_id == bundle.snapshot_id
    assert context.external_context_snapshot_id == external.context_snapshot_id
    assert tuple(kind for kind, _ in context.section_items()) == INTERPRETER_SECTION_ORDER
    assert len(context.section_items()) == 8
    assert len(context.selected_observations) <= 12
    assert len(context.interpretation_contracts) <= 16
    assert len(context.evidence_graph_subset) <= 12
    assert len(context.external_assertions) <= 12
    assert len(context.relevant_publications) <= 12
    assert any(
        item.source_badge.startswith("LITERATURE")
        for item in context.relevant_publications
    )
    assert any(
        item.source_identifier == external.context_snapshot_id
        for item in context.provenance
    )
    with pytest.raises(FrozenInstanceError):
        context.context_id = "changed"


def test_context_builder_never_materializes_raw_scientific_tables_or_mutates_bundle() -> None:
    bundle, _, external = _snapshot()
    report_json = bundle.snapshot.report.rows_json
    signed_json = bundle.snapshot.signed_redistribution.rows_json
    with patch.object(
        type(bundle.snapshot.report),
        "records",
        side_effect=AssertionError("raw table materialization forbidden"),
    ):
        context = build_interpreter_context(bundle, external_snapshot=external)
    assert bundle.snapshot.report.rows_json == report_json
    assert bundle.snapshot.signed_redistribution.rows_json == signed_json
    assert "rows_json" not in repr(context)
    assert "scientific_fingerprint" not in repr(context)
    assert not hasattr(context, "report")
    assert not hasattr(context, "score")


def test_context_is_provider_free_and_local_only_path_remains_complete() -> None:
    bundle, _, external = _snapshot()
    with patch("requests.Session.get", side_effect=AssertionError("network forbidden")):
        cached = build_interpreter_context(bundle, external_snapshot=external)
        local = build_interpreter_context(bundle)
        rendered = DeterministicInterpreter().render(local)
    assert cached.external_assertions
    assert local.external_context_snapshot_id is None
    assert local.external_assertions == ()
    assert local.relevant_publications == ()
    assert len(rendered.sections) == 8


def test_context_rejects_stale_or_mismatched_evidence() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle, external_snapshot=external)
    with pytest.raises(ValueError, match="stale"):
        build_interpreter_context(
            bundle,
            presenter_deck=replace(deck, simulation_snapshot_id="other"),
        )
    local_deck = build_presenter_deck(bundle)
    with pytest.raises(ValueError, match="do not match"):
        build_interpreter_context(
            bundle,
            presenter_deck=local_deck,
            external_snapshot=external,
        )


def test_context_item_rejects_raw_table_identifiers() -> None:
    for identifier in ("rows_json", "raw_report", "Raw Report"):
        with pytest.raises(ValueError, match="raw scientific"):
            InterpreterContextItem(identifier, "value", "SOPHIARK COMPUTED")


def test_simulation_context_uses_explicit_missing_labels_instead_of_guessing() -> None:
    bundle, _, _ = _snapshot()
    empty_scope_snapshot = replace(bundle.snapshot, tissue="", targets=())
    empty_scope_bundle = replace(bundle, snapshot=empty_scope_snapshot)
    context = build_interpreter_context(empty_scope_bundle)
    values = {item.identifier: item.text for item in context.simulation_context}
    assert values["tissue"] == "Not recorded"
    assert values["targets"] == "None recorded"


def test_interface_module_has_no_model_sdk_or_transport_imports() -> None:
    source = Path("src/research/interpreter.py").read_text(encoding="utf-8")
    import_lines = "\n".join(
        line.casefold()
        for line in source.splitlines()
        if line.startswith("import ") or line.startswith("from ")
    )
    for forbidden in (
        "openai", "anthropic", "gemini", "ollama", "requests", "httpx", "urllib"
    ):
        assert forbidden not in import_lines
