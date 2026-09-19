"""Phase 11 acceptance tests for separate structured Research Brief exports."""

from __future__ import annotations

from dataclasses import FrozenInstanceError, replace
from unittest.mock import patch

import pytest

from src.research.brief import (
    BRIEF_SECTION_ORDER,
    BriefFormat,
    ResearchBrief,
    build_research_brief,
    render_research_brief_download,
)
from src.research.presenter import build_presenter_deck
from tests.test_research_external_context_snapshot import _snapshot


class FakeUI:
    def __init__(self):
        self.calls = []

    def download_button(self, *args, **kwargs):
        self.calls.append((args, kwargs))


def test_brief_has_exact_structured_order_and_separate_snapshot_identity() -> None:
    bundle, _, external = _snapshot()
    brief = build_research_brief(bundle, external_snapshot=external)
    assert isinstance(brief, ResearchBrief)
    assert brief.brief_id.startswith("research-brief-")
    assert brief.simulation_id == bundle.simulation_id
    assert brief.simulation_snapshot_id == bundle.snapshot_id
    assert brief.external_context_snapshot_id == external.context_snapshot_id
    assert tuple(section.title for section in brief.sections) == BRIEF_SECTION_ORDER
    assert len(brief.sections) == 10
    assert all(section.entries for section in brief.sections)
    assert any(
        entry.source_badge == "SOPHIARK COMPUTED"
        for section in brief.sections
        for entry in section.entries
    )
    assert any(
        entry.source_badge == "LITERATURE · EUROPE PMC"
        for section in brief.sections
        for entry in section.entries
    )
    with pytest.raises(FrozenInstanceError):
        brief.title = "changed"


def test_brief_is_deterministic_and_does_not_mutate_or_embed_raw_scientific_tables() -> None:
    bundle, _, external = _snapshot()
    before_report = bundle.snapshot.report.rows_json
    before_signed = bundle.snapshot.signed_redistribution.rows_json
    first = build_research_brief(bundle, external_snapshot=external)
    second = build_research_brief(bundle, external_snapshot=external)
    assert first == second
    assert first.brief_id == second.brief_id
    assert bundle.snapshot.report.rows_json == before_report
    assert bundle.snapshot.signed_redistribution.rows_json == before_signed
    representation = repr(first)
    assert "rows_json" not in representation
    assert "scientific_fingerprint" not in representation
    assert not hasattr(first, "score")


def test_markdown_html_and_text_exports_preserve_boundaries_and_escape_html() -> None:
    bundle, _, external = _snapshot()
    brief = build_research_brief(bundle, external_snapshot=external)
    markdown = brief.render(BriefFormat.MARKDOWN)
    html = brief.render(BriefFormat.HTML)
    text = brief.render(BriefFormat.TEXT)
    assert markdown.startswith("# SOPHIARK Research Brief")
    assert "## Literature Context" in markdown
    assert "> Interpretation boundary:" in markdown
    assert html.startswith("<!doctype html>")
    assert '<meta charset="utf-8">' in html
    assert "<section><h2>Sources</h2>" in html
    assert "<script" not in html.casefold()
    assert "INTERPRETATION BOUNDARIES" in text
    assert external.context_snapshot_id in markdown


@pytest.mark.parametrize(
    ("format", "extension", "mime"),
    (
        (BriefFormat.MARKDOWN, ".md", "text/markdown"),
        (BriefFormat.HTML, ".html", "text/html"),
        (BriefFormat.TEXT, ".txt", "text/plain"),
    ),
)
def test_download_payload_is_utf8_and_separate_from_scientific_exports(format, extension, mime) -> None:
    bundle, _, external = _snapshot()
    brief = build_research_brief(bundle, external_snapshot=external)
    filename, content_type, data = brief.download_payload(format)
    assert filename.startswith("research-brief-")
    assert filename.endswith(extension)
    assert content_type.startswith(mime)
    assert data.decode("utf-8") == brief.render(format)
    assert not filename.endswith((".csv", ".xlsx"))


def test_brief_build_is_provider_free_and_local_only_export_remains_available() -> None:
    bundle, _, external = _snapshot()
    with patch("requests.Session.get", side_effect=AssertionError("network forbidden")):
        cached = build_research_brief(bundle, external_snapshot=external)
        local = build_research_brief(bundle)
    assert cached.external_context_snapshot_id == external.context_snapshot_id
    assert local.external_context_snapshot_id is None
    literature = next(section for section in local.sections if section.title == "Literature Context")
    assert literature.entries[0].value == "No cached literature context is attached."


def test_brief_rejects_stale_presenter_deck() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle, external_snapshot=external)
    stale = replace(deck, simulation_snapshot_id="other")
    with pytest.raises(ValueError, match="stale"):
        build_research_brief(bundle, presenter_deck=stale)


def test_brief_rejects_mismatched_presenter_and_external_snapshot() -> None:
    bundle, _, external = _snapshot()
    deck = build_presenter_deck(bundle)
    with pytest.raises(ValueError, match="do not match"):
        build_research_brief(
            bundle,
            presenter_deck=deck,
            external_snapshot=external,
        )


def test_download_renderer_uses_only_explicit_brief_payload() -> None:
    bundle, _, external = _snapshot()
    brief = build_research_brief(bundle, external_snapshot=external)
    ui = FakeUI()
    render_research_brief_download(brief, ui=ui, format=BriefFormat.HTML)
    assert len(ui.calls) == 1
    args, kwargs = ui.calls[0]
    assert args == ("Download Research Brief",)
    assert kwargs["file_name"].endswith(".html")
    assert kwargs["mime"].startswith("text/html")
    assert kwargs["data"] == brief.render(BriefFormat.HTML).encode("utf-8")
