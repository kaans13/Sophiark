"""Structured, deterministic Research Brief exports.

Briefs are separate artifacts.  They never modify the scientific CSV/XLSX
exports and require no PDF library, provider, network, LLM, or file write.
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
import hashlib
import html
import json
from typing import Any, Iterable

from .external_context import ExternalContextSnapshot
from .presenter import PresenterDeck, PresenterSectionKind, build_presenter_deck


class BriefFormat(str, Enum):
    MARKDOWN = "markdown"
    HTML = "html"
    TEXT = "text"


BRIEF_SECTION_ORDER = (
    "Analysis Context",
    "Network Response",
    "Detected Observations",
    "Family Context",
    "Functional Context",
    "Relationship Context",
    "Statistical Support",
    "Literature Context",
    "Interpretation Boundaries",
    "Sources",
)


@dataclass(frozen=True, slots=True)
class BriefEntry:
    label: str
    value: str
    source_badge: str | None = None
    source_identifier: str | None = None


@dataclass(frozen=True, slots=True)
class BriefSection:
    title: str
    entries: tuple[BriefEntry, ...]
    note: str | None = None

    def __post_init__(self) -> None:
        object.__setattr__(self, "entries", tuple(self.entries))


@dataclass(frozen=True, slots=True)
class ResearchBrief:
    brief_id: str
    simulation_id: str
    simulation_snapshot_id: str
    external_context_snapshot_id: str | None
    title: str
    sections: tuple[BriefSection, ...]

    def __post_init__(self) -> None:
        object.__setattr__(self, "sections", tuple(self.sections))
        if tuple(section.title for section in self.sections) != BRIEF_SECTION_ORDER:
            raise ValueError("Research Brief sections must follow the fixed structured order")
        if self.brief_id != _brief_id(
            self.simulation_id,
            self.simulation_snapshot_id,
            self.external_context_snapshot_id,
            self.title,
            self.sections,
        ):
            raise ValueError("Research Brief integrity check failed")

    def render(self, format: BriefFormat | str = BriefFormat.MARKDOWN) -> str:
        selected = BriefFormat(format)
        if selected is BriefFormat.MARKDOWN:
            return _markdown(self)
        if selected is BriefFormat.HTML:
            return _html(self)
        return _text(self)

    def download_payload(
        self,
        format: BriefFormat | str = BriefFormat.MARKDOWN,
    ) -> tuple[str, str, bytes]:
        selected = BriefFormat(format)
        extension, mime = {
            BriefFormat.MARKDOWN: ("md", "text/markdown; charset=utf-8"),
            BriefFormat.HTML: ("html", "text/html; charset=utf-8"),
            BriefFormat.TEXT: ("txt", "text/plain; charset=utf-8"),
        }[selected]
        safe_simulation = "".join(
            character if character.isalnum() or character in "-_" else "-"
            for character in self.simulation_id
        )
        return (
            f"research-brief-{safe_simulation}.{extension}",
            mime,
            self.render(selected).encode("utf-8"),
        )


def _entry_value(entry: BriefEntry) -> str:
    suffix = f" [{entry.source_badge}]" if entry.source_badge else ""
    identifier = f" ({entry.source_identifier})" if entry.source_identifier else ""
    return f"{entry.value}{suffix}{identifier}"


def _plain_sections(sections: Iterable[BriefSection]) -> list[dict[str, Any]]:
    return [
        {
            "title": section.title,
            "entries": [
                {
                    "label": entry.label,
                    "value": entry.value,
                    "source_badge": entry.source_badge,
                    "source_identifier": entry.source_identifier,
                }
                for entry in section.entries
            ],
            "note": section.note,
        }
        for section in sections
    ]


def _brief_id(
    simulation_id: str,
    simulation_snapshot_id: str,
    external_context_snapshot_id: str | None,
    title: str,
    sections: Iterable[BriefSection],
) -> str:
    encoded = json.dumps(
        {
            "simulation_id": simulation_id,
            "simulation_snapshot_id": simulation_snapshot_id,
            "external_context_snapshot_id": external_context_snapshot_id,
            "title": title,
            "sections": _plain_sections(sections),
        },
        ensure_ascii=False,
        allow_nan=False,
        sort_keys=True,
        separators=(",", ":"),
    ).encode("utf-8")
    return "research-brief-" + hashlib.sha256(encoded).hexdigest()[:32]


def _facts(section: Any) -> tuple[BriefEntry, ...]:
    return tuple(
        BriefEntry(
            fact.label,
            fact.value,
            fact.source_badge,
            fact.source_identifier,
        )
        for fact in section.sophiark_data
    )


def _context_entries(items: Iterable[Any], label: str) -> tuple[BriefEntry, ...]:
    return tuple(
        BriefEntry(label, item.text, item.source_badge, item.source_identifier)
        for item in items
    )


def build_research_brief(
    bundle: Any,
    *,
    presenter_deck: PresenterDeck | None = None,
    external_snapshot: ExternalContextSnapshot | None = None,
) -> ResearchBrief:
    deck = presenter_deck or build_presenter_deck(
        bundle,
        external_snapshot=external_snapshot,
    )
    if deck.simulation_id != bundle.snapshot.simulation_id:
        raise ValueError("presenter deck does not belong to the current simulation")
    if deck.simulation_snapshot_id != bundle.snapshot.snapshot_id:
        raise ValueError("presenter deck is stale")
    if (
        presenter_deck is not None
        and external_snapshot is not None
        and deck.external_context_snapshot_id != external_snapshot.context_snapshot_id
    ):
        raise ValueError("presenter deck and external context snapshot do not match")
    by_kind = {section.kind: section for section in deck.sections}
    perturbation = by_kind[PresenterSectionKind.PERTURBATION]
    network = by_kind[PresenterSectionKind.NETWORK_RESPONSE]
    structural = by_kind[PresenterSectionKind.STRUCTURAL_PATTERNS]
    functional = by_kind[PresenterSectionKind.FUNCTIONAL_CONTEXT]
    relationships = by_kind[PresenterSectionKind.RELATIONSHIPS]
    statistics = by_kind[PresenterSectionKind.STATISTICAL_SUPPORT]
    literature = by_kind[PresenterSectionKind.LITERATURE]
    limitations = by_kind[PresenterSectionKind.LIMITATIONS]

    family_entries = tuple(
        entry for entry in _facts(structural)
        if "family" in entry.label.casefold()
    )
    observation_entries = (
        _facts(structural)
        + _context_entries(structural.biological_context, "Biological context")
    )
    sections = (
        BriefSection(
            "Analysis Context",
            _facts(perturbation),
            perturbation.interpretation_boundary,
        ),
        BriefSection(
            "Network Response",
            _facts(network),
            network.interpretation_boundary,
        ),
        BriefSection(
            "Detected Observations",
            observation_entries,
            structural.interpretation_boundary,
        ),
        BriefSection(
            "Family Context",
            family_entries or (BriefEntry("Status", "No family pattern is included in this brief."),),
            "Family co-occurrence does not establish coordinated activation or causality.",
        ),
        BriefSection(
            "Functional Context",
            _facts(functional)
            + _context_entries(functional.biological_context, "Local context")
            + _context_entries(functional.external_context, "External context"),
            functional.interpretation_boundary,
        ),
        BriefSection(
            "Relationship Context",
            _facts(relationships)
            + _context_entries(relationships.external_context, "External context"),
            relationships.interpretation_boundary,
        ),
        BriefSection(
            "Statistical Support",
            _facts(statistics),
            statistics.interpretation_boundary,
        ),
        BriefSection(
            "Literature Context",
            _context_entries(literature.literature_context, "Literature context")
            or (BriefEntry("Status", literature.observation, "LITERATURE"),),
            literature.interpretation_boundary,
        ),
        BriefSection(
            "Interpretation Boundaries",
            tuple(
                BriefEntry(section.title, section.interpretation_boundary)
                for section in deck.sections
            ),
            limitations.interpretation_boundary,
        ),
        BriefSection(
            "Sources",
            tuple(BriefEntry("Source badge", source) for source in deck.source_badges),
            (
                f"External context snapshot: {deck.external_context_snapshot_id}"
                if deck.external_context_snapshot_id
                else "No external context snapshot is attached."
            ),
        ),
    )
    title = f"SOPHIARK Research Brief · {bundle.snapshot.tissue} · {bundle.snapshot.species}"
    return ResearchBrief(
        brief_id=_brief_id(
            bundle.snapshot.simulation_id,
            bundle.snapshot.snapshot_id,
            deck.external_context_snapshot_id,
            title,
            sections,
        ),
        simulation_id=bundle.snapshot.simulation_id,
        simulation_snapshot_id=bundle.snapshot.snapshot_id,
        external_context_snapshot_id=deck.external_context_snapshot_id,
        title=title,
        sections=sections,
    )


def _markdown(brief: ResearchBrief) -> str:
    lines = [f"# {brief.title}", "", f"Brief ID: `{brief.brief_id}`", ""]
    for section in brief.sections:
        lines.extend((f"## {section.title}", ""))
        for entry in section.entries:
            label = entry.label.replace("|", "\\|")
            value = _entry_value(entry).replace("|", "\\|")
            lines.append(f"- **{label}:** {value}")
        if section.note:
            lines.extend(("", f"> Interpretation boundary: {section.note}"))
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def _html(brief: ResearchBrief) -> str:
    parts = [
        "<!doctype html><html><head><meta charset=\"utf-8\">",
        f"<title>{html.escape(brief.title)}</title></head><body>",
        f"<h1>{html.escape(brief.title)}</h1>",
        f"<p>Brief ID: <code>{html.escape(brief.brief_id)}</code></p>",
    ]
    for section in brief.sections:
        parts.append(f"<section><h2>{html.escape(section.title)}</h2><ul>")
        for entry in section.entries:
            parts.append(
                f"<li><strong>{html.escape(entry.label)}:</strong> "
                f"{html.escape(_entry_value(entry))}</li>"
            )
        parts.append("</ul>")
        if section.note:
            parts.append(
                '<p class="interpretation-boundary"><strong>Interpretation boundary:</strong> '
                f"{html.escape(section.note)}</p>"
            )
        parts.append("</section>")
    parts.append("</body></html>")
    return "".join(parts)


def _text(brief: ResearchBrief) -> str:
    lines = [brief.title, f"Brief ID: {brief.brief_id}", ""]
    for section in brief.sections:
        lines.extend((section.title.upper(), "-" * len(section.title)))
        lines.extend(f"{entry.label}: {_entry_value(entry)}" for entry in section.entries)
        if section.note:
            lines.append(f"Interpretation boundary: {section.note}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def render_research_brief_download(
    brief: ResearchBrief,
    *,
    ui: Any,
    format: BriefFormat | str = BriefFormat.MARKDOWN,
) -> None:
    filename, mime, data = brief.download_payload(format)
    ui.download_button(
        "Download Research Brief",
        data=data,
        file_name=filename,
        mime=mime,
    )


__all__ = [
    "BRIEF_SECTION_ORDER",
    "BriefEntry",
    "BriefFormat",
    "BriefSection",
    "ResearchBrief",
    "build_research_brief",
    "render_research_brief_download",
]
