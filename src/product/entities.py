"""Canonical biological identity helpers outside calculation engines."""

from __future__ import annotations

from dataclasses import dataclass
import re
from typing import Iterable, Mapping


_ENSP = re.compile(r"^ENSP\d+$")
_STRING_ENSP = re.compile(r"^(?P<taxon>\d+)\.(?P<ensp>ENSP\d+)$")


def normalize_string_protein_id(value: object, *, taxon_id: int = 9606) -> str:
    """Return canonical ENSP, removing one matching STRING species prefix.

    The operation is deliberately strict and idempotent. A foreign taxon or an
    unknown identifier shape fails instead of being silently rewritten.
    """
    text = str(value).strip()
    if _ENSP.fullmatch(text):
        return text
    matched = _STRING_ENSP.fullmatch(text)
    if not matched:
        raise ValueError(f"unsupported STRING protein identifier: {text!r}")
    if int(matched.group("taxon")) != int(taxon_id):
        raise ValueError(
            f"STRING identifier taxon {matched.group('taxon')} does not match {taxon_id}"
        )
    return matched.group("ensp")


def normalize_ensp(value: object) -> str:
    return normalize_string_protein_id(value, taxon_id=9606)


@dataclass(frozen=True, slots=True)
class MappingReport:
    source_records: int
    mapped_records: int
    unmapped_records: int
    usable_records: int

    @property
    def coverage(self) -> float:
        return self.mapped_records / self.source_records if self.source_records else 0.0

    def as_dict(self) -> dict[str, int | float]:
        return {
            "source_records": self.source_records,
            "mapped_records": self.mapped_records,
            "unmapped_records": self.unmapped_records,
            "usable_records": self.usable_records,
            "mapping_coverage": self.coverage,
        }


def map_entities(values: Iterable[object], mapping: Mapping[str, str]) -> tuple[list[str | None], MappingReport]:
    resolved: list[str | None] = []
    mapped = 0
    usable = 0
    source = 0
    for raw in values:
        source += 1
        try:
            canonical = normalize_ensp(raw)
        except ValueError:
            resolved.append(None)
            continue
        usable += 1
        value = mapping.get(canonical)
        resolved.append(value)
        mapped += value is not None
    return resolved, MappingReport(source, mapped, source - mapped, usable)
