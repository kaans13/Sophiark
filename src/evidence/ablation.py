"""Frozen-topology, channel-level diagnostics for the Evidence BETA control.

Nothing in this module is selected by the application.  It is deliberately
limited to experimental ablations that retain the E1 edge set.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

from .scores import combine_string_channels


NON_TEXT_CHANNELS: tuple[str, ...] = (
    "neighborhood", "neighborhood_transferred", "fusion", "cooccurence",
    "coexpression", "coexpression_transferred", "experiments",
    "experiments_transferred", "database", "database_transferred",
)
DIRECT_CHANNELS: tuple[str, ...] = (
    "neighborhood", "fusion", "cooccurence", "coexpression", "experiments", "database",
)
TRANSFERRED_CHANNELS: tuple[str, ...] = (
    "neighborhood_transferred", "coexpression_transferred",
    "experiments_transferred", "database_transferred",
)


@dataclass(frozen=True, slots=True)
class AblationCondition:
    identifier: str
    description: str
    included_channels: tuple[str, ...]
    family: str

    @property
    def excluded_channels(self) -> tuple[str, ...]:
        return tuple(channel for channel in NON_TEXT_CHANNELS if channel not in self.included_channels)


def _without(*channels: str) -> tuple[str, ...]:
    removed = set(channels)
    return tuple(channel for channel in NON_TEXT_CHANNELS if channel not in removed)


def ablation_conditions() -> tuple[AblationCondition, ...]:
    """Pre-registered diagnostic panel; duplicate biological groups are retained as aliases."""
    all_channels = NON_TEXT_CHANNELS
    experimental = ("experiments", "experiments_transferred")
    coexpression = ("coexpression", "coexpression_transferred")
    database = ("database", "database_transferred")
    genomic = ("neighborhood", "neighborhood_transferred", "fusion", "cooccurence")
    rows = (
        ("C1", "Current Evidence: all non-text channels", all_channels, "control"),
        ("A1", "minus database", _without("database"), "single_drop"),
        ("A2", "minus database_transferred", _without("database_transferred"), "single_drop"),
        ("A3", "minus database + database_transferred", _without(*database), "paired_drop"),
        ("A4", "minus experiments", _without("experiments"), "single_drop"),
        ("A5", "minus experiments_transferred", _without("experiments_transferred"), "single_drop"),
        ("A6", "minus experiments + experiments_transferred", _without(*experimental), "paired_drop"),
        ("A7", "minus coexpression", _without("coexpression"), "single_drop"),
        ("A8", "minus coexpression_transferred", _without("coexpression_transferred"), "single_drop"),
        ("A9", "minus coexpression + coexpression_transferred", _without(*coexpression), "paired_drop"),
        ("A10", "minus neighborhood", _without("neighborhood"), "single_drop"),
        ("A11", "minus neighborhood_transferred", _without("neighborhood_transferred"), "single_drop"),
        ("A12", "minus both neighborhood forms", _without("neighborhood", "neighborhood_transferred"), "paired_drop"),
        ("A13", "minus fusion", _without("fusion"), "single_drop"),
        ("A14", "minus cooccurence", _without("cooccurence"), "single_drop"),
        ("B1", "direct non-text only", DIRECT_CHANNELS, "direct_vs_transferred"),
        ("B2", "transferred non-text only", TRANSFERRED_CHANNELS, "direct_vs_transferred"),
        ("G1", "remove curated database evidence", _without(*database), "family_drop"),
        ("G2", "remove experimental evidence", _without(*experimental), "family_drop"),
        ("G3", "remove coexpression evidence", _without(*coexpression), "family_drop"),
        ("G4", "remove evolutionary/genomic context", _without(*genomic), "family_drop"),
        ("P1", "experimental only", experimental, "retention_only"),
        ("P2", "coexpression only", coexpression, "retention_only"),
        ("P3", "database only", database, "retention_only"),
        ("P4", "experimental + coexpression", experimental + coexpression, "retention_only"),
        ("P5", "experimental + coexpression + genomic, excluding database", experimental + coexpression + genomic, "retention_only"),
        ("P6", "direct experimental + direct coexpression", ("experiments", "coexpression"), "retention_only"),
    )
    return tuple(AblationCondition(identifier, description, tuple(channels), family) for identifier, description, channels, family in rows)


def recombine_channels(profiles: pd.DataFrame, included_channels: Sequence[str], *, prior: float = .041) -> pd.Series:
    """Recombine channel inputs exactly as the validated STRING combiner does."""
    selected = tuple(included_channels)
    unknown = set(selected).difference(NON_TEXT_CHANNELS)
    if unknown:
        raise ValueError(f"Unknown non-text channel(s): {sorted(unknown)}")
    missing = set(selected).difference(profiles.columns)
    if missing:
        raise ValueError(f"Source profiles lack requested channel(s): {sorted(missing)}")
    if not selected:
        raise ValueError("An ablation must retain at least one evidence channel")
    numeric = profiles.loc[:, selected].apply(pd.to_numeric, errors="coerce").fillna(0.0).to_numpy(dtype=float) / 1000.0
    return pd.Series(
        [combine_string_channels(row, prior) for row in numeric], index=profiles.index, name="ablated_confidence",
    )


def apply_e1_ablation(graph, profiles: pd.DataFrame, condition: AblationCondition, *, prior: float = .041) -> pd.Series:
    """Replace only weights on an E1 copy; nodes and edges are never removed."""
    if len(profiles) != graph.ecount():
        raise ValueError("Profile/graph edge count mismatch")
    if "classic_weight" not in graph.es.attribute_names():
        raise ValueError("E1 ablation requires an Evidence graph with preserved classic weights")
    confidence = recombine_channels(profiles, condition.included_channels, prior=prior)
    official = pd.to_numeric(profiles["combined_score"], errors="coerce").to_numpy(dtype=float) / 1000.0
    classic_weight = np.asarray(graph.es["classic_weight"], dtype=float)
    tissue_multiplier = classic_weight / np.maximum(np.square(official), 1e-12)
    new_weight = np.maximum(np.square(confidence.to_numpy(dtype=float)) * tissue_multiplier, 1e-9)
    graph.es["ablated_confidence"] = confidence.tolist()
    graph.es["weight"] = new_weight.tolist()
    graph.es["distance"] = (1.0 / new_weight).tolist()
    # Diagnostic conditions must not accidentally reactivate optional modifiers.
    graph.es["text_bonus"] = [1.0] * graph.ecount()
    graph.es["physical_bonus"] = [1.0] * graph.ecount()
    graph["ablation_condition"] = condition.identifier
    graph["ablation_included_channels"] = list(condition.included_channels)
    return confidence

