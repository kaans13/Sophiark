from __future__ import annotations

import pandas as pd
import pytest

from src.evidence.ablation import (
    DIRECT_CHANNELS, NON_TEXT_CHANNELS, TRANSFERRED_CHANNELS,
    ablation_conditions, recombine_channels,
)
from src.evidence.scores import combine_string_channels


def _profiles() -> pd.DataFrame:
    return pd.DataFrame({channel: [0, 100, 700] for channel in NON_TEXT_CHANNELS})


def test_control_recombines_all_non_text_channels_with_validated_combiner():
    profiles = _profiles()
    value = recombine_channels(profiles, NON_TEXT_CHANNELS)
    assert value.iloc[2] == pytest.approx(combine_string_channels([.7] * len(NON_TEXT_CHANNELS)))


def test_drop_condition_excludes_only_its_requested_channel():
    profiles = _profiles()
    condition = next(row for row in ablation_conditions() if row.identifier == "A1")
    value = recombine_channels(profiles, condition.included_channels)
    expected = combine_string_channels([.7] * (len(NON_TEXT_CHANNELS) - 1))
    assert "database" not in condition.included_channels
    assert "experiments" in condition.included_channels
    assert value.iloc[2] == pytest.approx(expected)


def test_direct_and_transferred_diagnostics_are_disjoint_and_complete():
    assert set(DIRECT_CHANNELS).isdisjoint(TRANSFERRED_CHANNELS)
    assert set(DIRECT_CHANNELS) | set(TRANSFERRED_CHANNELS) == set(NON_TEXT_CHANNELS)


def test_conditions_never_include_text_or_physical_modifiers():
    forbidden = {"textmining", "textmining_transferred", "physical_score", "physical_status", "homology"}
    assert all(not (set(condition.included_channels) & forbidden) for condition in ablation_conditions())
