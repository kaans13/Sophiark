from __future__ import annotations

import numpy as np
import pandas as pd

from src.unified import clear_unified_cache, run_unified_analysis
from src.unified.service import PresentationConvention


def _classic() -> pd.DataFrame:
    return pd.DataFrame({
        "gene": ["ENSP_A", "ENSP_B", "ENSP_C"],
        "PageRank_Baseline": [.4, .3, .2], "PageRank_Perturbed": [.44, .285, .2],
        "Delta_PageRank_Pct": [10.0, -5.0, 0.0], "Hinterland_Skoru": [9.0, 7.0, 2.0],
    })


def _directed() -> pd.DataFrame:
    return pd.DataFrame({
        "entity": ["ENSP_A", "ENSP_B", "ENSP_D"],
        "Directed_PageRank_Baseline": [.4, .3, .1], "Directed_PageRank_Perturbed": [.436, .282, .11],
        "Directed_Redistribution_Pct": [9.0, -6.0, 10.0], "Directed_BC": [1.0, 2.0, 3.0],
    })


def _runner(**_):
    return _classic(), _directed(), {"status": "AVAILABLE", "direction_coverage": .5}, {
        "ENSP_A": "A", "ENSP_B": "B", "ENSP_C": "C", "ENSP_D": "D",
    }, {}


def test_unified_preserves_standalone_raw_engine_outputs_and_ranks():
    classic = _classic(); directed = _directed()
    report = run_unified_analysis(
        target="ENSP_A", tissue="Lung", convention=PresentationConvention(top_n=2),
        use_cache=False, engine_runner=_runner,
    )
    pd.testing.assert_frame_equal(report.classic_raw, classic)
    pd.testing.assert_frame_equal(report.directed_raw, directed)
    candidate = report.candidates.set_index("entity_id")
    assert candidate.loc["ENSP_A", "classic_absolute_response_rank"] == 1
    assert candidate.loc["ENSP_A", "directed_absolute_response_rank"] == 2
    assert candidate.loc["ENSP_A", "classic_absolute_rank_percentile"] == 1.0
    assert candidate.loc["ENSP_B", "classic_positive_rank"] is pd.NA
    assert candidate.loc["ENSP_B", "directed_negative_rank"] == 1
    assert candidate.loc["ENSP_A", "finding_class"] == "ROBUST_CROSS_MODEL"


def test_missing_entity_is_not_zero_and_engine_failure_is_not_substituted():
    report = run_unified_analysis(target="ENSP_A", tissue="Lung", use_cache=False, engine_runner=_runner)
    candidate = report.candidates.set_index("entity_id")
    assert bool(candidate.loc["ENSP_C", "classic_available"])
    assert not bool(candidate.loc["ENSP_C", "directed_available"])
    assert pd.isna(candidate.loc["ENSP_C", "directed_response"])
    assert candidate.loc["ENSP_C", "finding_class"] == "CLASSIC_ONLY_AVAILABLE"
    assert bool(candidate.loc["ENSP_D", "directed_available"])
    assert pd.isna(candidate.loc["ENSP_D", "classic_response"])


def test_context_layers_do_not_change_synthesis_ranking(monkeypatch):
    from src.unified import service

    def provenance(frame, *, target):
        changed = frame.copy(deep=True); changed["Evidence_Provenance_Available"] = False
        return changed, {"status": "AVAILABLE"}

    def complexes(*, candidates, **_):
        changed = candidates.copy(deep=True); changed["CORUM_Target_Complex_Member"] = False
        return changed, pd.DataFrame(), {"status": "AVAILABLE"}

    monkeypatch.setattr(service, "_target_edge_provenance", provenance)
    monkeypatch.setattr(service, "_complex_context", complexes)
    report = run_unified_analysis(target="ENSP_A", tissue="Lung", use_cache=False, engine_runner=_runner)
    values = report.candidates.set_index("entity_id")
    assert values.loc["ENSP_A", "classic_absolute_response_rank"] == 1
    assert values.loc["ENSP_A", "directed_absolute_response_rank"] == 2
    assert np.isclose(values.loc["ENSP_A", "classic_response"], 10.0)
    assert np.isclose(values.loc["ENSP_A", "directed_response"], 9.0)


def test_unified_cache_is_keyed_by_target_and_returns_detached_frames():
    from src.unified import service
    calls: list[str] = []

    def cached_runner(**kwargs):
        calls.append(kwargs["target"])
        return _runner()

    clear_unified_cache()
    original = service._default_engine_runs
    service._default_engine_runs = cached_runner
    try:
        first = run_unified_analysis(target="ENSP_A", tissue="Lung", use_cache=True)
        same_key = run_unified_analysis(target="ENSP_A", tissue="Lung", use_cache=True)
        second = run_unified_analysis(target="ENSP_B", tissue="Lung", use_cache=True)
    finally:
        service._default_engine_runs = original
        clear_unified_cache()
    first.candidates.loc[:, "finding_class"] = "MUTATED_TEST_VALUE"
    assert "MUTATED_TEST_VALUE" not in set(same_key.candidates["finding_class"])
    assert "MUTATED_TEST_VALUE" not in set(second.candidates["finding_class"])
    assert calls == ["ENSP_A", "ENSP_B"]
