from __future__ import annotations

import json
import numpy as np
import pandas as pd

from src.evidence.audit import audit_modifier_bias, hill_sensitivity, mid_evidence_audit
from src.evidence.config import EvidenceConfig, EvidenceMode, PhysicalPolicy
from src.evidence.graph import text_bonus
from src.evidence.provenance import attach_provenance
from src.evidence.scores import combine_string_channels
from src.evidence.comparison import compare_classic_evidence


def test_string_combiner_matches_known_independent_probability_math():
    prior = .041
    scores = [.4, .7]
    expected = prior + (1-prior) * (1-((1-(.4-prior)/(1-prior)) * (1-(.7-prior)/(1-prior))))
    assert combine_string_channels(scores, prior) == expected


def test_transferred_text_is_excluded_by_declared_policy():
    config = EvidenceConfig()
    assert "transferred" in config.s_nontext_policy
    assert "textmining" in config.s_nontext_policy


def test_text_bonus_is_bounded_and_near_zero_evidence_cannot_be_rescued():
    config = EvidenceConfig(lambda_text=.10)
    bonus = text_bonus(np.array([0.0, .5, 1.0]), np.array([1.0, 1.0, 1.0]), config)
    assert bonus[0] == 1.0
    assert np.all((bonus >= 1.0) & (bonus <= 1.1))


def test_physical_default_is_context_only_and_multiplier_off():
    config = EvidenceConfig()
    assert config.physical_policy is PhysicalPolicy.CONTEXT_ONLY
    assert config.lambda_physical == 0


def test_provenance_is_machine_readable_on_every_row():
    result = attach_provenance(pd.DataFrame({"gene": ["A", "B"]}), {"engine_version": "x", "tissue": "Lung"})
    assert result["Evidence_Run_Provenance_JSON"].nunique() == 1
    assert json.loads(result.iloc[0]["Evidence_Run_Provenance_JSON"])["tissue"] == "Lung"


def test_bias_gate_demotes_annotation_correlated_modifier_to_context():
    n = 100
    frame = pd.DataFrame({
        "modifier": np.arange(n), "Classic_Degree": np.arange(n),
        "Classic_PageRank": np.arange(n), "Annotation_Availability": np.arange(n),
    })
    decision = audit_modifier_bias(frame, modifier_column="modifier")
    assert not decision.usable_as_multiplier
    assert "context" in decision.reason


def test_mid_evidence_audit_reports_all_bands():
    frame = pd.DataFrame({"s_nontext": np.linspace(.04, 1, 90), "text_bonus": np.linspace(1, 1.05, 90)})
    assert set(mid_evidence_audit(frame)["Evidence_Band"].astype(str)) == {"LOW", "MID", "HIGH"}


def test_classic_evidence_comparison_preserves_missing_values():
    classic = pd.DataFrame({"gene": ["A"], "Delta_PageRank_Pct": [1.0]})
    evidence = pd.DataFrame({"gene": ["B"], "Delta_PageRank_Pct": [2.0]})
    compared = compare_classic_evidence(classic, evidence).set_index("gene")
    assert pd.isna(compared.at["A", "Evidence_Response"])
    assert pd.isna(compared.at["B", "Classic_Response"])


def test_hill_sensitivity_uses_only_predefined_compact_grid():
    result = hill_sensitivity(pd.Series(np.linspace(0, 1, 30)), pd.Series(np.linspace(1, 0, 30)))
    assert len(result) == 3 * 5 * 3
    assert result["max_bonus"].max() <= 1.1
