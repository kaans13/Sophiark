"""Disease reference comparison must remain post-hoc to the scientific engine."""

from __future__ import annotations

from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from src.research.disease_bank.benchmark import benchmark_disease_reference
from src.research.disease_bank.models import Disease, DiseaseReferenceSet, MappingSummary
from tests.test_forward_smoke import tiny_inputs


def test_disease_benchmark_leaves_deterministic_forward_output_bitwise_unchanged() -> None:
    from src import biology_logic

    graph, frame, target = tiny_inputs("ENSP")
    old_bonus, old_penalty = biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty
    biology_logic.apply_biological_bonus = lambda *args, **kwargs: None
    biology_logic.apply_functional_penalty = lambda *args, **kwargs: None
    reference = DiseaseReferenceSet(
        Disease("TEST:Disease", "Regression reference"), (), (), MappingSummary(0, 0, 0, 0), "AVAILABLE"
    )
    try:
        with TemporaryDirectory() as directory:
            before = biology_logic.run_infection_simulation(
                graph, frame, Path(directory) / "before", spesifik_hedefler=[target], redistribution_mode="top_n",
                exploratory_top_n=4, null_iterations=0, efficiency_sample_sources=None, local_efficiency_sample_nodes=None,
            )
            benchmark_disease_reference(reference, report=before, perturbed_proteins=(target,), population_size=5)
            after = biology_logic.run_infection_simulation(
                graph, frame, Path(directory) / "after", spesifik_hedefler=[target], redistribution_mode="top_n",
                exploratory_top_n=4, null_iterations=0, efficiency_sample_sources=None, local_efficiency_sample_nodes=None,
            )
        # Timestamp is operational provenance, not a scientific output.  All
        # computed response, loss, PageRank, BC, efficiency, Hinterland and
        # bottleneck fields must be bitwise identical.
        scientific_columns = [column for column in before.columns if column != "Analysis_Timestamp_UTC"]
        pd.testing.assert_frame_equal(before[scientific_columns], after[scientific_columns], check_exact=True)
    finally:
        biology_logic.apply_biological_bonus, biology_logic.apply_functional_penalty = old_bonus, old_penalty
