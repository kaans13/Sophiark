"""Deterministic, post-hoc comparison of a disease reference set and results."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Mapping, Sequence

import pandas as pd

from .models import DiseaseReferenceSet


_IDENTIFIER_COLUMNS = ("gene", "Protein", "ENSP", "protein")
_LOSS_COLUMNS = ("Network_Loss_Pct", "Network_Loss", "Network_Loss_Score")


def _ordered_ids(frame: pd.DataFrame | None) -> tuple[str, ...]:
    if not isinstance(frame, pd.DataFrame) or frame.empty:
        return ()
    column = next((name for name in _IDENTIFIER_COLUMNS if name in frame.columns), None)
    if column is None:
        return ()
    return tuple(dict.fromkeys(str(item).strip() for item in frame[column].dropna() if str(item).strip()))


def _loss_ids(report: pd.DataFrame | None) -> tuple[str, ...]:
    if not isinstance(report, pd.DataFrame) or report.empty:
        return ()
    loss_column = next((name for name in _LOSS_COLUMNS if name in report.columns), None)
    if loss_column is None:
        return ()
    numeric = pd.to_numeric(report[loss_column], errors="coerce")
    return _ordered_ids(report.loc[numeric.notna() & (numeric != 0)])


@dataclass(frozen=True, slots=True)
class RecallAtK:
    k: int
    returned_count: int
    overlap_count: int
    evaluation_reference_count: int
    recall: float | None
    expected_random_overlap: float | None
    enrichment_ratio: float | None


@dataclass(frozen=True, slots=True)
class DiseaseBenchmark:
    disease_id: str
    disease_name: str
    perturbed_disease_targets: tuple[str, ...]
    evaluation_reference: tuple[str, ...]
    known_response: tuple[str, ...]
    known_network_loss: tuple[str, ...]
    known_redistribution: tuple[str, ...]
    affected_entities: tuple[str, ...]
    network_candidates_outside_reference: tuple[str, ...]
    recalls: tuple[RecallAtK, ...]
    limitation: str = (
        "Disease association is external reference evidence; it does not demonstrate causality, activation, therapeutic relevance, or experimental validation."
    )

    def summary_rows(self) -> list[dict[str, object]]:
        return [{
            "Disease": self.disease_name, "Disease ID": self.disease_id,
            "Perturbed disease-associated targets": len(self.perturbed_disease_targets),
            "Recovered non-perturbed disease-associated genes": len(self.known_response),
            "Disease-associated network loss": len(self.known_network_loss),
            "Disease-associated redistribution": len(self.known_redistribution),
            "Network candidates outside the reference set": len(self.network_candidates_outside_reference),
        }]


def benchmark_disease_reference(
    reference: DiseaseReferenceSet,
    *,
    report: pd.DataFrame,
    redistribution: pd.DataFrame | None = None,
    perturbed_proteins: Iterable[str] = (),
    population_size: int | None = None,
    recall_ks: Sequence[int] = (50, 100, 200, 300),
) -> DiseaseBenchmark:
    """Compare immutable external reference IDs against already computed tables.

    The table order is retained for Recall@K.  Expected overlap is the
    deterministic uniform-reference expectation ``returned * eval / universe``;
    no p-value or significance claim is produced.
    """

    reference_ids = reference.canonical_proteins
    perturbed = frozenset(str(item).strip() for item in perturbed_proteins if str(item).strip())
    evaluation = reference_ids - perturbed
    response_order = _ordered_ids(redistribution if redistribution is not None else report)
    response = frozenset(response_order)
    affected = frozenset(_ordered_ids(report))
    loss = frozenset(_loss_ids(report))
    known_response = tuple(identifier for identifier in response_order if identifier in evaluation)
    known_redistribution = known_response if redistribution is not None else ()
    known_loss = tuple(identifier for identifier in _loss_ids(report) if identifier in evaluation)
    outside = tuple(identifier for identifier in response_order if identifier not in reference_ids)
    universe = int(population_size) if population_size is not None and population_size > 0 else len(affected | response | reference_ids)
    recalls: list[RecallAtK] = []
    for value in recall_ks:
        k = int(value)
        if k <= 0:
            raise ValueError("recall K values must be positive")
        selected = response_order[:k]
        overlap = sum(identifier in evaluation for identifier in selected)
        selected_count = len(selected)
        expected = (selected_count * len(evaluation) / universe) if universe else None
        recalls.append(RecallAtK(
            k=k, returned_count=selected_count, overlap_count=overlap, evaluation_reference_count=len(evaluation),
            recall=(overlap / len(evaluation)) if evaluation else None,
            expected_random_overlap=expected,
            enrichment_ratio=(overlap / expected) if expected and expected > 0 else None,
        ))
    return DiseaseBenchmark(
        disease_id=reference.disease.disease_id, disease_name=reference.disease.disease_name,
        perturbed_disease_targets=tuple(sorted(reference_ids & perturbed)), evaluation_reference=tuple(sorted(evaluation)),
        known_response=known_response, known_network_loss=known_loss, known_redistribution=known_redistribution,
        affected_entities=tuple(sorted(affected & evaluation)), network_candidates_outside_reference=outside,
        recalls=tuple(recalls),
    )


__all__ = ["DiseaseBenchmark", "RecallAtK", "benchmark_disease_reference"]
