"""UI katmanının bilimsel motorla paylaştığı hafif veri sözleşmeleri."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass(frozen=True)
class AnalysisContext:
    """Kalıcı snapshot değildir; aktif Streamlit çalıştırmasının bağlamıdır."""

    species: str
    tissue: str
    output_dir: Path
    report_path: Path
    protein_prefix: str


@dataclass
class PreparedAnalysis:
    """Motorun ürettiği mutable graf ve ona bağlı skor tablosu."""

    graph: Any
    scores: pd.DataFrame


@dataclass(frozen=True)
class PerturbationComparison:
    """Raw baseline-versus-perturbation values; not a biological activity score.

    ``local_stress_pct`` geriye dönük alan adıdır ve gerçek local efficiency
    değildir; Top_Positive_PageRank_Mean_Pct değerini taşır.
    """

    intervention_genes: tuple[str, ...]
    observed_genes: tuple[str, ...]
    system_shift_pct: float
    local_stress_pct: float
    blocked_edges: int
    observed_metrics: dict[str, dict[str, float]]
    top_pagerank_gains: tuple[dict[str, float | str], ...]
    report: pd.DataFrame


@dataclass
class PropagationTrace:
    target_gene: str
    nodes: pd.DataFrame
    routes: pd.DataFrame
    edges: pd.DataFrame = field(default_factory=pd.DataFrame)
    evidence: pd.DataFrame = field(default_factory=pd.DataFrame)
    metadata: dict[str, Any] = field(default_factory=dict)


@dataclass
class TissueDifferentialResult:
    summary: pd.DataFrame
    critical_genes: pd.DataFrame
    notices: list[str]
    similarity: pd.DataFrame


@dataclass(frozen=True)
class StructuralEfficiencyMetrics:
    global_baseline: float
    global_perturbed: float
    global_change: float
    global_change_pct: float
    mean_local_baseline: float
    mean_local_perturbed: float
    mean_local_change_pct: float
    target_local_baseline: float
    target_local_perturbed: float
    target_local_change_pct: float


@dataclass(frozen=True)
class ScientificRunMetadata:
    values: dict[str, Any]


@dataclass
class PredictedNetworkInteractionResult:
    targets: tuple[str, ...]
    metric_comparison: pd.DataFrame
    relative_tolerance: float
