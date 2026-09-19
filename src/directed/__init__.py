"""Optional hybrid directed calculation engine; classic remains authoritative by default."""

from .comparison import compare_classic_and_directed
from .engine import (
    compute_directed_betweenness, compute_directed_pagerank, run_directed_calculation,
    safe_run_directed_calculation, suppress_target_edges,
)
from .effective_graph import (
    EffectiveGraphSnapshot, effective_graph_fingerprint, graph_weight_parity,
    prepare_effective_directed_source_graph,
)
from .hybrid_graph import build_hybrid_directed_graph, directed_graph_cache_key
from .models import (
    DirectedCalculationResult, DirectedEngineStatus, EngineComparisonResult,
    HybridDirectedGraph, HybridGraphProvenance, PerturbationStrategy,
)

__all__ = [
    "DirectedCalculationResult", "DirectedEngineStatus", "EngineComparisonResult",
    "HybridDirectedGraph", "HybridGraphProvenance", "PerturbationStrategy",
    "build_hybrid_directed_graph", "compare_classic_and_directed",
    "compute_directed_pagerank", "directed_graph_cache_key",
    "compute_directed_betweenness",
    "EffectiveGraphSnapshot", "effective_graph_fingerprint", "graph_weight_parity",
    "prepare_effective_directed_source_graph",
    "run_directed_calculation", "safe_run_directed_calculation", "suppress_target_edges",
]
