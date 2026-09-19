"""Bilimsel motoru değiştirmeden UI için analiz orkestrasyonu sağlar."""

from __future__ import annotations

import inspect
from typing import Callable, Optional

import pandas as pd

from src.models import PreparedAnalysis
from src.product.runtime_resources import DetachedPreparedContextCache


_PREPARED_NETWORKS = DetachedPreparedContextCache(maxsize=3)


def prepare_network(motor_module, forced_genes: tuple | None, tissue: str, bc_sample_sources: int | None,
                    tissue_normalization_mode: str = "within_tissue",
                    cache_identity: object | None = None) -> PreparedAnalysis:
    identity = cache_identity if isinstance(cache_identity, tuple) else (repr(cache_identity),)
    graph, scores = _PREPARED_NETWORKS.get_or_build(
        (
            str(getattr(motor_module, "__name__", motor_module)), forced_genes,
            tissue, bc_sample_sources, tissue_normalization_mode, identity,
        ),
        lambda: motor_module.arayuz_icin_motoru_hazirla(
            forced_genes=list(forced_genes) if forced_genes else None,
            hedef_doku=None if tissue == "None" else tissue,
            bc_sample_sources=bc_sample_sources,
            tissue_normalization_mode=tissue_normalization_mode,
        ),
    )
    return PreparedAnalysis(graph=graph, scores=scores)


def execute_simulation(*, motor_module, graph, scores: pd.DataFrame, output_dir, targets: Optional[list], localization: str,
                       block_strength: float, damping: float, apply_druggability: Callable[[pd.DataFrame], pd.DataFrame],
                       ghost_filter: Callable[[pd.DataFrame, object], tuple[pd.DataFrame, list]],
                       scientific_options: dict | None = None,
                       signed_result_sink: Callable[[pd.DataFrame], None] | None = None) -> tuple[pd.DataFrame, list]:
    # prepare_network cache'i mutable bir igraph nesnesi döndürür. Simülasyon
    # motoru kenar ağırlıklarını geçici değiştirdiğinden, cache'teki grafı ve
    # aynı oturumdaki ağ görselleştirmesini korumak için çalışma kopyası kullan.
    working_graph = graph.copy()
    kwargs = {"spesifik_hedefler": targets} if targets else {"hedef_lokalizasyon": localization}
    simulation = motor_module.run_infection_simulation
    supported = inspect.signature(simulation).parameters
    optional_kwargs = {
        "block_weight_fraction": block_strength,
        "damping": damping,
        **(scientific_options or {}),
    }
    compatible_kwargs = {
        key: value for key, value in optional_kwargs.items() if key in supported
    }
    result = simulation(
        working_graph, scores, output_dir, **compatible_kwargs, **kwargs
    )
    if result is None or result.empty:
        return pd.DataFrame(), []
    # Signed sonuç Pandas attrs yerine motor modül belleğinden alınır: attrs
    # içinde DataFrame taşımak pandas concat/nlargest işlemlerini bozabilir.
    module_name = str(getattr(motor_module, "__name__", ""))
    if module_name.startswith("src.mouse"):
        from src.mouse.biology_logic import latest_signed_redistribution
    elif module_name.startswith("src.evidence"):
        from src.evidence.engine import latest_result
        evidence_result = latest_result()
        signed_redistribution = (
            evidence_result.signed_response if evidence_result is not None else None
        )
    else:
        from src.biology_logic import latest_signed_redistribution
    if not module_name.startswith("src.evidence"):
        signed_redistribution = latest_signed_redistribution()
    result = apply_druggability(result)
    filtered, ghosts = ghost_filter(result, working_graph)
    if isinstance(signed_redistribution, pd.DataFrame) and signed_result_sink is not None:
        signed_result_sink(signed_redistribution)
    return filtered, ghosts
