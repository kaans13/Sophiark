"""UI-facing orchestration adapters around unchanged Sophiark services.

This module intentionally owns no graph algorithm or scientific scoring.  It
is the stable seam between Streamlit pages and the existing motor/services.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import pandas as pd
import streamlit as st

from src.product.runtime_resources import file_identity, product_data_snapshot, symbol_mapping
from src.services.analysis_service import execute_simulation, prepare_network
from src.services.compensation import predicted_compensation_candidates
from src.services.target_stress import preselect_candidates, validate_candidates

HUMAN_TISSUES = ["None", "Adipose Tissue", "Brain", "Breast", "Heart Muscle", "Kidney", "Liver", "Lung", "Pancreas", "Skin", "Spleen", "Testis"]


@dataclass(frozen=True)
class RuntimeContext:
    """Non-persistent UI context; scientific engines remain unchanged."""

    motor_module: object
    output_dir: Path
    is_mouse: bool
    tissue: str


def runtime_for_selection(*, is_mouse: bool, tissue: str) -> RuntimeContext:
    """Resolve the existing human/mouse motor without changing its behaviour."""
    if is_mouse:
        from src.mouse import main as motor_module
        output_dir = Path("outputs_mouse")
    else:
        from src import main as motor_module
        output_dir = Path("outputs")
    return RuntimeContext(motor_module=motor_module, output_dir=output_dir, is_mouse=is_mouse, tissue=tissue)


def tissue_options(is_mouse: bool) -> list[str]:
    if not is_mouse:
        return HUMAN_TISSUES
    try:
        from src.mouse import config
        return ["None"] + sorted({str(value) for value in getattr(config, "MOUSE_TISSUES", [])})
    except Exception:
        return ["None"]


def symbol_options(is_mouse: bool) -> pd.DataFrame:
    """Read existing local symbol mappings; no network or scoring work happens here."""
    if is_mouse:
        from src.mouse import config
        mapping = getattr(config, "MOUSE_SYMBOL_MAP", {})
        return pd.DataFrame({"gene": list(mapping), "Symbol": list(mapping.values())}).dropna().sort_values("Symbol")
    path = Path("data") / "processed" / "ensp_with_symbols.csv"
    mapping = symbol_mapping(path)
    return pd.DataFrame({"gene": mapping.keys(), "Symbol": mapping.values()}).sort_values("Symbol")


def regulators_for_selection(is_mouse: bool) -> dict:
    """Expose the same species-specific TRRUST mapping used by the legacy UI."""
    if is_mouse:
        from src.mouse import config
        return getattr(config, "TARGET_TO_REGULATORS", {})
    from src import config
    return getattr(config, "TARGET_TO_REGULATORS", {})


def directed_evidence_for_selection(is_mouse: bool) -> tuple[pd.DataFrame, dict]:
    """TRRUST ve OmniPath'ı PPI grafına karıştırmadan ayrı overlay döndürür."""
    from src.scientific.regulatory import load_omnipath_overlay, trrust_records
    if is_mouse:
        from src.mouse import config
    else:
        from src import config
    trrust = trrust_records(getattr(config, "TARGET_TO_REGULATORS", {}))
    omnipath, metadata = load_omnipath_overlay(getattr(config, "OMNIPATH_SNAPSHOT_PATH", None))
    metadata["trrust_records"] = len(trrust)
    metadata["structural_network_evidence"] = "STRING undirected PPI"
    metadata["simulation_evidence"] = "Sophiark edge attenuation + measured network change"
    return pd.concat([trrust, omnipath], ignore_index=True), metadata


def prepare_active_network(context: RuntimeContext, *, forced_genes: tuple[str, ...] | None,
                           bc_sample_sources: int | None):
    if context.is_mouse:
        from src.mouse import config
        cache_identity = (
            "mouse",
            file_identity(getattr(config, "DB_PATH", "missing-mouse-db")),
            file_identity(getattr(config, "OUT_CSV", "missing-mouse-output")),
        )
    else:
        cache_identity = product_data_snapshot().identity
    return prepare_network(
        context.motor_module, forced_genes, context.tissue, bc_sample_sources,
        st.session_state.get("tissue_normalization_mode", "within_tissue"),
        cache_identity=cache_identity,
    )


def run_forward_simulation(*, context: RuntimeContext, graph, scores: pd.DataFrame, output_dir: Path,
                           targets: list[str] | None, localization: str, block_strength: float,
                           damping: float, apply_druggability: Callable, ghost_filter: Callable):
    output_dir.mkdir(parents=True, exist_ok=True)
    efficiency_full = st.session_state.get("efficiency_validation_mode", "Fast Approximate") == "Full Validation"
    return execute_simulation(
        motor_module=context.motor_module, graph=graph, scores=scores, output_dir=output_dir,
        targets=targets, localization=localization, block_strength=block_strength, damping=damping,
        apply_druggability=apply_druggability, ghost_filter=ghost_filter,
        scientific_options={
            "redistribution_mode": st.session_state.get("redistribution_mode", "null_fdr"),
            "exploratory_top_n": int(st.session_state.get("redistribution_top_n", 250)),
            "redistribution_percentile": float(st.session_state.get("redistribution_percentile", 99.0)),
            "null_iterations": int(st.session_state.get("null_iterations", 99)),
            "fdr_alpha": float(st.session_state.get("fdr_alpha", 0.05)),
            "random_seed": int(st.session_state.get("scientific_random_seed", 42)),
            "efficiency_sample_sources": None if efficiency_full else int(st.session_state.get("efficiency_sample_sources", 256)),
            "local_efficiency_sample_nodes": None if efficiency_full else int(st.session_state.get("local_efficiency_sample_nodes", 128)),
            "tissue": context.tissue,
            "tissue_normalization_mode": st.session_state.get("tissue_normalization_mode", "within_tissue"),
            "bc_mode": "fast_approximate" if st.session_state.get("bc_sample_sources", 1200) else "full_validation",
            "bc_sample_size": st.session_state.get("bc_sample_sources", 1200),
            "edge_evidence_weighting_mode": st.session_state.get("edge_evidence_weighting_mode", "structural_only"),
        },
        signed_result_sink=lambda frame: st.session_state.__setitem__(
            "signed_redistribution_result", frame
        ),
    )


def run_target_stress_search(*, context: RuntimeContext, graph, scores: pd.DataFrame, target_gene: str,
                             regulators: dict, gene_to_symbol: dict[str, str], mode: str, limit: int,
                             block_strength: float, damping: float, ghost_filter: Callable,
                             deep_validation_top_n: int = 0) -> pd.DataFrame:
    candidates = preselect_candidates(
        graph=graph, scores=scores, target_gene=target_gene, regulators=regulators,
        gene_to_symbol=gene_to_symbol, mode=mode, limit=limit,
    )
    if candidates.empty:
        return pd.DataFrame()
    return validate_candidates(
        candidates=candidates, motor_module=context.motor_module, graph=graph, scores=scores,
        target_gene=target_gene, output_root=context.output_dir, block_strength=block_strength,
        damping=damping, ghost_filter=ghost_filter,
        deep_validation_top_n=deep_validation_top_n,
    )


def run_compensation_analysis(*, context: RuntimeContext, graph, scores: pd.DataFrame, perturbed_gene: str,
                              block_strength: float, damping: float, ghost_filter: Callable) -> pd.DataFrame:
    symbols = symbol_options(context.is_mouse)
    gene_to_symbol = dict(zip(symbols.get("gene", []), symbols.get("Symbol", [])))
    return predicted_compensation_candidates(
        motor_module=context.motor_module, graph=graph, scores=scores, perturbed_gene=perturbed_gene,
        output_root=context.output_dir, block_strength=block_strength, damping=damping,
        ghost_filter=ghost_filter, gene_to_symbol=gene_to_symbol,
    )
