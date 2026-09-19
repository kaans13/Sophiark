"""Drop-in UI bridge for the separate Evidence calculation engine."""
from __future__ import annotations

from pathlib import Path
import pandas as pd

from src.config import *  # Existing UI constants; Classic source remains untouched.
from .config import EvidenceConfig
from .engine import prepare_evidence_network, run_evidence_simulation, latest_result
from .config import EVIDENCE_DB_PATH, PHYSICAL_DB_PATH, CORUM_PATH


_CONFIG = EvidenceConfig()


def configure(config: EvidenceConfig) -> None:
    global _CONFIG
    _CONFIG = config


def cache_identity() -> tuple:
    """Evidence-only cache namespace; invalidates when configuration or local evidence changes."""
    def stamp(path):
        stat = path.stat()
        return (str(path.resolve()), stat.st_size, stat.st_mtime_ns)
    return (_CONFIG.identity(), stamp(EVIDENCE_DB_PATH), stamp(PHYSICAL_DB_PATH), stamp(CORUM_PATH))


def arayuz_icin_motoru_hazirla(forced_genes=None, hedef_doku=None, bc_sample_sources=None,
                               tissue_normalization_mode="within_tissue"):
    graph, scores, _ = prepare_evidence_network(
        config=_CONFIG, forced_genes=forced_genes, tissue=hedef_doku,
        bc_sample_sources=bc_sample_sources,
    )
    return graph, scores


def run_infection_simulation(g, df: pd.DataFrame, out_dir: Path, hedef_lokalizasyon=None,
                             spesifik_hedefler=None, block_weight_fraction=BLOCK_WEIGHT_FRACTION,
                             damping=PAGERANK_DAMPING, tissue=None, **kwargs):
    if not spesifik_hedefler:
        gates = df[df["Gümrük_Kapisi"] == True] if "Gümrük_Kapisi" in df else df
        spesifik_hedefler = gates.loc[
            gates["Lokalizasyon"].astype(str).str.contains(str(hedef_lokalizasyon), na=False), "gene"
        ].astype(str).tolist()
    if not spesifik_hedefler:
        return pd.DataFrame()
    result = run_evidence_simulation(
        graph=g, scores=df, out_dir=Path(out_dir), targets=list(spesifik_hedefler),
        config=_CONFIG, tissue=tissue or g["evidence_tissue"],
        block_weight_fraction=block_weight_fraction, damping=damping,
        scientific_options=kwargs,
    )
    return result.report
