"""Extract normalized, read-only biological facts from completed results."""

from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
import pickle
from typing import Iterable, Mapping

import pandas as pd

from .family_detector import detect_family_complex_clusters
from .ontology_normalizer import row_values


_SYMBOL_COLUMNS = ("Symbol", "gene", "Gene", "Gen")
_ANNOTATION_COLUMNS: Mapping[str, tuple[str, ...]] = {
    "go_biological_process": ("GO Biyolojik Süreç", "GO_BP_Terimleri", "go_bp"),
    "go_molecular_function": ("GO Moleküler İşlev", "GO_MF_Terimleri", "go_mf"),
    "go_cellular_component": ("GO Hücresel Bileşen", "GO_CC_Terimleri", "go_cc"),
    "kegg": ("KEGG", "KEGG Yolakları", "kegg"),
    "localization": ("Lokalizasyon", "localization"),
    "protein_name": ("Protein_Adi", "MyGene Adı", "name"),
}
_CACHE_FIELDS: Mapping[str, str] = {
    "go_biological_process": "go_bp", "go_molecular_function": "go_mf",
    "go_cellular_component": "go_cc", "protein_name": "name",
}
_ROOT = Path(__file__).resolve().parents[2]
_HUMAN_MYGENE_CACHE = _ROOT / "data" / "processed" / "mygene_cache.pkl"


@dataclass(frozen=True)
class CandidateFact:
    symbol: str
    delta_pagerank_pct: float | None
    delta_pagerank: float | None
    pagerank_baseline: float | None
    pagerank_perturbed: float | None
    hinterland_baseline: float | None
    bc_baseline: float | None
    annotations: dict[str, tuple[str, ...]]
    annotation_sources: dict[str, str]

    def to_dict(self) -> dict[str, object]:
        value = asdict(self)
        value["annotations"] = {key: list(items) for key, items in self.annotations.items()}
        return value


def _symbol(row: pd.Series) -> str:
    for column in _SYMBOL_COLUMNS:
        if column in row.index and str(row[column]).strip() not in {"", "nan", "None"}:
            return str(row[column]).strip().upper()
    return ""


def _number(row: pd.Series, *columns: str) -> float | None:
    for column in columns:
        if column in row.index:
            value = pd.to_numeric(pd.Series([row[column]]), errors="coerce").iloc[0]
            if pd.notna(value):
                return float(value)
    return None


@lru_cache(maxsize=1)
def _local_human_mygene_cache() -> dict[str, dict]:
    """Read the pre-existing local cache once; analysis never calls MyGene."""
    if not _HUMAN_MYGENE_CACHE.exists():
        return {}
    try:
        with _HUMAN_MYGENE_CACHE.open("rb") as handle:
            raw = pickle.load(handle)
    except (OSError, pickle.UnpicklingError, EOFError):
        return {}
    return raw if isinstance(raw, dict) else {}


def _cached_record(gene: object, *, human: bool) -> dict:
    if not human or not gene:
        return {}
    cache = _local_human_mygene_cache()
    identifier = str(gene)
    return cache.get(identifier) or cache.get(identifier.removeprefix("9606.")) or {}


def extract_candidate_facts(candidates: pd.DataFrame | None, *, organism: str = "Homo sapiens") -> list[CandidateFact]:
    """Copy source facts into an explicit schema without mutating a DataFrame."""
    if candidates is None or candidates.empty:
        return []
    facts: list[CandidateFact] = []
    human = organism.casefold() in {"homo sapiens", "human", "insan"}
    for _, row in candidates.iterrows():
        symbol = _symbol(row)
        if not symbol:
            continue
        annotations = {name: row_values(row, columns) for name, columns in _ANNOTATION_COLUMNS.items()}
        sources = {name: "rapor" if values else "yok" for name, values in annotations.items()}
        cached = _cached_record(row.get("gene"), human=human)
        # The modern workspace receives the scientific report before its legacy
        # display adapter adds GO fields. Reuse the existing local cache only
        # for missing presentation facts, without writing the source frame.
        for name, cache_field in _CACHE_FIELDS.items():
            if not annotations[name] and cached.get(cache_field):
                annotations[name] = row_values(pd.Series({cache_field: cached[cache_field]}), (cache_field,))
                sources[name] = "yerel_mygene_cache" if annotations[name] else "yok"
        facts.append(CandidateFact(
            symbol=symbol,
            delta_pagerank_pct=_number(row, "Delta_PageRank_Pct", "delta_pagerank_pct"),
            delta_pagerank=_number(row, "Delta_PageRank", "delta_pagerank"),
            pagerank_baseline=_number(row, "PageRank_Baseline", "pagerank_baseline"),
            pagerank_perturbed=_number(row, "PageRank_Perturbed", "pagerank_perturbed"),
            hinterland_baseline=_number(row, "Hinterland_Skoru", "hinterland_baseline"),
            bc_baseline=_number(row, "BC_Skoru", "bc_baseline"),
            annotations=annotations,
            annotation_sources=sources,
        ))
    return facts


def facts_as_frame(facts: Iterable[CandidateFact]) -> pd.DataFrame:
    """Create a separate presentation frame, keeping scientific frames intact."""
    return pd.DataFrame([fact.to_dict() for fact in facts])


def extract_biological_facts(
    candidates: pd.DataFrame | None, *, direct_neighbors: Iterable[str] | None = None,
    organism: str = "Homo sapiens",
) -> dict[str, object]:
    """Return source-backed facts and clusters for later deterministic composition."""
    facts = extract_candidate_facts(candidates, organism=organism)
    symbols = [fact.symbol for fact in facts]
    human = organism.casefold() in {"homo sapiens", "human", "insan"}
    clusters = detect_family_complex_clusters(symbols, candidates=candidates, direct_neighbors=direct_neighbors) if human else []
    return {
        "candidates": [fact.to_dict() for fact in facts],
        "clusters": clusters,
        "organism": organism,
        "family_reference_available": human,
        "family_reference_warning": None if human else "Küratörlü HGNC aile referansı insan-geni odaklıdır; fare için aile kümesi üretilmedi.",
        "debug": {
            "input_columns": list(candidates.columns) if candidates is not None else [],
            "fact_count": len(facts),
            "annotation_source_counts": {
                source: sum(1 for fact in facts for value in fact.annotation_sources.values() if value == source)
                for source in ("rapor", "yerel_mygene_cache", "yok")
            },
        },
    }
