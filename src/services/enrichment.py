"""Tür ve analiz evreni bilinçli, yerel ORA enrichment servisi."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import time

import numpy as np
import pandas as pd
from scipy.stats import hypergeom

from src.scientific.statistics import benjamini_hochberg


@dataclass(frozen=True)
class EnrichmentLibrary:
    name: str
    family: str
    species: str


LIBRARIES = {
    "human": (
        EnrichmentLibrary("KEGG_2021_Human", "KEGG", "Homo sapiens"),
        EnrichmentLibrary("GO_Biological_Process_2021", "GO:BP", "Homo sapiens"),
    ),
    "mouse": (
        EnrichmentLibrary("KEGG_2019_Mouse", "KEGG", "Mus musculus"),
        EnrichmentLibrary("GO_Biological_Process_2021", "GO:BP", "Mus musculus"),
    ),
}


@lru_cache(maxsize=16)
def _get_gene_sets(library_name: str, organism: str) -> dict[str, list[str]]:
    import gseapy as gp
    last_error: Exception | None = None
    for attempt in range(2):
        try:
            return gp.get_library(name=library_name, organism=organism)
        except Exception as error:
            last_error = error
            if attempt == 0:
                time.sleep(0.35)
    assert last_error is not None
    raise last_error


def _clean(values: list[str]) -> set[str]:
    return {
        str(value).strip()
        for value in values
        if str(value).strip() not in {"", "—", "nan", "None"}
    }


def local_over_representation(
    input_genes: list[str],
    background_genes: list[str],
    gene_sets: dict[str, list[str]],
) -> pd.DataFrame:
    """Aktif ağ evreninde hypergeometric ORA + BH-FDR."""
    background = _clean(background_genes)
    selected = _clean(input_genes) & background
    if not selected or not background:
        return pd.DataFrame()
    population, draws = len(background), len(selected)
    rows: list[dict[str, object]] = []
    for term, members in gene_sets.items():
        pathway = _clean(list(members)) & background
        overlap = selected & pathway
        if not overlap:
            continue
        p_value = float(hypergeom.sf(len(overlap) - 1, population, len(pathway), draws))
        rows.append({
            "Term": str(term),
            "P-value": p_value,
            "Overlap": f"{len(overlap)}/{len(pathway)}",
            "Genes": ";".join(sorted(overlap)),
            "Matched_Genes": len(overlap),
            "Input_Genes": draws,
            "Background_Universe_Size": population,
            "Gene_Set_Size_In_Universe": len(pathway),
        })
    result = pd.DataFrame(rows)
    if result.empty:
        return result
    result["Adjusted P-value"] = benjamini_hochberg(result["P-value"].to_numpy())
    return result.sort_values(["Adjusted P-value", "P-value"], kind="mergesort").reset_index(drop=True)


def run_post_simulation_enrichment(
    symbols: list[str],
    *,
    is_mouse: bool,
    background_symbols: list[str] | None = None,
) -> tuple[pd.DataFrame, list[str]]:
    """KEGG ve GO-BP için tür-spesifik local ORA çalıştır.

    Enrichr'ın remote istatistiği kullanılmaz; gseapy yalnızca güncel gene-set
    üyeliklerini indirir. İstatistik aktif doku grafının custom background'u ile
    yerelde hesaplanır. Mouse kütüphanesi erişilemiyorsa insan orthologuna sessiz
    dönüşüm yapılmaz.
    """
    unique_symbols = sorted(_clean(symbols))
    universe = sorted(_clean(background_symbols or []))
    if is_mouse:
        # Mouse Enrichr gene-set üyeleri büyük harf sembol namespace'i kullanır.
        # Bu yalnızca eşleştirme canonicalizasyonudur; human ortholog dönüşümü değildir.
        unique_symbols = sorted({symbol.upper() for symbol in unique_symbols})
        universe = sorted({symbol.upper() for symbol in universe})
    if not unique_symbols:
        return pd.DataFrame(), ["Etkilenen genler için doğrulanmış sembol bulunamadı."]
    if not universe:
        return pd.DataFrame(), ["Aktif doku ağı sembol evreni bulunamadı; varsayılan genom evrenine sessizce geçilmedi."]

    key = "mouse" if is_mouse else "human"
    organism = "Mouse" if is_mouse else "Human"
    frames: list[pd.DataFrame] = []
    notices: list[str] = []
    for library in LIBRARIES[key]:
        try:
            gene_sets = _get_gene_sets(library.name, organism)
            result = local_over_representation(unique_symbols, universe, gene_sets)
            if result.empty:
                notices.append(f"{library.name}: aktif ağ evreninde eşleşen terim bulunamadı.")
                continue
            result["Kaynak"] = library.name
            result["Pathway_Family"] = library.family
            result["Species"] = library.species
            result["Source_Species"] = library.species
            result["Target_Species"] = library.species
            result["Symbol_Namespace"] = "HGNC symbol" if not is_mouse else "MGI symbol (uppercase Enrichr canonicalization)"
            result["Ortholog_Conversion_Applied"] = False
            result["Ortholog_Conversion_Rate"] = 1.0
            # FDR değeri ham sonuç olarak korunur. Sunum katmanı bu veriyi
            # "anlamlı/anlamsız" biçiminde yorumlamaz.
            frames.append(result.head(20))
        except Exception as error:
            notices.append(
                f"{library.name}: tür-spesifik gene set alınamadı ({type(error).__name__}); "
                "ortholog dönüşümü uygulanmadı. Eksik sonuç negatif evidence değildir."
            )
    return (pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()), notices
