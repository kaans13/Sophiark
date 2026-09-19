"""Presentation data adapters: symbol-first reports, expression and comparisons."""

from __future__ import annotations

import pickle
import sqlite3
from pathlib import Path

import pandas as pd


def symbol_map(symbols: pd.DataFrame) -> dict[str, str]:
    if symbols is None or symbols.empty:
        return {}
    return dict(zip(symbols["gene"].astype(str), symbols["Symbol"].astype(str)))


def symbol_first_report(report: pd.DataFrame, mapping: dict[str, str]) -> pd.DataFrame:
    result = report.copy()
    ids = result.get("gene", pd.Series("", index=result.index)).astype(str)
    resolved = ids.map(mapping)
    if "Symbol" in result:
        resolved = result["Symbol"].where(result["Symbol"].notna() & result["Symbol"].astype(str).ne(""), resolved)
    result.insert(0, "Gen Sembolü", resolved.fillna("Alias bulunamadı"))
    if "gene" in result:
        result = result.rename(columns={"gene": "Protein Kimliği"})
    return result


def critical_rows(report: pd.DataFrame) -> pd.DataFrame:
    if report is None or report.empty:
        return pd.DataFrame()
    if "Strong_Redistribution" in report:
        selected = report[report["Strong_Redistribution"].fillna(False).astype(bool)]
    elif "Hasar_Tipi" in report:
        selected = report[report["Hasar_Tipi"].astype(str).str.contains("Üçüncül|stres", case=False, regex=True, na=False)]
    else:
        selected = report
    order = "Abs_Delta_PageRank" if "Abs_Delta_PageRank" in selected else "Hinterland_Skoru"
    return selected.sort_values(order, ascending=False) if order in selected else selected


def query_expression(gene_id: str, *, is_mouse: bool) -> pd.DataFrame:
    if is_mouse:
        from src.mouse import config
        table = "mouse_tissue_expression"
    else:
        from src import config
        table = "tissue_expression"
    clean = str(gene_id).removeprefix("9606.").removeprefix("10090.")
    with sqlite3.connect(str(config.DB_PATH)) as connection:
        return pd.read_sql_query(
            f'SELECT tissue AS "Doku", MAX(expression_level) AS "İfade Seviyesi" FROM {table} '
            'WHERE protein_id = ? AND expression_level > 0 GROUP BY tissue ORDER BY "İfade Seviyesi" DESC',
            connection, params=(clean,),
        )


def cross_species_comparison(base_dir: Path = Path(".")) -> pd.DataFrame:
    human_path = base_dir / "outputs" / "reports" / "enfeksiyon_sok_dalgasi_raporu.csv"
    mouse_path = base_dir / "outputs_mouse" / "reports" / "mouse_enfeksiyon_sok_dalgasi_raporu.csv"
    ortholog_path = base_dir / "data" / "processed" / "mouse_ortholog_map.pkl"
    symbols_path = base_dir / "data" / "processed" / "ensp_with_symbols.csv"
    if not all(path.exists() for path in (human_path, mouse_path, ortholog_path, symbols_path)):
        return pd.DataFrame()
    human, mouse = pd.read_csv(human_path), pd.read_csv(mouse_path)
    symbols = pd.read_csv(symbols_path)
    with ortholog_path.open("rb") as handle:
        orthologs = pickle.load(handle)
    from src.mouse.config import MOUSE_SYMBOL_MAP
    human_map = dict(zip(symbols["gene"].astype(str), symbols["Symbol"].astype(str)))
    human = human.copy(); mouse = mouse.copy()
    human["İnsan ENSP"] = human["gene"].astype(str)
    mouse["Fare Protein Kimliği"] = mouse["gene"].astype(str)
    human["İnsan Sembolü"] = human["gene"].astype(str).map(human_map).fillna(human.get("Symbol"))
    human["Fare Sembolü"] = human["İnsan Sembolü"].map(orthologs)
    mouse["Fare Sembolü"] = mouse.get("Symbol", mouse["gene"].astype(str).map(MOUSE_SYMBOL_MAP))
    hcols = [c for c in ["İnsan ENSP", "İnsan Sembolü", "Fare Sembolü", "Hasar_Tipi", "Hinterland_Skoru", "BC_Skoru"] if c in human]
    mcols = [c for c in ["Fare Sembolü", "Fare Protein Kimliği", "Hasar_Tipi", "Hinterland_Skoru", "BC_Skoru"] if c in mouse]
    h = human[hcols].rename(columns={"Hasar_Tipi":"İnsan Etkisi","Hinterland_Skoru":"İnsan Hinterland","BC_Skoru":"İnsan BC"})
    m = mouse[mcols].rename(columns={"Hasar_Tipi":"Fare Etkisi","Hinterland_Skoru":"Fare Hinterland","BC_Skoru":"Fare BC"})
    result = h[h["Fare Sembolü"].notna()].merge(m, on="Fare Sembolü", how="outer", indicator=True)
    merge_state = result.pop("_merge")
    result["Karşılaştırma Durumu"] = merge_state.map({"both":"Her iki türde","left_only":"Yalnızca insan","right_only":"Yalnızca fare"})
    result["İlişki"] = merge_state.map({
        "both": "İnsan–fare ortolog eşleşmesi",
        "left_only": "İnsan kaydı; eşleşen fare sonuç satırı yok",
        "right_only": "Fare kaydı; insan eşleşmesi yok",
    })
    if {"İnsan Hinterland", "Fare Hinterland"} <= set(result):
        result["Hinterland Farkı"] = pd.to_numeric(result["İnsan Hinterland"], errors="coerce") - pd.to_numeric(result["Fare Hinterland"], errors="coerce")
        result = result.sort_values("Hinterland Farkı", key=lambda x: x.abs(), ascending=False)
    return result
