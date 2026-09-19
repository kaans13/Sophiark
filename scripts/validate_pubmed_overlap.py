"""Reproducible PubMed Literature Concordance validator.

`--pmids-json` ile Ã§alÄ±ÅŸtÄ±rÄ±ldÄ±ÄŸÄ±nda internet gerekmez; aynÄ± PMID snapshot
aynÄ± sonuÃ§u Ã¼retir. `--fetch-live` yalnÄ±zca snapshot oluÅŸturma amacÄ±yla vardÄ±r.
"""
from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests

PROJECT_ROOT = Path(__file__).resolve().parents[1]
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.scientific.metadata import SOPHIARK_VERSION
from src.scientific.pubmed_validation import calculate_concordance, extract_symbols_from_texts, write_validation_snapshot


EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"


def fetch_pmids(query: str, retmax: int) -> list[str]:
    response = requests.get(f"{EUTILS}/esearch.fcgi", params={"db": "pubmed", "term": query, "retmax": retmax, "retmode": "json"}, timeout=60)
    response.raise_for_status()
    return response.json()["esearchresult"].get("idlist", [])


def fetch_abstracts(pmids: list[str]) -> list[str]:
    if not pmids:
        return []
    response = requests.get(f"{EUTILS}/efetch.fcgi", params={"db": "pubmed", "id": ",".join(pmids), "rettype": "abstract", "retmode": "text"}, timeout=90)
    response.raise_for_status()
    return [response.text]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--predictions", type=Path, required=True)
    parser.add_argument("--mapping", type=Path, default=Path("data/processed/ensp_with_symbols.csv"))
    parser.add_argument("--target", required=True)
    parser.add_argument("--organism", default="Homo sapiens")
    parser.add_argument("--tissue", default="None")
    parser.add_argument("--query", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--pmids-json", type=Path)
    parser.add_argument("--fetch-live", action="store_true")
    parser.add_argument("--retmax", type=int, default=50)
    args = parser.parse_args()

    predictions = pd.read_csv(args.predictions)
    mapping = pd.read_csv(args.mapping, usecols=["gene", "Symbol"]).dropna()
    merged = predictions.merge(mapping, on="gene", how="left")
    predicted_symbols = set(merged["Symbol"].dropna().astype(str).str.upper())
    universe = set(mapping["Symbol"].astype(str).str.upper())
    if args.fetch_live:
        pmids = fetch_pmids(args.query, args.retmax)
        texts = fetch_abstracts(pmids)
    elif args.pmids_json:
        pmids = json.loads(args.pmids_json.read_text(encoding="utf-8"))
        texts_path = args.pmids_json.with_name("abstracts.txt")
        if not texts_path.exists():
            raise SystemExit("Offline validation iÃ§in pmids.json ile aynÄ± dizinde abstracts.txt gerekli.")
        texts = [texts_path.read_text(encoding="utf-8")]
    else:
        raise SystemExit("--fetch-live veya --pmids-json zorunludur.")
    literature, excluded = extract_symbols_from_texts(texts, universe)
    concordance = calculate_concordance(predicted_symbols, literature, universe)
    concordance["extraction_method"] = "exact-token match against local symbol universe; symbols <3 chars excluded"
    concordance["excluded_ambiguous_short_symbols"] = excluded
    query_metadata = {
        "target_gene": args.target, "organism": args.organism, "tissue": args.tissue,
        "perturbation_parameters": "read from prediction CSV metadata columns when present",
        "sophiark_version": SOPHIARK_VERSION, "selected_predicted_gene_set": sorted(predicted_symbols),
        "gene_selection_method": "provided prediction CSV", "pubmed_query": args.query,
        "query_date_utc": datetime.now(timezone.utc).isoformat(), "leakage_check": "PubMed set is read only after predictions; never passed to ranking.",
    }
    mapping_report = merged[["gene", "Symbol"]].rename(columns={"gene": "network_node_id", "Symbol": "resolved_symbol"})
    root = write_validation_snapshot(args.output, query=query_metadata, pmids=pmids, literature_symbols=literature, predictions=merged, concordance=concordance, mapping=mapping_report)
    if args.fetch_live:
        (root / "abstracts.txt").write_text("\n".join(texts), encoding="utf-8")
    print(json.dumps(concordance, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
