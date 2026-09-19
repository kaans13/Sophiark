from __future__ import annotations

import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import requests


BASE_DIR = Path(__file__).resolve().parents[1]
OUT_DIR = BASE_DIR / "data" / "orthology"
TSV_PATH = OUT_DIR / "ensembl_human_mouse_orthology.tsv"
METADATA_PATH = OUT_DIR / "ensembl_human_mouse_orthology.metadata.json"
# Pin the archive reached by the public Ensembl endpoint at snapshot time.
MART_URL = "https://jun2026.archive.ensembl.org/biomart/martservice"

QUERY = """<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE Query>
<Query virtualSchemaName="default" formatter="TSV" header="1" uniqueRows="1"
       count="" datasetConfigVersion="0.6">
  <Dataset name="hsapiens_gene_ensembl" interface="default">
    <Attribute name="ensembl_gene_id"/>
    <Attribute name="external_gene_name"/>
    <Attribute name="mmusculus_homolog_ensembl_gene"/>
    <Attribute name="mmusculus_homolog_associated_gene_name"/>
    <Attribute name="mmusculus_homolog_orthology_type"/>
    <Attribute name="mmusculus_homolog_orthology_confidence"/>
  </Dataset>
</Query>"""


def main() -> None:
    response = requests.get(MART_URL, params={"query": QUERY}, timeout=180)
    response.raise_for_status()
    content = response.content
    if not content.startswith(b"Gene stable ID\t"):
        raise RuntimeError("Unexpected Ensembl BioMart response; snapshot not written")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    TSV_PATH.write_bytes(content)
    lines = content.decode("utf-8").splitlines()
    mapped = sum(1 for line in lines[1:] if len(line.split("\t")) >= 4 and line.split("\t")[2])
    metadata = {
        "source": "Ensembl Compara via Ensembl BioMart",
        "requested_url": MART_URL,
        "resolved_url": response.url.split("?", 1)[0],
        "retrieved_at_utc": datetime.now(timezone.utc).isoformat(),
        "human_taxon_id": 9606,
        "mouse_taxon_id": 10090,
        "sha256": hashlib.sha256(content).hexdigest(),
        "row_count": max(0, len(lines) - 1),
        "mapped_row_count": mapped,
        "query": QUERY,
    }
    METADATA_PATH.write_text(
        json.dumps(metadata, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    print(json.dumps(metadata, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
