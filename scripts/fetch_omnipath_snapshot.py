"""OmniPath interaction endpointinden versioned, local snapshot indirir.

CanlÄ± API analiz sÄ±rasÄ±nda asla Ã§aÄŸrÄ±lmaz; bu komut yalnÄ±zca veri yenileme
iÃ§indir. Snapshot ve metadata birlikte saklanÄ±r.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd
import requests


ENDPOINT = "https://omnipathdb.org/interactions"
# Base interaction columns (is_directed/stimulation/inhibition) zaten endpoint
# tarafÄ±ndan gelir; `fields` yalnÄ±zca ek kanÄ±t alanlarÄ± iÃ§in kullanÄ±lÄ±r.
FIELDS = ["sources", "references", "curation_effort"]


def fetch_snapshot(organism_taxid: int, destination: Path) -> dict:
    params = {
        "format": "tsv", "genesymbols": "1", "organisms": str(organism_taxid),
        "fields": ",".join(FIELDS),
    }
    response = requests.get(ENDPOINT, params=params, timeout=180)
    response.raise_for_status()
    destination.mkdir(parents=True, exist_ok=True)
    snapshot = destination / f"omnipath_{organism_taxid}.tsv"
    snapshot.write_bytes(response.content)
    frame = pd.read_csv(snapshot, sep="\t")
    # Endpoint sÃ¼rÃ¼mÃ¼ ile explicit release her zaman yayÄ±nlanmadÄ±ÄŸÄ±ndan,
    # indirilen byte checksum'Ä± reproducibility kimliÄŸidir.
    metadata = {
        "source": ENDPOINT,
        "downloaded_at": datetime.now(timezone.utc).isoformat(),
        "organism": "Homo sapiens" if organism_taxid == 9606 else str(organism_taxid),
        "organism_taxid": organism_taxid,
        "release": response.headers.get("X-OmniPath-Version") or response.headers.get("ETag") or "endpoint_release_not_provided",
        "row_count": int(len(frame)),
        "checksum_sha256": hashlib.sha256(response.content).hexdigest(),
        "parser_version": "omnipath-tsv-v1",
        "request": {"endpoint": ENDPOINT, "params": params},
    }
    (destination / f"omnipath_{organism_taxid}_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return {"snapshot": str(snapshot), **metadata}


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--organism-taxid", type=int, default=9606)
    parser.add_argument("--destination", type=Path, default=Path("data/regulatory/omnipath"))
    print(json.dumps(fetch_snapshot(parser.parse_args().organism_taxid, parser.parse_args().destination), indent=2))
