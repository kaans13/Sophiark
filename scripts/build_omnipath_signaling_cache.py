"""Build the read-only normalized OmniPath signaling SQLite cache."""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from src.research.entity_resolution import EntityResolver
from src.signaling.ingestion import load_omnipath_signaling_dataset
from src.signaling.storage import write_signaling_cache


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--taxon-id", type=int, choices=(9606, 10090), default=9606)
    parser.add_argument("--project-root", type=Path, default=Path("."))
    arguments = parser.parse_args()
    root = arguments.project_root.resolve()
    resolver = EntityResolver.from_project_data(
        root, include_human=arguments.taxon_id == 9606, include_mouse=arguments.taxon_id == 10090,
    )
    snapshot = root / "data" / "regulatory" / "omnipath" / f"omnipath_{arguments.taxon_id}.tsv"
    destination = root / "data" / "processed" / f"omnipath_signaling_{arguments.taxon_id}.sqlite"
    dataset = load_omnipath_signaling_dataset(snapshot, resolver=resolver, taxon_id=arguments.taxon_id)
    write_signaling_cache(dataset, destination, source_snapshot=snapshot)
    print(json.dumps({
        "cache": str(destination), "fingerprint": dataset.fingerprint,
        "qa": asdict(dataset.qa), "mapping": asdict(dataset.mapping),
    }, indent=2))


if __name__ == "__main__":
    main()
