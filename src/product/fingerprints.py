"""Stable fingerprints for datasets, configurations and cache identities."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any, Mapping


def sha256_file(path: str | Path, *, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def stable_fingerprint(value: Mapping[str, Any]) -> str:
    encoded = json.dumps(
        value, ensure_ascii=False, sort_keys=True, separators=(",", ":"), default=str
    ).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def cache_key(
    *, dataset_fingerprints: Mapping[str, str], preprocessing_versions: Mapping[str, str],
    engine_version: str, engine_config: Mapping[str, Any], target: str,
    tissue: str, analysis_type: str,
) -> str:
    return stable_fingerprint({
        "datasets": dict(dataset_fingerprints),
        "preprocessing": dict(preprocessing_versions),
        "engine_version": engine_version,
        "engine_config": dict(engine_config),
        "target": target,
        "tissue": tissue,
        "analysis_type": analysis_type,
    })
