from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Mapping

import pandas as pd
from functools import lru_cache

from .config import CORUM_PATH, EVIDENCE_DB_PATH, FULL_LINKS_PATH, PHYSICAL_DB_PATH, EvidenceConfig, STRING_TAXON, STRING_VERSION
from .scores import index_metadata


@lru_cache(maxsize=8)
def _identity(path: Path) -> dict[str, object]:
    stat = path.stat()
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return {"path": str(path), "bytes": stat.st_size, "sha256": digest.hexdigest()}


def build_run_provenance(*, config: EvidenceConfig, tissue: str | None,
                         targets: list[str], perturbation: Mapping[str, object]) -> dict[str, object]:
    metadata = index_metadata()
    return {
        "engine": "Evidence · BETA", "engine_version": config.engine_version,
        "organism_taxon": STRING_TAXON, "tissue": tissue or "ALL",
        "targets": list(targets), "string_version": STRING_VERSION,
        "string_data_source": _identity(FULL_LINKS_PATH),
        "evidence_index": _identity(EVIDENCE_DB_PATH),
        "physical_data_source": _identity(PHYSICAL_DB_PATH),
        "corum_data_source": _identity(CORUM_PATH),
        "evidence_mode": config.mode.value, "s_nontext_policy": config.s_nontext_policy,
        "hill_parameters": {
            "gate": config.text_gate, "n": config.hill_n, "m": config.hill_m,
            "k_text": config.hill_k_text, "k_evidence": config.hill_k_evidence,
            "lambda_text": config.lambda_text,
        },
        "physical_policy": config.physical_policy.value,
        "lambda_physical": config.lambda_physical,
        "corum_policy": config.corum_policy.value,
        "string_combiner_validation": metadata,
        "perturbation": dict(perturbation), "configuration_identity": config.identity(),
    }


def attach_provenance(frame: pd.DataFrame, provenance: Mapping[str, object]) -> pd.DataFrame:
    result = frame.copy(deep=True)
    # Machine-readable provenance is intentionally repeated on every exported row.
    result["Evidence_Run_Provenance_JSON"] = json.dumps(dict(provenance), sort_keys=True, separators=(",", ":"))
    result.attrs["evidence_run_provenance"] = dict(provenance)
    return result


def provenance_frame(provenance: Mapping[str, object]) -> pd.DataFrame:
    return pd.DataFrame([{
        key: json.dumps(value, sort_keys=True) if isinstance(value, (dict, list)) else value
        for key, value in provenance.items()
    }])
