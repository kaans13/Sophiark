"""Application, engine-contract, data-processing and export versions."""

SOPHIARK_VERSION = "5.3.0"

# These are contract/provenance identifiers, not claims that engine code changed.
CLASSIC_ENGINE_VERSION = "classic-reference-v5.3"
DIRECTED_ENGINE_VERSION = "directed-beta-reference-v1"
EVIDENCE_ENGINE_VERSION = "evidence-experimental-reference-v1"
UNIFIED_SCHEMA_VERSION = "1.0"
EXPORT_SCHEMA_VERSION = "sophiark-export-v1"
BUNDLE_SCHEMA_VERSION = 1

PREPROCESSING_VERSIONS = {
    "string": "string-normalization-v1",
    "hpa": "hpa-tpm-to-ensp-v1",
    "omnipath": "omnipath-ingestion-v1",
    "corum": "corum-membership-v1",
    "trrust": "trrust-legacy-compatible-v1",
    "canonical_mapping": "ensp-mapping-v1",
}

ENGINE_VERSIONS = {
    "classic": CLASSIC_ENGINE_VERSION,
    "directed": DIRECTED_ENGINE_VERSION,
    "evidence": EVIDENCE_ENGINE_VERSION,
    "unified_schema": UNIFIED_SCHEMA_VERSION,
}
