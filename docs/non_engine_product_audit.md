# Non-engine product audit

## Final status

**NON_ENGINE_AUDIT_PARTIAL.** The audited product contracts around real Unified exports, bundle/history, registry health, capabilities, canonical identity and cross-species presentation pass. The audit did not complete live browser/session interaction, clean-environment installation, exhaustive visual rendering of every table family, or production builder lifecycle execution.

## Evidence-backed findings

- Engine freeze is 12/12 exact; no protected engine source changed in this audit.
- Product registry health and the runtime snapshot agree: all seven current capabilities are `READY`; start/end dataset fingerprints are identical.
- The real CFTR Unified fixture contains 100 unique `entity_id` rows, 16,954 unique Classic `gene` rows and 16,954 unique Directed `entity` rows. Required `Gen`/`ENSP` candidate fields are present.
- Reopened real `.sophiark` bundles preserve `COMPLETE` status and 100 candidates; history rebuild indexed two valid bundles with no corruption or duplicates. Bundle member scans found no personal absolute paths.
- Cross-species presentation exposes `İnsan ENSP`, `Fare Protein Kimliği` and `İlişki`. Current fixture rows are human-only because the existing Mouse report has no rows; this is presented explicitly.
- Important non-engine imports produced no output-file writes during the controlled import snapshot.
- The multipage pages delegate through service adapters. The legacy `app.py` remains a large active/legacy-referenced Streamlit controller with direct orchestration and session state, so UI/service separation is **PARTIAL**, not a clean architecture pass.

## Fixed product defect

The Unified JSON export serialized absolute workstation paths in provenance/context metadata, unlike the portable `.sophiark` serializer. `export_unified_report` now applies the same existing portable serializer before writing JSON. Numerical result tables and machine keys are unchanged. A focused regression test verifies the paths are reduced to filenames while canonical fields remain intact.

## Required gate summary

| Gates | Result |
| --- | --- |
| A–B, D–F, H–J, P, S–V, W–Y, AE–AH, AI–AK, AM–AN, AP, AR, AV–AW | YES |
| C, G, K–L, O, Q–R, U, AX | PARTIAL |
| Z, AA | NO observed in audited service/bundle paths |
| AC | NO for active runtime code; portable artifact scan passed |
| E | YES for controlled non-engine imports |
| M–N | PARTIAL: shared registry covers core fields; many candidate-specific fields remain intentionally local/unregistered |
| AO | PARTIAL: declared environments exist, but a clean installation was not executed |
| AQ | PARTIAL: logging was observed; full encoding/consumer matrix was not executed |
| AL | PARTIAL: fixture builder idempotence is covered by tests; production builders were not run |
| AB, AD, AT, AU | YES in audited paths |

## Remaining audit debt

- Live Streamlit browser flows, including stale-state and invalid-input UX, require an interactive session execution.
- A clean virtual-environment install was not performed.
- The audit did not manually render every active table family or execute production builders, avoiding unnecessary data writes.
- `app.py` remains substantial monolithic controller debt, although the newer multipage surface uses adapters.

## Commands verified

```powershell
.\.venv\Scripts\python.exe scripts\data_health.py --fingerprint
.\.venv\Scripts\python.exe -c "from src.product.engine_freeze import verify_engine_freeze; print(len(verify_engine_freeze()))"
.\.venv\Scripts\python.exe -m pytest -q
```

Final regression: **476 passed in 24.08s**.
