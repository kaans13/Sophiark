# Sophiark system acceptance audit

Final status: `SYSTEM_ACCEPTANCE_PARTIAL`.

The scientific freeze, production-data fingerprints, representative numerical parity, cache/reuse path, exact bundle roundtrips, and final regression passed. The acceptance remains partial because two legacy presentation surfaces cannot yet prove the `Gen` + `ENSP` contract from their source schemas without guessing identity data.

## Evidence index

| Check | Result | Evidence | Observed value |
| --- | --- | --- | --- |
| FREEZE_01 / FREEZE_02 | PASS | `outputs/system_acceptance/engine_freeze_start.txt`, `engine_freeze_final.txt` | 12/12 exact SHA-256 matches |
| DATA_01 / DATA_02 | PASS | `data_health_start.json`, `data_health_end.json`, `data_fingerprint_comparison.txt` | 9/9 fingerprints unchanged; 7 READY capabilities |
| PARITY_01 | PASS | `unified_panel/unified_validation_summary.json` | CFTR, TGFBR1, ANO2: 16,954 rows per Classic/Directed universe |
| RUNTIME_01 | PASS | `runtime_profile/profile.json` | cold 34.3843 s; warm 16.3902 s and 16.5522 s; cache hit 0.001543 s |
| BUNDLE_01 | PASS | `outputs/analysis_bundle_validation/summary.json` | 3/3 real live/bundle comparisons exact |
| PRESENTATION_01 | PASS | `presentation_contracts_final.txt` | 79 targeted contracts passed |
| ENTITY_01 | FIXED_AND_PASS | `identity_presentation_tests.txt` | 52 targeted contracts passed |
| REGRESSION_02 | PASS | `full_regression_final.txt` | 473 passed in 22.23 s |

## Scientific integrity

- Engine freeze at start and end: PASS (`12/12`).
- Protected engine modifications during this audit: none. All audit changes are presentation adapters, labels, and tests.
- Production data: PASS (`9/9` start/end fingerprints equal).
- Classic and Directed standalone/Unified parity: PASS. Maximum observed deltas were limited to floating-point scale.

| Target / tissue | Classic max delta | Directed max delta |
| --- | ---: | ---: |
| CFTR / Lung | 3.2975410596547405e-11 | 2.0604941364243956e-10 |
| TGFBR1 / Lung | 5.0903002299373945e-11 | 5.737327453403385e-11 |
| ANO2 / Lung | 9.639372221433018e-11 | 2.355099574559427e-10 |

Evidence and CORUM context remained available in every representative Unified run; the Unified validation recorded no errors.

## Runtime, cache, history, and bundle

- Same-tissue resource reuse: PASS. The cold CFTR run took 34.3843 s; warm TGFBR1 and ANO2 runs took 16.3902 s and 16.5522 s.
- Cache: PASS. CFTR cache miss took 32.7482 s and the immediately following cache hit took 0.001543 s. The profile records no measured engine boundary work for the hit.
- Bundle load: PASS. The runtime profile and the exact roundtrip script both record no engine rerun on load.
- History/index, malformed bundle, future schema, ID conflict, and v1 compatibility: PASS through `tests/test_analysis_bundle_history.py` and `tests/test_pre_ui_ux_hardening.py`.
- Headless/UI route: PASS by static audit: the UI page calls `run_unified_analysis`, renders its returned report, and saved-report opening uses `AnalysisBundleReader().load` before the same renderer. The renderer itself contains no Classic/Directed/Unified engine call.

## Presentation and identity fixes

The following presentation-only corrections were made and regression-tested:

- Central canonical identity headers: `symbol`/`Symbol` → `Gen`; `entity_id`/`gene` → `ENSP`.
- Directed and Classic-vs-Directed projections now carry adjacent `Gen` and `ENSP` fields rather than a symbol-only `Gene` field.
- Research Explorer member and export rows use `Gen` + `ENSP`.
- Active Classic table adapters and the forward response matrix display `Gen` + `ENSP` headings.
- The status registry now explicitly covers `UNSUPPORTED_BUNDLE_VERSION` and `ID_CONFLICT` in Turkish.

Shared metric formatting, percentage semantics, missing/false/unavailable distinctions, deterministic identity-first ordering, presentation non-mutation, and central label registry behavior passed the targeted presentation contracts.

## Remaining verified scope limits

1. `src/ui/biological_interpretation.py` contains candidate projections that are handed a single `gene` value rather than a separately proven preferred-symbol and canonical-ENSP pair. Adding a guessed second identity field would violate the identity contract; this is `UNVERIFIED`, not accepted as compliant.
2. The active human/mouse comparison presentation did not have a real cross-species result universe available in this run. Its organism/orthology display semantics therefore remain `UNVERIFIED` end-to-end.

These two items make the final readiness gate `PARTIAL`; they do not invalidate the passed scientific, cache, bundle, or regression gates.
