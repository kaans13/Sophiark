# SOPHIARK final live execution acceptance

Final result: **FULL_SYSTEM_EXECUTION_PARTIAL**.

All active scientific engines were run on real production data without changing frozen engine code, formulas, graph weights, or data fingerprints. The official product command `streamlit run app.py` started successfully at `http://127.0.0.1:8506`. `src/ui/pages/01_analiz.py` is an internal page module and is recorded as `NOT_SUPPORTED_ENTRYPOINT`, not as a product failure.

## Accepted live execution

| Area | Result | Live evidence |
| --- | --- | --- |
| Classic / Directed / Unified | PASS | CFTR, TGFBR1, ANO2: each 16,954 Classic and 16,954 Directed rows; all Unified reports COMPLETE |
| Multi-target | PASS | CFTR+TGFBR1 resolved, 22 rows, no ghost entities |
| Full dose response | PARTIAL | Six configured doses executed with Hill n=2 and valid six-row CSV; only optional PNG creation failed |
| Threshold / Null-FDR | PASS | Six thresholds; configured three null iterations with graph immutability preserved |
| Tissue | PASS | Real Lung and Liver CFTR runs; 11 critical rows total |
| Evidence Beta / CORUM | PASS | COMPLETE, 16,954 signed rows, three CFTR complex rows, provenance embedded |
| Research Explorer | PASS | Snapshot with 77 observations, 139 relationships, and honest source notices |
| Bundles | PASS | Save/open/rebuild preserved canonical tables and scientific provenance; 3 valid bundles |
| Mouse-native / cross-species | PASS at engine level | Eight real mouse-native targets on a distinct 17,627-node mouse graph; TP53 human–mouse direction concordance 95.7% across 12,855 one-to-one nodes |
| Docking | NOT TESTED | Explicitly excluded from scope |

## Explicit limitations

The dose-response calculation itself is complete. Its optional `doz_yanit_egrisi.png` was not written because the frozen runtime tries to save under a missing nested `reports` directory. The validated CSV contains all six dose rows and required finite numerical fields. This was not changed because the responsible runtime path is protected by the engine freeze.

Browser-only user journeys, downloads, and UI state leakage are **UNVERIFIED**, not failed or guessed: the official Streamlit server was live, but this environment exposed no controllable browser tab. The execution therefore cannot honestly claim end-user click-through acceptance.

## Integrity and regression

Engine freeze verification had 12 exact matches; all nine production data resources retained their expected SHA-256 fingerprints. Real-run graph weights retained digest `a16f8db90ac20daf90fa923af72e1ce663bdb41dba61f6694214463b6b4b844e`. Full regression: **485 passed**.

Machine-readable evidence is in `outputs/final_execution_acceptance/`, including the exact mode matrix, simulation results, export/state reports, and M1–M14 mouse capability record.
