# SOPHIARK UI information-preservation closure

Status: **UI_INFORMATION_PRESERVATION_READY** after the four verified non-engine blockers were repaired.

## Raw data

`render_dataframe` now supports a positional full-data inspector. It shows a page selector, visible row range, total rows, page count and a fixed 500-row page size. It uses `iloc` only: no filtering, re-sorting, rank recalculation, percentile recalculation, symbol collapse or value transformation is applied to canonical data. Classic and Directed raw panels, plus the root raw forward report, use this inspector.

For the existing CFTR/Lung bundle, Classic and Directed each contain 16,954 rows. Tests reach the first row, row 500, row 501 and final row, then concatenate every page and compare it exactly with the canonical frame. Downloads remain available but are no longer the only way to inspect rows after 500.

## Unified candidates

The Unified report has a dedicated **Düşük uyum** lens and an unfiltered, paginated **Tüm adaylar** inspector. Both retain the canonical candidate schema and export through the existing CSV/XLSX path. A deterministic two-row LOW_AGREEMENT fixture verifies canonical count, UI-accessible count, export count, ENSP identity, response values and ranks.

## Bundle state

Before the repair, an exception from `AnalysisBundleReader().load()` was caught while the previous `saved_unified_result` remained in session state. `load_saved_unified_result` now invalidates only that transient display state before every explicit open; it assigns a result only after successful validation. The invalid-after-valid test leaves no displayed result, preserves durable history, and then successfully loads a new valid bundle.

## Entrypoint contract

The supported product command is `streamlit run app.py`, documented in `readme.md`. `src/ui/pages/01_analiz.py` is an internal workspace page and is explicitly marked `NOT_SUPPORTED_ENTRYPOINT`; no `sys.path` or `PYTHONPATH` workaround was added. The supported root command started successfully from a controlled temporary working directory without manual `PYTHONPATH`.

Final verification passed: 485 tests, engine freeze 12/12 exact at start and end, and all 9 catalogued production-data fingerprints unchanged. The evidence is recorded in `outputs/ui_acceptance/closure_summary.json` and `outputs/ui_acceptance/final_summary.json`.
