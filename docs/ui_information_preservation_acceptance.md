# SOPHIARK UI information-preservation acceptance

Status: **UI_INFORMATION_PRESERVATION_BLOCKED**.

The audit loaded the existing `CFTR_Lung.sophiark` and `TGFBR1_Lung.sophiark` fixtures into the current comparison page without running an engine.  The live report showed identity (`target`, `tissue`, completion state), cross-engine agreement, evidence provenance, CORUM context, and the Classic/Directed raw panels.  CFTR's canonical bundle contains 100 candidates, 16,954 Classic rows, 16,954 Directed rows, and 3 CORUM rows.  The UI reported the matching totals and exposed full CSV/XLSX downloads.

This is not a READY result.  The shared table renderer displays `df.head(row_limit)`: 500 rows for raw reports and 100 for Unified finding/provenance/context views.  Its caption makes the limit explicit, so it is not a silent truncation, but the remaining raw rows are accessible only by export.  That fails the requested interactive raw-data condition.  The Unified report also renders only Robust, Direction Sensitive, Classic Dominant, and Directed Emergent candidate subsets; it has no full-candidate table for `LOW_AGREEMENT` rows.

The invalid-bundle stale-state contract is also broken in `src/ui/pages/03_karsilastirma.py`: a load error is caught but `saved_unified_result` is not cleared.  Thus a valid B followed by invalid C retains B on screen until a successful D overwrites it.  The browser session confirmed a valid TGFBR1 load; its file-input adapter could not set the malformed C/D files, so the exact A → B → C → D browser run is recorded as incomplete rather than claimed.  The stale-state failure is nevertheless deterministic from the executed UI control flow.

Directly launching `src/ui/pages/01_analiz.py` without a parent `PYTHONPATH` failed with `ModuleNotFoundError: src`.  The documented entrypoint remains `app.py`; the redesigned page files are conditional rather than registered under that route.

No protected source/data was changed.  The engine-freeze check passed for all 12 protected files at both audit start and finish; the full regression suite passed 479 tests.

Supporting machine-readable evidence is in `outputs/ui_acceptance/`.  The broader table inventory used as the baseline is `outputs/ui_inventory/table_inventory.json`.
