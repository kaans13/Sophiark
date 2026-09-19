# Sophiark desktop scientific workspace · 2026-09-05

Authoritative entry: `app.py`, launched in the existing project `.venv`, from this project directory. Baseline: 488 tests passed in 30.58 s. No dependencies installed or changed.

## Discovery and preservation

The pre-implementation, machine-readable inventory is `outputs/ui_reconstruction/coverage-manifest.json`: 745 source UI calls, including conditional renderers, inputs, outputs, explanatory text, errors and internal unregistered pages. `source-baseline.json` records 188 Python source hashes; `before/` preserves the original presentation sources. Source inventory does not imply browser verification.

Root workflows: Human/Mouse, TR/EN, tissue, Classic/Directed/Evidence/Compare, single/multiple genes, disease query and localization; tissue expression lookup; forward execution; target stress; compensation; threshold sweep; advanced parameters; diagnostics and recent history. Results: target and redistribution readout/cards; directed and evidence comparisons; 2D/3D network; signed loss selection; network matrix and candidate detail; baseline tissue expression and cellular localization; statistical support; GO/KEGG; interpretation; Research Explorer and disease reference; complete results/export; cross-species comparison. Every conditional branch remains in the inventory.

Internal page sources additionally contained Unified, bundle history/import, dose-response, propagation, tissue comparison and sensitivity that were not registered by the original application. The reconstructed host now exposes those existing implementations through the workspace selector and advanced-results tabs. Availability is separate from browser verification; no scientific implementation is replaced.

## Global interface rules

- One wide desktop workspace; persistent configuration rail, compact context header, six result destinations. No framework migration or multipage navigation.
- Reading order: observed response → candidate and supporting evidence → biological context → full data → further investigation.
- Use existing dark palette consistently, local system font fallback, readable 14 px body, 20 px section headings, restrained numeric readouts, 24–32 px main gutters. No external font request required.
- Configuration precedes execution. Target count and scientific parameters are available before the run. Secondary exploration and maintenance are grouped separately.
- Results have explicit context and stale-state notices. Values, defaults, scientific classifications and source order are unchanged.
- Candidate evidence joins only existing rows by exact protein identity. Enrichment memberships match exact gene tokens. Missing evidence is stated as unavailable, never inferred from absence.
- Computed response, baseline annotation, statistical support and hypothesis are labelled separately. PageRank does not imply expression or activity. No composite score or new ranking.
- All-results inspector preserves full rows/columns, source precision and full-data export. Search affects the displayed subset only. Technical field groups supplement, never replace, the complete table.
- Tabs and expanders provide navigation without triggering analysis. Candidate and table controls use Streamlit fragments to avoid re-entering the scientific execution path. No mutable graph is shared with an engine as a performance shortcut.

## Acceptance contract

Requirements 1–22 and scientific UX requirements in the two supplied design briefs remain binding. Added mandatory requirements:

23. Tests executed in the user's actual local Sophiark environment.
24. Actual local Streamlit application launched successfully after all changes.
25. Major browser workflows executed against that instance.
26. No critical local-only runtime failure remains.
27. No significant performance regression from reconstruction.
28. No unnecessary scientific recomputation from presentation interactions.
29. No alternative environment used as a substitute for final acceptance.

PASS requires verified evidence for all applicable requirements. UNVERIFIED remains PARTIAL; critical broken functionality is FAIL. Unit/AppTest checks supplement, never replace, actual browser acceptance. Record exact dependency/environment changes, source preservation, test results and browser limitations in the final report.

## Local performance finding and process boundary

The actual Lung/CFTR dose workflow reached the existing 3000-vertex congestion betweenness step and remained there for several minutes. A second browser navigation timed out while the original synchronous engine call occupied the local server. The run was interrupted by restarting only the verification server launched for this task; no dependency, cache, data source or scientific code was changed. This attempt does not establish dose acceptance.

`src/ui/dose_job.py` now runs that same call in a child process using `sys.executable`, the real project directory, serialized copies of the actual graph/scores, original arguments and Python/NumPy random state. No alternate environment or reduced dose series is used. The result is written atomically and shown only in the matching context. UI polling reads completion only, never re-executes the engine. Browser verification of this boundary and a completed scientific result are separate acceptance requirements.

The local browser run confirmed that this process boundary keeps the main workspace responsive while the unchanged dose engine executes. The worker reached the existing congestion step and was stopped from the UI; a completed full dose series was not obtained in the acceptance window, so that branch remains PARTIAL.

## Final local workflow evidence

- Fresh empty state: no target, disabled run action, question-first guidance and explicit target context.
- Classic, CFTR in Lung: 21 result rows (target plus 20 candidates), all six result destinations, full signed-response data and exports.
- Directed, CFTR in Lung: AVAILABLE, 20 candidates, separate signed directed metrics and comparison output.
- Unified, CFTR in Lung: COMPLETE, Classic complete, Directed available, 16,954 matched entities, 98.1% sign agreement, 46 robust rows; a 2,760,792-byte analysis bundle was written and indexed in local history.
- Research Explorer: 646 observations, 31 family observations, 493 functional observations, one bottleneck overlap and 677 relationships. Full export construction is explicit and snapshot-cached.
- Propagation trace: 40 affected nodes, four directed-evidence edges, full table and CSV; its wording identifies a predicted network propagation rather than biochemical causality.
- Presentation interactions produced zero additional startup, network preparation, Classic simulation, enrichment or presentation executions in the captured performance log.
- Concurrent preparation requests for the same context are coalesced; callers receive detached graph copies. This changes runtime coordination only and leaves the graph builder and calculations untouched.

Evidence and Compare calculation runs, every severity/tissue-differential combination, and a completed default dose-response series were not all exercised end to end in the browser. They remain available and covered by the local automated suite, but final browser acceptance is PARTIAL under the contract's strict rule.
