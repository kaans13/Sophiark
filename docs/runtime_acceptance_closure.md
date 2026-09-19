# Runtime acceptance closure

## Final status

**FINAL_RUNTIME_ACCEPTANCE_PARTIAL.** The active dose-response control-flow defect was proven and repaired with one localized removal. Its complete real six-dose output run exceeded the configured soft runtime budget during the first dose’s congestion calculation, so full dose/output/export/state-isolation acceptance remains incomplete.

## Control-flow repair

- Active entry: `src.biology_logic.run_pharmacological_dose_response`, called by `src/ui/analysis_results.py`.
- Pre-fix evidence: an invocation with three supplied doses returned a sentinel mapping `dict`, with zero dose rows; see `outputs/runtime_acceptance_closure/pre_fix_control_flow_probe.json`.
- Defect: the unconditional `return` directly after the sentinel log bypassed the existing `all_names` target resolution and `for sf in sorted(dozlar, reverse=True)` loop.
- Repair: only that return dictionary was removed. The existing loop, Hill equation, configured doses, graph weighting, PageRank and output construction were not edited.
- Targeted regression: `tests/test_dose_response_control_flow.py` verifies one result per dose, correct Hill metadata, each supplied dose entering the existing Hill call, and restored graph weights/distances.

## Real dose run

The CFTR/Lung run used the existing configured series `0.90, 0.75, 0.50, 0.25, 0.10, 0.01` and Hill `n=2`. It entered the real loop and applied `0.90`; the real congestion BC computation began for 3,000 nodes. It did not finish usefully within the 120-second soft mode budget and was stopped. See `outputs/runtime_acceptance_closure/dose_response_time_budget.json`.

## Bounded parity

Standalone CFTR/Lung Classic and Directed completed in 31.05 seconds. Compared with the previous Unified CFTR artifact, maximum differences were floating-point scale only: Classic response `4.56e-11`, Directed response `5.57e-11`, Directed BC `0.0`. The post-fix Unified synthesis used these standalone results and had exact zero deltas, with Evidence Provenance and CORUM both available. See `outputs/runtime_acceptance_closure/bounded_parity.json`.

## Safety

- Only `src/biology_logic.py` changed among protected files; its old SHA-256 was `d2ae2392ab11699aac3d9ecbbc640dad65de12dbe0178d7f311ada4f98017722` and new SHA-256 is `f831ab35bfa83e43a6094e75d78b8ad3832c2dd3ef30c6a9af00fbae556aaa73`.
- `config/engine-freeze.json` was updated only for that documented repair; freeze verification is 12/12 exact.
- Production capabilities and fingerprints are unchanged and remain READY.
- Final suite: **475 passed in 23.82s**.
