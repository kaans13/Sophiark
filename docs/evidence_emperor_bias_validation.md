# Evidence Engine BETA — Emperor / Hub Bias Validation

## Judgment: BIAS_WARNING

The current, unchanged Evidence default was tested: E1 topology, `S_nontext` edge weighting, text multiplier OFF, physical multiplier OFF, and CORUM restricted to post-calculation complex context. The panel reused the eight existing target/tissue conditions: CFTR, TP53, EGFR, and GLP1R in Lung; TP53, MYC, PTEN, and STAT3 in Liver.

Absolute-response promotion is consistently defined as `ClassicRank - EvidenceRank`; positive values indicate a higher Evidence rank. Positive and negative response ranks were audited separately in the CSV artifacts. The meaningful-response subset uses the existing `|Delta_PageRank_Pct| >= 0.05` threshold; it is reported alongside all-node results to prevent near-zero rank noise from driving the conclusion.

## Finding

In the meaningful absolute-response audit, D10 Classic-degree nodes had positive median promotion in 8/8 target conditions (aggregate median +240.25; 60.4% promoted). D10 Classic-baseline-PageRank nodes also had positive median promotion in 8/8 conditions (aggregate median +246.31; 60.9% promoted). In contrast, D1 baseline-PageRank nodes had aggregate median promotion -963.19 and only 48.1% promoted after the near-zero filter.

Promotion relationships with annotation/evidence coverage repeated across targets: median target-level Spearman rho was +0.1486 for channel coverage (7/8 positive), +0.2580 for database availability (8/8 positive), +0.2243 for experimental availability (8/8 positive), and +0.2756 for physical support (8/8 positive). These annotations are not calculation modifiers, but their consistent association with promotion remains a bias risk.

The coarse matched-control audit, stratified by Classic degree decile, Classic baseline-PageRank decile, and tissue, found positive high-annotation-minus-low-annotation median promotion for all 8 target conditions. This does not prove causality, but it makes a purely hub-structure explanation insufficient.

`S_nontext` itself remains centrality-associated. Mean incident `S_nontext` versus Classic degree/PageRank Spearman rho was approximately +0.42 in both tissues. Mean `S_nontext / official combined_score` increased from 0.42–0.44 in D1 to 0.68–0.69 in D10, meaning removal of text mining does not remove the centrality advantage in retained non-text support.

There were 397 Classic-low to Evidence-high absolute Top-100 transitions. Only 5 met the strict Emperor diagnostic group, and 309 were neither top-decile degree nor top-decile annotation coverage; therefore these transitions are not simply a list of famous genes. Still, their evidence profiles are heterogeneous, so they do not negate the global inflation signal.

No Evidence, Classic, or Directed calculation code was changed. The only corrected defects were in the audit: excluding the directly perturbed target from its otherwise signed-response universe, and removing duplicate annotation/matched-control summary rows. Full regression passed afterward.

## Artifacts

The machine-readable artifacts are under `outputs/evidence_validation/emperor_bias_final/`:

- `evidence_emperor_bias_summary.json`
- `evidence_degree_decile_audit.csv`
- `evidence_pagerank_decile_audit.csv`
- `evidence_annotation_bias.csv`
- `evidence_matched_controls.csv`
- `evidence_candidate_transitions.csv`

The appropriate action for this validation task is to report the warning, not retune the Evidence model.
