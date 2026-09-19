# OmniPath Directed & Signed Signaling — Engineering Checkpoint

This is an implementation handoff log, not user documentation.

## Phase 1 — Architecture Discovery

- **Status:** COMPLETE
- **Files changed:** `docs/omnipath_signaling_implementation.md`
- **Decisions:** The existing STRING/igraph graph remains Layer 1 and authoritative for all Sophiark mathematics. OmniPath will be an isolated, local-first Layer 2 post-processor. It will consume immutable simulation results and never mutate or augment the PPI graph.
- **Tests passed:** Read-only architecture inspection; no scientific code changed in this phase.
- **Unresolved issues:** Existing OmniPath normalization is an evidence overlay, not a traversable signed graph. Raw endpoints include atomic UniProt proteins and OmniPath `COMPLEX:` entities; complex symbols must not be split into guessed protein-level causal edges.
- **Next phase:** Normalize the existing versioned local OmniPath snapshots into a provenance-preserving data model, run dataset QA, and map atomic entities through the existing `EntityResolver`.
- **Assumptions:** Human and mouse local snapshots dated 2026-08-15 are the authoritative offline inputs for this implementation.

### Architecture map

| Concern | Existing authoritative implementation | Layer-2 integration seam |
| --- | --- | --- |
| Network construction | `src/graph_engine.py`, `src/mouse/graph_engine.py` | No integration; read-only separation |
| Perturbation/PageRank/BC/Hinterland | `src/biology_logic.py`, `src/metrics.py` and mouse equivalents | Completed report is input only |
| Efficiency/redistribution/network loss | `src/scientific/*`, `src/services/analysis_service.py` | No formula or ranking change |
| Tissue context | SQLite `tissue_expression` tables and `src/scientific/tissue.py` | Annotation only; OmniPath edges are not deleted |
| Existing OmniPath | `src/scientific/regulatory.py`, `scripts/fetch_omnipath_snapshot.py`, versioned TSV + metadata | Reuse snapshot and provenance; add graph-capable normalized layer |
| Entity resolution | `src/research/entity_resolution.py` | Reuse `EntityResolver`; preserve ambiguous/unresolved outcomes |
| Immutable result fixture | `src/research/snapshot.py`, `SimulationResultSnapshot` | Signaling post-processing consumes frozen report rows |
| Research orchestration | `src/research/service.py`, `OfflineResearchBundle` | Optional failure-isolated result attached after core bundle stages |
| Source/provenance | `ProvenanceRecord`, local source statuses, provider statuses | Preserve dataset checksum/version/resources/references per edge |
| Research observations | `src/research/observations.py` and immutable observation contracts | Add signaling-specific, limitation-bearing presentation records |
| UI/export | `src/research/ui.py` | One Directed Signaling Context table and existing CSV/XLSX export path |
| Cache/storage | Streamlit resource cache; research evidence SQLite for snapshots | Content-addressed dataset/traversal cache; no startup fetch |
| Regression | `tests/test_forward_smoke.py`, research snapshot/foundation tests, evidence tests | Reuse suite; compare core output with feature on/off |

### Discovered data

- Human snapshot: 86,189 rows, SHA-256 `0e5a94b40d07ca25f96ec88243ac96f0fd8b30967f571ac45acf462977e012a8`.
- Mouse snapshot: 31,384 rows, SHA-256 `bd8102ca7127ab9b5bb45eee973cd9e55bd8029d10ec26f753cebc7e538348e7`.
- Columns include stable OmniPath endpoints, gene symbols, direction, stimulation/inhibition and consensus flags, resources, references, and curation effort.
- Preliminary human Layer-1 mapping through the existing resolver: 7,760 / 8,855 unique displayed endpoint labels (87.63%); 32,782 / 86,189 raw rows have both endpoints mapped to Layer 1 (38.04%). The edge rate is depressed primarily by explicit OmniPath complexes such as `COMPLEX:*`; these must remain Layer-2-only rather than being expanded by guesswork. Phase 3 must report atomic and Layer-2-only coverage separately and must not hide the raw rate.

## Phase 2 — OmniPath Data Discovery & Ingestion

- **Status:** COMPLETE
- **Files changed:** `src/signaling/models.py`, `src/signaling/ingestion.py`, `src/signaling/storage.py`, `scripts/build_omnipath_signaling_cache.py`, `.gitignore`
- **Decisions:** Existing versioned TSV files remain raw authoritative snapshots. Normalized graph data is stored in separately generated SQLite caches. Application startup never fetches the network. Conflicts retain canonical sign `0`; raw row locators, resources, references and endpoint identifiers remain traceable.
- **Tests passed:** Small/corrupt dataset rejection, duplicate/conflict aggregation, SQLite round-trip, stale checksum rejection and required indexes.
- **Unresolved issues:** None for ingestion. Mouse source quality/mapping is recorded under Phase 3.
- **Next phase:** Map nodes through the existing resolver and classify Layer-2-only entities.
- **Assumptions:** Explicit OmniPath complex identifiers are biologically resolved Layer-2 entities, but are never expanded into guessed protein-level edges.

## Phase 3 — Entity Mapping

- **Status:** COMPLETE for human; PARTIAL for mouse
- **Files changed:** `src/signaling/ingestion.py`, generated `data/processed/omnipath_signaling_9606.sqlite`, generated `data/processed/omnipath_signaling_10090.sqlite`
- **Decisions:** `src.research.entity_resolution.EntityResolver` is authoritative. Mapping QA reports total raw Layer-1 coverage separately from atomic resolver-eligible coverage. Complex and stable OmniPath-only entities stay in Layer 2 and never receive Sophiark scores.
- **Tests passed:** Exact, ambiguous, unresolved and signaling-only mapping fixtures.
- **Human acceptance:** atomic entity 7,770 / 7,877 (98.64%); atomic edge 32,692 / 33,253 (98.31%), `PREFERRED`. Raw all-edge Layer-1 coverage remains 38.04% because complex-containing edges are intentionally not expanded.
- **Mouse acceptance:** atomic entity 5,698 / 7,185 (79.30%); atomic edge 20,807 / 31,384 (66.30%), `FAILURE`. Root cause: the local mouse OmniPath projection contains 1,487 endpoints that the authoritative local mouse STRING mapping cannot resolve; frequent examples include raw UniProt accessions in the gene-symbol field and nomenclature mismatches. No fuzzy or guessed ortholog mappings were introduced.
- **Unresolved issues:** Mouse signaling is correctly marked `PARTIAL` and is not production-ready under the requested 70% edge threshold. A curated UniProt→mouse ENSMUSP mapping snapshot is required to resolve this without sacrificing accuracy.
- **Next phase:** Build traversal over the accepted mapped and explicit Layer-2-only graph without touching Layer 1.
- **Assumptions:** Human is the production-accepted signaling dataset; mouse remains opt-in/failure-isolated and transparently partial.

## Phase 4 — Directed Graph

- **Status:** COMPLETE
- **Files changed:** `src/signaling/graph.py`
- **Decisions:** Separate adjacency and reverse-adjacency indexes; no STRING graph mutation. Bounded BFS provides ancestors, descendants and distance. Alias resolution refuses ambiguous node labels.
- **Tests passed:** Direction asymmetry, shortest directed distance, max-depth enforcement and cycle protection.
- **Unresolved issues:** None.
- **Next phase:** Signed bounded paths.
- **Assumptions:** Biological direction is source → target influence, not physical movement.

## Phase 5 — Signed Path Engine

- **Status:** COMPLETE
- **Files changed:** `src/signaling/graph.py`, `tests/test_omnipath_signaling.py`
- **Decisions:** Path sign is the product of deterministic edge signs. Unsigned/conflicting edges and mixed positive/negative route sets yield `uncertain`. Enumeration is hard-capped and reverse-distance-pruned.
- **Tests passed:** `+×+`, `+×-`, `-×-`, unsigned/conflicting uncertainty, multiple paths and bounded cycle behavior.
- **Unresolved issues:** Path counts are capped at 2,048 per pair by design; this is a safety bound, not a biological threshold.
- **Next phase:** Annotate immutable simulation results.
- **Assumptions:** Network-response sign is never converted into biochemical sign.

## Phase 6 — Perturbation Result Integration

- **Status:** COMPLETE
- **Files changed:** `src/signaling/analysis.py`, `src/research/config.py`, `src/research/service.py`, `app.py`
- **Decisions:** Signaling runs only after an immutable result snapshot and only when Research Explorer is explicitly opened. `safe_analyze_signaling_context` isolates missing/corrupt/mapping/path failures. Feature flag: `ResearchConfig.directed_signaling_enabled`; default false, app opt-in true.
- **Tests passed:** Per-target response annotation, feature-enabled/disabled exact core-field preservation, unavailable-layer core preservation.
- **Unresolved issues:** None for human.
- **Next phase:** Alternative routes and convergence.
- **Assumptions:** Tissue support annotates endpoints and never filters OmniPath edges.

## Phase 7 — Alternative Routes & Multi-target Analysis

- **Status:** COMPLETE
- **Files changed:** `src/signaling/analysis.py`
- **Decisions:** Alternative route candidates require bounded descendant-region overlap. Per-target rows remain separate. Convergence requires reachability from at least two targets and preserves distance/sign/evidence per target.
- **Tests passed:** Alternative bounded-region logic, per-target independence and directed convergence fixture.
- **Unresolved issues:** Labels are descriptive candidates, never resistance/compensation claims.
- **Next phase:** Research Explorer and export.
- **Assumptions:** Default maximum depth is 4.

## Phase 8 — Research Explorer / UI / Export

- **Status:** COMPLETE
- **Files changed:** `src/research/ui.py`, `tests/test_research_ui.py`, `tests/test_omnipath_signaling.py`
- **Decisions:** One `Directed Signaling Context` table is shown in Overview. Full path counts, references and provenance are exported through existing CSV ZIP/XLSX infrastructure. Directed convergence has a separate export table. Existing navigation and visual style remain unchanged.
- **Tests passed:** UI projection, source status, signaling/provenance export, full export round-trip suite.
- **Unresolved issues:** No full-graph visualization was added; rendering all OmniPath was explicitly out of scope and the compact table meets the minimal-complexity requirement.
- **Next phase:** Full regression and real analysis.
- **Assumptions:** Detailed references belong in export rather than first-paint UI.

## Phase 9 — Full Regression & Real Analysis

- **Status:** COMPLETE for human; overall feature PARTIAL because mouse Phase 3 is below threshold
- **Files changed:** `requirements.txt` (restored explicit pytest dependency), tests and this checkpoint
- **Tests passed:** `python -m pytest -q` → **356 passed**. Exact core preservation test covers PageRank, BC, Hinterland, redistribution delta, network shift/loss context, global/local efficiency and bottleneck classification with Layer 2 enabled/disabled. Existing golden fixtures were not regenerated.
- **Real analysis:** Existing human CFTR no-tissue report, 101 response rows; 86 resolved into Layer 2, 15 unresolved, 46 upstream-reachable, 0 downstream-reachable, 40 no directed path, 2 consistently signed, 44 uncertain, 0 alternative-route candidates. Example paths include `PRKACA → CFTR`, `PRKACB → TP53 → HSPA8 → CFTR`, and `GNB1 → PIK3CA → SYK → CFTR`. These are connectivity observations only.
- **Performance:** processed human cache + graph load 1.301 s; target forward/reverse traversal 0.0061 s; full 101-row annotation 0.305 s; measured process memory delta 167.8 MB. Signaling runs after, and separately from, the perturbation motor.
- **Unresolved issues:** Mouse mapping remains below the mandatory threshold. No safe local identifier source currently closes the gap.
- **Next phase:** Acquire/version a curated mouse UniProt→ENSMUSP mapping, add it to the existing resolver as an exact alias source, rebuild the mouse cache, and rerun the same acceptance suite.
- **Assumptions:** The existing CFTR production-shape report is a valid real-result sanity fixture; no hard-coded expected biological genes were used.

## Current handoff status

- **Overall:** PARTIAL
- **Completed phases:** 1, 2, 4, 5, 6, 7, 8; Phase 3 and Phase 9 complete for human.
- **Current blocker:** Mouse atomic edge mapping 66.30% is below the 70% production floor.
- **Next exact action:** Add an authoritative local mouse UniProt→ENSMUSP identifier snapshot to `EntityResolver`, rebuild `omnipath_signaling_10090.sqlite`, and require atomic edge mapping ≥70% before changing mouse status from `PARTIAL`.
- **Files to inspect first:** this checkpoint, `src/signaling/ingestion.py`, `src/research/entity_resolution.py`, `scripts/build_omnipath_signaling_cache.py`, `tests/test_omnipath_signaling.py`.

---

# Directed Calculation Engine v1

## Phase A — Calculation and perturbation discovery

- **Status:** COMPLETE
- **Files changed:** checkpoint only during discovery.
- **Decisions:** Classic Sophiark remains immutable and default. Its target-loss contract is weighted incident-edge attenuation. Because PageRank row-normalizes outgoing weights, uniformly scaling every edge leaving the target does not change that target's transition proportions. The classic effect is primarily reduced transition probability *into* the target from neighbors that have alternative edges.
- **Evidence:** `src/biology_logic.py::run_infection_simulation` scales each unique incident edge by `block_weight_fraction`, recomputes `distance=1/weight`, and runs undirected weighted PageRank before/after. The target remains present and continues to receive teleportation mass.
- **Tests:** Existing 356-test baseline was green before directed implementation.
- **Next:** Build the hybrid topology without adding OmniPath-only edges.

## Phase B — Hybrid directed graph

- **Status:** COMPLETE
- **Files changed:** `src/directed/models.py`, `src/directed/hybrid_graph.py`, `src/directed/__init__.py`.
- **Decisions:** STRING is the only calculation-edge universe. One supported OmniPath direction produces one arc; two supported directions produce two arcs; no supported direction produces two equal-weight `undirected_fallback` arcs. Sign is metadata only. Cache identity includes tissue, STRING version, OmniPath fingerprint, direction policy and fallback policy; perturbation strategy is not part of the graph cache because it does not change the baseline topology.
- **Tests:** single direction/no reverse, bidirectional fallback/equal weight, and OmniPath-only exclusion pass.
- **Next:** Validate directed PageRank and perturbation semantics.

## Phase C — Directed PageRank

- **Status:** COMPLETE
- **Files changed:** `src/directed/engine.py`.
- **Decisions:** python-igraph PRPACK is used with `directed=True`, the same damping and non-negative STRING weights as classic. PRPACK owns convergence/dangling handling and does not expose Python `max_iter`/`tolerance`; those values are retained and validated in the API/provenance, matching the classic motor's effective behavior (its `max_iter` argument is likewise not passed to PRPACK). No personalization is used.
- **Tests:** symmetric two-arc representation matches classic undirected PageRank within `rtol=1e-10, atol=1e-12`; reversing an edge changes the result.
- **Next:** Audit target-loss strategies rather than assuming edge scaling is sufficient.

## Phase D — Directed perturbation semantics decision gate

- **Status:** COMPLETE WITH EXPLICIT LIMITATION
- **Files changed:** `src/directed/engine.py`, `scripts/audit_directed_perturbation.py`, `docs/directed_perturbation_semantics_audit.json`, directed tests.
- **Problem:** Uniform outgoing-edge attenuation is canceled by weighted PageRank row normalization. Incident attenuation can also be ineffective in a chain/cycle where every upstream node has only one outgoing transition.
- **Evidence:** Four synthetic families were tested: one-way chain, competing route, alternative source and feedback cycle. Across all cases, outgoing-only was numerically identical to baseline and incident attenuation was numerically identical to incoming-only. In the competing-route graph, perturbing B or C changed target PR by approximately `-0.0571` and total L1 redistribution by approximately `0.1142`; chain/source/single-outdegree cycle cases produced zero redistribution. Node removal produced nonzero redistribution everywhere but changed the node universe and was substantially harsher.
- **Selected strategy:** `incident_edge_attenuation`. It is retained because it is the exact directed edge-level analogue of the immutable classic contract. Its effective PageRank mechanism is reduced accessibility from upstream nodes with alternative transitions; the reduced target mass then lowers absolute downstream flow even though target outgoing proportions remain unchanged.
- **Rejected/deferred:** outgoing-only (normalized away); incoming-only (PageRank-equivalent but breaks the classic incident-edge contract for future edge-sensitive calculations); node removal (different node universe and severity); custom transition leakage (would be a new non-standard propagation model). No teleport/personalization hack or magic constant was introduced.
- **Known limitation:** This is not a universal effective knockout. Pure sources and single-transition chains may show no response; target teleportation remains. Metadata now records `directed_perturbation_strategy`, `suppression_factor`, suppressed arc count and target baseline/perturbed PR.
- **Next:** Directed redistribution and BC, followed by real-network measurement of how frequently alternative transitions make attenuation effective.

## Phase E — Directed redistribution

- **Status:** COMPLETE
- **Files changed:** `src/directed/engine.py`.
- **Decisions:** `Directed Redistribution` is computed only as perturbed directed PageRank minus baseline directed PageRank, with percent change using the same baseline denominator convention. Targets are excluded from response rows but their before/after scores remain in metadata. Multi-target suppression uses one graph and a union of incident arc IDs, so shared arcs are changed once.
- **Tests:** competing-route response, simultaneous fork perturbation and final A—B/A→B calculation-chain proof pass.
- **Next:** Separate directed BC and comparison.

## Phase F — Directed BC and safe metric boundary

- **Status:** COMPLETE
- **Files changed:** `src/directed/engine.py`.
- **Decisions:** Directed BC uses `directed=True` and `distance=1/weight`, either exact or deterministically source-sampled. Directed Hinterland, efficiency, compartment bottleneck, signed PageRank and directed null significance are explicitly absent.
- **Tests:** directed report/provenance schema and graph direction tests pass.
- **Next:** Complete comparison/UI/export and real human benchmarks.

## Phase G — Classic vs directed comparison

- **Status:** COMPLETE
- **Files changed:** `src/directed/comparison.py`, `src/services/directed_analysis.py`.
- **Decisions:** Rows are ranked by absolute redistribution percent in both engines. The table preserves entity, classic response, directed response, raw response delta, both ranks, signed rank shift, Directed BC, direct target relation/support, directed distance, direct sign context, local direction-source coverage and an explicit `not_evaluated_in_calculation_engine` alternative-route marker. Validation includes Spearman rank correlation, Pearson response correlation, overlap and Jaccard at configurable K, median/mean absolute rank shift and largest increases/decreases.
- **Cache:** A bounded four-entry in-process cache reuses immutable hybrid baseline graphs. Its key contains graph object identity/shape, tissue, STRING version, OmniPath fingerprint and both direction/fallback policies. Perturbation strategy is result provenance, not a baseline graph-cache input.
- **Tests:** comparison schema, rank validation, config-sensitive cache key and safe failure preservation pass.
- **Next:** UI, Research Explorer provenance and export.

## Phase H — UI, Research Explorer and export

- **Status:** COMPLETE
- **Files changed:** `app.py`, `src/services/directed_analysis.py`.
- **Decisions:** Sidebar selector is `Classic` (default), `Directed`, `Compare`. Classic still runs through the unchanged path. Directed/Compare use the same prepared tissue graph; Compare does not build a second hybrid graph. Results are visibly labeled `Directed Redistribution` or `Classic vs Directed`. CSV export retains the comparison table; provenance export retains engine, strategy, suppression factor, dataset fingerprint, coverage, fallback count, timings and unsupported metrics. Research snapshots receive a separate `Sophiark Directed Calculation Engine` provenance record including comparison validation; the existing OmniPath Layer 2 context remains distinct.
- **Failure behavior:** Missing/corrupt cache or mouse selection returns `UNAVAILABLE`; the already-computed classic DataFrame is not mutated or discarded.
- **Scientific boundary:** UI text states that direction changes topology while activation/inhibition remains context only.
- **Next:** Full synthetic/regression validation.

## Phase I — Synthetic validation and regression

- **Status:** COMPLETE
- **Files changed:** `tests/test_directed_calculation_engine.py`, `scripts/audit_directed_perturbation.py`, `docs/directed_perturbation_semantics_audit.json`.
- **Synthetic coverage:** resolved A→B has no B→A; fallback is equal-weight bidirectional; OmniPath-only edges are absent; mixed resolved/fallback policy is preserved; symmetric directed PR equals classic; chain/source/sink/cycle behavior is measured; reverse-edge control differs; competing/fork and richer incoming+outgoing-alternative graphs change target and downstream responses; multi-target arcs are suppressed simultaneously once; result/export provenance is complete; directed failure preserves classic input exactly.
- **Critical proof:** `build_hybrid_directed_graph` creates only A→B for STRING A—B plus OmniPath A→B. `run_directed_calculation` passes that exact hybrid graph to baseline `compute_directed_pagerank`; `suppress_target_edges` copies and attenuates incident directed arc IDs; the same directed PageRank function consumes the perturbed copy; response rows are computed from `perturbed[entity] - baseline[entity]`.
- **Classic boundary:** No classic graph/PageRank/perturbation/redistribution/BC/Hinterland/efficiency/null-model formula was edited for this feature. Existing golden fixtures were not regenerated.
- **Final regression:** `python -m pytest -q` → **371 passed in 14.45 s**, zero failures (baseline before directed work: 356; continuation-session baseline: 366).
- **Next:** Human Lung benchmarks.

## Phase J — Real human Lung benchmarks

- **Status:** COMPLETE
- **Files changed:** `scripts/benchmark_directed_engine.py`, `docs/directed_engine_real_benchmark.json`.
- **Graph:** 16,955 nodes; 513,176 STRING edges; 1,013,437 directed arcs. 14,077 STRING pairs direction-resolved (2.7431%): 12,915 single-direction and 1,162 bidirectionally supported. 499,099 STRING pairs use bidirectional fallback. 5,164 unique nodes touch a direction-resolved pair; 15,069 mapped OmniPath Layer-1 pairs are absent from the Lung STRING graph and are not added; zero STRING edges are dropped.
- **Target-local coverage:** EGFR 165 / 1,909 incident arcs (8.6433%); TP53 246 / 2,872 (8.5655%); MYC 78 / 1,816 (4.2952%). One-hop resolved fractions are 3.2356%, 2.3457%, 2.6908%; two-hop fractions are approximately 1.52%. Fallback dominance therefore explains much of the high classic/directed similarity.
- **Comparison:** 

| Scenario | Spearman | Pearson | overlap@20/50/100/300 | Jaccard@50/100 | median / mean abs rank shift | mean abs response Δ (pp) | sign changes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| EGFR | 0.9747 | 0.8784 | 14 / 38 / 87 / 270 | 0.6129 / 0.7699 | 229 / 532.35 | 0.01969 | 117 |
| TP53 | 0.9274 | 0.9127 | 17 / 41 / 82 / 260 | 0.6949 / 0.6949 | 494.5 / 985.68 | 0.03704 | 268 |
| MYC | 0.9806 | 0.9453 | 17 / 43 / 93 / 279 | 0.7544 / 0.8692 | 181 / 420.24 | 0.01010 | 115 |
| TP53 + EGFR + MYC | 0.9432 | 0.9130 | 16 / 38 / 78 / 256 | 0.6129 / 0.6393 | 344 / 779.47 | 0.05708 | 271 |

- **Unforced examples:** EGFR classic/directed top is ERRFI1; largest amplified EPGN and attenuated SEMA4B. TP53 top is NOTCH2NLB; largest amplified PPP1R13L and attenuated PRELID3A. MYC classic top CDCA7L and directed top ZNF346; largest amplified MYCBP and attenuated USP36. Combined top is NOTCH2NLB; largest amplified EPGN and attenuated PRELID3A. These names are observations, not assertions or test expectations.
- **Target PR:** All targets materially lose PageRank under both models. EGFR directed `0.00082532 → 0.000009679`; TP53 `0.00117056 → 0.000010047`; MYC `0.00065639 → 0.000009504`. The combined run attenuates all three on one graph; it is not a sum of single-target outputs.
- **Layer 2 sanity:** Among top-20 directed responses, EGFR had 18 bidirectionally reachable annotations; TP53 had 12 bidirectionally reachable and 4 indirect downstream; MYC had 10 bidirectionally reachable, 5 indirect downstream, 1 upstream and 1 no-path. Nearly all bounded path signs were uncertain, so no biochemical-state claim is made. Combined bounded Layer-2 annotation remains the slow post-processing step.
- **Performance:** STRING graph build 5.037 s; hybrid build 3.267 s; shared four-source exploratory Directed BC 0.088 s. Core per-scenario directed calculation is 1.12–1.19 s, including ~0.05–0.06 s baseline PR and ~0.13–0.15 s perturbation+PR. Comparison is ~0.03 s. Single-target runs including top-20 Layer-2 sanity are 1.98–2.77 s; multi-target Layer-2 bounded path enumeration raises total to 84.0 s, while its calculation core remains 1.19 s.
- **Memory:** process delta is 566.6 MB total: 101.2 MB STRING graph, 140.9 MB loaded OmniPath dataset and 324.5 MB hybrid directed graph. This is a measured process-level delta and includes allocator effects.
- **Interpretation:** Direction has a limited incremental effect globally because 97.26% of STRING pairs fall back to symmetric arcs, but the effect is reproducible and nonzero. TP53 and the combined scenario are more topology-sensitive than MYC by rank correlation, rank turnover and response delta. High correlations are not treated as failure; they match fallback dominance.
- **Unresolved/technical debt:** Multi-target Layer-2 path annotation needs memoized/capped batch traversal; the 324.5 MB hybrid graph is material; the four-source BC is exploratory and its sampling mode must be considered when interpreting individual values. Mouse remains unavailable at 66.30% mapping.
- **Deferred:** signed calculation, directed null model/significance, Directed Hinterland, Directed Compartment Bottleneck and directed efficiency.

## Directed engine final status

- **Human Directed Calculation Engine v1:** COMPLETE.
- **Mouse:** UNAVAILABLE/experimental; the prior 66.30% exact atomic edge mapping remains below production acceptance and does not block human.
- **Final judgment:** YES — OmniPath direction materially participates in the calculation topology used by baseline PageRank, perturbation and redistribution. It is not a post-hoc annotation. The magnitude is incremental globally because only 2.7431% of Lung STRING pairs are direction-resolved; target-local coverage is higher but still fallback-dominated.
- **Boundary:** This validates direction-aware PageRank redistribution under incident-edge attenuation. It does not simulate activation/inhibition state, and the selected perturbation remains accessibility attenuation rather than a universal effective knockout.
