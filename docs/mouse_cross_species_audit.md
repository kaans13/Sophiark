# Mouse Audit

Status: implementation and validation checkpoint. Human calculation code is the immutable reference.

The critical pre-implementation gate passed for the default `score >= 500` Mouse Classic path: the local dataset is taxon 10090, IDs are `ENSMUSP`, `combined_score` is read from the correct column and scale, graph weights follow the Human formula, the calculation/cache namespaces are species-separated, and Mouse provenance reports `Mus musculus`.

Two pre-existing Mouse-only issues were found. The raw STRING snapshot contains both orientations of each undirected pair; graph `simplify(max weight)` already made PageRank topology correct, but degree SQL counted both orientations. Mouse degree SQL now counts distinct neighbors. The Mouse SQLite build retained only rows at or above 500, so rescue/sweep thresholds below 500 cannot access the raw lower-confidence rows. The latter remains a documented limitation; Human was not changed.

# Dataset

Source: local STRING v12.0 `10090.protein.links.v12.0.txt`; taxon 10090.

| Metric | Mouse raw | Mouse selected ≥500 | Human selected ≥500 |
|---|---:|---:|---:|
| Rows | 12,684,354 | 1,017,864 | 569,217 |
| Unique proteins | 21,645 | 20,495 | 19,220 |
| Unique undirected pairs | 6,342,177 | 508,932 | 569,217 |
| Duplicate undirected rows | 6,342,177 | 508,932 | 0 |
| Self loops | 0 | 0 | 0 |
| Malformed IDs | 0 | 0 | 0 |
| Missing score | 0 | 0 | 0 |
| Score min/max | 150 / 999 | 500 / 999 | 500 / 999 |
| Score mean/median | 264.199 / 207 | 689.149 / 656 | 698.366 / 657 |

All-proteome Mouse graph: 20,495 nodes, 508,932 edges, 53 components, largest-component ratio 99.395%, median degree 21, mean degree 49.664. Lung Mouse graph: 17,627 nodes, 427,169 edges, 44 components, largest-component ratio 99.444%.

# Combined Score

Both species use:

`base_weight = (combined_score / 1000.0) ** 2`

The tissue multiplier is then applied and distance is `1 / final_weight`. There is no Mouse-specific score multiplier, clipping rule, or second normalization. Synthetic parity covers scores 500, 777 and 999.

# Algorithm Parity

| Component | Human | Mouse | Contract |
|---|---|---|---|
| PageRank | igraph, undirected, weighted, damping .85 | same | SAME |
| Perturbation | unique incident edge × .001 | same | SAME |
| Redistribution | shared `classify_redistribution` | same | SAME |
| Candidate selection | shared Top-N/null-FDR contract | same | SAME |
| BC | weighted inverse-strength distance | same | SAME |
| Hinterland | PR rank 60% + degree rank 40% | same; degree duplicate fixed | SAME |
| Efficiency | shared scientific implementation | same | SAME |
| Null/random seed | shared degree-matched null, seed 42 | same | SAME |

An exact synthetic Human/Mouse calculation test verifies PageRank, perturbation percentage and positive/negative response selection after only the identifier prefix changes.

# Alias Resolution

Mouse direct entity data is STRING `protein.info`. Exact symbols are case-insensitive. Ensembl Compara one-to-one relationships may serve as explicit, provenance-bearing cross-species aliases; they are not fuzzy matches. Human `ENSP` identifiers are rejected in Mouse resolution.

Examples: CFTR/Cftr → Cftr/ENSMUSP00000049228; TP53/Trp53/tp53 → Trp53/ENSMUSP00000104298 through Ensembl one-to-one mapping for the Human-style names; EGFR → Egfr/ENSMUSP00000020329; MYC → Myc/ENSMUSP00000022971.

# Species Isolation

Human IDs are `ENSP`, Mouse IDs are `ENSMUSP`. Mouse DB, output, symbol, MyGene and graph caches use separate paths. Directed cache includes taxon. Human ENSP input is never accepted as a Mouse graph node.

# Orthology

Authoritative snapshot: Ensembl Compara through the pinned Ensembl BioMart June 2026 archive. SHA-256: `7afbb340eb3637b1a6a25e513cb5e0243428a87d03d89a09e0b21dadba7ca58b`.

The snapshot has 41,859 Human symbols and 55,572 relationships: 17,039 one-to-one, 1,560 one-to-many, 566 many-to-many, and 24,063 no-ortholog symbol records. Correlations use only symbol-level bijective one-to-one projections; ambiguous/non-bijective mappings are reported and excluded rather than collapsed.

# Functional Proxy

Functional proxy discovery is an explanatory annotation-overlap ranking. Family, pathway, GO-MF, GO-BP, compartment and tissue-support components remain separate. It produces no equivalence probability, does not use perturbation output, and `automatically_selected` is always false.

# Mouse Classic

Real Lung runs succeeded for Cftr, Trp53, Egfr, Myc, Pten, Stat3, Brca1 and Nfkb1. Full results are stored in `outputs_mouse/audits/mouse_cross_species_benchmark.json`.

# Mouse Directed

UNAVAILABLE for production. Re-audited local mapping: 5,698/7,185 atomic entities (79.30%) and 20,807/31,384 atomic edges (66.30%). The project edge threshold is 70%; no fuzzy mappings were introduced. UI offers Mouse Classic only and the Directed service fails closed for taxon 10090.

# Human–Mouse Comparison

The comparison requires equal STRING threshold, damping, attenuation, selection mode and algorithm version. It projects full signed responses through bijective one-to-one orthologs, calculates species-relative rank percentiles, direction concordance and Top-20/50/100 overlaps, and keeps raw ΔPageRank descriptive.

TP53/Trp53 Lung example: 12,855 comparable one-to-one responses, 2 non-bijective symbol mappings excluded, 11,713 concordant positive, 591 concordant negative, 551 discordant, direction concordance 95.71%, rank-percentile Spearman 0.7554, Top-20/50/100 overlap 4/18/39.

# Cross-Species Mathematical Limitations

Network size, density, degree, tissue expression scale and annotation coverage differ. A Mouse ΔPageRank twice a Human value is not interpreted as twice the biological effect. One-to-many/many-to-many mappings are not arbitrarily reduced for correlation. Functional proxies are not orthologs or validated replacements.

# Tests

Baseline before this work: 385 collected tests. Final: `398 passed, 29 subtests passed in 18.63s`. Cross-species tests cover parsing contracts, score/edge parity, calculation parity, BC/Hinterland parity, alias/species isolation, orthology classes, suitability, proxy safety, config mismatch, rank percentiles, direction concordance, Mouse UI exclusions and Mouse Directed fail-closed behavior.

# Performance

Human Lung + Mouse lung network preparation in the reproducible benchmark took 25.50 seconds total. Individual Mouse target simulations took approximately 1.37–1.49 seconds with sampled structural metrics.

# Remaining Limitations

- Mouse SQLite contains only `combined_score >= 500`; lower-threshold rescue/sweep requires a rebuilt, versioned full-score DB.
- Human HPA and Mouse Bgee expression scores are separately normalized within species and are not biological fold-change comparators.
- Mouse Directed remains unavailable below mapping acceptance.
- Functional proxy production candidates require complete authoritative annotation inputs; absence returns no candidate rather than guessing.
- Comparison requires both full signed species snapshots and refuses configuration mismatch.
