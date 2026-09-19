# Evidence Engine BETA data audit

- Human taxon: NCBI 9606.
- Production association source: local `data/raw/hinterland_core.db`, 6,857,702 unique undirected STRING v12 rows; production threshold 500; score scale 0–1000.
- Existing association schema before this work: `protein1`, `protein2`, `combined_score`. It cannot reconstruct channel-decomposed confidence.
- Existing physical source reused: `data/raw/hinterland_core_physical.db`, 584,794 STRING physical rows. HuRI was not downloaded because physical context is already locally available and the default physical multiplier is off.
- Existing complex source reused: `data/raw/coreComplexes.txt` plus the existing processed CORUM lookup. No replacement complex dataset was downloaded.
- One new resource was necessary: official human STRING v12 `9606.protein.links.full.v12.0.txt.gz`. The full variant was selected because it is the minimum v12 human resource preserving direct and transferred channel scores, including transferred text mining. The smaller detailed variant cannot implement the required transferred-text exclusion reliably.
- The generated `data/processed/string_evidence_v12_9606.sqlite` indexes only unique edges with official combined score ≥500. It contains 569,217 rows, exactly matching the production Classic edge universe at that threshold.
- STRING `homology` is retained as provenance but is not independently combined; transferred channel scores are combined explicitly. `S_nontext` excludes direct and transferred text-mining evidence.
- No STRING version upgrade, HuRI download, physical replacement, or CORUM replacement was performed.

The source file SHA-256, generated-index SHA-256, physical database SHA-256, CORUM SHA-256, tissue, Evidence mode, `S_nontext` policy, Hill parameters, physical policy, CORUM policy, perturbation settings, and engine version are embedded as machine-readable JSON in every Evidence result/export row.
