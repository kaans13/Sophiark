# Sophiark Parameter Reference

This reference is source-backed. `app.py` is the documented active entry point; the **conditional page layer** (`src/ui/pages`) is included separately because it contains real configurable controls but is not registered by the documented root launch. Values below are interface defaults/ranges, not biological recommendations.

## Summary

| Parameter | Default | Range/options | Controls | Practical effect |
| --- | --- | --- | --- | --- |
| Organism | Human | Human, Mouse | data/engine path | changes source network and tissue set |
| Tissue | first Human selectable tissue (`Adipose Tissue`) | supplied tissue list or `None` | active graph context | filters the graph by stored expression context |
| Engine | Classic | Classic; Human also Directed, Evidence, Compare | result flow | selects the presented calculation engine/mode |
| Targets | none | 1–10 in root gene search | attenuated incident edges | changes the modeled perturbation set |
| Number of genes to test | **100** | 20, 30, 50, 100, 200 | candidate scan scope | broader scan can take longer |
| Damping (`damping`) | 0.85 | 0.70–0.99, 0.01 | PageRank propagation | changes global/local propagation balance |
| Remaining edge weight (`block_strength`) | 0.001 | 0.0001–0.1, 0.0001 | attenuation severity | larger fraction retains more target-edge weight |
| Hill coefficient (`hill_n`) | 2.0 | 0.5–4.0, 0.5 | dose-response curve | changes curve steepness in associated tools |
| Paralog boost (`paralog_boost`) | 10% | 0–50%, 5% | paralog contribution model | raises modeled paralog boost |
| Betweenness accuracy (`bc_sample_sources`) | Fast / 1200 | Fast=1200, Balanced=3000, Full=None | BC computation scope | more/full sources increase computation |
| Real Network nodes (`harita_max_node`) | 600 | 150–1500, 50 | map display cap | changes displayed context, not result calculation |
| Real Network labels (`harita_max_etiket`) | 25; Restore defaults sets 40 | 10–150, 5 | map label cap | changes readability/display only |
| Target Stress mode | Network-mediated effect | Network-mediated, Direct, Combined | candidate evidence presentation | changes the search/evidence mode |
| Target Stress limit | 5 | 5, 10, 15 | tested candidate count | more candidate forward runs take longer |
| Negative loss threshold | 0.05% | 0–100%, 0.01 | signed-loss display selection | changes displayed negative candidates only |
| Negative candidate limit | 100 | 10–2000, 10 | signed-loss display length | changes display/export selection only |
| Language (`ui_locale`) | English | English, Turkish | UI text | presentation only |

## Active root workspace parameters

### Organism — `tur_secim_radio`

**Meaning/control:** Chooses Human or Mouse and thereby the engine module, symbol map, tissue list, output directory, and local network snapshot. **Unit:** categorical. **Default:** Human. **Options:** Human (*Homo sapiens*), Mouse (*Mus musculus*). **Increase/decrease:** not numeric. **Practical interpretation:** choose the organism matching the data context; do not compare values across organisms as if they share a graph. The root exposes only Classic for Mouse.

### Tissue — `analysis_tissue`

**Meaning/control:** Selects the stored tissue-expression context used to filter the PPI graph. **Unit:** categorical. **Default:** the second Human list entry, `Adipose Tissue`, when not already in session; `None` is also available. **Options:** the shipped Human or Mouse tissue list. **Increase/decrease:** not numeric. **Practical interpretation:** a tissue selection removes genes not produced in the selected stored context; `None` turns that filter off. Use the context you intend to study, and retain it when comparing analyses.

### Calculation Engine — `calculation_engine`

**Meaning/control:** Chooses the result flow. **Default:** Classic. **Options:** Human: Classic, Directed, Evidence, Compare; Mouse: Classic. **Unit:** categorical. **Practical interpretation:** Classic is the undirected PPI engine. Directed/Evidence are marked BETA and retain their own result/evidence semantics; Compare contrasts Classic and Directed. Changing it requires a new analysis run.

### Target type and targets — `mod`, `gen_sec`

**Meaning/control:** Selects Search gene name, Select disease, or Manual localization; gene search supplies the explicit target list. **Default:** Search gene name, no targets. **Range/options:** one to ten mapped choices in root gene search; disease results depend on the configured source; manual localization offers the eight stored locations. **Unit:** identifiers/category. **Practical interpretation:** selected target edges are attenuated. More targets expand the perturbation set; shared edges are weakened once. This is a computational target choice, not a claim of biological inhibition.

### Number of genes to test — `redistribution_candidate_limit`

**Meaning/control:** Limits the candidate scan used for redistribution-oriented output. **Unit:** count. **Default:** **100**. **Allowed values:** 20, 30, 50, 100, 200. **Increase:** broadens the candidate scan and can increase runtime. **Decrease:** narrows the scan and can finish sooner. **Practical interpretation:** use the displayed default as the baseline; record the count when comparing outputs because the tested candidate universe differs.

### Network-importance propagation — `damping`

**Meaning/control:** PageRank damping/propagation setting. **Unit:** coefficient. **Default:** 0.85. **Range:** 0.70–0.99 in steps of 0.01. **Increase:** retains more propagation through the graph, producing a more global influence pattern. **Decrease:** makes the modeled propagation more local. **Practical interpretation:** alter only for a defined sensitivity comparison and preserve the value alongside the result.

### Remaining edge weight — `block_strength`

**Meaning/control:** The fraction of target-incident edge weight remaining after attenuation. **Unit:** fraction. **Default:** 0.001. **Range:** 0.0001–0.1 in steps of 0.0001. **Increase:** leaves more edge weight and models a milder attenuation. **Decrease:** leaves less edge weight and models a stronger attenuation. **Practical interpretation:** it is a graph operation, not a calibrated knockout, drug dose, or biochemical inhibition value.

### Hill coefficient — `hill_n`

**Meaning/control:** Sets the steepness of the dose-response curve used by associated analysis tools. **Unit:** coefficient. **Default:** 2.0. **Range:** 0.5–4.0, step 0.5. **Increase:** makes the modeled curve steeper. **Decrease:** makes it more gradual. **Practical interpretation:** do not treat a chosen value as biologically validated unless supported externally; it is a model setting.

### Paralog boost — `paralog_boost`

**Meaning/control:** Percentage boost applied by the existing paralog contribution model. **Unit:** percent. **Default:** 10. **Range:** 0–50, step 5. **Increase/decrease:** raises/lowers that modeled boost; zero removes this added contribution. **Practical interpretation:** it changes model behavior and is not evidence that paralogs compensate in vivo.

### Betweenness accuracy — `bc_sample_sources`

**Meaning/control:** Chooses the source-node scope for betweenness calculation. **Unit:** count/mode. **Default:** Fast = 1200. **Options:** Fast=1200, Balanced=3000, Full=`None` (all nodes). **Increase:** Balanced uses more sampled sources; Full uses all nodes and can take longer. **Decrease:** Fast uses fewer sources and is quicker. **Practical interpretation:** keep the mode consistent for comparisons of betweenness-derived output.

### Real Network display limits — `harita_max_node`, `harita_max_etiket`, `tum_etiket_goster_chk`, `harita_modu_secim`

**Meaning/control:** These control presentation of the Real Network map. **Defaults:** 600 nodes, 25 labels (the Restore defaults button sets 40 labels), all-core-label toggle off, Relevant Genes (Fast) map mode. **Ranges/options:** nodes 150–1500 in 50s; labels 10–150 in 5s; modes Relevant Genes, Real Network 2D, Real Network 3D. **Increase:** shows more context/labels and can make the browser view denser or slower. **Decrease:** reduces visual density. **Practical interpretation:** these do not recalculate the network result; core target/stressed genes are included independently of the context cap.

### Target Stress settings — `target_stress_mode`, `target_stress_limit`

**Meaning/control:** Sets whether the search presents network-mediated effect, loaded direct TRRUST evidence, or combined evidence, and how many candidates are tested. **Defaults:** Network-mediated effect; 5. **Options:** three modes; limit 5/10/15. **Unit:** category/count. **Increase of limit:** performs more per-candidate forward simulations and may take longer. **Practical interpretation:** requires exactly one target. Direct mode depends on loaded directional evidence; missing evidence is not a negative biological finding.

### Signed-loss display controls — `negative_min_abs_delta_pct`, `negative_redistribution_top_n`

**Meaning/control:** Filter and limit the negative Network Importance-change display. **Units:** percent/count. **Defaults:** 0.05 and 100. **Ranges:** 0–100 by 0.01; 10–2000 by 10. **Increase:** threshold requires larger negative changes; limit shows more qualifying records. **Decrease:** threshold admits smaller negative changes; limit shows fewer records. **Practical interpretation:** source comments state these only affect result selection/display; they do not recompute network importance or alter the positive table.

### Language — `ui_locale`

**Meaning/control:** UI presentation language. **Default:** English. **Options:** English, Turkish. **Unit:** categorical. **Practical interpretation:** language changes visible labels only, not graph construction, analysis parameters, or results.

## Conditional page-layer parameters

These controls exist in `src/ui/pages` but are not top-level routes in the documented `streamlit run app.py` entry point.

| Parameter / config | Default | Options/range | Effect |
| --- | --- | --- | --- |
| Perturbation severity / `perturbation_severity`, `block_strength` | Near-complete / 0.001 | Mild=0.50, Moderate=0.20, Strong=0.05, Near-complete=0.001 | sets remaining edge fraction |
| BC sample sources | 1200 | 300, 600, 1200 | sampled BC scope |
| Efficiency validation | Fast Approximate | Fast Approximate, Full Validation | calculation mode |
| Global/local efficiency samples | 256 / 128 | 64/128/256/512; 32/64/128/256 | sampling scope |
| Redistribution mode | `null_fdr` | `null_fdr`, `percentile`, `top_n` | candidate-selection rule |
| Null sensitivity / `null_iterations` | Standard / 99 | Fast=19, Standard=99, Rigorous=499 | empirical-null iterations |
| Exploratory Top-N | 250 | 10–1000, step 10 | `top_n` selection count |
| Absolute change percentile | 99.0 | 90.0–100.0, 0.1 | percentile selection cutoff |
| Tissue normalization | `within_tissue` | `within_tissue`, `cross_tissue_comparable` | tissue normalization mode |
| Evidence edge weighting | `structural_only` | `structural_only`, `legacy_edge_modifiers` | evidence/edge policy |
| Random seed | 42 | 0–2,147,483,647 | reproducible sampled/random steps |
| Conditional map nodes | 350 | 50–800, 50 | display cap |
| Propagation nodes/routes | 40 / 12 | 20–100, 5 / 3–30, 1 | trace display scope |
| Discovery candidate count | 5 | 3–25 | Target Stress validations |
| Discovery deep validation | Off | Off, first 3, first 5 | Full-BC validation candidates |
| Tissue comparison | none | 2–8 tissues | tissues included in comparison |

For each conditional numeric control, increasing a sample/count/iteration/node/route setting broadens the respective calculation or display scope and can increase runtime; decreasing it narrows that scope. Increasing the percentile makes selection more restrictive. The selection modes change which rule chooses candidates; they do not mean one biological interpretation is “correct.” The random seed changes the deterministic random stream, not a biological variable.

## How should I choose parameters?

Start with the displayed root defaults to establish a reproducible baseline: Classic, the intended organism/tissue, damping 0.85, remaining edge weight 0.001, Hill coefficient 2.0, paralog boost 10%, Fast betweenness, and 100 genes to test. The settings most likely to change the result context or candidate universe are organism, tissue, target set, engine, attenuation, damping, candidate limit, and candidate-selection mode in the conditional layer. Sampling settings mainly trade compute cost against the scope used by approximate metrics; map settings affect presentation only.

For robustness work, vary one setting at a time and retain the complete context. Do not directly compare candidate ranks or percentages from runs whose organism, network snapshot, tissue filter, target list, engine, attenuation, damping, selection rule, candidate limit, or sampling scope differ without stating those differences. No value in this document is labeled optimal, best, or biologically correct: the repository defines computational behavior, not experimental calibration.
