# Sophiark Researcher Guide

## 1. What is Sophiark?

Sophiark is a network-topology research application for examining how a modeled perturbation of one or more gene/protein targets changes the structure of a tissue-aware protein-protein interaction (PPI) network. It uses local STRING-based interaction snapshots, optionally restricts the graph with stored tissue-expression context, attenuates edges connected to selected target(s), and compares baseline and perturbed network measures.

Its core approach is graph calculation, not machine learning. It calculates quantities such as PageRank, betweenness, Network Importance (the project’s Hinterland score), communities, and efficiency before and after attenuation. In Sophiark, a **network perturbation** is a computational reduction of target-edge weights. It is not a measured molecular perturbation, gene-expression model, or proof of a biological mechanism.

## 2. What does each component do?

### Analysis settings

This active root workspace defines the organism, tissue, calculation engine, target mode, and scenario. Use it before every forward run. Confirm the displayed species, tissue, engine, target(s), and settings before interpreting output.

### Target selection

**Search gene name** accepts up to ten mapped protein/gene choices. **Select disease** queries the configured disease source and maps returned symbols to available protein identifiers. **Manual localization** runs the existing localization-based targeting path. Use gene search for explicit targets, disease selection for mapped disease genes, and manual localization only when that is the intended model input. Check the resolved identifiers and any warning that a target is absent from the active tissue graph.

### Tissue and organism context

Human and Mouse use separate data paths. A selected tissue applies the stored expression filter; `None` disables it. This determines the graph passed to downstream calculations. Inspect the tissue-filter status and ghost/absent-target messages.

### Calculation Engine

**Classic** is the primary undirected PPI perturbation engine. In Human mode, **Directed · BETA** separately applies the directional policy, **Evidence · BETA** exposes the evidence-oriented mode, and **Compare** contrasts Classic and Directed outputs. Raw response values are not averaged. Mouse exposes Classic in the root workflow. Use only engines marked available and inspect their evidence/coverage notes.

### Network Response and visualization

After a run, Sophiark presents a response matrix, target and redistribution views, signed-loss controls, and a network map. The map has relevant-gene, capped 2D Real Network, and capped 3D Real Network views. It shows graph context, not physical cellular positions or omitted nodes.

### Biological Context, functional interpretation, and Research Context

Baseline tissue expression, cellular-localization annotations, and GO/KEGG results annotate the returned network records. They do not alter the core PPI calculation. The on-demand, read-only Research Context explorer has **Overview**, **Families**, **Functions**, and **Relationships** lenses; it builds local annotation/evidence views without changing simulation ranks or candidate classes.

### Target Stress and predicted compensation

Target Stress Search, for exactly one target, tests candidate interventions and compares target/local-network response. Its modes present network-mediated effects, loaded direct TRRUST evidence, or both. Predicted Compensation reports genes with increased relative network importance after perturbing the target. Both are computational candidate analyses, not evidence of regulation, compensation, or therapeutic effect.

### Advanced analyses

Where required data are available, the root/conditional interface provides propagation trace, tissue comparison, dose-response/threshold tools, Human–Mouse comparison, exports, and history. Availability is conditional on the active engine, species, and result context.

## 3. How Sophiark works

`target input → organism/tissue network context → baseline graph metrics → target-edge attenuation → perturbed graph metrics → network-response/redistribution presentation → optional annotations/evidence`.

The selected tissue removes records without the relevant stored expression context. Sophiark calculates baseline measures, reduces incident target-edge weights to the configured remaining fraction, recalculates graph measures, and reports relative topological differences. Candidate filtering, signed-loss display, enrichment, directional evidence, and research annotations are downstream views or separately identified modes; none turns a graph result into a causal biological claim.

## 4. How to use Sophiark

1. Launch with `streamlit run app.py`. The interface opens in English; use **System → Language** for Turkish.
2. Open **New analysis**. Select Human or Mouse, then choose a tissue or leave filtering off.
3. Select an engine available for that species.
4. Under **Target**, use Search gene name, Select disease, or Manual localization. Verify the resolved target IDs.
5. Open **Advanced settings** only if a documented comparison requires it. **Number of genes to test** defaults to 100.
6. Select **Start analysis** and wait for the preparation/simulation status.
7. Read the summary, target cards, Network Response, signed-loss/redistribution views, and visualization together. Confirm that species, tissue, targets, engine, and parameters match the intended experiment.
8. Review biological context, GO/KEGG, evidence, and Research Context as annotations. Use Target Stress or Predicted Compensation only with one selected target.
9. Export relevant tables together with their context and parameter metadata.

## 5. Understanding the results

**Network importance change (%) / Delta PageRank** is a relative change in the model after attenuation. Positive means a greater relative PageRank share in the perturbed graph; negative means lower relative share. It is not expression, protein abundance, activity, or functional loss.

**Baseline Network Importance / Hinterland score** is structural context in the baseline graph; it does not alone predict response. **Betweenness** and **Compartment Bottleneck** fields represent possible graph-bridging roles and are not proof of biological necessity or druggability. **Global/local efficiency changes** describe graph connectivity under the selected model, not measured tissue physiology.

**Redistribution candidates** are records selected because their relative structural role changes after target attenuation. A high value is a larger model-relative shift, not evidence of compensatory expression. **Endpoint**, **signed channel**, and directional fields, when supplied, describe the active engine/evidence representation. Missing direction evidence is not negative biological evidence.

GO/KEGG, localization, and baseline expression are stored annotations. Enrichment does not establish pathway activation or inhibition, and baseline expression is not a measured perturbation result. Outputs are network-level hypotheses that require independent biological and experimental validation.

## 6. Example workflow

The root interface suggests **CFTR**, **GLP1R**, or **PSEN1** as examples. Select Human, choose an available tissue such as the shipped `Adipose Tissue` option, use **Search gene name**, select the available `CFTR` mapping, and run Classic with displayed defaults. Inspect the target card, response matrix, Network Importance change, and tissue/localization annotations. Do not predetermine numerical outcomes: they depend on the installed snapshot and selected settings.

## 7. Biological terminology

| Term | Plain meaning | Sophiark context |
| --- | --- | --- |
| Gene | A DNA-encoded biological unit. | A selectable target mapped to a network record. |
| Protein | A molecule produced from gene information. | Usually the PPI-network node. |
| PPI | Protein-protein interaction. | A STRING-based weighted graph relationship. |
| Node / Edge | A graph item / connection. | A protein record / PPI relation. |
| Perturbation | A change applied to a system. | Reduction of selected target-edge weights. |
| Bottleneck | A bridge-like network position. | Structural/contextual candidate, not causality. |
| Redistribution | Relative reallocation of network role. | Structural change after attenuation. |
| Endpoint | A terminal/result-relevant graph entity. | Engine-specific response/evidence field when present. |
| Reachability | Whether a graph path connects nodes. | A graph property in relevant analysis views. |
| Biological context | Supporting entity annotation. | Stored tissue, localization, functional, or evidence information. |
| Pathway / GO / KEGG | Organized biological-process knowledge. | Interpretation/enrichment context, not demonstrated activation. |

## 8. Limitations

Sophiark calculates changes in a selected snapshot-based network model. It does not measure expression, protein activity, drug exposure, binding kinetics, cellular viability, or clinical outcome. Network coverage, interaction confidence, tissue records, identifier mapping, and directional evidence are limited by installed snapshots. Edge attenuation is not an experimental knockout or dose model. Results vary with tissue, species, snapshot, engine, targets, and settings; compare only compatible runs. Experimental validation remains necessary for mechanism and causality claims.

## 9. FAQ

**Is Sophiark machine learning?** No; the active workflow is graph calculation and source-backed annotation.

**What can I perturb?** Up to ten mapped gene/protein targets in root gene search, mapped disease-derived targets, or the available manual-localization path.

**What does a high redistribution value mean?** A larger relative structural change in this model, not proof of compensation or activation.

**Does Sophiark prove a biological mechanism?** No. It supplies network hypotheses and context.

**Can it replace experimental validation?** No.

**Who is it for?** Researchers investigating target-focused network-perturbation hypotheses in supplied Human or Mouse PPI/tissue contexts.
