# Sophiark

Sophiark is a Streamlit application for target-focused **edge-attenuation** simulations in tissue-filtered protein-protein interaction (PPI) networks. It weakens edges incident to selected targets, then compares network metrics before and after the perturbation.

> Sophiark is a network-topology research tool. Delta PageRank does not represent gene/protein expression, biological activation, causality, or clinical effect.

## Capabilities

- Separate Human (*Homo sapiens*) and Mouse (*Mus musculus*) network and tissue-data paths
- Tissue-specific PPI subnetworks through stored expression filtering
- Single- and multi-target edge attenuation
- PageRank, betweenness, Network Importance (Hinterland), community, and efficiency comparisons
- Redistribution, Target Stress, predicted compensation, propagation trace, and tissue-comparison tools when their context is available
- GO/KEGG enrichment, TRRUST/OmniPath evidence layers, exports, and reproducibility metadata

## Workflow

1. Load a STRING-based PPI network from a local SQLite snapshot.
2. Filter/weight it with the selected tissue context.
3. Calculate baseline graph metrics.
4. Reduce weights of edges connected to selected target(s).
5. Recalculate metrics and report topological differences.
6. Optionally add enrichment, regulatory evidence, and local Research Context; these do not change the core PPI calculation.

## Launch

Requirements: Python 3.10+ and the local data snapshots.

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
streamlit run app.py
```

The interface opens in English. Use **System → Language** to switch to Turkish.

## Required snapshots

| Path | Purpose |
| --- | --- |
| `data/raw/hinterland_core.db` | Human PPI interactions and tissue-expression table |
| `data/raw/mouse_hinterland_core.db` | Mouse PPI interactions and tissue-expression table |
| `data/processed/ensp_with_symbols.csv` | Human protein-ID to gene-symbol mapping |
| `data/raw/10090.protein.info.v12.0.txt.gz` | Mouse STRING-ID to symbol mapping |
| `data/regulatory/omnipath/omnipath_9606.tsv` | Human OmniPath directional-evidence snapshot |
| `data/regulatory/omnipath/omnipath_10090.tsv` | Mouse OmniPath directional-evidence snapshot |
| `data/processed/structural_binding.db` | Structural/binding-ligand lookup |

## Documentation

- [Researcher guide](docs/USER_GUIDE.md): workflow, components, interpretation, limitations, glossary, and FAQ.
- [Parameter reference](docs/PARAMETERS.md): source-backed defaults, ranges, options, and practical effects.
- [Product architecture](docs/product_architecture.md): product and data-layer contract.

## Interpretation boundaries

- STRING confidence is not causality or biological-effect strength.
- Edge attenuation is not a knockout, drug dose, or biophysical inhibition model.
- Centrality differences are not expression changes or clinical outcomes.
- Enrichment does not establish pathway activation/inhibition.
- Missing directional evidence is not negative biological evidence.
